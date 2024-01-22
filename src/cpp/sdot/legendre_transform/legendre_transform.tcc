#pragma once

#include "../WeightedPointSet_AABB.h"
#include "../cell_points.h"
#include "legendre_transform.h"
#include <eigen3/Eigen/LU>

#define DTP template<class Ci_nb_dims, class Ct_sval, class M_dirs, class M_offs, class B_dirs, class B_offs>
#define UTP LegendreTransform<Ci_nb_dims, Ct_sval, M_dirs, M_offs, B_dirs, B_offs>

DTP UTP::Ret UTP::transform() {
    using namespace Vfs;
    Vfs::VecImpl<Svec> nm_dirs;
    Vfs::VecImpl<Sval> nm_offs;
    Vfs::VecImpl<Svec> nb_dirs;
    Vfs::VecImpl<Sval> nb_offs;

    //
    if constexpr( nb_dims == 0 ) {
        nm_offs.append( m_offs );
        return { nm_dirs, nm_offs, nb_dirs, nb_offs };
    } else {
        // if there's an equality constraint, we're actually working at most on nb_dims - 1
        if ( Opt<std::pair<Svec,Svec>> p = first_eq_bnd() )
            return transform_without_dir( p->first, p->second, false );

        // else, get the power diagram vertices (needed for update of used_ms and used_bs)
        const auto coords_and_cuts = vertex_corr();
        const auto &vertex_coords = std::get<0>( coords_and_cuts );
        const auto &vertex_cuts = std::get<1>( coords_and_cuts );

        // update of used_ms and used_bs (needed for unused_dirs)
        update_used_ms_and_bs( vertex_coords, vertex_cuts );

        // if already lies in a sub-space...
        if ( Opt<std::pair<Svec,Svec>> p = unused_dir() )
            return transform_without_dir( p->first, p->second, true );

        // make the new affine functions
        make_new_affs( nm_dirs, nm_offs, vertex_coords, vertex_cuts);

        // make the new boundaries
        make_new_bnds( nb_dirs, nb_offs, vertex_coords, vertex_cuts );

        return { nm_dirs, nm_offs, nb_dirs, nb_offs };
    }
}

DTP auto UTP::vertex_corr() {
    using namespace Vfs;

    // acceleration structure
    VecImpl<Sval> weights( FromReservationSize(), m_offs.size() );
    for( PI i = 0; i < m_offs.size(); ++i )
        weights << norm_2_p2( m_dirs[ i ] ) - 2 * m_offs[ i ];
    auto wps = make_weighted_point_set_aabb( m_dirs, weights, ci_nb_dims );

    // correspondance
    return cell_points( wps, b_dirs, b_offs );
}

DTP void UTP::update_used_ms_and_bs( const auto &vertex_coords, const auto &vertex_cuts ) {
    used_ms = std::vector<bool>( m_dirs.size(), false );
    used_bs = std::vector<bool>( b_dirs.size(), false );
    for( PI num_vertex = 0, nb_vertices = vertex_coords.size() / nb_dims; num_vertex < nb_vertices; ++num_vertex ) {
        const auto *cuts = vertex_cuts.data() + num_vertex * ( nb_dims + 1 );
        for( PI r = 0; r <= nb_dims; ++r ) {
            PI cut = cuts[ r ];
            if ( cut < m_dirs.size() ) {
                used_ms[ cut ] = true;
                continue;
            }

            cut -= m_dirs.size();
            if ( cut < b_dirs.size() ) {
                used_bs[ cut ] = true;
                continue;
            }
        }
    }
}

DTP Vfs::Opt<std::pair<typename UTP::Svec,typename UTP::Svec>> UTP::unused_dir() {
    using namespace Vfs;

    // starting point
    const PI not_used = -1;
    PI first_used_m = not_used;
    for( PI n = 0; n < used_ms.size(); ++n ) {
        if ( used_ms[ n ] ) {
            first_used_m = n;
            break;
        }
    }

    // covariance matrix
    Svec init_dir( FromItemValue(), 0 );
    if ( first_used_m != not_used )
        init_dir = m_dirs[ first_used_m ];

    using CovMat = Eigen::Matrix<Sval,nb_dims,nb_dims>;
    CovMat cov_mat;
    for( PI r = 0; r < nb_dims; ++r )
        for( PI c = 0; c < nb_dims; ++c )
            cov_mat.coeffRef( r, c ) = 0;
    for( PI n = first_used_m + 1; n < used_ms.size(); ++n ) {
        if ( ! used_ms[ n ] )
            continue;
        Svec dir = Svec( m_dirs[ n ] ) - init_dir;
        Sval cdi = norm_2_p2( dir );
        for( PI r = 0; r < nb_dims; ++r )
            for( PI c = 0; c < nb_dims; ++c )
                cov_mat.coeffRef( r, c ) += dir[ r ] * dir[ c ] / cdi;
    }
    for( PI n = 0; n < used_bs.size(); ++n ) {
        if ( ! used_bs[ n ] )
            continue;
        Svec dir = b_dirs[ n ];
        Sval cdi = norm_2_p2( dir );
        for( PI r = 0; r < nb_dims; ++r )
            for( PI c = 0; c < nb_dims; ++c )
                cov_mat.coeffRef( r, c ) += dir[ r ] * dir[ c ] / cdi;
    }

    Eigen::FullPivLU<CovMat> lu( cov_mat );
    if ( lu.dimensionOfKernel() )
        return std::make_pair( init_dir, Svec( lu.kernel().col( 0 ) ) );

    return {};
}

DTP UTP::Ret UTP::transform_without_dir( Svec pos, Svec dir, bool add_bnd ) {
    using namespace Vfs;

    // normalization of dir
    dir /= norm_2( dir );

    // a simple grahm shmidt to find a base
    VecImpl<VecImpl<Sval,nb_dims>,nb_dims-1> base;
    PI dir_to_avoid = argmax( abs( dir ) );
    for( PI i = 0, j = 0; i < nb_dims; ++i ) {
        if ( i != dir_to_avoid ) {
            Svec b{ FromItemValue(), 0 };
            b[ i ] = 1;

            b = b - sp( b, dir ) * dir;
            for( PI k = 0; k < j; ++k )
                b = b - sp( b, base[ k ] ) * base[ k ];
            b /= norm_2( b );

            base[ j++ ] = b;
        }
    }

    // P( base, pos );

    // make new inputs from projection
    using Pvec = VecImpl<Sval,nb_dims-1>;
    Vfs::VecImpl<Pvec> pm_dirs( FromReservationSize(), m_offs.size() );
    Vfs::VecImpl<Sval> pm_offs( FromReservationSize(), m_offs.size() );
    for( PI i = 0; i < m_offs.size(); ++i ) {
        Svec m_dir = m_dirs[ i ];

        Pvec pm_dir;
        for( PI d = 0; d < nb_dims - 1; ++d )
            pm_dir[ d ] = sp( m_dir - pos, base[ d ] );
        pm_dirs << pm_dir;

        pm_offs << m_offs[ i ] - sp( m_dir, pos );
    }

    Vfs::VecImpl<Pvec> pb_dirs( FromReservationSize(), b_offs.size() );
    Vfs::VecImpl<Sval> pb_offs( FromReservationSize(), b_offs.size() );
    for( PI i = 0; i < b_offs.size(); ++i ) {
        Svec b_dir = b_dirs[ i ];
        if ( all( b_dir == dir ) || all( b_dir == -dir ) )
            continue;

        Pvec pb_dir;
        for( PI d = 0; d < nb_dims - 1; ++d )
            pb_dir[ d ] = sp( b_dir, base[ d ] );
        pb_dirs << pb_dir;

        pb_offs << b_offs[ i ] - sp( b_dir, pos );
    }

    // make a legendre transform with the new base
    ::LegendreTransform pnlt( Vfs::CtInt<nb_dims-1>(), ct_sval, pm_dirs, pm_offs, pb_dirs, pb_offs );
    auto plt = pnlt.transform();

    Vfs::VecImpl<Pvec> &nm_dirs = std::get<0>( plt );
    Vfs::VecImpl<Sval> &nm_offs = std::get<1>( plt );
    Vfs::VecImpl<Pvec> &nb_dirs = std::get<2>( plt );
    Vfs::VecImpl<Sval> &nb_offs = std::get<3>( plt );

    // inverse projection
    Vfs::VecImpl<Svec> im_dirs( FromReservationSize(), nm_offs.size() );
    Vfs::VecImpl<Sval> im_offs( FromReservationSize(), nm_offs.size() );
    for( PI i = 0; i < nm_offs.size(); ++i ) {
        Pvec nm_dir = nm_dirs[ i ];
        Sval im_off = nm_offs[ i ];
        Svec im_dir;
        for( PI e = 0; e < nb_dims; ++e ) {
            Sval v = pos[ e ];
            for( PI d = 0; d < nb_dims - 1; ++d ) {
                v += nm_dir[ d ] * base[ d ][ e ];
                im_off += nm_dir[ d ] * base[ d ][ e ] * pos[ e ];
            }
            im_dir[ e ] = v;
        }

        im_dirs << im_dir;
        im_offs << im_off;
    }

    Vfs::VecImpl<Svec> ib_dirs( FromReservationSize(), nb_offs.size() );
    Vfs::VecImpl<Sval> ib_offs( FromReservationSize(), nb_offs.size() );
    for( PI i = 0; i < nb_offs.size(); ++i ) {
        Pvec nb_dir = nb_dirs[ i ];
        Sval ib_off = nb_offs[ i ];
        Svec ib_dir;
        for( PI e = 0; e < nb_dims; ++e ) {
            Sval v = 0;
            for( PI d = 0; d < nb_dims - 1; ++d ) {
                v += nb_dir[ d ] * base[ d ][ e ];
                ib_off += nb_dir[ d ] * base[ d ][ e ] * pos[ e ];
            }
            ib_dir[ e ] = v;
        }

        ib_dirs << ib_dir;
        ib_offs << ib_off;
    }

    // add equalities for kernel vectors
    if ( add_bnd ) {
        const auto s = sp( dir, pos );

        ib_dirs << dir;
        ib_offs << s;

        ib_dirs << -dir;
        ib_offs << -s;
    }

    return { im_dirs, im_offs, ib_dirs, ib_offs };
}

DTP Vfs::Opt<std::pair<typename UTP::Svec,typename UTP::Svec>> UTP::first_eq_bnd() {
    for( PI i = 0; i < b_dirs.size(); ++i ) {
        Svec bdi = b_dirs[ i ];
        for( PI j = 0; j < i; ++j ) {
            Svec bdj = b_dirs[ j ];
            if ( all( bdi == -bdj ) && b_offs[ i ] == -b_offs[ j ] )
                return std::make_pair( b_offs[ i ] / Vfs::norm_2_p2( bdi ) * bdi, bdi );
        }
    }
    return {};
}

DTP Vfs::PI UTP::nb_inf_cuts( const auto *cuts ) {
    PI res = 0;
    for( PI i = 0; i <= nb_dims; ++i )
        if ( inf_cut( cuts[ i ] ) )
            ++res;
    return res;
}

DTP void UTP::make_new_affs( Vfs::VecImpl<Svec> &new_m_dirs, Vfs::VecImpl<Sval> &new_m_offs, const auto &vertex_coords, const auto &vertex_cuts ) {
    // system to solve for the new cell
    using EMat = Eigen::Matrix<Sval,nb_dims+1,nb_dims+1>;
    using EVec = Eigen::Matrix<Sval,nb_dims+1,1>;
    EMat M;
    EVec V;
    for( PI num_vertex = 0; num_vertex < vertex_coords.size() / nb_dims; ++num_vertex ) {
        const auto *coords = vertex_coords.data() + num_vertex * ( nb_dims + 0 );
        const auto *cuts = vertex_cuts.data() + num_vertex * ( nb_dims + 1 );
        if ( nb_inf_cuts( cuts ) )
            continue;

        for( PI r = 0; r <= nb_dims; ++r ) {
            PI cut_id = cuts[ r ];

            // cut coming from an affine function
            if ( auto ci = aff_cut( cut_id ) ) {
                Svec dir = m_dirs[ *ci ];
                for( PI c = 0; c < nb_dims; ++c )
                    M( r, c ) = dir[ c ];
                M( r, nb_dims ) = -1;
                V( r ) = m_offs[ *ci ];
                continue;
            }

            // cut coming from a boundary
            if ( auto ci = bnd_cut( cut_id ) ) {
                auto bnd = b_dirs[ *ci ];
                for( PI c = 0; c < nb_dims; ++c )
                    M( r, c ) = bnd[ c ];
                M( r, nb_dims ) = 0;

                Sval v = 0;
                for( PI c = 0; c < nb_dims; ++c )
                    v += bnd[ c ] * coords[ c ];
                V( r ) = v;
                continue;
            }

            // ?
            ERROR( "should not happen" );
        }

        // solve
        Eigen::FullPivLU<EMat> lu( M );
        EVec X = lu.solve( V );

        // data for the new cell
        new_m_dirs << Svec( Vfs::FromIterator(), X.data() );
        new_m_offs << X[ nb_dims ];
    }
}

DTP void UTP::make_new_bnds( Vfs::VecImpl<Svec> &new_b_dirs, Vfs::VecImpl<Sval> &new_b_offs, const auto &vertex_coords, const auto &vertex_cuts ) {
    using EMat = Eigen::Matrix<Sval,nb_dims+1,nb_dims+1>;
    using EVec = Eigen::Matrix<Sval,nb_dims+1,1>;
    using std::sqrt;

    EMat M;
    // EVec V;
    for( PI num_vertex = 0; num_vertex < vertex_coords.size() / nb_dims; ++num_vertex ) {
        const auto *coords = vertex_coords.data() + num_vertex * ( nb_dims + 0 );
        const auto *cuts = vertex_cuts.data() + num_vertex * ( nb_dims + 1 );
        if ( nb_inf_cuts( cuts ) != 1 )
            continue;

        for( PI r = 0; r <= nb_dims; ++r ) {
            // cell cut
            if ( auto ci = aff_cut( cuts[ r ] ) ) {
                Svec dir = m_dirs[ *ci ];
                for( PI c = 0; c < nb_dims; ++c )
                    M( r, c ) = dir[ c ];
                M( r, nb_dims ) = -1;
                // V( r ) = 0;
                continue;
            }

            // boundary cut
            if ( auto ci = bnd_cut( cuts[ r ] ) ) {
                Svec bnd = b_dirs[ *ci ];
                for( PI c = 0; c < nb_dims; ++c )
                    M( r, c ) = bnd[ c ];
                M( r, nb_dims ) = 0;
                // V( r ) = 0;
                continue;
            }

            // infnite cut => we say for now that the sum of the coefficients must be == 1
            //  the direction will be corrected in a second phase
            if ( auto ci = inf_cut( cuts[ r ] ) ) {
                //     for( PI c = 0; c < nb_dims; ++c )
                //         M( r, c ) = 1;
                //     M( r, nb_dims ) = 0;
                //     V( r ) = 1;
                for( PI c = 0; c <= nb_dims; ++c )
                    M( r, c ) = 0;
                continue;
            }

            // ?
            ERROR( "should not happen" );
        }


        // solve
        Eigen::FullPivLU<EMat> lu( M );
        ASSERT( lu.dimensionOfKernel() >= 1 );
        //std::cout << lu.kernel() << std::endl;
        // TODO;
        EVec X = lu.kernel().col( 0 );

        // coeff for normalization
        Sval cnorm = 0;
        for( PI r = 0; r < nb_dims; ++r )
            cnorm += X[ r ] * X[ r ];
        cnorm = 1 / sqrt( cnorm );

        // check orientation
        auto bad_new_bnd_orientation = [&]() {
            using namespace std;
            Sval best_sp = 0;

            // test with boundaries
            for( PI i = 0; i < b_dirs.size(); ++i ) {
                if ( ! used_bs[ i ] )
                    continue;
                Sval prop_sp = 0;
                Svec b_dir = b_dirs[ i ];
                for( PI d = 0; d < nb_dims; ++d )
                    prop_sp += b_dir[ d ] * X[ d ];
                if ( abs( best_sp ) < abs( prop_sp ) )
                    best_sp = prop_sp;
            }

            // test with other affine functions
            for( PI i = 0; i < m_dirs.size(); ++i ) {
                if ( ! used_ms[ i ] )
                    continue;
                Sval prop_sp = - X[ nb_dims ];
                Svec m_dir = m_dirs[ i ];
                for( PI d = 0; d < nb_dims; ++d )
                    prop_sp += m_dir[ d ] * X[ d ];
                if ( abs( best_sp ) < abs( prop_sp ) )
                    best_sp = prop_sp;
            }

            return best_sp > 0;
        };
        if ( bad_new_bnd_orientation() )
            cnorm = - cnorm;
        X *= cnorm;

        // data for the new cell
        new_b_dirs << Svec( Vfs::FromIterator(), X.data() );
        new_b_offs << X[ nb_dims ];
    }
}

DTP Vfs::Opt<Vfs::PI> UTP::aff_cut( Vfs::PI cut ) const {
    return cut < m_dirs.size() ? cut : Vfs::Opt<Vfs::PI>{};
}

DTP Vfs::Opt<Vfs::PI> UTP::bnd_cut( Vfs::PI cut ) const {
    return cut >= m_dirs.size() && cut - m_dirs.size() < b_dirs.size() ? cut - m_dirs.size() : Vfs::Opt<Vfs::PI>{};
}

DTP Vfs::Opt<Vfs::PI> UTP::inf_cut( Vfs::PI cut ) const {
    return cut >= m_dirs.size() + b_dirs.size() ? cut - ( m_dirs.size() + b_dirs.size() ) : Vfs::Opt<Vfs::PI>{};
}

#undef DTP
#undef UTP
