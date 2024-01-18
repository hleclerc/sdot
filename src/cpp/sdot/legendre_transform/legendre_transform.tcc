#pragma once

#include "legendre_transform.h"
// #include <gmpxx.h>
// #include <set>

// namespace Vfs {
//     auto *display( Displayer &ds, const mpq_class &q ) { return ds.number( q.get_str() ); }
// } // namespace Vfs

#define DTP template<class Vco,class Vcu,class Mdi,class Mof,class Bdi,class Bof,class CtSval,class CiNdim>
#define UTP LegendreTransformOnFilledSpace<Vco,Vcu,Mdi,Mof,Bdi,Bof,CtSval,CiNdim>

DTP Vfs::Opt<Vfs::PI> UTP::aff_cut( Vfs::PI cut ) const {
    return cut < m_dirs.size() ? cut : Vfs::Opt<Vfs::PI>{};
}

DTP Vfs::Opt<Vfs::PI> UTP::bnd_cut( Vfs::PI cut ) const {
    return cut >= m_dirs.size() && cut - m_dirs.size() < b_dirs.size() ? cut - m_dirs.size() : Vfs::Opt<Vfs::PI>{};
}

DTP Vfs::Opt<Vfs::PI> UTP::inf_cut( Vfs::PI cut ) const {
    return cut >= m_dirs.size() + b_dirs.size() ? cut - ( m_dirs.size() + b_dirs.size() ) : Vfs::Opt<Vfs::PI>{};
}

DTP void UTP::make_new_aff( Vfs::VecImpl<Sval> &new_m_dirs, Vfs::VecImpl<Sval> &new_m_offs, const auto *coords, const auto *cuts ) {
    // system to solve for the new cell
    using EMat = Eigen::Matrix<Sval,nb_dims+1,nb_dims+1>;
    using EVec = Eigen::Matrix<Sval,nb_dims+1,1>;
    EMat M;
    EVec V;
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
    for( PI r = 0; r < nb_dims; ++r )
        new_m_dirs << X[ r ];
    new_m_offs << X[ nb_dims ];
}

DTP void UTP::make_new_bnd( Vfs::VecImpl<Sval> &new_b_dirs, Vfs::VecImpl<Sval> &new_b_offs, const auto *coords, const auto *cuts ) {
    using EMat = Eigen::Matrix<Sval,nb_dims+1,nb_dims+1>;
    using EVec = Eigen::Matrix<Sval,nb_dims+1,1>;
    using std::sqrt;

    EMat M;
    EVec V;
    for( PI r = 0; r <= nb_dims; ++r ) {
        // cell cut
        if ( auto ci = aff_cut( cuts[ r ] ) ) {
            Svec dir = m_dirs[ *ci ];
            for( PI c = 0; c < nb_dims; ++c )
                M( r, c ) = dir[ c ];
            M( r, nb_dims ) = -1;
            V( r ) = 0;
            continue;
        }

        // boundary cut
        if ( auto ci = bnd_cut( cuts[ r ] ) ) {
            Svec bnd = b_dirs[ *ci ];
            for( PI c = 0; c < nb_dims; ++c )
                M( r, c ) = bnd[ c ];
            M( r, nb_dims ) = 0;
            V( r ) = 0;
            continue;
        }

        // infnite cut => we say for now that the sum of the coefficients must be == 1
        //  the direction will be corrected in a second phase
        if ( auto ci = inf_cut( cuts[ r ] ) ) {
            for( PI c = 0; c < nb_dims; ++c )
                M( r, c ) = 1;
            M( r, nb_dims ) = 0;
            V( r ) = 1;
            continue;
        }

        // ?
        ERROR( "should not happen" );
    }

    // solve
    Eigen::FullPivLU<EMat> lu( M );
    EVec X = lu.solve( V );

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

    // data for the new cell
    for( PI r = 0; r < nb_dims; ++r )
        new_b_dirs << cnorm * X[ r ];
    new_b_offs << cnorm * X[ nb_dims ];
}

DTP auto UTP::transform() {
    using namespace Vfs;

    //
    VecImpl<Sval> new_m_dirs;
    VecImpl<Sval> new_m_offs;
    VecImpl<Sval> new_b_dirs;
    VecImpl<Sval> new_b_offs;
    for( PI num_vertex = 0, nb_vertices = vertex_coords.size() / nb_dims; num_vertex < nb_vertices; ++num_vertex ) {
        const auto *coords = vertex_coords.data() + num_vertex * nb_dims;
        const auto *cuts = vertex_cuts.data() + num_vertex * ( nb_dims + 1 );

        // infinite vertices do not contribute to new affine funcs but may create boundaries
        PI nb_infinite_cuts = 0;
        for( PI r = 0; r <= nb_dims; ++r )
            if ( inf_cut( cuts[ r ] ) )
                ++nb_infinite_cuts;
        if ( nb_infinite_cuts ) {
            if ( nb_infinite_cuts == 1 )
                make_new_bnd( new_b_dirs, new_b_offs, coords, cuts );
            continue;
        }

        // else, we have a new affine func to create and add
        make_new_aff( new_m_dirs, new_m_offs, coords, cuts );
    }

    return std::tuple( std::move( new_m_dirs ), std::move( new_m_offs ), std::move( new_b_dirs ), std::move( new_b_offs ) );
}

#undef DTP
#undef UTP
