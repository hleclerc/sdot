#pragma once

#include <vfs/containers/EmptyArrayImpl.h>
#include <vfs/containers/VecImpl.h>
#include <vfs/support/Opt.h>


//auto legendre_transform_on_filled_space( const auto &vertex_coords, const auto &vertex_cuts, const auto &m_dirs, const auto &m_offs, const auto &b_dirs, const auto &b_offs, Vfs::CtType<Sval>, Vfs::CtInt<nb_dims> ) {
template<class Vco,class Vcu,class Mdi,class Mof,class Bdi,class Bof,class CtSval,class CiNdim>
struct LegendreTransformOnFilledSpace {
    static constexpr int nb_dims       = CiNdim::value;
    using                Sval          = CtSval::value;
    using                Svec          = Vfs::VecImpl<Sval,CiNdim::value>;
    using                Vb            = std::vector<bool>;
    using                PI            = Vfs::PI;

    auto                 transform     ();

    void                 make_new_aff  ( Vfs::VecImpl<Sval> &new_m_dirs, Vfs::VecImpl<Sval> &new_m_offs, const auto *coords, const auto *cuts );
    void                 make_new_bnd  ( Vfs::VecImpl<Sval> &new_b_dirs, Vfs::VecImpl<Sval> &new_b_offs, const auto *coords, const auto *cuts );
    Vfs::Opt<PI>         aff_cut       ( Vfs::PI cut ) const;
    Vfs::Opt<PI>         bnd_cut       ( Vfs::PI cut ) const;
    Vfs::Opt<PI>         inf_cut       ( Vfs::PI cut ) const;

    const Vco&           vertex_coords;
    const Vcu&           vertex_cuts;
    const Mdi&           m_dirs;
    const Mof&           m_offs;
    const Bof&           b_dirs;
    const Bdi&           b_offs;
    CtSval               sval;
    CiNdim               ndi;
    const Vb&            used_ms;
    const Vb&            used_bs;
};

auto legendre_transform( const auto &m_dirs, const auto &m_offs, const auto &b_dirs, const auto &b_offs, auto nb_dims ) -> std::tuple<
    Vfs::VecImpl<decltype( m_dirs[ 0 ][ 0 ] )>,
    Vfs::VecImpl<decltype( m_dirs[ 0 ][ 0 ] )>,
    Vfs::VecImpl<decltype( m_dirs[ 0 ][ 0 ] )>,
    Vfs::VecImpl<decltype( m_dirs[ 0 ][ 0 ] )>
>;

///
template<class Sval,int nb_dims>
std::tuple<Vfs::VecImpl<Sval>,Vfs::VecImpl<Sval>,Vfs::VecImpl<Sval>,Vfs::VecImpl<Sval>>
legendre_transform_on_filled_space( const auto &vertex_coords, const auto &vertex_cuts, const auto &m_dirs, const auto &m_offs, const auto &b_dirs, const auto &b_offs, Vfs::CtType<Sval> sval, Vfs::CtInt<nb_dims> ndims, const std::vector<bool> &used_ms, const std::vector<bool> &used_bs ) {
    LegendreTransformOnFilledSpace ltofs( vertex_coords, vertex_cuts, m_dirs, m_offs, b_dirs, b_offs, sval, ndims, used_ms, used_bs );
    return ltofs.transform();
}

template<class Sval>
std::tuple<Vfs::VecImpl<Sval>,Vfs::VecImpl<Sval>,Vfs::VecImpl<Sval>,Vfs::VecImpl<Sval>>
legendre_transform_with_proj( const auto &vertex_coords, const auto &vertex_cuts, const auto &weights, const auto &m_dirs, const auto &m_offs, const auto &b_dirs, const auto &b_offs, auto nb_dims_ori, auto nb_dims_img, const auto &init_dir, const auto &image, const auto &kernel, Vfs::CtType<Sval>, const std::vector<bool> &used_ms, const std::vector<bool> &used_bs ) {
    // check nb_dims_img (must be == image.size())
    using namespace Vfs;
    if constexpr ( nb_dims_img.value > 1 )
        if ( nb_dims_img.value != image.size() )
            return legendre_transform_with_proj( vertex_coords, vertex_cuts, weights, m_dirs, m_offs, b_dirs, b_offs, nb_dims_ori, CtInt<std::max(nb_dims_img-1,1)>(), init_dir, image, kernel, Vfs::CtType<Sval>(), used_ms, used_bs );

    // if space is already filled
    if ( nb_dims_img == nb_dims_ori )
        return legendre_transform_on_filled_space( vertex_coords, vertex_cuts, m_dirs, m_offs, b_dirs, b_offs, CtType<Sval>(), nb_dims_ori, used_ms, used_bs );

    // else, make a projection
    using Svec = VecImpl<Sval,nb_dims_ori>;
    using Pvec = VecImpl<Sval,nb_dims_img>;
    auto proj = [&]( const Svec &v ) {
        Pvec res;
        for( int r = 0; r < nb_dims_img; ++r )
            res[ r ] = sp( image[ r ], v );
        return res;
    };

    VecImpl<Pvec> pm_dirs( FromReservationSize(), m_dirs.size() );
    VecImpl<Sval> pm_offs( FromReservationSize(), m_offs.size() );
    for( PI i = 0; i < m_dirs.size(); ++i ) {
        pm_dirs << proj( m_dirs[ i ] );
        pm_offs << m_offs[ i ];
    }

    VecImpl<Pvec> pb_dirs( FromReservationSize(), b_dirs.size() );
    VecImpl<Sval> pb_offs( FromReservationSize(), b_offs.size() );
    for( PI i = 0; i < b_dirs.size(); ++i ) {
        pb_dirs << proj( b_dirs[ i ] );
        pb_offs << b_offs[ i ];
    }

    // make cell points again
    auto pwps = make_weighted_point_set_aabb( pm_dirs, weights, nb_dims_img );
    const auto pcoords_and_cuts = cell_points( pwps, pb_dirs, b_offs );
    const auto &pvertex_coords = std::get<0>( pcoords_and_cuts );
    const auto &pvertex_cuts = std::get<1>( pcoords_and_cuts );

    // legendre transform in nb_dims_img
    auto res = legendre_transform_on_filled_space( pvertex_coords, pvertex_cuts, pm_dirs, pm_offs, pb_dirs, pb_offs, CtType<Sval>(), nb_dims_img, used_ms, used_bs );

    // inv projection
    auto inv_proj = [&]( auto v ) {
        Svec res( FromItemValue(), 0 );
        for( PI r = 0; r < nb_dims_img; ++r )
            res = res + v[ r ] * image[ r ];
        return res;
    };

    VecImpl<Sval> nm_dirs( FromReservationSize(), nb_dims_ori * m_dirs.size() / nb_dims_img );
    auto &nm_offs = std::get<1>( res );
    for( PI i = 0; i < std::get<0>( res ).size(); i += nb_dims_img.value )
        nm_dirs.append( inv_proj( std::get<0>( res ).data() + i ) );


    VecImpl<Sval> nb_dirs( FromReservationSize(), nb_dims_ori * b_dirs.size() / nb_dims_img );
    auto &nb_offs = std::get<3>( res );
    for( PI i = 0, j = 0; i < std::get<2>( res ).size(); i += nb_dims_img.value, ++j ) {
        Svec b_dir = inv_proj( std::get<2>( res ).data() + i );
        //nb_offs[ j ] += sp( b_dir, init_dir );
        nb_dirs.append( b_dir );
    }

    // add equalities for kernel vectors
    for( PI k = 0; k < kernel.size(); ++k ) {
        const auto s = sp( kernel[ k ], init_dir );

        nb_dirs.append( kernel[ k ] );
        nb_offs << s;

        nb_dirs.append( -kernel[ k ] );
        nb_offs << -s;

    }

    return std::tuple( std::move( nm_dirs ), std::move( nm_offs ), std::move( nb_dirs ), std::move( nb_offs ) );
}

///
template<class Sval,int nb_dims>
auto legendre_transform_without_dir( const auto &m_dirs, const auto &m_offs, const auto &b_dirs, const auto &b_offs, const Vfs::VecImpl<Sval,nb_dims> &kernel_dir ) {
    using Svec = Vfs::VecImpl<Sval,nb_dims>;
    using namespace Vfs;

    // pb: il faut trouver une base. On peut faire un grahm shmidt
    VecImpl<VecImpl<Sval,nb_dims>,nb_dims-1> base;
    PI dir_to_avoid = argmin( abs( kernel_dir ) );
    Sval n2_kernel_dir = norm_2_p2( kernel_dir );
    for( PI i = 0, j = 0; i < nb_dims; ++i ) {
        if ( i != dir_to_avoid ) {
            Svec b{ FromItemValue(), 0 };
            b[ i ] = 1;

            b = b - sp( b, kernel_dir ) * kernel_dir / n2_kernel_dir;
            for( PI k = 0; k < j; ++k )
                b = b - sp( b, base[ k ] ) * base[ k ] / norm_2_p2( base[ k ] );

            base[ j++ ] = b;
        }
    }

    P( base );
    TODO;

    Vfs::VecImpl<Svec> pm_dirs( FromReservationSize(), m_offs.size() );
    Vfs::VecImpl<Sval> pm_offs( FromReservationSize(), m_offs.size() );
    for( PI i = 0; i < m_offs.size(); ++i ) {
        Svec m_dir = m_dirs[ i ];
        m_dir = m_dir - sp( m_dir, kernel_dir ) * kernel_dir / n2_kernel_dir;
        pm_dirs << m_dir;
    }


    Vfs::VecImpl<Svec> pb_dirs( FromReservationSize(), m_offs.size() );
    Vfs::VecImpl<Sval> pb_offs( FromReservationSize(), m_offs.size() );

    auto plt = legendre_transform( pm_dirs, pm_offs, pb_dirs, pb_offs, CtInt<nb_dims-1>() );

    Vfs::VecImpl<Svec> im_dirs( FromReservationSize(), m_offs.size() );
    Vfs::VecImpl<Sval> im_offs( FromReservationSize(), m_offs.size() );
    Vfs::VecImpl<Svec> ib_dirs( FromReservationSize(), m_offs.size() );
    Vfs::VecImpl<Sval> ib_offs( FromReservationSize(), m_offs.size() );

    return std::tuple( im_dirs, im_offs, ib_dirs, ib_offs );
}

/// return ( new_m_dirs, new_m_offs, new_b_dirs, new_b_offs )
auto legendre_transform( const auto &m_dirs, const auto &m_offs, const auto &b_dirs, const auto &b_offs, auto nb_dims ) -> std::tuple<
       Vfs::VecImpl<decltype( m_dirs[ 0 ][ 0 ] )>,
       Vfs::VecImpl<decltype( m_dirs[ 0 ][ 0 ] )>,
       Vfs::VecImpl<decltype( m_dirs[ 0 ][ 0 ] )>,
       Vfs::VecImpl<decltype( m_dirs[ 0 ][ 0 ] )>
> {
    using Sval = decltype( m_dirs[ 0 ][ 0 ] + m_offs[ 0 ] + b_dirs[ 0 ][ 0 ] + b_offs[ 0 ] );
    using Svec = Vfs::VecImpl<Sval,nb_dims>;
    using namespace Vfs;

    // if there's an equality constraint, we're actually working at most on nb_dims - 1. TODO: chose a O( n * log( n ) ) algorithm
    for( PI i = 0; i < b_dirs.size(); ++i ) {
        Svec bdi = b_dirs[ i ];
        for( PI j = 0; j < i; ++j ) {
            Svec bdj = b_dirs[ j ];
            if ( all( bdi == -bdj ) && b_offs[ i ] == -b_offs[ j ] )
                return legendre_transform_without_dir( m_dirs, m_offs, b_dirs, b_offs, bdi );
        }
    }

    // acceleration structure
    VecImpl<Sval> weights( FromReservationSize(), m_offs.size() );
    for( PI i = 0; i < m_offs.size(); ++i )
        weights << norm_2_p2( m_dirs[ i ] ) - 2 * m_offs[ i ];
    auto wps = make_weighted_point_set_aabb( m_dirs, weights, nb_dims );

    // get items that are really used in the diagram
    const auto coords_and_cuts = cell_points( wps, b_dirs, b_offs );
    const auto &vertex_coords = std::get<0>( coords_and_cuts );
    const auto &vertex_cuts = std::get<1>( coords_and_cuts );

    std::vector<bool> used_ms( m_dirs.size(), false );
    std::vector<bool> used_bs( b_dirs.size(), false );
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

    const PI not_used = -1;
    PI first_used_m = not_used;
    for( PI n = 0; n < used_ms.size(); ++n ) {
        if ( used_ms[ n ] ) {
            first_used_m = n;
            break;
        }
    }

    // make a covariance matrix from the directions
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

    Eigen::EigenSolver<CovMat> es( cov_mat );
    Eigen::Matrix<Sval,nb_dims,1> eva = es.eigenvalues().real();
    const Sval lim = eva.maxCoeff() * std::numeric_limits<Sval>::epsilon() * 100; // TODO: something more accurate :)
    VecImpl<Svec> image, kernel;
    for( PI i = 0; i < nb_dims; ++i )
        ( eva[ i ] > lim ? image : kernel ) << es.eigenvectors().col( i ).real();

    // make the legendre tranform in the right number of dimensions
    return legendre_transform_with_proj( vertex_coords, vertex_cuts, weights, m_dirs, m_offs, b_dirs, b_offs, nb_dims, nb_dims, init_dir, image, kernel, CtType<Sval>(), used_ms, used_bs );
}

#include "legendre_transform.tcc" // IWYU pragma: export
