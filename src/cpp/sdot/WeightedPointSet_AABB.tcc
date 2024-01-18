#pragma once

#include "WeightedPointSet_AABB.h"

#define DTP template<class Scalar,class Point,class Weight,int nb_dims>
#define UTP WeightedPointSet_AABB<Scalar,Point,Weight,nb_dims>

DTP void UTP::for_each_template_arg( auto &&f ) {
    f( Vfs::CtType<Scalar>() );
    f( Vfs::CtType<Point>() );
    f( Vfs::CtType<Weight>() );
    f( Vfs::CtInt<nb_dims>() );
}

DTP void UTP::get_compilation_flags( Vfs::CompilationFlags &cn ) {
    cn.add_inc_file( "sdot/WeightedPointSet_AABB.h" );
}

DTP auto UTP::template_type_name() {
    return "WeightedPointSet_AABB";
}

DTP void UTP::set_points_and_weights( const auto &points, const auto &weights ) {
    nb_points = points.size();
    auto indices = Vfs::VecImpl<PI>::range( nb_points );
    root = make_box_rec( points, weights, { indices.data(), indices.size() }, nullptr );
    root->sibling = root;
}

DTP UTP::Box *UTP::make_box_rec( const auto &points, const auto &weights, std::span<PI> indices, Box *parent ) {
    using std::min, std::max;

    Box *box = pool.create<Box>();
    box->parent = parent;
    if ( indices.empty() )
        return box;

    Point pt = points[ indices[ 0 ] ];
    Point mi = pt;
    Point ma = pt;
    for( PI i = 1; i < indices.size(); ++i ) {
        Point pt = points[ indices[ i ] ];
        for( PI d = 0; d < nb_dims; ++d ) {
            mi[ d ] = min( mi[ d ], pt[ d ] );
            ma[ d ] = max( ma[ d ], pt[ d ] );
        }
    }
    box->min_pos = mi;
    box->max_pos = ma;

    // final ?
    if ( indices.size() <= max_nb_points_per_cell ) {
        box->indices.reserve( indices.size() );
        box->weights.reserve( indices.size() );
        box->points.reserve( indices.size() );
        for( PI i = 0; i < indices.size(); ++i ) {
            box->indices.push_back( indices[ i ] );
            box->weights.push_back( weights[ indices[ i ] ] );
            box->points.push_back( points[ indices[ i ] ] );
        }

        leaves.push_back( box );
        return box;
    }

    // else, create a box with several children
    int max_d = 0;
    for( int d = 1; d < nb_dims; ++d )
        if ( ma[ max_d ] - mi[ max_d ] < ma[ d ] - mi[ d ] )
            max_d = d;

    auto mid = indices.begin() + indices.size() / 2;
    std::nth_element( indices.begin(), mid, indices.end(), [&]( PI a, PI b ) { return points[ a ][ max_d ] < points[ b ][ max_d ]; } );

    Box *c0 = make_box_rec( points, weights, { indices.begin(), mid }, box );
    Box *c1 = make_box_rec( points, weights, { mid, indices.end() }, box );
    box->children.push_back( c0 );
    box->children.push_back( c1 );
    c0->sibling = c1;
    c1->sibling = c0;

    return box;
}


DTP auto *UTP::display( Vfs::Displayer &ds ) const {
    return DS_OBJECT( WeightedPointSet_AABB, nb_points, root );
}

#undef DTP
#undef UTP

