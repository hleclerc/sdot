#pragma once

#include <eigen3/Eigen/LU>
#include "LaguerreCell.h"

#define DTP template<class Scalar,class Point,class Weight,int nb_dims>
#define UTP LaguerreCell<Scalar,Point,Weight,nb_dims>

DTP void UTP::init( PI orig_index, PI beg_inf_index, const Point &center, Scalar radius ) {
    make_the_initial_simplex( beg_inf_index, center, radius );
    this->orig_index = orig_index;
}

DTP TTI auto UTP::array_without_index( const std::array<T,i> &values, PI index ) {
    std::array<T,i-1> res;
    for( int d = 0, o = 0; d < i; ++d )
        if ( d != index )
            res[ o++ ] = values[ d ];
    return res;
}

DTP void UTP::make_the_initial_simplex( PI beg_inf_index, const Point &center, Scalar radius ) {
    vertices.clear();
    edges.clear();
    cuts.clear();

    // cuts
    for( int d = 0; d < nb_dims; ++d ) {
        Point dir( Vfs::FromItemValue(), 0 );
        dir[ d ] = -1;
        cuts.emplace_back( beg_inf_index++, dir, radius - center[ d ] );
    }

    Point dir( Vfs::FromItemValue(), 1 );
    Point vrd( Vfs::FromItemValue(), radius );
    cuts.emplace_back( beg_inf_index++, dir, sp( center + vrd, dir ) );

    // vertices
    for( int nc_0 = 0; nc_0 < nb_dims + 1; ++nc_0 ) {
        std::array<PI,nb_dims> num_cuts;
        for( int i = 0; i < nc_0; ++i )
            num_cuts[ i ] = i;
        for( int i = nc_0 + 1; i < nb_dims + 1; ++i )
            num_cuts[ i - 1 ] = i;

        vertices.emplace_back( num_cuts, compute_pos( num_cuts ) );
    }

    // edges
    for( int nc_0 = 0; nc_0 < nb_dims; ++nc_0 ) {
        std::array<PI,nb_dims-1> num_cuts;
        for( int i = 0; i < nc_0; ++i )
            num_cuts[ i ] = i;

        for( int nc_1 = nc_0 + 1; nc_1 < nb_dims + 1; ++nc_1 ) {
            for( int i = nc_0 + 1; i < nc_1; ++i )
                num_cuts[ i - 1 ] = i;
            for( int i = nc_1 + 1; i < nb_dims + 1; ++i )
                num_cuts[ i - 2 ] = i;

            edges.emplace_back( num_cuts, std::array<PI,2>{ PI( nc_0 ), PI( nc_1 ) } );
        }
    }
}

DTP void UTP::for_each_vertex( const std::function<void( const Vertex &v )> &f ) const {
    for( const Vertex &v : vertices )
        f( v );
}

DTP void UTP::for_each_face( const std::function<void( std::array<PI,nb_dims-2> num_cuts, std::span<const Vertex *> vertices )> &f ) const {
    // face => edge list
    std::map<std::array<PI,nb_dims-2>,std::vector<std::array<const Vertex *,2>>> map;
    for_each_edge( [&]( std::array<PI,nb_dims-1> bc, const Vertex *v0, const Vertex *v1 ) {
        for( PI nf = 0; nf < nb_dims - 1; ++nf ) {
            std::array<PI,nb_dims-2> bd = array_without_index( bc, nf );
            map[ bd ].push_back( { v0, v1 } );
        }
    } );
    if ( map.empty() )
        return;

    // for each face
    std::vector<const Vertex *> vertices;
    for( const auto &face : map ) {
        vertices.clear();
        ++curr_op_id;

        // helper to find an unused vertex
        auto find_unused_vertex = [&]() -> const Vertex * {
            for( const auto &a : face.second )
                if ( a[ 0 ]->op_id != curr_op_id )
                    return a[ 0 ];
            return nullptr;
        };

        // make edge sweep. TODO: avoid this O(n^2) behavior
        while ( const Vertex *p = find_unused_vertex() ) {
            vertices.push_back( p );
            p->op_id = curr_op_id;

            // helper to find the next vertex in the face
            auto find_next_vertex = [&]() -> const Vertex * {
                for( const auto &a : face.second ) {
                    if ( a[ 0 ] == p && a[ 1 ]->op_id != curr_op_id )
                        return a[ 1 ];
                    if ( a[ 1 ] == p && a[ 0 ]->op_id != curr_op_id )
                        return a[ 0 ];
                }
                return nullptr;
            };

            //
            while ( const Vertex *n = find_next_vertex() ) {
                vertices.push_back( n );
                n->op_id = curr_op_id;
                p = n;
            }

            //
            f( face.first, vertices );
        }
    }
}

DTP void UTP::for_each_edge( const std::function<void( std::array<PI,nb_dims-1> num_cuts, const Vertex *v0, const Vertex *v1 )> &f ) const {
    for( const Edge &e : edges )
        f( e.num_cuts, vertices.data() + e.vertices[ 0 ], vertices.data() + e.vertices[ 1 ] );
}

DTP Point UTP::compute_pos( std::array<PI,nb_dims> num_cuts ) const {
    using TM = Eigen::Matrix<Scalar,nb_dims,nb_dims>;
    using TV = Eigen::Matrix<Scalar,nb_dims,1>;

    TM m;
    TV v;
    for( PI i = 0; i < nb_dims; ++i ) {
        for( PI j = 0; j < nb_dims; ++j )
            m( i, j ) = cuts[ num_cuts[ i ] ].dir[ j ];
        v( i ) = cuts[ num_cuts[ i ] ].sp;
    }

    Eigen::PartialPivLU<TM> lu( m );
    TV x = lu.solve( v );

    Point res;
    for( PI i = 0; i < nb_dims; ++i )
        res[ i ] = x[ i ];

    return res;
}

DTP bool UTP::vertex_has_cut( const Vertex &vertex, const auto &cut_func ) const {
    for( auto num_cut : vertex.num_cuts )
        if ( cut_func( cuts[ num_cut ].n_index ) )
            return true;
    return false;
}

DTP void UTP::display_vtk( VtkOutput &vo, const auto &outside_cut ) const {
    auto to_vtk = [&]( const auto &pos ) {
        VtkOutput::Pt res;
        for( PI i = 0; i < min( PI( pos.size() ), res.size() ); ++i )
            res[ i ] = pos[ i ];
        for( PI i = PI( pos.size() ); i < res.size(); ++i )
            res[ i ] = 0;
        return res;
    };

    auto add_item = [&]( int vtk_id, std::span<const Vertex *> vertices ) {
        Vfs::VecImpl<VtkOutput::Pt> points;
        VtkOutput::VTF convex_function;
        VtkOutput::VTF is_outside;
        for( const Vertex *vertex : vertices ) {
            convex_function << Vfs::sp( vertex->pos, *orig_point ) - ( Vfs::norm_2_p2( *orig_point ) - *orig_weight ) / 2;
            is_outside << vertex_has_cut( *vertex, outside_cut );
            points << to_vtk( vertex->pos );
        }
        vo.add_polygon( points, { { "convex_function", convex_function }, { "is_outside", is_outside } } );
    };

    // edges
    if constexpr ( nb_dims >= 1 ) {
        for_each_edge( [&]( std::array<PI,nb_dims-1> num_cuts, const Vertex *v0, const Vertex *v1 ) {
            const Vertex *vs[] = { v0, v1 };
            add_item( VtkOutput::VtkLine, vs );
        } );
    }

    // faces
    if constexpr ( nb_dims >= 2 ) {
        for_each_face( [&]( std::array<PI,nb_dims-2> bc, std::span<const Vertex *> vertices ) {
            add_item( VtkOutput::VtkPolygon, vertices );
        } );
    }
}

DTP Vfs::DisplayItem *UTP::display( Vfs::Displayer &ds ) const {
    return DS_OBJECT( LaguerreCell, vertices, edges, cuts );
}

DTP bool UTP::is_ext( const Point &pos, const Point &dir, Scalar off ) {
    return sp( pos, dir ) > off;
}

DTP TT void UTP::apply_corr( std::vector<T> &vec, std::vector<int> &keep ) {
    int last_keep = keep.size();
    for( int i = 0; i < last_keep; ++i ) {
        if ( keep[ i ] ) {
            keep[ i ] = i;
            continue;
        }

        while( --last_keep > i && ! keep[ last_keep ] )
            keep[ last_keep ] = -1;

        vec[ i ] = std::move( vec[ last_keep ] );
        keep[ last_keep ] = i;
        keep[ i ] = -1;
    }

    vec.resize( last_keep );
}

DTP TTI auto UTP::array_with_value( const std::array<T,i> &a, T value ) {
    std::array<T,i+1> res;
    for( PI n = 0; n < i; ++n )
        res[ n ] = a[ n ];
    res[ i ] = value;
    return res;
}

DTP void UTP::cut( PI n_index, const Point &dir, Scalar off ) {
    // ext/int for each vertex
    vertex_corr.resize( vertices.size() );
    bool has_ext = false;
    for( PI num_vertex = 0; num_vertex < vertices.size(); ++num_vertex ) {
        bool ext = is_ext( vertices[ num_vertex ].pos, dir, off );
        vertex_corr[ num_vertex ] = ! ext;
        has_ext |= ext;
    }

    // all int ?
    if ( ! has_ext )
        return;

    // check dir and off are new. TODO: something more robust
    for( const Cut &cut : cuts ) {
        using namespace std;
        if ( Vfs::norm_2_p2( cut.dir - dir ) < 1e-10 && abs( cut.sp - off ) < 1e-10 )
            return;
    }


    // move vertex to the new positions
    apply_corr( vertices, vertex_corr );

    // add the new cut
    PI new_cut = cuts.size();
    cuts.emplace_back( n_index, dir, off );

    //
    if constexpr ( nb_dims >= 2 )
        waiting_vertices.init( cuts.size(), -1 );

    // for each edge
    edge_corr.resize( edges.size() );
    for( PI num_edge = 0, nb_edges = edges.size(); num_edge < nb_edges; ++num_edge ) {
        Edge *edge = &edges[ num_edge ];

        // helper to create the new edges (on faces that have a cut)
        auto add_to_waiting_vertices = [&]( auto face, PI vertex ) {
            int &wv = waiting_vertices[ face ];
            if ( wv >= 0 ) {
                edges.emplace_back( array_with_value( face, new_cut ), std::array<PI,2>{ PI( wv ), vertex } );
                edge = &edges[ num_edge ];
                wv = -1;
            } else
                wv = vertex;
        };

        // all ext => remove it
        bool e0 = vertex_corr[ edge->vertices[ 0 ] ] < 0;
        bool e1 = vertex_corr[ edge->vertices[ 1 ] ] < 0;
        if ( e0 && e1 ) {
            edge_corr[ num_edge ] = 0;
            continue;
        }

        // => we're going to keep this edge (potentially with a modification)
        edge_corr[ num_edge ] = 1;

        // only v0 is ext
        if ( e0 ) {
            edge->vertices[ 1 ] = vertex_corr[ edge->vertices[ 1 ] ];
            edge->vertices[ 0 ] = vertices.size();

            auto num_cuts = array_with_value( edge->num_cuts, new_cut );
            vertices.emplace_back( num_cuts, compute_pos( num_cuts ) );

            // add a waiting vertex for each face
            if constexpr ( nb_dims >= 2 )
                for( int d = 0; d < nb_dims - 1; ++d )
                    add_to_waiting_vertices( array_without_index( edge->num_cuts, d ), edge->vertices[ 0 ] );

            continue;
        }

        // only v1 is ext
        if ( e1 ) {
            edge->vertices[ 0 ] = vertex_corr[ edge->vertices[ 0 ] ];
            edge->vertices[ 1 ] = vertices.size();

            auto num_cuts = array_with_value( edge->num_cuts, new_cut );
            vertices.emplace_back( num_cuts, compute_pos( num_cuts ) );

            // add a waiting vertex for each face
            if constexpr ( nb_dims >= 2 )
                for( int d = 0; d < nb_dims - 1; ++d )
                    add_to_waiting_vertices( array_without_index( edge->num_cuts, d ), edge->vertices[ 1 ] );

            continue;
        }

        // => all int
        for( PI i = 0; i < 2; ++i )
            edge->vertices[ i ] = vertex_corr[ edge->vertices[ i ] ];
    }

    // move edges to the new positions
    while ( edge_corr.size() < edges.size() )
        edge_corr.push_back( 1 );
    apply_corr( edges, edge_corr );
}

DTP Vfs::DisplayItem *UTP::Edge::display( Vfs::Displayer &ds ) const {
    return DS_OBJECT( Edge, num_cuts, vertices );
}

DTP Vfs::DisplayItem *UTP::Cut::display( Vfs::Displayer &ds ) const {
    return DS_OBJECT( Cut, n_index, dir, sp );
}

#undef DTP
#undef UTP
