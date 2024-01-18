#pragma once

#include "for_each_laguerre_cell.h"

///
void display_vtk_laguerre_cells( VtkOutput &vo, auto &&weighted_point_set, auto &&b_dirs, auto &&b_offs, double fit_boundaries ) {
    using Wps = DECAYED_TYPE_OF( weighted_point_set );
    using namespace Vfs;
    using ST = Wps::Scalar;
    using Pt = VecImpl<ST,Wps::nb_dims>;

    const PI beg_bnd = weighted_point_set.size();
    const PI beg_inf = beg_bnd + b_dirs.size();

    const auto inf_cut = [&]( auto n_index ) {
        return n_index >= beg_inf;
    };

    // if we have to find a precise hull
    std::mutex m;
    VecImpl<Pt> dirs;
    VecImpl<ST> offs;
    if ( fit_boundaries >= 0 ) {
        // get some directions
        if ( Wps::nb_dims == 1 ) {
            dirs.push_back_br( -1 );
            dirs.push_back_br( +1 );
        } else if ( Wps::nb_dims == 2 ) {
            for( PI i = 0, n = 20; i < n; ++i ) {
                ST a = 2 * M_PI * i / n;
                dirs.push_back_br( cos( a ), sin( a ) );
            }
        } else {
            TODO;
        }

        // get max scalar product for each direction
        offs.resize( dirs.size(), std::numeric_limits<ST>::lowest() );
        for_each_laguerre_cell( weighted_point_set, b_dirs, b_offs, [&]( Vfs::PI nb_threads, auto &&func ) {
            func( [&]( const auto &cell, int ) {
                m.lock();
                cell.for_each_vertex( [&]( const auto &vertex ) {
                    if ( ! cell.vertex_has_cut( vertex, inf_cut ) )
                        for( PI i = 0; i < dirs.size(); ++i )
                            offs[ i ] = std::max( offs[ i ], sp( dirs[ i ], vertex.pos ) );
                } );
                m.unlock();
            } );
        }, 0, beg_bnd, beg_inf );

        if ( offs.size() && offs[ 0 ] == std::numeric_limits<ST>::lowest() ) {
            dirs.clear();
            offs.clear();
        }
    }

    //
    PI cur_inf = beg_inf;
    for_each_laguerre_cell( weighted_point_set, b_dirs, b_offs, [&]( Vfs::PI nb_threads, auto &&func ) {
        func( [&]( auto &cell, int ) {
            for( PI i = 0; i < offs.size(); ++i )
                cell.cut( cur_inf++, dirs[ i ], offs[ i ] + fit_boundaries );

            m.lock();
            cell.display_vtk( vo, inf_cut );
            m.unlock();
        } );
    } );
}
