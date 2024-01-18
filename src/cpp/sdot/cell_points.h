#pragma once

#include "for_each_laguerre_cell.h"
#include <map>

/// tuple( coord_list, id_list, inf_dirs );
auto cell_points( const auto &weighted_point_set, const auto &b_dirs, const auto &b_offs, Vfs::PI beg_index_wps, Vfs::PI beg_index_bnd, Vfs::PI beg_index_inf ) {
    using Wps = DECAYED_TYPE_OF( weighted_point_set );
    constexpr int nb_dims = Wps::nb_dims;
    using Scalar = Wps::Scalar;
    using namespace Vfs;

    std::map<std::array<PI,nb_dims+1>,PI> num_map;
    PI nb_points = 0;

    VecImpl<Scalar> coord_list;
    VecImpl<PI> id_list;

    std::mutex m;
    for_each_laguerre_cell( weighted_point_set, b_dirs, b_offs, [&]( Vfs::PI nb_threads, auto &&func ) {
        func( [&]( const auto &cell, int ) {
            m.lock();
            std::array<PI,nb_dims+1> ids;
            for( const auto &vertex : cell.vertices ) {
                for( PI d = 0; d < nb_dims; ++d )
                    ids[ d ] = cell.cuts[ vertex.num_cuts[ d ] ].n_index;
                ids[ nb_dims ] = cell.orig_index;
                std::sort( ids.begin(), ids.end() );

                auto iter = num_map.find( ids );
                if ( iter != num_map.end() )
                    continue;

                PI num = nb_points++;
                num_map.insert( iter, { ids, num } );

                for( PI d = 0; d < nb_dims; ++d )
                    coord_list.push_back( vertex.pos[ d ] );
                for( PI d = 0; d < nb_dims + 1; ++d )
                    id_list.push_back( ids[ d ] );
            }
            m.unlock();
        } );
    }, beg_index_wps, beg_index_bnd, beg_index_inf );

    return std::tuple( std::move( coord_list ), std::move( id_list ) );
}

auto cell_points( const auto &weighted_point_set, const auto &b_dirs, const auto &b_offs ) {
    Vfs::PI ws = weighted_point_set.size();
    return cell_points( weighted_point_set, b_dirs, b_offs, 0, ws, ws + b_dirs.size() );
}
