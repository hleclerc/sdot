#pragma once

#include "BoxForWeightedPointSet_AABB.h"
#include <vfs/containers/VecImpl.h>

///
template<class Scalar_,class Point_,class Weight_,int nb_dims_>
class WeightedPointSet_AABB {
public:
    static constexpr int nb_dims               = nb_dims_;
    using                Scalar                = Scalar_;
    using                Weight                = Weight_;
    using                Point                 = Point_;

    using                Box                   = BoxForWeightedPointSet_AABB<Scalar,Point,Weight,nb_dims>;
    using                PI                    = Vfs::PI;

    static void          for_each_template_arg ( auto &&f );
    static void          get_compilation_flags ( Vfs::CompilationFlags &cn );
    static auto          template_type_name    ();

    void                 set_points_and_weights( const auto &points, const auto &weights );
    Box*                 make_box_rec          ( const auto &points, const auto &weights, std::span<PI> indices, Box *parent );
    auto*                display               ( Vfs::Displayer &ds ) const;
    PI                   size                  () const { return nb_points; }

    // params
    bool                 get_axes_using_eigen_values = true;
    PI                   max_nb_points_per_cell = 30;

    Scalar               init_radius_mul = 10;
    PI                   nb_points = 0;
    Vfs::BumpPointerPool pool;

    Vfs::Vec<Box*>       leaves;
    Box*                 root;
};

auto make_weighted_point_set_aabb( const auto &points, const auto &weights, auto nb_dims ) {
    using namespace std;
    using namespace Vfs;
    using Scalar = decltype( exact_div( points[ 0 ][ 0 ] * weights[ 0 ], 1 + points[ 0 ][ 1 ] ) ); // TODO: use of exact_div
    using Weight = DECAYED_TYPE_OF( weights[ 0 ] );
    using TP = DECAYED_TYPE_OF( points[ 0 ][ 0 ] );
    constexpr int dim = GET_DT_VALUE( nb_dims );
    using Point = Vfs::VecImpl<TP,dim>;

    WeightedPointSet_AABB<Scalar,Point,Weight,dim> res;
    res.set_points_and_weights( points, weights );
    return res;
}

#include "WeightedPointSet_AABB.tcc" // IWYU pragma: export
