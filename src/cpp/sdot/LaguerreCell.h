#pragma once

#include "MapOfUniquePISortedArray.h"
#include "LaguerreCellVertex.h"
#include "VtkOutput.h"
#include <functional>

/**
*/
template<class Scalar,class Point,class Weight,int nb_dims>
class LaguerreCell {
public:
    using                PI                      = Vfs::PI;
    using                FaceToInt               = MapOfUniquePISortedArray<PI,nb_dims-2,int>;

    static constexpr PI  beg_index_outside_cut   = std::numeric_limits<Vfs::PI>::max() / 2; ///< for cut.n_index
    using                Vertex                  = LaguerreCellVertex<Scalar,nb_dims,PI>;
    struct               Edge                    { std::array<PI,std::max(nb_dims-1,0)> num_cuts; std::array<PI,2> vertices; Vfs::DisplayItem *display( Vfs::Displayer &ds ) const; };
    struct               Cut                     { PI n_index; Point dir; Scalar sp; Vfs::DisplayItem *display( Vfs::Displayer &ds ) const; };

    void                 make_the_initial_simplex( PI beg_inf_index, const Point &center, Scalar radius );
    void                 init                    ( PI orig_index, PI inf_index, const Point &center, Scalar radius );
    void                 cut                     ( PI n_index, const Point &dir, Scalar off );

    void                 for_each_vertex         ( const std::function<void( const Vertex &v )> &f ) const;
    void                 for_each_face           ( const std::function<void( std::array<PI,nb_dims-2> num_cuts, std::span<const Vertex *> vertices )> &f ) const;
    void                 for_each_edge           ( const std::function<void( std::array<PI,nb_dims-1> num_cuts, const Vertex *v0, const Vertex *v1 )> &f ) const;
    void                 display_vtk             ( VtkOutput &vo, const auto &outside_cut ) const;
    Scalar               measure                 () const { return orig_index; }
    Vfs::DisplayItem*    display                 ( Vfs::Displayer &ds ) const;

    bool                 vertex_has_cut          ( const Vertex &vertex, const auto &outside_cut ) const;
    Point                compute_pos             ( std::array<PI,nb_dims> num_cuts ) const;

    static void          get_compilation_flags   ( Vfs::CompilationFlags &cn ) { cn.add_inc_file( "sdot/LaguerreCell.h" ); }
    static void          for_each_template_arg   ( auto &&f ) { f( Vfs::CtType<Scalar>() ); f( Vfs::CtType<Point>() ); f( Vfs::CtType<Weight>() ); f( Vfs::CtInt<nb_dims>() ); }
    static auto          template_type_name      () { return "LaguerreCell"; }

    TTI static auto      array_without_index     ( const std::array<T,i> &values, PI index );
    TTI static auto      array_with_value        ( const std::array<T,i> &a, T value );
    TT static void       apply_corr              ( std::vector<T> &vec, std::vector<int> &keep );
    static bool          is_ext                  ( const Point &pos, const Point &dir, Scalar off );

    const Weight*        orig_weight;            ///<
    const Point*         orig_point;             ///<
    PI                   orig_index;             ///<

    std::vector<Vertex>  vertices;               ///<
    std::vector<Edge>    edges;                  ///<
    std::vector<Cut>     cuts;                   ///<

    FaceToInt            waiting_vertices;       ///<
    std::vector<int>     vertex_corr;            ///<
    std::vector<int>     edge_corr;              ///<
    std::vector<Scalar>  sps;                    ///<

    mutable PI           curr_op_id = 0;         ///<
};

#include "LaguerreCell.tcc" // IWYU pragma: export
