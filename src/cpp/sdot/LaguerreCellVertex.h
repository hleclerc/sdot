#pragma once

#include <vfs/containers/VecImpl.h>



///
template<class ST,int nb_dims,class PI=Vfs::PI>
struct LaguerreCellVertex {
    using                  VI         = Vfs::VecImpl<PI,nb_dims>;
    using                  Pt         = Vfs::VecImpl<ST,nb_dims>;

    Vfs::DisplayItem*      display    ( Vfs::Displayer &ds ) const{ return DS_OBJECT( Vertex, num_cuts, pos ); }

    VI                     num_cuts;  ///<
    Pt                     pos;       ///<

    mutable PI             op_id = 0; ///<
};

