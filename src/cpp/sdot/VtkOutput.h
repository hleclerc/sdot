#pragma once

#include <vfs/containers/VecImpl.h>
#include <map>

/***/
class VtkOutput {
public:
    enum {          VtkTriangle          = 5   };
    enum {          VtkPyramid           = 14  };
    enum {          VtkPolygon           = 7   };
    enum {          VtkWedge             = 13  };
    enum {          VtkPoint             = 1   };
    enum {          VtkTetra             = 10  };
    enum {          VtkHexa              = 12  };
    enum {          VtkLine              = 4   };
    enum {          VtkPoly              = 7   };
    enum {          VtkQuad              = 9   };

    using           TF                   = double;
    using           PI                   = Vfs::PI;
    using           Pt                   = Vfs::VecImpl<TF,3>;
    using           VTF                  = Vfs::VecImpl<TF>;
    using           FieldMap             = std::map<std::string,VTF>;

    /**/            VtkOutput            ();

    void            save                 ( std::string filename ) const;
    void            save                 ( std::ostream &os ) const;

    // fixed #pts
    void            add_triangle         ( std::array<Pt,3> pts, const std::map<std::string,VTF> &point_data = {}, const std::map<std::string,TF> &cell_data = {} );
    void            add_pyramid          ( std::array<Pt,5> pts, const std::map<std::string,VTF> &point_data = {}, const std::map<std::string,TF> &cell_data = {} );
    void            add_wedge            ( std::array<Pt,6> pts, const std::map<std::string,VTF> &point_data = {}, const std::map<std::string,TF> &cell_data = {} );
    void            add_tetra            ( std::array<Pt,4> pts, const std::map<std::string,VTF> &point_data = {}, const std::map<std::string,TF> &cell_data = {} );
    void            add_point            ( Pt pts, const std::map<std::string,VTF> &point_data = {}, const std::map<std::string,TF> &cell_data = {} );
    void            add_hexa             ( std::array<Pt,8> pts, const std::map<std::string,VTF> &point_data = {}, const std::map<std::string,TF> &cell_data = {} );
    void            add_quad             ( std::array<Pt,4> pts, const std::map<std::string,VTF> &point_data = {}, const std::map<std::string,TF> &cell_data = {} );
    void            add_edge             ( std::array<Pt,2> pts, const std::map<std::string,VTF> &point_data = {}, const std::map<std::string,TF> &cell_data = {} );

    // variable #pts
    void            add_polygon          ( std::span<Pt> pts, const std::map<std::string,VTF> &point_data = {}, const std::map<std::string,TF> &cell_data = {} );
    void            add_line             ( std::span<Pt> pts, const std::map<std::string,VTF> &point_data = {}, const std::map<std::string,TF> &cell_data = {} );

    // generic
    void            add_item             ( const Pt *pts_data, Vfs::PI pts_size, Vfs::PI vtk_type, const std::map<std::string,VTF> &point_data, const std::map<std::string,TF> &cell_data );

    // type info
    static void     get_compilation_flags( Vfs::CompilationFlags &cn ) { cn.add_inc_file( "sdot/VtkOutput.h" ); }
    static auto     type_name            () { return "VtkOutput"; }

    FieldMap        point_fields;        ///<
    FieldMap        cell_fields;         ///<
    std::vector<PI> cell_types;
    std::vector<PI> cell_items;
    std::vector<Pt> points;
};

