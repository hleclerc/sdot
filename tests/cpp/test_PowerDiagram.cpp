#include "../../src/sdot/PowerDiagram.h"
#include "../catch_main.h"

TEST_CASE( "PowerDiagram", "" ) {
    const int dim = 2;
    const int nb_points = 10;
    Vector v = Vector::randu( dim * nb_points, 0.0, 1.0 );
    auto points = List<Point>::from_size_and_external_item_values( nb_points, { .item_type = Point::type_for( dim ) }, v.data() );
    auto weights = Vector::ones( nb_points );

    WeightedPointSet wps = make_WeightedPointSet_PolytopRec( points, weights );
    PowerDiagram pd( &wps );

    VtkOutput vo;
    pd.for_each_cell_full( [&]( int nb_threads, const auto &f ) {
        f( [&]( const LaguerreCell &cell, int ) {
            cell.display_vtk( vo );
        });
    } );
    vo.save( "test.vtk" );

    // P( pd.volumes() );
}
