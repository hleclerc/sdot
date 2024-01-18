#include "../../src/sdot/ConvexFunction.h"
#include "../catch_main.h"

// List<Point> rand_points( int dim, int nb_points ) {
//     List<Point> res;
//     for( PI i = 0; i < nb_points; ++i )
//         res << Point::randu( dim, 0.0, 1.0 );
//     return res;
// }

TEST_CASE( "FacetedConvexFunction", "" ) {
    // auto cf = ConvexFunction::from_func( rand_points( 2, 10 ), [&]( Point p ) {
    //     return 0; // sum( p * p );
    // } );

    //    Tensor pts = Tensor::rand( { 3, 10 } );
    //    Vector val = Vector::map( pts.size( 1 ), [&]( const auto &index ) {
    //        return sum( pts( _, index ) );
    //    } );
    //    auto fc = FacetedConvexFunction::from_values( Tensor::rand( { 3, 10 } ) );

    //    VtkOutput vo;
    //    fc.display( vo );
    //    vo.save( "test.vtk" );

    // P( pd.volumes() );
}
