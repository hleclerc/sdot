#include "../src/sdot/VtkOutput.h"
#include "catch_main.h"

TEST_CASE( "VtkOutput", "" ) {
    VtkOutput vo;
    vo.add_edge( {
        VtkOutput::Pt{ 0.0, 0.0, 0.0 },
        VtkOutput::Pt{ 1.0, 0.0, 0.0 },
    } );
    vo.save( "test.vtk" );
}
