import numpy as np
import sdot

def test_2_cells():
    bnds = [ [ np.cos( a ), np.sin( a ), 4.0 ] for a in np.linspace( 0, 2 * np.pi, 4, endpoint = False ) ]
    pnts = [ [ 0.0, 0.0 ], [ 1.0, 0.0 ] ]
    wgts = [ 0.2, 0.0 ]

    pd = sdot.PowerDiagram( pnts, wgts, [ bnds ] )
    print( pd.cell_points() )
    # pd.write_vtk( "smurf.vtk" )

def test_rand( n ):
    bnds = [ [ np.cos( a ), np.sin( a ), 2.0 ] for a in np.linspace( 0, 2 * np.pi, 5, endpoint = False ) ]
    pnts = np.random.rand( n, 2 ) - 0.5
    wgts = np.ones( n )

    pd = sdot.PowerDiagram( pnts, wgts, [ bnds ] )
    pd.write_vtk( "smurf.vtk" )

def test_2_cells_with_pysdot():
    import pysdot

    pnts = [ [ 0.0, 0.0 ], [ 1.0, 0.0 ] ]
    wgts = [ 0.2, 0.0 ]

    p = pysdot.PowerDiagram( pnts, wgts )
    p.display_vtk( "yo.vtk" )

# test_2_cells_with_pysdot()
test_2_cells()
