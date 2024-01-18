import numpy as np
import vfs

# pre-loaded functions
make_weighted_point_set_aabb = vfs.function( 'make_weighted_point_set_aabb', [ f'inc_file:sdot/WeightedPointSet_AABB.h' ] )
rt_int = vfs.function( 'Vfs::RtInt', [ 'inc_file:vfs/vfs_system/RtInt.h' ] )

#
class PowerDiagram:
    def __init__( self, points_, weights_, b_dirs = None, b_offs = None ) -> None:
        points = np.asarray( points_ )
        weights = np.asarray( weights_ )
        assert( points.ndim == 2 )
        assert( weights.ndim == 1 )
        assert( points.shape[ 0 ] == weights.shape[ 0 ] )

        # points and weight are actually stored in a "weighted_point_set" structure
        self.wps = make_weighted_point_set_aabb( points, weights, rt_int( points.shape[ 1 ] ) )
        self.nb_dims = points.shape[ 1 ]

        self.b_dirs = np.asarray( b_dirs )
        self.b_offs = np.asarray( b_offs )
        assert( self.b_dirs.ndim == 2 )
        assert( self.b_offs.ndim == 1 )
        assert( self.b_dirs.shape[ 1 ] == self.nb_dims )
        assert( self.b_offs.shape[ 0 ] == self.b_dirs.shape[ 0 ] )

    def write_vtk( self, filename, fit_boundaries = 0.1 ):
        display_vtk_laguerre_cells = vfs.function( 'display_vtk_laguerre_cells', [ f'inc_file:sdot/display_vtk_laguerre_cells.h' ] )
        VtkOutput = vfs.function( 'VtkOutput', [ f'inc_file:sdot/VtkOutput.h' ] )
        save = vfs.method( 'save' )

        vo = VtkOutput()
        display_vtk_laguerre_cells( vo, self.wps, self.b_dirs, self.b_offs, fit_boundaries )
        save( vo, filename )

    def cell_points( self ):
        """ return a tuple with
             * a numpy array with the id of the diracs that made this point 
             * a numpy array with the coordinates 
        """
        func = vfs.function( 'cell_points', [ f'inc_file:sdot/cell_points.h' ] )
        
        coords, ids = func( self.wps, self.b_dirs, self.b_offs )
        return np.reshape( coords, [ -1, self.nb_dims ] ), np.reshape( ids, [ -1, self.nb_dims + 1 ] )

