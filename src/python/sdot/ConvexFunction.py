from .PowerDiagram import PowerDiagram, rt_int
import numpy as np
import vfs

class ConvexFunction:
    """ 
        Function =
          max( 
            -inf, 
            scalar_product( m_dirs[ 0 ], y ) - m_offs[ 0 ], 
            scalar_product( m_dirs[ 1 ], y ) - m_offs[ 1 ], 
            ...
          ) 
          for y such that 
            scalar_product( b_dirs[ 0 ], y ) - b_offs[ 0 ] <= 0,
            scalar_product( b_dirs[ 1 ], y ) - b_offs[ 1 ] <= 0,
            ...
    """

    def __init__( self, m_dirs, m_offs, b_dirs = None, b_offs = None, dtype = np.double ) -> None:
        # m_...
        self.m_dirs = np.asarray( m_dirs, dtype = dtype )
        self.m_offs = np.asarray( m_offs, dtype = dtype )

        assert( self.m_dirs.ndim == 2 )
        assert( self.m_offs.ndim == 1 )
        assert( self.m_offs.shape[ 0 ] == self.m_dirs.shape[ 0 ] )

        # b_...
        if b_dirs is None or not len( b_dirs ):
            assert( b_offs is None or not len( b_offs ) )
            b_dirs = np.zeros( [ 0, self.m_dirs.shape[ 1 ] ], dtype = dtype )
            b_offs = np.zeros( [ 0 ], dtype = dtype )
        self.b_dirs = np.asarray( b_dirs, dtype = dtype )
        self.b_offs = np.asarray( b_offs, dtype = dtype )

        assert( self.b_dirs.ndim == 2 )
        assert( self.b_offs.ndim == 1 )
        assert( self.b_offs.shape[ 0 ] == self.b_dirs.shape[ 0 ] )
        assert( self.b_dirs.shape[ 1 ] == self.m_dirs.shape[ 1 ] )

        # attributes by deduction
        self.nb_dims = self.m_dirs.shape[ 1 ]

    def make_approx_from_values_and_derivatives( sample_coords, f_val, f_der, b_dirs, b_offs ):
        m_dirs = []
        m_offs = []
        for num_point in range( sample_coords.shape[ 1 ] ):
            point = sample_coords[ :, num_point ]
            val = f_val( point )
            der = f_der( point )

            m_dirs.append( der )
            m_offs.append( np.dot( der, point ) - val )
        return ConvexFunction( m_dirs, m_offs, b_dirs, b_offs )

    def to_power_diagram( self ):
        weights = np.sum( self.m_dirs * self.m_dirs, axis = 1 ) - 2 * self.m_offs
        return PowerDiagram( self.m_dirs, weights, self.b_dirs, self.b_offs )

    def write_vtk( self, filename, fit_boundary = 1.0 ):
        pd = self.to_power_diagram()
        pd.write_vtk( filename, fit_boundary )

    def summary( self ):
        from sympy import Symbol, tensorproduct

        y = [ Symbol( f'y_{ d }' ) for d in range( self.nb_dims ) ]

        repr_max = lambda i: str( sum( tensorproduct( self.m_dirs[ i, : ], y ) ) -  self.m_offs[ i ] )
        repr_bnd = lambda i: str( sum( tensorproduct( self.b_dirs[ i, : ], y ) ) <= self.b_offs[ i ] )

        lst_max = [ repr_max( i ) for i in range( self.m_dirs.shape[ 0 ] ) ]
        lst_bnd = [ repr_bnd( i ) for i in range( self.b_dirs.shape[ 0 ] ) ]
        lst_max.sort()
        lst_bnd.sort()

        res = f"max({ ', '.join( lst_max ) })"
        if len( lst_bnd ):
            res += f" for { ', '.join( lst_bnd ) }"
        return res

    def legendre_transform( self ):
        # call the legendre_transform func
        func = vfs.function( 'legendre_transform', [ f'inc_file:sdot/legendre_transform/legendre_transform.h' ] )
        new_m_dirs, new_m_offs, new_b_dirs, new_b_offs = func( self.m_dirs, self.m_offs, self.b_dirs, self.b_offs, rt_int( self.nb_dims ) )
        new_m_dirs = np.reshape( new_m_dirs, [ -1, self.nb_dims ] )
        new_b_dirs = np.reshape( new_b_dirs, [ -1, self.nb_dims ] )

        # make a new ConvexFunction
        return ConvexFunction( new_m_dirs, new_m_offs, new_b_dirs, new_b_offs )

    def __repr__( self ):
        return self.summary()
