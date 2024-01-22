import numpy as np
import unittest
import sdot

class TestConvexFunction(unittest.TestCase):
    def test_1d( self ):
        # unbounbed
        self.check_transform(
            m_dirs = [ [ -2.0 ], [ +0.0 ], [ +3.0 ] ],
            m_offs = [   +5.0  ,   +0.0  ,   +4.0   ],
        )

        #
        self.check_transform(
            m_dirs = [ [ -3.0 ], [ +4.0 ] ],
            m_offs = [   +0.0  ,   +0.0   ],
        )

    def test_1d_add( self ):
        a = sdot.ConvexFunction( 
            m_dirs = [ [ 1 ] ],
            m_offs = [ 0 ],
            b_dirs = [ [ -1 ], [ 1 ] ],
            b_offs = [ 0, 2 ],
        )
        b = sdot.ConvexFunction( 
            m_dirs = [ [ 0 ] ],
            m_offs = [ 1 ],
        )
        print( "->", a + b )


    def test_2d( self ):
        self.check_transform( 
            m_dirs = [ [ 0, 0 ], [ -1, 0 ], [ 0, -1 ], [ +1, 0 ], [ 0, +1 ] ], 
            m_offs = [ 0, 1, 1, 1, 1 ]
        )

        # unbounbed 1d test case in a 2D space
        # self.check_transform(
        #     m_dirs = [ [ -2.0, 0.0 ], [ +0.0, 0.0 ], [ +3.0, 0.0 ] ],
        #     m_offs = [   +5.0  ,        +0.0  ,        +4.0   ],
        # )

        # self.check_transform(
        #     m_dirs = [ [ 0.0, -2.0 ], [ 0.0, +0.0 ], [ 0.0, +3.0 ] ],
        #     m_offs = [        +5.0  ,        +0.0  ,        +4.0   ],
        # )

        # eq-bounbed 1d test case in a 2D space
        # self.check_transform(
        #     m_dirs = [ [ -1.0, +0.0 ], [ +1.0, +0.0 ] ],
        #     m_offs = [   +5.0        ,   +2.0         ],
        #     b_dirs = [ [ +0.0, -1.0 ], [ +0.0, +1.0 ] ],
        #     b_offs = [         -1.0  ,         +1.0   ],
        # )

        # 0d test case in a 2D space
        # self.check_transform(
        #     m_dirs = [ [ 0.0, +0.0 ] ],
        #     m_offs = [   +0.0        ],
        # )

        # 1d test case in a 2D space
        # max(0) for -1.0*y_0 <= 3.0, 1.0*y_0 <= 4.0
        # self.check_transform(
        #     m_dirs = [ [ -3.0, 0.0 ], [ +4.0, 0.0 ] ],
        #     m_offs = [   +0.0  ,   +0.0   ],
        # )

        #
        # self.check_transform(
        #     m_dirs = [ [ +0.0, +0.0 ] ],
        #     m_offs = [        +0.0    ],
        #     b_dirs = [ [ +1.0, +0.0 ], [ -1.0, +0.0 ] ],
        #     b_offs = [         +0.0  ,         +0.0   ],
        # )

        # self.check_transform(
        #     m_dirs = [ [ -2.0, 0.0 ], [ +0.0, 0.0 ], [ +3.0, 0.0 ] ],
        #     m_offs = [        +5.0  ,        +0.0  ,        +4.0   ],
        # )

        #
        # self.check_transform(
        #     m_dirs = [ [ -3.0 ], [ +4.0 ] ],
        #     m_offs = [   +0.0  ,   +0.0   ],
        # )

    def test_2d_from_vals_and_ders( self, n = 5 ):
        pass
        # angles = np.linspace( 0, 2 * np.pi, 5, endpoint = False )
        # b_dirs = [ [ np.cos( a ), np.sin( a ) ] for a in angles ]
        # b_offs = [ 2.0 for a in angles ]

        # ca = sdot.ConvexFunction.make_approx_from_values_and_derivatives(
        #     sample_coords = np.random.rand( 2, n ) * 2 - 1,
        #     f_der = lambda p: [ 2 * p[ 0 ], 0.2 * 2 * p[ 1 ] ],
        #     f_val = lambda p: p[ 0 ]**2 + 0.2 * p[ 1 ]**2,
        #     # b_dirs = b_dirs,
        #     # b_offs = b_offs,
        # )
        # ca.write_vtk( "orig.vtk" )

        # lt = ca.legendre_transform()
        # lt.write_vtk( "lege.vtk" )
        # print( lt )

        # tl = lt.legendre_transform()
        # tl.write_vtk( "egel.vtk" )



    # def test_2d_unbounded( self ):
    #     #dir = np.ones( [ 2 ] ) / np.sqrt( 2 )
    #     dir = np.array( [ 0, 1 ] )
    #     off = np.array( [ 2, 5 ] )

    #     m_dirs = [ off - 3 * dir, off + 4 * dir ]
    #     m_offs = [        0.0   ,        0.0    ]

    #     ca = sdot.ConvexFunction( m_dirs, m_offs )
    #     # ca.write_vtk( "orig.vtk" )
    #     print( "ca:", ca )

    #     lt = ca.legendre_transform()
    #     print( "lt:", lt )
    #     # self.assertEqual( lt.summary(), "max(-2.5*y, 1.33*y) for y>=-2.0 and y<=3.0" )

    #     tl = lt.legendre_transform()
    #     print( "tl:", tl )
    #     # self.assertEqual( tl.summary(), "max(0.0, -2.0*y-5.0, 3.0*y-4.0)" )

    # def test_1d_bounded( self ):
    #     bnds = [ [ -1, 0 ], [ 1, 4 ] ]
    #     dirs = [ [ 0 ], [ 1 ] ]
    #     offs = [   0,     1   ]

    #     ca = sdot.ConvexFunction( dirs, offs, bnds )
    #     self.assertEqual( ca.summary(), "max(0.0, 1.0*y-1.0) for y>=-0.0 and y<=4.0" )

    #     lt = ca.legendre_transform()
    #     self.assertEqual( lt.summary(), "max(1.0*y, 0.0, 4.0*y-3.0)" )

    #     tl = lt.legendre_transform()
    #     self.assertEqual( tl.summary(), "max(0.0, 1.0*y-1.0) for y>=-0.0 and y<=4.0" )

    # def test_2d_unbounded( self ):
    #     bnds = []
    #     dirs = [ [ 0, 0 ], [ 3, -1 ], [ -2, 1 ] ]
    #     offs = [      0  ,       4  ,       5   ]

    #     ca = sdot.ConvexFunction( dirs, offs, bnds )
    #     self.assertEqual( ca.summary(), "max(0, 3.0*y_0 - 1.0*y_1 - 4.0, -2.0*y_0 + 1.0*y_1 - 5.0)" )
        
    #     lt = ca.legendre_transform()
    #     self.assertEqual( lt.summary(), "max(9.0*y_0 + 23.0*y_1) for -1.0*y_0 - 3.0*y_1 <= 0.0, -1.0*y_0 - 2.0*y_1 <= 0.0, 2.0*y_0 + 5.0*y_1 <= 1.0" )
    #     lt.write_vtk( "lega.vtk" )

    #     tl = lt.legendre_transform()
    #     self.assertEqual( tl.summary(), "max(0, 3.0*y_0 - 1.0*y_1 - 4.0, -2.0*y_0 + 1.0*y_1 - 5.0)" )

    # def test_2d_inf( self ):
    #     affs = np.zeros( [ 0, 3 ] )
    #     bnds = np.zeros( [ 0, 3 ] )

    #     ca = sdot.ConvexFunction( affs, bnds )
    #     print( ca )
        
    #     lt = ca.legendre_transform()
    #     print( lt )

    #     tl = lt.legendre_transform()
    #     print( tl )
    #     # self.assertEqual( tl.summary(), "max(0, 3.0*y_0 - 1.0*y_1 - 4.0, -2.0*y_0 + 1.0*y_1 - 5.0)" )

    # def test_2d_mix( self ):
    #     # bnds = [ [ -1, 0, 2 ], [ 0, -1, 2 ] ]
    #     bnds = []
    #     dirs = [ [ 0, 0 ], [ 1, 0 ] ] # , [ 0, 1 ], ]
    #     offs = [      0  ,     -1   ] # ,     -1  , ]

    #     ca = sdot.ConvexFunction( dirs, offs, bnds )
    #     ca.write_vtk( "orig.vtk" )
    #     print( ca )
        
    #     lt = ca.legendre_transform()
    #     lt.write_vtk( "lega.vtk" )
    #     print( lt )
    #     # self.assertEqual( lt.summary(), "max(9.0*y_0 + 23.0*y_1) for -1.0*y_0 - 3.0*y_1 <= 0.0, -1.0*y_0 - 2.0*y_1 <= 0.0, 2.0*y_0 + 5.0*y_1 <= 1.0" )

    #     # tl = lt.legendre_transform()
    #     # print( tl )
    #     # self.assertEqual( tl.summary(), "max(0, 3.0*y_0 - 1.0*y_1 - 4.0, -2.0*y_0 + 1.0*y_1 - 5.0)" )

    def check_transform( self, m_dirs, m_offs, b_dirs = None, b_offs = None, e_init = None, e_lege = None ):
        init = sdot.ConvexFunction( m_dirs, m_offs, b_dirs, b_offs )
        if e_init is not None:
            self.assertEqual( init.summary(), e_init )
        else:
            init.write_vtk( "orig.vtk" )
            print( init )
    
        lege = init.legendre_transform()
        if e_lege is not None:
            self.assertEqual( lege.summary(), e_lege )
        else:
            lege.write_vtk( "lege.vtk" )
            print( lege )

        egel = lege.legendre_transform()
        if e_init is not None:
            self.assertEqual( egel.summary(), e_init )
        else:
            egel.write_vtk( "egel.vtk" )
            print( egel )


if __name__ == '__main__':
    unittest.main()

