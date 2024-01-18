import numpy as np
import scipy.sparse as spr
import gurobipy as grb
import sympy
from sympy.solvers import solve
from sympy import *
import matplotlib.pyplot as plt
import tabulate as tb
from mec.ph import Polyhedral,polyhedral_from_strings

# u = Polyhedral(np.array([[0],[1]]),np.array([0,1]),np.array([[-1],[1]]),np.array([0,4]))
# u.plot1d()
# u = Polyhedral(np.array([[-2],[0],[3]]),np.array([5,0,4]) )
# print( u )
# u.plot1d()

# u = Polyhedral(np.array([[-2,1],[0,0],[3,-1]]),np.array([5,0,4]) )
def ex_mp_1():
    afds = np.array( [ [ 0, 0 ] ] )
    afos = np.array( [ 0 ] )
    bnds = np.array( [ [ -1, 0 ], [ 1, 0 ] ] )
    bnos = np.array( [ 1, 1 ] )

    u = Polyhedral( afds, afos, bnds, bnos )
    print( u )

    v = u.star()
    print( v )

afds = np.array( [ [ 0 ] ] )
afos = np.array( [ 0 ] )
bnds = np.array( [ [ -1 ], [ 1 ] ] )
bnos = np.array( [ 1, 1 ] )

u = Polyhedral( afds, afos, bnds, bnos )
print( u )

v = u.star()
print( v )

# v.plot1d()
