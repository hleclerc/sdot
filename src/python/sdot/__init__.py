# cpp source directory
def add_cpp_flags():
    from pathlib import Path
    import vfs

    vfs.add_global_cpp_flag( "-I" + str( Path( __file__ ).parents[ 3 ].joinpath( 'src', 'cpp' ) ) )
    vfs.add_global_cpp_flag( "--cpp-flag=-O3" )

add_cpp_flags()

# exports
from .ConvexFunction import ConvexFunction
from .PowerDiagram import PowerDiagram
