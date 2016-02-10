#!/usr/bin/env python

# takes an HDF5 file generated from the SMPM solver and produces an XDMF
# document representing its contents.

import getopt
import h5py
import sys

def usage( script_name ):
    """
    Takes a name of the script (full path, name, etc) and prints its usage to standard
    output.
    """

    print """Usage: %s [-h] <SMPM HDF5 file>

Takes an HDF5 file produced by the SMPM Incompressible Navier-Stokes solver and
generates an XDMF file on standard output that describes its contents.  The
XDMF file allows the HDF5 file to be visualized in VTK-based tools like ParaView,
Visit, and MayaVi.

The command line options above are described below:

  -h        Displays this help message and exits.
""" % script_name

def main( argv ):
    """
    Takes a list of arguments, including the name of the application calling this method,
    and generates an XDMF file on standard output.
    """

    # parse our command line options.
    try:
        opts, args = getopt.getopt( argv[1:], "h" )
    except getopt.GetoptError as error:
        sys.stderr.write( "Error processing option: %s\n" % str( error ) )
        sys.exit( 1 )

    # handle any valid options were were presented.
    for opt, arg in opts:
        if opt == '-h':
            usage( argv[0] )
            sys.exit()

    if len( args ) != 1:
        sys.stderr.write( "Expected 1 argument but received %d.\n" % len( args ) )
        print
        usage( argv[0] )
        sys.exit( 1 )

    # open the supplied file and pull out the variables we need to describe
    # its contents.
    file_name = args[0]
    h5 = h5py.File( file_name, "r" )

    n                 = h5['/grid/n'].value
    mx                = h5['/grid/mx'].value * n
    my                = h5['/grid/my'].value * n
    mz                = h5['/grid/mz'].value * n
    number_time_steps = h5['/field/number_steps'].value

    h5.close()

    # create a collection of temporal grids and define a default topology and
    # geometry.  each timestep's grid will reference this toplogy and geometry
    # via an xpointer path.
    format_string = """<?xml version="1.0" ?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>
<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">
 <Domain>
   <Grid Name="mesh1" CollectionType="Temporal" GridType="Collection">

     <Topology TopologyType="3DSMesh" NumberOfElements="%d %d %d"/>
     <Geometry GeometryType="X_Y_Z">
       <DataItem Dimensions="%d %d %d" NumberType="Float" Precision="8" Format="HDF">
         %s:/grid/x
       </DataItem>
       <DataItem Dimensions="%d %d %d" NumberType="Float" Precision="8" Format="HDF">
         %s:/grid/y
       </DataItem>
       <DataItem Dimensions="%d %d %d" NumberType="Float" Precision="8" Format="HDF">
         %s:/grid/z
       </DataItem>
     </Geometry>
"""

    print format_string % (mz, my, mx, mz, my, mx, file_name, mz, my, mx, file_name, mz, my, mx, file_name )

    # define a grid of scalar values (ux, uy, uz, and s) for each timestep written
    # by the solver.  the topology and geometry points back to the global
    # versions.
    format_string = """
     <Grid Name="time%d">
       <Time Value="%d" />

       <xi:include xpointer="element(/1/1/1/1)"/>
       <xi:include xpointer="element(/1/1/1/2)"/>

       <Attribute Name="ux" AttributeType="Scalar" Center="Cell">
       <DataItem DataType="Float" Precision="8" Dimensions="%d %d %d" Format="HDF">
         %s:/field/step%d/ux
       </DataItem>
       </Attribute>
       <Attribute Name="uy" AttributeType="Scalar" Center="Cell">
       <DataItem DataType="Float" Precision="8" Dimensions="%d %d %d" Format="HDF">
         %s:/field/step%d/uy
       </DataItem>
       </Attribute>
       <Attribute Name="uz" AttributeType="Scalar" Center="Cell">
       <DataItem DataType="Float" Precision="8" Dimensions="%d %d %d" Format="HDF">
         %s:/field/step%d/uz
       </DataItem>
       </Attribute>
       <Attribute Name="s" AttributeType="Scalar" Center="Cell">
       <DataItem DataType="Float" Precision="8" Dimensions="%d %d %d" Format="HDF">
         %s:/field/step%d/s
       </DataItem>
       </Attribute>
     </Grid>

"""

    for step_number in range( number_time_steps ):
        print format_string % (step_number + 1,
                               step_number + 1,
                               mz, my, mx,
                               file_name, step_number + 1,
                               mz, my, mx,
                               file_name, step_number + 1,
                               mz, my, mx,
                               file_name, step_number + 1,
                               mz, my, mx,
                               file_name, step_number + 1)

    # make the document well formed.
    print """
   </Grid>
 </Domain>
</Xdmf>
"""

if __name__ == "__main__":
    main( sys.argv )
