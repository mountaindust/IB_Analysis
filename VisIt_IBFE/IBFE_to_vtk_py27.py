#! /usr/bin/env python2

'''This script reads in IBFE Lagrangian data using VisIt and spits out VTK files.
NOTE: The resulting VTK files are of the UNSTRUCTURED_GRID type and describe
scalar values on a finite element mesh. While Python can technically read in all
of this data, making sense of the 3D mesh by hand will not be easy. Use this file
only if you have something specific in mind.

If you really must read in the data to Python, use the vtk Python library and
follow my example for an unstructured grid. The cells will be each mesh element,
with the cell type describing their shape.

Created on Tues June 14 2017

Author: Christopher Strickland
Email: cstric12@utk.edu
'''

import sys, os
import argparse
############# Edit with proper VisIt path!!! #############
sys.path.append(r"C:\Program Files\LLNL\VisIt 2.11.0\lib\site-packages")
#############
from visit import *

### Full default path to output.ex2 ###
path_to_data = r"IBFE_data/output.ex2"
### Default directory in which to save vtk files ###
out_directory = r"./vtk"
### Default variables to convert
var_list = ["X_0","X_1"]

parser = argparse.ArgumentParser(description="Convert IBFE Lagrangian data to VTK.")
parser.add_argument('-f', '--filename', type=str, default=path_to_data,
                    help="location of data to load")
parser.add_argument('-o', '--outdir', type=str, default=out_directory,
                    help="dir for output")
parser.add_argument('--vars', type=str, default=','.join(var_list),
                    help="variables in the dataset to convert.\n"+
                    "Specify as a comma separated list, no spaces!")

def main(path_to_data, out_directory, var_list):
    '''Read in Lagrangian data located at path_to_data and convert the variables
    in var_list to VTK files, storing them in out_directory.'''

    # Check for proper paths
    assert os.path.isfile(path_to_data), "Could not find %s." % path_to_data
    if not os.path.isdir(out_directory):
        os.mkdir(out_directory)

    # Begin spontaneous state
    View2DAtts = View2DAttributes()
    View2DAtts.windowCoords = (0, 1, 0, 1)
    View2DAtts.viewportCoords = (0.2, 0.95, 0.15, 0.95)
    View2DAtts.fullFrameActivationMode = View2DAtts.Auto  # On, Off, Auto
    View2DAtts.fullFrameAutoThreshold = 100
    View2DAtts.xScale = View2DAtts.LINEAR  # LINEAR, LOG
    View2DAtts.yScale = View2DAtts.LINEAR  # LINEAR, LOG
    View2DAtts.windowValid = 0
    SetView2D(View2DAtts)
    # End spontaneous state

    # Open the database
    OpenDatabase("localhost:"+path_to_data, 0)

    for var in var_list:
        # Add and draw the Pseudocolor plot
        AddPlot("Pseudocolor", var, 0, 1)
        DrawPlots()

        ExportDBAtts = ExportDBAttributes()
        ExportDBAtts.allTimes = 1
        ExportDBAtts.db_type = "VTK"
        ExportDBAtts.db_type_fullname = "VTK_1.0"
        ExportDBAtts.filename = var+"_"
        ExportDBAtts.dirname = out_directory
        ExportDBAtts.variables = (var)
        ExportDatabase(ExportDBAtts)
        DeleteAllPlots()
    CloseDatabase("localhost:"+path_to_data)



if __name__ == '__main__':
    args = parser.parse_args()
    Launch()
    main(args.filename, args.outdir, args.vars.split(','))