#! /usr/bin/env python2

''' This script interacts with VisIt to analyze Lagrangian Data from IBFE. It
automates the tutorial for analyzing IBFE Lagrangian data in time in VisIt found
in the Google Drive folder, IBAMR Tutorials - Public > IBFE > IBFE Tutorials.
Run this file directly through a Python 2.7 interpreter.'''

from __future__ import division # make float division the default
import sys
import os
import argparse
############# Edit with proper VisIt path!!! #############
sys.path.append("C:\Program Files\LLNL\VisIt 2.11.0\lib\site-packages")
#############
from visit import *
# visit_utils has several nice helper functions, among them "query" which just
#   returns the numerical output of Query.
from visit_utils import *
import numpy as np
import matplotlib.pyplot as plt

### Full default path to output.ex2 ###
path_to_data = r"IBFE_data/output.ex2"
### Default directory in which to save resulting plots ###
out_directory = r"."
###

# Arparse setup, to take options from the command line
parser = argparse.ArgumentParser(description="Analyze IBFE Lagrangian data.")
parser.add_argument('-f', '--filename', type=str, default=path_to_data,
                    help="location of data to load")
parser.add_argument('-o', '--outdir', type=str, default=out_directory,
                    help="dir for output")



def define_expressions(var_list):
    ''' Define any expressions (new variables, not in the dataset) that you want
    plots of here.
    '''

    ### Define additional variables to plot/analyze as Expressions ###
    DefineScalarExpression("radial", "sqrt(dX_0*dX_0+dX_1*dX_1)")
    ### Append their names to var_list ###
    var_list.append("radial")



def add_selection(var=None):
    ''' If you only want to analyze a subset of the Lagrangian mesh, specify
    the selection here. You can make this dependent upon the variable currently
    under consideration by using if-statements along with the name of the variable
    (which is passed in as the "var" argument, a string).
    '''

    ### Create a box selection to restrict the analysis ###
    AddOperator("Box", 0)
    BoxAtts = BoxAttributes()
    BoxAtts.amount = BoxAtts.Some  # Some, All
    BoxAtts.minx = -0.02
    BoxAtts.maxx = 0.02
    BoxAtts.miny = -0.02
    BoxAtts.maxy = 0.02
    BoxAtts.minz = 0
    BoxAtts.maxz = 0.004
    BoxAtts.inverse = 0
    SetOperatorOptions(BoxAtts, 0)
    ###


def apply_annotations(var=None):
    ''' Make the plot prettier here. If you want different effects for different
    variables, you can use if-statements along with the name of the variable
    (which is passed in as the "var" argument, a string).
    '''

    ### Record Annotations and put here! ###

    ###
    pass


def main(path_to_data=path_to_data, out_directory=out_directory,
         var_list=["X_0","X_1"]):
    '''Run the analysis'''

    # Check for proper paths
    assert os.path.isfile(path_to_data), "Could not find %s." % path_to_data
    if not os.path.isdir(out_directory):
        os.mkdir(out_directory)

    ### This bit just sets a default view of some sort... not really sure ###
    # (I just dumped this straight in from recording it in VisIt)
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
    ###

    # Open the database
    OpenDatabase("localhost:"+path_to_data, 0)

    ### Set the Query over time options ###
    QueryOverTimeAtts = GetQueryOverTimeAttributes()
    QueryOverTimeAtts.timeType = QueryOverTimeAtts.DTime  # Cycle, DTime, Timestep
    QueryOverTimeAtts.startTimeFlag = 0
    QueryOverTimeAtts.startTime = 0
    QueryOverTimeAtts.endTimeFlag = 0
    QueryOverTimeAtts.endTime = 1
    QueryOverTimeAtts.strideFlag = 0
    QueryOverTimeAtts.stride = 1
    QueryOverTimeAtts.createWindow = 1
    QueryOverTimeAtts.windowId = 2
    SetQueryOverTimeAttributes(QueryOverTimeAtts)
    ###

    # Define any additional variables to add to the analysis
    define_expressions(var_list)

    for n,var in enumerate(var_list):
        # Add and draw the Pseudocolor plot
        AddPlot("Pseudocolor", var, 1, 0)
        DrawPlots()

        ### Create a box selection to restrict the analysis ###
        add_selection()
        DrawPlots()
        ###

        SetActivePlots(n) # Highlight the current quantity
        ### Create the query, e.g. "Max", "Min", or "Average Value" ###
        QueryOverTime("Max", end_time=320, start_time=0, stride=1)

        SetActiveWindow(1)
        HideActivePlots()

    # make the plot pretty
    SetActiveWindow(2) # The analysis plot is window 2
    apply_annotations()

    ### Save finished plot ###
    SaveWindowAtts = SaveWindowAttributes()
    SaveWindowAtts.outputToCurrentDirectory = 0
    SaveWindowAtts.outputDirectory = out_directory
    SaveWindowAtts.fileName = "visit"
    SaveWindowAtts.family = 1
    SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
    # resolution of resulting figure
    SaveWindowAtts.width = 1024
    SaveWindowAtts.height = 1024
    #
    SaveWindowAtts.screenCapture = 0
    SaveWindowAtts.saveTiled = 0
    SaveWindowAtts.quality = 80
    SaveWindowAtts.progressive = 0
    SaveWindowAtts.binary = 0
    SaveWindowAtts.stereo = 0
    SaveWindowAtts.compression = SaveWindowAtts.PackBits  # None, PackBits, Jpeg, Deflate
    SaveWindowAtts.forceMerge = 0
    SaveWindowAtts.resConstraint = SaveWindowAtts.ScreenProportions  # NoConstraint, EqualWidthHeight, ScreenProportions
    SaveWindowAtts.advancedMultiWindowSave = 0
    SetSaveWindowAttributes(SaveWindowAtts)
    SaveWindow()
    ###

if __name__ == '__main__':
    Launch()
    args = parser.parse_args()
    main(args.filename, args.outdir)
