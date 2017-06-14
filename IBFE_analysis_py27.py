#! /usr/bin/env python2

''' This script interacts with VisIt to analyze Lagrangian Data from IBFE. It
automates the tutorial for analyzing IBFE Lagrangian data in time in VisIt found
in the Google Drive folder, IBAMR Tutorials - Public > IBFE > IBFE Tutorials.
Run this file directly through a Python 2.7 interpreter.'''

from __future__ import division # make float division the default
import sys
import os.path
############# Edit with proper VisIt path!!! #############
sys.path.append("C:\Program Files\LLNL\VisIt 2.11.0\lib\site-packages")
#############
from visit import *
# visit_utils has several nice helper functions, among them "query" which just
#   returns the numerical output of Query.
from visit_utils import *
import numpy as np
import matplotlib.pyplot as plt

##### Full path to output.ex2 #####
path_to_data = r"IBFE_data/output.ex2"
##### Directory in which to save resulting plots #####
out_directory = r"D:\Python\VisIt_IBFE"
#####
assert os.path.isfile(path_to_data)

def main(path_to_data=path_to_data, out_directory=out_directory):
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

    ### Define additional variables to plot as Expressions ###
    DefineScalarExpression("radial", "sqrt(dX_0*dX_0+dX_1*dX_1)")

    # Add and draw the Pseudocolor plot
    AddPlot("Pseudocolor", "X_0", 1, 0)
    DrawPlots()

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
    DrawPlots()
    ###

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

    SetActivePlots(0) # Highlight the current quantity
    # "Max", "Min", or "Average Value"
    QueryOverTime("Max", end_time=320, start_time=0, stride=1)
    # This is experimental
    print(GetQueryOutputValue())

    ### Record Annotations and put here! ###

    ###

    ### Save finished plot ###
    SetActiveWindow(2) # The analysis plot is window 2
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
    # If no arguments are passed on the command line, launch visit and run main()
    if len(sys.argv[1:]) == 0:
        Launch()
        #LaunchNowin()
        main()
    else:
        # Expect a string argument which gives the text filename of previously
        #   generated data from this file. Load and plot it.
        # data = np.loadtxt(sys.argv[1])
        # plot_avgs(data[0,:],data[1,:],data[2,:],data[3,:])
        pass