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
sys.path.append(r"C:\Program Files\LLNL\VisIt 2.11.0\lib\site-packages")
#############
from visit import *
# visit_utils has several nice helper functions, among them "query" which just
#   returns the numerical output of Query.
from visit_utils import *
import numpy as np
from matplotlib import cm

### Full default path to output.ex2 ###
path_to_data = r"IBFE_data/output.ex2"
### Default directory in which to save resulting plots ###
out_directory = r"."
### Default variables to analyze
var_list = ["X_0","X_1"]
### Query Over Time type (e.g., "Max", "Min", or "Average Value")
q_type = "Average Value"

# Arparse setup, to take options from the command line
parser = argparse.ArgumentParser(description="Analyze IBFE Lagrangian data.")
parser.add_argument('-f', '--filename', type=str, default=path_to_data,
                    help="location of data to load")
parser.add_argument('-o', '--outdir', type=str, default=out_directory,
                    help="dir for output")
parser.add_argument('--vars', type=str, default=','.join(var_list),
                    help="variables in the dataset to analyze.\n"+
                    "Specify as a comma separated list, no spaces!")
parser.add_argument('--type', type=str, default=q_type,
                    help='type of Query Over Time to run,\n'+
                    'e.g. Max, Min, or "Average Value". If the type\n'+
                    'has a space in it, be sure to enclose in quotations!')



def define_expressions(var_list):
    ''' Define any expressions (new variables, not already in the dataset) 
    that you want plots of here.
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



def apply_plot_attributes(var, clr_x):
    ''' Make the individual plot lines prettier here. If you want different 
    effects for different variables, you can use if-statements along with the 
    name of the variable (which is passed in as the "var" argument, a string).
    clr_x is a number between 0 and 1 to apply to a matplotlib colormap
    '''
    # Something about setting these curve attributes seems to kill the
    #   VisIt color cycler. So we do one better by replacing it with a
    #   matplotlib color cycler!
    # Good ones to try include Accent, Dark2, Set1, Set2, Vega10, spectral...

    this_color = tuple(int(val*255) for val in cm.Accent(clr_x))

    CurveAtts = CurveAttributes()
    ### set the curve color
    CurveAtts.curveColorSource = CurveAtts.Custom  # Cycle, Custom
    CurveAtts.curveColor = this_color # RGBalpha
    ### turn off/on legend and labels
    CurveAtts.showLegend = 1 # legend on/off
    CurveAtts.showLabels = 0 # in-plot labels on/off
    ### line properties
    CurveAtts.showLines = 1 # on/off
    CurveAtts.lineStyle = CurveAtts.SOLID  # SOLID, DASH, DOT, DOTDASH
    CurveAtts.lineWidth = 2
    ### point properties (plot data points along the line)
    CurveAtts.showPoints = 0 # on/off
    CurveAtts.symbol = CurveAtts.Point  # Point, TriangleUp, TriangleDown, Square, Circle, Plus, X
    CurveAtts.pointSize = 5
    CurveAtts.pointFillMode = CurveAtts.Static  # Static, Dynamic
    CurveAtts.pointStride = 1 # plot every n data-points
    CurveAtts.symbolDensity = 50
    ### fill below the plot line to the x-axis
    CurveAtts.fillMode = CurveAtts.NoFill  # NoFill, Solid, HorizontalGradient, VerticalGradient
    CurveAtts.fillColor1 = (0, 255, 0, 255)
    CurveAtts.fillColor2 = (100, 255, 100, 255) # used to make a gradient effect
    SetPlotOptions(CurveAtts)



def apply_annotations():
    ''' Make the plot as a whole prettier here.'''

    AnnotationAtts = AnnotationAttributes()
    AnnotationAtts.userInfoFlag = 0 # turn off user infomation
    AnnotationAtts.databaseInfoFlag = 0 # turn off database title
    AnnotationAtts.timeInfoFlag = 0
    AnnotationAtts.legendInfoFlag = 1 # turn on legend
    AnnotationAtts.backgroundColor = (255, 255, 255, 255) # RGBalpha
    AnnotationAtts.foregroundColor = (0, 0, 0, 255) # RGBalpha
    # Axes properties
    AnnotationAtts.axes2D.visible = 1
    AnnotationAtts.axes2D.autoSetTicks = 1
    AnnotationAtts.axes2D.autoSetScaling = 1
    AnnotationAtts.axes2D.lineWidth = 0 # for axes box
    AnnotationAtts.axes2D.tickLocation = AnnotationAtts.axes2D.Outside  # Inside, Outside, Both
    AnnotationAtts.axes2D.tickAxes = AnnotationAtts.axes2D.BottomLeft  # Off, Bottom, Left, BottomLeft, All
    # x-axis title properties
    AnnotationAtts.axes2D.xAxis.title.visible = 1
    AnnotationAtts.axes2D.xAxis.title.font.font = AnnotationAtts.axes2D.xAxis.title.font.Courier  # Arial, Courier, Times
    AnnotationAtts.axes2D.xAxis.title.font.scale = 1
    AnnotationAtts.axes2D.xAxis.title.font.useForegroundColor = 1
    AnnotationAtts.axes2D.xAxis.title.font.color = (0, 0, 0, 255)
    AnnotationAtts.axes2D.xAxis.title.font.bold = 1
    AnnotationAtts.axes2D.xAxis.title.font.italic = 1
    AnnotationAtts.axes2D.xAxis.title.userTitle = 0
    AnnotationAtts.axes2D.xAxis.title.title = "X-Axis" # Only used if userTitle = 1
    AnnotationAtts.axes2D.xAxis.title.userUnits = 1 # suppress uniformative "time" units
    AnnotationAtts.axes2D.xAxis.title.units = "" # Only used if userUnits = 1
    # x-label properties
    AnnotationAtts.axes2D.xAxis.label.visible = 1
    AnnotationAtts.axes2D.xAxis.label.font.font = AnnotationAtts.axes2D.xAxis.label.font.Courier  # Arial, Courier, Times
    AnnotationAtts.axes2D.xAxis.label.font.scale = 1
    AnnotationAtts.axes2D.xAxis.label.font.useForegroundColor = 1
    AnnotationAtts.axes2D.xAxis.label.font.color = (0, 0, 0, 255)
    AnnotationAtts.axes2D.xAxis.label.font.bold = 1
    AnnotationAtts.axes2D.xAxis.label.font.italic = 1
    AnnotationAtts.axes2D.xAxis.label.scaling = 0
    # x-tickmark properties
    AnnotationAtts.axes2D.xAxis.tickMarks.visible = 1
    AnnotationAtts.axes2D.xAxis.tickMarks.majorMinimum = 0
    AnnotationAtts.axes2D.xAxis.tickMarks.majorMaximum = 1
    AnnotationAtts.axes2D.xAxis.tickMarks.minorSpacing = 0.02
    AnnotationAtts.axes2D.xAxis.tickMarks.majorSpacing = 0.2
    AnnotationAtts.axes2D.xAxis.grid = 0 # x-grid on/off
    # y-axis title properties
    AnnotationAtts.axes2D.yAxis.title.visible = 1
    AnnotationAtts.axes2D.yAxis.title.font.font = AnnotationAtts.axes2D.yAxis.title.font.Courier  # Arial, Courier, Times
    AnnotationAtts.axes2D.yAxis.title.font.scale = 1
    AnnotationAtts.axes2D.yAxis.title.font.useForegroundColor = 1
    AnnotationAtts.axes2D.yAxis.title.font.color = (0, 0, 0, 255)
    AnnotationAtts.axes2D.yAxis.title.font.bold = 1
    AnnotationAtts.axes2D.yAxis.title.font.italic = 1
    AnnotationAtts.axes2D.yAxis.title.userTitle = 0
    AnnotationAtts.axes2D.yAxis.title.title = "Y-Axis" # Only used if userTitle = 1
    AnnotationAtts.axes2D.yAxis.title.userUnits = 0
    AnnotationAtts.axes2D.yAxis.title.units = "" # Only used if userUnits = 1
    # y-label properties
    AnnotationAtts.axes2D.yAxis.label.visible = 1
    AnnotationAtts.axes2D.yAxis.label.font.font = AnnotationAtts.axes2D.yAxis.label.font.Courier  # Arial, Courier, Times
    AnnotationAtts.axes2D.yAxis.label.font.scale = 1
    AnnotationAtts.axes2D.yAxis.label.font.useForegroundColor = 1
    AnnotationAtts.axes2D.yAxis.label.font.color = (0, 0, 0, 255)
    AnnotationAtts.axes2D.yAxis.label.font.bold = 1
    AnnotationAtts.axes2D.yAxis.label.font.italic = 1
    AnnotationAtts.axes2D.yAxis.label.scaling = 0
    # y-tickmark properties
    AnnotationAtts.axes2D.yAxis.tickMarks.visible = 1
    AnnotationAtts.axes2D.yAxis.tickMarks.majorMinimum = 0
    AnnotationAtts.axes2D.yAxis.tickMarks.majorMaximum = 1
    AnnotationAtts.axes2D.yAxis.tickMarks.minorSpacing = 0.02
    AnnotationAtts.axes2D.yAxis.tickMarks.majorSpacing = 0.2
    AnnotationAtts.axes2D.yAxis.grid = 0 # y-grid on/off
    SetAnnotationAttributes(AnnotationAtts)
    
    ### Create a plot title ###
    title = CreateAnnotationObject("Text2D")
    title.SetText('Data analysis')
    title.SetPosition(0.46,0.95)



def main(path_to_data=path_to_data, out_directory=out_directory,
         var_list=var_list,q_type=q_type):
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

    # Get spacing for matplotlib colors
    color_idx = np.linspace(0, 1, len(var_list))

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
        QueryOverTime(q_type, end_time=320, start_time=0, stride=1)

        ### Change how the output plot looks
        SetActiveWindow(2)
        SetActivePlots(n)
        apply_plot_attributes(var, color_idx[n])

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
    # Resolution of resulting figure. visit seems to mostly ignore what you
    #   put here if you are asking for resolutions bigger than your screen
    #   (this also happens in the GUI). If you need publication quality stuff, 
    #   you are best off trying postscript.
    SaveWindowAtts.resConstraint = SaveWindowAtts.NoConstraint  # NoConstraint, EqualWidthHeight, ScreenProportions 
    SaveWindowAtts.width = 3000
    SaveWindowAtts.height = 3000
    #
    SaveWindowAtts.screenCapture = 0
    SaveWindowAtts.saveTiled = 0
    SaveWindowAtts.quality = 100 # 0 - 100
    SaveWindowAtts.progressive = 0
    SaveWindowAtts.binary = 0
    SaveWindowAtts.stereo = 0
    SaveWindowAtts.compression = SaveWindowAtts.Jpeg  # None, PackBits, Jpeg, Deflate
    SaveWindowAtts.forceMerge = 0
    SaveWindowAtts.advancedMultiWindowSave = 0
    SetSaveWindowAttributes(SaveWindowAtts)
    SaveWindow()
    ###

    CloseDatabase("localhost:"+path_to_data)

if __name__ == '__main__':
    args = parser.parse_args()
    Launch()
    main(args.filename, args.outdir, args.vars.split(','), args.type)
