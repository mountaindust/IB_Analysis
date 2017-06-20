# Automate VisIt to analyze IBFE Lagrangian mesh data

Author: Christopher Strickland   
Email: cstric12@utk.edu   
Website: http://www.christopherstrickland.info   
GitHub: http://github.com/mountaindust   

## Automating VisIt for time-series plot creation

The relevant file is `IBFE_analysis_py27.py`. You will need to edit it with your 
particulars and run it under Python 2.7 (no later). An example of output is shown in `visit000.png`, which was generated using Laura's IBFE_jellyfish.

### Dependencies for automating plot creation in VisIt
- Python 2.7
- VisIt installed
- numpy
- matplotlib

### Instructions
1. Review the IBFE tutorial entitled "5-IBFE Tutorial: Visualizing Lagrangian Data" found in the Google Drive folder IBAMR Tutorials - Public > IBFE > IBFE Tutorials.
2. Edit the top of `IBFE_analysis_py27.py` with the path to your VisIt installation's lib/site-packages folder. This tells Python where to look for VisIt's Python interface.
3. You can specify the path to your data, the directory to put the results in, the variables to include in the analysis, and the type of time query to conduct at the top of `IBFE_analysis_py27.py`. Alternatively, you can also pass this information via command-line options when you run the program.
4. If you want to define any new variables to be included in the plot, do so in the function define_expressions, copying the example that is already there. These will automatically be added to the query over time plot by the code.
5. Edit the add_selection function to apply an operator to the Lagrangian mesh in order to exclude part of it from analysis. If you want to analyze the entire mesh, just comment all of the code inside this function and add a line with the command `pass` underneath (this tells Python that the function is intentionally empty and meant to do nothing).
6. Edit the apply\_plot\_attributes function if you would like the resulting plot curves to look different than they do in the example. You can change the color scheme by using different matplotlib colormaps in the top line.
7. Edit the apply\_annotations function to tweak the look of the plot as a whole, including the axes and plot title.

To run the program, just type `python IBFE_analysis_py27.py` into a terminal. Options are available to include different variables, change the query type, etc. without going back in and editing the code. Enter `python IBFE_analysis_py27.py -h` into a unix terminal to get help on what is available.

Example: `python IBFE_analysis_py27.py --vars X_0,X_1,X_2 --type "Max"` runs the "Max" query over time on variables "X\_0", "X\_1", and "X\_2", as well as any expressions defined in the function define\_expressions.

## Exporting Lagrangian data to VTK

The relevant file is `IBFE_to_vtk_py27.py`. You will need to edit it with your 
particulars and run in under Python 2.7 (no later). NOTE: The resulting VTK files 
are of the UNSTRUCTURED_GRID type and describe scalar values on a finite element mesh. 
While Python can technically read in all of this data, making sense of the 3D mesh 
by hand will not be easy. Use this file only if you have something specific in mind. 
Usage is similar to the instructions for `IBFE_analysis_py27.py` above.

If you really must read in the data to Python, use the vtk Python library and
follow my example for an unstructured grid. The cells will be each mesh element,
with the cell type describing their shape.
