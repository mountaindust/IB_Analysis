# Examples of reading and writing VTK data

### Dependencies (these vary amoung included files)
- Python 3.5+ and Python 2.7
- numpy/scipy
- matplotlib
- vtk (to get this under Anaconda, run `conda install -c menpo vtk`)
- VisIt installed.

You will need a Python 2.7 environment with numpy and VisIt installed to convert IBAMR data into vtk data.

## vtk\_RW\_example
Reads and writes back out a few VTK files. This is primarily meant as an example
of how to write out VTK data using the vtk Python library. For details on 
reading VTK files in Python, see data_IO.

## data_IO
This was pulled from my Planktos repo and provides examples of functions for
reading VTK data from a variety of sources using the vtk Python library.

## read\_IBAMR3d\_py27
This will read in samrai 3D data using VisIt, resample it, and convert it to VTK. Must be run under Python 2.7 with VisIt installed. You must edit this file with the path to your VisIt installation's lib/site-packages directory before running. Help is available by running this script from the terminal with the -h option.

## collect\_avgs\_py27
This is a script I used (after adding a for-loop to do this over many files) for averaging flow velocity in 3D samrai data over 2D slices perpendicular to the Z axis. We were examining 3D flow around one or more cylinders of varying height, Re, and spacing. The averages in the X, Y, and Z direction and dp/dx are saved to a text file for later reading and plotting. plot\_avgs.py was used to read these text files and plot the averaged flow profile over various simulations - it has been included for completion.