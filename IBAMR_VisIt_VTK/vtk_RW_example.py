import numpy as np
import vtk
from vtk.util import numpy_support

### STRUCTURED_POINTS - SCALARS ###
reader = vtk.vtkStructuredPointsReader()
reader.SetFileName("fMag.0001.vtk")
reader.Update()
vtk_data = reader.GetOutput()
scalar_data = numpy_support.vtk_to_numpy(vtk_data.GetPointData().GetScalars())
mesh_shape = vtk_data.GetDimensions()
origin = vtk_data.GetOrigin()
spacing = vtk_data.GetSpacing()
scalar_data_vtk = numpy_support.numpy_to_vtk(scalar_data)
scalar_data_vtk.SetName('fMag')
vtk_obj = vtk.vtkStructuredPoints()
vtk_obj.SetDimensions(mesh_shape)
vtk_obj.SetOrigin(origin)
vtk_obj.SetSpacing(spacing)
vtk_obj.GetPointData().SetScalars(scalar_data_vtk)
writer = vtk.vtkGenericDataObjectWriter()
writer.SetFileName("test_out.vtk")
writer.SetInputDataObject(vtk_obj)
writer.Update()
writer.Write() # returns 1 on success, 0 on failure

### STRUCTURED_POINTS - VECTORS ###
reader = vtk.vtkStructuredPointsReader()
reader.SetFileName("u.0001.vtk")
reader.Update()
vtk_data = reader.GetOutput()
# return an n x 3 array of vectors
vector_data = numpy_support.vtk_to_numpy(vtk_data.GetPointData().GetVectors())
mesh_shape = vtk_data.GetDimensions()
origin = vtk_data.GetOrigin()
spacing = vtk_data.GetSpacing()
vector_data_vtk = numpy_support.numpy_to_vtk(vector_data)
vector_data_vtk.SetName("u")
vtk_obj = vtk.vtkStructuredPoints()
vtk_obj.SetDimensions(mesh_shape)
vtk_obj.SetOrigin(origin)
vtk_obj.SetSpacing(spacing)
vtk_obj.GetPointData().SetVectors(vector_data_vtk)
writer = vtk.vtkGenericDataObjectWriter()
writer.SetFileName("test_vector_out.vtk")
writer.SetHeader("SpamAndEggs!")
writer.SetInputDataObject(vtk_obj)
writer.Update()
writer.Write()

### UNSTRUCTURED_GRID (lagsPts) ###
reader = vtk.vtkUnstructuredGridReader()
reader.SetFileName("lagsPts.0001.vtk")
reader.Update()
vtk_data = reader.GetOutput()
points = numpy_support.vtk_to_numpy(vtk_data.GetPoints().GetData()) # 2D array
cells = numpy_support.vtk_to_numpy(vtk_data.GetCells().GetData()) # 1D array
cell_types = numpy_support.vtk_to_numpy(vtk_data.GetCellTypesArray())
# We need new datatypes here
vtk_points = vtk.vtkPoints()
vtk_points.SetData(numpy_support.numpy_to_vtk(points))
vtk_cells = vtk.vtkCellArray()
vtk_cells.SetCells(cells.size//2, numpy_support.numpy_to_vtkIdTypeArray(cells))
vtk_obj = vtk.vtkUnstructuredGrid()
vtk_obj.SetCells(np.ones(cells.size//2, dtype=int), vtk_cells)
vtk_obj.SetPoints(vtk_points)
# write out
writer = vtk.vtkGenericDataObjectWriter()
writer.SetFileName("test_lagpts_out.vtk")
writer.SetHeader("SpamAndEggs!")
writer.SetInputDataObject(vtk_obj)
writer.Update()
writer.Write()
