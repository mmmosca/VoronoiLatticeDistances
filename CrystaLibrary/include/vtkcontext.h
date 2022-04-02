#ifndef _VTK_CONTEXT_H
#define _VTK_CONTEXT_H

#include <geom.h>
#include <linearcomplex.h>
#include <vector>
#include <Eigen/Dense>
#include <geometrystructure.h>

#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDataSetMapper.h>
#include <vtkIdList.h>
#include <vtkNamedColors.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkArrowSource.h>
#include <vtkMath.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolygon.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>
#include <vtkXMLPolyDataWriter.h>
//#include <vtkUnstructuredGrid.h>
//#include <vtkXMLUnstructuredGridWriter.h>

class VTKContext {
public:

	static vtkSmartPointer<vtkPolyData> getVTKPolygonalDataFromLinearComplexCombinatorialMap(Linear_Complex_Combinatorial_Map &lc_cm);

	static vtkSmartPointer<vtkPolyData> getVTKPolygonalData(std::vector<Point_CM> points, std::vector<std::vector<int>> faces);

	static vtkSmartPointer<vtkPolyData> getVTKPolygonalData(Polyhedron &P);

	static vtkSmartPointer<vtkPolyData> getVTKPolygonalData(Polyhedron_Ext &P);

	static vtkSmartPointer<vtkPolyData> getVTKPolygonalData(Nef_polyhedron_Ext &P);
	/*
	*	Draw with VTK (OpenGL2) the combinatorial map without the infinite vertex
	*/
	static void drawPolygonalData(vtkSmartPointer<vtkPolyData> polygonPolyData);

	/*
	*	Write the Combinatorial map into a .vtk file format
	*/
	static void writePolyDataInVTPFile(vtkSmartPointer<vtkPolyData> polyData, const char* filename);

	static void writePolyDataFilterInVTPFile(vtkSmartPointer<vtkTransformPolyDataFilter> polyDataFilter, const char* filename);

	static vtkSmartPointer<vtkPolyData> getPolyDataPoints(std::vector<Eigen::VectorXd> &points);

	static std::vector<vtkSmartPointer<vtkTransformPolyDataFilter>> getPolyDataAxis(std::vector<double> parameters);
};

#endif // _VTK_CONTEXT_H
