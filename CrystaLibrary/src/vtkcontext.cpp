/*
 * Permission is granted to copy, distribute and/or modify the documents
 * in this directory and its subdirectories unless otherwise stated under
 * the terms of the GNU Free Documentation License, Version 1.1 or any later version 
 * published by the Free Software Foundation; with no Invariant Sections, 
 * no Front-Cover Texts and no Back-Cover Texts. A copy of the license 
 * is available at the website of the GNU Project.
 * The programs and code snippets in this directory and its subdirectories
 * are free software; you can redistribute them and/or modify it under the 
 * terms of the GNU General Public License as published by the Free Software 
 * Foundation; either version 2 of the License, or (at your option) any later
 * version. This code is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * 
 * Author Marco M. Mosca, email: marcomichele.mosca@gmail.com
*/
#include "vtkcontext.h"

vtkSmartPointer<vtkPolyData> VTKContext::getVTKPolygonalDataFromLinearComplexCombinatorialMap(Linear_Complex_Combinatorial_Map &lc_cm) {
	Linear_Complex_Combinatorial_Map::size_type mark_todraw = lc_cm.get_new_mark();
	// Map Delaunay Point -> unique index
	std::map<Point_CM, int> point_to_index;
	std::map<Point_CM, int>::iterator it_point_to_index;
	std::vector<vtkIdType> pointIds;
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	std::vector<std::vector<vtkIdType>> facesIds;
	std::vector<vtkSmartPointer<vtkPolygon>> polygons;
	vtkSmartPointer<vtkCellArray> faces = vtkSmartPointer<vtkCellArray>::New();
	// It is the new index of the current vertex (chosen by me) 
	vtkIdType curr_v_index = 0;

	//std::cout << "Points from VTK context" << std::endl;
	// Extract the vertices, update the VTK context by marking the starting edge for every 2-cell
	// Iterate every face (2-cell): one dart for every face
	for (Linear_Complex_Combinatorial_Map::One_dart_per_cell_range<2>::iterator face_dart = lc_cm.one_dart_per_cell<2>().begin(),
		face_dart_end = lc_cm.one_dart_per_cell<2>().end(); face_dart != face_dart_end; ++face_dart)
	{
		vtkSmartPointer<vtkPolygon> face = vtkSmartPointer<vtkPolygon>::New();
		std::vector<vtkIdType> faceIds;
		Dart_handle dh_edges = face_dart;
		// Iterate every edge for every 2-cell
		while (!lc_cm.is_marked(dh_edges, mark_todraw)) {
			Point_CM p = lc_cm.point(dh_edges);
			// Check if the point exists in the map
			it_point_to_index = point_to_index.find(p);
			// Does not exist: The vertex has not been traversed before
			if (it_point_to_index == point_to_index.end()) {
				if (dh_edges == face_dart) {
					CGAL::mark_cell<Linear_Complex_Combinatorial_Map, 1>(lc_cm, dh_edges, mark_todraw);
				}
				// It will have the current index
				pointIds.push_back(curr_v_index);
				// Insert the vertex (point) and its index into the map
				point_to_index[p] = curr_v_index;
				// Insert the point in the VTKPoint structure
				points->InsertNextPoint(p[0], p[1], p[2]);
				//std::cout << "New [" << curr_v_index << "]: " << p[0] << " " << p[1] << " " << p[2] << std::endl;
				// Add the current vertex ID to the face IDs
				faceIds.push_back(curr_v_index);
				++curr_v_index;
				dh_edges = lc_cm.beta(dh_edges, 1);
			}
			// The vertex exists: Push its ID in the face indexes
			else {
				if (dh_edges == face_dart) {
					CGAL::mark_cell<Linear_Complex_Combinatorial_Map, 1>(lc_cm, dh_edges, mark_todraw);
				}
				// std::cout << "Old [" << it_point_to_index->second << "]: " << it_point_to_index->first[0] << " " << it_point_to_index->first[1] << " " << it_point_to_index->first[2] << std::endl;
				// Just push its index from the map to the face IDs
				faceIds.push_back(it_point_to_index->second);
				dh_edges = lc_cm.beta(dh_edges, 1);
			}
		}
		CGAL::unmark_cell<Linear_Complex_Combinatorial_Map, 1>(lc_cm, dh_edges, mark_todraw);
		// Set the Point IDs for the current face
		face->GetPointIds()->SetNumberOfIds(faceIds.size());
		//std::cout << "Face IDS: ";
		for (int i = 0; i < faceIds.size(); ++i) {
			face->GetPointIds()->SetId(i, faceIds[i]);
			//std::cout << faceIds[i] << " ";
		}
		//std::cout << std::endl;
		// Add the current face Ids to the container
		facesIds.push_back(faceIds);
		// Add the face to the poly data container
		polygons.push_back(face);
		//std::cout << "Number of Face points: " << faceIds.size() << std::endl;
	}

	// Add all faces to the Faces collection (Polygons array)
	for (auto &poly : polygons) {
		faces->InsertNextCell(poly);
	}
	//std::cout << "Number of Points: " << pointIds.size() << std::endl;
	//std::cout << "Number of faces: " << polygons.size() << std::endl;

	// Create a PolyData
	vtkSmartPointer<vtkPolyData> polygonPolyData = vtkSmartPointer<vtkPolyData>::New();
	polygonPolyData->SetPoints(points);
	polygonPolyData->SetPolys(faces);

	return polygonPolyData;
}

vtkSmartPointer<vtkPolyData> VTKContext::getVTKPolygonalData(std::vector<Point_CM> points, std::vector<std::vector<int>> faces)
{
	vtkSmartPointer<vtkPoints> vtk_points = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> vtk_faces = vtkSmartPointer<vtkCellArray>::New();

	for (Point_CM &p : points) {
		vtk_points->InsertNextPoint(CGAL::to_double(p[0]), CGAL::to_double(p[1]), CGAL::to_double(p[2]));
	}

	for (std::vector<int>& f : faces)
	{
		vtkSmartPointer<vtkPolygon> face = vtkSmartPointer<vtkPolygon>::New();
		// Set the number of Point IDs for the current face
		face->GetPointIds()->SetNumberOfIds(f.size());
		//std::cout << "Face IDS: ";
		for (int i = 0; i < f.size(); ++i) {
			face->GetPointIds()->SetId(i, f[i]);
		}
		// Add the face to face container
		vtk_faces->InsertNextCell(face);
	}

	// Create a PolyData
	vtkSmartPointer<vtkPolyData> polygonPolyData = vtkSmartPointer<vtkPolyData>::New();
	polygonPolyData->SetPoints(vtk_points);
	polygonPolyData->SetPolys(vtk_faces);

	return polygonPolyData;
}

vtkSmartPointer<vtkPolyData> VTKContext::getVTKPolygonalData(Polyhedron &P) {
	Polyhedron::Facet_iterator facet_it;
	// Map Poly Point -> unique index
	std::map<Kernel::Point_3, int> point_to_index;
	std::map<Kernel::Point_3, int>::iterator it_point_to_index;
	std::vector<vtkIdType> pointIds;
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	std::vector<std::vector<vtkIdType>> facesIds;
	std::vector<vtkSmartPointer<vtkPolygon>> polygons;
	vtkSmartPointer<vtkCellArray> faces = vtkSmartPointer<vtkCellArray>::New();
	// It is the new index of the current vertex (chosen by me) 
	vtkIdType curr_v_index = 0;

	//std::cout << "Points from VTK context" << std::endl;
	// Extract the vertices, update the VTK context
	// Iterate every face (2-cell)

	for (Polyhedron::Facet_iterator f_it = P.facets_begin(); f_it != P.facets_end(); ++f_it) {
		Polyhedron::Halfedge_around_facet_circulator v_circ;
		vtkSmartPointer<vtkPolygon> face = vtkSmartPointer<vtkPolygon>::New();
		std::vector<vtkIdType> faceIds;

		v_circ = f_it->facet_begin();
		do {
			Kernel::Point_3 p = v_circ->vertex()->point();
			// Check if the point exists in the map
			it_point_to_index = point_to_index.find(p);
			// Does not exist: The vertex has not been traversed before
			if (it_point_to_index == point_to_index.end()) {
				// It will have the current index
				pointIds.push_back(curr_v_index);
				// Insert the vertex (point) and its index into the map
				point_to_index[p] = curr_v_index;
				// Insert the point in the VTKPoint structure
				points->InsertNextPoint(CGAL::to_double(p[0]), CGAL::to_double(p[1]), CGAL::to_double(p[2]));
				//std::cout << "New [" << curr_v_index << "]: " << p[0] << " " << p[1] << " " << p[2] << std::endl;
				// Add the current vertex ID to the face IDs
				faceIds.push_back(curr_v_index);
				++curr_v_index;
			}
			// The vertex exists: Push its ID in the face indexes
			else {
				// std::cout << "Old [" << it_point_to_index->second << "]: " << it_point_to_index->first[0] << " " << it_point_to_index->first[1] << " " << it_point_to_index->first[2] << std::endl;
				// Just push its index from the map to the face IDs
				faceIds.push_back(it_point_to_index->second);
			}
		} while (++v_circ != f_it->facet_begin());

		// Set the Point IDs for the current face
		face->GetPointIds()->SetNumberOfIds(faceIds.size());
		//std::cout << "Face IDS: ";
		for (int i = 0; i < faceIds.size(); ++i) {
			face->GetPointIds()->SetId(i, faceIds[i]);
			//std::cout << faceIds[i] << " ";
		}
		//std::cout << std::endl;
		// Add the current face Ids to the container
		facesIds.push_back(faceIds);
		// Add the face to the poly data container
		polygons.push_back(face);
		//std::cout << "Number of Face points: " << faceIds.size() << std::endl;
	}

	// Add all faces to the Faces collection (Polygons array)
	for (auto &poly : polygons) {
		faces->InsertNextCell(poly);
	}
	//std::cout << "Number of Points: " << pointIds.size() << std::endl;
	//std::cout << "Number of faces: " << polygons.size() << std::endl;

	// Create a PolyData
	vtkSmartPointer<vtkPolyData> polygonPolyData = vtkSmartPointer<vtkPolyData>::New();
	polygonPolyData->SetPoints(points);
	polygonPolyData->SetPolys(faces);

	return polygonPolyData;
}

vtkSmartPointer<vtkPolyData> VTKContext::getVTKPolygonalData(Polyhedron_Ext &P) {
	Polyhedron_Ext::Facet_iterator facet_it;
	// Map Poly Point -> unique index
	std::map<Kernel_Ext::Point_3, int> point_to_index;
	std::map<Kernel_Ext::Point_3, int>::iterator it_point_to_index;
	std::vector<vtkIdType> pointIds;
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	std::vector<std::vector<vtkIdType>> facesIds;
	std::vector<vtkSmartPointer<vtkPolygon>> polygons;
	vtkSmartPointer<vtkCellArray> faces = vtkSmartPointer<vtkCellArray>::New();
	// It is the new index of the current vertex (chosen by me) 
	vtkIdType curr_v_index = 0;

	//std::cout << "Points from VTK context" << std::endl;
	// Extract the vertices, update the VTK context
	// Iterate every face (2-cell)

	for (Polyhedron_Ext::Facet_iterator f_it = P.facets_begin(); f_it != P.facets_end(); ++f_it) {
		Polyhedron_Ext::Halfedge_around_facet_circulator v_circ;
		vtkSmartPointer<vtkPolygon> face = vtkSmartPointer<vtkPolygon>::New();
		std::vector<vtkIdType> faceIds;

		v_circ = f_it->facet_begin();
		do {
			Kernel_Ext::Point_3 p = v_circ->vertex()->point();
			// Check if the point exists in the map
			it_point_to_index = point_to_index.find(p);
			// Does not exist: The vertex has not been traversed before
			if (it_point_to_index == point_to_index.end()) {
				// It will have the current index
				pointIds.push_back(curr_v_index);
				// Insert the vertex (point) and its index into the map
				point_to_index[p] = curr_v_index;
				// Insert the point in the VTKPoint structure
				points->InsertNextPoint(CGAL::to_double(p[0]), CGAL::to_double(p[1]), CGAL::to_double(p[2]));
				//std::cout << "New [" << curr_v_index << "]: " << p[0] << " " << p[1] << " " << p[2] << std::endl;
				// Add the current vertex ID to the face IDs
				faceIds.push_back(curr_v_index);
				++curr_v_index;
			}
			// The vertex exists: Push its ID in the face indexes
			else {
				// std::cout << "Old [" << it_point_to_index->second << "]: " << it_point_to_index->first[0] << " " << it_point_to_index->first[1] << " " << it_point_to_index->first[2] << std::endl;
				// Just push its index from the map to the face IDs
				faceIds.push_back(it_point_to_index->second);
			}
		} while (++v_circ != f_it->facet_begin());

		// Set the Point IDs for the current face
		face->GetPointIds()->SetNumberOfIds(faceIds.size());
		//std::cout << "Face IDS: ";
		for (int i = 0; i < faceIds.size(); ++i) {
			face->GetPointIds()->SetId(i, faceIds[i]);
			//std::cout << faceIds[i] << " ";
		}
		//std::cout << std::endl;
		// Add the current face Ids to the container
		facesIds.push_back(faceIds);
		// Add the face to the poly data container
		polygons.push_back(face);
		//std::cout << "Number of Face points: " << faceIds.size() << std::endl;
	}

	// Add all faces to the Faces collection (Polygons array)
	for (auto &poly : polygons) {
		faces->InsertNextCell(poly);
	}
	//std::cout << "Number of Points: " << pointIds.size() << std::endl;
	//std::cout << "Number of faces: " << polygons.size() << std::endl;

	// Create a PolyData
	vtkSmartPointer<vtkPolyData> polygonPolyData = vtkSmartPointer<vtkPolyData>::New();
	polygonPolyData->SetPoints(points);
	polygonPolyData->SetPolys(faces);

	return polygonPolyData;
}

vtkSmartPointer<vtkPolyData> VTKContext::getVTKPolygonalData(Nef_polyhedron_Ext &P) {
	// Map Poly Point -> unique index
	std::map<Kernel_Ext::Point_3, int> point_to_index;
	std::map<Kernel_Ext::Point_3, int>::iterator it_point_to_index;
	std::vector<vtkIdType> pointIds;
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	std::vector<std::vector<vtkIdType>> facesIds;
	std::vector<vtkSmartPointer<vtkPolygon>> polygons;
	vtkSmartPointer<vtkCellArray> faces = vtkSmartPointer<vtkCellArray>::New();
	// It is the new index of the current vertex (chosen by me) 
	vtkIdType curr_v_index = 0;

	//std::cout << "Points from VTK context" << std::endl;
	// Extract the vertices, update the VTK context
	// Iterate every face (2-cell)

	for (Nef_polyhedron_Ext::SFace_const_iterator f_it = P.sfaces_begin(); f_it != P.sfaces_end(); ++f_it) {
		vtkSmartPointer<vtkPolygon> face = vtkSmartPointer<vtkPolygon>::New();
		std::vector<vtkIdType> faceIds;

		for (Nef_polyhedron_Ext::SFace_cycle_const_iterator s_it = f_it->sface_cycles_begin(); s_it != f_it->sface_cycles_end(); s_it++) {
			if (s_it.is_shalfedge()) {
				Kernel_Ext::Point_3 p = Nef_polyhedron_Ext::SHalfedge_const_handle(s_it)->source()->point();
				//std::cout << p[0] << " " << p[1] << " " << p[2] << std::endl;
				// Check if the point exists in the map
				it_point_to_index = point_to_index.find(p);
				// Does not exist: The vertex has not been traversed before
				if (it_point_to_index == point_to_index.end()) {
					// It will have the current index
					pointIds.push_back(curr_v_index);
					// Insert the vertex (point) and its index into the map
					point_to_index[p] = curr_v_index;
					// Insert the point in the VTKPoint structure
					points->InsertNextPoint(CGAL::to_double(p[0]), CGAL::to_double(p[1]), CGAL::to_double(p[2]));
					//std::cout << "New [" << curr_v_index << "]: " << p[0] << " " << p[1] << " " << p[2] << std::endl;
					// Add the current vertex ID to the face IDs
					faceIds.push_back(curr_v_index);
					++curr_v_index;
				}
				// The vertex exists: Push its ID in the face indexes
				else {
					// std::cout << "Old [" << it_point_to_index->second << "]: " << it_point_to_index->first[0] << " " << it_point_to_index->first[1] << " " << it_point_to_index->first[2] << std::endl;
					// Just push its index from the map to the face IDs
					faceIds.push_back(it_point_to_index->second);
				}
			}
		}
		
		// Set the Point IDs for the current face
		face->GetPointIds()->SetNumberOfIds(faceIds.size());
		//std::cout << "Face IDS: ";
		for (int i = 0; i < faceIds.size(); ++i) {
			face->GetPointIds()->SetId(i, faceIds[i]);
			//std::cout << faceIds[i] << " ";
		}
		//std::cout << std::endl;
		// Add the current face Ids to the container
		facesIds.push_back(faceIds);
		// Add the face to the poly data container
		polygons.push_back(face);
		//std::cout << "Number of Face points: " << faceIds.size() << std::endl;
	}

	// Add all faces to the Faces collection (Polygons array)
	for (auto &poly : polygons) {
		faces->InsertNextCell(poly);
	}
	//std::cout << "Number of Points: " << pointIds.size() << std::endl;
	//std::cout << "Number of faces: " << polygons.size() << std::endl;

	// Create a PolyData
	vtkSmartPointer<vtkPolyData> polygonPolyData = vtkSmartPointer<vtkPolyData>::New();
	polygonPolyData->SetPoints(points);
	polygonPolyData->SetPolys(faces);

	return polygonPolyData;
}

void VTKContext::drawPolygonalData(vtkSmartPointer<vtkPolyData> polygonPolyData) {

	vtkSmartPointer<vtkNamedColors> colors = vtkSmartPointer<vtkNamedColors>::New();
	
	vtkSmartPointer<vtkDataSetMapper> mapper = vtkSmartPointer<vtkDataSetMapper>::New();
	mapper->SetInputData(polygonPolyData);

	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);
	actor->GetProperty()->SetColor(colors->GetColor3d("Silver").GetData());

	// Visualize the shape
	vtkSmartPointer<vtkRenderer> renderer =	vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderWindow> renderWindow =	vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->SetWindowName("Polygons");
	renderWindow->AddRenderer(renderer);
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =	vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renderWindowInteractor->SetRenderWindow(renderWindow);

	renderer->AddActor(actor);
	renderer->SetBackground(colors->GetColor3d("White").GetData());
//	renderer->ResetCamera();
//	renderer->GetActiveCamera()->Azimuth(30);
//	renderer->GetActiveCamera()->Elevation(30);
	renderWindow->Render();
	renderWindowInteractor->Start();

}

void VTKContext::writePolyDataInVTPFile(vtkSmartPointer<vtkPolyData> polyData, const char* filename) {
	vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer->SetFileName(filename);
	writer->SetInputData(polyData);
	writer->Write();
}

void VTKContext::writePolyDataFilterInVTPFile(vtkSmartPointer<vtkTransformPolyDataFilter> polyDataFilter, const char* filename) {
	vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer->SetFileName(filename);
	writer->SetInputConnection(polyDataFilter->GetOutputPort());
	writer->Write();
}

vtkSmartPointer<vtkPolyData> VTKContext::getPolyDataPoints(std::vector<Eigen::VectorXd> &points) {
	vtkSmartPointer<vtkPoints> polyPoints = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkCellArray> polyVertices = vtkSmartPointer<vtkCellArray>::New();
	std::vector<vtkIdType> pid;
	for (auto &p : points) {
		pid.push_back(polyPoints->InsertNextPoint(p.data()));
	}
	polyVertices->InsertNextCell(pid.size(), pid.data());
	polyData->SetPoints(polyPoints);
	polyData->SetVerts(polyVertices);
	return polyData;
}

std::vector<vtkSmartPointer<vtkTransformPolyDataFilter>> VTKContext::getPolyDataAxis(std::vector<double> cell_parameters) {
	double a_length = cell_parameters[0],
		b_length = cell_parameters[1],
		c_length = cell_parameters[2],
		alpha = cell_parameters[3],
		beta = cell_parameters[4],
		gamma = cell_parameters[5];
	double max_length = std::max({ a_length, b_length, c_length });
	Eigen::Matrix3d vec_comp = getTransformationMatrixFromFractionalToCartesian(cell_parameters);
	Eigen::Vector3d _a(vec_comp.col(0)), _b(vec_comp.col(1)), _c(vec_comp.col(2)), origin(0, 0, 0);
	_a.normalize();
	_b.normalize();
	_c.normalize();
	std::vector<vtkSmartPointer<vtkTransformPolyDataFilter>> arrows;

	vtkSmartPointer<vtkNamedColors> colors = vtkSmartPointer<vtkNamedColors>::New();
	// Create axis arrows.
	vtkSmartPointer<vtkArrowSource> arrowSourceA = vtkSmartPointer<vtkArrowSource>::New();
	vtkSmartPointer<vtkArrowSource> arrowSourceB = vtkSmartPointer<vtkArrowSource>::New();
	vtkSmartPointer<vtkArrowSource> arrowSourceC = vtkSmartPointer<vtkArrowSource>::New();
	arrowSourceA->SetShaftRadius(0.01);
	arrowSourceA->SetTipLength(0.1);
	arrowSourceA->SetTipRadius(0.03);
	arrowSourceA->SetTipResolution(20);
	arrowSourceB->SetShaftRadius(0.01);
	arrowSourceB->SetTipLength(0.1);
	arrowSourceB->SetTipRadius(0.03);
	arrowSourceB->SetTipResolution(20);
	arrowSourceC->SetShaftRadius(0.01);
	arrowSourceC->SetTipLength(0.1);
	arrowSourceC->SetTipRadius(0.03);
	arrowSourceC->SetTipResolution(20);
	arrowSourceA->Update();
	arrowSourceB->Update();
	arrowSourceC->Update();
	vtkSmartPointer<vtkMatrix4x4> matrixA = vtkSmartPointer<vtkMatrix4x4>::New();
	vtkSmartPointer<vtkMatrix4x4> matrixB = vtkSmartPointer<vtkMatrix4x4>::New();
	vtkSmartPointer<vtkMatrix4x4> matrixC = vtkSmartPointer<vtkMatrix4x4>::New();
	
	// Since A is fixed as X, transformation matrix for A will remain the Identity
	matrixA->Identity();
	matrixB->Identity();
	matrixC->Identity();

	for (auto i = 0; i < 3; i++) {
		matrixB->SetElement(i, 0, _b(i));
		matrixB->SetElement(i, 1, _a(i));
		matrixB->SetElement(i, 2, _c(i));
	}
	
	// Angles between C and X, Y, Z
	//std::cout << "ALPHA: " << alpha << " --- BETA: " << beta << " --- GAMMA: " << gamma << std::endl;

	for (auto i = 0; i < 3; i++) {
		matrixC->SetElement(i, 0, _c(i));
		matrixC->SetElement(i, 1, _b(i));
		matrixC->SetElement(i, 2, _a(i));
	}	
	
	// Apply the transforms
	vtkSmartPointer<vtkTransform> transformA = vtkSmartPointer<vtkTransform>::New();
	vtkSmartPointer<vtkTransform> transformB = vtkSmartPointer<vtkTransform>::New();
	vtkSmartPointer<vtkTransform> transformC = vtkSmartPointer<vtkTransform>::New();
	transformA->Translate(origin.data());
	transformA->Concatenate(matrixA);
	transformA->Scale(a_length, a_length, a_length);
	transformA->Update();
	transformB->Translate(origin.data());
	transformB->Concatenate(matrixB);
	transformB->Scale(b_length, b_length, b_length);
	transformB->Update();
	transformC->Translate(origin.data());
	transformC->Concatenate(matrixC);
	transformC->Scale(c_length, c_length, c_length);
	transformC->Update();
	
	// Transform the polydata
	vtkSmartPointer<vtkTransformPolyDataFilter> transformPDA = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
	vtkSmartPointer<vtkTransformPolyDataFilter> transformPDB = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
	vtkSmartPointer<vtkTransformPolyDataFilter> transformPDC = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
	transformPDA->SetTransform(transformA);
	transformPDA->SetInputConnection(arrowSourceA->GetOutputPort());
	transformPDA->Update();
	transformPDB->SetTransform(transformB);
	transformPDB->SetInputConnection(arrowSourceB->GetOutputPort());
	transformPDB->Update();
	transformPDC->SetTransform(transformC);
	transformPDC->SetInputConnection(arrowSourceC->GetOutputPort());
	transformPDC->Update();
	arrows.push_back(transformPDA);
	arrows.push_back(transformPDB);
	arrows.push_back(transformPDC);
	return arrows;
}
