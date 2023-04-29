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
#include <voronoicell.h>

/*****************************\
|******* VORONOICELL *********|
\*****************************/

std::vector<Point_CM> VoronoiCell::getPoints()
{
	return this->points;
}

void VoronoiCell::setPoints(std::vector<Point_CM> cell_points) {
	this->points.assign(cell_points.begin(), cell_points.end());
}

void VoronoiCell::setPoints(std::vector<std::vector<double>> cell_points) {
	if (!this->points.empty()) this->points.clear();
	for (int i = 0; i < cell_points.size(); ++i) {
		Point_CM p(cell_points[i][0], cell_points[i][1], cell_points[i][2]);
		this->points.push_back(p);
	}
}

std::vector<std::vector<int>> VoronoiCell::getFaces()
{
	return this->faces;
}

void VoronoiCell::setFaces(std::vector<std::vector<int>> cell_faces) {
	this->faces.assign(cell_faces.begin(), cell_faces.end());
}

void VoronoiCell::setSite(Eigen::VectorXd& voronoi_site) {
	this->site = voronoi_site;
}

void VoronoiCell::setSite(std::vector<double> voronoi_site) {
	this->site = Eigen::VectorXd(3);
	this->site << voronoi_site[0], voronoi_site[1], voronoi_site[2];
}

void VoronoiCell::setName(std::string& s) {
	this->name.assign(s.begin(), s.end());
}

void VoronoiCell::updatePolyhedron() {
	std::vector<Kernel::Point_3> ps_3;
	for (auto& p : this->points) {
		Kernel::Point_3 p3(p.x(), p.y(), p.z());
		ps_3.push_back(p3);
	}
	CGAL::convex_hull_3(ps_3.begin(), ps_3.end(), this->polyhedron);
}

bool VoronoiCell::hasCorrectGeometry() {
	/* If points are less then faces, the geometry will not be correct */
	if (this->getPoints().size() < this->getFaces().size()) return false;
	/* If empty vectors, geometry incorrect */
	if (this->getPoints().empty() || this->getFaces().empty()) return false;
	/* If a point is counted 2 times or less in total -> Incorrect geometry */
	for (auto& count = this->point_counter.begin(); count != this->point_counter.end(); ++count) {
		if (count->second < 3) return false;
	}
	
	return true;
}

struct classcomp {
	bool operator() (const int& lhs, const int& rhs) const
	{
		return lhs > rhs;
	}
};

void VoronoiCell::fixGeometry() {
	this->updatePointCounter();
	std::set<int, classcomp> toremove_indeces;
	// Delete incorrect points indeces in faces
	// Each point should appear in at least 3 faces
	for (int i = 0; i < this->faces.size(); ++i) {
		for (int j = 0; j < this->faces[i].size(); ++j) {
			if (this->point_counter[this->faces[i][j]] < 3) {
				toremove_indeces.insert(this->faces[i][j]);
			}
		}
		for (auto& toremove_p : toremove_indeces) {
			std::remove(this->faces[i].begin(), this->faces[i].end(), toremove_p);
		}
	}	
	// Delete incorrect point in point list
	std::set<int, classcomp>::iterator set_it;
	for (auto& toremove_p : toremove_indeces) {
		this->points.erase(this->points.begin() + toremove_p);
	}
	// Update indeces in faces
	for (int i = 0; i < this->faces.size(); ++i) {
		for (int j = 0; j < this->faces[i].size(); ++j) {
			for (auto& toremove_p : toremove_indeces) {
				// There should not be = number since they were deleted before
				if (toremove_p < this->faces[i][j]) {
					--this->faces[i][j];
				}
			}
		}
	}
	
	this->updatePointCounter();
}

void VoronoiCell::updatePointCounter() {
	std::map<int, int>::iterator map_it;
	this->point_counter.clear();
	/* Store point count in map */
	for (auto& face = this->faces.begin(); face != this->faces.end() ; ++face) {
		for (int p_index = 0; p_index < face->size(); ++p_index) {
			map_it = this->point_counter.find(face->at(p_index));
			if (map_it != this->point_counter.end()) {
				++this->point_counter[map_it->first];
			}
			else {
				this->point_counter.insert(std::pair<int, int>(face->at(p_index), 1));
			}
		}
	}
}

double VoronoiCell::getVolume()
{
	Triangulation_basic T;
	double volume = 0;
	std::vector<Kernel::Point_3> ps_3;

	for (auto& p : this->points) {
		Kernel::Point_3 p3(p.x(), p.y(), p.z());
		ps_3.push_back(p3);
	}
	T.insert(ps_3.begin(), ps_3.end());

	Tb_ACell_iterator t_cell;
	for (t_cell = T.all_cells_begin(); t_cell != T.all_cells_end(); ++t_cell) {
		if (!T.is_infinite(t_cell)) {
			// compute the volume of the tetrahedron
			volume += CGAL::to_double(CGAL::abs(CGAL::volume(t_cell->vertex(0)->point(), t_cell->vertex(1)->point(), t_cell->vertex(2)->point(), t_cell->vertex(3)->point())));
		}
	}
	return volume;
}

Point_CM getRotatedPoint_cm(Kernel_simple::Vector_3 p, Kernel_simple::Vector_3 a, double angle_rad) {
	double c = cos(angle_rad), s = sin(angle_rad);
	Kernel_simple::Vector_3 v = p * c + (1 - c) * (CGAL::to_double(a * p) * a) + s * (CGAL::cross_product(a, p));

	return Point_CM(CGAL::to_double(v.x()), CGAL::to_double(v.y()), CGAL::to_double(v.z()));
}

std::vector<Point_CM> VoronoiCell::getRotatedPointsAboutAngleAndAxis(std::vector<double> rotation_sample)
{
	double theta = rotation_sample[0],
		mu = rotation_sample[1],
		z = rotation_sample[2];
	std::vector<Point_CM> rotated_points;
	Kernel_simple::Vector_3 rotation_axis(cos(mu) * sqrt(1 - pow(z, 2)), sin(mu) * sqrt(1 - pow(z, 2)), z);
	for (int i = 0; i < this->points.size(); ++i) {
		Kernel_simple::Vector_3 v(CGAL::to_double(this->points[i].x()), CGAL::to_double(this->points[i].y()), CGAL::to_double(this->points[i].z()));
		rotated_points.push_back(getRotatedPoint_cm(v, rotation_axis, theta));
	}
	return rotated_points;
}

std::vector<Point_CM> VoronoiCell::getRotatedPointsAboutAngleAndAxis(std::vector<double> rotationaxis_sample, double anglerad_sample) {
	double theta = anglerad_sample;
	std::vector<Point_CM> rotated_points;
	Kernel_simple::Vector_3 rotation_axis(rotationaxis_sample[0], rotationaxis_sample[1], rotationaxis_sample[2]);
	for (int i = 0; i < this->points.size(); ++i) {
		Kernel_simple::Vector_3 v(CGAL::to_double(this->points[i].x()), CGAL::to_double(this->points[i].y()), CGAL::to_double(this->points[i].z()));
		rotated_points.push_back(getRotatedPoint_cm(v, rotation_axis, theta));
	}
	return rotated_points;
}

VoronoiCell VoronoiCell::getRotatedVoronoiDomain(std::vector<double> rotationaxis_sample, double anglerad_sample) {
	VoronoiCell v_rotated(*this);
	double theta = anglerad_sample;
	std::vector<Point_CM> rotated_points;
	Kernel_simple::Vector_3 rotation_axis(rotationaxis_sample[0], rotationaxis_sample[1], rotationaxis_sample[2]);
	for (int i = 0; i < this->points.size(); ++i) {
		Kernel_simple::Vector_3 v(CGAL::to_double(this->points[i].x()), CGAL::to_double(this->points[i].y()), CGAL::to_double(this->points[i].z()));
		rotated_points.push_back(getRotatedPoint_cm(v, rotation_axis, theta));
	}
	v_rotated.setPoints(rotated_points);

	return v_rotated;
}

bool VoronoiCell::isPointInside(Eigen::VectorXd& query) {
	AABB_Tree polyhedron_tree(this->polyhedron.facets_begin(), this->polyhedron.facets_end(), this->polyhedron);
	polyhedron_tree.accelerate_distance_queries();
	Point_inside_Poly3 point_inside_test(polyhedron_tree);
	std::vector<int> pointindeces_inside;

	Kernel::Point_3 query_p(query[0], query[1], query[2]);
	if (point_inside_test(query_p) == CGAL::ON_BOUNDED_SIDE) return true;
	else return false;
}

std::vector<int> VoronoiCell::getPointIndecesInside(std::vector<Eigen::VectorXd> &queries) {
	AABB_Tree polyhedron_tree(this->polyhedron.facets_begin(), this->polyhedron.facets_end(), this->polyhedron);
	polyhedron_tree.accelerate_distance_queries();
	Point_inside_Poly3 point_inside_test(polyhedron_tree);
	std::vector<int> pointindeces_inside;

	for (int i = 0; i < queries.size(); ++i) {
		Kernel::Point_3 query_p(queries[i][0], queries[i][1], queries[i][2]);
		if (point_inside_test(query_p) == CGAL::ON_BOUNDED_SIDE) pointindeces_inside.push_back(i);
	}
	return pointindeces_inside;
}

VoronoiCell VoronoiCell::getShiftedToOriginVoronoiDomain() {
	assert(! this->points.empty());
	VoronoiCell shifted_domain;
	std::vector<Point_CM> shifted_points;
	for (int i = 0; i < this->points.size(); ++i) {
		Point_CM shifted_point(CGAL::to_double(this->points[i][0]) - this->site[0], CGAL::to_double(this->points[i][1]) - this->site[1], CGAL::to_double(this->points[i][2]) - this->site[2]);
		shifted_points.push_back(shifted_point);
	}
	shifted_domain.setPoints(shifted_points);
	shifted_domain.setFaces(this->faces);
	//Eigen::VectorXd newsite(3); newsite << 0, 0, 0;
	Eigen::VectorXd newsite(this->site - this->site);
	shifted_domain.setSite(newsite);
	shifted_domain.updatePolyhedron();

	return shifted_domain;
}

VoronoiCell VoronoiCell::getVoronoiDomainShiftedBy(Eigen::VectorXd offset_vector) {
	assert(!this->points.empty());
	VoronoiCell shifted_domain;
	std::vector<Point_CM> shifted_points;
	for (int i = 0; i < this->points.size(); ++i) {
		Point_CM shifted_point(CGAL::to_double(this->points[i][0]) + offset_vector[0], CGAL::to_double(this->points[i][1]) + offset_vector[1], CGAL::to_double(this->points[i][2]) + offset_vector[2]);
		shifted_points.push_back(shifted_point);
	}
	shifted_domain.setPoints(shifted_points);
	shifted_domain.setFaces(this->faces);
	//Eigen::VectorXd newsite(3); newsite << 0, 0, 0;
	Eigen::VectorXd newsite(this->site + offset_vector);
	shifted_domain.setSite(newsite);
	shifted_domain.updatePolyhedron();

	return shifted_domain;
}

void VoronoiCell::writeToVTPFileFormat(const char* filename)
{
	vtkSmartPointer<vtkPolyData> polygons = VTKContext::getVTKPolygonalData(this->points, this->faces);
	VTKContext::writePolyDataInVTPFile(polygons, filename);
}

void VoronoiCell::writePointsToCSVFileFormat(const char* filename) {
	std::ofstream out_points;
	out_points.open(filename);
	out_points << "X,Y,Z" << std::endl;
	for (Point_CM &p : this->getPoints()) {
		for (int i = 0; i < 3; ++i) {
			if (i == 2) {
				out_points << p[i];
			}
			else {
				out_points << p[i] << ",";
			}
		}
		out_points << std::endl;
	}
	out_points.close();
}

void VoronoiCell::writeToOFFFileFormat(const char* filename) {
	std::ofstream out_points;
	out_points.open(filename);
	out_points << this->points.size() << "\t" << this->faces.size() << std::endl;
	for (Point_CM&p : this->points) {
		for (int i = 0; i < 3; ++i) {
			if (i == 2) {
				out_points << p[i];
			}
			else {
				out_points << p[i] << ",";
			}
		}
		out_points << std::endl;
	}
	for (std::vector<int> &f : this->faces) {
		out_points << f.size() << "\t";
		for (int i = 0; i < f.size(); ++i) {
			if (i == f.size()-1) {
				out_points << f[i];
			}
			else {
				out_points << f[i] << ",";
			}
		}
		out_points << std::endl;
	}

	out_points.close();
}

/************************************\
|******* VORONOICELLMETRICS *********|
\************************************/

void VoronoiCellMetrics::setThreads(int n_threads) {
	this->threads = n_threads;
}

void VoronoiCellMetrics::addVoronoiDomain(VoronoiCell& voronoi_domain, std::string name) {
	this->voronoi_vector.push_back(voronoi_domain);
	this->names.push_back(name);
}

void VoronoiCellMetrics::setSamplingMode(SamplingMode s_mode) {
	this->sampling_mode = s_mode;
}

/********************\
|***** UNIFORM ******|
\********************/
void VoronoiCellMetrics::getOffsetOverRotation(VoronoiCell v1, VoronoiCell v2, std::vector<double> rotation_sample, VoronoiMetricResult &vormetric_result)
{
	std::vector<Point_CM> rotated_points = v1.getRotatedPointsAboutAngleAndAxis(rotation_sample), points_lcc2 = v2.getPoints();
	std::vector<std::vector<int>> faces_lcc2 = v2.getFaces();
	//Maximum distance to cover all vertices
	double offset_overotation = -1;

	Kernel_simple::Point_3 origin(0.f, 0.f, 0.f);
	for (Point_CM& p : rotated_points) {
		const Kernel_simple::Point_3 p1(CGAL::to_double(p.x()), CGAL::to_double(p.y()), CGAL::to_double(p.z()));
		const Kernel_simple::Segment_3 segment = Kernel_simple::Segment_3(origin, p1);

		for (std::vector<int>& face : faces_lcc2) {
			//std::cout << "\tPlane Points: " << face[0] << ", " << face[1] << ", " << face[2] << std::endl;
			// A face is defined at least by 3 points
			Kernel_simple::Point_3 f1(CGAL::to_double(points_lcc2[face[0]].x()), CGAL::to_double(points_lcc2[face[0]].y()), CGAL::to_double(points_lcc2[face[0]].z()));
			Kernel_simple::Point_3 f2(CGAL::to_double(points_lcc2[face[1]].x()), CGAL::to_double(points_lcc2[face[1]].y()), CGAL::to_double(points_lcc2[face[1]].z()));
			Kernel_simple::Point_3 f3(CGAL::to_double(points_lcc2[face[2]].x()), CGAL::to_double(points_lcc2[face[2]].y()), CGAL::to_double(points_lcc2[face[2]].z()));
			const Kernel_simple::Plane_3 plane = Kernel_simple::Plane_3(f1, f2, f3);
			CGAL::cpp11::result_of<Kernel_simple::Intersect_3(Kernel_simple::Segment_3, Kernel_simple::Plane_3)>::type result = CGAL::intersection(segment, plane);
			if (result) {
				// Vertex p is outside the polyhedron

				double curr_distancetopoly_vertices = sqrt(CGAL::to_double(CGAL::squared_distance(p1, plane)));

				if (offset_overotation == -1) {
					offset_overotation = curr_distancetopoly_vertices;
				}
				else {
					if (curr_distancetopoly_vertices > offset_overotation) {
						offset_overotation = curr_distancetopoly_vertices;
					}
				}
			}
			else {
				// Vertex p is inside the polyhedron
				//std::cout << "\t\tNo Intersection..." << std::endl;
			}
		}
	}
	// Case: all v1 vertices inside v2 
	if (offset_overotation == -1) {
		offset_overotation = 0;
	}

	vormetric_result.distance_result = offset_overotation;
	vormetric_result.rotation_sample = rotation_sample;
}

void VoronoiCellMetrics::getOffsetOverRotationsSet(VoronoiCell v1, VoronoiCell v2, int thread_id, int start, int end, VoronoiMetricResult& vormetricresult_vector)
{
	// Add a check on current thread result and update vormetricresult_vector[i]
	for (int i = start; i < end; ++i) {
		VoronoiMetricResult offset_overRotationsSet;
		this->getOffsetOverRotation(v1, v2, this->rotation_samples[i], offset_overRotationsSet);
		if (i == start) vormetricresult_vector = offset_overRotationsSet;
		else {
			if (offset_overRotationsSet.distance_result < vormetricresult_vector.distance_result) {
				vormetricresult_vector.distance_result = offset_overRotationsSet.distance_result;
				vormetricresult_vector.rotation_sample = offset_overRotationsSet.rotation_sample;
			}
		}
	}
}

VoronoiMetricResult VoronoiCellMetrics::getOffset(VoronoiCell v1, VoronoiCell v2)
{
	VoronoiMetricResult offset_overAllRotations;
	std::vector<VoronoiMetricResult> offset_results;

	for (int i = 0; i < rotation_samples.size(); ++i) {
		VoronoiMetricResult voronoi_result;
		offset_results.push_back(voronoi_result);
	}
	// Get the offset over a rotation of every rotation sample
	for (int i = 0; i < rotation_samples.size(); ++i) {
		this->getOffsetOverRotation(v1, v2, this->rotation_samples[i], offset_results[i]);
	}

	// Get the minimum over all rotations: Offset(v1,v2)
	int min_offset_index = 0;
	for (int i = 1; i < offset_results.size(); ++i) {
		if (offset_results[i].distance_result < offset_results[min_offset_index].distance_result) {
			min_offset_index = i;
		}
	}

	offset_overAllRotations = offset_results[min_offset_index];
	return offset_overAllRotations;
}

VoronoiMetricResult VoronoiCellMetrics::getOffset_Parallel_std(VoronoiCell v1, VoronoiCell v2)
{
	VoronoiMetricResult offset_overAllRotations;
	std::vector<VoronoiMetricResult> offset_results(this->threads);
	std::vector<std::thread> pool(this->threads);
	int start=0, end=0, samples_per_thread = this->rotation_samples.size() / this->threads;
	//std::cout << "Rotation samples: " << this->rotation_samples.size() << " - Samples per thread: " << samples_per_thread << std::endl;

	for (int i = 0; i < this->threads; ++i) {
		start = (i * samples_per_thread);
		if (i == this->threads - 1) end = this->rotation_samples.size();
		else end = (start + samples_per_thread);
		//std::cout << "Start: " << start << " - End: " << end << std::endl;

		// Get the offset over a set of rotations
		auto f = std::bind(&VoronoiCellMetrics::getOffsetOverRotationsSet, this, v1, v2, i, start, end, std::ref(offset_results[i]));
		pool[i] = std::thread(f);
	}
	for (int i = 0; i < this->threads; ++i) {
		pool[i].join();
	}
	// Get the minimum over all rotations: Offset(v1,v2)
	int min_offset_index = 0;
	//std::cout << offset_results[0].distance_result << std::endl;
	for (int i = 1; i < offset_results.size(); ++i) {
		//std::cout << offset_results[i].distance_result << std::endl;
		if (offset_results[i].distance_result < offset_results[min_offset_index].distance_result) {
			min_offset_index = i;
		}
	}
	offset_overAllRotations = offset_results[min_offset_index];
	//std::cout << "Minimum: " << offset_overAllRotations.distance_result << std::endl;
	return offset_overAllRotations;
}

void VoronoiCellMetrics::getScaleOverRotation(VoronoiCell v1, VoronoiCell v2, std::vector<double> rotation_sample, VoronoiMetricResult &vormetric_result)
{
	std::vector<Point_CM> rotated_points = v1.getRotatedPointsAboutAngleAndAxis(rotation_sample), points_lcc2 = v2.getPoints();
	std::vector<std::vector<int>> faces_lcc2 = v2.getFaces();
	const Kernel_simple::Point_3 origin(0.f, 0.f, 0.f);
	//Maximum distance to cover all vertices
	double scale_overotation = -1;

	for (auto& p : rotated_points) {
		const Kernel_simple::Point_3 p1(CGAL::to_double(p.x()), CGAL::to_double(p.y()), CGAL::to_double(p.z()));
		Kernel_simple::Vector_3 min_distance_o_pint_vector;
		double distance_o_p = sqrt(CGAL::to_double(CGAL::squared_distance(origin, p1))), min_distance_o_pint = -1, curr_distance_o_pint;
		const Kernel_simple::Ray_3 ray = Kernel_simple::Ray_3(origin, p1);
		bool found = false;
		//std::cout << "Vertex: " << p_3 << std::endl;
		// Pick the intersection of the Ray_p with the closest face ( min distance(O, F_face intesected Ray_p) )

		for (auto& face : faces_lcc2) {

			//std::cout << "\tPlane Points: " << points_lcc2[face[0]] << ", " << points_lcc2[face[1]] << ", " << points_lcc2[face[2]] << std::endl;
			Kernel_simple::Point_3 f1(CGAL::to_double(points_lcc2[face[0]][0]), CGAL::to_double(points_lcc2[face[0]][1]), CGAL::to_double(points_lcc2[face[0]][2]));
			Kernel_simple::Point_3 f2(CGAL::to_double(points_lcc2[face[1]][0]), CGAL::to_double(points_lcc2[face[1]][1]), CGAL::to_double(points_lcc2[face[1]][2]));
			Kernel_simple::Point_3 f3(CGAL::to_double(points_lcc2[face[2]][0]), CGAL::to_double(points_lcc2[face[2]][1]), CGAL::to_double(points_lcc2[face[2]][2]));
			const Kernel_simple::Plane_3 plane = Kernel_simple::Plane_3(f1, f2, f3);
			CGAL::cpp11::result_of<Kernel_simple::Intersect_3(Kernel_simple::Ray_3, Kernel_simple::Plane_3)>::type result = CGAL::intersection(ray, plane);
			if (result) {
				// Vertex p is outside v2				
				// Get the vector with the minimum distance from origin to face intersection (closest intersection to origin)
				if (const Kernel_simple::Point_3 * p_int = boost::get<Kernel_simple::Point_3>(&*result)) {
					found = true;
					const Kernel_simple::Vector_3 v_int(CGAL::to_double(p_int->x()), CGAL::to_double(p_int->y()), CGAL::to_double(p_int->z()));
					curr_distance_o_pint = sqrt(CGAL::to_double(v_int.squared_length()));
					if (min_distance_o_pint == -1) {
						min_distance_o_pint = curr_distance_o_pint;
						min_distance_o_pint_vector = v_int;
					}
					else {
						if (curr_distance_o_pint < min_distance_o_pint) {
							min_distance_o_pint = curr_distance_o_pint;
							min_distance_o_pint_vector = v_int;
						}
					}
				}				
			}
			else {
				// Vertex p is inside the polyhedron
				//std::cout << "\t\tNo Intersection..." << std::endl;
			}
		}
		// If found a face intersection: Pick max (Scale(P1,P2))
		if (found) {
			double curr_scale = distance_o_p / min_distance_o_pint;
			if (scale_overotation == -1) {
				scale_overotation = curr_scale;
			}
			else {
				// Pick the maximum to cover all vertices
				if (curr_scale > scale_overotation) {
					scale_overotation = curr_scale;
				}
			}
			found = false;
		}
		else {
			// Never entered: Intersection always found
		}
	}

	vormetric_result.distance_result = scale_overotation;
	vormetric_result.rotation_sample = rotation_sample;
}

void VoronoiCellMetrics::getScaleOverRotationsSet(VoronoiCell v1, VoronoiCell v2, int thread_id, int start, int end, VoronoiMetricResult& vormetricresult_vector)
{
	// Add a check on current thread result and update vormetricresult_vector[i]
	for (int i = start; i < end; ++i) {
		VoronoiMetricResult scale_overRotationsSet;
		this->getScaleOverRotation(v1, v2, this->rotation_samples[i], scale_overRotationsSet);
		if (i == start) vormetricresult_vector = scale_overRotationsSet;
		else {
			if (scale_overRotationsSet.distance_result < vormetricresult_vector.distance_result) {
				vormetricresult_vector.distance_result = scale_overRotationsSet.distance_result;
				vormetricresult_vector.rotation_sample = scale_overRotationsSet.rotation_sample;
			}
		}
	}
}

VoronoiMetricResult VoronoiCellMetrics::getScale(VoronoiCell v1, VoronoiCell v2)
{
	VoronoiMetricResult scale_overAllRotations;
	std::vector<VoronoiMetricResult> scale_results;

	for (int i = 0; i < rotation_samples.size(); ++i) {
		VoronoiMetricResult voronoi_result;
		scale_results.push_back(voronoi_result);
	}
	// Get the scale over a rotation of every rotation sample
	for (int i = 0; i < rotation_samples.size(); ++i) {
		this->getScaleOverRotation(v1, v2, this->rotation_samples[i], scale_results[i]);
	}

	// Get the minimum over all rotations: Scale(v1,v2)
	int min_scale_index = 0;
	for (int i = 1; i < scale_results.size(); ++i) {
		if (scale_results[i].distance_result < scale_results[min_scale_index].distance_result) {
			min_scale_index = i;
		}
	}

	scale_overAllRotations = scale_results[min_scale_index];
	return scale_overAllRotations;
}

VoronoiMetricResult VoronoiCellMetrics::getScale_Parallel_std(VoronoiCell v1, VoronoiCell v2)
{
	VoronoiMetricResult scale_overAllRotations;
	std::vector<VoronoiMetricResult> scale_results(this->threads);

	std::vector<std::thread> pool(this->threads);
	int start = 0, end = 0, samples_per_thread = this->rotation_samples.size() / this->threads;
	//std::cout << "Rotation samples: " << this->rotation_samples.size() << " - Samples per thread: " << samples_per_thread << std::endl;

	for (int i = 0; i < this->threads; ++i) {
		start = (i * samples_per_thread);
		if (i == this->threads - 1) end = this->rotation_samples.size();
		else end = (start + samples_per_thread);
		//std::cout << "Start: " << start << " - End: " << end << std::endl;

		// Get the offset over a set of rotations
		auto f = std::bind(&VoronoiCellMetrics::getScaleOverRotationsSet, this, v1, v2, i, start, end, std::ref(scale_results[i]));
		pool[i] = std::thread(f);
	}
	for (int i = 0; i < this->threads; ++i) {
		pool[i].join();
	}
	// Get the minimum over all rotations: Scale(v1,v2)
	int min_scale_index = 0;
	//std::cout << scale_results[0].distance_result << std::endl;
	for (int i = 1; i < scale_results.size(); ++i) {
		//std::cout << scale_results[i].distance_result << std::endl;
		if (scale_results[i].distance_result < scale_results[min_scale_index].distance_result) {
			min_scale_index = i;
		}
	}

	scale_overAllRotations = scale_results[min_scale_index];
	//std::cout << "Minimum: " << scale_overAllRotations.distance_result << std::endl;
	return scale_overAllRotations;
}

void VoronoiCellMetrics::generateUniformRotationSamples(int intervals) {
	this->rotation_samples.clear();
	double angle_samples = 2 * M_PI * intervals;
	// Generate rotation samples
	for (double theta = 0; theta < 2 * M_PI; theta += 1.f/intervals) {
		this->rotation_samples.push_back(std::vector<double>{theta, 0, 1});
		for (double mu = 0; mu < 2 * M_PI; mu += 1.f / intervals) {
			for (double z = 1.f / (2 * intervals); z < 1; z += (1.f / intervals)) {
				std::vector<double> rotationSample(3);
				rotationSample[0] = theta;
				rotationSample[1] = mu;
				rotationSample[2] = z;
				this->rotation_samples.push_back(rotationSample);
			}
		}
	}
	//std::cout << "Number of rotation samples: " << this->rotation_samples.size() << std::endl;
}

VoronoiMetricResult VoronoiCellMetrics::getExtendedHausdorffDistance(VoronoiCell v1, VoronoiCell v2)
{
	VoronoiMetricResult r_v1, r_v2, r_final;

	r_v1 = this->getOffset_Parallel_std(v1, v2);
	r_final.distance_vector.push_back(r_v1.distance_result);
	r_final.sample_vector.push_back(std::vector<double>(r_v1.rotation_sample));

	r_v2 = this->getOffset_Parallel_std(v2, v1);
	r_final.distance_vector.push_back(r_v2.distance_result);
	r_final.sample_vector.push_back(std::vector<double>(r_v2.rotation_sample));

	// Get the MAX between two offsets
	if (r_v1.distance_result > r_v2.distance_result) {
		r_final.distance_result = r_v1.distance_result;
		r_final.rotation_sample = r_v1.rotation_sample;
	}
	else {
		r_final.distance_result = r_v2.distance_result;
		r_final.rotation_sample = r_v2.rotation_sample;
	}
	return r_final;
}

VoronoiMetricResult VoronoiCellMetrics::getScaleInvariantDistance(VoronoiCell v1, VoronoiCell v2)
{
	VoronoiMetricResult r_v1, r_v2, r_final;

	r_v1 = this->getScale_Parallel_std(v1, v2);
	r_final.distance_vector.push_back(r_v1.distance_result);
	r_final.sample_vector.push_back(std::vector<double>(r_v1.rotation_sample));
	
	r_v2 = this->getScale_Parallel_std(v2, v1);
	r_final.distance_vector.push_back(r_v2.distance_result);
	r_final.sample_vector.push_back(std::vector<double>(r_v2.rotation_sample));

	// Get the logarithm (ln) of the MAX between two scales
	if (r_v1.distance_result > r_v2.distance_result) {
		r_final.distance_result = std::log(r_v1.distance_result);
		r_final.rotation_sample = r_v1.rotation_sample;
	}
	else {
		r_final.distance_result = std::log(r_v2.distance_result);
		r_final.rotation_sample = r_v2.rotation_sample;
	}
	return r_final;
}

/*************************\
|*************************|
\*************************/

void VoronoiCellMetrics::updateExtendedHausdorffDistanceMatrix() {
	int voronoi_count = this->voronoi_vector.size(), pair_count=0;
	this->extHausdorffDistMatrix.resize(voronoi_count, voronoi_count);
	this->extHausdorffDistMatrix.setZero();
	this->optimalRotations.clear();

	for (int i = 0; i < voronoi_count; ++i) {
		this->optimalRotations.push_back(std::vector<VoronoiMetricResult>(voronoi_count));
	}

	for (int i = 0; i < voronoi_count - 1; ++i) {
		for (int j = i + 1; j < voronoi_count; ++j) {
			if (i != j) {
				// Extended Hausdorff distance
				VoronoiMetricResult result_EHD;
				switch (this->sampling_mode) {
					case SamplingMode::UNIFORM : {
						result_EHD = this->getExtendedHausdorffDistance(this->voronoi_vector[i], this->voronoi_vector[j]);
						break;
					}
					default: break;
				}

				this->extHausdorffDistMatrix(i, j) = result_EHD.distance_result;
				this->extHausdorffDistMatrix(j, i) = result_EHD.distance_result;
				this->optimalRotations[i][j] = result_EHD;
				this->optimalRotations[j][i] = result_EHD;
			}
			std::cout << '\r' << "Pair Count: " << ++pair_count;
		}
	}
	std::cout << std::endl;
}

void VoronoiCellMetrics::updateScaleInvariantDistanceMatrix() {
	int voronoi_count = this->voronoi_vector.size(), pair_count=0;
	this->scaleInvDistMatrix.resize(voronoi_count, voronoi_count);
	this->scaleInvDistMatrix.setZero();
	this->optimalRotations.clear();

	for (int i = 0; i < voronoi_count; ++i) {
		this->optimalRotations.push_back(std::vector<VoronoiMetricResult>(voronoi_count));
	}

	for (int i = 0; i < voronoi_count - 1; ++i) {
		for (int j = i + 1; j < voronoi_count; ++j) {
			if (i != j) {
				// Scale-Invariant distance
				VoronoiMetricResult result_SID;
				switch (this->sampling_mode) {
				case SamplingMode::UNIFORM: {
					result_SID = this->getScaleInvariantDistance(this->voronoi_vector[i], this->voronoi_vector[j]);
					break;
				}
				default: break;
				}
				
				this->scaleInvDistMatrix(i, j) = result_SID.distance_result;
				this->scaleInvDistMatrix(j, i) = result_SID.distance_result;
				this->optimalRotations[i][j] = result_SID;
				this->optimalRotations[j][i] = result_SID;
			}
			std::cout << '\r' <<  "Pair Count: "<< ++pair_count;
		}
	}
	std::cout << std::endl;

}

void VoronoiCellMetrics::writeExtendedHausdorffDistanceMatrixToCSV(std::string filename) {
	int voronoi_count = this->voronoi_vector.size();
	std::ofstream out_matrix_DH;
	out_matrix_DH.open(filename + ".csv");
	out_matrix_DH << "ID";
	for (int i = 0; i < voronoi_count; ++i) {
		out_matrix_DH << "," << this->voronoi_vector[i].name;
	}
	out_matrix_DH << '\n';
	for (int i = 0; i < voronoi_count; ++i) {
		out_matrix_DH << this->voronoi_vector[i].name;
		for (int j = 0; j < voronoi_count; ++j) {
			out_matrix_DH << ',' << this->extHausdorffDistMatrix(i, j);
		}
		out_matrix_DH << '\n';
	}
	out_matrix_DH.close();
}

void VoronoiCellMetrics::writeScaleInvariantDistanceMatrixToCSV(std::string filename) {
	int voronoi_count = this->voronoi_vector.size();
	std::ofstream out_matrix_DS;
	out_matrix_DS.open(filename + ".csv");
	out_matrix_DS << "ID";
	for (int i = 0; i < voronoi_count; ++i) {
		out_matrix_DS << "," << this->names[i];
	}
	out_matrix_DS << '\n';
	for (int i = 0; i < voronoi_count; ++i) {
		out_matrix_DS << this->names[i];
		for (int j = 0; j < voronoi_count; ++j) {
			out_matrix_DS << ',' << this->scaleInvDistMatrix(i, j);
		}
		out_matrix_DS << '\n';
	}
	out_matrix_DS.close();

}

void VoronoiCellMetrics::writeOptimalRotationsToVTP(std::string path) {
	int voronoi_vector_size = this->voronoi_vector.size();
	std::vector<VoronoiCell> vcell_vector(2);
	std::string final_name{};
	for (int row = 0; row < voronoi_vector_size - 1; ++row) {
		for (int col = row + 1; col < voronoi_vector_size; ++col) {
			VoronoiMetricResult res = this->optimalRotations[row][col];
			VoronoiCell v;
			if (res.rotated_index == 0) {
				v = this->voronoi_vector[row];
				vcell_vector[1] = this->voronoi_vector[col];
			}
			else {
				v = this->voronoi_vector[col];
				vcell_vector[1] = this->voronoi_vector[row];
			}
			VoronoiCell v_r1 = v.getRotatedVoronoiDomain(res.axis_r1_rotation_sample, res.angle_r1_rotation_sample);
			VoronoiCell v_r2_r1 = v_r1.getRotatedVoronoiDomain(res.axis_r2_rotation_sample, res.angle_r2_rotation_sample);
			vcell_vector[0] = v_r2_r1;
			for (int i = 0; i < vcell_vector.size(); ++i) {
				final_name = path + "/" + vcell_vector[0].name + "_" + vcell_vector[1].name + "_" + std::to_string(i) + ".vtp";
				vcell_vector[i].writeToVTPFileFormat(final_name.c_str());
			}
		}
	}
}