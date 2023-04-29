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
#ifndef _VORONOICELL_H
#define _VORONOICELL_H

#include <vtkcontext.h>
#include <crystalstructure.h>
#include <thread>
#include <map>

class VoronoiCell {
private:
	std::vector<Point_CM> points;
	std::vector<std::vector<int>> faces;
	std::map<int, int> point_counter;
	CGAL::Polyhedron_3<Kernel> polyhedron;
	
public:
	std::string name;
	Eigen::VectorXd site;
	VoronoiCell() {};

	std::vector<Point_CM> getPoints();
	void setPoints(std::vector<Point_CM> cell_points);
	void setPoints(std::vector<std::vector<double>> cell_points);
	std::vector<std::vector<int>> getFaces();
	void setFaces(std::vector<std::vector<int>> cell_faces);
	void setSite(Eigen::VectorXd& voronoi_site);
	void setSite(std::vector<double> voronoi_site);
	void setName(std::string& s);
	void updatePolyhedron();
	bool hasCorrectGeometry();
	void fixGeometry();
	/*	Every vertex must appear at least in 3 faces of the Voronoi Domain */
	void updatePointCounter();
	double getVolume();
	/* 
	 * rotation_sample: (theta, mu, z)
	 * theta: angle about axis of rotation
	 * mu: angle about Z axis
	 * z: height of rotation axis
	 */
	std::vector<Point_CM> getRotatedPointsAboutAngleAndAxis(std::vector<double> rotation_sample);
	std::vector<Point_CM> getRotatedPointsAboutAngleAndAxis(std::vector<double> rotationaxis_sample, double anglerad_sample);
	VoronoiCell getRotatedVoronoiDomain(std::vector<double> rotationaxis_sample, double anglerad_sample);
	bool isPointInside(Eigen::VectorXd& query);
	std::vector<int> getPointIndecesInside(std::vector<Eigen::VectorXd> &queries);
	VoronoiCell getShiftedToOriginVoronoiDomain();
	VoronoiCell getVoronoiDomainShiftedBy(Eigen::VectorXd offset_vector);
	void writeToVTPFileFormat(const char* filename);
	void writePointsToCSVFileFormat(const char* filename);
	void writeToOFFFileFormat(const char* filename);
};

struct VoronoiMetricResult {
	double distance_result;
	std::vector<double> rotation_sample;
	std::vector<double> axis_r1_rotation_sample;
	std::vector<double> axis_r2_rotation_sample;
	double angle_r1_rotation_sample;
	double angle_r2_rotation_sample;

	std::vector<std::string> filename_vector;
	std::vector<double> distance_vector;
	std::vector<std::vector<double>> sample_vector;
	int rotated_index;
};

enum class SamplingMode {UNIFORM};

class VoronoiCellMetrics {
private:
	int threads;
	SamplingMode sampling_mode;
	std::vector<std::string> names;
	std::vector<VoronoiCell> voronoi_vector;
	std::vector<std::vector<double>> rotation_samples;
	std::vector<std::vector<double>> axis_r1_rotation_samples;
	std::vector<std::vector<double>> axis_r2_rotation_samples;
	std::vector<double> angle_r1_rotation_samples;
	std::vector<double> angle_r2_rotation_samples;
	Eigen::MatrixXd extHausdorffDistMatrix, scaleInvDistMatrix;

	/*********************\
	|****** UNIFORM ******|
	\*********************/

	void getOffsetOverRotation(VoronoiCell v1, VoronoiCell v2, std::vector<double> rotation_sample, VoronoiMetricResult &vormetric_result);
	void getOffsetOverRotationsSet(VoronoiCell v1, VoronoiCell v2, int thread_id, int start, int end, VoronoiMetricResult& vormetric_result);
	VoronoiMetricResult getOffset(VoronoiCell v1, VoronoiCell v2);
	VoronoiMetricResult getOffset_Parallel_std(VoronoiCell v1, VoronoiCell v2);

	void getScaleOverRotation(VoronoiCell v1, VoronoiCell v2, std::vector<double> rotation_sample, VoronoiMetricResult &vormetric_result);
	void getScaleOverRotationsSet(VoronoiCell v1, VoronoiCell v2, int thread_id, int start, int end, VoronoiMetricResult& vormetricresult_vector);
	VoronoiMetricResult getScale(VoronoiCell v1, VoronoiCell v2);
	VoronoiMetricResult getScale_Parallel_std(VoronoiCell v1, VoronoiCell v2);

public:
	std::vector<std::vector<VoronoiMetricResult>> optimalRotations;

	void setThreads(int n_threads);
	void addVoronoiDomain(VoronoiCell& voronoi_domain, std::string name = "");
	void setSamplingMode(SamplingMode s_mode);
	void generateUniformRotationSamples(int intervals);

	VoronoiMetricResult getExtendedHausdorffDistance(VoronoiCell v1, VoronoiCell v2);
	VoronoiMetricResult getScaleInvariantDistance(VoronoiCell v1, VoronoiCell v2);

	VoronoiMetricResult getExtendedHausdorffDistance_Approximate(VoronoiCell v1, VoronoiCell v2);
	VoronoiMetricResult getScaleInvariantDistance_Approximate(VoronoiCell v1, VoronoiCell v2);

	Eigen::MatrixXd getExtendedHausdorffDistanceMatrix();
	void updateExtendedHausdorffDistanceMatrix();
	Eigen::MatrixXd getScaleInvariantDistanceMatrix();
	void updateScaleInvariantDistanceMatrix();	

	std::vector<VoronoiCell> getRIDOptimalRotations(int row, int col);

	void writeExtendedHausdorffDistanceMatrixToCSV(std::string filename);
	void writeScaleInvariantDistanceMatrixToCSV(std::string filename);

	void writeOptimalRotationsToVTP(std::string path);
};
#endif //!_VORONOICELL_H