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
#include <cmd.h>
#include <fileio.h>
#include <boost/filesystem.hpp>
#include <chrono>
#include <ctime>
#include "logging.h"

std::string usage{"|--- USAGE ---|:\nVoronoi Computation is performed on a 3x3x3 extended domain of Niggli's Reduced Cell\n"\
					"Required options:\n\t"\
					"-inputdir [Input Folder with CIF files]\n\t"\
					"-outputdir [Output Folder to write Results]\n"\
					"Output options:\n\t"\
					"-ds:\tOutputs .csv file with Scale Invariant Distance matrix (n x n) between all n crystal lattices\n\t"\
					"-dh:\tOutputs .csv file with Extended Hausdorff Distance matrix (n x n) between all n crystal lattices\n"\
					"Optional commands:\n\t"\
					"-intervals [integer n (default n=2)]:\tIt affects the number of rotation samples to be considered for metric computations (number of rotations: 4*pi^2*n^3)\n\t"\
					"-threads [integer t (default t=1)]:\tRotation samples are divided among t threads \n"\
					"-debug\tEnable debug message logging\n"};

int main(int argc, char* argv[]) {
	std::string input_dir, output_dir;
	CommandLine cmd;
	int nthreads = 1;
	int intervals = 2;
	bool INPUTDIR_OPT = false;
	bool OUTPUTDIR_OPT = false;
	bool CONFIG_OPT = false;
	bool EXTHAUSDORFF_OPT = false;
	bool SCALEINV_OPT = false;
	bool OUTPUT_OPT = false;
	bool DEBUG_OPT = false;

	char* w;
	while ((w = cmd.mygetoptW(argc, argv, "inputdir:|outputdir:|dh|ds|intervals:|threads:|debug|help|")) != NULL) {
		if (strcmp(w, "inputdir") == 0) {
			input_dir.assign(cmd.myoptarg);
			INPUTDIR_OPT = true;
			continue;
		}
		if (strcmp(w, "outputdir") == 0) {
			output_dir.assign(cmd.myoptarg);
			OUTPUTDIR_OPT = true;
			continue;
		}
		if (strcmp(w, "dh") == 0) {
			EXTHAUSDORFF_OPT = true;
			OUTPUT_OPT = OUTPUT_OPT || EXTHAUSDORFF_OPT;
			continue;
		}
		if (strcmp(w, "ds") == 0) {
			SCALEINV_OPT = true;
			OUTPUT_OPT = OUTPUT_OPT || SCALEINV_OPT;
			continue;
		}
		if (strcmp(w, "intervals") == 0) {
			intervals = stoi(cmd.myoptarg);
			continue;
		}
		if (strcmp(w, "threads") == 0) {
			nthreads = std::stoi(cmd.myoptarg);
			continue;
		}
		if (strcmp(w, "debug") == 0) {
			DEBUG_OPT = true;
			continue;
		}
		if (strcmp(w, "help") == 0) {
			Logger::info("USAGE\n" + usage);
			exit(EXIT_SUCCESS);
		}
	}

	if (!(INPUTDIR_OPT && OUTPUTDIR_OPT)) {
		Logger::error("No input and output directory specified\n" + usage);
		exit(EXIT_FAILURE);
	}

	if (!OUTPUT_OPT) {
		Logger::error("No option for results detected. Please choose an option!\n" + usage);
		exit(EXIT_FAILURE);
	}

	if (DEBUG_OPT) {
		Logger::setLevel(LogLevel::LOGDEBUG);
	} else {
		Logger::setLevel(LogLevel::LOGINFO);
	}

	Logger::info("Loading data from OFF files...");
	auto start = std::chrono::system_clock::now();

	VoronoiCellMetrics voronoi_metrics;
	int file_count = 0;
	for (boost::filesystem::directory_iterator itr(input_dir); itr != boost::filesystem::directory_iterator(); ++itr)
	{
		std::string filename(itr->path().string());
		filename.erase(0, filename.find_last_of("\\/") + 1);
		filename.erase(filename.find_last_of('.'), filename.length());
		if (strcmp(itr->path().extension().string().c_str(), ".off") == 0) {
			Logger::info(std::to_string(++file_count) + " - Loading OFF file: " + filename);
			OFFFile off_reader(itr->path().string());
			VoronoiCell v;
			v.setPoints(off_reader.getPoints());
			v.setFaces(off_reader.getFaces());
			v.setSite(std::vector<double>{0.0, 0.0, 0.0});
			v.setName(std::string(filename));
			voronoi_metrics.addVoronoiDomain(v, filename);
		}
	}
	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	int elapsed_s = elapsed_seconds.count();
	int hours = floor(elapsed_s / 3600), minutes = floor(fmod(elapsed_s, 3600) / 60), seconds = fmod(elapsed_s, 60);
	Logger::info("Elapsed time: " + std::to_string(hours) + ':' + std::to_string(minutes) + ':' + std::to_string(seconds));

	/******************************************************\
	|******** POST-PROCESSING: METRICS COMPUTATION ********|
	\******************************************************/

	if (!(SCALEINV_OPT || EXTHAUSDORFF_OPT)) return 0;
	Logger::info("Computing requested metrics...");
	start = std::chrono::system_clock::now();

	voronoi_metrics.setSamplingMode(SamplingMode::UNIFORM);
	voronoi_metrics.setThreads(nthreads);
	voronoi_metrics.generateUniformRotationSamples(intervals);

	if (EXTHAUSDORFF_OPT) {
		voronoi_metrics.updateExtendedHausdorffDistanceMatrix();
		voronoi_metrics.writeExtendedHausdorffDistanceMatrixToCSV(output_dir + std::string("/exthausdorffdistancematrix"));
	}

	if (SCALEINV_OPT) {
		voronoi_metrics.updateScaleInvariantDistanceMatrix();
		voronoi_metrics.writeScaleInvariantDistanceMatrixToCSV(output_dir + std::string("/scaleinvariantdistancematrix"));
	}

	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	elapsed_s = elapsed_seconds.count();
	hours = floor(elapsed_s / 3600);
	minutes = floor(fmod(elapsed_s, 3600) / 60);
	seconds = fmod(elapsed_s, 60);
	Logger::info("Elapsed time: " + std::to_string(hours) + ':' + std::to_string(minutes) + ':' + std::to_string(seconds));

	return(0);
}
