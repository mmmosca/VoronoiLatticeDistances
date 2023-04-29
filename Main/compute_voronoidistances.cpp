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
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>

#include <boost/filesystem.hpp>
#include <voronoidiagram.h>
#include <cmd.h>
#include <cifio.h>
#include <gemmi/dirwalk.hpp>
#include <gemmi/smcif.hpp>
#include <chrono>
#include <ctime>



#define PRINTUSAGE() printf("|--- USAGE ---|:\nVoronoi Computation is performed on a 3x3x3 extended domain of Niggli's Reduced Cell\n"\
							"Required options:\n\t"\
							"-inputdir [Input Folder with CIF files]\n\t"\
							"-outputdir [Output Folder to write Results]\n"\
							"Output options:\n\t"\
							"-vol:\tOutputs .csv file with Voronoi Cell volume\n\t"\
							"-vtp:\tOutputs .vtp file with Voronoi Cell of a Lattice (Paraview format file)\n\t"\
							"-csv:\tOutputs .csv file with Lattice and Voronoi Cell points\n\t"\
							"-off:\tOutputs .off file with Voronoi Cell verteces and faces\n\t"\
							"-ds:\tOutputs .csv file with Scale Invariant Distance matrix (n x n) between all n crystal lattices\n\t"\
							"-dh:\tOutputs .csv file with Extended Hausdorff Distance matrix (n x n) between all n crystal lattices\n"\
							"Optional commands:\n\t"\
							"-intervals [integer n (default n=2)]:\tIt affects the number of rotation samples to be considered for metric computations (number of rotations: 4*pi^2*n^3)\n\t"\
							"-threads [integer t (default t=1)]:\tRotation samples are divided among t threads \n");


int main(int argc, char* argv[]) {
	std::string input_dir, output_dir;
	CommandLine cmd;
	bool INPUTDIR_OPT = false;
	bool OUTPUTDIR_OPT = false;
	bool VOLUME_OPT = false;
	bool EXTHAUSDORFF_OPT = false;
	bool SCALEINV_OPT = false;
	bool VTP_OPT = false;
	bool CSV_OPT = false;
	bool OFF_OPT = false;
	bool OUTPUT_OPT = false;

	int intervals = 2;
	int n_threads = 1;

	char* w;
	while ((w = cmd.mygetoptW(argc, argv, "inputdir:|outputdir:|vol|dh|ds|vtp|csv|off|intervals:|threads:|")) != NULL) {
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
		if (strcmp(w, "vol") == 0) {
			VOLUME_OPT = true;
			OUTPUT_OPT = OUTPUT_OPT || VOLUME_OPT;
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
			n_threads = stoi(cmd.myoptarg);
			continue;
		}
		if (strcmp(w, "vtp") == 0) {
			VTP_OPT = true;
			OUTPUT_OPT = OUTPUT_OPT || VTP_OPT;
			continue;
		}
		if (strcmp(w, "csv") == 0) {
			CSV_OPT = true;
			OUTPUT_OPT = OUTPUT_OPT || CSV_OPT;
			continue;
		}
		if (strcmp(w, "off") == 0) {
			OFF_OPT = true;
			OUTPUT_OPT = OUTPUT_OPT || OFF_OPT;
			continue;
		}
	}

	if (!(INPUTDIR_OPT && OUTPUTDIR_OPT)) {
		PRINTUSAGE();
		exit(EXIT_SUCCESS);
	}

	if (!OUTPUT_OPT) {
		std::cout << "No option for results detected. Please choose an option!\n" << std::endl;
		PRINTUSAGE();
		exit(EXIT_SUCCESS);

	}

	std::string output_dir_vtp(output_dir + "/Voronoi_VTP"),
		output_dir_csv(output_dir + "/Voronoi_CSV"),
		output_dir_latticecsv(output_dir + "/Lattice_CSV"),
		output_dir_off(output_dir + "/Voronoi_OFF"),
		volume_file_name(output_dir + "/voronoi_volumes.csv");

	std::ofstream out_volume;

	if (VTP_OPT) {
		boost::filesystem::create_directory(boost::filesystem::path(output_dir_vtp));
	}
	if (CSV_OPT) {
		boost::filesystem::create_directory(boost::filesystem::path(output_dir_csv));
		boost::filesystem::create_directory(boost::filesystem::path(output_dir_latticecsv));
	}
	if (OFF_OPT) {
		boost::filesystem::create_directory(boost::filesystem::path(output_dir_off));
	}
	if (VOLUME_OPT) {
		out_volume.open(volume_file_name);
		out_volume << "ID,VORONOI_VOLUME\n";
	}
	
	VoronoiCellMetrics voronoi_metrics;

	auto start = std::chrono::system_clock::now();
	std::cout << "Computing Voronoi Cells..." << std::endl;

	int i = 0;
	for (std::string& cif_filepath : gemmi::CifWalk(input_dir)) {
		gemmi::cif::Document d = gemmi::cif::read_file(cif_filepath);
		int block_index;
		if ((block_index = CIFHandler::findValidBlockIndex(d)) == -1) {
			std::cerr << "File corrupted: " + cif_filepath << std::endl;
			continue;
		}

		gemmi::SmallStructure crystal_structure = gemmi::make_small_structure_from_block(d.blocks[block_index]);
		gemmi::UnitCell unit_cell = crystal_structure.cell;
		std::vector<double> cell_parameters{ unit_cell.a, unit_cell.b ,unit_cell.c, unit_cell.alpha, unit_cell.beta, unit_cell.gamma };

		// Filter filename
		std::string filename(cif_filepath);
		filename.erase(0, filename.find_last_of("\\/") + 1);
		filename.erase(filename.find_last_of('.'), filename.length());

		Lattice lattice;
		std::cout << "\r" << ++i << " - Processing file: " << filename << std::flush;
		lattice.setCellParameters(cell_parameters);

		VoronoiDiagram vd(lattice, 3);
		 
		std::vector<VoronoiCell> vcell_vector = vd.getVoronoiCellsFromPointsInside_Ord(lattice.getPointCloud());

		VoronoiCell v_cell = vcell_vector[0].getShiftedToOriginVoronoiDomain();
		v_cell.setName(filename);

		if (VOLUME_OPT) {
			out_volume << filename << ',' << cbrt(v_cell.getVolume()) << std::endl;
		}
		if (VTP_OPT) {
			std::string vtp_file(output_dir_vtp + "/" + filename + ".vtp");
			v_cell.writeToVTPFileFormat(vtp_file.c_str());
		}
		if (CSV_OPT) {
			std::string csv_file(output_dir_csv + "/" + filename + "_vorvertices.csv"),
				latticecsv_file(output_dir_latticecsv + "/" + filename + "_latticepoints.csv");
			lattice.writePointsToCSVFormatFile(latticecsv_file.c_str());
			v_cell.writePointsToCSVFileFormat(csv_file.c_str());
		}
		if (OFF_OPT) {
			std::string off_file(output_dir_off + "/" + filename + ".off");
			v_cell.writeToOFFFileFormat(off_file.c_str());
		}
		if (SCALEINV_OPT || EXTHAUSDORFF_OPT) {
			voronoi_metrics.addVoronoiDomain(v_cell, filename);
		}
	}

	if (VOLUME_OPT) {
		out_volume.close();
	}

	std::cout << std::endl;
	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	int elapsed_s = elapsed_seconds.count();
	int hours = floor(elapsed_s / 3600), minutes = floor(fmod(elapsed_s, 3600) / 60), seconds = fmod(elapsed_s, 60);
	std::cout << "Elapsed time: " << hours << ':' << minutes << ':' << seconds << "s\n";

	/***********************\
	|*** POST-PROCESSING ***|
	\***********************/

	if (!(SCALEINV_OPT || EXTHAUSDORFF_OPT)) return 0;

	voronoi_metrics.setThreads(n_threads);
	voronoi_metrics.setSamplingMode(SamplingMode::UNIFORM);
	voronoi_metrics.generateUniformRotationSamples(intervals);

	std::cout << "Computing requested metrics..." << std::endl;
	start = std::chrono::system_clock::now();

	if (EXTHAUSDORFF_OPT) {
		voronoi_metrics.updateExtendedHausdorffDistanceMatrix();
		voronoi_metrics.writeExtendedHausdorffDistanceMatrixToCSV(output_dir + std::string("/exthausdorffdistancematrix"));
	}

	if (SCALEINV_OPT) {
		voronoi_metrics.updateScaleInvariantDistanceMatrix();
		voronoi_metrics.writeScaleInvariantDistanceMatrixToCSV(output_dir + std::string("/scaleinvariantdistancematrix"));
	}

	std::cout << std::endl;
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	elapsed_s = elapsed_seconds.count();
	hours = floor(elapsed_s / 3600);
	minutes = floor(fmod(elapsed_s, 3600) / 60);
	seconds = fmod(elapsed_s, 60);
	std::cout << "Elapsed time: " << hours << ':' << minutes << ':' << seconds << "s\n";

	return 0;
}