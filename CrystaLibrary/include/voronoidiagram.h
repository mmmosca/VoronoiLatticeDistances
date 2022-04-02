#ifndef _VORONOIDIAGRAM_H
#define _VORONOIDIAGRAM_H

#include <voronoicell.h>

class VoronoiDiagram {
public:
	std::vector<VoronoiCell> voronoicells_vector;

	VoronoiDiagram(Lattice& lattice, int extend);
	std::vector<VoronoiCell> getVoronoiCellsFromPointsInside(std::vector<Eigen::VectorXd> points_inside);
	std::vector<VoronoiCell> getVoronoiCellsFromPointsInside_Ord(std::vector<Eigen::VectorXd> points_inside);
	void writeVoronoiDomainsToVTP(std::string filename);
	void writeVoronoiDomainsToOFF(std::string filename);

};
#endif // !_VORONOIDIAGRAM_H
