#include <gtest/gtest.h>
#include "voronoidiagram.h"
#include "crystalstructure.h"

struct VoronoiDiagramTestCase
{
    std::vector<double> cell_parameters;
    unsigned int extension;
    unsigned int voronoicells_count;
};

class VDParameterizedTest :
    public ::testing::TestWithParam<VoronoiDiagramTestCase>
{
};

INSTANTIATE_TEST_SUITE_P(
    VDTests,
    VDParameterizedTest,
    ::testing::Values(
        VoronoiDiagramTestCase{
            std::vector<double> { 
                10, 10, 10, 90, 90, 90 
            }, 3, 8
        },
        VoronoiDiagramTestCase{
            std::vector<double> { 
                10, 10, 10, 90, 90, 90
            }, 2, 1
        }
    )
);

TEST_P(VDParameterizedTest, CreateVoronoiDiagramFromALattice)
{
    Lattice lattice;
    const auto& testCase = GetParam();

    lattice.setCellParameters(testCase.cell_parameters);
    VoronoiDiagram vd(lattice, testCase.extension);

    EXPECT_EQ(vd.voronoicells_vector.size(), testCase.voronoicells_count);
}
