#include "tests.h"
#include "periodicpointcloud.h"

TEST(PeriodicPointCloudTest, Permutations)
{
    PeriodicPointCloud ppc(3);
    int extend = 3, expected_permutations;
    
    ppc.addGenerator(Eigen::Vector3d(1,0,0));
    ppc.addGenerator(Eigen::Vector3d(0,1,0));
    ppc.addGenerator(Eigen::Vector3d(0,0,1));
    ppc.updateCoefficientsForPositiveDirection(extend);
    ppc.updateCoefficientCombinations();
	ppc.updateGeneratorCombinations();

    expected_permutations = static_cast<int>(std::pow(extend+1, ppc.G.size()));

    EXPECT_EQ(ppc.C.size(), 4);
    EXPECT_EQ(ppc.coef_count, 4);
    EXPECT_EQ(ppc.C_combinations.size(), expected_permutations);
    EXPECT_EQ(ppc.c_combination_count, expected_permutations);
    EXPECT_EQ(ppc.G.size(), 3);
    EXPECT_EQ(ppc.generator_count, 3);
    EXPECT_EQ(ppc.G_combinations.size(), expected_permutations);
    EXPECT_EQ(ppc.g_combination_count, expected_permutations);
}
