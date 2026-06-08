#include <gtest/gtest.h>
#include "periodicpointcloud.h"

struct PPCTestCase
{
    std::vector<Eigen::VectorXd> generators;
    unsigned int extension;
    unsigned int expected_coefficients;
    unsigned int expected_c_permutations;
    unsigned int expected_generators;
    unsigned int expected_g_permutations;
};

class PPCParameterizedTest :
    public ::testing::TestWithParam<PPCTestCase>
{
};

INSTANTIATE_TEST_SUITE_P(
    PPCTests,
    PPCParameterizedTest,
    ::testing::Values(
        PPCTestCase{
            std::vector<Eigen::VectorXd> {
                Eigen::Vector3d(1,0,0),
                Eigen::Vector3d(0,1,0),
                Eigen::Vector3d(0,0,1)
            }, 3, 4, 64, 3, 64
        },
        PPCTestCase{
            std::vector<Eigen::VectorXd> {
                Eigen::Vector2d(1,0),
                Eigen::Vector2d(0,1)
            }, 2, 3, 9, 2, 9
        }
    )
);

TEST_P(PPCParameterizedTest, Permutations)
{
    const auto& testCase = GetParam();
    PeriodicPointCloud ppc(testCase.generators.size());
    
    for (auto& g : testCase.generators) {
        ppc.addGenerator(std::move(g));
    }
    ppc.updateCoefficientsForPositiveDirection(testCase.extension);
    ppc.updateCoefficientCombinations();
	ppc.updateGeneratorCombinations();

    EXPECT_EQ(ppc.C.size(), testCase.expected_coefficients);
    EXPECT_EQ(ppc.coef_count, testCase.expected_coefficients);
    EXPECT_EQ(ppc.C_combinations.size(), testCase.expected_c_permutations);
    EXPECT_EQ(ppc.c_combination_count,  testCase.expected_c_permutations);
    EXPECT_EQ(ppc.G.size(), testCase.expected_generators);
    EXPECT_EQ(ppc.generator_count, testCase.expected_generators);
    EXPECT_EQ(ppc.G_combinations.size(), testCase.expected_g_permutations);
    EXPECT_EQ(ppc.g_combination_count, testCase.expected_g_permutations);
}
