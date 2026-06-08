#include <gtest/gtest.h>
#include "crystalstructure.h"

TEST(LatticeTest, spanTheLattice)
{
    Lattice lattice;
    int extend = 3, expected_permutations;
    
    lattice.setCellParameters({10, 10, 10, 90, 90, 90});
    lattice.updateParametersToOrginalMode();
    lattice.spanTheLattice(extend);
    expected_permutations = static_cast<int>(std::pow(extend+1, lattice.G.size()));

    EXPECT_EQ(lattice.C.size(), 4);
    EXPECT_EQ(lattice.coef_count, 4);
    EXPECT_EQ(lattice.C_combinations.size(), expected_permutations);
    EXPECT_EQ(lattice.c_combination_count, expected_permutations);
    EXPECT_EQ(lattice.G.size(), 3);
    EXPECT_EQ(lattice.generator_count, 3);
    EXPECT_EQ(lattice.G_combinations.size(), expected_permutations);
    EXPECT_EQ(lattice.g_combination_count, expected_permutations);
}
