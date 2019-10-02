//   ______  _______   _______   ______  __
//  /      ||       \ |   ____| /      ||  |
// |  ,----'|  .--.  ||  |__   |  ,----'|  |
// |  |     |  |  |  ||   __|  |  |     |  |
// |  `----.|  '--'  ||  |     |  `----.|  |
//  \______||_______/ |__|      \______||__|
//
// Coordinate Descent Full Configuration Interaction (CDFCI) package in C++14
// https://github.com/quan-tum/CDFCI
//
// Copyright (c) 2019, Zhe Wang, Yingzhou Li and Jianfeng Lu
// All rights reserved.
//
// This source code is licensed under the BSD 3-Clause License found in the
// LICENSE file in the root directory of this source tree.

#include "test.h"
#include "../../src/run_wrapper.h"

bool test_system(std::string sys_name, std::string path)
{
    // Read input
    std::string input_file = path + "/" + sys_name + "/input.json";
    Option option = read_input(input_file);
    std::string fcidump_path = option["hamiltonian"]["fcidump_path"];
    option["hamiltonian"]["fcidump_path"] = path + "/" + sys_name + "/" + fcidump_path;
    // Run CDFCI
    double energy = run_cdfci_wrapper(option);
    // Read reference
    std::ifstream f(path + "/" + sys_name + "/" + "result.txt");
    double ref_energy, ref_error;
    f >> ref_energy >> ref_error;
    // Compare
    return (fabs(ref_energy - energy) < ref_error);
}

TEST_CASE("Test example systems")
{
    const std::string test_path = "example";

    SUBCASE("h2o_ccpvdz")
    {
        CHECK(test_system("h2o_ccpvdz", test_path));
    }
    SUBCASE("h2o_alternative")
    {
        CHECK(test_system("h2o_alternative", test_path));
    }
    SUBCASE("h2o_sto3g")
    {
        CHECK(test_system("h2o_sto3g", test_path));
    }
    SUBCASE("c2_ccpvdz")
    {
        CHECK(test_system("c2_ccpvdz", test_path));
    }
    SUBCASE("n2_ccpvdz")
    {
        CHECK(test_system("n2_ccpvdz", test_path));
    }
    SUBCASE("n2_ccpvdz_eps1e-2")
    {
        CHECK(test_system("n2_ccpvdz_eps1e-2", test_path));
    }
    SUBCASE("n2_ccpvdz_triplet")
    {
        CHECK(test_system("n2_ccpvdz_triplet", test_path));
    }
    SUBCASE("cr2_ahlrichs")
    {
        CHECK(test_system("cr2_ahlrichs", test_path));
    }
}