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

#include "run_wrapper.h"

int main(int argc, char *argv[])
{
    try
    {
        print_header();
        std::string input_file = check_argument(argc, argv);
        Option option = read_input(input_file);
        run_cdfci_wrapper(option);
        return 0;
    }
    catch(const std::exception& e)
    {
        std::cerr << std::endl << "ERROR: " << e.what() << std::endl;
        return 1;
    }
    
}
