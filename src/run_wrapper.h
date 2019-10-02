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

#ifndef CDFCI_RUN_WRAPPER_H
#define CDFCI_RUN_WRAPPER_H

#include "solver.h"

void print_header()
{
    std::cout << R"(  ______  _______   _______   ______  __  )" << std::endl;
    std::cout << R"( /      ||       \ |   ____| /      ||  | )" << std::endl;
    std::cout << R"(|  ,----'|  .--.  ||  |__   |  ,----'|  | )" << std::endl;
    std::cout << R"(|  |     |  |  |  ||   __|  |  |     |  | )" << std::endl;
    std::cout << R"(|  `----.|  '--'  ||  |     |  `----.|  | )" << std::endl;
    std::cout << R"( \______||_______/ |__|      \______||__|  )" << CDFCI_VERSION_MAJOR\
	      << '.' << CDFCI_VERSION_MINOR << '.' << CDFCI_VERSION_PATCH << std::endl;
    std::cout << std::endl;
    std::cout << "Compiled on " << __DATE__ << " at " << __TIME__ << std::endl;
    std::cout << std::endl;
    std::cout << "Authors:" << std::endl;
    std::cout << "    Yingzhou Li, Jianfeng Lu and Zhe Wang" << std::endl;
    std::cout << "    Department of Mathematics" << std::endl;
    std::cout << "    Duke University" << std::endl;
    std::cout << "Github:" << std::endl;
    std::cout << "    https://github.com/quan-tum/CDFCI/" << std::endl;
    std::cout << "=======================================================" << std::endl;
#ifdef _OPENMP
    std::cout << "OpenMP: enabled with " << omp_get_max_threads() << " threads" << std::endl;
#else
    std::cout << "OpenMP: disabled" << std::endl;
#endif
    std::cout << std::endl;
}

std::string check_argument(int argc, char *argv[])
{
    // The first argument is the path of the input json file.
    if (argc < 2)
    {
        throw std::invalid_argument("No input file specified. Please use the path "
         "of the input file as the first argument.");
    }
    std::string input_file = argv[1];
    return input_file;
}

Option read_input(const std::string input_file)
{
    // Default options.
    Option default_cdfci_option = {
        {"num_iterations", 10},
        {"report_interval", 1},
        {"coordinate_pick", "gcd_grad"},
        {"coordinate_update", "ls"},
        {"ref_det_occ", Option::array()},
        {"z_threshold", 0.0},
        {"z_threshold_search", false},
        {"max_wavefunction_size", 1},
        {"max_memory", 0.0}
    };
    Option default_hamiltonian_option = {
        {"type", "molecule"},
        {"threshold", 1e-13}
    };
    // Read the input file.
    Option option = option_from_file(input_file);

    default_cdfci_option.merge_patch(option["solver"]["cdfci"]);
    option["solver"]["cdfci"] = default_cdfci_option;
    default_hamiltonian_option.merge_patch(option["hamiltonian"]);
    option["hamiltonian"] = default_hamiltonian_option;

    // Validate the input file.
    if (option["hamiltonian"]["type"] != "molecule")
    {
        throw std::invalid_argument("hamiltonian:type is invalid. Currently the only supported "
         "Hamiltonian is \"molecule\".");
    }
    double max_memory = option["solver"]["cdfci"]["max_memory"];
    if (max_memory <= 0.0)
    {
        throw std::invalid_argument("solver:cdfci:max_memory is invalid. It is the maximum memory (GB)"
            " allowed for the wavefunction. Please set a positive number.");
    }
    bool z_threshold_search = option["solver"]["cdfci"]["z_threshold_search"];
    double z_threshold = option["solver"]["cdfci"]["z_threshold"];
    long report_interval = option["solver"]["cdfci"]["report_interval"];
    if (z_threshold_search)
    {
        std::cout << "Warning: z cut-off threshold auto search (z_threshold_search) is enabled. ";
        std::cout << "This feature is experimental and may not always work." << std::endl << std::endl;
        const double start_z_threshold = 1e-10;
        if (z_threshold < start_z_threshold)
        {
            z_threshold = start_z_threshold;
            std::cout << "Warning: The initial z_threshold may be too small. Set z_threshold = " << std::scientific \
             << std::setprecision(2) << start_z_threshold << std::endl << std::endl;
        }
    }
    if (z_threshold_search && (report_interval < 1000))
    {
        throw std::invalid_argument("Please set report_interval >= 1000 when z_threshold_search = true.");
    }

    return option;
}

template <int N>
double run_cdfci(Option& option, Fcidump& fci)
{
    // Step 1: Initialize the Hamiltonian.
    // Currently only molecules are supported.
    HamiltonianMolecule<N> h(fci, option["hamiltonian"]["threshold"]);

    // Step 2: Initialize the CDFCI solver.
    Option cdfci_option = option["solver"]["cdfci"];
    CDFCI<N> cdfci_solver(cdfci_option);
    CoordinatePickGcdGrad<N> coord_pick;
    CoordinateUpdateLS<N> coord_ls;
    Determinant<N>::constuct_masks();

    // Step 3: Solve.
    cdfci_solver.solve(h, coord_pick, coord_ls);
    return cdfci_solver.result_energy;
}

double run_cdfci_wrapper(Option& option)
{
    // Read FCIDUMP.
    Fcidump fci(option["hamiltonian"]["fcidump_path"]);

    // Calculate number of size_t to use to store the determinant.
    int n = 1 + (fci.norb - 1) / (sizeof(Determinant<1>::det_type) * CHAR_BIT);

    // Print basic information.
    std::cout << "FCIDUMP information" << std::endl;
    std::cout << "-------------------" << std::endl;
    std::cout << "Number of elections: " << fci.nelec << std::endl;
    std::cout << "Number of spin-orbitals: " << fci.norb << std::endl;
    std::cout << "Ms2: " << fci.ms2 << std::endl;
    std::string hf_type = fci.uhf ? "unrestricted" : "restricted";
    std::cout << "Hartree Fock: " + hf_type << std::endl << std::endl;
    std::cout << "Machine information" << std::endl;
    std::cout << "-------------------" << std::endl;
    std::cout << "sizeof(size_t): "  << sizeof(Determinant<1>::det_type) << " bytes or ";
    std::cout << (sizeof(Determinant<1>::det_type) * CHAR_BIT) << " bits" << std::endl;
    std::cout << "Number of size_t to represent one Slater determinant: " << n << std::endl;
#ifdef CDFCI_USE_LONG_DOUBLE
    std::cout << "__float128: disabled" << std::endl;
    std::cout << "Warning: __float128 may not be supported and long double is used instead. Long time iteration (> 1e9) may accumulate significant errors. Please use g++/icc or modify src/config.h to use the right definition." << std::endl << std::endl;
#else
    std::cout << "__float128: enabled" << std::endl << std::endl;
#endif

    // Compute max_wavefunction_size
    int bytes_per_det = n * sizeof(Determinant<1>::det_type) + sizeof(WaveFunction<1>::value_type);
    size_t max_wavefunction_size = 8;
    double max_load_factor = 0.79;
    double max_memory = option["solver"]["cdfci"]["max_memory"];
    while (max_wavefunction_size * bytes_per_det < max_memory * 1073741824)
    {
        max_wavefunction_size *= 2;
    }
    max_wavefunction_size /= 2;
    option["solver"]["cdfci"]["max_wavefunction_size"] = static_cast<size_t>(max_load_factor * max_wavefunction_size);

    std::cout << "Input option" << std::endl;
    std::cout << "------------" << std::endl;
    std::cout << option.dump(4) << std::endl;
    std::cout << "Note:" << std::endl;
    std::cout << "In the wavefunction, each determinant needs " << bytes_per_det << " bytes." << std::endl;
    std::cout << "The size of the wavefunction is 2^" << static_cast<int>(log2(max_wavefunction_size)) << " = " << max_wavefunction_size\
        << ", which can store " << static_cast<size_t>(max_load_factor * max_wavefunction_size) << " determinants at most." << std::endl;
    std::cout << "The wavefunction uses about " << std::fixed << std::setprecision(3) << \
        static_cast<double>(max_wavefunction_size) * bytes_per_det / 1073741824 << " GB memory." << std::endl << std::endl;

    // Run CDFCI
    switch (n)
    {
    case 1:
        return run_cdfci<1>(option, fci);
    case 2:
        return run_cdfci<2>(option, fci);
    case 3:
        return run_cdfci<3>(option, fci);
    case 4:
        return run_cdfci<4>(option, fci);
    default:
        throw std::invalid_argument("The CDFCI program only supports 1 - 4 size_t integers for a Slater determinant. \
         Please refer to the documentation, modify, recompile the code and run again.");
        break;
    }
}


#endif
