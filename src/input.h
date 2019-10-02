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

#ifndef CDFCI_INPUT_H
#define CDFCI_INPUT_H

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <regex>
#include <exception>
#include <type_traits>
#include <algorithm>
#include "lib/json/json.hpp"
#include "determinant.h"

using Option = nlohmann::json;

Option option_from_file(const std::string file_path)
{
    Option result;
    try
    {
        std::ifstream f(file_path);
        f >> result;
    }
    catch (...)
    {
        throw std::invalid_argument("Can not read the input file " + file_path);
    }
    return result;
}


struct pair_hash {
    template <class T>
    std::size_t operator () (const std::pair<T, T> &p) const {
        return std::hash<T>{}((p.first << 10) + p.second);
    }
};


class Fcidump
{
    public:

    using value_type = double;

    /* FCIDUMP header */

    int norb;  // Number of spin-orbitals (not orbitals, consistent with UHF).
    int nelec;  // Number of electrons.
    int ms2;  // Spin
    // orbsym is not used.
    // isym is not used.
    bool uhf; // Optional (Psi4/HANDE only). True if Unrestricted Hartree-Fock (UHF) is used. False if RHF is used.

    /* FCIDUMP body (integrals) */

    // Store two-body integrals in the format of ((i, j), ((a, b), integral), where i < j.
    std::unordered_map<OrbitalPair, std::unordered_map<OrbitalPair, value_type, pair_hash>, pair_hash> two_body_integral;

    // Store one-body integrals in a norb * norb dense matrix.
    std::vector<std::vector<value_type>> one_body_integral;

    // Store the orbital energy in the format of (i, energy).
    // Note this is optional, which is provided by Psi4/HANDE but not PySCF.
    // Currently, it is only used to find the Hartree-Fock state for initialization.
    std::unordered_map<Orbital, value_type> orbital_energy;

    // The core energy.
    value_type core_energy;

    Fcidump(const std::string file_path) 
    {
        from_file(file_path);
    }

    void from_file(const std::string file_path)
    {
        try
        {
            std::ifstream f(file_path);
            if (f.fail())
            {
                throw std::invalid_argument("Failed to open the FCIDUMP file " + file_path);
            }

            // Read FCIDUMP header.
            std::string header, line;
            while (std::getline(f, line) && line.find("&END") == std::string::npos &&
             line.find("/") == std::string::npos)
            {
                header += " ";
                header += line;
            }
            header += " &END";

            if (f.eof())
            {
                throw std::invalid_argument("FCIDUMP has the wrong header.");
            }

            const std::string int_regex =  R"([ ]*=[ ]*(\d+))";
            const std::string bool_regex = R"([ ]*=[ .]*(FALSE|TRUE))";

            norb = read_parameter<int>(header, "NORB", int_regex);
            nelec = read_parameter<int>(header, "NELEC", int_regex);
            ms2 = read_parameter<int>(header, "MS2", int_regex);
            uhf = read_parameter<bool>(header, "UHF", bool_regex, false);
            if (!uhf)
            {
                norb *= 2;  // Convert orbitals to spin-orbitals for RHF.
            }

            // Initialize one_body_integral
            one_body_integral.resize(norb);
            for (auto& v : one_body_integral)
            {
                v.resize(norb, 0);
            }

            // Read FCIDUMP body (integrals).
            if (uhf)
            {
                read_integral_uhf(f);
            }
            else
            {
                read_integral_rhf(f);
            }
            return;
        }
        catch (const std::exception& e)
        {
            throw std::invalid_argument("Can not read the FCIDUMP file " + file_path + ". " + e.what());
        }
    }

    // Read required parameters
    template<typename T>
    T read_parameter(const std::string& header, const std::string name, const std::string& regex_string) const
    {
        std::regex r(name + regex_string);
        std::smatch m;
        T parameter;
        if (std::regex_search(header, m, r))
        {
            std::string parameter_string = m[1];
            // Convert to type T
            if (std::is_same<T, bool>::value)
            {
                // Convert string "TRUE" or "FALSE" (case insensitive) to bool.
                std::transform(parameter_string.begin(), parameter_string.end(), parameter_string.begin(), ::tolower);
                std::istringstream(parameter_string) >> std::boolalpha >> parameter;
            }
            else
            {
                std::istringstream(parameter_string) >> parameter;
            }
        }
        else
        {
            throw std::invalid_argument(name + " is not found.");
        }
        
        return parameter;
    }

    // Read optional parameter. Default value shoule be provided.
    template<typename T>
    T read_parameter(const std::string& header, const std::string name, const std::string& regex_string, T default_value)
    {
        T parameter;
        try
        {
            parameter = read_parameter<T>(header, name, regex_string);
        }
        catch(const std::exception& e)
        {
            // Use default value and show a warning. Do not terminate the program.
            parameter = default_value;
            std::cerr << "Warning: FCIDUMP optional parameter " << e.what() << " Use default value " << default_value\
                 << "." << std::endl << std::endl;
        }
        return parameter;   
    }

    void read_integral_uhf(std::ifstream& f)
    {
        Orbital i, j, a, b;
        value_type integral;
        while (f >> integral >> i >> a >> j >> b)
        {
            // Two-body integrals
            if (i && a && j && b)
            {
                // [symmetry]
                // <ij|g|ab>
                insert_two_body_(i - 1, j - 1, a - 1, b - 1, integral);
                // <ib|g|aj>
                insert_two_body_(i - 1, b - 1, a - 1, j - 1, integral);
                // <aj|g|ib>
                insert_two_body_(a - 1, j - 1, i - 1, b - 1, integral);
                // <ab|g|ij>
                insert_two_body_(a - 1, b - 1, i - 1, j - 1, integral);
            }
            // One-body integrals <i|f|a>
            else if (a)
            {
                one_body_integral[i - 1][a - 1] = integral;
                one_body_integral[a - 1][i - 1] = integral;
            }
            // Orbital energy
            else if (i)
            {
                orbital_energy.insert({i - 1, integral});
            }
            // Core energy
            else
            {
                core_energy = integral;
            }
        }
        return;
    }

    Orbital alpha(const Orbital i) const {return 2 * i - 2;}
    Orbital beta(const Orbital i) const {return 2 * i - 1;}

    void read_integral_rhf(std::ifstream& f)
    {
        Orbital i, j, a, b;
        value_type integral;
        while (f >> integral >> i >> a >> j >> b)
        {
            // Two-body integrals
            if (i && a && j && b)
            {
                // [symmetry]
                // <ij|g|ab>
                rhf_insert_two_body_(i, j, a, b, integral);
                // <ib|g|aj>
                rhf_insert_two_body_(i, b, a, j, integral);
                // <aj|g|ib>
                rhf_insert_two_body_(a, j, i, b, integral);
                // <ab|g|ij>
                rhf_insert_two_body_(a, b, i, j, integral);
            }
            // One-body integrals <i|f|a>
            else if (a)
            {
                one_body_integral[alpha(i)][alpha(a)] = integral;
                one_body_integral[beta(i)][beta(a)] = integral;
                one_body_integral[alpha(a)][alpha(i)] = integral;
                one_body_integral[beta(a)][beta(i)] = integral;
            }
            // Orbital energy
            else if (i)
            {
                orbital_energy.insert({alpha(i), integral});
                orbital_energy.insert({beta(i), integral});
            }
            // Core energy
            else
            {
                core_energy = integral;
            }
        }
        return;
    }

    // Helper function. Insert four combinations of spin-orbitals from orbitals. [symmetry]
    void rhf_insert_two_body_ (const Orbital i, const Orbital j, const Orbital a, const Orbital b, const value_type integral)
    {
        // Insert four spin-orbitals combination.
        // alpha-alpha
        insert_two_body_(alpha(i), alpha(j), alpha(a), alpha(b), integral);
        // beta-beta
        insert_two_body_(beta(i), beta(j), beta(a), beta(b), integral);
        // alpha-beta
        insert_two_body_(alpha(i), beta(j), alpha(a), beta(b), integral);
        // beta-alpha
        insert_two_body_(beta(i), alpha(j), beta(a), alpha(b), integral);
        return;
    }

    // Helper function. Keep i <= j.
    // Note: when i = j, <ii|ab> and <ii|ba> are the same and both inserted. [symmetry] (bug fixed)
    void insert_two_body_ (const Orbital i, const Orbital j, const Orbital a, const Orbital b, const value_type integral)
    {
        if (i <= j)
        {
            two_body_integral[{i, j}].insert({{a, b}, integral});
        }
        // Exchange i, j (and a, b accordingly)
        if (i >= j)
        {
            two_body_integral[{j, i}].insert({{b, a}, integral});
        }
        return;
    }

    // Return the two_body_integral. Read-only.
    value_type get_two_body_integral(Orbital i, Orbital j, Orbital a, Orbital b) const
    {
        try
        {
            return two_body_integral.at({i, j}).at({a, b});
        }
        catch(const std::out_of_range& e)
        {
            return 0;
        }

    }

    // Return the one_body_integral. Read-only.
    value_type get_one_body_integral(Orbital i, Orbital a) const
    {
        return one_body_integral[i][a];
    }

};

#endif
