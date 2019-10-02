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

#ifndef CDFCI_HAMILTONIAN_H
#define CDFCI_HAMILTONIAN_H 1

#include <algorithm>
#include "input.h"
#include "determinant.h"

template<int N = 1>
class Hamiltonian
{
    public:

    using value_type = double;
    using Column = std::vector<std::pair<Determinant<N>, value_type>>;

    int norb;
    int nelec;
    int ms2;
    bool uhf;
    
    Hamiltonian() : norb{0}, nelec{0}, ms2{0}, uhf{false} {};

    virtual value_type get_diagonal(DeterminantDecoded<N>& det) const = 0;
    virtual Column get_column(DeterminantDecoded<N>& det) const = 0;
    virtual Determinant<N> get_hartree_fock() const = 0;
    virtual Column get_column_diagonal_single(DeterminantDecoded<N>& det) const = 0;
    virtual void get_double_excitation(DeterminantDecoded<N>& det, Column& result, int idx_i, int idx_j) const = 0;

    virtual ~Hamiltonian() {};  // TODO: how to use the destructor?
};

template<int N = 1>
class HamiltonianMolecule : public Hamiltonian<N>
{
    using OrbitalList = std::vector<Orbital>;
    using typename Hamiltonian<N>::value_type;
    using typename Hamiltonian<N>::Column;
    using Hamiltonian<N>::norb;
    using Hamiltonian<N>::nelec;
    using Hamiltonian<N>::ms2;
    using Hamiltonian<N>::uhf;

    private:

    // Double excitation data structure.
    using DoubleExcitationEntry = std::pair<OrbitalPair, value_type>;
    std::vector<std::vector<DoubleExcitationEntry>> double_excitation;  // norb * norb

    // Single excitation data structure.
    struct SingleExcitationEntry
    {
        Orbital a;
        value_type one_body_int;  // <i|f|a>
        std::vector<value_type> two_body_int;  // <ik||ak>

        SingleExcitationEntry(int norb) : a{0}, one_body_int{0}
        {
            two_body_int.resize(norb);
        }
    };
    std::vector<std::vector<SingleExcitationEntry>> single_excitation;

    // Diagonal data structure.
    struct DiagonalData
    {
        std::vector<value_type> one_body_integral;  // <i|f|i>, size = norb
        std::vector<value_type> two_body_integral;  // <ij||ij>, size = norb * norb
    };
    DiagonalData diagonal;

    // Other members    
    value_type core_energy;
    std::unordered_map<Orbital, value_type> orbital_energy;

    public:

    using Hamiltonian<N>::Hamiltonian;

    HamiltonianMolecule(Fcidump& fci, value_type threshold = 0.0)
    {
        norb = fci.norb;
        nelec = fci.nelec;
        ms2 = fci.ms2;
        uhf = fci.uhf;
        core_energy = fci.core_energy;

        construct_double_excitation(fci, threshold);
        construct_single_excitation(fci, threshold);
        construct_diagonal(fci);
        orbital_energy = fci.orbital_energy;
    }

    int index(const Orbital i, const Orbital j) const {return i * norb + j;}

    // Construct double excitation <ij||ab> = <ij|g|ab> - <ij|g|ba>.
    // Only store i < j. (In Fcidump, i <= j).
    void construct_double_excitation(const Fcidump& fci, double threshold)
    {
        // Initialize
        double_excitation.resize(norb * norb);

        for (const auto& ij_entry : fci.two_body_integral)
        {
            auto ij = ij_entry.first;
            auto i = ij.first;
            auto j = ij.second;
            if (i == j)
            {
                continue;  // <ii||ab> is always zero.
            }

            for (const auto& ab_entry: ij_entry.second)  // key = (a, b), value = integral
            {
                auto ab = ab_entry.first;
                auto a = ab.first;
                auto b = ab.second;
                auto integral_ijab = ab_entry.second;

                // opposite string, simply <ij|g|ab>. Note a > b and a < b are both allowed.
                // i % 2 = 0 if i is alpha. i % 2 = 1 if i is beta.
                if ((i % 2) ^ (j % 2))
                {
                    if (fabs(integral_ijab) > threshold)
                    {
                        double_excitation[index(i, j)].push_back(ab_entry);
                    }
                }
                // same string, <ij||ab> = <ij|g|ab> - <ij|g|ba>
                // only store a < b
                else if (a < b)
                {
                    auto integral = integral_ijab - fci.get_two_body_integral(i, j, b, a);  // Q: Why can not const fci? A: fixed. Cannot use [] for unordered_map.
                    std::pair<OrbitalPair, value_type> ab_int {ab, integral};
                    if (fabs(integral) > threshold)
                    {
                        double_excitation[index(i, j)].push_back(ab_int);
                    }
                }
                else if ((b < a) && (fci.two_body_integral.at(ij).find({b, a}) == fci.two_body_integral.at(ij).end()))
                {
                    auto integral = -integral_ijab;  // <ij|g|min(a,b) max(a,b)> = 0 but <ij|g|max(a,b) min(a,b)> != 0
                    std::pair<OrbitalPair, value_type> ab_int {{b, a}, integral};
                    if (fabs(integral) > threshold)
                    {
                        double_excitation[index(i, j)].push_back(ab_int);
                    }
                }
            }
        }
        // Sort
        const auto comp = [](auto& x, auto& y) {return fabs(x.second) > fabs(y.second);};
        for (auto& ij_column : double_excitation)
            std::sort(ij_column.begin(), ij_column.end(), comp);

        return;
    }

    // Construct single excitation
    void construct_single_excitation(const Fcidump& fci, double threshold)
    {
        // Initialization
        single_excitation.resize(norb);

        // for (const auto& i_entry : fci.one_body_integral)
        for (auto i = 0; i < norb; ++i)
        {
            for(auto a = 0; a < norb; ++a)
            {
                SingleExcitationEntry entry(norb);
                entry.a = a;
                entry.one_body_int = fci.get_one_body_integral(i, a);

                auto max_val = fabs(entry.one_body_int);

                // Find all double integrals <ik|g|ak> and <ik||ak> = <ik|g|ak> - <ik|g|ka>
                for (auto k = 0; k < norb; ++k)
                {
                    value_type integral;
                    if (i <= k)
                    {
                        integral = fci.get_two_body_integral(i, k, a, k) - fci.get_two_body_integral(i, k, k, a);
                    }
                    else
                    {
                        integral = fci.get_two_body_integral(k, i, k, a) - fci.get_two_body_integral(k, i, a, k);
                    }
                    
                   entry.two_body_int[k] = integral;
                   max_val = std::max(max_val, fabs(integral));
                }

                // Store the excitation if large enough.
                if (max_val > threshold)
                {
                    single_excitation[i].push_back(entry);
                }
            }
        }
        // Sort
        const auto comp = [](auto& x, auto& y) {return fabs(x.one_body_int) > fabs(y.one_body_int);};
        for (auto& i_column : single_excitation)
            std::sort(i_column.begin(), i_column.end(), comp);

        return;
    }

    // Size of one_body_integral = norb.
    // one_body_integral[i] = <i|f|i>
    // Size of two_body_integral = norb * norb.
    // two_body_integral{{i, j}} = <ij||ij>, i < j.
    // The lower diagonal is not used.
    void construct_diagonal(const Fcidump& fci)
    {
        // Reserve space
        diagonal.one_body_integral.resize(norb);
        diagonal.two_body_integral.resize(norb * norb);

        // One-body integral <i|f|i>, diagonal
        for (auto i = 0; i < norb; ++i)
        {
            diagonal.one_body_integral[i] = fci.get_one_body_integral(i, i);
        }

        // Two-body integral <ij||ij> (same string) and <ij|g|ij> (opposite string)
        for (auto i = 0; i < norb; ++i)
            for (auto j = i + 1; j < norb; ++j)
            {
                diagonal.two_body_integral[index(i, j)] = fci.get_two_body_integral(i, j, i, j)\
                 - fci.get_two_body_integral(i, j, j, i);
            }

        return;
    }

    void get_diagonal(DeterminantDecoded<N>& det, Column& result) const
    {
        auto value = get_diagonal(det);
        result.push_back({det, value});
    }

    value_type get_diagonal(DeterminantDecoded<N>& det) const
    {
        auto occupied_orbitals = det.get_occupied_orbitals();
        return get_diagonal(occupied_orbitals);
    }

    value_type get_diagonal(OrbitalList& occupied_orbitals) const
    {
        value_type value = 0;
        
        for (auto idx_i = 0; idx_i < occupied_orbitals.size(); ++idx_i)
        {
            auto i = occupied_orbitals[idx_i];
            // One body integral
            value += diagonal.one_body_integral[i];
            // Two body integral, i < j
            for (auto idx_j = idx_i + 1; idx_j < occupied_orbitals.size(); ++idx_j)
            {
                auto j = occupied_orbitals[idx_j];
                value += diagonal.two_body_integral[index(i, j)];
            }
        }
        value += core_energy;
        return value;
    }

    Column get_single_excitation(DeterminantDecoded<N>& det) const
    {
        Column result;
        get_single_excitation(det, result);
        return result;
    }

    void get_single_excitation(DeterminantDecoded<N>& det, Column& result) const
    {
        // Single excitation from orbital i
        for (auto i : det.occupied_orbitals)
        {
            // Loop over all possible excitation orbital a
            for (auto& entry_a : single_excitation[i])
            {
                auto a = entry_a.a;
                // The excitation exists if a is not occupied.
                if (!det.is_occupied(a))
                {
                    // <i|f|a>
                    value_type value = entry_a.one_body_int;
                    // <ik||ak>
                    for (auto k : det.occupied_orbitals) value += entry_a.two_body_int[k];
                    // value is the Hamiltonian. Construct the new Determinant.
                    Determinant<N> new_det(det);
                    new_det.clear_orbital(i);
                    // Calculate sign
                    auto sign = new_det.parity(i, a);
                    if (sign) value = -value;

                    new_det.set_orbital(a);
                    result.push_back({new_det, value});
                }
            }
        }
        return;
    }

    Column get_double_excitation(DeterminantDecoded<N>& det) const
    {
        Column result;
        get_double_excitation(det, result);
        return result;
    }

    void get_double_excitation(DeterminantDecoded<N>& det, Column& result) const
    {
        // Double excitation from orbitals (i, j).
        for (auto idx_i = 0; idx_i < det.occupied_orbitals.size(); ++idx_i)
        {
            // auto i = det.occupied_orbitals[idx_i];
            for (auto idx_j = idx_i + 1; idx_j < det.occupied_orbitals.size(); ++idx_j)
            {
                // auto j = det.occupied_orbitals[idx_j];
                get_double_excitation(det, result, idx_i, idx_j);
            }
        }

        return;
    }

    // Get all double excitations from occupied orbital index (idx_i, idx_j).
    void get_double_excitation(DeterminantDecoded<N>& det, Column& result, int idx_i, int idx_j) const
    {
        auto i = det.occupied_orbitals[idx_i];
        auto j = det.occupied_orbitals[idx_j];
        // Loop over all possible orbitals a, b
        for (auto &entry : double_excitation[index(i, j)])
        {
            auto a = entry.first.first;
            auto b = entry.first.second;
            // The excitation exists if a, b are not occupied
            if (!det.is_occupied(a) & !det.is_occupied(b))
            {
                value_type value = entry.second;
                Determinant<N> new_det(det);
                // Calculate sign
                new_det.clear_orbital(i);
                auto sign1 = new_det.parity(i, a);
                new_det.set_orbital(a);
                new_det.clear_orbital(j);
                auto sign2 = new_det.parity(j, b);
                new_det.set_orbital(b);
                auto sign = sign1 ^ sign2;
                if (sign)  value = -value;

                result.push_back({new_det, value});
            }
        }
        return;
    }

    Column get_column(DeterminantDecoded<N>& det) const
    {
        Column result;
        get_diagonal(det, result);
        get_single_excitation(det, result);
        get_double_excitation(det, result);
        return result;
    }

    Column get_column_diagonal_single(DeterminantDecoded<N>& det) const
    {
        Column result;
        get_diagonal(det, result);
        get_single_excitation(det, result);
        return result;
    }

    std::vector<std::pair<Orbital, value_type>> get_sorted_orbital() const
    {
        std::vector<std::pair<Orbital, value_type>> result;
        if (orbital_energy.empty())
        {
            for (auto i = 0; i < norb; ++i)
            {
                result.push_back({i, i});
            }
            std::cout << "Warning: no orbital energy found. Assume orbitals are sorted." << std::endl << std::endl;
        }
        // Sort orbitals by the orbital energy.
        else
        {
            for (auto &entry : orbital_energy)
            {
                result.push_back({entry.first, entry.second});
            }
        }
        const auto comp = [](auto& x, auto& y) {return x.second < y.second;};
        std::sort(result.begin(), result.end(), comp);
        return result;
    }

    Determinant<N> get_hartree_fock() const
    {
        Determinant<N> hf;
        auto sorted_orbital = get_sorted_orbital();

        // Number of occupied alpha orbials and occupied beta orbitals.
        auto n_alpha = (nelec + ms2) / 2;
        auto n_beta = (nelec - ms2) / 2;

        for (auto i_alpha = 0, i_beta = 0, i = 0; i < nelec; ++i)
        {
            auto orb = sorted_orbital[i].first;
            // Alpha orbital
            if ((orb % 2 == 0) && i_alpha < n_alpha)
            {
                hf.set_orbital(orb);
                ++i_alpha;
            }
            // Beta orbital
            else if ((orb % 2 == 1) && i_beta < n_beta)
            {
                hf.set_orbital(orb);
                ++i_beta;
            }
        }
        
        return hf;
    }
};

#endif
