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

#ifndef CDFCI_DETERMINANT_H
#define CDFCI_DETERMINANT_H 1

#include <vector>
#include <climits>
#include "config.h"

template<int N = 1>
class Determinant
{
    public:

        using det_type = size_t;  // Use of size_t is consistent with the macro definition above.
        
        constexpr static int size = N;
        constexpr static int det_size = sizeof(det_type) * CHAR_BIT;

        // Masks for parity calculation.
        static std::array<std::array<det_type, size>, det_size * size> masks;

        // Bit representation of the Slater determinant, using N numbers of det_type integers.
        std::array<det_type, size> repr;
        
        Determinant()
        {
            // Initialize repr to 0.
            std::fill(repr.begin(), repr.end(), 0);
        }
        Determinant(std::array<det_type, N> input_repr) : repr{input_repr} {}
        Determinant(const Determinant<N>& det)
        {
            repr = det.repr;
        }

        // Need to be called explicitly to construct the masks.
        static void constuct_masks()
        {
            for (auto i = 0; i < masks.size(); ++i)
            {
                masks[i][i / det_size] = (static_cast<det_type>(1) << (i % det_size)) - 1;
                for (auto j = 0; j < i / det_size; ++j)
                {
                    masks[i][j] = ~static_cast<det_type>(0);  // all bits are 1.
                }
                for (auto j = i / det_size + 1; j < size; ++j)
                {
                    masks[i][j] = 0;
                }
            }
            return;
        }

        // FLAG: efficiency needed.
        // Return the occupied orbitals in a vector, from the smallest to the largest one.
        std::vector<Orbital> get_occupied_orbitals() const
        {
            std::vector<Orbital> orbitals;
            for (auto i = 0; i < size; ++i)
            {
                det_type det = repr[i];
                while (det)
                {
                    int orb = CTZ(det);  // the index of the least significant bit
                    orbitals.push_back(orb + det_size * i);
                    det &= (det - 1);  // unset the least significant bit
                }
            }
            return orbitals;
        }

        // FLAG: efficiency needed.
        bool is_occupied(Orbital i) const
        {
            return repr[i / det_size] & (static_cast<det_type>(1) << (i % det_size));
        }

        void set_orbital(Orbital i)
        {
            repr[i / det_size] |= (static_cast<det_type>(1) << (i % det_size));
        }

        void clear_orbital(Orbital i)
        {
            repr[i / det_size] &= ~(static_cast<det_type>(1) << (i % det_size));
        }

        // Return the number of 1's mod 2 in bits min(i, a) <= k < max(i, a).
        int parity(Orbital i, Orbital a) const
        {
            int result = 0;
            for (auto t = 0; t < size; ++t)
            {
                result ^= PARITY(repr[t] & (masks[i][t] ^ masks[a][t]));
            }
            return result;
        }
};

template<>
int Determinant<1>::parity(Orbital i, Orbital a) const
{
    return PARITY(repr[0] & (((static_cast<det_type>(1) << i) - 1) ^ ((static_cast<det_type>(1) << a) - 1)));
}

template<int N>
std::array<std::array<typename Determinant<N>::det_type, Determinant<N>::size>, Determinant<N>::det_size * Determinant<N>::size> Determinant<N>::masks;

template<int N = 1>
class DeterminantDecoded : public Determinant<N>
{
    public:

        std::vector<Orbital> occupied_orbitals;

        DeterminantDecoded(Determinant<N>& det)
        {
            this->repr = det.repr;
            occupied_orbitals = det.get_occupied_orbitals();
        }
};

#endif
