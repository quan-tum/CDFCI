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

#ifndef CDFCI_CONFIG_H
#define CDFCI_CONFIG_H

/* version */
#define CDFCI_VERSION_MAJOR 0 // for incompatible API changes
#define CDFCI_VERSION_MINOR 1 // for adding functionality in a backwards-compatible manner
#define CDFCI_VERSION_PATCH 0 // for backwards-compatible bug fixes

/* macros */

// Count Trailing Zero
// GCC builtin function
#ifdef __x86_64__
#   define CTZ __builtin_ctzll
#   define PARITY __builtin_parityll
#else
#   define CTZ __builtin_ctz
#   define PARITY __builtin_parity
#endif

// OpenMP
#ifdef _OPENMP
#   include <omp.h>
#endif

// float128 support
// temporary solution for clang: use long double instead
#if defined(__clang__)
#   define CDFCI_USE_LONG_DOUBLE
#   define QUAD_PRECISION long double
#else
#   define QUAD_PRECISION __float128
#endif



/* Global type definition */
using Orbital = int;  // spin-orbitals
using OrbitalPair = std::pair<Orbital, Orbital>;

#endif
