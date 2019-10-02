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

#ifndef CDFCI_SOLVER_H
#define CDFCI_SOLVER_H 1

#include <chrono>
#include <iomanip>
#include "wavefunction.h"
#include "hamiltonian.h"
#include "lib/robin_hood/robin_hood.h"

/* Coordinate Pick functors*/
template<int N = 1>
class CoordinatePick
{
    public:

    virtual typename WaveFunctionVector<N>::value_type operator () (const WaveFunctionVector<N>& det_list) const = 0;

    virtual ~CoordinatePick() {};
};

template<int N = 1>
class CoordinatePickGcdGrad : public CoordinatePick<N>
{
    public:

    ~CoordinatePickGcdGrad() {};

    typename WaveFunctionVector<N>::value_type operator () (const WaveFunctionVector<N>& det_list) const
    {
        typename WaveFunctionVector<N>::value_type result;
        // Warning: assume det_list is non-empty.
        double max_abs_grad = 0;
        for (auto& det : det_list.data)
        {
            auto x = det.second[0];
            auto z = det.second[1];
            double abs_grad = fabs(x * det_list.xx() + z);
            if (abs_grad >= max_abs_grad)
            {
                max_abs_grad = abs_grad;
                result = det;  // FLAG: may be inefficient.
            }
        }
        return result;
    }
};

/* Coordinate Update functors*/
template<int N = 1>
class CoordinateUpdate
{
    public:

    virtual double operator () (Hamiltonian<N>& h, WaveFunction<N>& xz, 
        typename WaveFunctionVector<N>::value_type& det_picked, DeterminantDecoded<N>& det_decoded) const = 0;

    virtual ~CoordinateUpdate() {};
};

template<int N = 1>
class CoordinateUpdateLS : public CoordinateUpdate<N>
{
    public:

    ~CoordinateUpdateLS() {};

    double operator () (Hamiltonian<N>& h, WaveFunction<N>& xz, 
        typename WaveFunctionVector<N>::value_type& det_picked, DeterminantDecoded<N>& det_decoded) const
    {
        double dx = 0;

        auto det = det_picked.first;
        auto x = det_picked.second[0];
        auto z = det_picked.second[1];
        auto xx = xz.xx();

        // Diagonal dA
        auto dA = h.get_diagonal(det_decoded);
        dA = -dA;

        // Line Search
        auto p1 = xx - x * x - dA;
        auto q = z + dA * x;  //
        auto p3 = p1/3;
        auto q2 = q/2;
        auto d = p3 * p3 * p3 + q2 * q2;
        double rt = 0;

        const double pi = atan(1.0) * 4;

        if (d >= 0)
        {
            auto qrtd = sqrt(d);
            rt = cbrt(-q2 + qrtd) + cbrt(-q2 - qrtd);
        }
        else
        {
            auto qrtd = sqrt(-d);
            if (q2 >= 0)
            {
                rt = 2 * sqrt(-p3) * cos((atan2(-qrtd, -q2) - 2*pi)/3.0);
            }
            else
            {
                rt = 2 * sqrt(-p3) * cos(atan2(qrtd, -q2)/3.0);
            }
            
        }

        dx = rt - x;

        // Newton iteration to improve accuracy
        auto dxn = dx - (dx*(dx*(dx + 3*x) + (3*x*x + p1)) + p1*x + q + x*x*x)
             / (dx*(3*dx + 6*x) + 3*x*x + p1);

        const double depsilon = 1e-12;
        const int max_iter = 10;
        int iter = 0;
        while (fabs((dxn - dx)/dx) > depsilon && iter < max_iter)
        {
            dx = dxn;
            dxn = dx - (dx*(dx*(dx + 3*x) + (3*x*x + p1)) + p1*x + q + x*x*x)
                / (dx*(3*dx + 6*x) + 3*x*x + p1);
            ++iter;
        }

        return dx;

    }

};

/* CDFCI solver */
class Solver
{
    public:

    // Store the energy solved.
    double result_energy = 0;

    virtual ~Solver() {};
};

template<int N = 1>
class CDFCI : public Solver
{
    public:

    Option option = {
        {"num_iterations", 10},
        {"report_interval", 1},
        {"coordinate_pick", "gcd_grad"},
        {"coordinate_update", "ls"},
        {"ref_det_occ", Option::array()},
        {"z_threshold", 0.0},
        {"z_threshold_search", false},
        {"max_wavefunction_size", 1000000}
    };

    CDFCI() {};

    CDFCI(Option& opt)
    {
        option.merge_patch(opt);
    }

    ~CDFCI() {};

    int solve(Hamiltonian<N>& h, CoordinatePick<N>& coord_pick, CoordinateUpdate<N>& coord_update, WaveFunction<N>& xz, WaveFunctionVector<N>& sub_xz)
    {
        double energy = 0;

        int num_iter = option["num_iterations"];
        int report_interval = option["report_interval"];
        double z_threshold = option["z_threshold"];
        bool z_threshold_search = option["z_threshold_search"];
        typename Hamiltonian<N>::Column column;
        size_t last_xz_size = 0;

        auto time_start = std::chrono::high_resolution_clock::now();
        double dx = 0;

        std::cout << std::setw(13) << std::left << "Iteration";
        std::cout << std::setw(18) << std::right << "Energy";
        std::cout << std::setw(18) << "dx";
        std::cout << std::setw(15) << "|z|_0";
        std::cout << std::setw(9) << "|H_i|_0";
        std::cout << std::setw(10) << "Time";
        std::cout << std::endl;

        for (auto i = 0; i < num_iter / report_interval; ++i)
        {
            for (auto j = 0; j < report_interval; ++j)
            {
                // Coordinate Pick
                auto det_picked = coord_pick(sub_xz);
                DeterminantDecoded<N> det_decoded(det_picked.first);

                // Coordinate Update
                dx = coord_update(h, xz, det_picked, det_decoded);

                // Update x
                xz.update_x(det_picked.first, dx);

#ifdef _OPENMP
                // Use serial version for diagonal and single excitation
                column = h.get_column_diagonal_single(det_decoded);
                xz.update_z_and_get_sub(column, dx, sub_xz, z_threshold);
                // Recalculate z_i
                double new_z = sub_xz.new_z;
                // Update size
                xz.update_size(sub_xz.n_new_element);

                // Use OpenMP for double excitation
                #pragma omp parallel shared(det_decoded, dx, column, h, sub_xz, xz)
                {
                typename Hamiltonian<N>::Column column_double;
                WaveFunctionVector<N> sub_xz_double;
                // get double excitation
                #pragma omp for
                for (auto idx = 0; idx < h.nelec * (h.nelec - 1) / 2; ++idx)
                {
                    int idx_j = 0.5 + sqrt(0.25 + 2 * idx);
                    int idx_i = idx - (idx_j - 1) * idx_j / 2;
                    h.get_double_excitation(det_decoded, column_double, idx_i, idx_j);
                }
                xz.update_z_only(column_double, dx, sub_xz_double, z_threshold);

                // Reduce sub_xz_double to sub_xz and xz
                #pragma omp critical
                {
                    // Update xz
                    xz.update_xz(sub_xz_double.xz());
                    xz.update_size(sub_xz_double.n_new_element);  // update size

                    // Reduce recalculated z_i
                    new_z += sub_xz_double.new_z;
                    // Update sub_xz
                    sub_xz.append(sub_xz_double);

                }
                #pragma omp barrier
                #pragma omp master
                {
                    sub_xz.set_xx(xz.xx());
                    sub_xz.set_xz(xz.xz());
                    // Update recalculated z_i
                    xz.reinsert_z(det_picked.first, new_z);
                }

                }  // end omp parallel
#else
                // Get H(:, det_picked)
                column = h.get_column(det_decoded);

                // Update z
                xz.update_z_and_get_sub(column, dx, sub_xz, z_threshold);
#endif
            }
            // xz.refresh();
            energy = xz.get_variational_energy();
            auto xz_size = xz.size();

            auto time_end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> time_elpse = time_end - time_start;

            // Output
            std::cout << std::setw(13) << std::left << (i+1)*report_interval;
            std::cout << std::setw(18) << std::right << std::fixed << std::setprecision(10) << energy;
            std::cout << std::setw(18) << std::scientific <<std::setprecision(4) << dx;
            std::cout << std::setw(15) << xz_size;
            std::cout << std::setw(9) << sub_xz.size();
            std::cout << std::setw(10) << std::fixed << std::setprecision(2) << time_elpse.count();
            std::cout << std::endl;

            // Auto adjust z_threshold
            if (z_threshold_search)
            {
                auto inc_ratio = 1000.0 * (xz_size - last_xz_size) / sub_xz.size() / report_interval;
                if ((inc_ratio >= 1) && (xz_size> 0.7 * xz.max_size))
                {
                    return -1;
                }
            }
            last_xz_size = xz.size();

            // Check overflow. May be unnecessary but cheap.
            if (xz_size > 0.79 * xz.max_size)
            {
                throw std::overflow_error("The hash table is full. Please increase max_memory or z_threshold.");
            }
        }

        // Store the final computed energy in the Solver class.
        result_energy = energy;
        return 0;
    }

    int solve(Hamiltonian<N>& h, CoordinatePick<N>& coord_pick, CoordinateUpdate<N>& coord_update)
    {
        bool z_threshold_search = option["z_threshold_search"];
        double z_threshold = option["z_threshold"];

        // Print information
        std::cout << "CDFCI calculation" << std::endl;
        std::cout << "-----------------" << std::endl;
        // Find the initial value.
        auto hf = get_ref_det(h);
        DeterminantDecoded<N> hf_decoded(hf);
        std::cout << "Hartree Fock determinant occupied spin-orbitals:" << std::endl;
        for (auto& orb : hf_decoded.occupied_orbitals)
        {
            std::cout << orb << "  ";
        }
        std::cout << std::endl;
        std::cout << "Hartree Fock energy: " << std::fixed << std::setprecision(10) << \
             h.get_diagonal(hf_decoded) << std::endl << std::endl;
        
        while (true)
        {
            // Initialize xz
#ifdef _OPENMP
            WaveFunctionCuckoo<N> xz(option["max_wavefunction_size"]);
#else
            WaveFunctionStd<N> xz(option["max_wavefunction_size"]);
#endif
            // Initialize sub_xz
            WaveFunctionVector<N> sub_xz;
            std::array<double, 2> val{0, 0};
            sub_xz.data.push_back({hf, val});

            auto ierr = solve(h, coord_pick, coord_update, xz, sub_xz);
            if (z_threshold_search && (ierr == -1))
            {
                z_threshold = z_threshold * 10;
                option["z_threshold"] = z_threshold;
                std::cout << "z_threshold_search: z_threshold is too small. Increase z_threhold to ";
                std::cout << std::scientific << std::setprecision(2) << option["z_threshold"] << " and restart." << std::endl;
            }
            else
            {
                return ierr;
            }
            
        }
    }

    Determinant<N> get_ref_det(Hamiltonian<N>& h) const
    {
        std::vector<Orbital> occ_list = option["ref_det_occ"];
        Determinant<N> hf;
        if (occ_list.empty())
        {
            std::cout << "Reference determinant is not provided. Use Hartree Fock state from FCIDUMP." << std::endl;
            return h.get_hartree_fock();
        }
        else
        {
            for (auto orb : occ_list)
            {
                if ((orb >= 0) && (orb < h.norb))
                {
                    hf.set_orbital(orb);
                }
            }
            if (hf.get_occupied_orbitals().size() != h.nelec)
            {
                throw std::invalid_argument("Reference determinant is invalid. The number of valid orbitals does not match "
                      "the number of electrons. Please provide " + std::to_string(h.nelec) + " orbitals from 0 to "
                       + std::to_string(h.norb - 1) + " or remove the argument.");
            }
        }
        return hf;
    }
};

#endif
