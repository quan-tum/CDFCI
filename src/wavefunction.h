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

#ifndef CDFCI_WAVEFUNCTION_H
#define CDFCI_WAVEFUNCTION_H 1

#include <unordered_map>
#include <cmath>
#include "lib/libcuckoo/cuckoohash_map.hh"
#include "lib/robin_hood/robin_hood.h"
#include "determinant.h"

template<int N = 1>
class WaveFunctionVector
{
    public:

    using value_type = std::pair<Determinant<N>, std::array<double, 2>>;
    std::vector<value_type> data;  // TODO: use iterator, and move to private.

    std::array<QUAD_PRECISION, 2> norm_square_ = {0, 0};
    QUAD_PRECISION dot_product_ = 0;
    // Temporary data when updating xz.
    size_t n_new_element = 0;  // Number of elements not in xz. Used to update xz.size() for Cuckoo hash table.
    double new_z = 0;  // Store the recalculated z_i.

    public:

    WaveFunctionVector(long capacity = 0)
    {
        data.reserve(capacity);
    }

    ~WaveFunctionVector() {};

    double xx() const
    {
        return norm_square_[0];
    }

    void set_xx(double xx)
    {
        norm_square_[0] = xx;
    }

    double xz() const
    {
        return dot_product_;
    }

    void set_xz(double xz)
    {
        dot_product_ = xz;
    }

    size_t size() const
    {
        return data.size();
    }

    void clear()
    {
        data.clear();
        norm_square_ = {0, 0};
        dot_product_ = 0;
    }

    void append(WaveFunctionVector<N>& sub_xz)
    {
        data.insert(data.end(), sub_xz.data.begin(), sub_xz.data.end());
    }

};

template<int N = 1>
class WaveFunction
{
    protected:

    std::array<QUAD_PRECISION, 2> norm_square_ = {0, 0};
    QUAD_PRECISION dot_product_ = 0;

    public:

    using key_type = Determinant<N>;
    using value_type = std::array<double, 2>;

    size_t max_size = 0;

    double xx() const
    {
        return norm_square_[0];
    }

    double xz() const
    {
        return dot_product_;
    }

    void update_xz(QUAD_PRECISION inc)
    {
        dot_product_ += inc;
    }

    double get_variational_energy() const
    {
        // xz/xx
        return dot_product_ / norm_square_[0];
    }

    virtual void update_x(key_type& key, double dx) {};
    virtual void update_z_and_get_sub(std::vector<std::pair<key_type, double>>& column, 
        double dx, WaveFunctionVector<N>& sub_xz, double z_threshold) {};
    virtual void update_z_only(std::vector<std::pair<key_type, double>>& column, 
        double dx, WaveFunctionVector<N>& sub_xz, double z_threshold) = 0;
    virtual size_t size() = 0;
    virtual size_t update_size(size_t n_new_element) {return 0;}
    virtual void reinsert_z(key_type& key, double new_z) {};

    virtual ~WaveFunction() {};
};

/* Hash function mapping determinants to size_t */
// Warning: The following hash functions are naive and just okay to use. They should be replaced with
// better hash functions.
// The determinants are represented as an array of size_t integers.
template<int N>
struct DeterminantHash {
    size_t operator() (const Determinant<N>& det) const {
        throw std::invalid_argument("The DeterminantHash does not support N > 4. Please modify the "
            "\"DeterminantHash\" class in \"wavefunction.h\" and provide a DeterminantHash.");
        return 0;
    }
};

template<>
size_t DeterminantHash<1>::operator() (const Determinant<1>& det) const
{
    return det.repr[0];
}

template<>
size_t DeterminantHash<2>::operator() (const Determinant<2>& det) const
{
    return det.repr[0] * 2038076783 + det.repr[1] * 179426549;
}

template<>
size_t DeterminantHash<3>::operator() (const Determinant<3>& det) const
{
    return det.repr[0] * 2038076783 + det.repr[1] * 179426549 + det.repr[2] * 500002577;
}

template<>
size_t DeterminantHash<4>::operator() (const Determinant<4>& det) const
{
    return det.repr[0] * 2038076783 + det.repr[1] * 179426549 + det.repr[2] * 500002577 + det.repr[3] * 255477023;
}

template <int N>
struct DeterminantHashRobinhood
{
    size_t operator()(const Determinant<N> &det) const
    {
        DeterminantHash<N> hash1;
        return (robin_hood::hash_int(hash1(det)));
    }
};

template<int N>
struct DeterminantEqual {
    bool operator () (const Determinant<N>& det1, const Determinant<N>& det2) const {
        bool result = true;
        for (auto i = 0; i < N; ++i)
        {
            result = result && (det1.repr[i] == det2.repr[i]);
        }
        return result;
    }
};

template<int N = 1>
class WaveFunctionStd : public WaveFunction<N>
{
    using typename WaveFunction<N>::key_type;
    using typename WaveFunction<N>::value_type;

    using WaveFunction<N>::norm_square_;
    using WaveFunction<N>::dot_product_;
    using WaveFunction<N>::max_size;

    private:

    robin_hood::unordered_flat_map<key_type, value_type, DeterminantHash<N>, DeterminantEqual<N>> data_;

    public:

    WaveFunctionStd(size_t capacity)
    {
        data_.reserve(capacity);
        max_size = capacity;
    }

    ~WaveFunctionStd() {};

    void update_x(key_type& key, double dx)
    {
        auto iter = data_.find(key);
        if (iter != data_.end())
        {
            // Update xx and xz
            norm_square_[0] += 2 * iter->second[0] * dx + dx * dx;
            dot_product_ += iter->second[1] * dx;
            // Update x
            iter->second[0] += dx;
        }
        else
        {
            value_type val {dx, 0};
            data_.insert({key, val});
            // Update xx and xz
            norm_square_[0] += dx * dx;
        }
        return;
    }

    void update_z_only(std::vector<std::pair<key_type, double>>& column, 
        double dx, WaveFunctionVector<N>& sub_xz, double z_threshold = 0.0)
    {
        // Used to recalculate z_i exactly. Assume the first entry of column is the one to be updated.
        // And column is the full corresponding column of the Hamiltonian.
        double new_z = 0.0;

        // Loop over column, update z, xz and get sub_xz
        for (auto& entry : column)
        {
            auto& det = entry.first;
            auto h = entry.second;
            auto dz = dx * h;

            auto iter = data_.find(det);
            if (iter != data_.end())
            {
                // Update z
                iter->second[1] += dz;
                // Update xz. Note: zz is not used and updated.
                sub_xz.dot_product_ += dz * iter->second[0];
                // Update sub_xz
                value_type val {iter->second[0], iter->second[1]};
                sub_xz.data.push_back({det, val});
                // Recaculate new_z
                new_z += h * iter->second[0];
            }
            else
            {
                // Only update if above threold
                if (fabs(dz) > z_threshold)
                {
                    value_type val {0, dz};
                    data_.insert({det, val});
                    // Update xz. Note d(xz) = 0;
                    // Update sub_xz
                    sub_xz.data.push_back({det, val});
                }
            }
        }
        // Update z_i. Assume column[0] is the coordinate updated.
        auto iter = data_.find(column[0].first);
        iter->second[1] = new_z;
        
        return;
    }

    void update_z_and_get_sub(std::vector<std::pair<key_type, double>>& column, 
        double dx, WaveFunctionVector<N>& sub_xz, double z_threshold = 0.0)
    {
        // Clear sub_xz
        sub_xz.clear();
        
        update_z_only(column, dx, sub_xz, z_threshold);

        // Update dot product. (atomic update)
        dot_product_ += sub_xz.dot_product_;
        // Update sub_xz
        sub_xz.set_xx(norm_square_[0]);
        sub_xz.set_xz(dot_product_); 
        return;
    }

    size_t size()
    {
        return data_.size();
    }
};

template<int N = 1>
class WaveFunctionCuckoo : public WaveFunction<N>
{
    using typename WaveFunction<N>::key_type;
    using typename WaveFunction<N>::value_type;

    using WaveFunction<N>::norm_square_;
    using WaveFunction<N>::dot_product_;
    using WaveFunction<N>::max_size;
    
    private:

    cuckoohash_map<key_type, value_type, DeterminantHashRobinhood<N>, DeterminantEqual<N>,\
         std::allocator<std::pair<const key_type, value_type>>, 8> data_;
    size_t size_ = 0;

    public:

    WaveFunctionCuckoo(size_t capacity)
    {
        data_.reserve(capacity);
        max_size = capacity;
    }

    ~WaveFunctionCuckoo() {};

    void update_x(key_type& key, double dx)
    {
        value_type xz {0, 0};
        value_type val_new {dx, 0};
        auto update_fn = [dx, &xz](value_type& val){xz[0] = val[0]; xz[1] = val[1];
            val[0] += dx; return false;};
        data_.upsert(key, update_fn, val_new);

        norm_square_[0] += 2 * xz[0] * dx + dx * dx;
        dot_product_ += xz[1] * dx;
        return;
    }

    void update_z_only(std::vector<std::pair<key_type, double>>& column, 
		       double dx, WaveFunctionVector<N>& sub_xz, double z_threshold = 0.0)
    {
        // Recalculate z_i
        sub_xz.new_z = 0;
        // Update size
        sub_xz.n_new_element = 0;

        // Loop over column, update z, xz and get sub_xz
        for (auto& entry : column)
        {
            auto& det = entry.first;
            auto h = entry.second;
            auto dz = dx * h;

            value_type val_new{0, 0};
            auto update_functor = [dz, &val_new](value_type &val) {
                // Update z
                val[1] += dz;
                // Get the value after update
                val_new[0] = val[0];
                val_new[1] = val[1];
                return false;  // always return false
            };
            // Call update_functor if det is found
            auto flag_found = data_.update_fn(det, update_functor);
            if (flag_found)
            {
                // Update xz. Note: zz is not used and updated.
                sub_xz.dot_product_ += dz * val_new[0];
                // Update sub_xz
                sub_xz.data.push_back({det, val_new});
                // Recaculate z_i
                sub_xz.new_z += h * val_new[0];
            }
            else
            {
                // Only update if above threold
                if (fabs(dz) > z_threshold)
                {
                    value_type val{0, dz};
                    data_.insert(det, val);
                    // Update xz. Note d(xz) = 0;
                    // Update sub_xz
                    sub_xz.data.push_back({det, val});
                    // Update size
                    sub_xz.n_new_element += 1;
                }
            }
        }
        // Check overflow
        if (size() > 0.79 * max_size)
        {
            throw std::overflow_error("The hash table is full. Please increase max_wavefunction_size or z_threshold.");
        }
        return;
    }

    void update_z_and_get_sub(std::vector<std::pair<key_type, double>>& column, 
			      double dx, WaveFunctionVector<N>& sub_xz, double z_threshold = 0.0)
    {
        // Clear sub_xz
        sub_xz.clear();

        update_z_only(column, dx, sub_xz, z_threshold);
        
        dot_product_ += sub_xz.dot_product_;
        // Update sub_xz
        sub_xz.set_xx(norm_square_[0]);
        sub_xz.set_xz(dot_product_); 
        return;
    }

    void reinsert_z(key_type& key, double new_z)
    {
        value_type val_new {0, new_z};

        auto update_fn = [new_z](value_type& val){
            val[1] = new_z;
            return false;
        };
        data_.upsert(key, update_fn, val_new);

        return;
    }

    size_t size()
    {
        return size_;
    }

    size_t update_size(size_t n_new_element)
    {
        size_ += n_new_element;
        return size_;
    }
};

#endif
