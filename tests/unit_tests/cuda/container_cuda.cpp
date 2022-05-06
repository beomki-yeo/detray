/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Detray include(s)
#include "container_cuda_kernel.hpp"

// Vecmem include(s)
#include <vecmem/memory/cuda/managed_memory_resource.hpp>

// Thrust include(s)
#include <thrust/tuple.h>

// GTest include(s)
#include <gtest/gtest.h>

// System include(s)
#include <numeric>

using namespace detray;

TEST(container_cuda, tuple_vector_container) {

    // Vecmem memory resource
    vecmem::cuda::managed_memory_resource resource;

    // Create tuple vector container
    tuple_vector_container_type container(resource);

    // Base container function check
    EXPECT_EQ(container.size(), 3);
    EXPECT_EQ(container.empty<0>(), true);
    EXPECT_EQ(container.empty<1>(), true);
    EXPECT_EQ(container.empty<2>(), true);

    EXPECT_EQ(container.to_id<>(0), 0);
    EXPECT_EQ(container.to_id<1>(2), 2);
    EXPECT_EQ(container.to_id<1>(0), 3);

    // Add elements to the container
    container.add_element<0>(1);
    container.add_element<0>(2);
    container.add_element<1>(3.1);
    container.add_element<1>(4.5);
    container.add_element<2>(5.5);
    container.add_element<2>(6.);

    vecmem::vector<int> int_vec{3, 4, 5};
    container.add_vector(int_vec);

    vecmem::vector<float> float_vec{12.1, 5.6};
    container.add_vector(float_vec);

    container.add_vector(vecmem::vector<double>{10.5, 7.6});

    // CPU sum check
    double cpu_sum = 0;
    cpu_sum = std::accumulate(container.group<0>().begin(),
                              container.group<0>().end(), cpu_sum);
    cpu_sum = std::accumulate(container.group<1>().begin(),
                              container.group<1>().end(), cpu_sum);
    cpu_sum = std::accumulate(container.group<2>().begin(),
                              container.group<2>().end(), cpu_sum);
    EXPECT_FLOAT_EQ(cpu_sum, 69.9);

    // CUDA sum check
    auto container_data = get_data(container);

    vecmem::vector<double> cuda_sum(&resource);
    cuda_sum.push_back(0);
    auto sum_data = vecmem::get_data(cuda_sum);

    get_sum(container_data, sum_data);

    EXPECT_FLOAT_EQ(cpu_sum, cuda_sum[0]);
}

TEST(container_cuda, tuple_array_container) {

    // Vecmem memory resource
    vecmem::cuda::managed_memory_resource resource;

    // Create tuple array container
    tuple_array_container_type container(resource);

    // Populate the elements
    auto& gr0 = container.group<0>();
    gr0 = std::array<vecmem::vector<int>, 1>{
        vecmem::vector<int>{{3, 4, 5}, &resource}};

    auto& gr1 = container.group<1>();
    gr1 = std::array<vecmem::vector<float>, 2>{
        vecmem::vector<float>{{1.1, 4.2}, &resource},
        vecmem::vector<float>{{5.1, 2.3, 1.7}, &resource}};

    // CPU sum check
    double cpu_sum = 0;

    cpu_sum = std::accumulate(gr0[0].begin(), gr0[0].end(), cpu_sum);
    cpu_sum = std::accumulate(gr1[0].begin(), gr1[0].end(), cpu_sum);
    cpu_sum = std::accumulate(gr1[1].begin(), gr1[1].end(), cpu_sum);
    EXPECT_FLOAT_EQ(cpu_sum, 26.4);

    // CUDA sum check
    auto container_data = get_data(container, resource);

    vecmem::vector<double> cuda_sum(&resource);
    cuda_sum.push_back(0);
    auto sum_data = vecmem::get_data(cuda_sum);

    get_sum(container_data, sum_data);

    EXPECT_FLOAT_EQ(cpu_sum, cuda_sum[0]);
}
