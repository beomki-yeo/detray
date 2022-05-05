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

using namespace detray;

TEST(container_cuda, tuple_vector_container) {

    // Vecmem memory resource
    vecmem::cuda::managed_memory_resource resource;

    // Create tuple vector container
    tuple_vector_container<thrust::tuple, vecmem::vector, std::size_t, int,
                           float, double>
        container(resource);
}

TEST(container_cuda, tuple_array_container) {}
