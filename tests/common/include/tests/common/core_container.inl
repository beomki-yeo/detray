/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Detray include(s)
#include "detray/core/detail/tuple_array_container.hpp"
#include "detray/core/detail/tuple_vector_container.hpp"

// Vecmem include(s)
#include <vecmem/containers/vector.hpp>
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s)
#include <gtest/gtest.h>

// System include(s)
#include <array>
#include <tuple>
#include <vector>

using namespace detray;

TEST(container, tuple_vector_container) {

    // Vecmem memory resource
    vecmem::host_memory_resource resource;

    // Create tuple vector container
    tuple_vector_container<std::tuple, vecmem::vector, std::size_t, int, float,
                           double>
        container(resource);

    container.add_elements<0>(1);
    container.add_elements<0>(2);
    container.add_elements<1>(3.1);
    container.add_elements<1>(4.5);
    container.add_elements<2>(5.5);
    container.add_elements<2>(6.);
}

TEST(container, tuple_array_container) {}