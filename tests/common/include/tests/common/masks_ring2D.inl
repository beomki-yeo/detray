/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <gtest/gtest.h>

#include "detray/definitions/units.hpp"
#include "detray/masks/masks.hpp"

using namespace detray;

constexpr scalar tol{1e-7f};

/// This tests the basic functionality of a ring
TEST(mask, ring2D) {
    using point_t = typename mask<ring2D<>>::loc_point_t;

    point_t p2_pl_in = {0.5f, -2.f};
    point_t p2_pl_edge = {0.f, 3.5f};
    point_t p2_pl_out = {3.6f, 5.f};

    constexpr scalar inner_r{0.f * unit<scalar>::mm};
    constexpr scalar outer_r{3.5f * unit<scalar>::mm};

    mask<ring2D<>> r2{0u, inner_r, outer_r};

    ASSERT_NEAR(r2[ring2D<>::e_inner_r], 0.f, tol);
    ASSERT_NEAR(r2[ring2D<>::e_outer_r], 3.5f, tol);

    ASSERT_TRUE(r2.is_inside(p2_pl_in) == intersection::status::e_inside);
    ASSERT_TRUE(r2.is_inside(p2_pl_edge) == intersection::status::e_inside);
    ASSERT_TRUE(r2.is_inside(p2_pl_out) == intersection::status::e_outside);
    // Move outside point inside using a tolerance
    ASSERT_TRUE(r2.is_inside(p2_pl_out, 1.2f) ==
                intersection::status::e_inside);

    // Check projection matrix
    const auto proj = r2.projection_matrix<e_bound_size>();
    for (unsigned int i = 0u; i < 2u; i++) {
        for (unsigned int j = 0u; j < e_bound_size; j++) {
            if (i == j && i < decltype(r2)::shape::meas_dim) {
                ASSERT_EQ(getter::element(proj, i, j), 1u);
            } else {
                ASSERT_EQ(getter::element(proj, i, j), 0u);
            }
        }
    }

    // Check bounding box
    constexpr scalar envelope{0.01f};
    const auto loc_bounds = r2.local_min_bounds(envelope);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_min_x], -(outer_r + envelope), tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_min_y], -(outer_r + envelope), tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_min_z], -envelope, tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_max_x], (outer_r + envelope), tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_max_y], (outer_r + envelope), tol);
    ASSERT_NEAR(loc_bounds[cuboid3D<>::e_max_z], envelope, tol);
}