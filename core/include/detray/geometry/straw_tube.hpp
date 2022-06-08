/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/geometry/geometry_unit.hpp"
#include "detray/intersection/planar_intersector.hpp"
#include "detray/masks/masks.hpp"
#include "detray/materials/material_rod.hpp"

namespace detray {

template <typename scalar_t, typename wire_mat_t, typename gas_mat_t,
          typename cylinder_mat_t>
class straw_tube final
    : public geometry_unit<1, scalar_t, wire_mat_t, gas_mat_t, cylinder_mat_t> {

    public:
    using base_type =
        geometry_unit<1, scalar_t, wire_mat_t, gas_mat_t, cylinder_mat_t>;
    static constexpr int DIM = base_type::DIM;
    using scalar_type = typename base_type::scalar_type;

    //// Left empty for the moment
};

}  // namespace detray