/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/intersection/planar_intersector.hpp"
#include "detray/masks/masks.hpp"
#include "detray/materials/material_slab.hpp"
#include "detray/measurements/geometry_unit.hpp"

namespace detray {

template <typename scalar_t, typename material_t>
class strip final : public geometry_unit<1, scalar_t, material_t> {

    public:
    using base_type = geometry_unit<1, scalar_t, material_t>;
    static constexpr int DIM = base_type::DIM;
    using scalar_type = typename base_type::scalar_type;

    using mask_type = rectangle2<>;
    using material_structure_type = material_slab<material_t>;

    strip(const mask_type& mask, const material_structure_type& slab)
        : m_mask(mask), m_slab(slab) {}

    private:
    mask_type m_mask;
    material_structure_type m_slab;
};

}  // namespace detray