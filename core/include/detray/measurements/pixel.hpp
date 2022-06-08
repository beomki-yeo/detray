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
#include "detray/materials/homogenous_surface_material.hpp"
#include "detray/measurement/tuple_container.hpp"

namespace detray {

class pixel final
    : public measurement_unit<2, rectangle2,
                              homogenous_surface_material<material_slab>> {

    pixel() {}
};

}  // namespace detray