/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

namespace detray {

template <int D, typename scalar_t, typename... material_t>
class geometry_unit {

    static constexpr int DIM = D;
    using scalar_type = scalar_t;

    public:
};

}  // namespace detray