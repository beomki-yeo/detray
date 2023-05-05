/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

namespace detray {

struct mask_property {

    enum class mask_type {
        e_sensitive = 0,
        e_portal = 1,
        e_passive = 2,
    };
};

struct portal : mask_property {
    inline static constexpr const mask_type P = mask_type::e_portal;
};

struct passive : mask_property {
    inline static constexpr const mask_type P = mask_type::e_passive;
};

template <unsigned int kMeasDim, bool kNormalOrdering>
struct sensitive : mask_property {
    inline static constexpr const mask_type P = mask_type::e_sensitive;
    inline static constexpr const unsigned int meas_dim{kMeasDim};
    inline static constexpr const bool normal_ordering{kNormalOrdering};
};

}  // namespace detray