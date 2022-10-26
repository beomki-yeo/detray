/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

namespace detray {

template <typename algebra_t>
struct gain_matrix_smoother {

    using output_type bool;
    using matrix_operator = typename algebra_t::matrix_operator;

    template <typename mask_group_t, typename index_t, typename surface_t,
              typename track_state_t>
    DETRAY_HOST_DEVICE inline output_type operator()(
        const mask_group_t& mask_group, const index_t& /*index*/,
        const surface_t& surface, track_state_t& trk_state) const {

        // Mask associated with the track state
        const auto& mask = mask_group[surface.mask().index()];
    }
};

}  // namespace detray