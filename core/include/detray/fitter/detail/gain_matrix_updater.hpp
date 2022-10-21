/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).

namespace detray {

template <typename algebra_t>
struct gain_matrix_updater {
    using matrix_operator = typename algebra_t::matrix_operator;

    template <typename track_state_t, typename propagation_state_t>
    bool operator()(track_state_t& trk_state,
                    const propagation_state_t& propagation) const {


    }
};

}  // namespace detray