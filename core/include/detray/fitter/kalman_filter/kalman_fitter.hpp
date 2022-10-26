/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/qualifiers.hpp"
#include "detray/propagator/base_actor.hpp"

namespace detray {

template <typename algebra_t>
struct kalman_fitter : actor {

    struct state {};

    struct kernel {};

    template <typename propagator_state_t>
    DETRAY_HOST_DEVICE void operator()(state& /*resetter_state*/,
                                       propagator_state_t& propagation) const {}
};

}  // namespace detray