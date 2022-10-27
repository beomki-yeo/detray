/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/tracks/bound_track_parameters.hpp"

namespace detray {

/// Track state definition for kalman filtering
template <typename algebra_t>
struct track_state{

    DETRAY_HOST_DEVICE 
    inline bound_track_parameters<algebra_t> filtered() const {
        return m_filtered;
    }

    DETRAY_HOST_DEVICE 
    inline bound_track_parameters<algebra_t> smoothed() const {
        return m_smoothed;
    }

    bound_track_parameters<algebra_t> m_filtered;
    bound_track_parameters<algebra_t> m_smoothed;
};

}  // namespace detray