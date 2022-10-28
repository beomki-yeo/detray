/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/propagator/navigator.hpp"
#include "detray/tracks/bound_track_parameters.hpp"

namespace detray {

/// Track state definition for kalman filtering
template <typename algebra_t, typename measurement_t>
struct track_state {

    using matrix_operator = typename algebra_t::matrix_operator;
    using size_type = typename matrix_operator::size_ty;
    template <size_type ROWS, size_type COLS>
    using matrix_type =
        typename matrix_operator::template matrix_type<ROWS, COLS>;

    DETRAY_HOST_DEVICE
    track_state(const dindex& surface_link, const measurement_t& meas)
        : m_measurement(meas), m_surface_link(surface_link) {}

    DETRAY_HOST_DEVICE
    inline dindex surface_link() const { return m_surface_link; }

    // Get the local position of measurement
    // FIXME: There is an inefficient conversion from vector to matrix
    DETRAY_HOST_DEVICE
    inline matrix_type<2, 1> measurement_local() const {
        matrix_type<2, 1> ret = matrix_operator().template zero<2, 1>();
        matrix_operator().element(ret, 0, 0) = m_measurement.local[0];
        matrix_operator().element(ret, 1, 0) = m_measurement.local[1];
        return ret;
    }

    // Get the local covariance of measurement
    DETRAY_HOST_DEVICE
    inline matrix_type<2, 2> measurement_covariance() const {
        matrix_type<2, 2> ret = matrix_operator().template zero<2, 2>();
        matrix_operator().element(ret, 0, 0) = m_measurement.variance[0];
        matrix_operator().element(ret, 1, 1) = m_measurement.variance[1];
        return ret;
    }

    DETRAY_HOST_DEVICE inline bound_track_parameters<algebra_t>& filtered(
        const navigation::direction dir) {
        return dir = navigation::direction::e_forward ? m_forward_filtered
                                                      : m_backward_filtered;
    }

    DETRAY_HOST_DEVICE
    inline const bound_track_parameters<algebra_t>& filtered(
        const navigation::direction dir) const {
        return dir = navigation::direction::e_forward ? m_forward_filtered
                                                      : m_backward_filtered;
    }

    DETRAY_HOST_DEVICE
    inline bound_track_parameters<algebra_t>& smoothed() { return m_smoothed; }

    DETRAY_HOST_DEVICE
    inline const bound_track_parameters<algebra_t>& smoothed() const {
        return m_smoothed;
    }

    dindex m_surface_link;
    measurement_t m_measurement;
    bound_track_parameters<algebra_t> m_forward_filtered;
    bound_track_parameters<algebra_t> m_backward_filtered;
    bound_track_parameters<algebra_t> m_smoothed;
};

}  // namespace detray