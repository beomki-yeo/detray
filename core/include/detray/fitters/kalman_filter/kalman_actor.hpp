/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/qualifiers.hpp"
#include "detray/fitters/kalman_filter/gain_matrix_updater.hpp"
#include "detray/fitters/kalman_filter/track_state.hpp"
#include "detray/propagator/base_actor.hpp"

namespace detray {

template <typename algebra_t, typename measurement_t,
          template <typename...> class vector_t>
struct kalman_actor : actor {

    using track_state_type = track_state<algebra_t, measurement_t>;

    struct state {
        DETRAY_HOST_DEVICE
        state(vector_t<track_state_type>&& track_states)
            : m_track_states(std::move(track_states)) {
            m_forward_it = m_track_states.begin();
            m_backward_it = m_track_states.rbegin();
        }

        DETRAY_HOST_DEVICE
        track_state_type& operator()(const navigation::direction dir) {
            if (dir == navigation::direction::e_forward) {
                return *m_forward_it;
            } else {
                return *m_backward_it;
            }
        }

        DETRAY_HOST_DEVICE
        void next(const navigation::direction dir) {
            if (dir == navigation::direction::e_forward) {
                m_forward_it++;
            } else {
                m_backward_it++;
            }
        }

        DETRAY_HOST_DEVICE
        bool is_complete(const navigation::direction dir) const {
            if (dir == navigation::direction::e_forward &&
                m_forward_it == m_track_states.end()) {
                return true;
            } else if (dir == navigation::direction::e_backward &&
                       m_backward_it == m_track_states.rend()) {
                return true;
            }
            return false;
        }

        vector_t<track_state_type> m_track_states;
        // iterator for forward filtering
        typename vector_t<track_state_type>::iterator m_forward_it;
        // iterator for backward filtering
        typename vector_t<track_state_type>::reverse_iterator m_backward_it;
    };

    template <typename propagator_state_t>
    DETRAY_HOST_DEVICE void operator()(state& actor_state,
                                       propagator_state_t& propagation) const {

        auto& navigation = propagation._navigation;
        auto& stepping = propagation._stepping;

        // If the iterator reaches the end, terminate the propagation
        if (actor_state.is_complete()) {
            propagation._heartbeat &= navigation.abort();
        }

        // triggered only for sensitive surfaces
        if (navigation.is_on_module() &&
            navigation.current()->sf_id == surface_id::e_sensitive) {

            auto& trk_state = actor_state(navigation.direction());

            // Abort if the propagator fails to find the next measurement
            if (navigation.current_object() != trk_state.surface_link()) {
                propagation._heartbeat &= navigation.abort();
            }

            auto det = navigation.detector();
            const auto& mask_store = det->mask_store();

            // Surface
            const auto& surface =
                det->surface_by_index(trk_state.surface_link());

            // Run kalman updater
            mask_store.template call<gain_matrix_updater<algebra_t>>(
                surface.mask(), surface, trk_state);

            // Update iterator
            actor_state.next(navigation.direction());
        }
    }
};

}  // namespace detray