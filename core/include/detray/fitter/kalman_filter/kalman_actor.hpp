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

template <typename algebra_t, template <typename...> class vector_t>
struct kalman_actor : actor {

    struct state {
        state() {}

        vector_t<track_state> m_track_states;
        std::size_t current_state_id = 0;
    };

    struct kernel {
        using output_type bool;
        using matrix_operator = typename algebra_t::matrix_operator;

        template <typename mask_group_t, typename index_t, typename surface_t,
                  typename track_state_t>
        DETRAY_HOST_DEVICE inline output_type operator()(
            const mask_group_t& mask_group, const index_t& /*index*/,
            const surface_t& surface, track_state_t& trk_state) const {

            // Reference: Application of Kalman filtering to track and vertex
            // fitting, R.Frühwirth, NIM A

            // Mask associated with the track state
            const auto& mask = mask_group[surface.mask().index()];

            // Some identity matrices
            const auto I66 =
                matrix_operator().identity<e_bound_size, e_bound_size>();
            const auto I22 = matrix_operator().identity<2, 2>();

            // projection matrix
            // Dimension: (2 X 6)
            const auto H = mask.projection_matrix();

            // Measurement data on surface
            // Dimension: (2 X 1)
            const auto& meas = trk_state.measurement();

            // Predicted vector of bound track parameters
            // Dimension: (6 X 1)
            const auto& predicted_vec = trk_state.predicted().vector();

            // Predicted covaraince of bound track parameters
            // Dimension: (6 X 6)
            const auto& predicted_cov = trk_state.predicted().covaraince();

            // Spatial resolution (Measurement error)
            // Dimension: (2 X 2)
            const auto V = trk_state.measurement_error();

            // Dimension: (2 X 6) * (6 X 6) * (6 X 2)
            const auto M =
                H * predicted_cov * matrix_operator().transpose(H) + V;

            // Kalman gain matrix
            // Dimension: (6 X 6) * (6 X 2) * (2 X 2)
            const auto K = predicted_cov * matrix_operator().transpose(H) *
                           matrix_operator().inverse(M);

            // Update filtered track parameters with Kalman gain
            auto& filtered_vec = trk_state.filtered().vector();
            auto& filtered_cov = trk_state.filtered().covariance();

            // Dimension: (6 X 1) + (6 X 2) * [ (2 X 1) - (2 X 6) * (6 X 1) ]
            filtered_vec = predicted_vec + K * (meas - H * predicted_vec);
            // Dimension: [ (6 X 6) - (6 X 2) * (2 X 6) ] * (6 X 6)
            filtered_cov = (I66 - K * H) * predicted_cov;

            // Residual between measurement and (projected) filtered vector
            // Dimension: (2 X 1) - (2 X 6) * (6 X 1);
            const auto residual = meas - H * filtered_vec;

            // Calculate the chi square
            auto& chi2 = trk_state.chi2();
            // Dimension: [ (2 X 2) - (2 X 6) * (6 X 2) ] * (2 X 2)
            const auto R = (I22 - H * K) * V;
            // Dimension: (1 X 2) * ( 2 X 2 ) * (2 X 1)
            chi2 = matrix_operator().transpose(residual) *
                   matrix_operator().inverse(R) * residual;

            return true;
        }
    };

    template <typename propagator_state_t>
    DETRAY_HOST_DEVICE void operator()(state& s,
                                       propagator_state_t& propagation) const {

        auto& navigation = propagation._navigation;
        auto& stepping = propagation._stepping;

        // triggered only for sensitive surfaces
        if (navigation.is_on_module() &&
            navigation.current()->sf_id == surface_id::e_sensitive) {

            if (navigation.current_object() !=
                s.m_track_states[s.current_id].surface_id()) {
                propagation._heartbeat &= navigation.abort();
            }

            auto det = navigation.detector();
            const auto& mask_store = det->mask_store();

            // Intersection
            const auto& is = navigation.current();

            // Surface
            const auto& surface = det->surface_by_index(is->index);

            mask_store.template call<kernel>(surface.mask(), surface,
                                             s.m_track_states[s.current_id]);

            s.current_id++;
        }
    };

}  // namespace detray