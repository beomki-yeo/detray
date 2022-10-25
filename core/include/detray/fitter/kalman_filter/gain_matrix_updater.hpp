/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

namespace detray {

template <typename algebra_t>
struct gain_matrix_updater {
    using matrix_operator = typename algebra_t::matrix_operator;
    using scalar_type = typename algebra_t::size_type;

    struct functor {
        output_type bool;

        template <typename mask_group_t, typename index_t,
                  typename track_state_t, typename propagator_state_t>
        DETRAY_HOST_DEVICE inline output_type operator()(
            const mask_group_t& mask_group, const index_t& /*index*/,
            track_state_t& trk_state, const propagation_state_t& propagation) {

            // projection matrix
            const auto H = mask.projection_matrix();

            // Measurement data on surface
            const auto& meas = trk_state.measurement();

            // Predicted vector of bound track parameters
            const auto& predicted_vec = trk_state.predicted().vector();

            // Predicted covaraince of bound track parameters
            const auto& predicted_cov = trk_state.predicted().covaraince();

            // Spatial resolution (Measurement error)
            // Dimension: (Meas_Dim X Meas_Dim)
            const auto R = trk_state.measurement_error();

            // Mask for track state
            const auto& mask = mask_group[surface.mask().index()];

            // Dimension: (Meas_Dim X 6) * (6 X 6) * (6 X Meas_Dim)
            const auto M =
                H * predicted_cov * matrix_operator().transpose(H) + R;

            // Kalman gain matrix
            // Dimension: (6 X 6) * (6 X Meas_Dim) * ()
            const auto K = predicted_cov * matrix_operator().transpose(H) *
                           matrix_operator().inverse(M);

            // Update filtered track parameters with Kalman gain
            auto& filtered_vec = trk_state.filtered().vector();
            auto& filtered_cov = trk_state.filtered().covariance();

            filtered_vec = predicted_vec + K * (meas - H * predicted_vec);
            filtered_cov =
                (matrix_operator().identity<>() - K * H) * predicted_cov;
            /*
            filtered = predicted + K * (calibrated - H * predicted);
            filteredCovariance =
                (BoundSymMatrix::Identity() - K * H) * predictedCovariance;
            */
        }
    };

    template <typename track_state_t, typename propagation_state_t>
    bool operator()(track_state_t& trk_state,
                    const propagation_state_t& propagation) const {

        /*
        /// Measurement data on surface
        const auto& meas = trk_state.measurement();

        /// Predicted bound track parameters
        const auto& predicted = trk_state.predicted();

        /// Projection matrix
        */
    }
};

}  // namespace detray