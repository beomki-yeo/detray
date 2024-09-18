/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/units.hpp"
#include "detray/materials/interaction.hpp"
#include "detray/materials/material.hpp"
#include "detray/materials/predefined_materials.hpp"
#include "detray/navigation/policies.hpp"
#include "detray/propagator/base_stepper.hpp"
#include "detray/tracks/tracks.hpp"
#include "detray/utils/matrix_helper.hpp"

namespace detray {

/// Runge-Kutta-Nystrom 4th order stepper implementation
///
/// @tparam magnetic_field_t the type of magnetic field
/// @tparam track_t the type of track that is being advanced by the stepper
/// @tparam constraint_ the type of constraints on the stepper
template <typename magnetic_field_t, typename algebra_t,
          typename constraint_t = unconstrained_step,
          typename policy_t = stepper_rk_policy,
          typename inspector_t = stepping::void_inspector,
          template <typename, std::size_t> class array_t = darray>
class rk_stepper final
    : public base_stepper<algebra_t, constraint_t, policy_t, inspector_t> {

    public:
    using base_type =
        base_stepper<algebra_t, constraint_t, policy_t, inspector_t>;

    using algebra_type = algebra_t;
    using scalar_type = dscalar<algebra_t>;
    using point3_type = dpoint3D<algebra_t>;
    using vector3_type = dvector3D<algebra_t>;
    using transform3_type = dtransform3D<algebra_t>;
    using matrix_operator = dmatrix_operator<algebra_t>;
    using mat_helper = matrix_helper<matrix_operator>;
    using free_track_parameters_type =
        typename base_type::free_track_parameters_type;
    using bound_track_parameters_type =
        typename base_type::bound_track_parameters_type;
    using magnetic_field_type = magnetic_field_t;
    template <std::size_t ROWS, std::size_t COLS>
    using matrix_type = dmatrix<algebra_t, ROWS, COLS>;

    rk_stepper() = default;

    struct state : public base_type::state {

        static constexpr const stepping::id id = stepping::id::e_rk;

        DETRAY_HOST_DEVICE
        state(const free_track_parameters_type& t,
              const magnetic_field_t& mag_field)
            : base_type::state(t), _magnetic_field(mag_field) {}

        template <typename detector_t>
        DETRAY_HOST_DEVICE state(
            const bound_track_parameters_type& bound_params,
            const magnetic_field_t& mag_field, const detector_t& det)
            : base_type::state(bound_params, det), _magnetic_field(mag_field) {}

        /// stepping data required for RKN4
        struct {
            vector3_type b_first{0.f, 0.f, 0.f};
            vector3_type b_middle{0.f, 0.f, 0.f};
            vector3_type b_last{0.f, 0.f, 0.f};
            // t = tangential direction = dr/ds
            std::array<vector3_type, 4u> t;
            // dt/ds = d^2r/ds^2 = q/p ( t X B )
            std::array<vector3_type, 4u> dtds;
            // dqopds at the last RK point
            scalar_type dqopds_last = 0.f;
        } _step_data;

        /// Additional data with volume material
        struct qop_data {
            // q/p
            std::array<scalar_type, 4u> qop{0.f,0.f,0.f,0.f};
            // d(q/p)/ds
            std::array<scalar_type, 4u> dqopds{0.f,0.f,0.f,0.f};
        };

        /// Magnetic field view
        const magnetic_field_t _magnetic_field;

        /// Material that track is passing through. Usually a volume material
        const detray::material<scalar_type>* _mat{nullptr};

        /// Access the current volume material
        DETRAY_HOST_DEVICE
        const auto& volume_material() const {
            assert(_mat != nullptr);
            return *_mat;
        }

        /// Update the track state w/o volume material
        DETRAY_HOST_DEVICE
        inline void advance_track();

        /// Update the track state w/ volume material
        DETRAY_HOST_DEVICE
        inline void advance_track(const qop_data& qd);

        /// Update the jacobian transport from free propagation
        DETRAY_HOST_DEVICE
        inline void advance_jacobian(const stepping::config& cfg = {});

        DETRAY_HOST_DEVICE
        inline void advance_jacobian(const qop_data& qd,
                                     const stepping::config& cfg = {});

        /// evaulate dqopds for a given step size and material
        DETRAY_HOST_DEVICE
        inline scalar_type evaluate_dqopds(const std::size_t i,
                                           const scalar_type h,
                                           const scalar_type dqopds_prev,
                                           const detray::stepping::config& cfg,
                                           qop_data& qd);

        /// evaulate dtds for runge kutta stepping
        DETRAY_HOST_DEVICE
        inline vector3_type evaluate_dtds(const vector3_type& b_field,
                                          const std::size_t i,
                                          const scalar_type h,
                                          const vector3_type& dtds_prev,
                                          const scalar_type qop);

        DETRAY_HOST_DEVICE
        inline matrix_type<3, 3> evaluate_field_gradient(
            const point3_type& pos);

        scalar_type scale_step_size(auto&& error_estimator,
                                    const stepping::config& cfg) {

            scalar_type error{1e20f};

            // Whenever navigator::init() is called the step size is set to
            // navigation path length (navigation()). We need to reduce it down
            // to make error small enough
            for (unsigned int i_t = 0u; i_t < cfg.max_rk_updates; i_t++) {
                this->count_trials();

                error = math::max(error_estimator(this->_step_size),
                                  static_cast<scalar_type>(1e-20));

                // Error is small enough
                // ---> break and advance track
                if (error <= 4.f * cfg.rk_error_tol) {
                    break;
                }
                // Error estimate is too big
                // ---> Make step size smaller and esimate error again
                else {

                    scalar_type step_size_scaling =
                        math::sqrt(math::sqrt(cfg.rk_error_tol / error));

                    this->_step_size *= step_size_scaling;

                    // Run inspection while the stepsize is getting adjusted
                    this->run_inspector(cfg, "Adjust stepsize: ", i_t + 1,
                                        step_size_scaling);
                }
            }

            return error;
        }

        void pre_step_update(const stepping::config& cfg) {

            // Update navigation direction
            const step::direction step_dir = this->_step_size >= 0.f
                                                 ? step::direction::e_forward
                                                 : step::direction::e_backward;
            this->set_direction(step_dir);

            // Check constraints
            if (math::fabs(this->_step_size) >
                math::fabs(
                    this->constraints().template size<>(this->direction()))) {

                // Run inspection before step size is cut
                this->run_inspector(cfg, "Before constraint: ");

                this->set_step_size(
                    this->constraints().template size<>(this->direction()));
            }
        }

        void post_step_update(const scalar_type error,
                              const stepping::config& cfg) {

            const scalar_type step_size_scaling =
                static_cast<scalar_type>(math::min(
                    math::max(math::sqrt(math::sqrt(cfg.rk_error_tol / error)),
                              static_cast<scalar_type>(0.25)),
                    static_cast<scalar_type>(4.)));

            // Save the current step size
            this->_prev_step_size = this->_step_size;

            // Update the step size
            this->_step_size *= step_size_scaling;
        }

        /// Evaluate dtds, where t is the unit tangential direction
        DETRAY_HOST_DEVICE
        inline vector3_type dtds() const;

        /// Evaulate d(qop)/ds
        DETRAY_HOST_DEVICE
        inline scalar_type dqopds() const;

        DETRAY_HOST_DEVICE
        inline scalar_type dqopds(const scalar_type qop) const;

        /// Evaulate d(d(qop)/ds)dqop
        DETRAY_HOST_DEVICE
        inline scalar_type d2qopdsdqop(const scalar_type qop) const;

        /// Call the stepping inspector
        template <typename... Args>
        DETRAY_HOST_DEVICE inline void run_inspector(
            [[maybe_unused]] const stepping::config& cfg,
            [[maybe_unused]] const char* message,
            [[maybe_unused]] Args&&... args) {
            if constexpr (!std::is_same_v<inspector_t,
                                          stepping::void_inspector>) {
                this->_inspector(*this, cfg, message,
                                 std::forward<Args>(args)...);
            }
        }
    };

    /// Take a step, using an adaptive Runge-Kutta algorithm.
    ///
    /// @return returning the heartbeat, indicating if the stepping is alive
    template <typename propagation_state_t>
    DETRAY_HOST_DEVICE bool step(propagation_state_t& propagation,
                                 const stepping::config& cfg) const;
};

}  // namespace detray

#include "detray/propagator/rk_stepper.ipp"
