/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "detray/definitions/detail/macros.hpp"
#include "detray/geometry/detector_volume.hpp"
#include "detray/materials/interaction.hpp"
#include "detray/propagator/actors/random_scatterer.hpp"

template <typename magnetic_field_t, typename transform3_t,
          typename constraint_t, typename policy_t, typename random_device_t,
          typename inspector_t, template <typename, std::size_t> class array_t>
DETRAY_HOST_DEVICE void
detray::rk_stepper<magnetic_field_t, transform3_t, constraint_t, policy_t,
                   random_device_t, inspector_t, array_t>::state::
    advance_track(const detray::stepping::config<scalar_type>& cfg) {

    const auto& sd = this->_step_data;
    const scalar_type h{this->_step_size};
    const scalar_type h_6{h * static_cast<scalar_type>(1. / 6.)};
    auto& track = this->_track;
    auto pos = track.pos();
    auto dir = track.dir();

    // Update the track parameters according to the equations of motion
    // Reference: Eq (82) of https://doi.org/10.1016/0029-554X(81)90063-X
    pos = pos + h * (sd.t[0u] + h_6 * (sd.dtds[0] + sd.dtds[1] + sd.dtds[2]));
    track.set_pos(pos);

    // Reference: Eq (82) of https://doi.org/10.1016/0029-554X(81)90063-X
    dir =
        dir + h_6 * (sd.dtds[0] + 2.f * (sd.dtds[1] + sd.dtds[2]) + sd.dtds[3]);
    dir = vector::normalize(dir);
    track.set_dir(dir);

    auto qop = track.qop();
    const auto mat = this->_mat;
    if (!(mat == vacuum<scalar_type>())) {

        // only for simulation with valid random device
        if constexpr (!std::is_same<random_device_t,
                                    stepping::void_random_device>::value) {

            // Scatter if it is for the simulation
            if (cfg.is_simulation) {
                const scalar_type projected_scattering_angle =
                    interaction_type().compute_multiple_scattering_theta0(
                        h / mat.X0(), this->_pdg, this->_mass, qop,
                        track.charge());

                dir = random_scatterer<transform3_t>().scatter(
                    dir, projected_scattering_angle,
                    this->rand_device.generator);
            }
        }

        // Reference: Eq (82) of https://doi.org/10.1016/0029-554X(81)90063-X
        if (!cfg.is_simulation) {
            qop = qop +
                  h_6 * (sd.dqopds[0u] + 2.f * (sd.dqopds[1u] + sd.dqopds[2u]) +
                         sd.dqopds[3u]);
        } else {
            // only for simulation with valid random device
            if constexpr (!std::is_same<random_device_t,
                                        stepping::void_random_device>::value) {
                // Apply mean probable energy loss instead
                const scalar_type e_loss_mpv =
                    interaction_type().compute_energy_loss_landau(
                        h, mat, this->_pdg, this->_mass, qop, track.charge());

                const scalar_type e_loss_sigma =
                    interaction_type().compute_energy_loss_landau_sigma(
                        h, mat, this->_pdg, this->_mass, qop, track.charge());

                // Get the new momentum
                const auto new_mom = random_scatterer<transform3_t>().attenuate(
                    e_loss_mpv, e_loss_sigma, this->mass, track.charge() / qop,
                    this->rand_device.generator);

                qop = track.charge() / new_mom;
            }
        }
    }
    track.set_qop(qop);

    // Update path length
    this->_path_length += h;
    this->_path_length_per_surface += h;
}

template <typename magnetic_field_t, typename transform3_t,
          typename constraint_t, typename policy_t, typename random_device_t,
          typename inspector_t, template <typename, std::size_t> class array_t>
DETRAY_HOST_DEVICE void
detray::rk_stepper<magnetic_field_t, transform3_t, constraint_t, policy_t,
                   random_device_t, inspector_t, array_t>::state::
    advance_jacobian(const detray::stepping::config<scalar_type>& cfg) {

    // Immediately exit if it is for simulation
    if (cfg.is_simulation) {
        return;
    }

    /// The calculations are based on ATL-SOFT-PUB-2009-002. The update of the
    /// Jacobian matrix is requires only the calculation of eq. 17 and 18.
    /// Since the terms of eq. 18 are currently 0, this matrix is not needed
    /// in the calculation. The matrix A from eq. 17 consists out of 3
    /// different parts. The first one is given by the upper left 3x3 matrix
    /// that are calculated by the derivatives dF/dT (called dFdT) and dG/dT
    /// (calles dGdT). The second is given by the top 3 lines of the rightmost
    /// column. This is calculated by dFdqop and dGdqop. The remaining non-zero
    /// term is calculated directly. The naming of the variables is explained in
    /// eq. 11 and are directly related to the initial problem in eq. 7. The
    /// evaluation is based by propagating the parameters T and lambda as given
    /// in eq. 16 and evaluating the derivations for matrix A.
    /// @note The translation for u_{n+1} in eq. 7 is in this case a
    /// 3-dimensional vector without a dependency of Lambda or lambda neither in
    /// u_n nor in u_n'. The second and fourth eq. in eq. 14 have the constant
    /// offset matrices h * Id and Id respectively. This involves that the
    /// constant offset does not exist for rectangular matrix dGdu' (due to the
    /// missing Lambda part) and only exists for dFdu' in dlambda/dlambda.

    // Set transport matrix (D) and update Jacobian transport
    //( JacTransport = D * JacTransport )
    auto D = matrix_operator().template identity<e_free_size, e_free_size>();

    const auto& sd = this->_step_data;
    const scalar_type h{this->_step_size};
    // const auto& mass = this->_mass;
    auto& track = this->_track;

    // Half step length
    const scalar_type h2{h * h};
    const scalar_type half_h{h * 0.5f};
    const scalar_type h_6{h * static_cast<scalar_type>(1. / 6.)};

    // 3X3 Identity matrix
    const matrix_type<3, 3> I33 = matrix_operator().template identity<3, 3>();

    // Initialize derivatives
    std::array<matrix_type<3u, 3u>, 4u> dkndt{I33, I33, I33, I33};
    std::array<vector3, 4u> dkndqop;
    std::array<matrix_type<3u, 3u>, 4u> dkndr;
    std::array<scalar, 4u> dqopn_dqop{1.f, 1.f, 1.f, 1.f};

    /*---------------------------------------------------------------------------
     *  dk_n/dt1
     *    = qop_n * (dt_n/dt1 X B_n)
     *      + qop_n * ( t_n X dB_n/dt1 ),
     *  where dB_n/dt1 == dB_n/dr_n * dr_n/dt1.
     *
     *  The second term is non-zero only for inhomogeneous magnetic fields
     *
     *  Note that [ t_n = t1 + h * d(t_{n-1})/ds) ] as indicated by
     *  Eq (84) of https://doi.org/10.1016/0029-554X(81)90063-X

     *  [ Table for dt_n/dt1 ]
     *  dt1/dt1 = I
     *  dt2/dt1 = d( t1 + h/2 * dt1/ds ) / dt1 = I + h/2 * dk1/dt1
     *  dt3/dt1 = d( t1 + h/2 * dt2/ds ) / dt1 = I + h/2 * dk2/dt1
     *  dt4/dt1 = d( t1 + h * dt3/ds ) / dt1 = I + h * dk3/dt1
     *
     *  [ Table for dr_n/dt1 ]
     *  dr1/dt1 = 0
     *  dr2/dt1 = d(r1 + h/2 * t1 + h^2/8 dt1/ds)/dt1 = h/2 * I + h^2/8 dk1/dt1
     *  dr3/dt1 = d(r1 + h/2 * t1 + h^2/8 dt1/ds)/dt1 = h/2 * I + h^2/8 dk1/dt1
     *  dr4/dt1 = d(r1 + h * t1 + h^2/2 dt3/ds)/dt1 = h * I + h^2/2 dk3/dt1
     *
     *  Note that
     *  d/dr [ F(T) X B ]  = dF(T)/dr (X) B, where (X) means the column wise
     *  cross product
    ---------------------------------------------------------------------------*/

    /*---------------------------------------------------------------------------
     *  dk_n/dqop_1
     *    = dqop_n/dqop1 * ( t_n X B_n )
     *      + qop_n * ( dt_n/dqop1 X B_n )
     *      + qop_n * ( t_n X dB_n/dqop1 ),
     *  where dB_n/dqop1 = dB_n/dr_n * dr_n/dqop1
     *
     *  Note that [ qop_n = qop1 + h * dqop_{n-1}/ds) ] as indicated by
     *  Eq (84) of https://doi.org/10.1016/0029-554X(81)90063-X
     *
     *  [ Table for dqop_n/dqop1 ]
     *  dqop1/dqop1 = 1
     *  dqop2/dqop1 = 1 + h/2 * d(dqop1/ds)/dqop1
     *  dqop3/dqop1 = 1 + h/2 * d(dqop1/ds)/dqop1 + h^2/4 d(d^2qop1/ds^2)/dqop1
     *  dqop4/dqop1 = 1 + h * d(dqop1/ds)/dqop1 + h^2/2 d(d^2qop1/ds^2)/dqop1
     *                + h^3/4 d(d^3qop1/ds^3)/dqop1
     *
     *  [ Table for dt_n/dqop1 ]
     *  dt1/dqop1 = 0
     *  dt2/dqop1 = d(t1 + h/2 dt1/ds)/dqop1 = h/2 * dk1/dqop1
     *  dt3/dqop1 = d(t1 + h/2 dt2/ds)/dqop1 = h/2 * dk2/dqop1
     *  dt4/dqop1 = d(t1 + h dt3/ds)/dqop1 = h * dk3/dqop1
     *
     *  [ Table for dr_n/dqop1 ]
     *  dr1/dqop1 = 0
     *  dr2/dqop1 = d(r1 + h/2 * t1 + h^2/8 dt1/ds)/dqop1 = h^2/8 * dk1/dqop1
     *  dr3/dqop1 = d(r1 + h/2 * t1 + h^2/8 dt1/ds)/dqop1 = h^2/8 * dk1/dqop1
     *  dr4/dqop1 = d(r1 + h * t1 + h^2/2 dt3/ds)/dqop1 = h^2/2 dk3/dqop1
    ---------------------------------------------------------------------------*/

    /*---------------------------------------------------------------------------
     *  dk_n/dr1
     *    = qop_n * ( dt_n/dr1 X B_n )
     *      + qop_n * ( t_n X dB_n/dr1 ),
     *  where dB_n/dr1 = dB_n/dr_n * dr_n/dr1
     *
     *  [ Table for dt_n/dr1 ]
     *  dt1/dr1 = 0
     *  dt2/dr1 = d(t1 + h/2 * dt1/ds)/dr1 = h/2 * dk1/dr1
     *  dt2/dr1 = d(t1 + h/2 * dt2/ds)/dr1 = h/2 * dk2/dr1
     *  dt3/dr1 = d(t1 + h * dt3/ds)/dr1 = h * dk3/dr1
     *
     *  [ Table for dr_n/dr1 ]
     *  dr1/dr1 = I
     *  dr2/dr1 = (r1 + h/2 * t1 + h^2/8 dt1/ds ) / dr1 = I + h^2/8 dk1/dr1
     *  dr3/dr1 = (r1 + h/2 * t1 + h^2/8 dt1/ds ) / dr1 = I + h^2/8 dk1/dr1
     *  dr4/dr1 = (r1 + h * t1 + h^2/2 dt3/ds ) / dr1 = I + h^2/2 dk3/dr1
    ---------------------------------------------------------------------------*/

    /*---------------------------------------------------------------------------
     *  d(dqop_n/ds)/dqop1
     *
     *  [ Table for dqop_n/ds ]
     *  dqop1/ds = qop1^3 * E * (-dE/ds) / q^2
     *  dqop2/ds = d(qop1 + h/2 * dqop1/ds)/ds = dqop1/ds + h/2 * d^2(qop1)/ds^2
     *  dqop3/ds = d(qop1 + h/2 * dqop2/ds)/ds = dqop1/ds + h/2 * d^2(qop2)/ds^2
     *  dqop4/ds = d(qop1 + h * dqop3/ds)/ds = dqop1/ds + h * d^2(qop3)/ds^2
    ---------------------------------------------------------------------------*/

    if (!cfg.use_eloss_gradient) {
        getter::element(D, e_free_qoverp, e_free_qoverp) = 1.f;
    } else {
        // Pre-calculate dqop_n/dqop1
        // Note that terms with d^2/ds^2 or d^3/ds^3 are ignored
        const scalar_type d2qop1dsdqop1 = this->d2qopdsdqop(sd.qop[0u]);
        dqopn_dqop[0u] = 1.f;
        dqopn_dqop[1u] = 1.f + half_h * d2qop1dsdqop1;
        dqopn_dqop[2u] = dqopn_dqop[1u];
        dqopn_dqop[3u] = 1.f + h * d2qop1dsdqop1;

        /*-----------------------------------------------------------------
         * Calculate the first terms of d(dqop_n/ds)/dqop1
        -------------------------------------------------------------------*/

        // Note that terms with d^2/ds^2 are ignored
        getter::element(D, e_free_qoverp, e_free_qoverp) =
            1.f + d2qop1dsdqop1 * h;
    }

    // Calculate in the case of not considering B field gradient
    if (!cfg.use_field_gradient) {

        /*-----------------------------------------------------------------
         * Calculate the first terms of dk_n/dt1
        -------------------------------------------------------------------*/
        // dk1/dt1
        dkndt[0u] =
            sd.qop[0u] * mat_helper().column_wise_cross(dkndt[0u], sd.b_first);

        // dk2/dt1
        dkndt[1u] = dkndt[1u] + half_h * dkndt[0u];
        dkndt[1u] =
            sd.qop[1u] * mat_helper().column_wise_cross(dkndt[1u], sd.b_middle);

        // dk3/dt1
        dkndt[2u] = dkndt[2u] + half_h * dkndt[1u];
        dkndt[2u] =
            sd.qop[2u] * mat_helper().column_wise_cross(dkndt[2u], sd.b_middle);

        // dk4/dt1
        dkndt[3u] = dkndt[3u] + h * dkndt[2u];
        dkndt[3u] =
            sd.qop[3u] * mat_helper().column_wise_cross(dkndt[3u], sd.b_last);

        /*-----------------------------------------------------------------
         * Calculate the first and second terms of dk_n/dqop1
        -------------------------------------------------------------------*/
        // dk1/dqop1
        dkndqop[0u] = dqopn_dqop[0u] * vector::cross(sd.t[0u], sd.b_first);

        // dk2/dqop1
        dkndqop[1u] =
            dqopn_dqop[1u] * vector::cross(sd.t[1u], sd.b_middle) +
            sd.qop[1u] * half_h * vector::cross(dkndqop[0u], sd.b_middle);

        // dk3/dqop1
        dkndqop[2u] =
            dqopn_dqop[2u] * vector::cross(sd.t[2u], sd.b_middle) +
            sd.qop[2u] * half_h * vector::cross(dkndqop[1u], sd.b_middle);

        // dk4/dqop1
        dkndqop[3u] = dqopn_dqop[3u] * vector::cross(sd.t[3u], sd.b_last) +
                      sd.qop[3u] * h * vector::cross(dkndqop[2u], sd.b_last);
    } else {

        // Positions at four stages
        std::array<vector3, 4u> r;
        r[0u] = track.pos();
        r[1u] = r[0u] + half_h * sd.t[0u] + h2 * 0.125f * sd.dtds[0u];
        r[2u] = r[1u];
        r[3u] = r[0u] + h * sd.t[0u] + h2 * 0.5f * sd.dtds[2u];

        // Field gradients at four stages
        std::array<matrix_type<3, 3>, 4u> dBdr;
        dBdr[0u] = evaluate_field_gradient(r[0u]);
        dBdr[1u] = evaluate_field_gradient(r[1u]);
        dBdr[2u] = dBdr[1u];
        dBdr[3u] = evaluate_field_gradient(r[3u]);

        // Temporary variable for dBdt and dBdr
        matrix_type<3u, 3u> dBdt_tmp;
        matrix_type<3u, 3u> dBdr_tmp;

        /*-----------------------------------------------------------------
         * Calculate all terms of dk_n/dt1
        -------------------------------------------------------------------*/
        // dk1/dt1
        dkndt[0u] =
            sd.qop[0u] * mat_helper().column_wise_cross(dkndt[0u], sd.b_first);

        // dk2/dt1
        dkndt[1u] = dkndt[1u] + half_h * dkndt[0u];
        dkndt[1u] =
            sd.qop[1u] * mat_helper().column_wise_cross(dkndt[1u], sd.b_middle);
        dBdt_tmp = dBdr[1u] * (half_h * I33 + h2 * 0.125f * dkndt[0u]);
        dkndt[1u] = dkndt[1u] - sd.qop[1u] * mat_helper().column_wise_cross(
                                                 dBdt_tmp, sd.t[1u]);

        // dk3/dt1
        dkndt[2u] = dkndt[2u] + half_h * dkndt[1u];
        dkndt[2u] =
            sd.qop[2u] * mat_helper().column_wise_cross(dkndt[2u], sd.b_middle);
        dBdt_tmp = dBdr[2u] * (half_h * I33 + h2 * 0.125f * dkndt[0u]);
        dkndt[2u] = dkndt[2u] - sd.qop[2u] * mat_helper().column_wise_cross(
                                                 dBdt_tmp, sd.t[2u]);

        // dk4/dt1
        dkndt[3u] = dkndt[3u] + h * dkndt[2u];
        dkndt[3u] =
            sd.qop[3u] * mat_helper().column_wise_cross(dkndt[3u], sd.b_last);
        dBdt_tmp = dBdr[3u] * (h * I33 + h2 * 0.5f * dkndt[2u]);
        dkndt[3u] = dkndt[3u] - sd.qop[3u] * mat_helper().column_wise_cross(
                                                 dBdt_tmp, sd.t[3u]);

        /*-----------------------------------------------------------------
         * Calculate all terms of dk_n/dqop1
        -------------------------------------------------------------------*/
        // dk1/dqop1
        dkndqop[0u] = dqopn_dqop[0u] * vector::cross(sd.t[0u], sd.b_first);

        // dk2/dqop1
        dkndqop[1u] =
            dqopn_dqop[1u] * vector::cross(sd.t[1u], sd.b_middle) +
            sd.qop[1u] * half_h * vector::cross(dkndqop[0u], sd.b_middle);
        dkndqop[1u] =
            dkndqop[1u] +
            sd.qop[1u] *
                vector::cross(sd.t[1u], h2 * 0.125f * dBdr[1u] * dkndqop[0u]);

        // dk3/dqop1
        dkndqop[2u] =
            dqopn_dqop[2u] * vector::cross(sd.t[2u], sd.b_middle) +
            sd.qop[2u] * half_h * vector::cross(dkndqop[1u], sd.b_middle);
        dkndqop[2u] =
            dkndqop[2u] +
            sd.qop[2u] *
                vector::cross(sd.t[2u], h2 * 0.125f * dBdr[2u] * dkndqop[0u]);

        // dk4/dqop1
        dkndqop[3u] = dqopn_dqop[3u] * vector::cross(sd.t[3u], sd.b_last) +
                      sd.qop[3u] * h * vector::cross(dkndqop[2u], sd.b_last);
        dkndqop[3u] =
            dkndqop[3u] +
            sd.qop[3u] *
                vector::cross(sd.t[3u], h2 * 0.5f * dBdr[3u] * dkndqop[2u]);

        /*-----------------------------------------------------------------
         * Calculate all terms of dk_n/dr1
        -------------------------------------------------------------------*/
        // dk1/dr1
        dkndr[0u] =
            -sd.qop[0u] * mat_helper().column_wise_cross(dBdr[0u], sd.t[0u]);

        // dk2/dr1
        dkndr[1u] = sd.qop[1u] * mat_helper().column_wise_cross(
                                     half_h * dkndr[0u], sd.b_middle);
        dBdr_tmp = dBdr[1u] * (I33 + h2 * 0.125 * dkndr[0u]);
        dkndr[1u] = dkndr[1u] - sd.qop[1u] * mat_helper().column_wise_cross(
                                                 dBdr_tmp, sd.t[1u]);

        // dk3/dr1
        dkndr[2u] = sd.qop[2u] * mat_helper().column_wise_cross(
                                     half_h * dkndr[1u], sd.b_middle);
        dBdr_tmp = dBdr[2u] * (I33 + h2 * 0.125 * dkndr[0u]);
        dkndr[2u] = dkndr[2u] - sd.qop[2u] * mat_helper().column_wise_cross(
                                                 dBdr_tmp, sd.t[2u]);

        // dk4/dr1
        dkndr[3u] = sd.qop[3u] *
                    mat_helper().column_wise_cross(h * dkndr[2u], sd.b_last);
        dBdr_tmp = dBdr[3u] * (I33 + h2 * 0.5 * dkndr[2u]);
        dkndr[3u] = dkndr[3u] - sd.qop[3u] * mat_helper().column_wise_cross(
                                                 dBdr_tmp, sd.t[3u]);

        // Set dF/dr1 and dG/dr1
        auto dFdr = matrix_operator().template identity<3, 3>();
        auto dGdr = matrix_operator().template identity<3, 3>();
        dFdr = dFdr + h * h_6 * (dkndr[0u] + dkndr[1u] + dkndr[2u]);
        dGdr = h_6 * (dkndr[0u] + 2.f * (dkndr[1u] + dkndr[2u]) + dkndr[3u]);

        matrix_operator().set_block(D, dFdr, 0u, 0u);
        matrix_operator().set_block(D, dGdr, 4u, 0u);
    }

    // Set dF/dt1 and dG/dt1
    auto dFdt = matrix_operator().template identity<3, 3>();
    auto dGdt = matrix_operator().template identity<3, 3>();
    dFdt = dFdt + h_6 * (dkndt[0u] + dkndt[1u] + dkndt[2u]);
    dFdt = h * dFdt;
    dGdt = dGdt + h_6 * (dkndt[0u] + 2.f * (dkndt[1u] + dkndt[2u]) + dkndt[3u]);

    matrix_operator().set_block(D, dFdt, 0u, 4u);
    matrix_operator().set_block(D, dGdt, 4u, 4u);

    // Set dF/dqop1 and dG/dqop1
    vector3 dFdqop = h * h_6 * (dkndqop[0u] + dkndqop[1u] + dkndqop[2u]);
    vector3 dGdqop =
        h_6 * (dkndqop[0u] + 2.f * (dkndqop[1u] + dkndqop[2u]) + dkndqop[3u]);
    matrix_operator().set_block(D, dFdqop, 0u, 7u);
    matrix_operator().set_block(D, dGdqop, 4u, 7u);

    this->_jac_transport = D * this->_jac_transport;
}

template <typename magnetic_field_t, typename transform3_t,
          typename constraint_t, typename policy_t, typename random_device_t,
          typename inspector_t, template <typename, std::size_t> class array_t>
DETRAY_HOST_DEVICE void
detray::rk_stepper<magnetic_field_t, transform3_t, constraint_t, policy_t,
                   random_device_t, inspector_t, array_t>::state::
    add_multiple_scattering_covariance(
        const detray::stepping::config<scalar_type>& cfg) {

    // Return if there is no material
    if (this->mat == vacuum<scalar_type>()) {
        return;
    }

    // Return if it is for simulation
    if (cfg.is_simulation) {
        return;
    }

    // Implement the thick scatterer method of Eq (4.103 - 111) of "Pattern
    // Recognition, Tracking and Vertex Reconstruction in Particle Detectors"

    // Variance per unit length of the projected scattering angle
    const scalar s1 = interaction_type().compute_multiple_scattering_theta0(
        1.f / this->mat.X0(), this->_pdg, this->_mass, this->_qop_i,
        this->_track.charge());
    const scalar s2 = interaction_type().compute_multiple_scattering_theta0(
        1.f / this->mat.X0(), this->_pdg, this->_mass, this->_track.qop(),
        this->_track.charge());
    const scalar variance = s1 * s2;

    // 4X4 Identity matrix
    matrix_type<4, 4> E = matrix_operator().template zero<4, 4>();

    const scalar L = this->_path_length_per_surface;
    const scalar half_L2 = L * L / 2.f;
    const scalar third_L3 = L * L * L / 3.f;

    // Set E (Eq. 4.107). Note that we flip the index of (0,1) <-> (2,3)
    getter::element(E, 0, 0) = third_L3;
    getter::element(E, 1, 1) = third_L3;
    getter::element(E, 2, 2) = L;
    getter::element(E, 3, 3) = L;
    getter::element(E, 0, 2) = half_L2;
    getter::element(E, 1, 3) = half_L2;
    getter::element(E, 2, 0) = half_L2;
    getter::element(E, 3, 1) = half_L2;

    const matrix_type<3, 3> C2G =
        mat_helper().curvilinear_to_global(this->_track.dir());
    // Take the bottom-right 2x2 block
    const matrix_type<2, 2> block =
        matrix_operator().template block<2, 2>(C2G, 0, 0);

    // Jacobian for (Eq. 4.108).  Note that we flip the index of (0,1) <-> (2,3)
    matrix_type<4, 4> T = matrix_operator().template identity<4, 4>();
    matrix_operator().template set_block<2, 2>(T, block, 2, 2);

    E = T * E * matrix_operator().transpose(T);

    // Jacobian d(phi,theta)/d(tx,ty). Note that it is different from
    // (Eq. 4.109) due to the different coordinate system
    matrix_type<4, 4> J = matrix_operator().template identity<4, 4>();

    const auto t = this->_track.dir();
    const scalar_type phi = getter::phi(t);
    const scalar_type theta = getter::theta(t);
    const scalar_type sin_phi = math::sin(phi);
    const scalar_type cos_phi = math::cos(phi);
    const scalar_type sin_theta = math::sin(theta);
    const scalar_type cos_theta = math::cos(theta);

    // -sin(phi)/sin(theta)
    getter::element(J, 2, 2) = -sin_phi / sin_theta;
    getter::element(J, 2, 3) = -cos_phi / sin_theta;
    getter::element(J, 3, 2) = -cos_phi * cos_theta;
    getter::element(J, 3, 3) = -sin_phi * cos_theta;

    E = J * E * matrix_operator().transpose(J);

    // Set the joint covariance
    matrix_operator().set_block<4, 4>(this->_joint_cov, E, 0, 0);
}

template <typename magnetic_field_t, typename transform3_t,
          typename constraint_t, typename policy_t, typename random_device_t,
          typename inspector_t, template <typename, std::size_t> class array_t>
DETRAY_HOST_DEVICE auto detray::rk_stepper<
    magnetic_field_t, transform3_t, constraint_t, policy_t, random_device_t,
    inspector_t,
    array_t>::state::evaluate_dqopds(const std::size_t i,
                                     const typename transform3_t::scalar_type h,
                                     const scalar dqopds_prev) ->
    typename transform3_t::scalar_type {

    const auto& track = this->_track;
    const scalar_type qop = track.qop();
    auto& sd = this->_step_data;

    if (this->_mat == detray::vacuum<scalar_type>()) {
        sd.qop[i] = qop;
        return 0.f;
    } else {

        if (i == 0u) {
            sd.qop[i] = qop;
        }
        // For reconstruction with mean energy loss
        else {

            // qop_n is calculated recursively like the direction of
            // evaluate_dtds.
            //
            // https://doi.org/10.1016/0029-554X(81)90063-X says:
            // "For y  we  have  similar  formulae  as  for x, for y' and
            // \lambda similar  formulae as for  x'"
            sd.qop[i] = qop + h * dqopds_prev;
        }

        return this->dqopds(sd.qop[i]);
    }
}

template <typename magnetic_field_t, typename transform3_t,
          typename constraint_t, typename policy_t, typename random_device_t,
          typename inspector_t, template <typename, std::size_t> class array_t>
DETRAY_HOST_DEVICE auto detray::rk_stepper<
    magnetic_field_t, transform3_t, constraint_t, policy_t, random_device_t,
    inspector_t,
    array_t>::state::evaluate_dtds(const vector3& b_field, const std::size_t i,
                                   const typename transform3_t::scalar_type h,
                                   const vector3& dtds_prev,
                                   const typename transform3_t::scalar_type qop)
    -> vector3 {
    auto& track = this->_track;
    const auto dir = track.dir();
    auto& sd = this->_step_data;

    if (i == 0u) {
        sd.t[i] = dir;
    } else {
        // Eq (84) of https://doi.org/10.1016/0029-554X(81)90063-X
        sd.t[i] = dir + h * dtds_prev;
    }

    // dtds = qop * (t X B) from Lorentz force
    return qop * vector::cross(sd.t[i], b_field);
}

template <typename magnetic_field_t, typename transform3_t,
          typename constraint_t, typename policy_t, typename random_device_t,
          typename inspector_t, template <typename, std::size_t> class array_t>
DETRAY_HOST_DEVICE auto detray::rk_stepper<
    magnetic_field_t, transform3_t, constraint_t, policy_t, random_device_t,
    inspector_t, array_t>::state::evaluate_field_gradient(const vector3& pos)
    -> matrix_type<3, 3> {

    matrix_type<3, 3> dBdr = matrix_operator().template zero<3, 3>();

    constexpr auto delta{1e-1f * unit<scalar_type>::mm};

    for (unsigned int i = 0; i < 3; i++) {

        vector3 dpos1 = pos;
        dpos1[i] += delta;
        const auto bvec1_tmp =
            this->_magnetic_field.at(dpos1[0], dpos1[1], dpos1[2]);
        vector3 bvec1;
        bvec1[0u] = bvec1_tmp[0u];
        bvec1[1u] = bvec1_tmp[1u];
        bvec1[2u] = bvec1_tmp[2u];

        vector3 dpos2 = pos;
        dpos2[i] -= delta;
        const auto bvec2_tmp =
            this->_magnetic_field.at(dpos2[0], dpos2[1], dpos2[2]);
        vector3 bvec2;
        bvec2[0u] = bvec2_tmp[0u];
        bvec2[1u] = bvec2_tmp[1u];
        bvec2[2u] = bvec2_tmp[2u];

        const vector3 gradient = (bvec1 - bvec2) * (1.f / (2.f * delta));

        getter::element(dBdr, 0u, i) = gradient[0u];
        getter::element(dBdr, 1u, i) = gradient[1u];
        getter::element(dBdr, 2u, i) = gradient[2u];
    }

    return dBdr;
}

template <typename magnetic_field_t, typename transform3_t,
          typename constraint_t, typename policy_t, typename random_device_t,
          typename inspector_t, template <typename, std::size_t> class array_t>
DETRAY_HOST_DEVICE auto
detray::rk_stepper<magnetic_field_t, transform3_t, constraint_t, policy_t,
                   random_device_t, inspector_t, array_t>::state::dqopds() const
    -> typename transform3_t::scalar_type {
    return this->_step_data.dqopds[3u];
}

template <typename magnetic_field_t, typename transform3_t,
          typename constraint_t, typename policy_t, typename random_device_t,
          typename inspector_t, template <typename, std::size_t> class array_t>
DETRAY_HOST_DEVICE auto detray::rk_stepper<
    magnetic_field_t, transform3_t, constraint_t, policy_t, random_device_t,
    inspector_t, array_t>::state::dqopds(const scalar_type qop) const ->
    typename transform3_t::scalar_type {

    const auto& mat = this->_mat;

    // d(qop)ds is zero for empty space
    if (mat == detray::vacuum<scalar_type>()) {
        return 0.f;
    }

    const auto pdg = this->_pdg;
    const scalar_type q = this->_track.charge();
    const scalar_type p = q / qop;
    const scalar_type mass = this->_mass;
    const scalar_type E = math::sqrt(p * p + mass * mass);

    // Compute stopping power
    const scalar_type stopping_power =
        interaction<scalar_type>().compute_stopping_power(mat, pdg,
                                                          {mass, qop, q});

    // Assert that a momentum is a positive value
    assert(p >= 0.f);

    // d(qop)ds, which is equal to (qop) * E * (-dE/ds) / p^2
    // or equal to (qop)^3 * E * (-dE/ds) / q^2
    return qop * qop * qop * E * stopping_power / (q * q);
}

template <typename magnetic_field_t, typename transform3_t,
          typename constraint_t, typename policy_t, typename random_device_t,
          typename inspector_t, template <typename, std::size_t> class array_t>
DETRAY_HOST_DEVICE auto detray::rk_stepper<
    magnetic_field_t, transform3_t, constraint_t, policy_t, random_device_t,
    inspector_t, array_t>::state::d2qopdsdqop(const scalar_type qop) const ->
    typename transform3_t::scalar_type {

    using scalar_t = typename transform3_t::scalar_type;

    if (this->_mat == vacuum<scalar_t>()) {
        return 0.f;
    }

    auto& track = this->_track;
    const scalar_t q = track.charge();
    const scalar_t p = q / qop;
    const scalar_t p2 = p * p;

    const auto& mass = this->_mass;
    const scalar_t E2 = p2 + mass * mass;

    // Interaction object
    interaction<scalar_t> I;

    // g = dE/ds = -1 * (-dE/ds) = -1 * stopping power
    const detail::relativistic_quantities<scalar_t> rq(mass, qop, q);
    // We assume that stopping power ~ mean ionization eloss per pathlength
    const scalar_type bethe = I.compute_bethe(this->_mat, rq);
    const scalar_type g = -1.f * bethe;

    // dg/d(qop) = -1 * derivation of stopping power
    const scalar_t dgdqop =
        -1.f * interaction<scalar_t>().derive_bethe(this->_mat, rq, bethe);

    // d(qop)/ds = - qop^3 * E * g / q^2
    const scalar_t dqopds = this->dqopds(qop);

    // Check Eq 3.12 of
    // (https://iopscience.iop.org/article/10.1088/1748-0221/4/04/P04016/meta)
    return dqopds * (1.f / qop * (3.f - p2 / E2) + 1.f / g * dgdqop);
}

template <typename magnetic_field_t, typename transform3_t,
          typename constraint_t, typename policy_t, typename random_device_t,
          typename inspector_t, template <typename, std::size_t> class array_t>
template <typename propagation_state_t>
DETRAY_HOST_DEVICE bool detray::rk_stepper<
    magnetic_field_t, transform3_t, constraint_t, policy_t, random_device_t,
    inspector_t, array_t>::step(propagation_state_t& propagation,
                                const detray::stepping::config<scalar_type>&
                                    cfg) {

    // Get stepper and navigator states
    state& stepping = propagation._stepping;
    auto& magnetic_field = stepping._magnetic_field;
    auto& navigation = propagation._navigation;

    auto vol = detector_volume{*navigation.detector(), navigation.volume()};
    stepping._mat = vol.material();

    auto& sd = stepping._step_data;

    scalar_type error_estimate{0.f};

    // First Runge-Kutta point
    const vector3 pos = stepping().pos();
    const auto bvec = magnetic_field.at(pos[0], pos[1], pos[2]);
    sd.b_first[0] = bvec[0];
    sd.b_first[1] = bvec[1];
    sd.b_first[2] = bvec[2];

    // qop should be recalcuated at every point
    // Reference: Eq (84) of https://doi.org/10.1016/0029-554X(81)90063-X
    sd.dqopds[0u] = stepping.evaluate_dqopds(0u, 0.f, 0.f);
    sd.dtds[0u] = stepping.evaluate_dtds(sd.b_first, 0u, 0.f,
                                         vector3{0.f, 0.f, 0.f}, sd.qop[0u]);

    const auto try_rk4 = [&](const scalar_type& h) -> bool {
        // State the square and half of the step size
        const scalar_type h2{h * h};
        const scalar_type half_h{h * 0.5f};

        // Second Runge-Kutta point
        // qop should be recalcuated at every point
        // Eq (84) of https://doi.org/10.1016/0029-554X(81)90063-X
        const vector3 pos1 =
            pos + half_h * sd.t[0u] + h2 * 0.125f * sd.dtds[0u];
        const auto bvec1 = magnetic_field.at(pos1[0], pos1[1], pos1[2]);
        sd.b_middle[0] = bvec1[0];
        sd.b_middle[1] = bvec1[1];
        sd.b_middle[2] = bvec1[2];

        sd.dqopds[1u] = stepping.evaluate_dqopds(1u, half_h, sd.dqopds[0u]);
        sd.dtds[1u] = stepping.evaluate_dtds(sd.b_middle, 1u, half_h,
                                             sd.dtds[0u], sd.qop[1u]);

        // Third Runge-Kutta point
        // qop should be recalcuated at every point
        // Reference: Eq (84) of https://doi.org/10.1016/0029-554X(81)90063-X
        sd.dqopds[2u] = stepping.evaluate_dqopds(2u, half_h, sd.dqopds[1u]);
        sd.dtds[2u] = stepping.evaluate_dtds(sd.b_middle, 2u, half_h,
                                             sd.dtds[1u], sd.qop[2u]);

        // Last Runge-Kutta point
        // qop should be recalcuated at every point
        // Eq (84) of https://doi.org/10.1016/0029-554X(81)90063-X
        const vector3 pos2 = pos + h * sd.t[0u] + h2 * 0.5f * sd.dtds[2u];
        const auto bvec2 = magnetic_field.at(pos2[0], pos2[1], pos2[2]);
        sd.b_last[0] = bvec2[0];
        sd.b_last[1] = bvec2[1];
        sd.b_last[2] = bvec2[2];

        sd.dqopds[3u] = stepping.evaluate_dqopds(3u, h, sd.dqopds[2u]);
        sd.dtds[3u] =
            stepping.evaluate_dtds(sd.b_last, 3u, h, sd.dtds[2u], sd.qop[3u]);

        // Compute and check the local integration error estimate
        // @Todo
        constexpr const scalar_type one_sixth{
            static_cast<scalar_type>(1. / 6.)};
        const vector3 err_vec =
            one_sixth * h2 *
            (sd.dtds[0u] - sd.dtds[1u] - sd.dtds[2u] + sd.dtds[3u]);
        error_estimate =
            math::max(getter::norm(err_vec), static_cast<scalar_type>(1e-20));

        return (error_estimate <= cfg.rk_error_tol);
    };

    // Initial step size estimate
    stepping.set_step_size(navigation());

    scalar_type step_size_scaling{1.f};
    std::size_t n_step_trials{0u};

    // Adjust initial step size to integration error
    while (!try_rk4(stepping._step_size)) {

        step_size_scaling = math::min(
            math::max(static_cast<scalar>(0.25),
                      math::sqrt(math::sqrt(
                          (cfg.rk_error_tol / math::abs(error_estimate))))),
            static_cast<scalar_type>(4));

        // Only step size reduction is allowed so that we don't overstep
        assert(step_size_scaling <= 1.f);

        stepping._step_size *= step_size_scaling;

        // If step size becomes too small the particle remains at the
        // initial place
        if (math::abs(stepping._step_size) < math::abs(cfg.min_stepsize)) {
            // Not moving due to too low momentum needs an aborter
            return navigation.abort();
        }

        // If the parameter is off track too much or given step_size is not
        // appropriate
        if (n_step_trials > cfg.max_rk_updates) {
            // Too many trials, have to abort
            return navigation.abort();
        }
        n_step_trials++;

        // Run inspection while the stepsize is getting adjusted
        stepping.run_inspector(cfg, "Adjust stepsize: ", n_step_trials,
                               step_size_scaling);
    }

    // Update navigation direction
    const step::direction step_dir = stepping._step_size >= 0.f
                                         ? step::direction::e_forward
                                         : step::direction::e_backward;
    stepping.set_direction(step_dir);

    // Check constraints
    if (math::abs(stepping.step_size()) >
        math::abs(
            stepping.constraints().template size<>(stepping.direction()))) {

        // Run inspection before step size is cut
        stepping.run_inspector(cfg, "Before constraint: ");

        stepping.set_step_size(
            stepping.constraints().template size<>(stepping.direction()));
    }

    // Advance track state
    stepping.advance_track(cfg);

    // Advance jacobian transport
    stepping.advance_jacobian(cfg);

    // Call navigation update policy
    typename rk_stepper::policy_type{}(stepping.policy_state(), propagation);

    // Run final inspection
    stepping.run_inspector(cfg, "Step complete: ");

    return true;
}
