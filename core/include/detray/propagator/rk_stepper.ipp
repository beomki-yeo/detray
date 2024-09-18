/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "detray/geometry/tracking_volume.hpp"

// System include(s).
#include <type_traits>

template <typename magnetic_field_t, typename algebra_t, typename constraint_t,
          typename policy_t, typename inspector_t,
          template <typename, std::size_t> class array_t>
DETRAY_HOST_DEVICE void
detray::rk_stepper<magnetic_field_t, algebra_t, constraint_t, policy_t,
                   inspector_t, array_t>::state::advance_track() {

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

    // Update path length
    this->_path_length += h;
    this->_abs_path_length += math::fabs(h);
    this->_s += h;
}

template <typename magnetic_field_t, typename algebra_t, typename constraint_t,
          typename policy_t, typename inspector_t,
          template <typename, std::size_t> class array_t>
DETRAY_HOST_DEVICE void detray::rk_stepper<
    magnetic_field_t, algebra_t, constraint_t, policy_t, inspector_t,
    array_t>::state::advance_track(const qop_data& qd) {

    this->advance_track();

    const scalar_type h_6{this->_step_size * static_cast<scalar_type>(1. / 6.)};

    auto qop = this->_track.qop();
    // Reference: Eq (82) of https://doi.org/10.1016/0029-554X(81)90063-X
    qop = qop + h_6 * (qd.dqopds[0u] + 2.f * (qd.dqopds[1u] + qd.dqopds[2u]) +
                       qd.dqopds[3u]);

    this->_track.set_qop(qop);
}

template <typename magnetic_field_t, typename algebra_t, typename constraint_t,
          typename policy_t, typename inspector_t,
          template <typename, std::size_t> class array_t>
DETRAY_HOST_DEVICE void detray::rk_stepper<
    magnetic_field_t, algebra_t, constraint_t, policy_t, inspector_t,
    array_t>::state::advance_jacobian(const detray::stepping::config& cfg) {

    // Reference DOI: NIM A 1068 (2024) 169734

    auto D = matrix_operator().template identity<e_free_size, e_free_size>();

    const auto& sd = this->_step_data;
    const scalar_type h{this->_step_size};
    auto& track = this->_track;

    // Half step length
    const scalar_type h2{h * h};
    const scalar_type half_h{h * 0.5f};
    const scalar_type h_6{h * static_cast<scalar_type>(1. / 6.)};

    // 3X3 Identity matrix
    const matrix_type<3, 3> I33 = matrix_operator().template identity<3, 3>();

    // Initialize derivatives
    std::array<matrix_type<3u, 3u>, 4u> dkndt{I33, I33, I33, I33};
    std::array<matrix_type<3u, 3u>, 4u> dkndr;

    const auto qop = track.qop();

    // Calculate in the case of not considering B field gradient
    if (!cfg.use_field_gradient) {

        /*-----------------------------------------------------------------
         * Calculate the first terms of dk_n/dt1
        -------------------------------------------------------------------*/
        // dk1/dt1
        dkndt[0u] = qop * mat_helper().column_wise_cross(dkndt[0u], sd.b_first);

        // dk2/dt1
        dkndt[1u] = dkndt[1u] + half_h * dkndt[0u];
        dkndt[1u] =
            qop * mat_helper().column_wise_cross(dkndt[1u], sd.b_middle);

        // dk3/dt1
        dkndt[2u] = dkndt[2u] + half_h * dkndt[1u];
        dkndt[2u] =
            qop * mat_helper().column_wise_cross(dkndt[2u], sd.b_middle);

        // dk4/dt1
        dkndt[3u] = dkndt[3u] + h * dkndt[2u];
        dkndt[3u] = qop * mat_helper().column_wise_cross(dkndt[3u], sd.b_last);

    } else {

        // Positions at four stages
        std::array<vector3_type, 4u> r;
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
        dkndt[0u] = qop * mat_helper().column_wise_cross(dkndt[0u], sd.b_first);

        // dk2/dt1
        dkndt[1u] = dkndt[1u] + half_h * dkndt[0u];
        dkndt[1u] =
            qop * mat_helper().column_wise_cross(dkndt[1u], sd.b_middle);
        dBdt_tmp = dBdr[1u] * (half_h * I33 + h2 * 0.125f * dkndt[0u]);
        dkndt[1u] = dkndt[1u] -
                    qop * mat_helper().column_wise_cross(dBdt_tmp, sd.t[1u]);

        // dk3/dt1
        dkndt[2u] = dkndt[2u] + half_h * dkndt[1u];
        dkndt[2u] =
            qop * mat_helper().column_wise_cross(dkndt[2u], sd.b_middle);
        dBdt_tmp = dBdr[2u] * (half_h * I33 + h2 * 0.125f * dkndt[0u]);
        dkndt[2u] = dkndt[2u] -
                    qop * mat_helper().column_wise_cross(dBdt_tmp, sd.t[2u]);

        // dk4/dt1
        dkndt[3u] = dkndt[3u] + h * dkndt[2u];
        dkndt[3u] = qop * mat_helper().column_wise_cross(dkndt[3u], sd.b_last);
        dBdt_tmp = dBdr[3u] * (h * I33 + h2 * 0.5f * dkndt[2u]);
        dkndt[3u] = dkndt[3u] -
                    qop * mat_helper().column_wise_cross(dBdt_tmp, sd.t[3u]);

        /*-----------------------------------------------------------------
         * Calculate all terms of dk_n/dr1
        -------------------------------------------------------------------*/
        // dk1/dr1
        dkndr[0u] = -qop * mat_helper().column_wise_cross(dBdr[0u], sd.t[0u]);

        // dk2/dr1
        dkndr[1u] = qop * mat_helper().column_wise_cross(half_h * dkndr[0u],
                                                         sd.b_middle);
        dBdr_tmp = dBdr[1u] * (I33 + h2 * 0.125 * dkndr[0u]);
        dkndr[1u] = dkndr[1u] -
                    qop * mat_helper().column_wise_cross(dBdr_tmp, sd.t[1u]);

        // dk3/dr1
        dkndr[2u] = qop * mat_helper().column_wise_cross(half_h * dkndr[1u],
                                                         sd.b_middle);
        dBdr_tmp = dBdr[2u] * (I33 + h2 * 0.125 * dkndr[0u]);
        dkndr[2u] = dkndr[2u] -
                    qop * mat_helper().column_wise_cross(dBdr_tmp, sd.t[2u]);

        // dk4/dr1
        dkndr[3u] =
            qop * mat_helper().column_wise_cross(h * dkndr[2u], sd.b_last);
        dBdr_tmp = dBdr[3u] * (I33 + h2 * 0.5 * dkndr[2u]);
        dkndr[3u] = dkndr[3u] -
                    qop * mat_helper().column_wise_cross(dBdr_tmp, sd.t[3u]);

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

    this->_jac_transport = D * this->_jac_transport;
}

template <typename magnetic_field_t, typename algebra_t, typename constraint_t,
          typename policy_t, typename inspector_t,
          template <typename, std::size_t> class array_t>
DETRAY_HOST_DEVICE void detray::rk_stepper<
    magnetic_field_t, algebra_t, constraint_t, policy_t, inspector_t,
    array_t>::state::advance_jacobian(const qop_data& qd,
                                      const detray::stepping::config& cfg) {

    auto D = matrix_operator().template identity<e_free_size, e_free_size>();

    const auto& sd = this->_step_data;
    const scalar_type h{this->_step_size};
    auto& track = this->_track;

    // Half step length
    const scalar_type h2{h * h};
    const scalar_type half_h{h * 0.5f};
    const scalar_type h_6{h * static_cast<scalar_type>(1. / 6.)};

    // 3X3 Identity matrix
    const matrix_type<3, 3> I33 = matrix_operator().template identity<3, 3>();

    // Initialize derivatives
    std::array<matrix_type<3u, 3u>, 4u> dkndt{I33, I33, I33, I33};
    std::array<matrix_type<3u, 3u>, 4u> dkndr;

    if (cfg.use_eloss_gradient) {

        std::array<scalar_type, 4u> dqopn_dqop{1.f, 1.f, 1.f, 1.f};

        // Pre-calculate dqop_n/dqop1
        const scalar_type d2qop1dsdqop1 = this->d2qopdsdqop(qd.qop[0u]);

        dqopn_dqop[0u] = 1.f;
        dqopn_dqop[1u] = 1.f + half_h * d2qop1dsdqop1;

        const scalar_type d2qop2dsdqop1 =
            this->d2qopdsdqop(qd.qop[1u]) * dqopn_dqop[1u];
        dqopn_dqop[2u] = 1.f + half_h * d2qop2dsdqop1;

        const scalar_type d2qop3dsdqop1 =
            this->d2qopdsdqop(qd.qop[2u]) * dqopn_dqop[2u];
        dqopn_dqop[3u] = 1.f + h * d2qop3dsdqop1;

        const scalar_type d2qop4dsdqop1 =
            this->d2qopdsdqop(qd.qop[3u]) * dqopn_dqop[3u];

        // Calculate dkndqop
        std::array<vector3_type, 4u> dkndqop;

        /*-----------------------------------------------------------------
         * Calculate the first and second terms of dk_n/dqop1
         * NOTE: We ignore the field gradient term
        -------------------------------------------------------------------*/

        // dk1/dqop1
        dkndqop[0u] = dqopn_dqop[0u] * vector::cross(sd.t[0u], sd.b_first);

        // dk2/dqop1
        dkndqop[1u] =
            dqopn_dqop[1u] * vector::cross(sd.t[1u], sd.b_middle) +
            qd.qop[1u] * half_h * vector::cross(dkndqop[0u], sd.b_middle);

        // dk3/dqop1
        dkndqop[2u] =
            dqopn_dqop[2u] * vector::cross(sd.t[2u], sd.b_middle) +
            qd.qop[2u] * half_h * vector::cross(dkndqop[1u], sd.b_middle);

        // dk4/dqop1
        dkndqop[3u] = dqopn_dqop[3u] * vector::cross(sd.t[3u], sd.b_last) +
                      qd.qop[3u] * h * vector::cross(dkndqop[2u], sd.b_last);

        // Set dF/dqop1 and dG/dqop1
        vector3_type dFdqop =
            h * h_6 * (dkndqop[0u] + dkndqop[1u] + dkndqop[2u]);
        vector3_type dGdqop =
            h_6 *
            (dkndqop[0u] + 2.f * (dkndqop[1u] + dkndqop[2u]) + dkndqop[3u]);
        matrix_operator().set_block(D, dFdqop, 0u, 7u);
        matrix_operator().set_block(D, dGdqop, 4u, 7u);

        /*-----------------------------------------------------------------
         * Calculate the first terms of d(dqop_n/ds)/dqop1
        -------------------------------------------------------------------*/

        getter::element(D, e_free_qoverp, e_free_qoverp) =
            1.f + h_6 * (d2qop1dsdqop1 + 2.f * (d2qop2dsdqop1 + d2qop3dsdqop1) +
                         d2qop4dsdqop1);
    }

    // Calculate in the case of not considering B field gradient
    if (!cfg.use_field_gradient) {

        /*-----------------------------------------------------------------
         * Calculate the first terms of dk_n/dt1
        -------------------------------------------------------------------*/
        // dk1/dt1
        dkndt[0u] =
            qd.qop[0u] * mat_helper().column_wise_cross(dkndt[0u], sd.b_first);

        // dk2/dt1
        dkndt[1u] = dkndt[1u] + half_h * dkndt[0u];
        dkndt[1u] =
            qd.qop[1u] * mat_helper().column_wise_cross(dkndt[1u], sd.b_middle);

        // dk3/dt1
        dkndt[2u] = dkndt[2u] + half_h * dkndt[1u];
        dkndt[2u] =
            qd.qop[2u] * mat_helper().column_wise_cross(dkndt[2u], sd.b_middle);

        // dk4/dt1
        dkndt[3u] = dkndt[3u] + h * dkndt[2u];
        dkndt[3u] =
            qd.qop[3u] * mat_helper().column_wise_cross(dkndt[3u], sd.b_last);

    } else {

        // Positions at four stages
        std::array<vector3_type, 4u> r;
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
            qd.qop[0u] * mat_helper().column_wise_cross(dkndt[0u], sd.b_first);

        // dk2/dt1
        dkndt[1u] = dkndt[1u] + half_h * dkndt[0u];
        dkndt[1u] =
            qd.qop[1u] * mat_helper().column_wise_cross(dkndt[1u], sd.b_middle);
        dBdt_tmp = dBdr[1u] * (half_h * I33 + h2 * 0.125f * dkndt[0u]);
        dkndt[1u] = dkndt[1u] - qd.qop[1u] * mat_helper().column_wise_cross(
                                                 dBdt_tmp, sd.t[1u]);

        // dk3/dt1
        dkndt[2u] = dkndt[2u] + half_h * dkndt[1u];
        dkndt[2u] =
            qd.qop[2u] * mat_helper().column_wise_cross(dkndt[2u], sd.b_middle);
        dBdt_tmp = dBdr[2u] * (half_h * I33 + h2 * 0.125f * dkndt[0u]);
        dkndt[2u] = dkndt[2u] - qd.qop[2u] * mat_helper().column_wise_cross(
                                                 dBdt_tmp, sd.t[2u]);

        // dk4/dt1
        dkndt[3u] = dkndt[3u] + h * dkndt[2u];
        dkndt[3u] =
            qd.qop[3u] * mat_helper().column_wise_cross(dkndt[3u], sd.b_last);
        dBdt_tmp = dBdr[3u] * (h * I33 + h2 * 0.5f * dkndt[2u]);
        dkndt[3u] = dkndt[3u] - qd.qop[3u] * mat_helper().column_wise_cross(
                                                 dBdt_tmp, sd.t[3u]);

        /*-----------------------------------------------------------------
         * Calculate all terms of dk_n/dr1
        -------------------------------------------------------------------*/
        // dk1/dr1
        dkndr[0u] =
            -qd.qop[0u] * mat_helper().column_wise_cross(dBdr[0u], sd.t[0u]);

        // dk2/dr1
        dkndr[1u] = qd.qop[1u] * mat_helper().column_wise_cross(
                                     half_h * dkndr[0u], sd.b_middle);
        dBdr_tmp = dBdr[1u] * (I33 + h2 * 0.125 * dkndr[0u]);
        dkndr[1u] = dkndr[1u] - qd.qop[1u] * mat_helper().column_wise_cross(
                                                 dBdr_tmp, sd.t[1u]);

        // dk3/dr1
        dkndr[2u] = qd.qop[2u] * mat_helper().column_wise_cross(
                                     half_h * dkndr[1u], sd.b_middle);
        dBdr_tmp = dBdr[2u] * (I33 + h2 * 0.125 * dkndr[0u]);
        dkndr[2u] = dkndr[2u] - qd.qop[2u] * mat_helper().column_wise_cross(
                                                 dBdr_tmp, sd.t[2u]);

        // dk4/dr1
        dkndr[3u] = qd.qop[3u] *
                    mat_helper().column_wise_cross(h * dkndr[2u], sd.b_last);
        dBdr_tmp = dBdr[3u] * (I33 + h2 * 0.5 * dkndr[2u]);
        dkndr[3u] = dkndr[3u] - qd.qop[3u] * mat_helper().column_wise_cross(
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

    this->_jac_transport = D * this->_jac_transport;
}

template <typename magnetic_field_t, typename algebra_t, typename constraint_t,
          typename policy_t, typename inspector_t,
          template <typename, std::size_t> class array_t>
DETRAY_HOST_DEVICE auto detray::rk_stepper<
    magnetic_field_t, algebra_t, constraint_t, policy_t, inspector_t,
    array_t>::state::evaluate_dqopds(const std::size_t i, const scalar_type h,
                                     const scalar_type dqopds_prev,
                                     const detray::stepping::config& cfg,
                                     qop_data& qd) -> scalar_type {

    const auto& track = this->_track;
    const scalar_type qop = track.qop();

    if (cfg.use_mean_loss) {
        if (i == 0u) {
            qd.qop[i] = qop;
        } else {

            // qop_n is calculated recursively like the direction of
            // evaluate_dtds.
            //
            // https://doi.org/10.1016/0029-554X(81)90063-X says:
            // "For y  we  have  similar  formulae  as  for x, for y' and
            // \lambda similar  formulae as for  x'"
            qd.qop[i] = qop + h * dqopds_prev;
        }
    }
    return this->dqopds(qd.qop[i]);
}

template <typename magnetic_field_t, typename algebra_t, typename constraint_t,
          typename policy_t, typename inspector_t,
          template <typename, std::size_t> class array_t>
DETRAY_HOST_DEVICE auto detray::rk_stepper<
    magnetic_field_t, algebra_t, constraint_t, policy_t, inspector_t,
    array_t>::state::evaluate_dtds(const vector3_type& b_field,
                                   const std::size_t i, const scalar_type h,
                                   const vector3_type& dtds_prev,
                                   const scalar_type qop) -> vector3_type {
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

template <typename magnetic_field_t, typename algebra_t, typename constraint_t,
          typename policy_t, typename inspector_t,
          template <typename, std::size_t> class array_t>
DETRAY_HOST_DEVICE auto detray::rk_stepper<
    magnetic_field_t, algebra_t, constraint_t, policy_t, inspector_t,
    array_t>::state::evaluate_field_gradient(const point3_type& pos)
    -> matrix_type<3, 3> {

    matrix_type<3, 3> dBdr = matrix_operator().template zero<3, 3>();

    constexpr auto delta{1e-1f * unit<scalar_type>::mm};

    for (unsigned int i = 0; i < 3; i++) {

        point3_type dpos1 = pos;
        dpos1[i] += delta;
        const auto bvec1_tmp =
            this->_magnetic_field.at(dpos1[0], dpos1[1], dpos1[2]);
        vector3_type bvec1;
        bvec1[0u] = bvec1_tmp[0u];
        bvec1[1u] = bvec1_tmp[1u];
        bvec1[2u] = bvec1_tmp[2u];

        point3_type dpos2 = pos;
        dpos2[i] -= delta;
        const auto bvec2_tmp =
            this->_magnetic_field.at(dpos2[0], dpos2[1], dpos2[2]);
        vector3_type bvec2;
        bvec2[0u] = bvec2_tmp[0u];
        bvec2[1u] = bvec2_tmp[1u];
        bvec2[2u] = bvec2_tmp[2u];

        const vector3_type gradient = (bvec1 - bvec2) * (1.f / (2.f * delta));

        getter::element(dBdr, 0u, i) = gradient[0u];
        getter::element(dBdr, 1u, i) = gradient[1u];
        getter::element(dBdr, 2u, i) = gradient[2u];
    }

    return dBdr;
}

template <typename magnetic_field_t, typename algebra_t, typename constraint_t,
          typename policy_t, typename inspector_t,
          template <typename, std::size_t> class array_t>
DETRAY_HOST_DEVICE auto
detray::rk_stepper<magnetic_field_t, algebra_t, constraint_t, policy_t,
                   inspector_t, array_t>::state::dtds() const -> vector3_type {

    // In case there was no step before
    if (this->_path_length == 0.f) {
        const point3_type pos = this->_track.pos();

        const auto bvec_tmp = this->_magnetic_field.at(pos[0], pos[1], pos[2]);
        vector3_type bvec;
        bvec[0u] = bvec_tmp[0u];
        bvec[1u] = bvec_tmp[1u];
        bvec[2u] = bvec_tmp[2u];

        return this->_track.qop() * vector::cross(this->_track.dir(), bvec);
    }
    return this->_step_data.dtds[3u];
}

template <typename magnetic_field_t, typename algebra_t, typename constraint_t,
          typename policy_t, typename inspector_t,
          template <typename, std::size_t> class array_t>
DETRAY_HOST_DEVICE auto
detray::rk_stepper<magnetic_field_t, algebra_t, constraint_t, policy_t,
                   inspector_t, array_t>::state::dqopds() const -> scalar_type {

    // In case there was no step before
    if (this->_path_length == 0.f) {
        return this->dqopds(this->_track.qop());
    }

    return this->_step_data.dqopds_last;
}

template <typename magnetic_field_t, typename algebra_t, typename constraint_t,
          typename policy_t, typename inspector_t,
          template <typename, std::size_t> class array_t>
DETRAY_HOST_DEVICE auto detray::rk_stepper<
    magnetic_field_t, algebra_t, constraint_t, policy_t, inspector_t,
    array_t>::state::dqopds(const scalar_type qop) const -> scalar_type {

    // d(qop)ds is zero for empty space
    if (this->_mat == nullptr) {
        return 0.f;
    }

    const auto& mat = this->volume_material();

    const scalar_type q = this->_ptc.charge();
    const scalar_type p = q / qop;
    const scalar_type mass = this->_ptc.mass();
    const scalar_type E = math::sqrt(p * p + mass * mass);

    // Compute stopping power
    const scalar_type stopping_power =
        interaction<scalar_type>().compute_stopping_power(mat, this->_ptc,
                                                          {mass, qop, q});

    // Assert that a momentum is a positive value
    assert(p >= 0.f);

    // d(qop)ds, which is equal to (qop) * E * (-dE/ds) / p^2
    // or equal to (qop)^3 * E * (-dE/ds) / q^2
    return qop * qop * qop * E * stopping_power / (q * q);
}

template <typename magnetic_field_t, typename algebra_t, typename constraint_t,
          typename policy_t, typename inspector_t,
          template <typename, std::size_t> class array_t>
DETRAY_HOST_DEVICE auto detray::rk_stepper<
    magnetic_field_t, algebra_t, constraint_t, policy_t, inspector_t,
    array_t>::state::d2qopdsdqop(const scalar_type qop) const -> scalar_type {

    if (this->_mat == nullptr) {
        return 0.f;
    }

    const auto& mat = this->volume_material();
    const scalar_type q = this->_ptc.charge();
    const scalar_type p = q / qop;
    const scalar_type p2 = p * p;

    const auto& mass = this->_ptc.mass();
    const scalar_type E2 = p2 + mass * mass;

    // Interaction object
    interaction<scalar_type> I;

    // g = dE/ds = -1 * (-dE/ds) = -1 * stopping power
    const detail::relativistic_quantities<scalar_type> rq(mass, qop, q);
    const scalar_type g = -1.f * I.compute_stopping_power(mat, this->_ptc, rq);

    // dg/d(qop) = -1 * derivation of stopping power
    const scalar_type dgdqop =
        -1.f * I.derive_stopping_power(mat, this->_ptc, rq);

    // d(qop)/ds = - qop^3 * E * g / q^2
    const scalar_type dqopds = this->dqopds(qop);

    // Check Eq 3.12 of
    // (https://iopscience.iop.org/article/10.1088/1748-0221/4/04/P04016/meta)
    return dqopds * (1.f / qop * (3.f - p2 / E2) + 1.f / g * dgdqop);
}

template <typename magnetic_field_t, typename algebra_t, typename constraint_t,
          typename policy_t, typename inspector_t,
          template <typename, std::size_t> class array_t>
template <typename propagation_state_t>
DETRAY_HOST_DEVICE bool detray::rk_stepper<
    magnetic_field_t, algebra_t, constraint_t, policy_t, inspector_t,
    array_t>::step(propagation_state_t& propagation,
                   const detray::stepping::config& cfg) const {

    // Get stepper and navigator states
    state& stepping = propagation._stepping;
    auto& magnetic_field = stepping._magnetic_field;
    auto& navigation = propagation._navigation;

    if (stepping._step_size == 0.f) {
        stepping._step_size = cfg.min_stepsize;
    } else if (stepping._step_size > 0) {
        stepping._step_size = math::min(stepping._step_size, navigation());
    } else {
        stepping._step_size = math::max(stepping._step_size, navigation());
    }

    const point3_type pos = stepping().pos();

    auto vol = navigation.get_volume();
    if (vol.has_material()) {
        stepping._mat = vol.material_parameters(pos);
    } else {
        stepping._mat = nullptr;
    }

    auto& sd = stepping._step_data;

    // First Runge-Kutta point
    const auto bvec = magnetic_field.at(pos[0], pos[1], pos[2]);
    sd.b_first[0] = bvec[0];
    sd.b_first[1] = bvec[1];
    sd.b_first[2] = bvec[2];

    scalar_type error_estimate{0.f};
    scalar_type error{1e20f};

    // Calculate for vacuum
    if (stepping._mat == nullptr) {

        sd.dqopds_last = 0.f;

        const auto estimate_error = [&](const scalar_type& h) -> scalar {
            // qop
            const scalar_type qop = stepping._track.qop();

            // State the square and half of the step size
            const scalar_type h2{h * h};
            const scalar_type half_h{h * 0.5f};

            // Second Runge-Kutta point
            const point3_type pos1 =
                pos + half_h * sd.t[0u] + h2 * 0.125f * sd.dtds[0u];
            const auto bvec1 = magnetic_field.at(pos1[0], pos1[1], pos1[2]);
            sd.b_middle[0] = bvec1[0];
            sd.b_middle[1] = bvec1[1];
            sd.b_middle[2] = bvec1[2];

            sd.dtds[1u] = stepping.evaluate_dtds(sd.b_middle, 1u, half_h,
                                                 sd.dtds[0u], qop);

            sd.dtds[2u] = stepping.evaluate_dtds(sd.b_middle, 2u, half_h,
                                                 sd.dtds[1u], qop);

            // Last Runge-Kutta point
            const point3_type pos2 =
                pos + h * sd.t[0u] + h2 * 0.5f * sd.dtds[2u];
            const auto bvec2 = magnetic_field.at(pos2[0], pos2[1], pos2[2]);
            sd.b_last[0] = bvec2[0];
            sd.b_last[1] = bvec2[1];
            sd.b_last[2] = bvec2[2];

            sd.dtds[3u] =
                stepping.evaluate_dtds(sd.b_last, 3u, h, sd.dtds[2u], qop);

            // Compute and check the local integration error estimate
            // @Todo
            constexpr const scalar_type one_sixth{
                static_cast<scalar_type>(1. / 6.)};
            const vector3_type err_vec =
                one_sixth * h2 *
                (sd.dtds[0u] - sd.dtds[1u] - sd.dtds[2u] + sd.dtds[3u]);
            error_estimate = getter::norm(err_vec);

            return error_estimate;
        };

        error = stepping.scale_step_size(estimate_error, cfg);
        stepping.pre_step_update(cfg);

        // Advance track state
        stepping.advance_track();

        // Advance jacobian transport
        if (cfg.do_covariance_transport) {
            stepping.advance_jacobian(cfg);
        }

        // Call navigation update policy
        typename rk_stepper::policy_type{}(stepping.policy_state(),
                                           propagation);

        stepping.post_step_update(error, cfg);

    }
    // Calculate for volume material
    else {
        typename state::qop_data qd;
        const scalar_type qop = stepping._track.qop();
        qd.qop[0u] = qop;

        // qop should be recalcuated at every point
        // Reference: Eq (84) of https://doi.org/10.1016/0029-554X(81)90063-X
        const auto estimate_error = [&](const scalar_type& h) -> scalar {
            // State the square and half of the step size
            const scalar_type h2{h * h};
            const scalar_type half_h{h * 0.5f};

            // Second Runge-Kutta point
            // qop should be recalcuated at every point
            // Eq (84) of https://doi.org/10.1016/0029-554X(81)90063-X
            const point3_type pos1 =
                pos + half_h * sd.t[0u] + h2 * 0.125f * sd.dtds[0u];
            const auto bvec1 = magnetic_field.at(pos1[0], pos1[1], pos1[2]);
            sd.b_middle[0] = bvec1[0];
            sd.b_middle[1] = bvec1[1];
            sd.b_middle[2] = bvec1[2];

            qd.dqopds[1u] =
                stepping.evaluate_dqopds(1u, half_h, qd.dqopds[0u], cfg, qd);
            sd.dtds[1u] = stepping.evaluate_dtds(sd.b_middle, 1u, half_h,
                                                 sd.dtds[0u], qd.qop[1u]);

            // Third Runge-Kutta point
            // qop should be recalcuated at every point
            // Reference: Eq (84) of
            // https://doi.org/10.1016/0029-554X(81)90063-X
            qd.dqopds[2u] =
                stepping.evaluate_dqopds(2u, half_h, qd.dqopds[1u], cfg, qd);
            sd.dtds[2u] = stepping.evaluate_dtds(sd.b_middle, 2u, half_h,
                                                 sd.dtds[1u], qd.qop[2u]);

            // Last Runge-Kutta point
            // qop should be recalcuated at every point
            // Eq (84) of https://doi.org/10.1016/0029-554X(81)90063-X
            const point3_type pos2 =
                pos + h * sd.t[0u] + h2 * 0.5f * sd.dtds[2u];
            const auto bvec2 = magnetic_field.at(pos2[0], pos2[1], pos2[2]);
            sd.b_last[0] = bvec2[0];
            sd.b_last[1] = bvec2[1];
            sd.b_last[2] = bvec2[2];

            qd.dqopds[3u] =
                stepping.evaluate_dqopds(3u, h, qd.dqopds[2u], cfg, qd);
            sd.dtds[3u] = stepping.evaluate_dtds(sd.b_last, 3u, h, sd.dtds[2u],
                                                 qd.qop[3u]);

            // Compute and check the local integration error estimate
            // @Todo
            constexpr const scalar_type one_sixth{
                static_cast<scalar_type>(1. / 6.)};
            const vector3_type err_vec =
                one_sixth * h2 *
                (sd.dtds[0u] - sd.dtds[1u] - sd.dtds[2u] + sd.dtds[3u]);
            error_estimate = getter::norm(err_vec);

            return error_estimate;
        };

        error = stepping.scale_step_size(estimate_error, cfg);
        stepping.pre_step_update(cfg);

        // Advance track state
        stepping.advance_track(qd);

        // Advance jacobian transport
        if (cfg.do_covariance_transport) {
            stepping.advance_jacobian(qd, cfg);
        }
    }

    // Call navigation update policy
    typename rk_stepper::policy_type{}(stepping.policy_state(), propagation);

    stepping.post_step_update(error, cfg);

    // Run final inspection
    stepping.run_inspector(cfg, "Step complete: ");

    return true;
}
