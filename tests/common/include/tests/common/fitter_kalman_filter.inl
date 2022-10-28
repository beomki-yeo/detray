/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "detray/fitters/kalman_filter/kalman_actor.hpp"
#include "detray/fitters/kalman_filter/track_state.hpp"
#include "detray/propagator/actors/aborters.hpp"
#include "detray/propagator/actors/parameter_resetter.hpp"
#include "detray/propagator/actors/parameter_transporter.hpp"
#include "detray/propagator/actors/pointwise_material_interactor.hpp"
#include "detray/propagator/propagator.hpp"
#include "tests/common/tools/create_telescope_detector.hpp"

// Vecmem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

// google-test include(s).
#include <gtest/gtest.h>

using namespace detray;
using transform3 = __plugin::transform3<detray::scalar>;

// measurement type
struct measurement {
    /// Local 2D coordinates for a measurement on a detector module
    std::array<scalar, 2> local{0., 0.};
    /// Variance on the 2D coordinates of the measurement
    std::array<scalar, 2> variance{0., 0.};
};

// @note: Lots of lines are overlapped with ones in check_simulation.inl.
TEST(fitter, kalman_filter_telescope_geometry) {

    // Memory resource
    vecmem::host_memory_resource host_mr;

    // Create B field
    const vector3 B{0, 0, 2 * unit_constants::T};

    // Build telescope geometry from given module positions
    detail::ray<transform3> traj{
        {0, 0, 0}, 0, {0, 0, 1}, -1};  // aligned to z-axis
    std::vector<scalar> positions = {0.,   50., 100., 150., 200., 250.,
                                     300., 350, 400,  450., 500.};

    const auto mat = silicon_tml<scalar>();
    const scalar thickness = 0.17 * unit_constants::cm;

    using b_field_t = decltype(create_telescope_detector(host_mr))::bfield_type;
    const auto det = create_telescope_detector(
        host_mr,
        b_field_t(b_field_t::backend_t::configuration_t{B[0], B[1], B[2]}),
        positions, traj, 1000. * unit_constants::mm, 1000. * unit_constants::mm,
        mat, thickness);

    // Generate simulation data

    // Read simulation data

    // Run Kalman filtering

    // Some tests
}