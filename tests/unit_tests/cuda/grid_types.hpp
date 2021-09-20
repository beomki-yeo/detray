/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include "tests/common/test_defs.hpp"
#include "grids/axis.hpp"
#include "grids/grid2.hpp"
#include "grids/serializer2.hpp"
#include "grids/populator.hpp"
#include "utils/indexing.hpp"

#pragma once

namespace detray{

    using populator_t = complete_populator<3, false, test::point3>;
    using grid2r_complete = grid2<complete_populator<3, false, test::point3>,
				  axis::regular<>,
				  axis::regular<>,
				  serializer2>;

    using grid2r_complete_data = grid2_data<complete_populator<3, false, test::point3>,
					    axis::regular<>,
					    axis::regular<>,
					    serializer2>;
    
} // namespace
