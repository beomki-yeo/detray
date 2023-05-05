/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "detray/masks/detail/mask_property.hpp"

// Google test include(s).
#include <gtest/gtest.h>

#include <iostream>
#include <type_traits>
using namespace detray;

class A {};
class B : A {};

template <typename T>
struct test {};

// Test mask propery
TEST(masks, mask_property) {

    passive a;
    (void)a;
    sensitive<2u, false> b;
    (void)b;

    static_assert(std::is_base_of_v<mask_property, passive> == true);
    // mask_property<true, 2u, false> b;
}