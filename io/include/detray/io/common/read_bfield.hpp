/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/io/common/detail/file_handle.hpp"

// Covfie include(s)
#include <covfie/core/utility/binary_io.hpp>

// System include(s)
#include <ios>
#include <string>

namespace detray::io {

/// @brief Function that reads the first 4 bytes of a potential bfield file and
/// checks that it contains data for a covfie field
inline bool check_covfie_file(io::detail::file_handle& file) {

    // See "covfie/lib/core/utility/binary_io.hpp"
    std::uint32_t hdr = covfie::utility::read_binary<std::uint32_t>(*file);

    // Compare to magic bytes
    return (hdr == covfie::utility::MAGIC_HEADER);
}

/// @brief function that reads a covfie field from file
template <typename bfield_t>
bfield_t read_bfield(const std::string& file_name) {

    // Open binary file
    io::detail::file_handle file{file_name,
                                 std::ios_base::in | std::ios_base::binary};

    check_covfie_file(file);

    return bfield_t(*file);
}

}  // namespace detray::io
