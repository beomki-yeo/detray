# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

# Set the CMake minimum version.
cmake_minimum_required( VERSION 3.13 )

# Set up the library.
detray_add_library( detray_io_csv io_csv
  "include/detray/io/csv_io_types.hpp"
  "include/detray/io/csv_io.hpp"
  "include/detray/io/utils.hpp" )
target_link_libraries( detray_io_csv INTERFACE dfelibs::dfelibs vecmem::core )

# Make the "main I/O library" link against it.
target_link_libraries( detray_io INTERFACE detray::io_csv )
