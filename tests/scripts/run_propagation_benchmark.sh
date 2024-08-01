#!/bin/bash
#
# Detray library, part of the ACTS project (R&D line)
#
# (c) 2021-2024 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

echo "===> CI Benchmark running script for detray"

# Set the script path
BASEDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Set workspace directory
if [ -z "${GITHUB_WORKSPACE}" ]; then
    WORKSPACE=${BASEDIR}/../../
else
    WORKSPACE=${GITHUB_WORKSPACE}
fi

git config --global --add safe.directory $(pwd)
