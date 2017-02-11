# Copyright (c) 2017 Matthew J. Smith and Overkit contributors
# License: MIT (http://opensource.org/licenses/MIT)

#===============
# Overkit tests
#===============

# Coverage report
add_custom_target(coverage
  COMMAND "${CMAKE_COMMAND}" -DSOURCE_DIR="${CMAKE_SOURCE_DIR}" -DBINARY_DIR="${CMAKE_BINARY_DIR}"
    -DCOVERAGE_DIR="src/core" -P "${CMAKE_SOURCE_DIR}/config/scripts/gen-coverage-data.cmake"
  COMMAND "${CMAKE_COMMAND}" -DSOURCE_DIR="${CMAKE_SOURCE_DIR}" -DBINARY_DIR="${CMAKE_BINARY_DIR}"
    -DCOVERAGE_DIR="src/support" -P "${CMAKE_SOURCE_DIR}/config/scripts/gen-coverage-data.cmake"
  COMMAND "${CMAKE_COMMAND}" -DSOURCE_DIR="${CMAKE_SOURCE_DIR}" -DBINARY_DIR="${CMAKE_BINARY_DIR}"
    -DCOVERAGE_DIR="src/utils" -P "${CMAKE_SOURCE_DIR}/config/scripts/gen-coverage-data.cmake"
  COMMAND "${CMAKE_COMMAND}" -DSOURCE_DIR="${CMAKE_SOURCE_DIR}" -DBINARY_DIR="${CMAKE_BINARY_DIR}"
    -P "${CMAKE_SOURCE_DIR}/config/scripts/gen-coverage-report.cmake"
)

# Reset coverage data
add_custom_target(coverageclean
  COMMAND "${CMAKE_COMMAND}" -DBINARY_DIR="${CMAKE_BINARY_DIR}"
    -P "${CMAKE_SOURCE_DIR}/config/scripts/coverageclean.cmake"
)

#================
# Subdirectories
#================

add_subdirectory(system)
add_subdirectory(unit)
