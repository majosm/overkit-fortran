# Copyright (c) 2017 Matthew J. Smith and Overkit contributors
# License: MIT (http://opensource.org/licenses/MIT)

#==================
# Overkit examples
#==================

if(EXAMPLES)

# Link to overkit and overkit-support
include_directories("${CMAKE_BINARY_DIR}/src/core/mod")
include_directories("${CMAKE_BINARY_DIR}/src/support/mod")
link_libraries(overkit-support overkit)

add_executable(Blobs Blobs.F90)
#add_executable(Airfoil Airfoil.F90)
#add_executable(Sphere Sphere.F90)

endif()
