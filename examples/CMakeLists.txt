# Copyright (c) 2017 Matthew J. Smith and Overkit contributors
# License: MIT (http://opensource.org/licenses/MIT)

#==================
# Overkit examples
#==================

if(EXAMPLES)

include_directories("${CMAKE_BINARY_DIR}/src/core")
include_directories("${CMAKE_BINARY_DIR}/src/support")
link_libraries(overkit overkit-support)

add_executable(BoxInBox BoxInBox.F90)
add_executable(Blobs Blobs.F90)
#add_executable(Airfoil Airfoil.F90)
#add_executable(Sphere Sphere.F90)

endif()
