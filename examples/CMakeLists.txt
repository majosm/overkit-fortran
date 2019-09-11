# Copyright (c) 2018 Matthew J. Smith and Overkit contributors
# License: MIT (http://opensource.org/licenses/MIT)

#==================
# Overkit examples
#==================

if(EXAMPLES)

# Link to overkit and overkit-support
include_directories("${PROJECT_BINARY_DIR}/src/core/mod")
include_directories("${PROJECT_BINARY_DIR}/src/support/mod")
link_libraries(overkit-support overkit)

add_executable(Blobs Blobs.F90)
add_executable(Bump Bump.F90)
add_executable(Inlet Inlet.F90)
add_executable(Circle Circle.F90)

endif()
