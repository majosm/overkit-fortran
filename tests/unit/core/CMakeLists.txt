# Copyright (c) 2017 Matthew J. Smith and Overkit contributors
# License: MIT (http://opensource.org/licenses/MIT)

#============
# Unit tests
#============

if(TESTS)

  find_package(pFUnit)

  if(PFUNIT_FOUND)

    set(UNIT_TEST_SOURCES
      BoundingBoxTest.pf
      CartTest.pf
      DonorAccelTest.pf
      DonorsTest.pf
      FieldTest.pf
      GeometryTest.pf
      GridTest.pf
      HashGridTest.pf
      InterpTest.pf
      MaskTest.pf
      OversetTest.pf
      PegasusTest.pf
    )

    # Create the test suites file
    set(TEST_SUITES_FILE_PATH "${CMAKE_CURRENT_BINARY_DIR}/testSuites.inc")
    file(WRITE "${TEST_SUITES_FILE_PATH}" "")

    foreach(TEST_FILE ${UNIT_TEST_SOURCES})
      # Process the .pf file into .F90 and place it in the build directory
      get_filename_component(TEST_NAME ${TEST_FILE} NAME_WE)
      set(TEST_PATH "${CMAKE_CURRENT_SOURCE_DIR}/${TEST_NAME}.pf")
      set(PROCESSED_TEST_PATH "${CMAKE_CURRENT_BINARY_DIR}/${TEST_NAME}.F90")
      add_custom_command(
        OUTPUT "${PROCESSED_TEST_PATH}"
        COMMAND python "${PFUNIT_PARSER}" "${TEST_PATH}" "${PROCESSED_TEST_PATH}"
        # Strip out #line directives (they seem to mess up gfortran's error output)
        COMMAND "${CMAKE_COMMAND}" -DFILE="${PROCESSED_TEST_PATH}"
          -P "${CMAKE_SOURCE_DIR}/config/scripts/strip-line.cmake"
        DEPENDS "${TEST_PATH}"
      )
      # Append the test to the test suites file
      set(PROCESSED_TEST_SOURCES ${PROCESSED_TEST_SOURCES} "${PROCESSED_TEST_PATH}")
      file(APPEND "${TEST_SUITES_FILE_PATH}" "ADD_TEST_SUITE(${TEST_NAME}_suite)\n")
    endforeach()

    # Find pFUnit modules, etc.
    include_directories("${PFUNIT_INCLUDE_DIR}" "${PFUNIT_MODULE_DIR}")

    # Find Overkit module file
    include_directories("${CMAKE_BINARY_DIR}/src/core/mod")
    include_directories("${CMAKE_BINARY_DIR}/src/support/mod")

    # Find testSuites.inc file
    include_directories("${CMAKE_CURRENT_BINARY_DIR}")

    add_executable(runall "${PFUNIT_DRIVER}" ${PROCESSED_TEST_SOURCES})
    target_link_libraries(runall overkit overkit-support "${PFUNIT_LIBRARY}")

    add_test(unit-tests "${CMAKE_CURRENT_BINARY_DIR}/runall" -v)

  else()

    message(WARNING "Could not find pFUnit. Unit tests will not be built.")
    set(GLOBAL_WARNING_FLAG TRUE CACHE INTERNAL "")

  endif()

endif()
