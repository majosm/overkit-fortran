# Copyright (c) 2018 Matthew J. Smith and Overkit contributors
# License: MIT (http://opensource.org/licenses/MIT)

#=======================
# Overkit library tests
#=======================

if(TESTS)

  # Unit tests
  if(PFUNIT_FOUND)

    set(UNIT_TEST_SOURCES
      CommandLineTest.pf
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
          -P "${PROJECT_SOURCE_DIR}/config/scripts/strip-line.cmake"
        DEPENDS "${TEST_PATH}"
      )
      # Append the test to the test suites file
      set(PROCESSED_TEST_SOURCES ${PROCESSED_TEST_SOURCES} "${PROCESSED_TEST_PATH}")
      file(APPEND "${TEST_SUITES_FILE_PATH}" "ADD_TEST_SUITE(${TEST_NAME}_suite)\n")
    endforeach()

    # Find pFUnit modules, etc.
    include_directories("${PFUNIT_INCLUDE_DIR}" "${PFUNIT_MODULE_DIR}")

    # Find Overkit module file
    include_directories("${PROJECT_BINARY_DIR}/src/core/mod")
    include_directories("${PROJECT_BINARY_DIR}/src/support/mod")

    # Find testSuites.inc file
    include_directories("${CMAKE_CURRENT_BINARY_DIR}")

    add_executable(runtests-support "${PFUNIT_DRIVER}" ${PROCESSED_TEST_SOURCES})
    target_link_libraries(runtests-support overkit-fortran overkit-support "${PFUNIT_LIBRARY}")

    add_test(NAME support-unit-tests
      COMMAND "${CMAKE_COMMAND}" -DTEST_EXECUTABLE=${CMAKE_CURRENT_BINARY_DIR}/runtests-support
        -P "${PROJECT_SOURCE_DIR}/config/scripts/pfunit-passfail.cmake")

  endif()

endif()
