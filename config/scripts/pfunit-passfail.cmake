# Copyright (c) 2017 Matthew J. Smith and Overkit contributors
# License: MIT (http://opensource.org/licenses/MIT)

# pFUnit driver doesn't seem to provide exit codes for pass/fail; this script remedies that
execute_process(COMMAND ${TEST_EXECUTABLE} -v
  OUTPUT_VARIABLE TEST_OUTPUT
  ERROR_VARIABLE TEST_OUTPUT
  RESULT_VARIABLE TEST_RESULT
)
message(STATUS ${TEST_OUTPUT})
if(NOT TEST_RESULT EQUAL 0 OR
  TEST_OUTPUT MATCHES "FAILURE" OR
  TEST_OUTPUT MATCHES "Failure" OR
  TEST_OUTPUT MATCHES "failure")
  message(FATAL_ERROR "")
endif()
