# Copyright (c) 2017 Matthew J. Smith and Overkit contributors
# License: MIT (http://opensource.org/licenses/MIT)

file(READ "${FILE}" FILE_CONTENTS)
string(REGEX REPLACE "#line[^\n]*" "" STRIPPED_FILE_CONTENTS "${FILE_CONTENTS}")
file(WRITE "${FILE}" "${STRIPPED_FILE_CONTENTS}")
