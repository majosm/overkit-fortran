# Copyright (c) 2017 Matthew J. Smith and Overkit contributors
# License: MIT (http://opensource.org/licenses/MIT)

if(NOT DEFINED C_ENABLED)
  message(FATAL_ERROR "Must set C_ENABLED for ReleaseBuildType")
endif()
if(NOT DEFINED CXX_ENABLED)
  message(FATAL_ERROR "Must set CXX_ENABLED for ReleaseBuildType")
endif()
if(NOT DEFINED Fortran_ENABLED)
  message(FATAL_ERROR "Must set Fortran_ENABLED for ReleaseBuildType")
endif()

if(${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")

  if (C_ENABLED)
    set(CMAKE_C_DIALECT "-std=c99")
    set(CMAKE_C_OPTS "")
  endif()

  if (CXX_ENABLED)
    set(CMAKE_CXX_DIALECT "-std=c++03")
    set(CMAKE_CXX_OPTS "")
  endif()

  if (Fortran_ENABLED)
    set(CMAKE_Fortran_DIALECT "-std=gnu -ffree-line-length-0")
    set(CMAKE_Fortran_OPTS "-fomit-frame-pointer -ffast-math -funroll-loops")
  endif()

elseif(${CMAKE_Fortran_COMPILER_ID} MATCHES "Intel")

  if(C_ENABLED)
    set(CMAKE_C_DIALECT "")
    set(CMAKE_C_OPTS "")
  endif()

  if(CXX_ENABLED)
    set(CMAKE_CXX_DIALECT "")
    set(CMAKE_CXX_OPTS "")
  endif()

  if(Fortran_ENABLED)
    set(CMAKE_Fortran_DIALECT "")
#    set(CMAKE_Fortran_OPTS "-ip -ipo -xHost")
    set(CMAKE_Fortran_OPTS "-ip -xHost")
  endif()

else()

  if(C_ENABLED)
    set(CMAKE_C_DIALECT "")
    set(CMAKE_C_OPTS "")
  endif()

  if(CXX_ENABLED)
    set(CMAKE_CXX_DIALECT "")
    set(CMAKE_CXX_OPTS "")
  endif()

  if(Fortran_ENABLED)
    set(CMAKE_Fortran_DIALECT "")
    set(CMAKE_Fortran_OPTS "")
  endif()

endif()

if(C_ENABLED)
  set(CMAKE_C_FLAGS_RELEASE_EXTRA "${CMAKE_C_DIALECT} ${CMAKE_C_OPTS}")
endif()
if(CXX_ENABLED)
  set(CMAKE_CXX_FLAGS_RELEASE_EXTRA "${CMAKE_CXX_DIALECT} ${CMAKE_CXX_OPTS}")
endif()
if(Fortran_ENABLED)
  set(CMAKE_Fortran_FLAGS_RELEASE_EXTRA "${CMAKE_Fortran_DIALECT} ${CMAKE_Fortran_OPTS}")
endif()

if(NOT DEFINED RELEASEBUILDTYPE)
  if(C_ENABLED)
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${CMAKE_C_FLAGS_RELEASE_EXTRA}"
      CACHE STRING "Flags used by the C compiler during release builds." FORCE)
  endif()
  if(CXX_ENABLED)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${CMAKE_CXX_FLAGS_RELEASE_EXTRA}"
      CACHE STRING "Flags used by the C++ compiler during release builds." FORCE)
  endif()
  if(Fortran_ENABLED)
    set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${CMAKE_Fortran_FLAGS_RELEASE_EXTRA}"
      CACHE STRING "Flags used by the Fortran compiler during release builds." FORCE)
  endif()
  set(RELEASEBUILDTYPE TRUE CACHE INTERNAL "")
endif()
