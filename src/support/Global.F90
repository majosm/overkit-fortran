! Copyright (c) 2017 Matthew J. Smith and Overkit contributors
! License: MIT (http://opensource.org/licenses/MIT)

module ovsGlobal

#ifdef f2003
  use, intrinsic iso_fortran_env, only : INPUT_UNIT, OUTPUT_UNIT, ERROR_UNIT
#endif

  implicit none

  private

  public :: rk
  public :: lk
  public :: bk
  public :: INPUT_UNIT, OUTPUT_UNIT, ERROR_UNIT
  public :: DEBUG
  public :: MAX_DIMS
  public :: PATH_LENGTH
  public :: Pi

  integer, parameter :: rk = selected_real_kind(15, 307)
  integer, parameter :: lk = selected_int_kind(18)
  integer, parameter :: bk = selected_int_kind(1)

#ifndef f2003
  integer, parameter :: INPUT_UNIT = 5
  integer, parameter :: OUTPUT_UNIT = 6
  integer, parameter :: ERROR_UNIT = 0
#endif

#ifdef OVERKIT_DEBUG
  logical, parameter :: DEBUG = .true.
#else
  logical, parameter :: DEBUG = .false.
#endif

  integer, parameter :: MAX_DIMS = 3

  integer, parameter :: PATH_LENGTH = 256

  real(rk), parameter :: Pi = 3.1415926535897932385_rk

end module ovsGlobal
