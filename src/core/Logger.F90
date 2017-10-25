! Copyright (c) 2017 Matthew J. Smith and Overkit contributors
! License: MIT (http://opensource.org/licenses/MIT)

module ovkLogger

  use ovkGlobal
  implicit none

  private

  ! Internal
  public :: t_logger
  public :: t_logger_

  type t_logger
    type(t_noconstruct) :: noconstruct
    logical :: verbose
  end type t_logger

  ! Trailing _ added for compatibility with compilers that don't support F2003 constructors
  interface t_logger_
    module procedure t_logger_Default
    module procedure t_logger_Assigned
  end interface t_logger_

contains

  function t_logger_Default() result(Logger)

    type(t_logger) :: Logger

    Logger%verbose = .false.

  end function t_logger_Default

  function t_logger_Assigned(Verbose) result(Logger)

    type(t_logger) :: Logger
    logical, intent(in) :: Verbose

    Logger%verbose = Verbose

  end function t_logger_Assigned

end module ovkLogger
