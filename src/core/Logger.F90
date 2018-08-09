! Copyright (c) 2018 Matthew J. Smith and Overkit contributors
! License: MIT (http://opensource.org/licenses/MIT)

module ovkLogger

  use ovkGlobal
  implicit none

  private

  ! Internal API
  public :: t_logger
  public :: t_logger_
  public :: operator (==)
  public :: operator (/=)
  public :: EnableStatusLog
  public :: DisableStatusLog
  public :: EnableErrorLog
  public :: DisableErrorLog

  type t_logger
    type(t_noconstruct) :: noconstruct
    logical :: log_status
    logical :: log_errors
    integer :: status_file
    integer :: error_file
  end type t_logger

  ! Trailing _ added for compatibility with compilers that don't support F2003 constructors
  interface t_logger_
    module procedure t_logger_Default
  end interface t_logger_

  interface operator (==)
    module procedure t_logger_Equal
  end interface operator (==)

  interface operator (/=)
    module procedure t_logger_NotEqual
  end interface operator (/=)

contains

  pure function t_logger_Default() result(Logger)

    type(t_logger) :: Logger

    Logger%log_status = .false.
    Logger%log_errors = .false.
    Logger%status_file = -1
    Logger%error_file = -1

  end function t_logger_Default

  pure function t_logger_Equal(LeftLogger, RightLogger) result(Equal)

    type(t_logger), intent(in) :: LeftLogger, RightLogger
    logical :: Equal

    Equal = &
      (LeftLogger%log_status .eqv. RightLogger%log_status) .and. &
      (LeftLogger%log_errors .eqv. RightLogger%log_errors)

    if (Equal) then
      if (LeftLogger%log_status) then
        Equal = Equal .and. LeftLogger%status_file == RightLogger%status_file
      end if
      if (LeftLogger%log_errors) then
        Equal = Equal .and. LeftLogger%error_file == RightLogger%error_file
      end if
    end if

  end function t_logger_Equal

  pure function t_logger_NotEqual(LeftLogger, RightLogger) result(NotEqual)

    type(t_logger), intent(in) :: LeftLogger, RightLogger
    logical :: NotEqual

    NotEqual = .not. t_logger_Equal(LeftLogger, RightLogger)

  end function t_logger_NotEqual

  subroutine EnableStatusLog(Logger, StatusFile)

    type(t_logger), intent(inout) :: Logger
    integer, intent(in) :: StatusFile

    Logger%log_status = .true.
    Logger%status_file = StatusFile

  end subroutine EnableStatusLog

  subroutine DisableStatusLog(Logger)

    type(t_logger), intent(inout) :: Logger

    Logger%log_status = .false.
    Logger%status_file = -1

  end subroutine DisableStatusLog

  subroutine EnableErrorLog(Logger, ErrorFile)

    type(t_logger), intent(inout) :: Logger
    integer, intent(in) :: ErrorFile

    Logger%log_errors = .true.
    Logger%error_file = ErrorFile

  end subroutine EnableErrorLog

  subroutine DisableErrorLog(Logger)

    type(t_logger), intent(inout) :: Logger

    Logger%log_errors = .false.
    Logger%error_file = -1

  end subroutine DisableErrorLog

end module ovkLogger
