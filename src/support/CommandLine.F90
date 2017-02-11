! Copyright (c) 2017 Matthew J. Smith and Overkit contributors
! License: MIT (http://opensource.org/licenses/MIT)

module ovsCommandLine

  use ovsGlobal
  implicit none

  private

  public :: t_cmd_opt
  public :: t_cmd_opt_
  public :: ParseArguments
  public :: OptionPresent
  public :: GetOptionValue
  public :: CMD_OPT_NONE
  public :: CMD_OPT_INTEGER
  public :: CMD_OPT_REAL
  public :: CMD_OPT_STRING

  integer, parameter :: CMD_ARG_LENGTH = 256

  type t_cmd_opt
    character(len=CMD_ARG_LENGTH) :: long_name
    character(len=1) :: short_name
    integer :: value_type
    character(len=4096) :: description
    logical :: is_present
    character(len=CMD_ARG_LENGTH) :: value
  end type t_cmd_opt

  ! Trailing _ added for compatibility with compilers that don't support F2003 constructors
  interface t_cmd_opt_
    module procedure t_cmd_opt_Default
    module procedure t_cmd_opt_Specified
  end interface t_cmd_opt_

  interface GetOptionValue
    module procedure GetOptionValue_Int
    module procedure GetOptionValue_Int_Default
    module procedure GetOptionValue_Real
    module procedure GetOptionValue_Real_Default
    module procedure GetOptionValue_String
    module procedure GetOptionValue_String_Default
  end interface GetOptionValue

  integer, parameter :: CMD_OPT_NONE = 1
  integer, parameter :: CMD_OPT_INTEGER = 2
  integer, parameter :: CMD_OPT_REAL = 3
  integer, parameter :: CMD_OPT_STRING = 4

contains

  function t_cmd_opt_Default() result(Option)

    type(t_cmd_opt) :: Option

    Option%long_name = ""
    Option%short_name = ""
    Option%value_type = CMD_OPT_NONE
    Option%description = ""
    Option%is_present = .false.
    Option%value = ""

  end function t_cmd_opt_Default

  function t_cmd_opt_Specified(LongName, ShortName, ValueType, Description) result(Option)

    character(len=*), intent(in) :: LongName
    character(len=*), intent(in) :: ShortName
    integer, intent(in) :: ValueType
    character(len=*), intent(in) :: Description
    type(t_cmd_opt) :: Option

    Option%long_name = LongName
    Option%short_name = ShortName
    Option%value_type = ValueType
    Option%description = Description
    Option%is_present = .false.
    Option%value = ""

  end function t_cmd_opt_Specified

  subroutine ParseArguments(RawArguments, Usage, Description, LongDescription, Options, &
    Arguments, MinArguments, MaxArguments)

    character(len=*), dimension(:), intent(in) :: RawArguments
    character(len=*), intent(in), optional :: Usage
    character(len=*), intent(in), optional :: Description
    character(len=*), intent(in), optional :: LongDescription
    type(t_cmd_opt), dimension(:), intent(inout), optional :: Options
    character(len=*), dimension(:), allocatable, intent(out), optional :: Arguments
    integer, intent(in), optional :: MinArguments, MaxArguments

    integer :: i, j, k, l
    logical :: Help
    character(len=CMD_ARG_LENGTH) :: RawArgument, NextRawArgument
    logical :: SkipNext
    character(len=CMD_ARG_LENGTH) :: OptionName, OptionValue
    logical :: ValidOptionValue
    character(len=CMD_ARG_LENGTH), dimension(:), allocatable :: ArgumentsTemp
    integer :: nArguments

    if (present(Arguments)) then
      allocate(ArgumentsTemp(size(RawArguments)))
      nArguments = 0
    end if

    Help = .false.

    i = 1
    argloop: do while (i <= size(RawArguments))

      RawArgument = RawArguments(i)
      if (i < size(RawArguments)) then
        NextRawArgument = RawArguments(i+1)
      else
        NextRawArgument = ""
      end if

      SkipNext = .false.

      if (RawArgument(1:2) == "--") then

        j = index(RawArgument, "=")
        j = merge(j, len_trim(RawArgument)+1, j /= 0)
        OptionName = RawArgument(3:j-1)
        OptionValue = RawArgument(j+1:)

        if (OptionName == "help") then
          Help = .true.
          exit argloop
        end if

        l = 0
        if (present(Options)) then
          do k = 1, size(Options)
            if (Options(k)%long_name == OptionName) then
              l = k
              exit
            end if
          end do
        end if
        if (l == 0) then
          write (ERROR_UNIT, '(3a)') "ERROR: Unrecognized option '", trim(OptionName), "'."
          stop 1
        end if
        if (OptionValue == "" .and. Options(l)%value_type /= CMD_OPT_NONE) then
          if (NextRawArgument /= "") then
            OptionValue = NextRawArgument
            SkipNext = .true.
          else
            write (ERROR_UNIT, '(3a)') "ERROR: Missing value for option '", trim(OptionName), &
              "'."
            stop 1
          end if
        end if
        ValidOptionValue = ValidateOptionValue(OptionValue, Options(l)%value_type)
        if (ValidOptionValue) then
          Options(l)%is_present = .true.
          Options(l)%value = OptionValue
        else
          write (ERROR_UNIT, '(3a)') "ERROR: Unrecognized value for option '", &
            trim(OptionName), "'."
          stop 1
        end if

      else if (RawArgument(1:1) == "-") then

        do j = 2, len_trim(RawArgument)
          OptionName = RawArgument(j:j)
          if (OptionName == "h") then
            Help = .true.
            exit argloop
          end if
          l = 0
          if (present(Options)) then
            do k = 1, size(Options)
              if (Options(k)%short_name == OptionName) then
                l = k
                exit
              end if
            end do
          end if
          if (l == 0) then
            write (ERROR_UNIT, '(3a)') "ERROR: Unrecognized option '", trim(OptionName), "'."
            stop 1
          end if
          if (Options(l)%value_type == CMD_OPT_NONE) then
            Options(l)%is_present = .true.
            Options(l)%value = ""
          else
            OptionValue = RawArgument(j+1:)
            if (OptionValue == "") then
              if (NextRawArgument /= "") then
                OptionValue = NextRawArgument
                SkipNext = .true.
              else
                write (ERROR_UNIT, '(3a)') "ERROR: Missing value for option '", &
                  trim(OptionName), "'."
                stop 1
              end if
            end if
            ValidOptionValue = ValidateOptionValue(OptionValue, Options(l)%value_type)
            if (ValidOptionValue) then
              Options(l)%is_present = .true.
              Options(l)%value = OptionValue
            else
              write (ERROR_UNIT, '(3a)') "ERROR: Unrecognized value for option '", &
                trim(OptionName), "'."
              stop 1
            end if
            exit
          end if
        end do

      else

        if (present(Arguments)) then
          ArgumentsTemp(nArguments+1) = RawArgument
          nArguments = nArguments + 1
        else
          write (ERROR_UNIT, '(3a)') "ERROR: Unrecognized argument '", trim(RawArgument), "'."
          stop 1
        end if

      end if

      i = merge(i+2, i+1, SkipNext)

    end do argloop

    if (Help) then
      if (present(Usage)) then
        write (*, '(2a)') "Usage: ", trim(Usage)
        write (*, '(a)') ""
      end if
      if (present(Description)) then
        write (*, '(2a)') "Description: ", trim(Description)
        write (*, '(a)') ""
      end if
      if (present(LongDescription)) then
        write (*, '(a)') trim(LongDescription)
        write (*, '(a)') ""
      end if
      if (present(Options)) then
        write (*, '(a)') "Options:"
        do i = 1, size(Options)
          write (*, '(6a)') "  --", trim(Options(i)%long_name), " (-", Options(i)%short_name, ") -- ", &
            trim(Options(i)%description)
        end do
        write (*, '(a)') ""
      end if
      stop
    end if

    if (present(Arguments)) then

      if (present(MinArguments)) then
        if (nArguments < MinArguments) then
          write (ERROR_UNIT, '(a)') "ERROR: Not enough command-line arguments."
          stop 1
        end if
      end if

      if (present(MaxArguments)) then
        if (nArguments > MaxArguments) then
          write (ERROR_UNIT, '(a)') "ERROR: Too many command-line arguments."
          stop 1
        end if
      end if

      allocate(Arguments(nArguments))
      Arguments = ArgumentsTemp(:nArguments)

    end if

  end subroutine ParseArguments

  function ValidateOptionValue(OptionValue, OptionValueType) result(ValidOptionValue)

    character(len=CMD_ARG_LENGTH), intent(in) :: OptionValue
    integer, intent(in) :: OptionValueType
    logical :: ValidOptionValue

    integer :: Error
    integer :: IntegerValue
    real(rk) :: RealValue

    select case (OptionValueType)
    case (CMD_OPT_NONE)
      ValidOptionValue = OptionValue == ""
    case (CMD_OPT_INTEGER)
      read (OptionValue, *, iostat=Error) IntegerValue
      ValidOptionValue = Error == 0
    case (CMD_OPT_REAL)
      read (OptionValue, *, iostat=Error) RealValue
      ValidOptionValue = Error == 0
    case (CMD_OPT_STRING)
      ValidOptionValue = OptionValue /= ""
    end select

  end function ValidateOptionValue

  pure function OptionPresent(Option) result(IsPresent)

    type(t_cmd_opt), intent(in) :: Option
    logical :: IsPresent

    IsPresent = Option%is_present

  end function OptionPresent

  subroutine GetOptionValue_Int(Option, Value)

    type(t_cmd_opt), intent(in) :: Option
    integer, intent(out) :: Value

    if (DEBUG) then
      if (.not. Option%is_present) then
        write (ERROR_UNIT, '(3a)') "ERROR: Attempted to read value of command line option '", &
          trim(Option%long_name), "' which is not present."
        stop 1
      end if
      if (Option%value_type == CMD_OPT_NONE) then
        write (ERROR_UNIT, '(3a)') "ERROR: Attempted to read value of command line option '", &
          trim(Option%long_name), "' which has no value type specified."
        stop 1
      end if
    end if

    read (Option%value, *) Value

  end subroutine GetOptionValue_Int

  subroutine GetOptionValue_Int_Default(Option, Value, DefaultValue)

    type(t_cmd_opt), intent(in) :: Option
    integer, intent(out) :: Value
    integer, intent(in) :: DefaultValue

    if (DEBUG) then
      if (Option%value_type == CMD_OPT_NONE) then
        write (ERROR_UNIT, '(3a)') "ERROR: Attempted to read value of command line option '", &
          trim(Option%long_name), "' which has no value type specified."
        stop 1
      end if
    end if

    if (Option%is_present) then
      read (Option%value, *) Value
    else
      Value = DefaultValue
    end if

  end subroutine GetOptionValue_Int_Default

  subroutine GetOptionValue_Real(Option, Value)

    type(t_cmd_opt), intent(in) :: Option
    real(rk), intent(out) :: Value

    if (DEBUG) then
      if (.not. Option%is_present) then
        write (ERROR_UNIT, '(3a)') "ERROR: Attempted to read value of command line option '", &
          trim(Option%long_name), "' which is not present."
        stop 1
      end if
      if (Option%value_type == CMD_OPT_NONE) then
        write (ERROR_UNIT, '(3a)') "ERROR: Attempted to read value of command line option '", &
          trim(Option%long_name), "' which has no value type specified."
        stop 1
      end if
    end if

    read (Option%value, *) Value

  end subroutine GetOptionValue_Real

  subroutine GetOptionValue_Real_Default(Option, Value, DefaultValue)

    type(t_cmd_opt), intent(in) :: Option
    real(rk), intent(out) :: Value
    real(rk), intent(in) :: DefaultValue

    if (DEBUG) then
      if (Option%value_type == CMD_OPT_NONE) then
        write (ERROR_UNIT, '(3a)') "ERROR: Attempted to read value of command line option '", &
          trim(Option%long_name), "' which has no value type specified."
        stop 1
      end if
    end if

    if (Option%is_present) then
      read (Option%value, *) Value
    else
      Value = DefaultValue
    end if

  end subroutine GetOptionValue_Real_Default

  subroutine GetOptionValue_String(Option, Value)

    type(t_cmd_opt), intent(in) :: Option
    character(len=*), intent(out) :: Value

    if (DEBUG) then
      if (.not. Option%is_present) then
        write (ERROR_UNIT, '(3a)') "ERROR: Attempted to read value of command line option '", &
          trim(Option%long_name), "' which is not present."
        stop 1
      end if
      if (Option%value_type == CMD_OPT_NONE) then
        write (ERROR_UNIT, '(3a)') "ERROR: Attempted to read value of command line option '", &
          trim(Option%long_name), "' which has no value type specified."
        stop 1
      end if
    end if

    Value = trim(Option%value)

  end subroutine GetOptionValue_String

  subroutine GetOptionValue_String_Default(Option, Value, DefaultValue)

    type(t_cmd_opt), intent(in) :: Option
    character(len=*), intent(out):: Value
    character(len=*), intent(in) :: DefaultValue

    if (DEBUG) then
      if (Option%value_type == CMD_OPT_NONE) then
        write (ERROR_UNIT, '(3a)') "ERROR: Attempted to read value of command line option '", &
          trim(Option%long_name), "' which has no value type specified."
        stop 1
      end if
    end if

    if (Option%is_present) then
      Value = trim(Option%value)
    else
      Value = trim(DefaultValue)
    end if

  end subroutine GetOptionValue_String_Default

end module ovsCommandLine
