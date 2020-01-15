
!-------------------
! Count number of lines in a file

subroutine row_count(in_unit, N)


implicit none

integer     :: in_unit, N, io
integer :: id

  rewind(in_unit)
  N = 0
  do
    read(in_unit,*,iostat=io) id
    if (io /= 0) exit
    N = N+1
  enddo

end subroutine row_count


subroutine colcount(in_unit, nargs)

implicit none

integer :: in_unit, nargs
integer :: n_non_space, n_char, n_word, io
logical :: prev_space
character(len=1) :: c

  n_non_space = 0
  n_char = 0
  n_word = 0
  prev_space = .TRUE. ! if 1st character is space, word count is increased
  do
    read(unit=in_unit, FMT='(A)', ADVANCE='NO', iostat=io) c
    if (io /= 0) exit
    n_char = n_char + 1
    if ((c == ' ').OR.(ichar(c) < 32)) then ! space or invisible character?
    !if (c == ' ') then
      prev_space = .TRUE.
    else
      if (prev_space) n_word = n_word + 1
      n_non_space = n_non_space + 1
      prev_space = .FALSE.
    endif
  enddo

  !print *,'Number of characters on 1st line=',n_char
  !print *,'Number of non-spaces on 1st line=',n_non_space
  !print *,'Number of words      on 1st line=',n_word

  nargs = n_word

end subroutine colcount
!-------------------
