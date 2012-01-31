subroutine restart(length)
  use commonarrays, only: c, nword
  use mcci_in
  use precision

  implicit none

  integer, intent(out) ::  length

  character (len=120)  :: line
  character (len=15)   :: numbers = '0123456789+-.Ee'
  character            :: char_tmp
  integer              :: stat = 0
  integer              :: number_of_numbers
  integer              :: iscan
  logical              :: in_number
  logical              :: wavefunction_is_real
  integer              :: lines_in_file, i, n
  integer              :: csf_number
  real (kind=pr)       :: c_real
  real (kind=pr)       :: c_imag

  open(unit=60,file='civ_in', form='formatted', status='old')
  rewind(60)
  lines_in_file = 0

  ! Check the integrity of the file
  do while (stat >= 0)
    read(60, *, IOStat=stat)
    if (stat == 0) then
      lines_in_file = lines_in_file + 1
    else if (stat > 0) then
      write(*,*) "Error reading civ_in"
      write(*,*) "iostat", stat
      write(*,*) "Stopping..."
      stop
    endif
  enddo

  if (mod(lines_in_file, nword) /= 0) then
    write(*,*) "Error reading civ_in"
    write(*,*) "Number of lines is not a multiple of nword:"
    write(*,*) "nlines: ", lines_in_file
    write(*,*) "nword: ", nword
    write(*,*) "Stopping..."
    stop
  endif

  length = lines_in_file / nword
  
  ! Find out if the wave function is real or complex (it is stored as
  ! complex(kind=pr) either way)

  rewind(60)
  read(60, '(a120)') line

  in_number = .false.
  number_of_numbers = 0

  do iscan = 1, len(line)
    char_tmp = line(iscan:iscan)
    if ((scan(char_tmp, numbers) == 0) .and. in_number) then
      in_number = .false.
      cycle
    endif

    if ((scan(char_tmp, numbers) > 0) .and. .not. in_number) then
      in_number = .true.
      number_of_numbers = number_of_numbers + 1
    end if
  end do

  select case (number_of_numbers)
    case (4)
      wavefunction_is_real = .true.
    case (5)
      wavefunction_is_real = .false.
    case default
      write (*,*) "Error reading civ_in"
      write (*,*) "Found "
      write (*,*) number_of_numbers
      write (*,*) "numbers in the first line of civ_in."
      write (*,*) "This must be 4 or 5, for the real and complex case, resp."
      write (*,*) "Stopping..."
      stop
  end select

  ! Now read in the actual contents of the file
  icij = 0
  c    = (0.0_pr, 0.0_pr)

  rewind(60)
  i = 0
  stat = 0

  do while (stat >= 0)
    i = i + 1
    do n = 1, nword
      if (n == 1) then
        if (wavefunction_is_real) then
          read(60, *, IOStat=stat) csf_number, c_real, icij(1,n,i), icij(2,n,i)
          c(i) = c_real
        else
          read(60, *, IOStat=stat) csf_number, c_real, c_imag, icij(1,n,i), icij(2,n,i)
          c(i) = cmplx(c_real, c_imag)
        endif
      else
        read(60, *, IOStat=stat) icij(1,n,i), icij(2,n,i)
      endif
    enddo
  enddo

  close(60)

end subroutine restart
