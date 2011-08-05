subroutine get_int_W(filename)
  use commonarrays
  use dyn_par
  use precision

  implicit none

  character (len=*), intent(inout)  :: filename  ! inout because we trim!

  integer           :: iline
  integer           :: num_oneints
  integer           :: i, j, ij
  complex(kind=pr)  :: eint

  !stop "succesfully tried to read CAP mo file"

  filename = trim(filename)
  open(11, file=filename, form='formatted')
  rewind(11)

  do iline = 1, 2     ! "CAP molecular ...", "Number of data lines..."
     read(11,*)
  end do

  read(11,*) num_oneints

  if (num_oneints > max1) then
     write(50,*)
     write(50,*) 'Problem in get_int_W.f'
     write(50,*) 'We are told there are ',num_oneints
     write(50,*) 'one-electron integrals in the Wints file'
     write(50,*) 'But max1 has been set at ',max1
     write(50,*)
     write(50,*) 'Increase max1.'
     write(50,*) 'Stopping...'
     stop
  end if

  do iline = 1, 2     ! "Format of data line:",  "i_mo, ..."
     read(11,*)
  end do

  if (ipoint(2) .ne. 1) then
    stop "get_int_W: ipoint array not initialized" 
  end if


  do iline = 1, num_oneints

     read(11,*) i,j,eint
     if (j>i) then
        write(50,*) 
        write(50,*) 'Problem in get_int_W.f'
        write(50,*) 'Found j>i while reading in the oneints.'
        write(50,*) 'i,j =',i,j
        write(50,*) 'Stopping...'
        stop
     end if
     ij = ipoint(i)+j
     e1ints(ij) = e1ints(ij) + eint

  end do

  close(11)

end subroutine get_int_W
