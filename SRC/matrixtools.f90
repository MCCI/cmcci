module matrixtools
  use precision

  ! Subroutines mxv_mreal and mxv_mcompl are nearly identical, the only
  ! difference being the type of the sparse matrix factor; real in
  ! one case, complex in the other. The two are overloaded (see the
  ! interface) so other parts of the code can call mxv_sparse with
  ! real (s) or complex (h) arguments.
  ! If, in the future, we want to make an all real mcci again, a
  ! third version will be needed, in which all arguments are reals.
  ! MS

  ! The module also contains the unpack_matrix subroutine, which converts
  ! a packed matrix into conventional form.

  interface mxv_sparse
    module procedure mxv_mreal, mxv_mcompl
  end interface

  contains

  subroutine mxv_mreal(length,m,ijm,v,mv)
  ! real matrix (in our case, S) but complex vectors
  implicit none

  integer,             intent(in)   :: length
  real(kind=pr),       intent(in)   :: m(*)
  integer,             intent(in)   :: ijm(*)
  complex(kind=pr),    intent(in)   :: v(*)
  complex(kind=pr),    intent(out)  :: mv(*)

  integer  :: ici, jci, k
  if( ijm(1) .ne. length + 2) STOP 'mxv_mreal: length mismatch'  

  do ici=1, length
     !diagonal
     mv(ici) = m(ici)*v(ici)
  enddo

  do ici=1, length
     do k=ijm(ici), ijm(ici+1)-1
        !lower off diagonals
        mv(ici) = mv(ici) + m(k)*v(ijm(k))
        !upper off diagonals
        jci =ijm(k)
        mv(jci) = mv(jci) + m(k)*v(ici)
     enddo
  enddo

  return
 
  end subroutine

  subroutine mxv_mcompl(length,m,ijm,v,mv)
  ! both complex matrix and complex vectors
  implicit none

  integer,             intent(in)   :: length
  complex(kind=pr),    intent(in)   :: m(*)
  integer,             intent(in)   :: ijm(*)
  complex(kind=pr),    intent(in)   :: v(*)
  complex(kind=pr),    intent(out)  :: mv(*)

  integer  :: ici, jci, k
  if( ijm(1) .ne. length + 2) STOP 'mxv_mcompl: length mismatch'  

  do ici=1, length
     !diagonal
     mv(ici) = m(ici)*v(ici)
  enddo

  do ici=1, length
     do k=ijm(ici), ijm(ici+1)-1
        !lower off diagonals
        mv(ici) = mv(ici) + m(k)*v(ijm(k))
        !upper off diagonals
        jci =ijm(k)
        mv(jci) = mv(jci) + m(k)*v(ici)
     enddo
  enddo

  return
  end subroutine

  subroutine unpack_matrix(kk, packed_matrix, unpacked_matrix)
    ! Converts a matrix in upper-packed form (UPLO = U in LAPACK subroutines
    ! such as dspgv) to the more conventional, two-dimensional form.
    implicit none

    integer,          intent(in)   :: kk
    complex(kind=pr), intent(in)   :: packed_matrix(:)
    complex(kind=pr), intent(out)  :: unpacked_matrix(:,:)

    integer  :: i, j

    do j = 1, kk
      do i = 1, j
        unpacked_matrix(i, j) = packed_matrix(i + (j-1)*j/2)
        unpacked_matrix(j, i) = unpacked_matrix(i, j)
      end do
    end do
  end subroutine

end module matrixtools
