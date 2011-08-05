module zggev_wrapper
use precision

contains

subroutine diagonalize(N, A, B, EVALS, EVECS)
  implicit none

  integer,          intent(in)  :: N
  complex(kind=pr), intent(in)  :: A(N,N)
  complex(kind=pr), intent(in)  :: B(N,N)
  complex(kind=pr), intent(out) :: EVALS(N)
  complex(kind=pr), intent(out) :: EVECS(N,N)

  complex(kind=pr)               :: alpha(N)
  complex(kind=pr)               :: beta(N)
  complex(kind=pr)               :: vr(N,N)
  complex(kind=pr)               :: vl(N,N)
  complex(kind=pr), allocatable  :: work(:)
  integer                        :: lwork
  real(kind=pr)                  :: rwork(8*N)
  integer                        :: info

  integer           :: k, j
  complex(kind=pr)  :: temp
  complex(kind=pr)  :: tempvec(N)
  logical           :: switch

  lwork = -1
  allocate(work(1))
  ! initial call to fetch optimal value of lwork
  call zggev ('V', 'V', N, A, N, B, N, alpha, beta, vl, N, vr, N, &
               work, lwork, rwork, info)
  deallocate(work)
  
  lwork = int(work(1))
  allocate(work(lwork))
  call zggev ('V', 'V', N, A, N, B, N, alpha, beta, vl, N, vr, N, &
               work, lwork, rwork, info)

  if(info.ne.0) then
    print*, 'info = ', info
    call dumparrays()
    STOP 'Problems with zggev '
  end if

  ! zggev returns the generalized eigenvalues as ratios of alpha(i) and
  ! beta(i). We need to check the conditioning and put the reults in 
  ! the e array.
  do k = 1, N
    if(abs(beta(k)) .lt.1.0d-4) then
      call dumparrays()
      STOP 'davidson: Dangerously low beta(k)'
    else
      EVALS(k) = alpha(k) / beta(k)
    end if
  end do

  ! Now we need to sort the eigenvalues, to make sure we grab the ieig
  ! most significant ones.
  do k = 2, N
    do j = k, 2, -1
      switch = real(EVALS(j)) < real(EVALS(j-1))      &
      .or. (real(EVALS(j)) == real(EVALS(j-1)) &
      .and. aimag(EVALS(j)) < aimag(EVALS(j-1)))
      if(switch) then
        temp       = EVALS(j)
        EVALS(j)   = EVALS(j-1)
        EVALS(j-1) = temp
        tempvec(1:N)  = vr(:,j)
        vr(:,j)       = vr(:,j-1)
        vr(:,j-1)     = tempvec(1:N)
        tempvec(1:N)  = vl(:,j)
        vl(:,j)       = vl(:,j-1)
        vl(:,j-1)     = tempvec(1:N)
        !print*, 'davidson: some swapping occurred'
      end if
    end do
  end do

  deallocate(work)

  contains
  subroutine dumparrays()
 
    print*, 'zggev_wrapper: dumping all variables and arrays'
    print*, 'N: ',        N
    print*, 'A: ',        A
    print*, 'EVALS: ',    EVALS
    print*, 'EVECS: ',    EVECS
    print*, 'alpha: ',    alpha
    print*, 'beta: ',     beta
    print*, 'vr: ',       vr
    print*, 'vl: ',       vl
    print*, 'work: ',     work
    print*, 'lwork: ',    lwork
    print*, 'rwork: ',    rwork
    print*, 'info: ',     info
    print*, 'k: ',  k, 'j: ', j
    print*, 'temp: ',     temp
    print*, 'tempvec: ',  tempvec
    print*, 'switch: ',   switch
   
  end subroutine dumparrays

  end subroutine diagonalize

end module
