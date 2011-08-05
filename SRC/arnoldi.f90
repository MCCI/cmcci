subroutine arnoldi(h, s, length, ieig, idiag)
  use precision
  implicit none
  complex, dimension(:), intent(inout)  :: h
  real,    dimension(:)  intent(inout)  :: s
  integer, intent(in)                   :: length
  integer, intent(in)                   :: ieig
  integer, intent(in)                   :: idiag
  
  integer, parameter :: maxn   = 100000 !maximum dimension of matrix
  integer, parameter :: maxnev = 25     !max num of eigenvalues requested
  integer, parameter :: maxncv = 75     !max num of Lanczos basis vectors
  integer, parameter :: ldv    = maxn

  !Local arrays
  integer            :: iparam(11)
  integer            :: ipntr(14)
  integer            :: ipix(maxn)
  logical            :: select(maxncv)
  complex(kind=pr)   :: ax(maxn)
  complex(kind=pr)   :: d(maxncv)
  complex(kind=pr)   :: resid(maxn)
  complex(kind=pr)   :: v(ldv, maxncv)
  complex(kind=pr)   :: workd(3*maxn)
  complex(kind=pr)   :: workev(2*maxncv)
  complex(kind=pr)   :: workl(3*maxncv*maxncv+5*maxncv)
  complex(kind=pr)   :: dd(maxn), dl(maxn), du(maxn), du2(maxn)
  real(kind=pr)      :: rwork(maxn)    !specifically required for znaupd
  real(kind=pr)      :: rd(maxncv, 3)

  !Local scalars
  character          :: bmat*1, which*2
  integer            :: ido, n, nev, ncv
  integer            :: lworkl, info, ierr
  integer            :: j, nconv, maxitr, ishifts, mode
  complex(kind=pr)   :: sigma, hh
  real(kind=pr)      :: tol
  logical            :: rvec

  !Parameters
  complex(kind=pr), parameter   :: zero = (0.0_pr, 0.0_pr)
  complex(kind=pr), parameter   ::  one = (1.0_pr, 0.0_pr)

  !I dont think I need to declare the lapack rountines used,
  !but otherwise will declare them here.

  n   = length   !dimension of the matrix to diagonalise
  nev = ieig     !number of eigenvalue pairs we want computed
  ncv = 2*nev    !as indicated in ARPACK user's guide, p30.

  if ( n .gt. maxn )     stop "Arnoldi: ERROR: n is greater than maxn"
  if ( nev .gt. maxnev ) stop "Arnoldi: ERROR: nev is greater than maxnev"
  if ( ncv .gt. maxncv ) stop "Arnoldi: ERROR: ncv is greater than maxncv"

  bmat   = 'G'     !G for generalized eigenvalue problem
  which  = 'SM'    !Smallest Magnitude eigenvalues wanted, AUG p35
  sigma  = zero    !not sure we'll need this, as we won't be shift-inverting

  lworkl = 3*ncv**2 + 5*ncv  !Dimension of workl array
  tol    = 0.0_pr            !tol = 0.0 => tolerance = machine epsilon
  ido    = 0                 !reverse communication parameter; must be set
                             !to zero during initialization
  info   = 0                 !Request random starting vector
  
  ishifts = 1       !set regular (not shift-invert) mode
  maxitr  = 300     !max num of Arnoldi iterations  
  mode    = 2       !mode 2 is ??????????

  

end subroutine
