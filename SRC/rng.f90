module rng
use precision
use commonarrays, only: me

implicit none

logical :: fixed_seed = .false.

contains

subroutine rng_init_clock()
  implicit none

  integer :: i
  ! i will be set to 0, but the system clock will be started and thus will
  ! return a nonzero value when next called (in rng_gen_seed).
  ! Accordingly, this call would ideally be placed as far as possible ahead of
  ! the call to rng_gen_seed(), to increase the variance in the returned values
  ! of nseed.
  ! This is only strictly necessary for some compilers (notably pgi); others
  ! will just return a clock count offset from an arbitrary time in the past,
  ! each time system_clock is called, and for those this call is
  ! inconsequential.
  call system_clock( count=i )

end subroutine


subroutine rng_gen_seed(seed)
  implicit none

  ! Must be called before generating any random numbers, and after a call to
  ! rng_init_clock (see above).

  real(kind=pr), intent(out) :: seed
  integer                    :: nseed, nfrac

  if ( fixed_seed ) then
    nseed= 18800925
    if(me.eq.0) write(50,*) 'Warning: seed set by hand'
  else
    call system_clock( count=nseed )
    print*, 'nseed : ', nseed
  end if

  nseed = nseed / ( me + 1)
  nfrac = nseed/199017
  nseed = nseed - nfrac*199017
  seed  = float(nseed)

  return
end subroutine rng_gen_seed

subroutine rng_random(seed,rand)
  implicit none

  real(kind=pr), intent(inout)  ::  seed
  real(kind=pr), intent(out)    ::  rand

  real(kind=pr), parameter  :: m = 199017.0
  real(kind=pr), parameter  :: a = 24298.0
  real(kind=pr), parameter  :: c = 1.0

  real(kind=pr)  :: k

  seed=a*seed+c
  k    = dint(seed/m)
  seed = seed-m*k
  rand = seed/(m+1.0)
  return
end subroutine rng_random


end module
