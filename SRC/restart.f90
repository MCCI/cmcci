subroutine restart(length)
  use commonarrays, only: c, icij, nword
  use dyn_par
  use precision

  implicit none

  integer, intent(out) ::  length

  integer        :: idummy, i, n
  real(kind=pr)  :: temp

  open(unit=60,file='civ_in')
  i = 0
11 continue
  i = i + 1
  do n=1, nword
     !        read(60,*,end=22) idummy,icij(1,n,i),icij(2,n,i)
     if(n.eq.1) then
        !           read(60,*,end=22) c(i),icij(1,n,i),icij(2,n,i)
        read(60,*,end=22) idummy,temp,icij(1,n,i),icij(2,n,i)
        c(i) = cmplx(temp, 0.0, pr)
     else
        read(60,*,end=22) icij(1,n,i),icij(2,n,i)
     endif
  enddo
  goto 11
22 continue
  length = i-1

  return
end subroutine restart
