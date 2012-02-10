subroutine energy(length,eval,dnorm)
  use commonarrays, only: c, h, ijh, s, ijs
  use mcci_in
  use precision

  implicit none

  integer,           intent(in)    :: length
  complex(kind=pr),  intent(out)   :: eval
  complex(kind=pr),     intent(out)   :: dnorm

  integer   :: ici, jci
  if(ijh(1) .ne. length + 2) STOP 'energy: h array mismatch'
  if(ijs(1) .ne. length + 2) STOP 'energy: s array mismatch'

  eval  = (0.0_pr, 0.0_pr)
  dnorm = 0.0_pr
  do ici=1, length
     eval  = eval  + c(ici)*h(ici)*c(ici)
     dnorm = dnorm + c(ici)*s(ici)*c(ici)
     do jci = ijh(ici), ijh(ici+1) - 1
        eval  =  eval  + 2.0*c(ici)*h(jci)*c(ijh(jci))
     enddo
     do jci = ijs(ici), ijs(ici+1) - 1
        dnorm =  dnorm + 2.0*c(ici)*s(jci)*c(ijs(jci))
     enddo
  enddo

  eval = eval/dnorm

  return
end subroutine energy
