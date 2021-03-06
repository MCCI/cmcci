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
  dnorm = (0.0_pr, 0.0_pr)
  do ici=1, length
     eval  = eval  + conjg(c(ici))*h(ici)*c(ici)
     dnorm = dnorm + conjg(c(ici))*s(ici)*c(ici)
     do jci = ijh(ici), ijh(ici+1) - 1
        !eval  =  eval  + 2.0*real(conjg(c(ici))*h(jci)*c(ijh(jci)), kind=pr)
        eval  =  eval  + conjg(c(ici))*h(jci)*c(ijh(jci))  &
                       + c(ici)*h(jci)*conjg(c(ijh(jci)))
     enddo
     do jci = ijs(ici), ijs(ici+1) - 1
        dnorm =  dnorm + conjg(c(ici))*s(jci)*c(ijs(jci))  &
                       + c(ici)*s(jci)*conjg(c(ijs(jci)))
     enddo
  enddo

  eval = eval/dnorm

  return
end subroutine energy
