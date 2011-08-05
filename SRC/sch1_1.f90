subroutine sch1_1(element,i_am_mu,i_am_nu,nu_doubly,idiff1,kck,n_2p,ep)
  use commonarrays, only: nbft, e1ints, e2ints, ntotal, i_sx2, ipoint, list
  use mcci_in
  implicit none


  complex(kind=pr), intent(out)    :: element
  integer,          intent(in)     :: i_am_mu, i_am_nu, nu_doubly
  integer,          intent(in)     :: idiff1
  integer,          intent(in)     :: kck, n_2p
  integer,          intent(in)     :: ep

  integer         :: n, lorb, iorb, jorb 
  real(kind=pr)   :: sc
  real(kind=pr)   :: ck
  integer         :: ispin, jspin, kspin, lspin
  integer         :: jl, ii, iijl, il, ji, ilji


  !     one orbital difference Slater Condon Harris rules

  jorb = list(i_am_mu,idiff1)
  lorb = list(i_am_nu,idiff1)

  jspin = 1
  if(jorb.gt.nbft) then
     jorb  = jorb - nbft
     jspin = 0
  endif

  if(lorb.gt.nbft) lorb  = lorb - nbft

  if(jorb .ge. lorb) then
     jl = ipoint(jorb) + lorb
  else
     jl = ipoint(lorb) + jorb
  endif

  element = e1ints(jl)

  do n = 2, nu_doubly, 2

     iorb = list(i_am_nu,n)
     if(iorb.gt.nbft) iorb  = iorb-nbft

     ii = ipoint(iorb) + iorb

     if(ii.ge.jl) then
        iijl = ipoint(ii) + jl
     else 
        iijl = ipoint(jl) + ii
     endif

     element = element + 2.0*e2ints(iijl)          ! "Coulomb terms"

     if(iorb .ge. lorb) then
        il = ipoint(iorb) + lorb
     else
        il = ipoint(lorb) + iorb
     endif

     if(jorb .ge. iorb) then
        ji = ipoint(jorb) + iorb
     else
        ji = ipoint(iorb) + jorb
     endif

     if(il .ge. ji ) then
        ilji = ipoint(il) + ji
     else
        ilji = ipoint(ji) + il
     endif

     element = element - e2ints(ilji)        ! "Exchange terms"

  enddo

  do n=nu_doubly+1, ntotal

     iorb  = list(i_am_mu,n)
     ispin = 1
     if(iorb.gt.nbft) then
        iorb  = iorb-nbft
        ispin = 0
     endif

     ii = ipoint(iorb) + iorb

     if(ii.ge.jl) then
        iijl = ipoint(ii) + jl
     else 
        iijl = ipoint(jl) + ii
     endif

     element = element + e2ints(iijl)     ! single "Coulomb terms"

     !        sum is spin restricted
     if(ispin.eq.jspin) then

        if(iorb .ge. lorb) then
           il = ipoint(iorb) + lorb
        else
           il = ipoint(lorb) + iorb
        endif

        if(jorb .ge. iorb) then
           ji = ipoint(jorb) + iorb
        else
           ji = ipoint(iorb) + jorb
        endif

        if(il.ge.ji) then
           ilji = ipoint(il) + ji
        else 
           ilji = ipoint(ji) + il
        endif

        element = element - e2ints(ilji)

     endif

  enddo

  sc = ck(i_sx2,n_2p,kck)
  element =  ep*sc*element

  return
end subroutine sch1_1
