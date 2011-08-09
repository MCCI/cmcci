subroutine sch1_2(element,i_am_mu,i_am_nu,nu_doubly,idiff1,kck,n_2p,ep)
  use commonarrays, only: nbft, e1ints, e2ints, ipoint, &
                          list, my_pair
  use mcci_in
  implicit none


  complex(kind=pr), intent(out)    :: element
  integer,          intent(in)     :: i_am_mu, i_am_nu, nu_doubly
  integer,          intent(in)     :: idiff1
  integer,          intent(in)     :: kck, n_2p
  integer,          intent(in)     :: ep

  integer         :: n, lorb, iorb, jorb 
  real(kind=pr)   :: sc1, sc2
  real(kind=pr)   :: ck
  integer         :: ispin, jspin, kspin, lspin
  integer         :: jl, ii, iijl, il, ji, ilji
  integer         :: idiff1_spin 


  !     one orbital difference Slater Condon Harris rules

  jorb = list(i_am_mu,idiff1)
  lorb = list(i_am_nu,idiff1)

  if(jorb.gt.nbft) jorb  = jorb - nbft
  if(lorb.gt.nbft) lorb  = lorb - nbft

  if(jorb.ge.lorb) then
     jl = ipoint(jorb) + lorb
  else
     jl = ipoint(lorb) + jorb
  endif

  element = e1ints(jl)

  do n = 2, nu_doubly, 2

     iorb = list(i_am_nu,n)
     if(iorb.gt.nbft) iorb = iorb-nbft

     ii = ipoint(iorb) + iorb

     if(ii.ge.jl) then
        iijl = ipoint(ii) + jl
     else 
        iijl = ipoint(jl) + ii
     endif

     element = element + 2.0*e2ints(iijl)          ! "Coulomb terms"

     if(iorb.ge.lorb) then
        il = ipoint(iorb) + lorb
     else
        il = ipoint(lorb) + iorb
     endif

     if(jorb.ge.iorb) then
        ji = ipoint(jorb) + iorb
     else
        ji = ipoint(iorb) + jorb
     endif

     if(il.ge.ji) then
        ilji = ipoint(il) + ji
     else
        ilji = ipoint(ji) + il
     endif

     element = element - e2ints(ilji)        ! "Exchange terms"

  enddo

  sc1 = ck(i_sx2,n_2p,kck)
  element  = sc1*element

  ispin = 1
  if(list(i_am_mu,my_pair(i_am_mu,idiff1)).gt.nbft) ispin = 0

  kspin = 1
  if(list(i_am_nu,my_pair(i_am_mu,idiff1)).gt.nbft) kspin = 0

  idiff1_spin = 1
  if(list(i_am_mu,idiff1).gt.nbft) idiff1_spin = 0

  do n=nu_doubly+1, ntotal

     iorb = list(i_am_nu,n)
     if(iorb.gt.nbft) iorb = iorb-nbft

     jspin = 1
     if(list(i_am_mu,n).gt.nbft) jspin = 0

     lspin = 1
     if(list(i_am_nu,n).gt.nbft) lspin = 0

     ii = ipoint(iorb) + iorb

     if(ii.ge.jl) then
        iijl = ipoint(ii) + jl
     else 
        iijl = ipoint(jl) + ii
     endif

     element = element + sc1*e2ints(iijl)    ! single "Coulomb terms"

     if(jspin.eq.idiff1_spin) then

        if( iorb .ge. lorb) then
           il = ipoint(iorb) + lorb
        else
           il = ipoint(lorb) + iorb
        endif

        if( jorb .ge. iorb) then
           ji = ipoint(jorb) + iorb
        else
           ji = ipoint(iorb) + jorb
        endif

        if(il.ge.ji) then
           ilji = ipoint(il) + ji
        else 
           ilji = ipoint(ji) + il
        endif

        if(kspin.eq.lspin) then
           sc2= 0.0                       ! coefs cancel  
        elseif(ispin.eq.kspin.and.jspin.eq.lspin) then
           sc2= ck(i_sx2,n_2p,kck+1) - sc1
        elseif(ispin.eq.lspin.and.jspin.eq.kspin) then
           sc2= ck(i_sx2,n_2p,kck-1) - sc1
        endif

        element = element + sc2*e2ints(ilji)

     endif
  enddo

  element =  ep*element

  return
end subroutine sch1_2
