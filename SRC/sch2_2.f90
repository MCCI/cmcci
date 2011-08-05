subroutine sch2_2(element,i_am_mu,i_am_nu,idiff1,idiff2,kck,n_2p,ep)
  use commonarrays, only: list, nbft, e2ints, ipoint, i_sx2
  use dyn_par
  use precision
  implicit none

  complex(kind=pr),  intent(out)     :: element
  integer,           intent(in)      :: i_am_mu, i_am_nu
  integer,           intent(in)      :: idiff1, idiff2
  integer,           intent(in)      :: kck, n_2p
  integer,           intent(in)      :: ep

  integer        :: iorb, jorb, korb, lorb
  integer        :: ispin, jspin, kspin, lspin
  integer        :: ik, jl, ikjl, il, jk, iljk
  real(kind=pr)  :: ck, sc1, sc2

  iorb    = list(i_am_mu,idiff1)
  jorb    = list(i_am_mu,idiff2)
  korb    = list(i_am_nu,idiff1)
  lorb    = list(i_am_nu,idiff2)

  ispin = 1
  jspin = 1
  kspin = 1
  lspin = 1
  if(iorb.gt.nbft) then
     iorb = iorb - nbft
     ispin = 0
  endif
  if(jorb.gt.nbft) then
     jorb  = jorb - nbft
     jspin = 0
  endif
  if(korb.gt.nbft) then
     korb  = korb - nbft
     kspin = 0
  endif
  if(lorb.gt.nbft) then
     lorb  = lorb - nbft
     lspin = 0
  endif

  if(iorb.ge.korb) then
     ik = ipoint(iorb) + korb
  else
     ik = ipoint(korb) + iorb
  endif

  if(jorb.ge.lorb) then
     jl = ipoint(jorb) + lorb
  else
     jl = ipoint(lorb) + jorb
  endif

  if(ik.ge.jl) then
     ikjl = ipoint(ik) + jl
  else
     ikjl = ipoint(jl) + ik
  endif

  sc1 = ck(i_sx2,n_2p,kck)

  element = cmplx(sc1*e2ints(ikjl), 0.0, pr)

  if(iorb.ge.lorb) then
     il = ipoint(iorb) + lorb
  else
     il = ipoint(lorb) + iorb
  endif

  if(jorb.ge.korb) then
     jk = ipoint(jorb) + korb
  else
     jk = ipoint(korb) + jorb
  endif

  if(il.ge.jk) then
     iljk = ipoint(il) + jk
  else
     iljk = ipoint(jk) + il
  endif

  if( (ispin.eq.jspin) .or. (kspin.eq.lspin) ) then
     sc2 = sc1
  elseif( (ispin.eq.kspin) .and. (jspin.eq.lspin) ) then
     sc2 = ck(i_sx2,n_2p,kck+1)
  elseif( (ispin.eq.lspin) .and. (jspin.eq.kspin) ) then
     sc2 = ck(i_sx2,n_2p,kck-1)
  endif

  element = element - cmplx(sc2*e2ints(iljk), 0.0, pr)

  element = ep*element

  return
end subroutine sch2_2
