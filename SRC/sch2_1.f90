subroutine sch2_1(element,i_am_mu,i_am_nu,idiff1,idiff2,kck,n_2p,ep)
  use commonarrays, only: list, nbft, e2ints, ipoint
  use mcci_in
  use precision
  implicit none

  complex(kind=pr),  intent(out)     :: element
  integer,           intent(in)      :: i_am_mu, i_am_nu
  integer,           intent(in)      :: idiff1, idiff2
  integer,           intent(in)      :: kck, n_2p
  integer,           intent(in)      :: ep

  integer        :: iorb, jorb, korb, lorb
  integer        :: ik, jl, ikjl, il, jk, iljk
  real(kind=pr)  :: ck, sc


  iorb    = list(i_am_mu,idiff1)
  jorb    = list(i_am_mu,idiff2)
  korb    = list(i_am_nu,idiff1)
  lorb    = list(i_am_nu,idiff2)

  if(iorb.gt.nbft) iorb = iorb - nbft
  if(jorb.gt.nbft) jorb = jorb - nbft
  if(korb.gt.nbft) korb = korb - nbft
  if(lorb.gt.nbft) lorb = lorb - nbft

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

  element = cmplx(e2ints(ikjl), 0.0, pr)

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

  element = element - cmplx(e2ints(iljk), 0.0, pr)

  sc = ck(i_sx2,n_2p,kck)

  element = ep*sc*element

  return
end subroutine sch2_1
