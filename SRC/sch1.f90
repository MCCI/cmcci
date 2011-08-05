subroutine sch1(element,i_am_mu,i_am_nu,nu_doubly,idiff1,kck,n_2p,ep) 
  !use commonarrays, only: list, nsym, nbft, nbpsy !none of these are used ?
  use dyn_par
  implicit none

  complex(kind=pr),  intent(out)     :: element
  integer,           intent(in)      :: i_am_mu, i_am_nu, nu_doubly
  integer,           intent(in)      :: idiff1
  integer,           intent(in)      :: kck, n_2p
  integer,           intent(in)      :: ep

  integer                            :: my_case 

!     determine what singles Slater-Condon-Harris rules apply
  call get_case1(idiff1,i_am_mu,i_am_nu,my_case)

  if(my_case.eq.0) then
   call sch1_0(element,i_am_mu,i_am_nu,nu_doubly,idiff1,kck,n_2p,ep)
  elseif(my_case.eq.1) then
   call sch1_1(element,i_am_mu,i_am_nu,nu_doubly,idiff1,kck,n_2p,ep)
  elseif(my_case.eq.2) then
   call sch1_2(element,i_am_mu,i_am_nu,nu_doubly,idiff1,kck,n_2p,ep)
  else
   write(*,*)'my_case',my_case
   call ldump(i_am_nu,i_am_mu,nu_doubly)
   STOP 'sch1: no conditions were met'
  endif

return
end
