subroutine sch2(element,i_am_mu,i_am_nu,nu_doubly,idiff1,idiff2,kck,n_2p,ep)
  use mcci_in
  implicit none

  complex(kind=pr),  intent(out)     :: element
  integer,           intent(in)      :: i_am_mu, i_am_nu, nu_doubly
  integer,           intent(in)      :: idiff1, idiff2
  integer,           intent(in)      :: kck, n_2p
  integer,           intent(in)      :: ep

  integer                            :: my_case 


  !     determine what doubles Slater-Condon-Harris rules apply
  call get_case2(idiff1,idiff2,nu_doubly,i_am_mu,i_am_nu,&
       ep,kck,my_case)

  if(my_case.eq.0) then
     call sch2_0(element,i_am_mu,i_am_nu,idiff1,idiff2,kck,n_2p,ep)
  elseif(my_case.eq.1) then
     call sch2_1(element,i_am_mu,i_am_nu,idiff1,idiff2,kck,n_2p,ep)
  elseif(my_case.eq.2) then
     call sch2_2(element,i_am_mu,i_am_nu,idiff1,idiff2,kck,n_2p,ep)
  elseif(my_case.eq.3) then
     call sch2_3(element,i_am_mu,i_am_nu,&
          idiff1,idiff2,nu_doubly,kck,n_2p,ep)
  elseif(my_case.eq.4) then
     call sch2_4(element,i_am_mu,i_am_nu,idiff1,idiff2,kck,n_2p,ep)
  elseif(my_case.eq.5) then
     call sch2_5(element,i_am_mu,i_am_nu,idiff1,idiff2,kck,n_2p,ep)
  else
     write(*,*)'my_case',my_case
     call ldump(i_am_nu,i_am_mu,nu_doubly)
     STOP 'sch2: no conditions were met'
  endif

  return 
end subroutine sch2
