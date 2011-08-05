subroutine ldump(i_am_nu,i_am_mu,nu_doubly)
  use commonarrays
  use mcci_in

  implicit none

  integer, intent(in) :: i_am_nu, i_am_mu, nu_doubly

  integer             :: nn

  !     dump maximum coincidence lists to standard output

  do nn=1,ntotal

     if(   list(i_am_nu,nn).le.nbft.and.list(i_am_mu,nn).le.nbft)then

        write(*,'(i4,a,i4,a,i5,i5)')&
             &list(i_am_mu,nn),'u',list(i_am_nu,nn),'u'&
             &,my_pair(i_am_mu,nn),my_pair(i_am_nu,nn)

     elseif(   list(i_am_nu,nn).gt.nbft.and.list(i_am_mu,nn).le.nbft)then


        write(*,'(i4,a,i4,a,i5,i5)')&
             &list(i_am_mu,nn),'u',list(i_am_nu,nn)-nbft,'d'&
             &,my_pair(i_am_mu,nn),my_pair(i_am_nu,nn)


     elseif(list(i_am_nu,nn).le.nbft.and.list(i_am_mu,nn).gt.nbft)then


        write(*,'(i4,a,i4,a,i5,i5)')&
             &list(i_am_mu,nn)-nbft,'d',list(i_am_nu,nn),'u'&
             &,my_pair(i_am_mu,nn),my_pair(i_am_nu,nn)


     else

        write(*,'(i4,a,i4,a,i5,i5)')&
             &list(i_am_mu,nn)-nbft,'d',list(i_am_nu,nn)-nbft,'d'&
             &,my_pair(i_am_mu,nn),my_pair(i_am_nu,nn)

     endif

     if(nn.eq.nu_doubly) write(*,*)
  enddo

  return
end subroutine ldump
