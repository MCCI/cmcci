subroutine exc(l_end,l_begin,l_to,iexc)
  use commonarrays, only: c, icij, me
  use mcci_in
  use precision

  implicit none
  integer, intent(in)     :: l_end
  integer, intent(in)     :: l_begin
  integer, intent(inout)  :: l_to
  integer, intent(in)     :: iexc

  integer             :: ntype, ibyte, ntime, lbuf
  integer             :: lpts, lptr
  integer             :: ls, lr
  integer, parameter  :: MAXBUF = 1048576 ! must be multiple of 8
  integer             :: ierr
  integer             :: idiff
  integer             :: k, ici
  include          'mpif.h' 

  !     l_end       end of buffer to broadcast from node iexc
  !     l_begin+1   beginning of buffer to broadcast from node iexc
  !     l_to        end of buffer on receiving nodes, upon return l_to=l_to+lr

  ntype = 1

  ls = l_end - l_begin      ! buffer to send   / 8
  if(ls.lt.0) ls = 0

  ntype = ntype + 1
  if(me.eq.iexc) then
     !        call brdcst(ntype,ls,4,iexc)                          ! TCGMSG
     call MPI_BCAST(ls,4,MPI_BYTE,iexc,MPI_COMM_WORLD,ierr)! MPI
     ibyte = 2*nbyte*iword*ls
  else
     !        call brdcst(ntype,lr,4,iexc)		               ! TCGMSG
     call MPI_BCAST(lr,4,MPI_BYTE,iexc,MPI_COMM_WORLD,ierr)! MPI
     ibyte = 2*nbyte*iword*lr
  endif

  ntime = ibyte/MAXBUF + 1

  idiff = ibyte - (ntime-1)*MAXBUF

  lpts = l_begin
  lptr = l_to

  do k = 1, ntime
     lbuf = MAXBUF
     if(k .eq. ntime) lbuf = idiff
     ntype = ntype + 1
     if(me.eq.iexc) then
        !           call brdcst(ntype,icij(1,1,lpts+1),lbuf,iexc)         ! TCGMSG
        call MPI_BCAST(icij(1,1,lpts+1),lbuf,MPI_BYTE,iexc,MPI_COMM_WORLD,ierr)     ! MPI
     else
        !           call brdcst(ntype,icij(1,1,lptr+1),lbuf,iexc)         ! TCGMSG
        call MPI_BCAST(icij(1,1,lptr+1),lbuf,MPI_BYTE,iexc,MPI_COMM_WORLD,ierr)   ! MPI
     endif
     lpts = lpts + MAXBUF/2/nbyte/iword
     lptr = lptr + MAXBUF/2/nbyte/iword
  enddo

  if(me.ne.iexc) then 
     do ici = l_to+1, l_to+lr+1
        c(ici) = 0.0
     enddo
     l_to = l_to + lr
  endif

  if(l_to.gt.maxc-1) STOP 'error exc: exceeded max. array length'

  return
end subroutine exc
