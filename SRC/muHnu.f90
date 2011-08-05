subroutine muHnu(ici,jci,element,n_2p,s_overlap)
  use commonarrays, only: nbft, icij, nword, i_sx2
  use mcci_in
  use precision
  implicit none

  integer,           intent(in)     :: ici
  integer,           intent(in)     :: jci
  complex(kind=pr),  intent(out)    :: element
  integer,           intent(out)    :: n_2p
  real(kind=pr),     intent(out)    :: s_overlap

  integer        :: isingles, jsingles, idoubles, jdoubles
  integer        :: ijsingles(iword), ijdoubles(iword),ij_sd(iword)
  integer        :: n, ndiff, j, jshift, itest
  logical        :: btest
  integer        :: i_am_mu, i_am_nu, nu_doubly, idiff1, idiff2
  integer        :: kck
  integer        :: ep
  real(kind=pr)  :: ck

  element = cmplx(0.0, 0.0, pr)
  s_overlap = 0.0_pr

  !     find number of differences and if matrix element = 0, return
  do n=1,nword 

     isingles = ieor(icij(1,n,ici),icij(2,n,ici))
     jsingles = ieor(icij(1,n,jci),icij(2,n,jci))

     idoubles = iand(icij(1,n,ici),icij(2,n,ici))
     jdoubles = iand(icij(1,n,jci),icij(2,n,jci))

     !        all singles
     ijsingles(n) = ieor(isingles,jsingles)
     !        doubles only in one configuration
     ijdoubles(n) = ieor(idoubles,jdoubles)
     ij_sd(n)     = iand(ijsingles(n),ijdoubles(n))

  enddo

  !     sum number of difference orbitals
  ndiff = 0
  do j = 0, nbft-1
     n = j/int_bits+1
     jshift = j - (n-1)*int_bits
     itest = 0
     if(btest(ijdoubles(n),jshift)) itest = 1
     ndiff = ndiff+2*itest
     itest = 0
     if(btest(ijsingles(n),jshift)) itest = 1
     ndiff = ndiff+  itest
     itest = 0
     if(btest(ij_sd(n),jshift)) itest = 1
     ndiff = ndiff-2*itest
  enddo
  if(ndiff.gt.4) return                     ! matrix element = 0

  !     note ndiff can change due to restrictions on the orderings,
  !     in any event, 4 --> 2, 2 --> 1 , 0 --> 0
  call reorder(ici,jci,i_am_mu,i_am_nu,nu_doubly,ndiff,idiff1,idiff2,  &
               kck,n_2p,ep)

  if(ndiff.gt.2) return

  !     Slater Condon Harris rules
  if( ndiff .eq. 0) then
     s_overlap = ck(i_sx2,n_2p,kck)
     ! no orbital differences
     call sch0(element,i_am_mu,i_am_nu,nu_doubly,kck,n_2p,ep)
  elseif( ndiff .eq. 1) then
     ! one orbital difference
     call sch1(element,i_am_mu,i_am_nu,nu_doubly,idiff1,kck,n_2p,ep)
  elseif( ndiff .eq. 2) then
     ! two orbital differences
     call sch2(element,i_am_mu,i_am_nu,nu_doubly,idiff1,idiff2,kck,n_2p,ep)
  else
     write(*,*) 'ndiff', ndiff
     STOP 'muHnu: no conditions were met'
  endif

      return
      end
