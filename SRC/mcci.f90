program mcci
  use commonarrays
  use mcci_in
  use precision
  use matrixtools
  use rng
  !use mpi

  implicit none
  include          'mpif.h'
  character(len=20)     :: fdate
  real                  :: twb,twg1,twp0,twp,twi,twf !timing variables
  real                  :: tub,tug1,tup0,tup,tui,tuf
  real                  :: tsb,tsg1,tsp0,tsp,tsi,tsf
  real(kind=pr)                 :: dwsum, dusum, dssum
  real(kind=pr)                 :: t_cumulative
  logical                       :: prune_all
  real(kind=pr)                 :: maxde,maxdn
  complex(kind=pr), allocatable :: w(:)
  real(kind=pr)                 :: state(7),stateave(2),de(7)
  real(kind=pr)                 :: statnave(2),dn(7)
  real(kind=pr)                 :: dnorm, vnorm
  integer                       :: statn(7)
  integer                       :: ialloc(20)
  integer                       :: ierr
  integer                       :: length
  real(kind=pr)                 :: ecore
  logical                       :: ecnvrgd,ncnvrgd,crit
  complex(kind=pr)              :: eref
  integer                       :: lb4prune, last_ok
  real(kind=pr)                 :: f_boost
  real(kind=pr)                 :: seed
  integer                       :: i_try
  integer                       :: nu_configs
  integer                       :: ltarget, inxcss
  integer                       :: lb4brnch
  real(kind=pr)                 :: compare
  integer                       :: idiag
  complex(kind=pr)              :: eval
  complex(kind=pr)              :: entot
  integer                       :: i_got_hit
  integer                       :: lconv
  integer                       :: lb4exc
  integer                       :: jj, i, n

  call rng_init_clock()

  ! call pbeginf                            ! TCGMSG
  call mpi_init(ierr)                       ! MPI

  ! who am i and how many others are there?
  ! me = nodeid()                             ! TCGMSG
  ! nproc = nnodes()                          ! TCGMSG
  call mpi_comm_rank(MPI_COMM_WORLD, me, ierr)    ! MPI
  call mpi_comm_size(MPI_COMM_WORLD, nproc, ierr) ! MPI

  if(me.eq.0) then
     open(50,file="e_summary")
     write(50,'(24a)') fdate()
     call header
  end if

  call read_params()
  call allocate_memory
  call init(seed,ecore)
  ! initialize variables for diagonalization steps
  c(1)   = 1.0_pr
  length = 1
  dnorm  = 1.0_pr
  if (inflg.eq.0) then
     call h_s_sparse(1,0)
     ! call dump('hdump',h,ijh)
     ! call dump('sdump',s,ijs)
     eref   = h(1)/s(1)
     vnorm  = sqrt(s(1))
     c(1)   = c(1)/vnorm
  else
     call restart(length)
     call chk_list(length,0)
     call h_s_sparse(length,0)
     eref   = h(1)/s(1)
  endif
  write(50,*) "ecore = ", ecore

  if(me.eq.0) write(*,*) '***initialization has been done***'


  lb4prune = length
  last_ok =  length

  !Write info to e_summary
  if(me.eq.0) call initial_info

  !for branching with nproc>1, how many extra configs before an exchange?
  f_boost = bmin

  t_cumulative = 0.0_pr

  call init_bk(length)

  !Divide frac among processors
  frac = frac/nproc

  !This is the loop responsible for iterating the cycles in mcci.  In each cycle,
  !new configurations are appended onto the current CI vector.  The resultant
  !hamiltonian is diagonalised, and poor configurations are pruned from the vector.
  i_try = 1
  do while (i_try.le.maxtry)

     !Timing data is will be written to e_summary
     if(me.eq.0.and.time) call timer(twi,tui,tsi)

     if(me.eq.0) then
        open(50,file='e_summary',status='old',position='append')
        write(50,*)
        write(50,*)
        write(50,*)'Diagonalization',i_try
        write(50,*)
     endif

     if(i_try .lt. maxtry) then

        !The number of new configurations to be appended
        nu_configs = int(frac*real(length, kind=pr))
        if(nu_configs.lt.1) nu_configs=0
        !The number of new configurations must be large enough to
        !satisfy lmin
        if(nu_configs.lt.lmin-length) nu_configs=lmin-length

        if(nproc .gt. 1) then
           ltarget    = length + nu_configs    ! desired length
           inxcss     = int(f_boost*dble(length))
           nu_configs = nu_configs + inxcss    ! configs+excess
        else
           ltarget    = length + nu_configs    ! desired length 
        endif

        !call gen_seed(seed) !not sure why this was put here..  MS

        if(nobrnch_first.and.(i_try.eq.1)) then
           nu_configs=0
           ltarget=length
        endif

        call branch(nu_configs,length,lb4brnch,seed)

        call genealogy(length,lb4brnch,seed)

        call chk_list(length,lb4brnch)

        if(nproc.gt.1) then

           if(me.eq.0) then
              write(50,*)'Boosted CI vector length =',length
              write(50,*)'f boost= ',f_boost
           endif

           compare = real(ltarget-length, kind=pr)/real(lb4brnch, kind=pr)

           ! adjust boost keeping within bounds
           f_boost = min(bmax,max(bmin,f_boost+compare))

           if(length.gt.ltarget) length = ltarget

        endif

        if(me.eq.0 .and. time_all) call timer(twb,tub,tsb)

     endif

     call h_move(length,lb4prune)
     call s_move(length,lb4prune)

     if(length.gt.last_ok) then 
        call h_s_sparse(length,last_ok)
     else
        call h_s_sparse(length,0)
     endif
     !write(50,*) "main: h(1) = ", h(1)
     !write(50,*) "main: h(100) = ", h(100)
     ! call dump('hdump',h,ijh)
     ! call dump('sdump',s,ijs)
     !print*, 'h matrix'
     !do i = 1, 16  !ijh(1)-2
     !   print*,  h(i) ,ijh(i)
     !enddo
     !print*, 's matrix'
     !do i = 1, 120  !ijh(1)-2
     !   print*,  s(i) ,ijs(i)
     !enddo

     if(me.eq.0 .and. time_all) call timer(twg1,tug1,tsg1)

     dwsum  = 0.0
     dusum  = 0.0
     dssum  = 0.0
     if(nodiag) goto 10

     !The davidson routine is responsible for diagonalising the hamiltonian
     print*, 'iteration: ', i_try
     call davidson(length,ieig,idiag)

     ! write(50,*)'edavidson', e(ieig) + ecore
     eval = e(ieig)

     if(me.eq.0) then
        write(50,*)'Iterations to convergence',idiag
        write(50,*)'E   =  ', eval + cmplx(ecore, kind=pr), 'from iter. diag.'
        !write(50,1002) eval+ecore
        write(50,*)'Branched  CI vector length=',length
     endif

     entot = eval+ecore

     call energy(length,eval,dnorm)

     call mxv_sparse(length,h,ijh,c,w)

10   continue

     if(me.eq.0) call energy_eval 

     if(i_try.lt.maxtry)then
        prune_all = .false. 
        if(me.eq.0 .and. time_all) call timer(twp0,tup0,tsp0)
        if( mod(i_try,npfull).eq.0.or.i_try.eq.(maxtry-1).or.i_try.eq.1) prune_all = .true.

        lb4prune = length

        call prune(length,lb4brnch,i_got_hit,prune_all,ieig)

        if(me.eq.0) then
           if(prune_all) then
              write(50,*)'Pruned    CI vector length=',length,&
                   '   *** Full pruning ***'              
              lconv = length
           else
              write(50,*)'Pruned    CI vector length=',length
           endif
        endif


        if(i_got_hit.ne.maxc) then
           last_ok = i_got_hit-1
        else
           last_ok = length                  ! no pruning took place
        endif

        if(nproc.gt.1) then
           lb4exc = length
           do jj=0, nproc-1
              call exc(lb4exc,lb4brnch,length,jj)    !parallel
           enddo
           if(me.eq.0) write(50,*)'Exchanged CI vector length=',length
           call chk_list(length,lb4exc)
        endif

        if(me.eq.0 .and. time_all) call timer(twp,tup,tsp)

        if(me.eq.0) then
           if(time_all) then 
              write(50,*)'BRANCHING                                    total'
              write(50,2004) twb - twi
              write(50,2005) tub - tui
              write(50,2006) tsb - tsi
              write(50,*)'H GENERATION FOR BRANCHED CONFIGURATIONS     total'
              write(50,2004) twg1 - twb
              write(50,2005) tug1 - tub
              write(50,2006) tsg1 - tsb
              write(50,*)'H DIAGONALIZATION average                    total'
              write(50,2001) dwsum/idiag,dwsum
              write(50,2002) dusum/idiag,dusum
              write(50,2003) dssum/idiag,dssum
              write(50,*)'PRUNING                                      total'
              write(50,2004) twp - twp0
              write(50,2005) tup - tup0
              write(50,2006) tsp - tsp0
           endif
        endif
     endif

     if(me.eq.0) then
        if(time) then 
           call timer(twf,tuf,tsf)
           twf = twf - twi
           tuf = tuf - tui
           tsf = tsf - tsi
           write(50,2007) twf
           write(50,2008) tuf
           write(50,2009) tsf
           t_cumulative = t_cumulative + twf
           write(50,2010) t_cumulative
           write(50,'(24a)') fdate()
        endif
        close(50)
     endif

     !The root processor checks for convergence
     if(me.eq.0 .and. i_want_conv) call convergence_test

     if(nproc.gt.1) call MPI_BCAST(maxtry,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

     i_try = i_try + 1

  enddo

1000 format(1x,'Eref= ',f20.14)
1001 format(1x,'Ecor= ',f20.14,' from Davidson diag.')
1002 format(1x,'E   = ',f20.14,' from iter. diag.')
1003 format(1x,'Egnd= ',f20.14,' end of iteration ')
1004 format(1x,'Ecor= ',f20.14,' best guess on a node from cHc',' after pruning') 
1005 format(1x,'Ecor= ',f20.14,' starting guess from cHc')

2001 format(1x,'wall time',f12.3,' s      wall time',f12.3,' s')
2002 format(1x,'user time',f12.3,' s      user time',f12.3,' s')
2003 format(1x,'sys  time',f12.3,' s      sys  time',f12.3,' s')

2004 format(1x,'         ', 12x ,'        wall time',f12.3,' s')
2005 format(1x,'         ', 12x ,'        user time',f12.3,' s')
2006 format(1x,'         ', 12x ,'        sys  time',f12.3,' s')

2007 format(1x,'Cycle WALL time',f12.3,' s')
2008 format(1x,'      USER time',f12.3,' s')
2009 format(1x,'      SYS  time',f12.3,' s')

2010 format(1x,'Cumul. run time',f12.3,' s')

  !  call pend                              ! TCGMSG
  call mpi_finalize(ierr)                   ! MPI


  !End of mcci.f90


  contains

  subroutine allocate_memory
     allocate(e(kmax),             stat=ialloc(1))
     allocate(icij(2,iword,maxc),  stat=ialloc(2))
     allocate(h(maxh),             stat=ialloc(3))
     allocate(s(maxs),             stat=ialloc(4))
     allocate(ijh(maxh),           stat=ialloc(5))
     allocate(ijs(maxs),           stat=ialloc(6))
     allocate(c(maxc),             stat=ialloc(7))
     allocate(hf(kmax*(kmax+1)/2), stat=ialloc(8))
     allocate(sf(kmax*(kmax+1)/2), stat=ialloc(9))
     allocate(b(maxc,kmax),        stat=ialloc(10))
     allocate(irrep(0:irmax-1),    stat=ialloc(11))
     allocate(list(2,maxocc),      stat=ialloc(12))
     allocate(my_pair(2,maxocc),   stat=ialloc(13))
     allocate(icase(16),           stat=ialloc(14))
     allocate(ipoint(max2),        stat=ialloc(15))
     allocate(e1ints(max1),        stat=ialloc(16))
     allocate(e2ints(max2),        stat=ialloc(17))
     allocate(nbpsy(irmax),        stat=ialloc(18))
     allocate(cnorm(maxc),         stat=ialloc(19))
     allocate(w(maxc),             stat=ialloc(20))
     if(any(ialloc /=0)) STOP "Error allocating memory in mcci.f90"
  end subroutine allocate_memory

  !Initial information regarding the branching factor etc. is
  !written to e_summary
  subroutine initial_info
     write(50,*)'Calculating',ieig,'th state in this irrep'
     write(50,*)'Running on',nproc,' nodes'
     write(50,*)'Branching factor      f =',frac
     write(50,*)'Davidson tolerance stop =',davidson_stop
     write(50,*)'Coef.    tolerance cmin =',cmin
     write(50,*)'H        tolerance hmin =',hmin
     if(inflg.ne.0)write(50,*)'RESTARTED'
     close(50)
  end subroutine initial_info

  subroutine energy_eval
     complex(kind=pr)  :: x
     open(40,file='civ_out',form='formatted')
        open(60,file='weight',form='formatted')
        write(60,*) 'Energy = ', eval+ecore
        eval = 0.0
        do i=1,length
           do n=1,nword
              if(n.eq.1) then
                 write(40,*)&
                 !write(40,'(i6,2x,e23.17,2x,i11,2x,i11)')&
                      i,c(i)/sqrt(dnorm),icij(1,n,i),icij(2,n,i)
              else
                 write(40,*)&
                 !write(40,'(33x,i11,2x,i11)')&
                      icij(1,n,i),icij(2,n,i)
              endif
           enddo
           x = c(i)*w(i)/dnorm
           write(60,*) i, x, x/e(ieig)
           eval = eval + x
        enddo
        write(60,*) 'Energy sum',eval, 'Core energy',ecore
        write(60,*) 'Energy check', eval+ecore
        close(40,err=111)
111     continue
        close(60)
  end subroutine energy_eval

  subroutine convergence_test
    implicit none

    integer  :: i, nnn

     if(i_try.eq.1) then
        do i=1,conv_average
           state(i)=0.0
           statn(i)=0
        enddo
        do i=1,2
           stateave(i)=0.0
           statnave(i)=0.0
        enddo
        do i=1,conv_history
           de(i)=0.0
           dn(i)=0.0
        enddo
        ecnvrgd  = .false.
        ncnvrgd = .false.
        open(70,file='convergence',form='formatted')
        write(70,*) '************CONVERGENCE TEST IS RUNNING***********'
        close(70)
     endif

     if(npfull_conv) then
        crit = (mod(i_try,npfull).eq.1).and.(.not.(ecnvrgd.and.ncnvrgd).and.i_try.gt.1)
        nnn = npfull
     else
        crit = (.not.(ecnvrgd.and.ncnvrgd)).and.(i_try.gt.1)
        nnn = 1
     endif

     if(i_try .eq. 1) then
        open(70,file='convergence',status='old',position='append')
        write(70,*) 'Convergence checking will begin after ',&
                    (conv_average+conv_history)*nnn, ' cycles.'
        write(70,*) '**************************************************'
        close(70)
     endif
        
     if (crit) then

        do i=1,conv_average-1
           state(i)=state(i+1)
        enddo
        state(conv_average)= entot
        stateave(1)=stateave(2)
        stateave(2)=0.0
        do i=1,conv_average
           stateave(2)=stateave(2)+state(i)
        enddo
        stateave(2)=stateave(2)/real(conv_average, kind=pr)
        do i=1,conv_history-1
           de(i)=de(i+1)
        enddo
        de(conv_history)=stateave(2)-stateave(1)
        maxde=0.0
        do i=1,conv_history
           if(abs(de(i)).gt.maxde) maxde=abs(de(i))
        enddo

        do i=1,conv_average-1
           statn(i)=statn(i+1)
        enddo
        statn(conv_average)=lconv
        statnave(1)=statnave(2)
        statnave(2)=0.0
        do i=1,conv_average
           statnave(2)=statnave(2)+real(statn(i), kind=pr)
        enddo
        statnave(2)=statnave(2)/real(conv_average, kind=pr)
        do i=1,conv_history-1
           dn(i)=dn(i+1)
        enddo
        if (statnave(1).ne.0) then
           dn(conv_history)=statnave(2)/statnave(1) - 1.0
        else
           dn(conv_history)=999.
        endif
        maxdn=0.0
        do i=1,conv_history
           if(abs(dn(i)).gt.maxdn) maxdn=abs(dn(i))
        enddo

        if (i_try.ge.((conv_average+conv_history)*nnn)) then 
           ecnvrgd = .false.
           if(maxde.lt.conv_thresh_e) ecnvrgd=.true.
           ncnvrgd = .false.
           if(maxdn.lt.conv_thresh_l) ncnvrgd=.true.
           open(70,file='convergence',status='old',position='append')
           write(70,*) 'i_try:',i_try,'ecnvrgd:',ecnvrgd,'ncnvrgd:',ncnvrgd
           write(70,*) 'de:',de(conv_history),'dn:',dn(conv_history)
           ! write(70,*) 'maxde:',maxde,'maxdn:',maxdn
           close(70,err=112)
112        continue
           if(ecnvrgd.and.ncnvrgd) maxtry = i_try + 2
        endif

     endif

  end subroutine convergence_test

end program mcci

subroutine header
  write(50,*)
  write(50,*)'                   m c c i  3.0'
  write(50,*) 
  write(50,*)'              written by J.C. Greer'
  write(50,*)
  write(50,*)'   ============================================='
  write(50,*)'   J.C. Greer, J. Chem. Phys. 103 (1995) p. 1821'
  write(50,*)'   J.C. Greer, J. Comp. Phys. 146 (1998) p. 181 '
  write(50,*)'   ============================================='
  write(50,*)
  return
end subroutine header

