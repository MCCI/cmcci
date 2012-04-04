module davidson
  use commonarrays, only: hf, sf, b, h, ijh, s, ijs, c, e
  use mcci_in
  use precision
  use matrixtools

  implicit none

contains

subroutine davidson_main(length,ieig,idiag)

  implicit none

  integer, intent(in)   :: length
  integer, intent(in)   :: ieig
  integer, intent(out)  :: idiag

  real(kind=pr), parameter         :: eps = 1.0d-18
  integer,       parameter         :: maxit = 100
  !integer,       parameter         :: maxit = 10
  integer                          :: info
  integer                          :: i, j, ici, k, kk, kkk, klim
  complex(kind=pr), allocatable    :: vl(:,:), vr(:,:)
  complex(kind=pr), allocatable    :: hf_unpacked(:,:)
  complex(kind=pr), allocatable    :: sf_unpacked(:,:)
  complex(kind=pr), allocatable    :: alpha(:), beta(:)
  complex(kind=pr)                 :: btemp(maxc,kmax)
  complex(kind=pr)                 :: vnorm, rnorm
  complex(kind=pr)                 :: a(maxc), d(maxc), r(maxc)
  !real(kind=pr)                   :: work(3*kmax)
  complex(kind=pr), allocatable    :: work(:)
  integer                          :: lwork
  real(kind=pr), allocatable       :: rwork(:)
  complex(kind=pr)                 :: dot, dot2, temp
  logical                          :: switch
  complex(kind=pr)                 :: tempvec(kmax)

  !     write(34,*)
  !     write(34,*) 'entering davidson...'

  idiag = 1

  if (ieig > length) then
    write(*,*) "Davidson: Error: ieig cannot be greater than ci vector length"
    write(*,*) "ieig: ", ieig, "    length: ", length
    write(*,*) "Stopping..."
    stop
  end if

111 continue

  klim = min(kmax,length)
  do kk=ieig, klim

    allocate(hf_unpacked(kk,kk), sf_unpacked(kk,kk))
    hf_unpacked = (0.0, 0.0)
    sf_unpacked = (0.0, 0.0)

    call hs_red(kk, hf_unpacked, sf_unpacked, length)

    allocate(alpha(kk), beta(kk))
    alpha = (0.0, 0.0)
    beta  = (0.0, 0.0)

    if (allocated(vl)) then   !! These ugly constructs are needed because
      deallocate(vl)          !! of the use of goto statements further
    end if                    !! down, which should, in time be eliminated.
    allocate(vl(kk,kk))
    vl = (0.0, 0.0)

    if (allocated(vr)) then
      deallocate(vr)
    end if
    allocate(vr(kk,kk))
    vr = (0.0, 0.0)

    lwork = 4*kk
    allocate(work(lwork), rwork(8*kk))
    work = (0.0, 0.0)
    rwork = 0.0

    lwork = -1
    call zggev('V', 'V', kk, hf_unpacked, kk, sf_unpacked, kk, alpha, &
               beta, vl, kk, vr, kk, work, lwork, rwork, info)
    lwork = int(work(1))
    deallocate(work)
    allocate(work(lwork))
    work = (0.0, 0.0)
    call zggev('V', 'V', kk, hf_unpacked, kk, sf_unpacked, kk, alpha, &
                beta, vl, kk, vr, kk, work, lwork, rwork, info)
     
    if(info.ne.0) then
      write(*,*) 'info = ', info
      write(*,*) 'N (kk) = ', kk
      write(*,*) 'hf_unpacked = ', hf_unpacked
      write(*,*) 'sf_unpacked = ', sf_unpacked
      write(*,*) 'idiag = ', idiag
      STOP 'Problems with zggev '
    end if
    deallocate(hf_unpacked, sf_unpacked, work, rwork)

    ! zggev returns the generalized eigenvalues as ratios of alpha(i) and
    ! beta(i). We need to check the conditioning and put the reults in 
    ! the e array.
    do k = 1, kk
      if(abs(beta(k)) .lt.1.0d-10) then
        write(0,*) 'Error in davidson iteration ', kk
        write(0,*) 'k: ', k
        write(0,*) 'alpha(k): ', alpha(k)
        write(0,*) 'beta(k): ', beta(k)
        STOP 'davidson: Dangerously low beta(k)'
      else
        e(k) = alpha(k) / beta(k)
      end if
    end do
    deallocate(alpha, beta)

    ! Now we need to sort the eigenvalues, to make sure we grab the ieig
    ! most significant ones.
    do k = 2, kk
      do j = k, 2, -1
        switch = real(e(j)) < real(e(j-1))      &
               .or. (real(e(j)) == real(e(j-1)) &
               .and. aimag(e(j)) < aimag(e(j-1)))
        if(switch) then
          temp   = e(j)
          e(j)   = e(j-1)
          e(j-1) = temp
          tempvec(1:kk) = vr(:,j)
          vr(:,j)       = vr(:,j-1)
          vr(:,j-1)     = tempvec(1:kk)
          tempvec(1:kk) = vl(:,j)
          vl(:,j)       = vl(:,j-1)
          vl(:,j-1)     = tempvec(1:kk)
          !write(*,*) 'davidson: some swapping happened'
        end if
      end do
    end do

    ! Now we calculate the residual, with Davidson's correction
    r = (0.0_pr, 0.0_pr)

    do k=1, kk
       !     note that a and d are calculated in h_s_reduced,
       !     but are currently not being saved for use 
       !     a = H * b
       call mxv_sparse(length,h,ijh,b(:,k),a)
       !     d = S * b
       call mxv_sparse(length,s,ijs,b(:,k),d)
       do ici=1, length
          r(ici) = r(ici)+vr(k,ieig)*(a(ici)-e(ieig)*d(ici))
       enddo
    enddo

     !     Undo the premultiplication by S: S_inv might be hard to calculate, so...
     do ici=1,length
        r(ici)=r(ici)/s(ici) 
     end do
     !     approximate correction to get closer to the true residual: r is now the
     !     exact residual if S is diagonal which should be nearly true
     !     (S diagonally dominant) - PD. Does not seem to make much difference to
     !     convergence anyway.

     !     residual norm
     vnorm = (0.0_pr, 0.0_pr)
     do ici = 1, length
         ! approximate norm:diag. dom.
        vnorm = vnorm + s(ici)*r(ici)*r(ici)
        !vnorm = vnorm + real(conjg(r(ici))*r(ici), kind=pr)
     enddo
         rnorm = sqrt(vnorm)
!    write(34,*) 'idiag,kk,rnorm=',idiag,kk,rnorm

     if(sqrt(real(rnorm)**2 + aimag(rnorm)**2).lt.davidson_stop) goto 222 ! finished?
     
!    form correction vector and orthogonalize to b_ks
     if(kk.lt.klim) then
        do ici=1, length
           if(abs(r(ici)) .lt. eps) then
              b(ici,kk+1) = (0.0_pr, 0.0_pr)
           else
              b(ici,kk+1) = r(ici)/(e(ieig)*s(ici)-h(ici))
           endif
        enddo
        vnorm = sqrt(dot_product(b(:,kk+1),b(:,kk+1)))
        do ici= 1,length
           b(ici,kk+1) = b(ici,kk+1)/vnorm
        enddo

! orthogonalize b_k vectors               
        do k=kk,1, -1
           dot = dot_product(b(:,k), b(:,kk+1))
           b(:,kk+1) = b(:,kk+1) - dot*b(:,k)
        enddo
        vnorm = sqrt(dot_product(b(:,kk+1), b(:,kk+1)))
        do ici= 1,length
           b(ici,kk+1) = b(ici,kk+1)/vnorm
        enddo
     endif
         
  enddo                         ! if allowed, go back again with one more
                                ! b_k vector and rediagonalise in this new
                                ! space.
     !stop 
!     If we get to here, we've brought in as many new b_k vectors as we
!     can, and we're still not converged. So, form the first ieig 
!     eigenvectors from the last LAPACK diagonalisation, make these the
!     new b_k for k=1...ieig, and forget the other b_k vectors.
      
!     generate new b_k vectors from eigenvectors of hred
      do i=1, maxc
         do k=1, kmax
            btemp(i,k) = (0.0_pr, 0.0_pr)
         enddo
      enddo
      do i= 1, ieig
         do ici=1, length
            do k=1, klim        ! we've just expanded in this many b_k's
               btemp(ici,i) = btemp(ici,i) + vr(k,i)*b(ici,k)
            enddo
         enddo
      enddo
!     We guess these eigenvectors are linearly independent: true if they
!     come from different eigenvalues, and likely to be true other wise
!     if the LAPACK routine is doing its job.
      
!     normalize b_k vectors ! don't know why, or the orthogonalisation either
!     - no difference in space spanned. In fact it means the z matrix, the
!     next times around will not be as close to the identity as it could be
!     (on that part where k1,k2 <= ieig).
      do k=1, ieig
         vnorm = (0.0_pr, 0.0_pr)
         do ici=1, length
            b(ici,k) = btemp(ici,k)
            vnorm  = vnorm + b(ici,k)*b(ici,k)
         enddo
         vnorm = sqrt(vnorm)
         do ici=1, length
            b(ici,k) = b(ici,k)/ vnorm
         enddo
         do ici=length+1, maxc
            b(ici,k) = (0.0_pr, 0.0_pr)
         enddo
      enddo
!     orthogonalize the new b_k vectors
      do kkk=2, ieig
         do k=kkk-1, 1, -1
            dot2 = dot_product(b(:,k), b(:,kkk))
            b(:,kkk) = b(:,kkk) - dot2*b(:,k)
         enddo
         vnorm = sqrt(dot_product(b(:,kkk), b(:,kkk)))
         do ici = 1, length
            b(ici,kkk) = b(ici,kkk)/vnorm
         enddo
      enddo
      
!     We know we're not converged.
      idiag = idiag + 1 
      if(idiag.gt.maxit) STOP 'TOO MANY DAVIDSON ITERATIONS'
      goto 111                  ! start afresh with this older, wiser set
                                ! of 1..ieig b_k's. 
      
 222  continue                  ! jump to here when converged.
      
!     We guess that the next call to davidson will be with a hamiltonian
!     not too different from this one: so the final b_k, k=1..ieig will
!     still be good approximate eigenvectors of the new projected hamiltonian.
!     So, save these to b(ici,k),k=1..ieig, so the first LAPACK
!     diagonalisation next Hamiltonian will (if the Hamiltonian were exactly
!     the same) would give us the same approximate Psi, and the same
!     residual as we have now (ie. ~0 residual). This stops us redoing work
!     next call that's done already. 
!     As it is, if we've run over klim LAPACK calls, we've already saved
!     the b_k, k=1..ieig, to memory at this point, and not changed them
!     since, so we have a not too bad space already, and any subsequent calls
!     will start with this space, so with the residual set to the value at
!     the point where we saved the (Gram-Schmidt orthonormalised) b_k. We
!     could do better, though, and save the b_k (k=1..ieig) at the very
!     end, and so start next time round with this residual (modulo the
!     Hamiltonian change) which should be smaller.

!     generate new b_k vectors from eigenvectors of hred
      do i=1, maxc
         do k=1, kmax
            btemp(i,k) = (0.0_pr, 0.0_pr)
         enddo
      enddo
      do i= 1, ieig
         do ici=1, length
            do k=1, kk        ! we've just expanded in this many b_k's
               btemp(ici,i) = btemp(ici,i) + vr(k,i)*b(ici,k)
            enddo
         enddo
      enddo
!     We guess these eigenvectors are linearly independent: true if they
!     come from different eigenvalues, and likely to be true otherwise
!     if the LAPACK routine is doing its job.



      do k=1, ieig
         do ici=1, maxc         ! better than length, for branching
            b(ici,k) = btemp(ici,k)
         enddo
      enddo
 
!     update c vector            
      do ici=1, maxc
         c(ici) = (0.0_pr, 0.0_pr)
      enddo

      do ici= 1, length
         c(ici) = b(ici,ieig)
      enddo

     if (allocated(vl)) then
       deallocate(vl)
     end if
     if (allocated(vr)) then
       deallocate(vr)
     end if
      
      return
  end subroutine davidson_main

  subroutine hs_red(kk, hred, sred, length)
  implicit none
 
  ! Construct the reduced matrices, as explicit 2D arrays.
  ! Technically, sred will still be symmetric (unlike hred, see
  ! Sommerfeld & Tarantelli, J. Chem. Phys. 112, 2106 (2000)!), but since
  ! LAPACK's zggev only takes explicit 2D arrays, it's less messy to just
  ! construct it alongside hred.
  integer,          intent(in)  :: kk
  integer,          intent(in)  :: length
  complex(kind=pr), intent(out) :: hred(kk, kk)
  complex(kind=pr), intent(out) :: sred(kk, kk)
  
  complex(kind=pr)  :: ai(maxc), di(maxc)
  integer           :: i, j, k

  hred = (0.0_pr, 0.0_pr)
  sred = (0.0_pr, 0.0_pr)
 
  do j=1, kk
    call mxv_sparse(length,h,ijh,b(:,j),ai) !ai now holds the jth column of H*B
    call mxv_sparse(length,s,ijs,b(:,j),di)
    do i = 1, kk
      hred(i,j) = dot_product(b(:,i), ai)
      sred(i,j) = dot_product(b(:,i), di)
    end do
  end do
  end subroutine hs_red
 
end module
