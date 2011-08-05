subroutine prune(length,llast,i_got_hit,prune_all,ieig)
  use commonarrays, only: icij, nword, c, b, s
  use mcci_in
  use precision

  implicit none

  integer, intent(inout)  :: length
  integer, intent(inout)  :: llast
  integer, intent(out)    :: i_got_hit
  logical, intent(in)     :: prune_all
  integer, intent(in)     :: ieig

  integer           :: lshift
  real(kind=pr)     :: sum_diag_s
  real(kind=pr)     :: cutoff
  integer           :: lower_limit
  integer           :: n, i, k, kkk, ici
  complex(kind=pr)  :: dot
  real(kind=pr)     :: vnorm

  !     find the small contributions in list and write array contigous 

  lshift = 0
  i_got_hit = maxc

  sum_diag_s = 0.0
  do,i=1,length
     sum_diag_s = sum_diag_s + real(conjg(c(i))*c(i))*s(i)
  enddo

  !     cutoff = cmin*dsqrt(dnorm)
  cutoff = cmin*sqrt(sum_diag_s)

  if (lref.eq.0) then
     lower_limit = llast+1
     if (prune_all) lower_limit = lkeep+1   
  else
     lower_limit = llast+1
     if (prune_all) lower_limit = lref+1   
  endif

  do ici = lower_limit, length
     !        if(abs(c(ici)).lt.cmin) then         ! compare unnormalized
     !        if(abs(c(ici)).lt.cutoff) then       ! compare normalized
     if(abs(c(ici))*sqrt(s(ici)).lt.cutoff) then  ! new
        if(i_got_hit.eq.maxc) i_got_hit = ici
        lshift = lshift + 1
        if(ici.le.llast) llast = llast - 1
     else
        if(lshift .gt. 0) then
           do n=1,nword
              icij(1,n,ici-lshift) = icij(1,n,ici)
              icij(2,n,ici-lshift) = icij(2,n,ici)
           enddo
           c(ici-lshift)      = c(ici)
           do k=1, ieig
              b(ici-lshift,k) = b(ici,k)
           enddo
        endif
     endif
  enddo

  !     set end buffer elements to zero
  do ici = length - lshift + 1, length
     c(ici) = (0.0_pr, 0.0_pr)
     do k=1, kmax
        b(ici,k) = (0.0_pr, 0.0_pr)
     enddo
     do n=1,nword
        icij(1,n,ici) = 0
        icij(2,n,ici) = 0
     enddo
  enddo

  !     re-orthogonalize the pruned b_k vectors
  do kkk=2, ieig
     do k=kkk-1, 1, -1
        dot = dot_product(b(:,k), b(:,kkk))
        do ici = 1, length
           b(ici,kkk) = b(ici,kkk) - dot*b(ici,k)
        enddo
     enddo
     vnorm = sqrt(real(dot_product(b(:,kkk),b(:,kkk)), kind=pr))
     do ici = 1, length
        b(ici,kkk) = b(ici,kkk)/vnorm
     enddo
  enddo

  length = length - lshift

  return
end subroutine prune
