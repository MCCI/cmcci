subroutine init_bk(ieig,length,inflg)
  use commonarrays, only: c, b
  use dyn_par
  use precision

  implicit none

  integer, intent(in)  :: ieig
  integer, intent(in)  :: length
  integer, intent(in)  :: inflg

  integer          :: i, j
  integer          :: k, kk
  complex(kind=pr) :: dot
  real(kind=pr)    :: vnorm

  do j=1, kmax
     do i= 1, maxc
        b(i,j) = (0.0_pr, 0.0_pr)
        if (i.eq.j) b(i,j) = (1.0_pr, 0.0_pr)
     enddo
  enddo

  if (inflg .ne. 0) then
     do i= 1, length
        b(i,ieig) = c(i)     ! copy input civ_in vector to b(:,ieig)
        !            b(i,1) = c(i)        ! copy input civ_in vector to b(:,1)
     enddo
  end if

  ! Orthogonalize b_k vectors, keeping in mind they are complex
  ! and that the dot_product intrinsic automatically conjugates
  ! its first argument.
  
  do kk = 2, kmax
    do k = kk-1, 1, -1
      dot = dot_product(b(:,k), b(:,kk))
      b(:,kk) = b(:,kk) - dot*b(:,k)
    end do
    vnorm = sqrt(real(dot_product(b(:,kk),b(:,kk)), kind=pr))
    b(:,kk) = b(:,kk)/vnorm
  end do


  return
end subroutine init_bk
