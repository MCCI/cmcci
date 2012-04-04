subroutine init_bk(length)
  use commonarrays, only: c, b
  use mcci_in
  use precision

  implicit none

  integer, intent(in)  :: length

  integer          :: i, j
  integer          :: k, kk
  complex(kind=pr) :: dot
  complex(kind=pr)    :: vnorm

  do j=1, kmax
     do i= 1, maxc
        b(i,j) = (0.0_pr, 0.0_pr)
        if (i.eq.j) b(i,j) = (1.0_pr, 0.0_pr)
     enddo
  enddo

  if (inflg .ne. 0) then
     if (ieig == 1) then
       do i= 1, length
          b(i,1) = c(i)     ! copy input civ_in vector to b(:,ieig)
       end do
     end if
  end if

  ! Orthogonalize b_k vectors, keeping in mind they are complex
  ! and that the dot_product intrinsic automatically conjugates
  ! its first argument.
  
  do kk = 2, kmax
    do k = kk-1, 1, -1
      dot = dot_product(b(:,k), b(:,kk))
      b(:,kk) = b(:,kk) - dot*b(:,k)
    end do
    vnorm = sqrt(dot_product(b(:,kk),b(:,kk)))
    b(:,kk) = b(:,kk)/vnorm
  end do

  return
end subroutine init_bk
