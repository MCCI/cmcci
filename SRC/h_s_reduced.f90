subroutine h_s_reduced(length,kl,ku)
  use commonarrays, only: h, ijh, s, ijs, hf, sf, b
  use mcci_in
  use precision
  use matrixtools

  implicit   none

  integer, intent(in)     :: length
  integer, intent(inout)  :: kl
  integer, intent(in)     :: ku

  integer   :: i, j, k
  !real(kind=pr)       :: ai(maxc), di(maxc) 
  complex(kind=pr)     :: ai(maxc), di(maxc)

  !     matrices are UPPER packed
  !if(kl.le.0) kl=l !! original makes no sense??
  if(kl.le.0) kl=0  !! think this is what they meant
  do j = kl, ku
     do i=1, j

        call mxv_sparse(length,h,ijh,b(:,j),ai)
        call mxv_sparse(length,s,ijs,b(:,j),di)

        hf(i+(j-1)*j/2) = (0.0_pr, 0.0_pr)
        sf(i+(j-1)*j/2) = (0.0_pr, 0.0_pr)
        do k = 1, length   
           hf(i+(j-1)*j/2) = hf(i+(j-1)*j/2) + b(k,i)*ai(k)
           sf(i+(j-1)*j/2) = sf(i+(j-1)*j/2) + b(k,i)*di(k)
        enddo

     enddo
  enddo

  return
end subroutine h_s_reduced
