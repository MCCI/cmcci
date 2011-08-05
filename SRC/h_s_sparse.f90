subroutine h_s_sparse(length,llast)
  use commonarrays
  use dyn_par
  use precision
  implicit none

  integer, intent(in)  :: length
  integer, intent(in)  :: llast

  complex(kind=pr)     :: eij
  integer              :: ici, jci
  integer              :: n_2p
  real(kind=pr)        :: s_overlap, ck
  integer              :: lnew, k, l


  !     diagonals
  do ici = llast+1, length
     call muHnu(ici,ici,eij,n_2p,s_overlap)
     cnorm(ici) = sqrt( abs(ck(i_sx2,n_2p,0)) )
     s(ici) = ck(i_sx2,n_2p,0)
     h(ici) = eij
  enddo
  ijh(1) = length + 2
  ijs(1) = length + 2

  if(llast.gt.0) then 
     !        number of new configurations
     lnew = length - llast
     !        set index to last location to resume generating h
     k = ijh(llast+1) - 1
     l = ijs(llast+1) - 1
  else
     k = length + 1
     l = length + 1
  endif

  do ici = llast+1, length
     do jci = 1, ici-1

        call muHnu(ici,jci,eij,n_2p,s_overlap)

        if( abs( eij ).ge.hmin ) then
           k = k+1
           if(k.gt.maxh) STOP 'gen_h_s: exceeded matrix storage h'
           h(k) = eij
           ijh(k) = jci
        endif

        if( s_overlap .ne. 0.0 ) then
           l = l+1
           if(l.gt.maxs) STOP 'gen_h_s: exceeded matrix storage s'
           s(l) = s_overlap
           ijs(l) = jci
        endif



     enddo

     ijh(ici+1) = k+1
     ijs(ici+1) = l+1

  enddo
  return
end subroutine h_s_sparse
