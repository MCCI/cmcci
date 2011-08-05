module commonarrays
  use precision
  use dyn_par
  
  integer ,allocatable     ::   icij(:,:,:)
  integer ,allocatable     ::   ijh(:), ijs(:)
  integer ,allocatable     ::   ifreeze(:)
  integer ,allocatable     ::   iactive(:)
  integer ,allocatable     ::   nbpsy(:)
  integer ,allocatable     ::   irrep(:)
  integer ,allocatable     ::   list(:,:)
  integer ,allocatable     ::   my_pair(:,:)
  integer ,allocatable     ::   icase(:)
  integer ,allocatable     ::   ipoint(:)
  complex  (kind=pr)   ,allocatable     ::   e1ints(:)
  real     (kind=pr)   ,allocatable     ::   e2ints(:)
  complex  (kind=pr)   ,allocatable     ::   wints(:)
  real     (kind=pr)   ,allocatable     ::   cnorm(:)
  complex  (kind=pr)   ,allocatable     ::   e(:)
  complex  (kind=pr)   ,allocatable     ::   hf(:), sf(:)
  complex  (kind=pr)   ,allocatable     ::   b(:,:)
  complex  (kind=pr)   ,allocatable     ::   c(:)
  complex  (kind=pr)   ,allocatable     ::   h(:)
  real     (kind=pr)   ,allocatable     ::   s(:)
  integer                  ::   nword,nfreeze,nactive
  integer                  ::   me, nproc
  integer                  ::   ntotal, n_alpha, n_beta, i_sx2
  integer                  ::   nsym, nbft

end module commonarrays
