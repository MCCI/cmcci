subroutine init(seed,ecore)
  use precision
  use commonarrays 
  use rng
  use mcci_in

  implicit none

  real(kind=pr), intent(inout) :: seed
  real(kind=pr), intent(out)   :: ecore

  integer              :: icount, i
  logical              :: is_frozen(maxbfs)

  call sym_init

  if (me.eq.0) write(50,*)

  if (SCF_integral_filename .eq. 'moints.ascii') then
     if (me.eq.0) write(50,*) 'Getting the one- and two-electron integrals from moints.ascii'
     call get_int(ecore)
  end if

  if (SCF_integral_filename .eq. 'moints.TM') then
     if (me.eq.0) write(50,*) 'Getting the one- and two-electron integrals from moints.TM'
     call get_int_TM(ecore)
  end if

  if (SCF_integral_filename .eq. 'moints.bTM') then
     if (me.eq.0) write(50,*) 'Getting the one- and two-electron integrals from moints.bTM'
     call get_int_bTM(ecore)
  end if

  if (me.eq.0 .and. caps) then
    write(50,*) 'Getting the one-electron CAP integrals from ', wints_filename
    call get_int_w(wints_filename)
  end if

  !     now we know nbft, set up the defaults for the complete active space part,
  !     where if nactive has been set to 0 or left unset, all nonfrozen orbitals are active

  if (nactive == 0) then

     nactive = nbft-nfreeze

     do i = 1,maxbfs
        is_frozen(i) = .FALSE.
     end do

     do i=1,nfreeze
        is_frozen(ifreeze(i)) = .TRUE.
     end do

     icount = 0
     do i = 1,nbft
        if (.not. is_frozen(i)) then
           icount=icount+1
           iactive(icount) = i
        end if
     end do

  end if

  call rng_gen_seed(seed) ! Moved to after the disk reads to increase
                          ! variance of generated sys_clock value, see rng.f90

  return  
end subroutine init
