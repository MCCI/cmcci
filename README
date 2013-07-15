!================================!
!README file for MCCI!           !
!================================!

!==========================!
!Original MCCI Documentation!
!==========================!

Monte Carlo configuration generation computer program for the calculation of
electronic states of atoms, molecules, and quantum dots
(2000) Computer Physics Communications, 131 (1), pp. 142-163.

!====================!
!Compile Instructions!
!====================!


Firstly set path for mcci, $MCCI = whatever_path
This directory contains the source code, a Makefile for
a serial machine, a Makefile for a parallel calculation and
the parameter file, params, details of which may be found in
the accompanying paper.

In the Makefiles, the paths for the LAPACK and BLAS libraries
will be changed, depending on where the user has them installed.
Also, the compiler name and associated compiler
flags may be changed. If the compiler is changed, please note
that the system clock called in the subroutine gen_seed.f may also
need to be changed.

!=============================!
!Keywords for the mcci.in file!
!=============================!

  ieig              =  Eigenvalue to be calculated. Each eigenvalue in each irrep begins with ieig=1 
                       i.e. ieig=1 is lowest energy eigenvalue in an irrep, ieig=2 is first excitation in an irrep, ... 
  n_up              =  Number of spin up (alpha) electrons
  n_dn              =  Number of spin down (beta) electrons
  cmin              =  Coefficient threshold:  Defines the minimum value of CSF's coefficient to be kept following
                       a matrix diagonalization, i.e. the "pruning threshold"
  s                 =  Spin ( units of h_bar ) (principle case M_s=S)
  maxtry            =  The maximum number of cycles (branch, diagonalization, prune) that should be performed
  mo_up(:)          =  Molecular orbital indices for up (alpha) electrons
  mo_dn(:)          =  Molecular orbital indices for down (beta) electrons
  lmin              =  Minimum vector length- branch will always generate lmin CSFs
  npfull            =  A "full prune" step will be performed every npfull cycles
  lkeep             =  Configurations to keep throughout a calculation. The first lkeep CSFs will not be pruned.
  hmin              =  H matrix threshold- elements of the H matrix with absolute values below hmin are not stored
  davidson_stop     =  Maximum number of iterations before stopping for convergence in the Davidson algorithm
  bmin              =  Minimum vector boost (used to adjust Monte Carlo sample size in early stages of a calculation)
  bmax              =  Maximum vector boost (used to adjust Monte Carlo sample size in early stages of a calculation)
  frac              =  Branching ratio: number of new CSFs / number of CSFs in current CI vector
  conv_thresh_e     =  Energy convergence: change in energy between npfull cycles is less that conv_thresh_e after ncheck tests
  conv_thresh_l     =  Vector convergencer: change in vector length between npfull cycles is less than conv_thresh_l for ncheck tests 
                       (Note, the number of tests, n, is defined by the keyword conv_history)
  restart           = .TRUE. or .FALSE. If this is set to true, then, a calculation will be 
                       restarted from a previous calculation and civ_in must be present (i.e. previous runs civ_out file 
                       should be renamed to civ_in).  If set to false, a calculation will begin from the CSF with
                       occupations define in mcci.in
  test              = .TRUE. or .FALSE. Debugging flag
  time              = .TRUE. or .FALSE. Timing information, general
  time_all          = .TRUE. or .FALSE. Timing information, detailed
  nobrnch_first     = .TRUE. or .FALSE. No branching in the first step
  nodiag            = .TRUE. or .FALSE. Used for collecting CSFs only (No matrix diagonalisation is peformed during cycles)
  i_want_conv       = .TRUE. or .FALSE. Stop mcci by using convergence criteria.
  npfull_conv       = .TRUE. or .FALSE. Convergence test only in npfull steps
  conv_average      =  Number of npfull steps to be averaged for convergence
  conv_history      =  Window of npfull steps to be include in convergence checking
  frozen_doubly     =  Indices of orbitals to be frozen (i.e. no excitations are allowed from these molecular orbitals 
                       during the CI calculation). Note frozen orbitals must be doubly occupied.
  mo_active           =  Orbital lists of active orbitals, includes all occupied and (initially) unoccupied orbitals to include 
                       in the CI calculations. Note excludes frozen_doubly and all inactive virtual orbitals
  lref              =  Branch (i.e. generate new CSFs) relative to only the first lref CSFs in a CI vector


!=====================================================================!
! Filling in orbital labels in mcci for D2h symmetry using mpgrad_ints!
!=====================================================================!

From the original control file that is obtained from the geometry relaxed
system, note the numbers for the occupied symmetries, e.g.
A1 1-7
A2 1-3
B1 1-9
B2 1-2

Then run eiger -a and go as high up in energy as required and note the numbers
of orbitals that you want to run for mpgrad. Make sure to leave no gaps
between representations.

Then change A1 - B2 for the new symmetries and run mpgrad_ints twice.

To fill this into mcci, we only need to look at the latest version of the
control file.

Fill in the n_up and n_down as required. Then we must label the symmetry
representations in correct order, so if have:

MPGRAD control     MCCI counting
==============     =============
A1 1-10            1-10
A2 1-5             11-16
B1 1-15            17-...
B2 1-4             .......

The general rule is: for the first representation it stays as is, for the
following :

#orbitals + #electrons -1

And when filling in the mcci.in file , the same rule applies.

