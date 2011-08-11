#
#=======================================#
#   Makefile to build  m c c i   4      #
#=======================================#
#
       DEST   = ./
 
#------gfortran on mark's desktop----------------
#        FC     = mpif90
#        FFLAGS = -g
#	LFLAGS = -llapack -lblas

#------PGF90 on sfi400.tyndall.ie----------------

        FC     = /usr/local/pgi108/linux86-64/2010/mpi/mpich/bin/mpif90
#        FFLAGS =  -i4 -r8 -O4 
        FFLAGS = -g -C -i4
        LIBS   = -llapack -lblas
   
  
#------PATHF90 on sfi400.tyndall.ie--------------
#         FC     = /opt/mpi/bin/mpif90
#        FFLAGS =  -i4 -r8 -O4 -OPT:Ofast
#         FFLAGS = -g -i4
#        LFLAGS = -llapack -lblas
#         LIBS   = /opt/acml4.0.0/pathscale64/lib/libacml.a
 

OBJS =  mcci.o get_int.o branch.o singles.o doubles.o reorder.o swap.o\
        h_s_sparse.o chk_list.o sym_init.o prune.o init.o\
        energy.o sch0.o sch1.o sch1_0.o sch1_1.o sch1_2.o\
        sch2.o sch2_0.o sch2_1.o sch2_2.o sch2_3.o sch2_4.o sch2_5.o\
        get_case1.o get_case2.o restart.o muHnu.o matrixtools.o \
        wxmxv_sparse.o h_move.o s_move.o timer.o ck.o dump.o   \
        ldump.o genealogy.o davidson.o h_s_reduced.o init_bk.o exc.o \
        get_int_TM.o get_int_bTM.o get_int_W.o precision.o mcci_in.o \
        commonarrays.o rng.o

mcci:   $(OBJS)
	$(FC) $(FFLAGS) $(LFLAGS) -o $(DEST)$@ $(OBJS) $(LIBS) 


clean:
	rm -f *.o *.mod mcci
 
mcci.o:	      mcci_in.o commonarrays.o precision.o matrixtools.o rng.o
get_int.o:              mcci_in.o commonarrays.o               precision.o
init.o:              mcci_in.o mcci_in.o commonarrays.o               precision.o
h_s_sparse.o:           mcci_in.o commonarrays.o               precision.o
h_s_reduced.o:                    commonarrays.o matrixtools.o precision.o
muHnu.o:		mcci_in.o commonarrays.o               precision.o
reorder.o:		mcci_in.o commonarrays.o               precision.o
swap.o:   		          commonarrays.o               precision.o
matrixtools.o:		          commonarrays.o               precision.o
wxmxv_sparse.o:		          commonarrays.o               precision.o
h_move.o:		          commonarrays.o               precision.o
s_move.o:		          commonarrays.o               precision.o
branch.o:		mcci_in.o commonarrays.o               precision.o  rng.o
genealogy.o:		mcci_in.o commonarrays.o               precision.o  rng.o
prune.o: 		mcci_in.o commonarrays.o               precision.o
singles.o:      	mcci_in.o commonarrays.o               precision.o  rng.o
doubles.o:	        mcci_in.o commonarrays.o               precision.o  rng.o
chk_list.o:		          commonarrays.o               precision.o
sym_init.o:		          commonarrays.o               precision.o
energy.o:  		          commonarrays.o               precision.o
sch0.o:			          commonarrays.o               precision.o
sch1.o:			          commonarrays.o               precision.o
sch1_0.o:		          commonarrays.o               precision.o
sch1_1.o:		          commonarrays.o               precision.o
sch1_2.o:		          commonarrays.o               precision.o
sch2.o:	   		          commonarrays.o               precision.o
sch2_0.o:		          commonarrays.o               precision.o
sch2_1.o:		          commonarrays.o               precision.o
sch2_2.o:		          commonarrays.o               precision.o
sch2_3.o:		          commonarrays.o               precision.o
sch2_4.o:		          commonarrays.o               precision.o
sch2_5.o:		          commonarrays.o               precision.o
get_case1.o:		          commonarrays.o               precision.o
get_case2.o:		          commonarrays.o               precision.o
restart.o:		          commonarrays.o               precision.o
dump.o:  		          commonarrays.o               precision.o
ldump.o:  		          commonarrays.o               precision.o
davidson.o:		mcci_in.o commonarrays.o matrixtools.o  precision.o
init_bk.o:		          commonarrays.o               precision.o
get_int_TM.o:		mcci_in.o commonarrays.o               precision.o
get_int_bTM.o:		mcci_in.o commonarrays.o               precision.o
get_int_W.o:		mcci_in.o commonarrays.o               precision.o
precision.o: 
exc.o:                  mcci_in.o commonarrays.o               precision.o
mcci_in.o:	                  commonarrays.o               precision.o 
commonarrays.o:                                                precision.o
timer.o:                                                       precision.o
ck.o:                                                          precision.o
rng.o:     		          commonarrays.o               precision.o

.SUFFIXES:
.SUFFIXES: .o .f .f90

.f.o:
	$(FC) $(FFLAGS) -c $*.f

.f90.o:
	$(FC) $(FFLAGS) -c $*.f90



