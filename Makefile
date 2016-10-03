SHELL=/bin/csh
GAU_DIR = $(gdvroot)/gdv
GAU_DIRL = $(GAU_DIR)
GAU_DIRA = $(GAU_DIR)
UTIL_NAME = util.a

FLC = flc
#BLAS = -lcxml
#BLAS = -lblas
#BLAS = -lesslp2
BLAS = -lscs_i8
IBM_FC = xlf_r -q64 -qsmp=noauto -qextname -qintsize=8 -qrealsize=8
SGI_FC = f77 -w -i8 -r8 -r8const -mips4 -64 -mp -r10000 -align64	\
-trapuv -Wl,-Xlocal,tstampc_ -Wl,-Xlocal,savedt_
SGI_FC2 = -lfastm
SUN_FC = f95 -fast -xtypemap=real:64,double:64,integer:64	\
-xtarget=native -xarch=v9 -xcache=generic
ALP_FC = f90 -O5 -transform_loops -omp -automatic -i8 -r8 -align	\
dcommons -tune host -trapuv -assume noaccuracy_sensitive -math_library	\
fast -reentrancy threaded
LIN_FC = pgf77 -mp -Mnostdlib
PGILIBS = $(PGI)/linux86/lib/libpgthread.a $(PGI)/linux86/lib/libpgc.a	\
$(PGI)/linux86/lib/libpgftnrtl.a $(PGI)/linux86/lib/libpgc.a		\
$(PGI)/linux86/lib/libpgmp.a $(PGI)/linux86/lib/libpgc.a
LIN_FC2 = /usr/local/lib/libf77blas.a /usr/local/lib/libatlas.a $(PGILIBS)
CRY_FC = f90 
CRY_FC2 = -Wl"-M /dev/null -D DUPENTRY=NOTE -D FORCE=OFF -f indef"
HP_FC = f77 +U77
HP_FC2 = -lveclib -lcps -lpthread -lm
PROFFLAG =
#NUTILM = ../nutil/nutil-mark
#NUTIL = ../nutil/*.o
#NUTILL = ../nutill/*.lo ../nutill/*.o
#NUTILLM = ../nutill

FC0 = pgf90 -g $(FC) $(PROFFLAG)
FC2 = $(CRY_FC2)
FC1 = $(GAU_DIRL)/$(UTIL_NAME) $(FC2)

LAPext = ./lapack_LINUX.a

.SUFFIXES:
.SUFFIXES: .lo .o .F

.F.o:
	rm -f $*.o
	$(MAKE) -f $(GAU_DIRL)/bsd/gdv.make MAKE='$(MAKE)' \
	PROFFLAG='$(PROFFLAG)' $*.o

.F.lo:
	$(MAKE) -f $(GAU_DIRL)/bsd/gdv.make MAKE='$(MAKE)' \
	PROFFLAG='$(PROFFLAG)' $*.lo

all: l535.exe

bdrys.o:
	gau-get bdrys utilam
	make -f $(GAU_DIRL)/bsd/gdv.make MAKE='$(MAKE)' \
	PROFFLAG='$(PROFFLAG)' $*.o
	rm -f bdrys.F

bdam1.o:
	gau-get bdam1 utilam
	make -f $(GAU_DIRL)/bsd/gdv.make MAKE='$(MAKE)' \
	PROFFLAG='$(PROFFLAG)' $*.o
	rm -f bdam1.F

OBJUT = 
OBJPARUT =
OBJLIN = $(GAU_DIRL)/mdutil-linda.o $(GAU_DIRL)/drum.lo	\
	 $(GAU_DIRL)/caldsu.lo $(GAU_DIRL)/coulsu.lo	\
	 $(GAU_DIRL)/prnxtv.lo $(GAU_DIRL)/lindev.lo	\
	 $(GAU_DIRL)/chain.lo $(GAU_DIRL)/linint.lo	\
	 $(GAU_DIRL)/glini1.lo $(GAU_DIRL)/prlin3.lo	\
	 $(GAU_DIRL)/glinlo.lo $(GAU_DIRL)/prrfsu.lo	\
	 $(GAU_DIRL)/glinre.lo $(GAU_DIRL)/prsmsu.lo

$(COBJS):
	$(MAKE) -f $(GAU_DIRL)/bsd/gdv.make MAKE='$(MAKE)' \
	PROFFLAG='$(PROFFLAG)' $*.o

# = link 535 =

MAIN535 = ml535.o

OBJ535 = aobstf.o  binom.o   derictr.o  detmat.o detmt3.o \
	 detmt2.o  dmblck.o  erictr.o evaldg.o evalgg.o evallg.o evalng.o  evalpg.o  \
	 frc1e.o   frc2e.o  fixgaufmt.o fixmat.o  fixmt2.o  \
	 formden.o gaufmt.o  \
	 gauleg.o  indgrd.o  intgrd.o intdrv.o mltdg.o \
	 oneectr.o phfgrad.o rtfact.o\
	 setgrd.o  sptblk.o  sptmat.o\
	 trpgrd.o  trprct.o  wigarr.o  wignerd.o

l535.exe: $(MAIN535) $(OBJ535) $(NUTILM)
	$(FC0) -o l535.exe $(MAIN535) $(OBJ535) $(NUTIL) \
	$(LAPext) $(FC1) $(BLAS)
	chmod o-rx l535.exe
