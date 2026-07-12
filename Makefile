# MoCHII Makefile — adapted from MoCafe_v2.00/Makefile (SEDust dropped;
# G0 module set).  Build:  make            ->  MoCHII.x
#                          make DEBUG=1    ->  bounds/traceback build
#************* Choose Compiler *************
UNAME = $(shell uname)

ifeq ($(UNAME), Linux)
    FC      = mpiifort
    ifeq (, $(shell which ifort))
        FC  = mpiifx
    endif
    ifeq (, $(shell which mpiifx))
        FC  = mpif90
    endif
endif
ifeq ($(F90), mpif90)
    FC      = mpif90
endif
ifeq ($(F90), mpiifx)
    FC      = mpiifx
endif
ifeq ($(F90), mpiifort)
    FC      = mpiifort
endif

#---------------------------
SRCDIR  = src
MAIN    = main
FLAGS   = -cpp -DMPI -ISEDust_lib
DEBUG   = 0
HDF5    = 1
SEDUST_LIB = SEDust_lib/libsedust.a

# HDF5 installation prefix (set when HDF5=1).
HDF5_PREFIX ?= /data/opt/hdf5_intel

ifneq ($(HDF5), 0)
   FLAGS    += -DHDF5 -I$(HDF5_PREFIX)/include
   HDF5_LIBS = $(HDF5_PREFIX)/lib/libhdf5_fortran.a $(HDF5_PREFIX)/lib/libhdf5.a -lsz -ldl -lz -lm
endif

ifeq ($(FC), $(filter $(FC), mpiifort mpiifx))
    FFLAGS  = -ipo -O3 -no-prec-div -fp-model fast=2 $(FLAGS) $(Fextra)
    MODFLAG = -module $(SRCDIR)
else ifeq ($(FC), mpif90)
    FFLAGS  = -O3 -ffpe-summary=none $(FLAGS) $(EXTRAFLAG)
    MODFLAG = -J$(SRCDIR) -I$(SRCDIR)
endif

ifeq ($(DEBUG), 1)
   ifeq ($(FC), $(filter $(FC), mpiifort mpiifx))
      FFLAGS  = -check all,noarg_temp_created -fpe0 -debug all -traceback -g -O0 $(FLAGS)
   else
      FFLAGS  = -O0 -fimplicit-none -Wall -Wextra -fcheck=all -fbacktrace $(FLAGS)
   endif
endif

LDFLAGS = $(extra) $(FFLAGS) -lcfitsio -L/usr/local/lib $(HDF5_LIBS) $(SEDUST_LIB) -qopenmp
#*********************************************************************
.SUFFIXES: .f .f90 .o

$(SRCDIR)/%.o: $(SRCDIR)/%.f90
	$(FC) $(FFLAGS) $(MODFLAG) -c -o $@ $<

OBJS	= \
	define.o \
	random_mt.o \
	utility.o \
	memory_mod_mpi.o \
	mathlib.o \
	fitsio_mod.o \
	hdf5io_mod.o \
	iofile_mod.o \
	read_mod.o \
	physics_amr_mod.o \
	octree_mod.o \
	read_generic_amr.o \
	sed_mod.o \
	cellinfo_mod.o \
	jtally_mod.o \
	photo_xsec.o \
	recomb_mod.o \
	gaunt.o \
	cooling_mod.o \
	gas_state_mod.o \
	ion_band_mod.o \
	species_mod.o \
	gas_opacity_mod.o \
	raytrace_amr.o \
	ion_peel_mod.o \
	grid_mod_car.o \
	grid_mod_amr.o \
	amr_refine_mod.o \
	dust_temp_mod.o \
	sedust_mod.o \
	diffuse_mod.o \
	nebcont_mod.o \
	nlevel_mod.o \
	sh95_mod.o \
	lines_mod.o \
	gas_rates_mod.o \
	ion_balance_mod.o \
	thermal_mod.o \
	setup.o \
	main.o

OBJECTS = $(patsubst %.o, $(SRCDIR)/%.o, $(OBJS))

# ion_peel_mod: precise FP (the fast=2 build produced Inf pixels through
# the extreme-dynamic-range peel accumulation; -O0 -check is clean).
$(SRCDIR)/ion_peel_mod.o: FFLAGS += -fp-model precise


# MoCafe policy: no inter-module dependency lists — always a full rebuild
# (module .mod staleness otherwise bites after editing define.f90 etc.).
default: clean
	$(MAKE) MoCHII.x

MoCHII.x: $(OBJECTS)
	$(FC) -o $@ $(OBJECTS) $(LDFLAGS)

clean:
	rm -f $(SRCDIR)/*.o $(SRCDIR)/*.mod MoCHII.x

.PHONY: default clean
