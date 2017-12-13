# This is the makefile for diablo.
# To compile the code, just type make.  Such an approach makes
# recompilation of the code easy, recompiling only as necessary
# to account for recent changes to the code.
#
# As the user, set the following definitions:

#**********   User Defined Variables Below *********
# Fortran 90 complier to use:
COMPILER = ifort

# Any desired options for the compiler (e.g. -O2, -g, etc.)
#USEROPTS = -O3 -fpp -check bounds -mcmodel medium -shared-intel -heap-arrays -g -traceback
USEROPTS = -O3 -fpp -check bounds -fp-model precise -heap-arrays -mcmodel medium -shared-intel
#USEROPTS = -O3 -fpp -check bounds -mcmodel medium -shared-intel
#USEROPTS = -O3 -fpp -check bounds

# Location where fftw and netcdf (optional) libraries are installed
#LINKDIR=/glade/p/work/dwhitt/fftw2_local/lib

# Location where the optional netcdf include file (netcdf.inc) is installed
#INCLUDEDIR=/glade/p/work/dwhitt/fftw2_local/include

# Option to compile with MPI libraries
PARALLEL = TRUE

# Option to enable the LES model (loads the required variables into memory)
LES = TRUE

# Option to complile with the HDF5 libarary
HDF5 = TRUE

# Option to complile with the flux limited vertical advection
FLUX = FALSE
 
# Option to complile with biogeochemistry 
BIO = TRUE

# **********    END of user definitions ************

ifeq ($(ARCH),ifort)
COMPILER = ifort
USEROPTS = -O3 -fpp
else ifeq ($(ARCH),gfortran)
COMPILER = gfortran
USEROPTS = -O3 -cpp
endif

ifeq ($(LES),TRUE)
LES_o = les.o
else
LES_o = no_les.o
endif

ifeq ($(BIO),TRUE)
BIOOPTS=-DBIO
endif

ifeq ($(FLUX),TRUE)
FLUXOPTS=-DFLUX
endif

ALL2ALL=1

COMPILER = mpif90 
MPI = mpi.o
MPIF = mpi.f
#ifeq ($(HDF5),TRUE)
HDF5_o = hdf5.o 
COMPILER = mpif90 
HDF5OPTS=-DHDF5
#endif

MAIN = diablo.f
#HEADER = header

# Use the parameters to set flags
#ifeq ($(NETCDF),TRUE)
#COMPOPTS = $(USEROPTS) $(HDF5OPTS) -I$(INCLUDEDIR)
#LINKOPTS = -L$(LINKDIR) -lrfftw -lfftw -lnetcdf
#NETCDF_o = netcdf.o
#else
COMPOPTS = $(USEROPTS) $(HDF5OPTS) $(FLUXOPTS) $(BIOOPTS)
LINKOPTS = -L$(LINKDIR) -lfftw2xf_intel
#NETCDF_o = no_netcdf.o
#endif

diablo: $(MAIN) diablo_io.o periodic.o channel.o courant.o $(LES_o) \
	duct.o cavity.o fft.o rand.o set_ics.o filter.o \
        save_stats.o user_rhs.o set_bcs.o bio.o $(MPI) \
	header $(HDF5_o)
	$(COMPILER) $(COMPOPTS) \
        $(MAIN) -o diablo \
	diablo_io.o periodic.o channel.o courant.o $(LES_o) \
	duct.o cavity.o fft.o rand.o set_ics.o filter.o \
        save_stats.o user_rhs.o set_bcs.o bio.o \
	$(MPI) $(LINKOPTS) $(HDF5_o)

diablo_io.o: diablo_io.f bio.o \
        header
	$(COMPILER) $(COMPOPTS) \
         -c diablo_io.f

periodic.o: periodic.f fft.o header
	$(COMPILER) $(COMPOPTS) \
        -c periodic.f

channel.o: channel.f fft.o user_rhs.o \
           $(MPI) header
	$(COMPILER) $(COMPOPTS) \
        -c channel.f

set_ics.o: set_ics.f channel.o fft.o \
           $(MPI) header
	$(COMPILER) $(COMPOPTS) \
        -c set_ics.f

filter.o: filter.f header
	$(COMPILER) $(COMPOPTS) \
        -c filter.f

save_stats.o: save_stats.f channel.o fft.o $(MPI) header
	$(COMPILER) $(COMPOPTS) -c save_stats.f

user_rhs.o: user_rhs.f bio.o fft.o header
	$(COMPILER) $(COMPOPTS) \
        -c user_rhs.f

set_bcs.o: set_bcs.f header
	$(COMPILER) $(COMPOPTS) \
        -c set_bcs.f

bio.o: bio.f fft.o header
	$(COMPILER) $(COMPOPTS) \
        -c bio.f

courant.o: courant.f $(MPI) header
	$(COMPILER) $(COMPOPTS) \
        -c courant.f

ifeq ($(LES),TRUE) 
les.o: les.f fft.o channel.o $(MPI) \
       header header_les grid_def
	$(COMPILER) $(COMPOPTS) \
        -c les.f
else
no_les.o: dummy_code/no_les.f
	$(COMPILER) $(COMPOPTS) \
        -c dummy_code/no_les.f
endif

ifeq ($(PARALLEL),TRUE)
mpi.o: $(MPIF) header
	$(COMPILER) $(COMPOPTS) \
        -c -o mpi.o $(MPIF) 
else
mpi_serial.o: dummy_code/mpi_serial.f header grid_def
	$(COMPILER) $(COMPOPTS) \
        -c dummy_code/mpi_serial.f
endif

header : header_mpi grid_def

header_mpi : grid_mpi

hdf5.o : hdf5.f channel.o $(MPI) courant.o fft.o \
         header
	$(COMPILER) $(COMPOPTS)  \
        -c hdf5.f

hdf5s.o : hdf5s.f fft.o channel.o $(MPI) courant.o \
          header
	$(COMPILER) $(COMPOPTS) \
        -c hdf5s.f

duct.o: duct.f header grid_def
	$(COMPILER) $(COMPOPTS) \
        -c duct.f

cavity.o: cavity.f header grid_def
	$(COMPILER) $(COMPOPTS) \
        -c cavity.f

fft.o:  fft.f header
	$(COMPILER) $(COMPOPTS) \
        -c fft.f

rand.o:  rand.f
	$(COMPILER) $(COMPOPTS) \
        -c rand.f

clean:
	rm -f *.o fort.* *~ diablo core

# Compiler specific notes:
#
# Compilation with Absoft Linux Fortran 77 appears to be impossible, as it
# cannot handle the INTEGER*8 option required by FFTW.  If someone finds
# a way around this, please let me know.
# 
# Compilation with Absoft Linux Fortran 90 is possible, but the option
# -YEXT_NAMES=LCS must be used as one of the link options so the compiler
# can find the lowercase external library function names.
#
# Compilation with Lahey Fortran 95 (lf95) is possible, but there is an
# underscore incompatability with the FFTW libraries, which are compiled
# with g77.  To get around this, you need to go into fft.f and add 
# trailing underscores to the name of every fftw function where they
# appear throughout the code.

