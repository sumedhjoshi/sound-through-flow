# parameters governing how we build the model.  by default, assume that we're
# not explicitly debugging and want things to run as fast as possible.
DEBUG          ?= yes
OPTIMIZE       ?= yes

# set the default Fortran 90 compiler, if nothing else has been defined.
#
# NOTE: this should be the compiler that is wrapped, not the wrapper itself.
#       that is, use "gfortran", "ifort", "pgf90", etc and not "mpif90",
#       "scorep", etc.
FC90           ?= gfortran

# set the default Fortran 90 compiler wrapper, if nothing else has been
# defined.
#
# NOTE: this should be the wrapper around the Fortran 90 compiler, and not
#       the Fortran 90 compiler itself.  that is, use "mpif90", "scorep", etc,
#       and not "gfortran", "ifort", "pgf90", etc.
FC90_WRAPPER   ?= mpif90

# directory where ACML is installed.  a multi-threaded implementation is
# preferred over the serial given that faster is typically better for
# development.
ACML_ROOT      ?= /opt/acml5.3.1/gfortran64_mp

# directory where OpenBLAS' libraries are located.
#
# NOTE: this is the location of the shared and static libraries, not the
# 	    parent directory to the libraries.
OPENBLAS_ROOT ?= /usr/lib/openblas-base

# use a a custom directory to search for the HDF5.mod file if requested,
# otherwise, assume it is in the system path.
ifeq (,$(HDF5_INCLUDE))
HDF5_INCLUDE := /usr/include
endif

F90FLAGS += -I$(HDF5_INCLUDE)
LIBS     += -lhdf5hl_fortran -lhdf5_fortran -lhdf5

# source code for the model.
SEMA_SOURCE     = apply_dirichlet_conditions.f90 \
                  apply_wave_operator.f90 \
                  check_consistency.f90 \
                  compute_divergence.f90 \
                  compute_dxi.f90 \
                  compute_gradient.f90 \
                  compute_laplacian.f90 \
                  derv.f90 \
                  display_inputs.f90 \
                  enforce_continuity.f90 \
                  estimate_memory_footprint.f90 \
                  gll.f90 \
                  internal_mesher.f90 \
                  main.f90 \
                  mod_derivatives.f90 \
                  mod_fields.f90 \
                  mod_geom.f90 \
                  mod_io.f90 \
                  mod_constants.f90 \
                  mod_legendre.f90 \
                  mod_mesh_deformation_maps.f90 \
                  mod_options.f90 \
                  mod_parallel_linalg.f90 \
                  mod_timing.f90 \
                  notify.f90 \
                  permutations.f90 \
                  quad.f90 \
                  read_initfile_data.f90 \
                  read_input_file.f90 \
                  seconds_since_epoch.f90 \
                  setup_deformed_derivatives.f90 \
                  sync_ranks.f90 \
                  write_fieldfile.f90

SEMA_OBJECTS = $(patsubst %.f90, %.o, $(SEMA_SOURCE))

# build a list of Fortran module files from the source file list.  take care
# to translate all of the file names into lower case since the module files
# generated with only lower case.
SEMA_MODULES := $(patsubst mod_%.f90, %.mod, $(filter mod_%.f90, \
                                                 $(shell echo $(SEMA_SOURCE) | tr 'A-Z' 'a-z')))

# the woodbury model is built from everything in the project.
SEMA_CMD     = sem_acoustics

# global, generic variables that contain each of the major targets in this
# Makefile.
SOURCE                   = $(SEMA_SOURCE) $(VTK_SOURCE)
OBJECTS                  = $(SEMA_OBJECTS) $(VTK_OBJECTS)
EXECUTABLES              = $(SEMA_CMD)
LIBRARIES                =
MODULES                  = $(SEMA_MODULES) $(VTK_MODULES)

# list of things that are created during development that can be removed, and
# easily regenerated, if necessary.
CLEAN_LIST               = $(OBJECTS) $(EXECUTABLES) $(LIBRARIES) $(MODULES) \
                           $(IFORT_GENMOD_FILES)

# determine the F90 compilation and linking flags.
#
# NOTE: we use the $(filter) function throughout the compiler, and wrapper,
#       selections so that custom versions of individual tools do not require
#       special support.  in particular, things like 'gfortran-4.8.2' will
#       match 'gfortran'.
ifneq (,$(filter gfortran, $(FC90)))

# gfortran does not have BLAS/LAPACK implementations, so select which versions
# to link against.
#
# NOTE: if ACML or OpenBLAS are requested, ensure that it can be found both at
#       link-time as well as run-time via -rpath.
ifeq ($(BLAS), ACML)
LIBS += -L$(ACML_ROOT)/lib -lacml -Wl,-rpath,$(ACML_ROOT)/lib
else ifeq ($(BLAS), ACML_OMP)
LIBS += -L$(ACML_ROOT)/lib -lacml_mp -Wl,-rpath,$(ACML_ROOT)/lib
else ifeq ($(BLAS), OPENBLAS)
LIBS += -L$(OPENBLAS_ROOT) -Wl,-rpath,$(OPENBLAS_ROOT) -lopenblas -llapack
else
# default to the system's implementation which is likely to be slow.
LIBS += -lblas -llapack
endif

# compile all source assuming 64-bit floating point and integers.  don't set a
# maximum line length.  turn on all warnings so that questionable code is noticed
# early on.
#F90FLAGS  = -c -fdefault-real-8 -ffixed-line-length-none -Wall -Wno-maybe-uninitialized #-fdefault-integer-4
F90FLAGS  += -c -fdefault-real-8 -ffixed-line-length-none -std=f2003 -fall-intrinsics -Wall -Wno-maybe-uninitialized #-fdefault-integer-4

ifeq ($(PEDANTIC),yes)
# if we're grooming the code for nits, turn on a slew of warnings that aren't
# otherwise enabled because they generate a lot of noise with bad code.
F90FLAGS += -Wconversion -pedantic -Warray-temporaries -Wcharacter-truncation -Wsurprising -Wintrinsic-shadow -Wextra -pedantic
endif

# generate the most verbose debugging information we can in compiled objects,
# regardless of the request for debugging.  this allows for improved profiling
# at the cost of some extra diskspace.
F90FLAGS += -ggdb -g3

ifneq ($(DEBUG),no)
# enable array bounds checking and ensure that invalid floating point
# operations signal an exception.
F90FLAGS += -fbounds-check -ffpe-trap=invalid -fbacktrace
endif

ifeq ($(OPTIMIZE),yes)
# enable fast math optimizations (which may not be strictly IEEE-754 compliant),
# unroll loops where it makes sense.
F90FLAGS += -O3 -funroll-loops -march=native -mtune=native
endif

ifeq ($(OPTIMIZE),mac)
# enable fast math optimizations (which may not be strictly IEEE-754 compliant),
# unroll loops where it makes sense.
# this optimization flag is for mac only, and it removes the march and native
# flags which break compilation on mac os.
F90FLAGS += -O3 -funroll-loops
endif

ifneq ($(OPENMP),no)
# if we haven't been requested to disable OpenMP, enable it during compilation
# and linking.
F90FLAGS += -fopenmp
LIBS     += -fopenmp
endif

ifneq ($(PLUGIN),)
# the GCC compiler suite supports plugins used for compilation.  pass it onto
# the compilation flags if the user has supplied one.
F90FLAGS += $(PLUGIN)
endif

else # FC90 !~ gfortran

ifneq (,$(filter ifort, $(FC90)))
# BLAS and LAPACK are provided by Intel's Math Kernel Library (MKL) when using
# Intel's Fortran compiler.
LIBS = -r8 -i8 -mkl=sequential

# compile all source assuming 64-bit floating point and integers.  assume ANSI
# aliasing rules, and warn about everything the compiler is concerned about.
F90FLAGS = -c -double-size 64 -real-size 64 -integer-size 64 -r8 -i8 -ansi-alias -warn all -mkl=sequential

ifeq ($(PEDANTIC),yes)
# if we're grooming the code for nits, turn on additional warnings that aren't
# otherwise enabled because they generate a lot of noise with bad code.
F90FLAGS += -debug all
endif

# generate the most verbose debugging information we can in compiled objects,
# regardless of the request for debugging.  this allows for improved profiling
# at the cost of some extra diskspace.
F90FLAGS += -g

ifneq ($(DEBUG),no)
# enable array bounds checking and ensure that invalid floating point
# operations signal an exception.
F90FLAGS += -fpe0 -debug full -check bounds
endif

ifeq ($(OPTIMIZE),yes)
# optimize for the host architecture and align arrays to 32-byte boundaries for
# maximum performance wrt vectorization.  enable inter-procedural optimizations
# across all translation units, unroll loops with a compiler-detected trip
# count, and use OpenMP for threading.
F90FLAGS += -xHost -align array32byte -ipo -unroll
endif

ifneq ($(OPENMP),no)
# if we haven't been requested to disable OpenMP, enable it during compilation
# and linking.
F90FLAGS += -openmp
LIBS     += -openmp
endif

# ifort explicitly builds interfaces for each sub-routine and function it
# encounters so that they can be checked.  attempt to keep track of them.
#
# NOTE: we cannot (easily) build a list of the "*__genmod.*" files since there
#       is one per sub-routine/function.
IFORT_GENMOD_FILES := *__genmod.f90 *__genmod.mod

else # FC90 !~ ifort
$(warning Unknown F90 compiler - $(FC90)!)
endif
endif

# determine which Fortran 90 wrapper we're using.
ifneq (,$(filter mpif90, $(FC90_WRAPPER)))
# nothing is needed for MPI.  all compilation flags are passed directly to
# $(FC90).
else # FC90_WRAPPER !~ mpif90
ifneq (,$(filter scorep, $(FC90_WRAPPER)))
override FC90_WRAPPER += --mpi --openmp mpif90

# Score-P assumes that OpenMP directives always generate OpenMP-based code.
# since we only enable OpenMP when requested, we need to turn it on in the
# non-requested case so that we do not get cryptic link errors about
# mismatches in thread-level storage (TLS)-related symbols.
ifneq ($(OPENMP),yes)

# NOTE: unfortunately, we have to handle each of the different types of
#       compilers so that we ensure OpenMP is properly enabled.
ifneq (,$(filter gfortran, $(FC90)))
F90FLAGS += -fopenmp
LIBS     += -fopenmp
else
ifneq (,$(filter ifort, $(FC90)))
F90FLAGS += -openmp
LIBS     += -openmp
endif
endif
endif

else
$(warning Unknown F90 wrapper - $(FC90_WRAPPER))
endif
endif

# by default, build the model.
all: $(SEMA_CMD)

# link the executable.
$(SEMA_CMD): $(SEMA_OBJECTS)
	$(FC90_WRAPPER) -o $@ $^ $(LIBS)

%.o: %.f90
	$(FC90_WRAPPER) -o $@ $(F90FLAGS) $<

%.o: %.F90
	$(FC90_WRAPPER) -o $@ $(F90FLAGS) $<

clean clena:
	-rm -rf $(CLEAN_LIST)

distclean distclena:
	-rm -rf $(DISTCLEAN_LIST)

# explicit dependencies between objects and the modules they require for
# compilation.
apply_dirichlet_conditions.o: mod_constants.o
apply_wave_operator.o: mod_constants.o mod_fields.o mod_derivatives.o
check_consistency.o: mod_constants.o mod_parallel_linalg.o mod_geom.o
compute_divergence.o: mod_constants.o mod_legendre.o mod_mesh_deformation_maps.o mod_derivatives.o
compute_dxi.o: mod_constants.o mod_legendre.o mod_timing.o
compute_gradient.o: mod_constants.o mod_legendre.o mod_mesh_deformation_maps.o mod_derivatives.o
compute_laplacian.o: mod_constants.o mod_legendre.o mod_mesh_deformation_maps.o mod_derivatives.o
derv.o:
display_inputs.o: mod_constants.o
enforce_continuity.o: mod_constants.o mod_timing.o
estimate_memory_footprint.o: mod_constants.o
gll.o:
internal_mesher.o: mod_constants.o mod_legendre.o mod_geom.o
main.o: mod_constants.o mod_geom.o mod_legendre.o mod_mesh_deformation_maps.o mod_options.o mod_parallel_linalg.o mod_fields.o mod_timing.o
mod_constants.o:
mod_derivatives.o: mod_legendre.o mod_mesh_deformation_maps.o mod_constants.o
mod_fields.o:
mod_geom.o:
mod_io.o:
mod_legendre.o:
mod_mesh_deformation_maps.o:
mod_options.o:
notify.o:
permutations.o:
quad.o:
read_initfile_data.o: mod_constants.o mod_fields.o mod_geom.o
read_input_file.o: mod_constants.o mod_options.o
seconds_since_epoch.o:
setup_deformed_derivatoes.o: mod_mesh_deformation_maps.o mod_constants.o mod_geom.o
sync_ranks.o: mod_constants.o
write_fieldfile.o: mod_constants.o mod_io.o mod_geom.o mod_fields.o
