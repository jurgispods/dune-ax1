noinst_PROGRAMS = test_multidomain \
						test_channels \
						test_hdf5 \
						test_stuff \
						test_geometrywrapper \
						test_iterators \
						test_acme2_cyl_par \
						test_mpi

test_multidomain_SOURCES = test_multidomain.cc
test_channels_SOURCES = test_channels.cc
test_hdf5_SOURCES = test_hdf5.cc
test_stuff_SOURCES = test_stuff.cc
test_geometrywrapper_SOURCES = test_geometrywrapper.cc
test_iterators_SOURCES = test_iterators.cc
test_acme2_cyl_par_SOURCES = test_acme2_cyl_par.cc
test_mpi_SOURCES = test_mpi.cc

# Uncomment this to switch from clang++ to g++
#CXX=g++ -std=c++11
CXX=clang++ -std=c++11

AM_CPPFLAGS = @AM_CPPFLAGS@ $(SUPERLU_CPPFLAGS)
AM_LDFLAGS = @AM_LDFLAGS@ $(SUPERLU_LIBS)

test_hdf5_CPPFLAGS = $(AM_CPPFLAGS) $(HDF5_INCLUDES)
test_hdf5_LDFLAGS = $(AM_LDFLAGS) $(HDF5_LIBS)

test_geometrywrapper_LDADD = $(SUPERLU_LIBS) $(UG_LDFLAGS) $(UG_LIBS) $(DUNE_LDFLAGS) $(DUNE_LIBS)
test_geometrywrapper_CXXFLAGS = $(AM_CPPFLAGS) $(UG_CPPFLAGS)

test_acme2_cyl_par_LDADD = $(HDF5_LIBS) $(DUNE_LIBS) $(DUNEMPILIBS) $(SUPERLU_LIBS) $(UG_LIBS) $(LDADD) 
test_acme2_cyl_par_LDFLAGS = $(AM_LDFLAGS) $(DUNEMPILDFLAGS) $(HDF5_LDFLAGS)
test_acme2_cyl_par_CXXFLAGS = $(AM_CPPFLAGS) $(UG_CPPFLAGS)

test_mpi_LDADD = $(HDF5_LIBS) $(DUNE_LIBS) $(DUNEMPILIBS) $(SUPERLU_LIBS) $(UG_LIBS) $(LDADD) 
test_mpi_LDFLAGS = $(AM_LDFLAGS) $(DUNEMPILDFLAGS) $(HDF5_LDFLAGS)
test_mpi_CXXFLAGS = $(AM_CPPFLAGS) $(UG_CPPFLAGS)

test_acme2_cyl_par_CPPFLAGS = $(AM_CPPFLAGS) $(DUNEMPICPPFLAGS) $(UG_CPPFLAGS) $(HDF5_CPPFLAGS) -DUSE_CYLINDER_COORDS -DAX1_PARALLEL #-DNDEBUG

# Override default CXXFLAGS
CXXFLAGS = -g -O3 -Wfatal-errors -Wno-tautological-compare #-Wall

#CXX = gccfilter -c -a g++ -std=c++0x
#CC  = "/home/jpods/bin/gccfilter/gccfilter -c -r -w gcc"

# don't follow the full GNU-standard
# we need automake 1.5
AUTOMAKE_OPTIONS = foreign 1.5

# pass most important options when "make distcheck" is used
DISTCHECK_CONFIGURE_FLAGS = --with-dune-common=$(DUNE_COMMON_ROOT) --with-dune-istl=$(DUNE_ISTL_ROOT) --with-dune-grid=$(DUNE_GRID_ROOT) --with-dune-localfunctions=$(DUNE_LOCALFUNCTIONS_ROOT) --with-dune-pdelab=$(DUNE_PDELAB_ROOT)  CXX="$(CXX)" CC="$(CC)"

include $(top_srcdir)/am/global-rules
-include $(top_srcdir)/git.mk
