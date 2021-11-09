noinst_PROGRAMS = acme2_cyl_par \
						acme2_cyl_par_myelin

#noinst_PROGRAMS = acme2_2d \
#						acme2_cyl_par \
#						acme2_cyl_par_mme \
#						acme2_cyl_par_laplace \
#						acme2_cyl_par_mori \
#						aggregate_fluxes \
#						eval_grids

#example01_SOURCES = example01.cc
#example02_SOURCES = example02.cc
#example03_SOURCES = example03.cc
#membrane1_SOURCES = membrane1.cc
#membrane2_SOURCES = membrane2.cc
#membrane3_SOURCES = membrane3.cc
#membrane4_SOURCES = membrane4.cc
#membrane5_SOURCES = membrane5.cc
#acme0_operator_split_SOURCES  = acme0.cc
#acme0_fully_implicit_SOURCES 	= acme0.cc
#acme1_SOURCES    = acme1.cc
#nemo1_SOURCES    = nemo1.cc
#acme1MD_SOURCES  = acme1MD.cc
#acme2_SOURCES    = acme2.cc
acme2_2d_SOURCES    = acme2_cyl_par.cc
acme2_cyl_par_SOURCES    = acme2_cyl_par.cc
acme2_cyl_par_mme_SOURCES    = acme2_cyl_par.cc
acme2_cyl_par_myelin_SOURCES    = acme2_cyl_par.cc
acme2_cyl_par_laplace_SOURCES    = acme2_cyl_par.cc
acme2_cyl_par_mori_SOURCES    = acme2_cyl_par.cc
aggregate_fluxes_SOURCES    = aggregate_fluxes.cc
eval_grids_SOURCES    = eval_grids.cc


# Choose your compiler here (g++ or clang++)
CXX=g++ -std=c++11 
#CXX=clang++ -std=c++11
#CXX=/opt/local/bin/mpic++ -std=c++11

# Use clang++ with gcc stdlib (Wrapper by Steffen); note: had to deactivate port 'hdf5-18' from MacPorts do get this to work, as otherwise there are two hdf5.hh in the include path and the one from MacPorts takes preference!
#CXX=clang++-gcc -std=c++11

#CXX = gccfilter -c -a g++ -std=c++0x
#CC  = "/home/jpods/bin/gccfilter/gccfilter -c -r -w gcc"

# Uncomment this for gSTLfilt'ered output!
#CXX  = "/home/jpods/bin/gfilt" -std=c++0x
#CC  = "/home/jpods/bin/gfilt" -std=c++0x


#AM_CPPFLAGS = @AM_CPPFLAGS@ $(SUPERLU_CPPFLAGS)
# Activate ISTL container checks
AM_CPPFLAGS = @AM_CPPFLAGS@ $(SUPERLU_CPPFLAGS) -DDUNE_ISTL_WITH_CHECKING
# Uncomment this for productive runs
#AM_CPPFLAGS = @AM_CPPFLAGS@ $(SUPERLU_CPPFLAGS) -DNDEBUG -UDUNE_ISTL_WITH_CHECKING

AM_LDFLAGS = @AM_LDFLAGS@ $(DUNE_LDFLAGS) $(SUPERLU_LDFLAGS)
#AM_LDFLAGS = @AM_LDFLAGS@ $(DUNE_LDFLAGS) $(SUPERLU_LDFLAGS) -rdynamic

#acme2_LDADD = $(SUPERLU_LIBS) $(UG_LIBS) $(DUNE_LIBS)
#acme2_LDFLAGS = $(AM_LDFLAGS) $(UG_LDFLAGS)
#acme2_CXXFLAGS = $(AM_CPPFLAGS) $(UG_CPPFLAGS)

acme2_2d_LDADD = $(HDF5_LIBS) $(DUNE_LIBS) $(DUNEMPILIBS) $(SUPERLU_LIBS) $(LDADD)
acme2_2d_LDFLAGS = $(AM_LDFLAGS) $(DUNEMPILDFLAGS) $(HDF5_LDFLAGS)
acme2_2d_CXXFLAGS = $(AM_CPPFLAGS) $(UG_CPPFLAGS)

# Commented version with UG flags
#acme2_cyl_par_LDADD = $(HDF5_LIBS) $(DUNE_LIBS) $(DUNEMPILIBS) $(SUPERLU_LIBS) $(UG_LIBS) $(LDADD) 
#acme2_cyl_par_LDFLAGS = $(AM_LDFLAGS) $(UG_LDFLAGS) $(DUNEMPILDFLAGS) $(HDF5_LDFLAGS)
#acme2_cyl_par_CXXFLAGS = $(AM_CPPFLAGS) $(UG_CPPFLAGS)

# Remove UG flags as long as we don't need it; if you need it, use above flags
acme2_cyl_par_LDADD = $(HDF5_LIBS) $(DUNE_LIBS) $(DUNEMPILIBS) $(SUPERLU_LIBS) $(LDADD) 
acme2_cyl_par_LDFLAGS = $(AM_LDFLAGS) $(DUNEMPILDFLAGS) $(HDF5_LDFLAGS)
acme2_cyl_par_CXXFLAGS = $(AM_CPPFLAGS)

acme2_cyl_par_mme_LDADD = $(acme2_cyl_par_LDADD)
acme2_cyl_par_mme_LDFLAGS = $(acme2_cyl_par_LDFLAGS)
acme2_cyl_par_mme_CXXFLAGS = $(acme2_cyl_par_CXXFLAGS)

acme2_cyl_par_myelin_LDADD = $(acme2_cyl_par_LDADD)
acme2_cyl_par_myelin_LDFLAGS = $(acme2_cyl_par_LDFLAGS)
acme2_cyl_par_myelin_CXXFLAGS = $(acme2_cyl_par_CXXFLAGS)

acme2_cyl_par_laplace_LDADD = $(acme2_cyl_par_LDADD)
acme2_cyl_par_laplace_LDFLAGS = $(acme2_cyl_par_LDFLAGS)
acme2_cyl_par_laplace_CXXFLAGS = $(acme2_cyl_par_CXXFLAGS)

acme2_cyl_par_mori_LDADD = $(acme2_cyl_par_LDADD)
acme2_cyl_par_mori_LDFLAGS = $(acme2_cyl_par_LDFLAGS)
acme2_cyl_par_mori_CXXFLAGS = $(acme2_cyl_par_CXXFLAGS)

aggregate_fluxes_LDADD = $(acme2_cyl_par_LDADD)
aggregate_fluxes_LDFLAGS = $(acme2_cyl_par_LDFLAGS)
aggregate_fluxes_CXXFLAGS = $(acme2_cyl_par_CXXFLAGS)

eval_grids_LDADD = $(acme2_cyl_par_LDADD)
eval_grids_LDFLAGS = $(acme2_cyl_par_LDFLAGS)
eval_grids_CXXFLAGS = $(acme2_cyl_par_CXXFLAGS)


# Override default CXXFLAGS
CXXFLAGS = -g -O3 -march=native -Wfatal-errors -Wno-tautological-compare #-ftemplate-backtrace-limit=0 #-Wall

# acme2 in plain 2D
acme2_2d_CPPFLAGS = $(AM_CPPFLAGS) $(UG_CPPFLAGS) $(DUNEMPICPPFLAGS) $(HDF5_CPPFLAGS) -UUSE_CYLINDER_COORDS -DAX1_PARALLEL -UMULTIPLE_MEMBRANE_ELEMENTS -DAX1_USE_JACOBIAN_METHODS_IN_OPERATOR #-DNDEBUG

# Commented version with UG flags
#acme2_cyl_par_CPPFLAGS = $(AM_CPPFLAGS) $(UG_CPPFLAGS) $(DUNEMPICPPFLAGS) $(HDF5_CPPFLAGS) -DUSE_CYLINDER_COORDS -DAX1_PARALLEL -UMULTIPLE_MEMBRANE_ELEMENTS -DAX1_USE_JACOBIAN_METHODS_IN_OPERATOR #-DNDEBUG

# Parallel version of acme2_cyl
acme2_cyl_par_CPPFLAGS = $(AM_CPPFLAGS) $(DUNEMPICPPFLAGS) $(HDF5_CPPFLAGS) -DUSE_CYLINDER_COORDS -DAX1_PARALLEL -UMULTIPLE_MEMBRANE_ELEMENTS -DAX1_USE_JACOBIAN_METHODS_IN_OPERATOR -DHAVE_INITIALIZER_LIST #-DNDEBUG

# Multiple membrane elements in acme2_cyl_par
acme2_cyl_par_mme_CPPFLAGS = $(AM_CPPFLAGS) $(DUNEMPICPPFLAGS) $(HDF5_CPPFLAGS) -DUSE_CYLINDER_COORDS -DAX1_PARALLEL -DMULTIPLE_MEMBRANE_ELEMENTS #-DNDEBUG

# Parallel version of acme2_cyl with myelination
acme2_cyl_par_myelin_CPPFLAGS = $(AM_CPPFLAGS) $(DUNEMPICPPFLAGS) $(HDF5_CPPFLAGS) -DUSE_CYLINDER_COORDS -DAX1_PARALLEL -UMULTIPLE_MEMBRANE_ELEMENTS -DAX1_USE_JACOBIAN_METHODS_IN_OPERATOR -DHAVE_INITIALIZER_LIST -DUSE_GMRES_SOLVER #-DNDEBUG

# Use plain Laplace operator instead of PNP in acme2_cyl_par
acme2_cyl_par_laplace_CPPFLAGS = $(AM_CPPFLAGS) $(DUNEMPICPPFLAGS) $(HDF5_CPPFLAGS) -DUSE_CYLINDER_COORDS -DAX1_PARALLEL -DUSE_LAPLACE_OPERATOR #-DNDEBUG

# Use Mori operator-splitting instead of PNP in acme2_cyl_par
acme2_cyl_par_mori_CPPFLAGS = $(AM_CPPFLAGS) $(DUNEMPICPPFLAGS) $(HDF5_CPPFLAGS) -DUSE_CYLINDER_COORDS -DAX1_PARALLEL -DUSE_MORI_OPERATOR_SPLIT #-DNDEBUG

# Setup for multiple grids (grid convergence evaluations)
eval_grids_CPPFLAGS = $(AM_CPPFLAGS) $(DUNEMPICPPFLAGS) $(HDF5_CPPFLAGS) -DUSE_CYLINDER_COORDS -DAX1_PARALLEL -DUSE_EVAL_GRIDS #-DNDEBUG

# Helper program to sum up all transmembrane fluxes into one single file
aggregate_fluxes_CPPFLAGS = $(AM_CPPFLAGS) $(DUNEMPICPPFLAGS) #-DNDEBUG

# don't follow the full GNU-standard
# we need automake 1.5
AUTOMAKE_OPTIONS = foreign 1.5

# pass most important options when "make distcheck" is used
DISTCHECK_CONFIGURE_FLAGS = --with-dune-common=$(DUNE_COMMON_ROOT) --with-dune-istl=$(DUNE_ISTL_ROOT) --with-dune-grid=$(DUNE_GRID_ROOT) --with-dune-localfunctions=$(DUNE_LOCALFUNCTIONS_ROOT) --with-dune-pdelab=$(DUNE_PDELAB_ROOT)  CXX="$(CXX)" CC="$(CC)"

include $(top_srcdir)/am/global-rules
-include $(top_srcdir)/git.mk
