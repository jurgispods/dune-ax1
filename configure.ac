# -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ(2.50)
DUNE_AC_INIT # gets module version from dune.module file
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src/example03.cc])
AM_CONFIG_HEADER([config.h])


# we need no more than the standard DE-stuff
# this module depends on dune-common dune-istl dune-grid dune-localfunctions dune-pdelab
# this implies checking for [dune-common], [dune-istl], [dune-grid], [dune-localfunctions], [dune-pdelab]
DUNE_CHECK_ALL

# implicitly set the Dune-flags everywhere
AC_SUBST(AM_CPPFLAGS, $DUNE_CPPFLAGS)
AC_SUBST(AM_LDFLAGS, $DUNE_LDFLAGS)
LIBS="$DUNE_LIBS"

AC_CONFIG_FILES([
  Makefile
  src/Makefile
  test/Makefile
  doc/Makefile
  doc/doxygen/Makefile
  doc/doxygen/Doxyfile
  dune/Makefile
  dune/ax1/Makefile
  dune/ax1/acme0/Makefile
  dune/ax1/acme0/common/Makefile
  dune/ax1/acme0/configurations/Makefile
  dune/ax1/acme0/configurations/default/Makefile
  dune/ax1/acme0/configurations/hamburger/Makefile
  dune/ax1/acme0/fully-implicit/Makefile
  dune/ax1/acme0/operator-split/Makefile
  dune/ax1/acme1/Makefile
  dune/ax1/acme1/common/Makefile
  dune/ax1/acme1/configurations/Makefile
  dune/ax1/acme1/configurations/default/Makefile
  dune/ax1/acme1/configurations/hamburger/Makefile
  dune/ax1/acme1/operator-split/Makefile
  dune/ax1/acme1MD/Makefile
  dune/ax1/acme1MD/common/Makefile
  dune/ax1/acme1MD/configurations/Makefile
  dune/ax1/acme1MD/configurations/bigmac/Makefile
  dune/ax1/acme1MD/configurations/cheeseburger/Makefile
  dune/ax1/acme1MD/configurations/default/Makefile
  dune/ax1/acme1MD/configurations/hamburger/Makefile
  dune/ax1/acme1MD/configurations/step/Makefile
  dune/ax1/acme1MD/fully-implicit/Makefile
  dune/ax1/acme1MD/operator-split/Makefile
  dune/ax1/acme2/Makefile
  dune/ax1/acme2/common/Makefile
  dune/ax1/acme2/configurations/Makefile
  dune/ax1/acme2/configurations/default/Makefile
  dune/ax1/acme2/configurations/step/Makefile
  dune/ax1/acme2/configurations/test_scales/Makefile
  dune/ax1/acme2/operator/Makefile
  dune/ax1/acme2_cyl/Makefile
  dune/ax1/acme2_cyl/common/Makefile
  dune/ax1/acme2_cyl/configurations/Makefile
  dune/ax1/acme2_cyl/configurations/default/Makefile
  dune/ax1/acme2_cyl/configurations/step/Makefile
  dune/ax1/acme2_cyl/configurations/test_scales/Makefile
  dune/ax1/acme2_cyl/operator/Makefile
  dune/ax1/channels/Makefile
  dune/ax1/common/Makefile
  m4/Makefile
  dune-ax1.pc
])
AC_OUTPUT
# finally print the summary information
DUNE_SUMMARY_ALL
