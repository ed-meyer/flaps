# Process this file with autoconf to produce a configure script.
AC_INIT(fem, [m4_esyscmd_s([cat $FROOT/Version])], flapsconsulting@gmail.com)
AC_CONFIG_HEADERS([config.h])
AC_LANG(C++)
AM_INIT_AUTOMAKE()
AC_PREFIX_DEFAULT($FROOT)

# Checks for programs.
AC_PROG_CC
AC_PROG_CXX
# ... and libraries
AC_CHECK_LIB([blas],scopy_)
AC_CHECK_LIB([lapack],dposv_)

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
