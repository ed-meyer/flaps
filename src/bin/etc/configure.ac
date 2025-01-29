# Process this file with autoconf to produce a configure script.
AC_INIT(etc, [m4_esyscmd_s([cat $FROOT/Version])], flapsconsulting@gmail.com)
AC_CONFIG_HEADERS(config.h)
AC_LANG(C++)
AC_CONFIG_SRCDIR([dtree.cpp])
AM_INIT_AUTOMAKE()
AC_PREFIX_DEFAULT($FROOT)

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
# ... and libraries
AC_CHECK_LIB([blas],scopy_)
AC_CHECK_LIB([lapack],dposv_)

# Checks for header files.
AC_HEADER_DIRENT
AC_CHECK_HEADERS([limits.h stdlib.h string.h unistd.h vnode.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_SIZE_T
AC_TYPE_LONG_DOUBLE
AC_TYPE_LONG_LONG_INT

# Checks for library functions.
AC_FUNC_LSTAT

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
