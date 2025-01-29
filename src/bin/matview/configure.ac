#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_INIT(matview, [m4_esyscmd_s([cat $FROOT/Version])], flapsconsulting@gmail.com)
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_SRCDIR([matview.c])
AM_INIT_AUTOMAKE()
AC_PREFIX_DEFAULT($FROOT)

# Checks for programs.
AC_PROG_CC
AC_PROG_CXX

# Checks for X window include files and libraries
# sets X_CFLAGS and X_PRE_LIBS, X_EXTRA_LIBS, and X_LIBS
AC_PATH_XTRA

# Checks for typedefs, structures, and compiler characteristics.
AC_LANG(C++)

# Checks for library functions.
AC_CHECK_FUNCS([feenableexcept])
# Check for header files
AC_CHECK_HEADERS([tk.h tcl.h])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
