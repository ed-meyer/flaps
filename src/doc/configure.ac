# Process this file with autoconf to produce a configure script.
AC_INIT(doc, [m4_esyscmd_s([cat $FROOT/Version])], flapsconsulting@gmail.com)
AM_INIT_AUTOMAKE()
AC_CONFIG_SRCDIR([manual.tex])
AC_CONFIG_HEADERS([config.h])

AC_PREFIX_DEFAULT($FROOT)
AC_CONFIG_SUBDIRS(figures)

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
