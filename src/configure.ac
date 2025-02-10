# Process this file with autoconf to produce a configure script.
AC_INIT(src, [m4_esyscmd_s([cat $FROOT/Version])], flapsconsulting@gmail.com)
AM_INIT_AUTOMAKE(no-dependencies)
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_SUBDIRS([scripts] [lib] [bin])
AC_PREFIX_DEFAULT($FROOT)

AC_CONFIG_FILES([Makefile])
AC_OUTPUT

AC_PROG_MAKE_SET
