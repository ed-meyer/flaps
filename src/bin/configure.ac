# Process this file with autoconf to produce a configure script.
AC_INIT(bin, [m4_esyscmd_s([cat $FROOT/Version])], flapsconsulting@gmail.com)
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE()
AC_PREFIX_DEFAULT($FROOT)

AC_CONFIG_SUBDIRS(call display etc fem flaps flut matview octlab parameters pz viz)

AC_CONFIG_FILES([Makefile])
AC_OUTPUT

AC_PROG_MAKE_SET
