# Process this file with autoconf to produce a configure script.
AC_INIT(flaps, [m4_esyscmd_s([cat $FROOT/Version])], flapsconsulting@gmail.com)
AM_INIT_AUTOMAKE()
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_SUBDIRS([src] [demo])
AC_PREFIX_DEFAULT($FROOT)


AC_ARG_ENABLE(matview, [ --enable-matview   build matview, a program for viewing matrices (requires tk and tcl)],
	[case "${enableval}" in
		yes) matview=true;;
		no) matview=false;;
		*) AC_MSG_ERROR(bad value ${enableval} for --enable-matview) ;;
		esac],[matview=false])
AM_CONDITIONAL(MATVIEW, test x$matview = xtrue)

AC_CONFIG_FILES([Makefile])
AC_OUTPUT

AC_PROG_MAKE_SET
