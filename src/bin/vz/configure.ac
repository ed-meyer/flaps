#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_INIT(vz, [m4_esyscmd_s([cat $FROOT/Version])], flapsconsulting@gmail.com)
AC_CONFIG_HEADERS([config.h])
AC_LANG(C++)
AC_CONFIG_SRCDIR([amvz.cpp])
AM_INIT_AUTOMAKE()
AC_PREFIX_DEFAULT($FROOT)

# Checks for X window include files and libraries
# sets X_CFLAGS and X_PRE_LIBS, X_EXTRA_LIBS, and X_LIBS
AC_PATH_XTRA

# Checks for programs.
AC_PROG_CC
AC_PROG_CXX
# ... and libraries
AC_CHECK_LIB([blas],scopy_)
AC_CHECK_LIB([lapack],dposv_)


AC_ARG_WITH(wx_config,
[ --with-wx_config=path path to wx-config ],
wx_config="$withval")

if test -n "$wx_config"; then
	WX_CONFIG="$wx_config"
	AC_MSG_RESULT($wx_config)
else
	WX_CONFIG="wx-config"
	AC_MSG_RESULT(default)
fi

AC_SUBST(WX_CONFIG)

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
