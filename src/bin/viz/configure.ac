#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_INIT(viz, [m4_esyscmd_s([cat $FROOT/Version])], flapsconsulting@gmail.com)
AC_CONFIG_HEADERS([config.h])
AC_LANG(C++)
AC_CONFIG_SRCDIR([viz.cpp])
AM_INIT_AUTOMAKE()
AC_PREFIX_DEFAULT($FROOT)

# Checks for X window include files and libraries
# sets X_CFLAGS and X_PRE_LIBS, X_EXTRA_LIBS, and X_LIBS
AC_PATH_XTRA

# Checks for programs...
AC_PROG_CC
AC_PROG_CXX
# ... and libraries
AC_CHECK_LIB([blas],scopy_)
AC_CHECK_LIB([lapack],dposv_)

# check for Posix message queue
AC_CHECK_HEADERS([mqueue.h])
# AC_CHECK_FUNC([mq_send],
# 	[AC_DEFINE([HAVE_MQ_SEND], [1], [Define to 1 if you have the mq_send function])],
# 	[AC_MSG_ERROR([mq_send is required but was not found])])

# check for wxEntry
# AC_CHECK_HEADERS([wx/wx.h])
# AC_CHECK_FUNC([wxEntry(argc,argv)],
# 	[AC_DEFINE([HAVE_WXENTRY], [1], [Define to 1 if you have the wxEntry function])],
# 	[AC_MSG_ERROR([wxEntry is required but was not found])])

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
