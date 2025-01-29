# Process this file with autoconf to produce a configure script.
AC_INIT(parameters, [m4_esyscmd_s([cat $FROOT/Version])], flapsconsulting@gmail.com)
AC_CONFIG_HEADERS([config.h])
AC_LANG(C++)
AC_CONFIG_SRCDIR([main.cpp])
AM_INIT_AUTOMAKE()
AC_PREFIX_DEFAULT($FROOT)

# Checks for programs...
AC_PROG_CXX
# ... and libraries
AC_CHECK_LIB([blas],scopy_)
AC_CHECK_LIB([lapack],dposv_)

# Check for headers
AC_CHECK_HEADERS([cfenv])

# Check for functions
AC_CHECK_FUNCS([feenableexcept])

# Compile and load with profiling
AC_ARG_ENABLE(profiling,
[ --enable-profiling    Turn on profiling],
[case "${enableval}" in
   yes) profiling=true
	;;
	no)  profiling=false ;;
	*) AC_MSG_ERROR(bad value ${enableval} for --enable-profiling) ;;
esac],[profiling=false])
AM_CONDITIONAL(PROFILING, test x$profiling = xtrue)


AC_CONFIG_FILES([Makefile])
AC_OUTPUT
