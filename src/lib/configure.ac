# Process this file with autoconf to produce a configure script.
AC_INIT([libflaps],[m4_esyscmd_s([cat $FROOT/Version])],[flapsconsulting@gmail.com])
AC_CONFIG_HEADERS([config.h])
# AC_LANG sets the compiler to use for subsequent tests
AC_LANG(C++)
LT_INIT
AM_INIT_AUTOMAKE()
AC_PREFIX_DEFAULT($FROOT)
AC_PROG_CC
AC_PROG_CXX
AC_CHECK_LIB([blas],scopy_)
AC_CHECK_LIB([lapack],dposv_)

# recommended by libtoolize:
AC_CONFIG_MACRO_DIRS([m4])

# AC_ENABLE_SHARED
# AC_DISABLE_STATIC
LT_INIT([libtool])

AC_C_BIGENDIAN

AC_TYPE_LONG_DOUBLE

# Checks for header files.
AC_CHECK_HEADERS([ftw.h numbers boost/core/demangle.hpp fftw3.h])

AC_CHECK_DECL(nftw, includes=ftw.h)

# Checks for library functions.
# AC_CHECK_FUNCS([hstrerror ftw nftw isfdtype readline])
AC_CHECK_FUNCS([ftw nftw])
AC_CHECK_LIB(fftw3,fftw_execute)

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

# Compile with -fanalyze
AC_ARG_ENABLE(analyze,
[ --enable-analyze    Turn on analyze],
[case "${enableval}" in
   yes) analyze=true
	;;
	no)  analyze=false ;;
	*) AC_MSG_ERROR(bad value ${enableval} for --enable-analyze) ;;
esac],[analyze=false])
AM_CONDITIONAL(ANALYZE, test x$analyze = xtrue)

AC_CONFIG_FILES([Makefile])

AC_OUTPUT
