# pcalg configure.ac by Alain Hauser <alain.hauser@biology.unibe.ch>
#
# Based on the example provided in "Writing R extensions", and the Autoconf
# script of the R package "RQuantLib"
AC_INIT([pcalg], 2.0.0) # package name, version
AC_CONFIG_SRCDIR([DESCRIPTION])

# We are using C++
AC_LANG(C++)
AC_REQUIRE_CPP

# Determine compile flags from R
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
    AC_MSG_ERROR([Could not determine R_HOME.])   
fi
CXX=`${R_HOME}/bin/R CMD config CXX`
CXXFLAGS=`"${R_HOME}/bin/R" CMD config CXXFLAGS`
LDFLAGS=`"${R_HOME}/bin/R" CMD config LDFLAGS`

# Look for Rscript, but use the one found via R_HOME to allow for multiple installations
# Block taken from RQuantLib's configure.in
AC_DEFUN(AC_PROG_RSCRIPT, [AC_CHECK_PROG(RSCRIPT,Rscript,yes)])
AC_PROG_RSCRIPT
#AC_MSG_CHECKING([for Rscript])
if test x"${RSCRIPT}" == x"yes" ; then
#     AC_MSG_RESULT([yes])    
     ## Rcpp compile flag providing header directory containing Rcpp.h
     rcpp_cxxflags=`${R_HOME}/bin/Rscript -e 'Rcpp:::CxxFlags()'`
     ## link flag providing libary as well as path to library, and optionally rpath
     rcpp_ldflags=`${R_HOME}/bin/Rscript -e 'Rcpp:::LdFlags()'`
     ## now use all these
     AC_SUBST([CXXFLAGS],["${CXXFLAGS} $rcpp_cxxflags"])
     AC_SUBST([LDFLAGS],["${LDFLAGS} $rcpp_ldflags"])
else
#    AC_MSG_RESULT([no])    
    echo "
  Your installation does not appear to have Rscript installed.

  Please make sure that you have a working and complete R installation.
"
    exit 1
fi

# Check for boost, this segment is from QuantLib's configure.ac
AC_ARG_WITH([boost-include],
            AC_HELP_STRING([--with-boost-include=INCLUDE_PATH],
                           [Supply the location of Boost header files]),
            [boost_include_path=$withval],
            [boost_include_path=""])
if test [ -n "$boost_include_path" ] ; then
   AC_SUBST([CXXFLAGS],["${CXXFLAGS} -I${boost_include_path}"])
fi

# Check for Boost headers
AC_DEFUN([CHECK_BOOST_DEVEL],
[AC_MSG_CHECKING([for Boost development files])
 AC_TRY_COMPILE(
    [@%:@include <boost/version.hpp>
     @%:@include <boost/lambda/lambda.hpp>
     @%:@include <boost/lambda/bind.hpp>
     @%:@include <boost/graph/adjacency_list.hpp>
     @%:@include <boost/dynamic_bitset.hpp>],
    [],
    [AC_MSG_RESULT([yes])],
    [AC_MSG_RESULT([no])
     AC_MSG_ERROR([Boost development files not found])
    ])
])

# Check for Boost version
AC_DEFUN([CHECK_BOOST_VERSION],
[AC_MSG_CHECKING([for minimal Boost version])
 AC_REQUIRE([CHECK_BOOST_DEVEL])
 AC_TRY_COMPILE(
    [@%:@include <boost/version.hpp>],
    [@%:@if BOOST_VERSION < 103400
     @%:@error too old
     @%:@endif],
    [AC_MSG_RESULT([yes])],
    [AC_MSG_RESULT([no])
     AC_MSG_ERROR([outdated Boost installation])
    ])
])

AC_DEFUN([CHECK_BOOST],
[AC_REQUIRE([CHECK_BOOST_DEVEL])
 AC_REQUIRE([CHECK_BOOST_VERSION])
])
CHECK_BOOST

AC_SUBST(CXXFLAGS)
AC_SUBST(LDFLAGS)
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
echo "Completed configuration and ready to build."
