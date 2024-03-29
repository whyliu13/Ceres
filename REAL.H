
#ifndef BL_REAL_H
#define BL_REAL_H

#ifdef BL_USE_FLOAT
#  undef BL_USE_DOUBLE
#  undef BL_USE_FLOAT
/*
  The macro BL_USE_FLOAT indicates that C++ floating-point calculations should
  use "float" variables and Fortran floating-point calculations should use
  "real*4" variables.  One of BL_USE_FLOAT or BL_USE_DOUBLE must always be
  defined when compiling and using BoxLib.
*/
#  define BL_USE_FLOAT 1
#else
#  undef  BL_USE_FLOAT
#  undef  BL_USE_DOUBLE
/*
  The macro BL_USE_DOUBLE indicates that C++ floating-point calculations
  should use "double" variables and Fortran floating-point calculations
  should use "real*8" variables.  One of BL_USE_FLOAT or BL_USE_DOUBLE must
  always be defined when compiling and using BoxLib.
*/
#  define BL_USE_DOUBLE 1
#endif

#if !defined(BL_LANG_FORT)

/* 
  Real is a typedef specifying the precision of the floating-point
  calculations in C++ code.  It will be either float or double
  depending upon which of the macros BL_USE_FLOAT or
  BL_USE_DOUBLE, respectively, is defined during compilations.  For
  portability, you should write floating-point code in terms of this
  typedef, instead of using float or double directly.

  Note that exactly one of these macros must be defined
  when compiling any module that uses floating-point.
*/
#ifdef BL_USE_FLOAT
    typedef float Real;
#else
    typedef double Real;
#endif

#else

/*
  The REAL_T macro specifies the precision of the floating-point
  calculations in Fortran code.  It will be either real*4 or
  real*8 depending upon which of the symbols BL_USE_FLOAT or
  BL_USE_DOUBLE, respectively, is defined during compilations.  For
  portability, you should write floating-point code in terms of this
  macro, instead of using real*4 or real*8 directly.

  Note that exactly one of these macros must be defined
  when compiling any module that uses floating-point.
*/

#    define INTEGER_T integer

#ifdef BL_USE_FLOAT

#ifdef BL_USE_FORT_STAR_PRECISION
#    define REAL_T real*4
#else
#    define REAL_T REAL
#endif

#ifdef BLANKCONCAT
#        define BL_REAL(a) a/**/E0
#        define BL_REAL_E(a,b) a/**/E/**/b
#else

#if __STDC__==1
#        define BL_REAL(a) a##E0
#        define BL_REAL_E(a,b) a##E##b
#else
#        define BL_REAL(a) a/**/E0
#        define BL_REAL_E(a,b) a/**/E/**/b
#endif

#endif

#else
#ifdef BL_USE_FORT_STAR_PRECISION
#    define REAL_T real*8
#else
#    define REAL_T DOUBLE PRECISION
#endif

#ifdef BLANKCONCAT
#        define BL_REAL(a) a/**/D0
#        define BL_REAL_E(a,b) a/**/D/**/b
#else

#if __STDC__==1
#        define BL_REAL(a) a##D0
#        define BL_REAL_E(a,b) a##D##b
#else
#        define BL_REAL(a) a/**/D0
#        define BL_REAL_E(a,b) a/**/D/**/b
#endif

#endif

#endif

#endif /* !BL_LANG_FORT */

#endif /*BL_REAL_H*/
