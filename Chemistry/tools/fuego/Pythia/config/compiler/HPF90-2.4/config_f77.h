/*
// -*- C -*-
//--------------------------------------------------------------------------
//
//
//                       Michael A.G. Aivazis
//                California Institute of Technology
//                   (C) 1999 All Rights Reserved
//
//--------------------------------------------------------------------------

// $Log: config_f77.h,v $
// Revision 1.3  2000/09/15 01:18:59  cummings
// Changed translation type for Fortran symbols to
// F77EXTERNS_LOWERCASE_TRAILINGBAR, since the HP
// F90 compiler does indeed append an underbar to
// its symbols.  (The online instructions given by
// "f90 +usage" claim that no trailing underbar is
// the default, but I checked that this is not true!
// At least, not in 64-bit mode...)
//
// Revision 1.2  2000/08/15 16:49:44  rapa
// fixed bug
//
// Revision 1.1  2000/08/14 23:05:21  aivazis
// Original source
//
// Revision 1.1  2000/04/27 19:06:09  aivazis
// Mixed language improvements
//
*/

#if !defined(__config_f77_h__)
#define __config_f77_h__

#define NEEDS_F77_TRANSLATION
#define F77EXTERNS_LOWERCASE_TRAILINGBAR

#endif

/* End of file */
