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
// Revision 1.2  2002/01/28 23:43:40  aivazis
// Defined the external symbol translation arlgorithm
//
// Revision 1.1  2001/07/11 02:07:26  cummings
// Build procedure files for new gcc 3.0 compiler, borrowed from gcc-2.95.2.  I have removed the *using namespace std* hack from the config_compiler.h file here.
//
//
*/

#if !defined(__config_f77_h__)
#define __config_f77_h__

#define NEEDS_F77_TRANSLATION
#define F77EXTERNS_LOWERCASE_TRAILINGBAR

#endif

/* End of file */
