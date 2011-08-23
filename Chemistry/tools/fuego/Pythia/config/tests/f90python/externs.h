/*
// -*- C++ -*-
//
//--------------------------------------------------------------------------------
//
//                              Michael A.G. Aivazis
//                       California Institute of Technology
//                          (C) 2000 All Rights Reserved
//
// <LicenseText>
//
//--------------------------------------------------------------------------------
//
// $Log: externs.h,v $
// Revision 1.2  2000/09/12 03:43:48  aivazis
// Added F77EXTERNS_LOWERCASE_TRAILINGBAR case
//
// Revision 1.1  2000/09/11 02:15:13  aivazis
// Original source
//
*/

#if !defined(externs_h)
#define externs_h

#if defined(F77EXTERNS_LOWERCASE_TRAILINGBAR)

/* This is the default. No translation */

#elif defined(F77EXTERNS_NOTRAILINGBAR)

#define hello_ hello

#elif defined(F77EXTERNS_EXTRATRAILINGBAR)

#define hello_ hello__

#elif defined(F77EXTERNS_UPPERCASE_NOTRAILINGBAR)

#define hello_ HELLO

#elif defined(F77EXTERNS_COMPAQ_F90)

/*
    This symbol does not contain "_" in its name and therefore requires no translation
     hello_
*/

#else
#error Uknown translation for FORTRAN external symbols
#endif

/* The declarations */
#if defined(__cplusplus)
extern "C" {
#endif

void hello_();

#if defined(__cplusplus)
}
#endif

#endif /* EXTERNS_H */

/* End of file */

