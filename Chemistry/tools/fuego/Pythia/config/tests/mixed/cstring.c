/*
//
// -*- C -*-
//--------------------------------------------------------------------------
//
//
//                       Michael A.G. Aivazis
//                California Institute of Technology
//                   (C) 2000 All Rights Reserved
//
//--------------------------------------------------------------------------

// $Log: cstring.c,v $
// Revision 1.4  2000/08/14 23:51:16  cummings
// Added F77 to C symbol translation for the cases of all uppercase letters
// with or without a trailing underbar.
//
// Revision 1.3  2000/08/14 22:53:12  aivazis
// Populated the symbol translation section a bit more
//
// Revision 1.2  2000/08/10 00:55:49  aivazis
// Added the translation section
// Improved the trailing info dump
//
// Revision 1.1  2000/08/09 22:57:15  aivazis
// Original source
//
 
*/
#include <portinfo.h>
#include <stdio.h>
#include <string.h>

/* Symbol translation */
#if defined(NEEDS_F77_TRANSLATION)

#if defined(F77EXTERNS_LOWERCASE_TRAILINGBAR)
/*
 * This is the default naming convention for us.
 * No translation necessary
 */
#elif defined(F77EXTERNS_LOWERCASE_NOTRAILINGBAR)
#define cget_string_ cget_string
#define fsend_string_ fsend_string

#elif defined(F77EXTERNS_UPPERCASE_TRAILINGBAR)
#define cget_string_ CGET_STRING_
#define fsend_string_ FSEND_STRING_

#elif defined(F77EXTERNS_UPPERCASE_NOTRAILINGBAR)
#define cget_string_ CGET_STRING
#define fsend_string_ FSEND_STRING

#else
#error Unknown translation for FORTRAN external symbols
#endif

#endif

/* Prototype for the FORTRAN string */
void fsend_string_();

void cget_string_(void * s)
{
    char buffer[8];
    int int_bytes;
    
    printf("Got it!\n");
    
/* The string is supposed to be "FORTRAN" */
/* Try to see whether it starts at the passed address */
    
    memcpy(buffer, s, 8);
    buffer[7] = 0;
    printf("The string is: %s\n", buffer);
    
/*
 * print the first four bytes at the string address as an int to
 * see whether it's the string length
 */
    
    memcpy(&int_bytes, s, sizeof(int)); 
    printf("The first four bytes are: %d\n", int_bytes);
    
/*
 * print the first four bytes after the string as an int to
 * see whether it's the string length
 */
    
    printf("sizeof(\"FORTRAN\") = %d\n", (int)sizeof("FORTRAN"));
    memcpy(&int_bytes, (char *)s + sizeof("FORTRAN"), sizeof(int));
    printf("The trailing four bytes are: %d\n", int_bytes);
    
    return;
}

int main()
{
    fsend_string_();
    
    return 0;
}

/* End of file */
