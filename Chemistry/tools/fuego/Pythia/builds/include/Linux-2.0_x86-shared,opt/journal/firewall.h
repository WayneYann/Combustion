/*
 * -*- C++ -*-
 * 
 *  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * 
 *                               Julian C. Cummings
 *                       California Institute of Technology
 *                        (C) 1998-2003 All Rights Reserved
 * 
 *  <LicenseText>
 * 
 *  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * 
 */


#if !defined(journal_firewall_h)
#define journal_firewall_h


/* get definition of __HERE__ macros */
#include "macros.h"


#ifdef __cplusplus
extern "C"
#endif
void firewall_hit(__HERE_DECL__, const char * fmt, ...);

#ifdef __cplusplus
extern "C"
#endif
void firewall_affirm(int condition, __HERE_DECL__, const char * fmt, ...);

#endif


/* version
 * $Id$
 */

/* End of file */

