// -*- C++ -*-
//
//--------------------------------------------------------------------------------
//
//                              Michael A.G. Aivazis
//                       California Institute of Technology
//                       (C) 1998-2003  All Rights Reserved
//
// <LicenseText>
//
//--------------------------------------------------------------------------------
//

#if !defined(journal_journal_h)
#define journal_journal_h

#include <string>
#include <sstream>

#include "Category.h"
#include "Diagnostic.h"
#include "Firewall.h"
#include "Debug.h"
#include "Info.h"
#include "Warning.h"
#include "Error.h"

namespace journal {

    typedef Firewall firewall_t;
    typedef Debug debug_t;
    typedef Info info_t;
    typedef Warning warning_t;
    typedef Error error_t;

}

#endif

// version
// $Id$

// End of file
