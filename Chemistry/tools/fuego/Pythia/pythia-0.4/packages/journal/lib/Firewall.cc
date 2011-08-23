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

#include <portinfo>
#include <map>
#include <string>
#include <sstream>

#include "Category.h"
#include "Diagnostic.h"
#include "Firewall.h"
#include "Index.h"
#include "LocalIndex.h"

using namespace journal;

// statics
static Firewall::index_t * _index = new LocalIndex(true);

// interface
void Firewall::newIndex(index_t * newIndex) {
    delete _index;
    _index = newIndex;
    return;
}

// meta-methods
Firewall::~Firewall() {}

Firewall::Firewall(string_t name):
    Diagnostic("firewall", name, _index->category(name))
{}

// version
// $Id$

// End of file
