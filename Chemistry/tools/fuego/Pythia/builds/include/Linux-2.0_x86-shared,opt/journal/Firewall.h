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

#if !defined(journal_Firewall_h)
#define journal_Firewall_h

// forward declarations
namespace journal {

    class Firewall;
    class Index;
    class Diagnostic;

}

// 
class journal::Firewall : public journal::Diagnostic {
// types
public:
    typedef Index index_t;

// interface
public:
    static void newIndex(index_t *);

    bool affirm(bool expression) { return !(expression); }

// meta-methods
public:
    Firewall(string_t name);
    virtual ~Firewall();

// disable these
private:
    Firewall(const Firewall &);
    const Firewall & operator=(const Firewall &);
};

#endif

// version
// $Id$

// End of file
