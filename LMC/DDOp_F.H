#ifndef _DDOp_F_H_
#define _DDOp_F_H_

#ifdef BL_LANG_FORT
#    define FORT_DDFORCE         ddforce
#    define FORT_DDC2E           ddc2e
#    define FORT_DDFLUX          ddflux
#    define FORT_DDETC           ddetc
#    define FORT_CRSNCCBND       crsnccbnd
#    define FORT_DDCCAVG         ddccavg
#    define FORT_DDCCINT         ddccint
#else
#  if defined(BL_FORT_USE_UPPERCASE)
#    define FORT_DDFORCE         DDFORCE
#    define FORT_DDC2E           DDC2E
#    define FORT_DDFLUX          DDFLUX
#    define FORT_DDETC           DDETC
#    define FORT_CRSNCCBND       CRSNCCBND
#    define FORT_DDCCAVG         DDCCAVG
#    define FORT_DDCCINT         DDCCINT
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define FORT_DDFORCE         ddforce
#    define FORT_DDC2E           ddc2e
#    define FORT_DDFLUX          ddflux
#    define FORT_DDETC           ddetc
#    define FORT_CRSNCCBND       crsnccbnd
#    define FORT_DDCCAVG         ddccavg
#    define FORT_DDCCINT         ddccint
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define FORT_DDFORCE         ddforce_
#    define FORT_DDC2E           ddc2e_
#    define FORT_DDFLUX          ddflux_
#    define FORT_DDETC           ddetc_
#    define FORT_CRSNCCBND       crsnccbnd_
#    define FORT_DDCCAVG         ddccavg_
#    define FORT_DDCCINT         ddccint_
#  endif

#include <ArrayLim.H>

extern "C" {

    void FORT_DDFORCE(const int* lo, const int* hi,
                      Real*       f, ARLIM_P(flo), ARLIM_P(fhi),
                      const Real* X, ARLIM_P(Xlo), ARLIM_P(Xhi),
                      const Real* T, ARLIM_P(Tlo), ARLIM_P(Thi),
                      const Real* dx, const int* dir);

    void FORT_DDC2E(const int* lo, const int* hi,
                    const Real* c,ARLIM_P(clo),ARLIM_P(chi),
                    Real* e,ARLIM_P(elo),ARLIM_P(ehi),
                    const int* nComp, const int* dir);

    void FORT_DDFLUX(const int* lo, const int* hi, const Real* dx, const int* dir,
                     Real*   fluxYH,  ARLIM_P(fluxYHlo), ARLIM_P(fluxYHhi),
                     const Real* YTc, ARLIM_P(YTclo),    ARLIM_P(TYchi),
                     const Real* cpe,  ARLIM_P(ceplo),     ARLIM_P(cpehi),
                     const int* for_T_or_H, const int* fillA,
                     const Real* a,   ARLIM_P(alo),      ARLIM_P(ahi));

    void FORT_DDETC(const int* lo, const int* hi,
                    Real*     div, ARLIM_P(divlo), ARLIM_P(divhi),
                    const Real* V, ARLIM_P(Vlo),   ARLIM_P(Vhi),
                    const Real* a, const int* dir, const int* nc, const int* diff0_avg1);

    void FORT_CRSNCCBND(const int* lo, const int* hi,
                        const Real* fine, ARLIM_P(flo), ARLIM_P(fhi),
                        const Real* crse, ARLIM_P(clo), ARLIM_P(chi),
                        const int* nc, const int* face, const int* ratio);

    void FORT_DDCCAVG(Real*       crse, ARLIM_P(clo), ARLIM_P(chi),
                      const Real* fine, ARLIM_P(flo), ARLIM_P(fhi),
                      const int* lo, const int* hi, const int* nc, const int* ratio);

    void FORT_DDCCINT(Real*       fine, ARLIM_P(flo), ARLIM_P(fhi),
                      const Real* crse, ARLIM_P(clo), ARLIM_P(chi),
                      const int* lo, const int* hi, const int* nc, const int* ratio);

}

#endif
#endif

