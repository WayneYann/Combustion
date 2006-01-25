#ifndef _DDOp_F_H_
#define _DDOp_F_H_

#ifdef BL_LANG_FORT
#    define FORT_DIFFFORCE       diffforce
#    define FORT_DDC2E           ddc2e
#    define FORT_FLUX            getflux
#    define FORT_INCRDIV         incrdiv
#    define FORT_CPSCALE         cpscale
#    define FORT_THERM           therm
#    define FORT_CRSNCCBND       crsnccbnd
#    define FORT_DDCCAVG         ddccavg
#    define FORT_DDCCINT         ddccint
#else
#  if defined(BL_FORT_USE_UPPERCASE)
#    define FORT_DIFFFORCE       DIFFFORCE
#    define FORT_DDC2E           DDC2E
#    define FORT_FLUX            GETFLUX
#    define FORT_INCRDIV         INCRDIV
#    define FORT_CPSCALE         CPSCALE
#    define FORT_THERM           THERM
#    define FORT_CRSNCCBND       CRSNCCBND
#    define FORT_DDCCAVG         DDCCAVG
#    define FORT_DDCCINT         DDCCINT
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define FORT_DIFFFORCE       diffforce
#    define FORT_DDC2E           ddc2e
#    define FORT_FLUX            getflux
#    define FORT_INCRDIV         incrdiv
#    define FORT_CPSCALE         cpscale
#    define FORT_THERM           therm
#    define FORT_CRSNCCBND       crsnccbnd
#    define FORT_DDCCAVG         ddccavg
#    define FORT_DDCCINT         ddccint
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define FORT_DIFFFORCE       diffforce_
#    define FORT_DDC2E           ddc2e_
#    define FORT_FLUX            getflux_
#    define FORT_INCRDIV         incrdiv_
#    define FORT_CPSCALE         cpscale_
#    define FORT_THERM           therm_
#    define FORT_CRSNCCBND       crsnccbnd_
#    define FORT_DDCCAVG         ddccavg_
#    define FORT_DDCCINT         ddccint_
#  endif

#include <ArrayLim.H>

extern "C" {

    void FORT_DIFFFORCE(const int* lo, const int* hi,
                        Real*       f, ARLIM_P(flo), ARLIM_P(fhi),
                        const Real* X, ARLIM_P(Xlo), ARLIM_P(Xhi),
                        const Real* T, ARLIM_P(Tlo), ARLIM_P(Thi),
                        const Real* dx, const int* dir);

    void FORT_DDC2E(const int* lo, const int* hi,
                    ARLIM_P(clo),ARLIM_P(chi),const Real* c,
                    ARLIM_P(elo),ARLIM_P(ehi),Real* e,
                    const int* nComp, const int* dir);

    void FORT_FLUX(const int* lo, const int* hi,
                   Real*   fluxH,  ARLIM_P(flxHlo), ARLIM_P(flxHhi),
                   Real*   fluxY,  ARLIM_P(flxYlo), ARLIM_P(flxYhi),
                   const Real* de, ARLIM_P(delo),   ARLIM_P(dehi),
                   const Real* Ye, ARLIM_P(Yelo),   ARLIM_P(Yehi),
                   const Real* Te, ARLIM_P(Telo),   ARLIM_P(Tehi),
                   const Real* He, ARLIM_P(Helo),   ARLIM_P(Hehi),
                   const int* i);
    
    void FORT_INCRDIV(const int* lo, const int* hi,
                      Real*     div, ARLIM_P(divlo), ARLIM_P(divhi),
                      const Real* V, ARLIM_P(Vlo),   ARLIM_P(Vhi),
                      const int* dir, const int* nc);

    void FORT_THERM(const int* lo, const int* hi,
                   const Real* lam, ARLIM_P(lamlo),   ARLIM_P(lamhi),
                   const Real* Ye, ARLIM_P(Yelo),   ARLIM_P(Yehi),
                   const Real* Te, ARLIM_P(Telo),   ARLIM_P(Tehi),
                   const int* i);

    void FORT_CPSCALE(const int* lo, const int* hi,
                        Real*       lam, ARLIM_P(lamlo), ARLIM_P(lamhi),
                        const Real* X, ARLIM_P(Xlo), ARLIM_P(Xhi),
                        const Real* T, ARLIM_P(Tlo), ARLIM_P(Thi));

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

