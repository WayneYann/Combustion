#ifndef _DDBndry_F_H_
#define _DDBndry_F_H_

/*
** $Id: DDBndry_F.H,v 1.1 2011-01-25 22:41:36 marc Exp $
*/

#ifdef BL_LANG_FORT
#    define FORT_DDBNDCCAVG   ddbndccavg
#    define FORT_CRSNCCBND    crsnccbnd
#else
#  if  defined(BL_FORT_USE_UPPERCASE)
#    define FORT_DDBNDCCAVG   DDBNDCCAVG
#    define FORT_CRSNCCBND    CRSNCCBND
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define FORT_DDBNDCCAVG   ddbndccavg
#    define FORT_CRSNCCBND    crsnccbnd
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define FORT_DDBNDCCAVG   ddbndccavg_
#    define FORT_CRSNCCBND    crsnccbnd_
# endif
#include <ArrayLim.H>

extern "C"
{
    void FORT_DDBNDCCAVG(Real*       crse, ARLIM_P(clo), ARLIM_P(chi),
                         const Real* fine, ARLIM_P(flo), ARLIM_P(fhi),
                         const int* lo, const int* hi, const int* nc, const int* ratio);

    void FORT_CRSNCCBND(const int* lo, const int* hi,
                        const Real* fine, ARLIM_P(flo), ARLIM_P(fhi),
                        const Real* crse, ARLIM_P(clo), ARLIM_P(chi),
                        const int* nc, const int* face, const int* ratio);

}
#endif

#endif /*_DDBndry_F_H_*/
