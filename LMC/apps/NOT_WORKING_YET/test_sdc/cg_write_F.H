#ifndef _cg_write_F_H_
#define _cg_write_F_H_

#ifdef BL_LANG_FORT
#    define FORT_CG_WRITE       dcg_write
#else
#  if defined(BL_FORT_USE_UPPERCASE)
#    define FORT_CG_WRITE       DCG_WRITE
#  elif defined(BL_FORT_USE_LOWERCASE)
#    define FORT_CG_WRITE       dcg_write
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define FORT_CG_WRITE       dcg_write_
# endif

#include <ArrayLim.H>

extern "C" {
  void FORT_CG_WRITE(const int* lo, const int* hi, const int* ncomp,
		     const Real* fb, ARLIM_P(fblo), ARLIM_P(fbhi),
		     const char* fname, const int* ix);
}
#endif
#endif
