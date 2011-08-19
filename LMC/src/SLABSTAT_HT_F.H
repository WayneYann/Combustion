#ifndef _SLABSTAT_HT_F_H_
#define _SLABSTAT_HT_F_H_

#ifdef BL_LANG_FORT
#    define FORT_HT_BASICSTATS_NCTRAC    ht_basicstats_nctrac
#    define FORT_HT_BASICSTATS_CTRAC     ht_basicstats_ctrac
#else
#  ifdef BL_FORT_USE_UPPERCASE
#    define FORT_HT_BASICSTATS_NCTRAC    HT_BASICSTATS_NCTRAC
#    define FORT_HT_BASICSTATS_CTRAC     HT_BASICSTATS_CTRAC
#  else
#  ifdef BL_FORT_USE_LOWERCASE
#    define FORT_HT_BASICSTATS_NCTRAC    ht_basicstats_nctrac
#    define FORT_HT_BASICSTATS_CTRAC     ht_basicstats_ctrac
#  else
#    define FORT_HT_BASICSTATS_NCTRAC    ht_basicstats_nctrac_
#    define FORT_HT_BASICSTATS_CTRAC     ht_basicstats_ctrac_
#  endif
#  endif

extern "C"
{
     void FORT_HT_BASICSTATS_NCTRAC (const Real* src, 
                                       ARLIM_P(src_lo), ARLIM_P(src_hi),
                                     const int* nsrc,
                                     Real* dst, 
                                       ARLIM_P(dst_lo), ARLIM_P(dst_hi),
                                     const int* ndst,
                                     const Real* dt,
                                     const Real* dx);

     void FORT_HT_BASICSTATS_CTRAC (const Real* src, 
                                      ARLIM_P(src_lo), ARLIM_P(src_hi),
                                    const int* nsrc,
                                    Real* dst, 
                                      ARLIM_P(dst_lo), ARLIM_P(dst_hi),
                                    const int* ndst,
                                    const Real* dt,
                                    const Real* dx);
}

#endif

#endif /*_SLABSTAT_HT_F_H_*/
