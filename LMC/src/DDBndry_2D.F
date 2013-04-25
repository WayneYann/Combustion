
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include <BC_TYPES.H>
#include <REAL.H>
#include <CONSTANTS.H>
#include <DDBndry_F.H>
#include <ArrayLim.H>

#define SDIM BL_SPACEDIM
      
      subroutine FORT_DDBNDCCAVG(c, DIMS(c), f, DIMS(f), lo, hi, nc, ratio)
      integer nc
      integer DIMDEC(f)
      integer DIMDEC(c)
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      REAL_T f(DIMV(f),nc)
      REAL_T c(DIMV(c),nc)
      integer ratio(SDIM)

      integer i,j,n,ii,jj
      REAL_T fac

      fac = 1.d0 / (ratio(1)*ratio(2))

c     NOTE: generalize to rz
      do n = 1, nc
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               c(i,j,n) = 0.d0
               do jj = 1,ratio(2)
                  do ii = 1,ratio(1)
                     c(i,j,n) = c(i,j,n) + f(2*i+ii-1,2*j+jj-1,n)
                  enddo
               enddo
               c(i,j,n) = c(i,j,n)*fac
            enddo
         enddo
      enddo
      end
      
