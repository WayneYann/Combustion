c
c $Id: DDBndry_3D.F,v 1.1 2011-01-25 22:41:36 marc Exp $
c
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "BC_TYPES.H"
#include "REAL.H"
#include "CONSTANTS.H"
#include "DDBndry_F.H"
#include "ArrayLim.H"

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
c
      integer i,j,k,n,ii,jj,kk
      REAL_T fac

      fac = 1.d0/(ratio(1)*ratio(2)*ratio(3))
c     NOTE: generalize to rz
      do n = 1, nc
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  c(i,j,k,n) = 0.d0
                  do kk = 1,ratio(3)
                     do jj = 1,ratio(2)
                         do ii = 1,ratio(1)
                            c(i,j,k,n) = c(i,j,k,n) + f(2*i+ii-1,2*j+jj-1,2*k+kk-1,n)
                         enddo
                     enddo
                  enddo
                  c(i,j,k,n) = c(i,j,k,n)*fac
               enddo
            enddo
         enddo
      enddo
      end

