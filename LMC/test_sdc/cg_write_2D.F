#include "cg_write_F.H"
#include "ArrayLim.H"

#   if   BL_SPACEDIM==1
#       define  ARLIM(x)  x(1)
#   elif BL_SPACEDIM==2
#       define  ARLIM(x)  x(1),x(2)
#   elif BL_SPACEDIM==3
#       define  ARLIM(x)  x(1),x(2),x(3)
#   endif

#define SDIM 2

      subroutine FORT_CG_WRITE (lo, hi, ncomp, fb, DIMS(fb), fname, ix)

      integer lo(SDIM), hi(SDIM)
      integer ncomp
      integer DIMDEC(fb)
      REAL_T fb(DIMV(fb),*)
      character*(20) fname
      integer ix

      integer j,n


      open(UNIT=11, FILE=fname, STATUS = 'OLD', ACCESS='APPEND')      

      do j=lo(2),hi(2)
         write(11,'(I5,1X)',advance='no')j  
         do n = 1,ncomp
            write(11,'(E22.15,1X)',advance='no') fb(ix,j,n) 
         enddo
         write(11,*)
      enddo

      close(11)

      end
