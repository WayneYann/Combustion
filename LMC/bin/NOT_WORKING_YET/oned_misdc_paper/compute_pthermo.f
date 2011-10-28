      subroutine compute_pthermo(scal,pthermo)
      include 'spec.h'
      real*8  scal(-1:nx,nscal)
      real*8 pthermo(-1:nx)
      
      real*8 Y(maxspec), RWRK
      integer i,n
      integer ispec, IWRK

c     Define thermodynamic pressure
      do i = 0,nx-1
         do n = 1,nspec
            ispec = FirstSpec-1+n
            Y(n) = scal(i,ispec) / scal(i,Density)
         end do
         CALL CKPY(scal(i,Density),scal(i,Temp),Y,IWRK,RWRK,pthermo(i))
         if (i.eq.0) then
            pthermo_min = pthermo(i)
            pthermo_max = pthermo(i)
         else
            pthermo_min = MIN(pthermo_min,pthermo(i))
            pthermo_max = MAX(pthermo_max,pthermo(i))
         endif
      end do
      pthermo(-1) = pthermo(0)
      pthermo(nx) = pthermo(nx-1)
      
      print *,'PTHERMO MIN/MAX ',pthermo_min,pthermo_max
      
      end
