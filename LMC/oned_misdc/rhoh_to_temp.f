      subroutine rhoh_to_temp(scal)
      implicit none
      include 'spec.h'
      real*8 scal(-1:nx,nscal)
      real*8 rho, Y(maxspec), hmix, hmixTYP
      integer i,n
      integer cnt,itemp,is

      integer NiterMAX, Niter
      parameter (NiterMAX = 30)
      double precision res(NiterMAX), errMAX, epsHtoTemp
      parameter (epsHtoTemp = 1.e-12)

      integer IWRK
      double precision RWRK

      
      hmixTYP = scal(0,RhoH) / scal(0,Density)
      do i = 1,nx-1
         hmixTYP = MAX(ABS(scal(i,RhoH) / scal(i,Density)),hmixTYP)
      enddo
      errMAX = hmixTYP * epsHtoTemp
 
      do i = 0,nx-1
         rho = 0.d0
         do n=1,Nspec
            is = FirstSpec + n - 1
            rho = rho + scal(i,is)
         enddo
         do n = 1,Nspec
            is = FirstSpec - 1 + n
            Y(n) = scal(i,is) / rho
         enddo
         hmix = scal(i,RhoH) / scal(i,Density)

         call FORT_TfromHYpt(scal(i,Temp),hmix,Y,Nspec,
     &                       errMax,NiterMAX,res,Niter)

         if (Niter.lt.0) then
            print *,'H to T solve failed, Niter,i=',Niter,i
            print *,'T:',scal(i,Temp)
            do n=1,NiterMAX
               print *,'res(',n,')',res(n)
            enddo
            stop
         endif
      enddo

      end
