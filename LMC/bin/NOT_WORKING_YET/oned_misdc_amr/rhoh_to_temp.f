      subroutine rhoh_to_temp(scal,lo,hi)
      implicit none
      include 'spec.h'
      real*8 scal(-2:nfine+1,nscal)
      integer lo,hi

      real*8 rho, Y(Nspec), hmix
      integer i,n
      integer is

      integer NiterMAX, Niter
      parameter (NiterMAX = 30)
      double precision res(NiterMAX), errMAX, epsHtoTemp
      parameter (epsHtoTemp = 1.e-12)

C CEG:: LMC just sets errMAX to 1.d-8
      errMAX = hmix_TYP*1.d-20

      do i=lo,hi
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
