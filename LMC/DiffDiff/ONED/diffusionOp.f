      subroutine neg_divF_H(LofS,S,PTCec,rhoTDec,rhoDijec,dx)
      implicit none
      include 'spec.h'
      real*8 LofS(maxspec+1,1:nx), S(maxscal,0:nx+1)
      real*8 PTCec(1:nx+1)
      real*8 rhoTDec(maxspec,1:nx+1),rhoDijec(maxspec,maxspec,1:nx+1),dx
      real*8 coef(maxspec+1,1:nx)
      real*8 Y(maxspec,0:nx+1), X(maxspec,0:nx+1), WWe, PTC
      real*8 de(maxspec+1,1:nx+1), q(1:nx+1), F(maxspec,1:nx+1)
      real*8 Ye(maxspec), Te, dxInv, He(maxspec)
      real*8 rhoe, enthe
      integer i,n,m,Niter,maxIter
      real*8 res(NiterMAX)

      dxInv = 1.d0/dx
      do i=0,nx+1
         do n=1,Nspec
            Y(n,i) = S(FirstSpec+n-1,i) / S(Density,i) 
         enddo
         call CKYTX(Y(1,i),IWRK,RWRK,X(1,i))
      enddo

      maxiter=0
      do i=1,nx+1
         do n=1,Nspec
            de(n,i) = (X(n,i)-X(n,i-1)) / dx
         enddo
         de(Nspec+1,i) = (S(Temp,i)-S(Temp,i-1)) / dx

         do n=1,Nspec
            Ye(n) = 0.5d0*(Y(n,i)+Y(n,i-1))
         enddo
         call CKMMWY(Ye,IWRK,RWRK,WWe)

         Te = 0.5d0*(S(Temp,i)+S(Temp,i-1))
         rhoe = 0.5d0*(S(Density,i)+S(Density,i-1))
         if (setTfromH.eq.1) then
            enthe = 0.5d0*(S(RhoH,i)+S(RhoH,i-1))/rhoe
            call FORT_TfromHYpt(Te,enthe,Ye,Nspec,errMax,NiterMAX,res,Niter)
            if (Niter.lt.0) then
               print *,'RhoH->T failed at i=',i
               stop
            endif
            maxiter=MAX(maxiter,Niter)
         else if (setTfromH.eq.0) then
            Te = Pcgs * WWe / (rhoe * RU)
         endif
         call CKHMS(Te,IWRK,RWRK,He) 

         q(i) = 0.d0
         do n = 1,Nspec
            F(n,i) = - Ye(n)*rhoTDec(n,i)*de(Nspec+1,i)/Te
            do m = 1,Nspec
               F(n,i) = F(n,i) - rhoDijec(n,m,i)*de(m,i)
            enddo
            q(i) = q(i) - (RU*Te/WWe)*rhoTDec(n,i)*de(n,i) + He(n)*F(n,i)
         enddo
         q(i) = q(i) - PTCec(i)*de(Nspec+1,i)
      enddo

      do i=1,nx
         do n = 1,Nspec
            LofS(n,i) = - dxInv*(F(n,i+1) - F(n,i))
         enddo
         LofS(Nspec+1,i) = - dxInv*(q(i+1) - q(i))
      enddo

      end

      subroutine neg_divF_H_Approx(LofS,S,PTCec,rhoDiec,dx)
      implicit none
      include 'spec.h'
      real*8 LofS(maxspec+1,1:nx), S(maxscal,0:nx+1)
      real*8 PTCec(1:nx+1)
      real*8 rhoDiec(maxspec,1:nx+1),cpicc(maxspec,0:nx+1),dx
      real*8 Y(maxspec,0:nx+1), X(maxspec,0:nx+1), WWe, PTC
      real*8 de(maxspec+1,1:nx+1), q(1:nx+1), F(maxspec,1:nx+1)
      real*8 Ye(maxspec), Te, dxInv2, He(maxspec)
      real*8 rhoe, enthe, cpb
      integer i,n,m,Niter,maxIter
      real*8 res(NiterMAX), lamOverCp

      dxInv2 = 1.d0/(dx*dx)
      do i=0,nx+1
         do n=1,Nspec
            Y(n,i) = S(FirstSpec+n-1,i) / S(Density,i) 
         enddo
      enddo

      maxiter=0
      do i=1,nx+1
         do n=1,Nspec
            de(n,i) = (Y(n,i)-Y(n,i-1)) / dx
         enddo
         de(Nspec+1,i) = (S(RhoH,i)-S(RhoH,i-1)) / dx

         Te = 0.5d0*(S(Temp,i)+S(Temp,i-1))
         rhoe = 0.5d0*(S(Density,i)+S(Density,i-1))
         if (setTfromH.eq.1) then
            enthe = 0.5d0*(S(RhoH,i)+S(RhoH,i-1))/rhoe
            call FORT_TfromHYpt(Te,enthe,Ye,Nspec,errMax,NiterMAX,res,Niter)
            if (Niter.lt.0) then
               print *,'RhoH->T failed at i=',i,
     &              'in neg_divF_H_Approx'
               stop
            endif
            maxiter=MAX(maxiter,Niter)
         endif
         call CKHMS(Te,IWRK,RWRK,He) 
         cpb = 0.d0
         do n=1,Nspec
            cpb = cpb + 0.5d0*(cpicc(n,i-1)+cpicc(n,i))
         enddo

         lamOverCp = PTCec(i) / cpb

         q(i) = 0.d0
         do n = 1,Nspec
            F(n,i) = - rhoDiec(n,i)*de(n,i)
            q(i) = q(i) + He(n)*(rhoDiec(n,i) - lamOverCp)*de(n,i)
         enddo
         q(i) = q(i) - PTCec(i)*de(Nspec+1,i)
      enddo

      do i=1,nx
         do n = 1,Nspec
            LofS(n,i) = - dxInv2*(F(n,i+1) - F(n,i))
         enddo
         LofS(Nspec+1,i) = - dxInv2*(q(i+1) - q(i))
      enddo

      end

      subroutine neg_divF_T(LofS,S,PTCec,rhoTDec,rhoDijec,cpicc,dx)
      implicit none
      include 'spec.h'
      real*8 LofS(maxspec+1,1:nx), S(maxscal,0:nx+1)
      real*8 PTCec(1:nx+1)
      real*8 rhoTDec(maxspec,1:nx+1),rhoDijec(maxspec,maxspec,1:nx+1),dx
      real*8 cpicc(maxspec,0:nx+1)
      real*8 Y(maxspec,0:nx+1), X(maxspec,0:nx+1), WWe, PTC
      real*8 de(maxspec+1,1:nx+1), q(1:nx+1), F(maxspec,1:nx+1)
      real*8 Ye(maxspec), Te, dxInv, sum
      real*8 rhoe, enthe, cpb
      integer i,n,m,Niter,maxIter
      real*8 res(NiterMAX)

      dxInv = 1.d0/dx
      do i=0,nx+1
         do n=1,Nspec
            Y(n,i) = S(FirstSpec+n-1,i) / S(Density,i) 
         enddo
         call CKYTX(Y(1,i),IWRK,RWRK,X(1,i))
      enddo

      maxiter=0
      do i=1,nx+1
         do n=1,Nspec
            de(n,i) = dxInv * (X(n,i)-X(n,i-1))
         enddo
         de(Nspec+1,i) = dxInv * (S(Temp,i)-S(Temp,i-1))

         do n=1,Nspec
            Ye(n) = 0.5d0*(Y(n,i)+Y(n,i-1))
         enddo
         call CKMMWY(Ye,IWRK,RWRK,WWe)

         Te = 0.5d0*(S(Temp,i)+S(Temp,i-1))
         rhoe = 0.5d0*(S(Density,i)+S(Density,i-1))
         if (setTfromH.eq.1) then
            enthe = 0.5d0*(S(RhoH,i)+S(RhoH,i-1))/rhoe
            call FORT_TfromHYpt(Te,enthe,Ye,Nspec,errMax,NiterMAX,res,Niter)
            if (Niter.lt.0) then
               print *,'RhoH->T failed at i=',i
               stop
            endif
            maxiter=MAX(maxiter,Niter)
         else if (setTfromH.eq.0) then
            Te = Pcgs * WWe / (rhoe * RU)
         endif

         q(i) = 0.d0
         do n = 1,Nspec
            F(n,i) = - Ye(n)*rhoTDec(n,i)*de(Nspec+1,i)/Te
            do m = 1,Nspec
               F(n,i) = F(n,i) - rhoDijec(n,m,i)*de(m,i)
            enddo
            q(i) = q(i) - (RU*Te/WWe)*rhoTDec(n,i)*de(n,i)
         enddo
         q(i) = q(i) - PTCec(i)*de(Nspec+1,i)
      enddo

      do i=1,nx
         do n = 1,Nspec
            LofS(n,i) = - dxInv*(F(n,i+1) - F(n,i))
         enddo
         sum = 0.d0
         cpb = 0.d0
         do n = 1,Nspec
            cpb = cpb + Y(n,i)*cpicc(n,i)
            sum = sum - 0.25d0*dxInv*(
     &           + F(n,i+1)*(cpicc(n,i  )+cpicc(n,i+1))*(S(Temp,i+1)-S(Temp,i  ) )
     &           + F(n,i  )*(cpicc(n,i-1)+cpicc(n,i  ))*(S(Temp,i  )-S(Temp,i-1) ) )
         enddo
         LofS(Nspec+1,i) = -(dxInv * (q(i+1) - q(i))  -  sum)/(S(Density,i) * cpb)
c         LofS(Nspec+1,i) = -(dxInv * (q(i+1) - q(i))) 
      enddo
      end

      subroutine neg_divF_T_Approx(LofS,S,PTCec,rhoDiec,cpicc,dx)
      implicit none
      include 'spec.h'
      real*8 LofS(maxspec+1,1:nx), S(maxscal,0:nx+1)
      real*8 PTCec(1:nx+1), rhoDiec(maxspec,1:nx+1),cpicc(maxspec,0:nx+1),dx
      real*8 q(1:nx+1), F(maxspec,1:nx+1),dxInv,sum,cpb
      integer i,n,m
      dxInv = 1.d0/dx
      do i=1,nx+1
         do n = 1,Nspec
            F(n,i) = -rhoDiec(n,i)*dxInv
     &           * ( S(FirstSpec+n-1,i)/S(Density,i) - S(FirstSpec+n-1,i-1)/S(Density,i-1) )
         enddo
         q(i) = -PTCec(i) * dxInv * ( S(Temp,i) - S(Temp,i-1) )
      enddo

      do i=1,nx
         do n = 1,Nspec
            LofS(n,i) = -dxInv*(F(n,i+1) - F(n,i))
         enddo
         sum = 0.d0         
         cpb = 0.d0
         do n = 1,Nspec
            cpb = cpb + cpicc(n,i)*S(FirstSpec+n-1,i)/S(Density,i  )
            sum = sum - 0.25d0*dxInv*(
     &           + F(n,i+1)*(cpicc(n,i  )+cpicc(n,i+1))*(S(Temp,i+1)-S(Temp,i  ) )
     &           + F(n,i  )*(cpicc(n,i-1)+cpicc(n,i  ))*(S(Temp,i  )-S(Temp,i-1) ) )
         enddo
         LofS(Nspec+1,i) = -(dxInv * (q(i+1) - q(i))  -  sum)/(S(Density,i) * cpb)
      enddo

      end

