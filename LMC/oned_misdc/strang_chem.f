      subroutine strang_chem(scal_old,scal_new,
     $                       const_src,lin_src_old,lin_src_new,
     $                       I_R,dt)
      implicit none
      include 'spec.h'
      real*8     scal_old(-1:nx  ,nscal)
      real*8     scal_new(-1:nx  ,nscal)
      real*8    const_src( 0:nx-1,nscal)
      real*8  lin_src_old( 0:nx-1,nscal)
      real*8  lin_src_new( 0:nx-1,nscal)
      real*8          I_R( 0:nx-1,0:maxspec)
      real*8  dt
      
      integer i,is,n,ifail
      real*8 RYold(maxspec), RYnew(maxspec), Told, Tnew
      real*8 linSrcOLD(nscal),linSrcNEW(nscal), rho_half, rho_old
      integer FuncCount, do_diag, IWRK
      real*8 diag(maxreac),hmix,RWRK,TSAVE
      real*8 Yhalf(maxspec),Y(maxspec), cp
      real*8 HK(maxspec),HK_new(maxspec),HK_old(maxspec)
      
      integer NiterMAX, Niter
      parameter (NiterMAX = 30)
      double precision res(NiterMAX), errMAX

C debugging remove me
      real*8 sum, diff, rhoh_intg


c     Shut off diagnostics
      do_diag = 0

      if (nochem_hack) then
         write(*,*)'WARNING! nochem_hack--skipping reactions'
         return
      endif

      print *,'... chemistry'
c     Evolve chem over grid
      do i = 0,nx-1

C CEG:: done to mimick john's SNe code.  firstSpec should be the fuel (i think)
C I don't think this is appropriate for the more complex reactions network here
C$$$         if (scal_new(i,FirstSpec) .le. 0.d0) then
C$$$            do n = 1, nscal
C$$$               scal_new(i,n) = scal_old(i,n)+
C$$$     &              + dt*const_src(i,n)
C$$$     &              + 0.5d0*dt*(lin_src_old(i,n)+lin_src_new(i,n))
C$$$            enddo
C$$$            do n = 0,Nspec
C$$$               I_R(i,n) = 0.d0
C$$$            enddo
C$$$         else
         if (use_strang) then
C integrating Y_m's not rhoY_m
            do n = 1,Nspec
               RYold(n) = scal_old(i,FirstSpec+n-1)/scal_old(i,Density)
            enddo
            Told = scal_old(i,Temp)            

         else 

            rho_old = 0.d0
            do n = 1,Nspec
               RYold(n) = scal_old(i,FirstSpec+n-1)
               rho_old = rho_old + RYold(n)
            enddo
            Told = scal_old(i,Temp)
C     FIXME!!!
C            Told = scal_old(i,RhoH)

c     Set linear source terms in common for ode integrators access
            do n = 1,Nspec
               is = FirstSpec + n - 1
               c_0(n) = const_src(i,is) + lin_src_old(i,is)
               c_1(n) = (lin_src_new(i,is) - lin_src_old(i,is))/dt
            enddo
            c_0(0) = const_src(i,RhoH) + lin_src_old(i,RhoH)
            c_1(0) = (lin_src_new(i,RhoH) - lin_src_old(i,RhoH))/dt
            rhoh_INIT = scal_old(i,RhoH)

         endif

C         write(*,*)i
         call chemsolve(RYnew, Tnew, RYold, Told, FuncCount, dt,
     &                  diag, do_diag, ifail)
         if (ifail.ne.0) then
            print *,'solve failed, i=',i
            stop
         endif
CCC CEWG debugging fixme
C         if (i .eq. 125) then
C            write(*,*) 'hit 125'
C            stop
C         endif

         if(use_strang) then

            do n = 1,Nspec
               scal_new(i,FirstSpec+n-1) = RYnew(n)*scal_old(i,Density)
            enddo
            scal_new(i,Temp) = Tnew
            do n = 1,Nspec
               is = FirstSpec + n - 1
               I_R(i,n) = (scal_new(i,is)-scal_old(i,is)) / dt
            enddo

         else

            scal_new(i,Density) = 0.d0
            do n = 1,Nspec
               scal_new(i,Density) = scal_new(i,Density) + RYnew(n)
            enddo
            do n = 1,Nspec
               scal_new(i,FirstSpec+n-1) = RYnew(n)
               Y(n) = RYnew(n)/scal_new(i,Density)
            enddo

            if (use_temp_eqn) then

               scal_new(i,RhoH) = scal_old(i,RhoH)+
     &              + dt*const_src(i,RhoH)
     &              + 0.5d0*dt*(lin_src_old(i,RhoH)+lin_src_new(i,RhoH))
               hmix = scal_new(i,RhoH) / scal_new(i,Density)

               TSAVE = Tnew
               errMax = hmix_TYP * 1.e-10
               call FORT_TfromHYpt(Tnew,hmix,Y,
     &              Nspec,errMax,NiterMAX,res,Niter)
               if (Niter.lt.0) then
                  print *,'SC: H to T solve failed in F, Niter=',Niter
                  print *,'hmix_TYP:',hmix_TYP
                  stop
               endif
               scal_new(i,Temp) = Tnew

c     Define change in state due to chemistry.
               call CKHMS(0.5d0*(Told+Tnew),IWRK,RWRK,HK)
               call CKHMS(Told,IWRK,RWRK,HK_old)
               call CKHMS(Tnew,IWRK,RWRK,HK_new)

               I_R(i,0) = 0.d0
               sum = 0.d0
               do n = 1,Nspec
                  is = FirstSpec + n - 1
                  I_R(i,n) =
     $                 (scal_new(i,is)-scal_old(i,is)) / dt
     $                 - const_src(i,is)
     $                 - 0.5d0*(lin_src_old(i,is)+lin_src_new(i,is))
                  Yhalf(n) = 0.5d0*(scal_old(i,is)/rho_old
     &                 +           Y(n)/scal_new(i,Density))
                  I_R(i,0) = I_R(i,0) - ( HK(n)*const_src(i,is)
     &                 + 0.5d0*(HK_old(n)*lin_src_old(i,is)
     &                 + HK_new(n)*lin_src_new(i,is)) )
C     sum = sum - I_R(i,n) * HK(n)
C     write(19,*)i,HK_old(n),HK(n),HK_new(n)
               enddo
               rho_half = 0.5d0*(rho_old + scal_new(i,Density))
               CALL CKCPBS(0.5d0*(scal_old(i,Temp)+Tnew),Yhalf,
     &              IWRK,RWRK,cp)
C     I_R(i,0) = -I_R(i,0)/(cp*rho_half)

               I_R(i,0) =  
     $              (scal_new(i,Temp)-scal_old(i,Temp)) / dt
     $              - ( I_R(i,0) + const_src(i,RhoH)
     $              + 0.5d0*(lin_src_old(i,RhoH)+lin_src_new(i,RhoH)) )/
     $              (rho_half*cp)

C     diff = sum - I_R(i,0)
C     write(18,*)i,I_R(i,0),sum,diff

            else

               scal_new(i,Temp) =  Tnew
               CALL CKHBMS(Tnew,Y,IWRK,RWRK,scal_new(i,RhoH))
               scal_new(i,RhoH) = scal_new(i,RhoH)*scal_new(i,Density)
CCCCCCCCCCCCCCC
C               scal_new(i,RhoH) = 
C$$$               scal_new(i,RhoH) = scal_old(i,RhoH)+
C$$$     &              + dt*const_src(i,RhoH)
C$$$     &              + 0.5d0*dt*(lin_src_old(i,RhoH)+lin_src_new(i,RhoH))
C$$$               hmix = scal_new(i,RhoH) / scal_new(i,Density)
C$$$               errMax = hmix_TYP * 1.e-10
C$$$C               call FORT_TfromHYpt(scal_new(i,Temp),hmix,Y,
C$$$C     &              Nspec,errMax,NiterMAX,res,Niter)
C$$$               call FORT_TfromHYpt(Tnew,hmix,Y,
C$$$     &              Nspec,errMax,NiterMAX,res,Niter)
C$$$               if (Niter.lt.0) then
C$$$                  print *,'SC: H to T solve failed in F, Niter=',Niter
C$$$                  print *,'hmix_TYP:',hmix_TYP
C$$$                  stop
C$$$               endif
C$$$               scal_new(i,Temp) = Tnew
CCCCCCCCCCc

C CEG just chekcing FIXME-- biggest difference ~5/10^5 
C$$$               rhoh_intg = scal_old(i,RhoH)+
C$$$     &              + dt*const_src(i,RhoH)
C$$$     &              + 0.5d0*dt*(lin_src_old(i,RhoH)+lin_src_new(i,RhoH))
C$$$               write(*,*) scal_new(i,RhoH), rhoh_intg
C$$$               write(*,*) scal_new(i,RhoH)-rhoh_intg

               do n = 1,Nspec
                  is = FirstSpec + n - 1
                  I_R(i,n) =
     $                 (scal_new(i,is)-scal_old(i,is)) / dt
     $                 - const_src(i,is)
     $                 - 0.5d0*(lin_src_old(i,is)+lin_src_new(i,is))
               enddo

            endif

         endif
C         endif
      enddo

      end

      subroutine chem(scal_old,scal_new,const_src,lin_src,I_R,dt)
      implicit none
      include 'spec.h'
      real*8     scal_old(-1:nx  ,nscal)
      real*8     scal_new(-1:nx  ,nscal)
      real*8    const_src( 0:nx-1,nscal)
      real*8      lin_src( 0:nx-1,nscal)
      real*8          I_R( 0:nx-1,0:maxspec)
      real*8  dt
      
      integer i,is,n,ifail
      real*8 RYold(maxspec), RYnew(maxspec), Told, Tnew
      real*8 linSrcOLD(nscal),linSrcNEW(nscal), rho_half, rho_old
      integer FuncCount, do_diag, IWRK
      real*8 diag(maxreac),hmix,RWRK,TSAVE
      real*8 Yhalf(maxspec),Y(maxspec), cp
      real*8 HK(maxspec),HK_new(maxspec),HK_old(maxspec)
      
      integer NiterMAX, Niter
      parameter (NiterMAX = 30)
      double precision res(NiterMAX), errMAX

C debugging remove me
      real*8 sum, diff, rhoh_intg


c     Shut off diagnostics
      do_diag = 0

      if (nochem_hack) then
         write(*,*)'WARNING! nochem_hack--skipping reactions'
         return
      endif

C      open(UNIT=11, FILE='ZP(0).dat', STATUS = 'REPLACE')

      print *,'... chemistry'
c     Evolve chem over grid
      do i = 0,nx-1
         
C         write(11,*)'i = ', i

C I don't think this is appropriate for the more complex reactions network here
C CEG:: done to mimick john's SNe code.  firstSpec should be the fuel (i think)
C$$$         if (scal_new(i,FirstSpec) .le. 0.d0) then
C$$$            do n = 1, nscal
C$$$               scal_new(i,n) = scal_old(i,n) + dt*const_src(i,n)
C$$$     &              + 0.5d0*dt*dt*lin_src(i,n)               
C$$$            enddo
C$$$         else
         if (use_strang) then

C integrating Y_m's not rhoY_m
            do n = 1,Nspec
               RYold(n) = scal_old(i,FirstSpec+n-1)/scal_old(i,Density)
            enddo
            Told = scal_old(i,Temp)            

         else 

            rho_old = 0.d0
            do n = 1,Nspec
               RYold(n) = scal_old(i,FirstSpec+n-1)
               rho_old = rho_old + RYold(n)
            enddo
            Told = scal_old(i,Temp)

c     Set linear source terms in common for ode integrators access
            do n = 1,Nspec
               is = FirstSpec + n - 1
               c_0(n) = const_src(i,is)
               c_1(n) = lin_src(i,is)
            enddo
            c_0(0) = const_src(i,RhoH) 
            c_1(0) = lin_src(i,RhoH)
            rhoh_INIT = scal_old(i,RhoH)

         endif

         call chemsolve(RYnew, Tnew, RYold, Told, FuncCount, dt,
     &                  diag, do_diag, ifail)
         if (ifail.ne.0) then
            print *,'solve failed, i=',i
            stop
         endif

         if(use_strang) then

            do n = 1,Nspec
               scal_new(i,FirstSpec+n-1) = RYnew(n)*scal_old(i,Density)
            enddo
            scal_new(i,Temp) = Tnew
            do n = 1,Nspec
               is = FirstSpec + n - 1
               I_R(i,n) = (scal_new(i,is)-scal_old(i,is)) / dt
            enddo

         else

            scal_new(i,Density) = 0.d0
            do n = 1,Nspec
               scal_new(i,Density) = scal_new(i,Density) + RYnew(n)
            enddo
            do n = 1,Nspec
               scal_new(i,FirstSpec+n-1) = RYnew(n)
               Y(n) = RYnew(n)/scal_new(i,Density)
            enddo

C CEG this was not a good idea--leads to bigger kinks
C            call CKRHOPY(scal_new(i,Density),Pcgs,Y,IWRK,RWRK,Tnew)
            scal_new(i,Temp) =  Tnew
            CALL CKHBMS(Tnew,Y,IWRK,RWRK,scal_new(i,RhoH))
            scal_new(i,RhoH) = scal_new(i,RhoH)*scal_new(i,Density)

CCCCCCCCCCCC FIXME 
C$$$            scal_new(i,RhoH) = scal_old(i,RhoH)
C$$$     &           + dt*const_src(i,RhoH)
C$$$     &           + 0.5d0*dt*dt*lin_src(i,RhoH)

C$$$            hmix = scal_new(i,RhoH) / scal_new(i,Density)
C$$$            errMax = hmix_TYP * 1.e-10
C$$$            do n = 1,Nspec
C$$$               Y(n) = scal_new(i,FirstSpec+n-1)/scal_new(i,Density)
C$$$            enddo
C$$$            call FORT_TfromHYpt(scal_new(i,Temp),hmix,Y,
C$$$     &           Nspec,errMax,NiterMAX,res,Niter)
C$$$            if (Niter.lt.0) then
C$$$               print *,'SC: H to T solve failed in F, Niter=',Niter
C$$$               print *,'hmix_TYP:',hmix_TYP
C$$$               stop
C$$$            endif
CCCCCCCCCCCCccCCC

            do n = 1,Nspec
               is = FirstSpec + n - 1
               I_R(i,n) =
     $              (scal_new(i,is)-scal_old(i,is)) / dt
     $              - const_src(i,is)
     $              - 0.5d0*dt*lin_src(i,is)
            enddo
            
         endif
C         endif
C         write(11,*)
C         write(11,*)
      enddo
C      close (11)
      end
