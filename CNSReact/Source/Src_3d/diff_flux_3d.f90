      subroutine diffFlux(lo,hi,BL_FARG(q), &
                          BL_FARG(flux1),BL_FARG(flux2),BL_FARG(flux3), &
                          dx,dy,dz,dt,bcx,bcy,bcz,initFlux)

      ! Get diffusion flux on faces surrounding Box(lo,hi), requires 
      ! valid data on region Box(lo,hi).grow(1)

      implicit none

      integer lo(SDIM),hi(SDIM)
      BL_FBOUNDS(q)
      BL_FBOUNDS(flux1)
      BL_FBOUNDS(flux2)
      BL_FBOUNDS(flux3)
      BL_FARRAY(q,QVAR)
      BL_FARRAY(flux2,NVAR)
      BL_FARRAY(flux1,NVAR)
      BL_FARRAY(flux3,NVAR)
      REAL_T dx,dy,dz,dt
      integer bcx, bcy, bcz, initFlux
      integer i,j,k

#define DVAR (3+NSPECMAX)
      REAL_T  D(QSTATE*DVAR)
      REAL_T TEMP(QSTATE)
      REAL_T CP(QSTATE)
     
      integer loD(3), hiD(3)
      integer n

      do n = 1,3
         loD(n) = lo(n) - 1
         hiD(n) = hi(n) + 1
      end do

      if (initFlux.eq.1) then
      do n=1,NVAR
      do k=lo(3),hi(3)
      do j=lo(2),hi(2)
      do i=lo(1),hi(1)+1
          flux1(i,j,k,n) = 0.
      enddo
      enddo
      enddo
      enddo

      do n=1,NVAR
      do k=lo(3),hi(3)
      do j=lo(2),hi(2)+1
      do i=lo(1),hi(1)
          flux2(i,j,k,n) = 0.
      enddo
      enddo
      enddo
      enddo

      do n=1,NVAR
      do k=lo(3),hi(3)+1
      do j=lo(2),hi(2)
      do i=lo(1),hi(1)
          flux3(i,j,k,n) = 0.
      enddo
      enddo
      enddo
      enddo
      endif

      call transCoef(loD,hiD,BL_FARG(q),TEMP,CP, &
                     D,loD(1),loD(2),loD(3),hiD(1),hiD(2),hiD(3))

      call myViscFlux(lo,hi,BL_FARG(q), &
                      D,loD(1),loD(2),loD(3),hiD(1),hiD(2),hiD(3), &
                      BL_FARG(flux1),BL_FARG(flux2),BL_FARG(flux3), &
                      dx,dy,dz,dt,bcx,bcy,bcz)

      end subroutine diffFlux
#endif /* NSPECMAX > 0 */

#if NSPECMAX > 0
      subroutine transCoef(lo,hi,BL_FARG(q),TEMP,CP,BL_FARG(D))

      use cdwrk_module
      implicit none
      integer lo(SDIM), hi(SDIM)
      BL_FBOUNDS(q)
      BL_FBOUNDS(D)
      BL_FARRAY(q,QVAR)
      BL_FARRAY(D,DVAR)
      REAL_T TEMP(DIMV(q))
      REAL_T CP(DIMV(q))

      integer i,j,k,n,doVelVisc,doTemp
      REAL_T Dmax
      
#if 1
      integer sizeq,sizeD
      sizeq = (ARG_H1(q)-ARG_L1(q)+1)*
     $        (ARG_H2(q)-ARG_L2(q)+1)*
     $        (ARG_H3(q)-ARG_L3(q)+1)
      if (sizeq>QSTATE) then
         write(*,*) 'not enough work space for q'
         stop
      end if
      if (sizeq>QSTATE) then
         write(*,*) 'not enough work space for TEMP'
         stop
      end if
      sizeD = (ARG_H1(D)-ARG_L1(D)+1)*
     $        (ARG_H2(D)-ARG_L2(D)+1)*
     $        (ARG_H3(D)-ARG_L3(D)+1)
      if (sizeD>QSTATE) then
         write(*,*) 'not enough work space for D'
         stop
      end if
#endif

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               TEMP(i,j,k) = MIN(MAX(q(i,j,k,QTHERM+1),300.d0),6500.d0)
            end do
         end do
      end do
      
      doTemp = 1
      doVelVisc = 1
      call FORT_MIXAVG_RHODIFF_TEMP_PRES(lo, hi,
     &     D,   DIMS(D),
     &     TEMP,DIMS(q),
     &     q(BL_ARGL1(q),BL_ARGL2(q),BL_ARGL3(q),QTHERM+NADV+1),DIMS(q),
     &     q(BL_ARGL1(q),BL_ARGL2(q),BL_ARGL3(q),QP),           DIMS(q),
     &     doTemp, doVelVisc)

      call FORT_CPMIXfromTY(lo, hi,
     &     CP,  DIMS(q),
     &     TEMP,DIMS(q),
     &     q(BL_ARGL1(q),BL_ARGL2(q),BL_ARGL3(q),QTHERM+NADV+1),DIMS(q))

c ::: Replace lambda with lambda/Cp
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               D(i,j,k,Nspec+1) = D(i,j,k,Nspec+1)/CP(i,j,k)
            end do
         end do
      end do
#if 0
      Dmax = -1.d0
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               do n=1,Nspec
                  Dmax = MAX(Dmax,D(i,j,k,n) / q(i,j,k,QR))
               enddo
            end do
         end do
      end do
      print*, '********************      Dmax: ', Dmax
#endif
      end

#if NSPECMAX > 0
      subroutine myViscFlux(lo,hi,BL_FARG(q),BL_FARG(D), &
                               BL_FARG(fluxx),BL_FARG(fluxy),BL_FARG(fluxz), &
                               dx,dy,dz,dt,bcx,bcy,bcz)
      implicit none

      use cdwrk_module

      integer lo(SDIM), hi(SDIM)
      BL_FBOUNDS(q)
      BL_FBOUNDS(D)
      BL_FBOUNDS(fluxx)
      BL_FBOUNDS(fluxy)
      BL_FBOUNDS(fluxz)
      BL_FARRAY(q,QVAR)
      BL_FARRAY(D,DVAR)
      BL_FARRAY(fluxx,NVAR)
      BL_FARRAY(fluxy,NVAR)
      BL_FARRAY(fluxz,NVAR)
      integer bcx,bcy,bcz
      REAL_T dx,dy,dz,dt

c     Local variables
      integer i,j,k, iv
      integer is,ie,js,je,ks,ke,istart,kstart
      REAL_T ux,uy,uz,vx,vy,vz,wx,wy,wz,tx,ty,tz
      REAL_T tauxxp,tauxxm,tauxyp,tauxym,tauxzp,tauxzm
      REAL_T tauyyp,tauyym,tauyxp,tauyxm,tauyzp,tauyzm
      REAL_T tauzzp,tauzzm,tauzxp,tauzxm,tauzyp,tauzym
      REAL_T divxp, divxm,divyp,divym,divzp,divzm
      REAL_T muxp,muxm,muyp,muym,muzp,muzm
      REAL_T lamxp,lamxm,lamyp,lamym,lamzp,lamzm
      REAL_T kxp,kxm,kyp,kym,kzp,kzm
      REAL_T phiflx
#if 1
      integer sizeD,sizeq,sizefx,sizefy,sizefz
#endif
#if NSPECMAX > 0
      REAL_T rhoDm(maxspec)
      integer n
      integer ifirstSp
      ifirstSp = NTHERM+NADV+1
#endif /* NSPECMAX > 0 */

#if 1
      sizeD = (ARG_H1(D)-ARG_L1(D)+1)*
     $        (ARG_H2(D)-ARG_L2(D)+1)*
     $        (ARG_H3(D)-ARG_L3(D)+1)
      if (sizeD>QSTATE) then
         write(*,*) 'not enough work space for D'
         stop
      end if
      sizeq = (ARG_H1(q)-ARG_L1(q)+1)*
     $        (ARG_H2(q)-ARG_L2(q)+1)*
     $        (ARG_H3(q)-ARG_L3(q)+1)
      if (sizeq>QSTATE) then
         write(*,*) 'not enough work space for q'
         stop
      end if
      sizefx = (ARG_H1(fluxx)-ARG_L1(fluxx)+1)*
     $         (ARG_H2(fluxx)-ARG_L2(fluxx)+1)*
     $         (ARG_H3(fluxx)-ARG_L3(fluxx)+1)
      if (sizefx>DFLUX) then
         write(*,*) 'not enough work space for fluxx'
         stop
      end if
      sizefy = (ARG_H1(fluxy)-ARG_L1(fluxy)+1)*
     $         (ARG_H2(fluxy)-ARG_L2(fluxy)+1)*
     $         (ARG_H3(fluxy)-ARG_L3(fluxy)+1)
      if (sizefy>DFLUX) then
         write(*,*) 'not enough work space for fluxy'
         stop
      end if
      sizefz = (ARG_H1(fluxz)-ARG_L1(fluxz)+1)*
     $         (ARG_H2(fluxz)-ARG_L2(fluxz)+1)*
     $         (ARG_H3(fluxz)-ARG_L3(fluxz)+1)
      if (sizefz>DFLUX) then
         write(*,*) 'not enough work space for fluxz'
         stop
      end if
#endif

      is = lo(1) 
      ie = hi(1)
      js = lo(2)
      je = hi(2)
      ks = lo(3)
      ke = hi(3)

c     compute x fluxes
      do k = ks,ke
         do j = js,je
            do i = is,ie+1
               
               muxm =  half*(D(i-1,j,k,Nspec+2) + D(i,j,k,Nspec+2))
               lamxm = -two3rd*muxm
               kxm =   half*(D(i-1,j,k,Nspec+1) + D(i,j,k,Nspec+1))
               do n=1,Nspec
                  rhoDm(n) = half*(D(i-1,j,k,n) + D(i,j,k,n))
               end do
               
               tauxxm = two*muxm*(q(i,j,k,QU) - q(i-1,j,k,QU))/dx 
               
               divxm = lamxm*(q(i,j,k,QU) - q(i-1,j,k,QU))/dx 
     &              +.25d0*lamxm*(q(i-1,j+1,k,QV)+q(i,j+1,k,QV)
     &              -             q(i-1,j-1,k,QV)-q(i,j-1,k,QV))/dy
     &              +.25d0*lamxm*(q(i-1,j,k+1,QW)+q(i,j,k+1,QW)
     &              -             q(i-1,j,k-1,QW)-q(i,j,k-1,QW))/dz
               
               tauyxm =  muxm*(q(i,j,k,QV) - q(i-1,j,k,QV))/dx 
     &              +.25d0*muxm*(q(i,j+1,k,QU)+q(i-1,j+1,k,QU)
     &              -            q(i,j-1,k,QU)-q(i-1,j-1,k,QU))/dy
               
               tauzxm =  muxm*(q(i,j,k,QW) - q(i-1,j,k,QW))/dx 
     &              +.25d0*muxm*(q(i,j,k+1,QU)+q(i-1,j,k+1,QU)
     &              -            q(i,j,k-1,QU)-q(i-1,j,k-1,QU))/dz

               phiflx =  tauxxm*(q(i,j,k,QU)+q(i-1,j,k,QU))
     &              +     divxm*(q(i,j,k,QU)+q(i-1,j,k,QU))
     &              +    tauyxm*(q(i,j,k,QV)+q(i-1,j,k,QV))
     &              +    tauzxm*(q(i,j,k,QW)+q(i-1,j,k,QW))
               
               fluxx(i,j,k,2) = fluxx(i,j,k,2) - dt*dy*dz*(tauxxm+divxm)
               fluxx(i,j,k,3) = fluxx(i,j,k,3) - dt*dy*dz*tauyxm
               fluxx(i,j,k,4) = fluxx(i,j,k,4) - dt*dy*dz*tauzxm
               fluxx(i,j,k,5) = fluxx(i,j,k,5) - dt*dy*dz*(half*phiflx 
     &              + kxm*( (q(i  ,j,k,QRE)+q(i  ,j,k,QP))/q(i  ,j,k,QR)
     &              -       (q(i-1,j,k,QRE)+q(i-1,j,k,QP))/q(i-1,j,k,QR) )/dx)

               do n=1,Nspec
                  fluxx(i,j,k,n+ifirstSp-1) = fluxx(i,j,k,n+ifirstSp-1)
     &                 - dt*dy*dz*rhoDm(n)*
     &                 (q(i,j,k,QTHERM+NADV+n)-q(i-1,j,k,QTHERM+NADV+n))/dx
               end do
               
            enddo
         enddo
      enddo

c     compute y fluxes
      do k = ks,ke
         do j = js,je+1
            do i = is,ie
               
               muym =  half*(D(i,j-1,k,Nspec+2) + D(i,j,k,Nspec+2))
               lamym = -two3rd*muym
               kym =   half*(D(i,j-1,k,Nspec+1) + D(i,j,k,Nspec+1))
               do n=1,Nspec
                  rhoDm(n) = half*(D(i,j-1,k,n) + D(i,j,k,n))
               end do
               
               tauxym =  muym*(q(i,j,k,QU) - q(i,j-1,k,QU))/dy 
     &              +.25d0*muym*(q(i+1,j,k,QV)+q(i+1,j-1,k,QV)
     &              -            q(i-1,j,k,QV)-q(i-1,j-1,k,QV))/dx
               
               tauyym = two*muym*(q(i,j,k,QV) - q(i,j-1,k,QV))/dy 
               
               divym = lamym*(q(i,j,k,QV) - q(i,j-1,k,QV))/dy 
     &              +.25d0*lamym*(q(i+1,j-1,k,QU)+q(i+1,j,k,QU)
     &              -             q(i-1,j-1,k,QU)-q(i-1,j,k,QU))/dx
     &              +.25d0*lamym*(q(i,j-1,k+1,QW)+q(i,j,k+1,QW)
     &              -             q(i,j-1,k-1,QW)-q(i,j,k-1,QW))/dz
               
               tauzym =  muym*(q(i,j,k,QW) - q(i,j-1,k,QW))/dy 
     &              +.25d0*muym*(q(i,j,k+1,QV)+q(i,j-1,k+1,QV)
     &              -            q(i,j,k-1,QV)-q(i,j-1,k-1,QV))/dz

               phiflx = tauxym*(q(i,j,k,QU)+q(i,j-1,k,QU))
     &              +   tauyym*(q(i,j,k,QV)+q(i,j-1,k,QV))
     &              +    divym*(q(i,j,k,QV)+q(i,j-1,k,QV))
     &              +   tauzym*(q(i,j,k,QW)+q(i,j-1,k,QW))
               
               fluxy(i,j,k,2) = fluxy(i,j,k,2) - dt*dx*dz*tauxym
               fluxy(i,j,k,3) = fluxy(i,j,k,3) - dt*dx*dz*(tauyym+divym)
               fluxy(i,j,k,4) = fluxy(i,j,k,4) - dt*dx*dz*tauzym
               fluxy(i,j,k,5) = fluxy(i,j,k,5) - dt*dy*dz*(half*phiflx 
     &              + two*kym*( (q(i,j  ,k,QRE)+q(i,j  ,k,QP))/q(i,j  ,k,QR)
     &              -           (q(i,j-1,k,QRE)+q(i,j-1,k,QP))/q(i,j-1,k,QR) )/dy)

               do n=1,Nspec
                  fluxy(i,j,k,n+ifirstSp-1) = fluxy(i,j,k,n+ifirstSp-1)
     &                 - dt*dx*dz*rhoDm(n)*
     &                 (q(i,j,k,QTHERM+NADV+n)-q(i,j-1,k,QTHERM+NADV+n))/dy
               end do
               
            enddo
         enddo
      enddo

c     compute z fluxes
      do k = ks,ke+1
         do j = js,je
            do i = is,ie
               
               muzm =  half*(D(i,j,k-1,Nspec+2) + D(i,j,k,Nspec+2))
               lamzm = -two3rd*muzm
               kzm =   half*(D(i,j,k-1,Nspec+1) + D(i,j,k,Nspec+1))
               do n=1,Nspec
                  rhoDm(n) = half*(D(i,j,k-1,n) + D(i,j,k,n))
               end do
               
               tauxzm =  muzm*(q(i,j,k,QU) - q(i,j,k-1,QU))/dz 
     &              +.25d0*muzm*(q(i+1,j,k,QW)+q(i+1,j,k-1,QW)
     &              -            q(i-1,j,k,QW)-q(i-1,j,k-1,QW))/dx

               tauyzm =  muzm*(q(i,j,k,QV) - q(i,j,k-1,QV))/dz 
     &              +.25d0*muzm*(q(i,j+1,k,QW)+q(i,j+1,k-1,QW)
     &              -            q(i,j-1,k,QW)-q(i,j-1,k-1,QW))/dy
               
               tauzzm = two*muzm*(q(i,j,k,QW) - q(i,j,k-1,QW))/dz 

               divzm = lamzm*(q(i,j,k,QW) - q(i,j,k-1,QW))/dz 
     &              +.25d0*lamzm*(q(i+1,j,k-1,QU)+q(i+1,j,k,QU)
     &              -             q(i-1,j,k-1,QU)-q(i-1,j,k,QU))/dx
     &              +.25d0*lamzm*(q(i,j+1,k-1,QV)+q(i,j+1,k,QV)
     &              -             q(i,j-1,k-1,QV)-q(i,j-1,k,QV))/dy
               
               phiflx =tauxzm*(q(i,j,k,QU)+q(i,j,k-1,QU))
     &              +  tauyzm*(q(i,j,k,QV)+q(i,j,k-1,QV))
     &              +  tauzzm*(q(i,j,k,QW)+q(i,j,k-1,QW))
     &              +   divzm*(q(i,j,k,QW)+q(i,j,k-1,QW))

               fluxz(i,j,k,2) = fluxz(i,j,k,2) - dt*dx*dy*tauxzm
               fluxz(i,j,k,3) = fluxz(i,j,k,3) - dt*dx*dy*tauyzm
               fluxz(i,j,k,4) = fluxz(i,j,k,4) - dt*dx*dy*(tauzzm+divzm)
               fluxz(i,j,k,5) = fluxz(i,j,k,5) - dt*dx*dy*(half*phiflx
     &              + two*kzm*( (q(i,j,k  ,QRE)+q(i,j,k  ,QP))/q(i,j,k  ,QR)
     &              -           (q(i,j,k-1,QRE)+q(i,j,k-1,QP))/q(i,j,k-1,QR) )/dz)

               do n=1,Nspec
                  fluxz(i,j,k,n+ifirstSp -1) = fluxz(i,j,k,n+ifirstSp -1)
     &                 - dt*dx*dy*rhoDm(n)*
     &                 (q(i,j,k,QTHERM+NADV+n)-q(i,j,k-1,QTHERM+NADV+n))/dz
               end do
               
            enddo
         enddo
      enddo
      end

#endif /* NSPECMAX > 0 */
