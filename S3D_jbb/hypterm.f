      module constants

      double precision, parameter :: alp = .8d0
      double precision, parameter :: bet = -.2d0
      double precision, parameter :: gam = 4.d0/105.d0
      double precision, parameter :: del = -1.d0/280.d0

      integer, parameter :: irho = 1
      integer, parameter :: imx  = 2
      integer, parameter :: imy  = 3
      integer, parameter :: imz  = 4
      integer, parameter :: iene = 5
      integer, parameter :: isp  = 6

      end module constants

      subroutine hypterm(lo,hi,ng,nspec,dx,cons,pres,flux)

      use constants

      implicit none

c     inputs:  lo,hi,ng,nspec,cons,pres
c     outputs: flux

      integer lo(3),hi(3),ng,nspec

      double precision dx(3)

      double precision pres(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,
     1             -ng+lo(3):hi(3)+ng)
      double precision cons(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,
     1             -ng+lo(3):hi(3)+ng,nspec+5)
      double precision flux(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),nspec+5)

      integer i,j,k,nsp

      double precision unp1,unp2,unp3,unp4,unm1,unm2,unm3,unm4

      do k=lo(3),hi(3)
      do j=lo(2),hi(2)
      do i=lo(1),hi(1)

      unp1 = cons(i+1,j,k,imx)/cons(i+1,j,k,irho)
      unp2 = cons(i+2,j,k,imx)/cons(i+2,j,k,irho)
      unp3 = cons(i+3,j,k,imx)/cons(i+3,j,k,irho)
      unp4 = cons(i+4,j,k,imx)/cons(i+4,j,k,irho)

      unm1 = cons(i-1,j,k,imx)/cons(i-1,j,k,irho)
      unm2 = cons(i-2,j,k,imx)/cons(i-2,j,k,irho)
      unm3 = cons(i-3,j,k,imx)/cons(i-3,j,k,irho)
      unm4 = cons(i-4,j,k,imx)/cons(i-4,j,k,irho)

      flux(i,j,k,irho)=flux(i,j,k,irho) -
     1   (alp*(cons(i+1,j,k,imx)-cons(i-1,j,k,imx))
     1  + bet*(cons(i+2,j,k,imx)-cons(i-2,j,k,imx))
     1  + gam*(cons(i+3,j,k,imx)-cons(i-3,j,k,imx))
     1  + del*(cons(i+4,j,k,imx)-cons(i-4,j,k,imx)))/dx(1)

      flux(i,j,k,imx)=flux(i,j,k,imx) -
     1   (alp*(cons(i+1,j,k,imx)*unp1-cons(i-1,j,k,imx)*unm1
     1      + (pres(i+1,j,k)-pres(i-1,j,k)))
     1  + bet*(cons(i+2,j,k,imx)*unp2-cons(i-2,j,k,imx)*unm2
     1      + (pres(i+2,j,k)-pres(i-2,j,k)))
     1  + gam*(cons(i+3,j,k,imx)*unp3-cons(i-3,j,k,imx)*unm3
     1      + (pres(i+3,j,k)-pres(i-3,j,k)))
     1  + del*(cons(i+4,j,k,imx)*unp4-cons(i-4,j,k,imx)*unm4
     1      + (pres(i+4,j,k)-pres(i-4,j,k)))) / dx(1)

      flux(i,j,k,imy)=flux(i,j,k,imy) -
     1   (alp*(cons(i+1,j,k,imy)*unp1-cons(i-1,j,k,imy)*unm1)
     1  + bet*(cons(i+2,j,k,imy)*unp2-cons(i-2,j,k,imy)*unm2)
     1  + gam*(cons(i+3,j,k,imy)*unp3-cons(i-3,j,k,imy)*unm3)
     1  + del*(cons(i+4,j,k,imy)*unp4-cons(i-4,j,k,imy)*unm4))/dx(1)

      flux(i,j,k,imz)=flux(i,j,k,imz) -
     1   (alp*(cons(i+1,j,k,imz)*unp1-cons(i-1,j,k,imz)*unm1)
     1  + bet*(cons(i+2,j,k,imz)*unp2-cons(i-2,j,k,imz)*unm2)
     1  + gam*(cons(i+3,j,k,imz)*unp3-cons(i-3,j,k,imz)*unm3)
     1  + del*(cons(i+4,j,k,imz)*unp4-cons(i-4,j,k,imz)*unm4))/dx(1)

      flux(i,j,k,iene)=flux(i,j,k,iene) -
     1   (alp*(cons(i+1,j,k,iene)*unp1-cons(i-1,j,k,iene)*unm1
     1      + (pres(i+1,j,k)*unp1-pres(i-1,j,k)*unm1))
     1  + bet*(cons(i+2,j,k,iene)*unp2-cons(i-2,j,k,iene)*unm2
     1      + (pres(i+2,j,k)*unp2-pres(i-2,j,k)*unm2))
     1  + gam*(cons(i+3,j,k,iene)*unp3-cons(i-3,j,k,iene)*unm3
     1      + (pres(i+3,j,k)*unp3-pres(i-3,j,k)*unm3))
     1  + del*(cons(i+4,j,k,iene)*unp4-cons(i-4,j,k,iene)*unm4
     1      + (pres(i+4,j,k)*unp4-pres(i-4,j,k)*unm4))) / dx(1)

      do nsp = 0,nspec-1

         flux(i,j,k,isp+nsp)=flux(i,j,k,isp+nsp) -
     1      (alp*(cons(i+1,j,k,isp+nsp)*unp1-cons(i-1,j,k,isp+nsp)*unm1)
     1     + bet*(cons(i+2,j,k,isp+nsp)*unp2-cons(i-2,j,k,isp+nsp)*unm2)
     1     + gam*(cons(i+3,j,k,isp+nsp)*unp3-cons(i-3,j,k,isp+nsp)*unm3)
     1     + del*(cons(i+4,j,k,isp+nsp)*unp4-cons(i-4,j,k,isp+nsp)*unm4)
     1       )/dx(1)

      enddo

      enddo
      enddo
      enddo

      do k=lo(3),hi(3)
      do j=lo(2),hi(2)
      do i=lo(1),hi(1)

      unp1 = cons(i,j+1,k,imy)/cons(i,j+1,k,irho)
      unp2 = cons(i,j+2,k,imy)/cons(i,j+2,k,irho)
      unp3 = cons(i,j+3,k,imy)/cons(i,j+3,k,irho)
      unp4 = cons(i,j+4,k,imy)/cons(i,j+4,k,irho)

      unm1 = cons(i,j-1,k,imy)/cons(i,j-1,k,irho)
      unm2 = cons(i,j-2,k,imy)/cons(i,j-2,k,irho)
      unm3 = cons(i,j-3,k,imy)/cons(i,j-3,k,irho)
      unm4 = cons(i,j-4,k,imy)/cons(i,j-4,k,irho)

      flux(i,j,k,irho)=flux(i,j,k,irho) -
     1   (alp*(cons(i,j+1,k,imy)-cons(i,j-1,k,imy))
     1  + bet*(cons(i,j+2,k,imy)-cons(i,j-2,k,imy))
     1  + gam*(cons(i,j+3,k,imy)-cons(i,j-3,k,imy))
     1  + del*(cons(i,j+4,k,imy)-cons(i,j-4,k,imy)))/dx(2)

      flux(i,j,k,imx)=flux(i,j,k,imx) -
     1   (alp*(cons(i,j+1,k,imx)*unp1-cons(i,j-1,k,imx)*unm1)
     1  + bet*(cons(i,j+2,k,imx)*unp2-cons(i,j-2,k,imx)*unm2)
     1  + gam*(cons(i,j+3,k,imx)*unp3-cons(i,j-3,k,imx)*unm3)
     1  + del*(cons(i,j+4,k,imx)*unp4-cons(i,j-4,k,imx)*unm4))/dx(2)

      flux(i,j,k,imy)=flux(i,j,k,imy) -
     1   (alp*(cons(i,j+1,k,imy)*unp1-cons(i,j-1,k,imy)*unm1
     1      + (pres(i,j+1,k)-pres(i,j-1,k)))
     1  + bet*(cons(i,j+2,k,imy)*unp2-cons(i,j-2,k,imy)*unm2
     1      + (pres(i,j+2,k)-pres(i,j-2,k)))
     1  + gam*(cons(i,j+3,k,imy)*unp3-cons(i,j-3,k,imy)*unm3
     1      + (pres(i,j+3,k)-pres(i,j-3,k)))
     1  + del*(cons(i,j+4,k,imy)*unp4-cons(i,j-4,k,imy)*unm4
     1      + (pres(i,j+4,k)-pres(i,j-4,k)))) / dx(2)

      flux(i,j,k,imz)=flux(i,j,k,imz) -
     1   (alp*(cons(i,j+1,k,imz)*unp1-cons(i,j-1,k,imz)*unm1)
     1  + bet*(cons(i,j+2,k,imz)*unp2-cons(i,j-2,k,imz)*unm2)
     1  + gam*(cons(i,j+3,k,imz)*unp3-cons(i,j-3,k,imz)*unm3)
     1  + del*(cons(i,j+4,k,imz)*unp4-cons(i,j-4,k,imz)*unm4))/dx(2)

      flux(i,j,k,iene)=flux(i,j,k,iene) -
     1   (alp*(cons(i,j+1,k,iene)*unp1-cons(i,j-1,k,iene)*unm1
     1      + (pres(i,j+1,k)*unp1-pres(i,j-1,k)*unm1))
     1  + bet*(cons(i,j+2,k,iene)*unp2-cons(i,j-2,k,iene)*unm2
     1      + (pres(i,j+2,k)*unp2-pres(i,j-2,k)*unm2))
     1  + gam*(cons(i,j+3,k,iene)*unp3-cons(i,j-3,k,iene)*unm3
     1      + (pres(i,j+3,k)*unp3-pres(i,j-3,k)*unm3))
     1  + del*(cons(i,j+4,k,iene)*unp4-cons(i,j-4,k,iene)*unm4
     1      + (pres(i,j+4,k)*unp4-pres(i,j-4,k)*unm4))) / dx(2)

      do nsp = 0,nspec-1

         flux(i,j,k,isp+nsp)=flux(i,j,k,isp+nsp) -
     1      (alp*(cons(i,j+1,k,isp+nsp)*unp1-cons(i,j-1,k,isp+nsp)*unm1)
     1     + bet*(cons(i,j+2,k,isp+nsp)*unp2-cons(i,j-2,k,isp+nsp)*unm2)
     1     + gam*(cons(i,j+3,k,isp+nsp)*unp3-cons(i,j-3,k,isp+nsp)*unm3)
     1     + del*(cons(i,j+4,k,isp+nsp)*unp4-cons(i,j-4,k,isp+nsp)*unm4)
     1       )/dx(2)

      enddo

      enddo
      enddo
      enddo


      do k=lo(3),hi(3)
      do j=lo(2),hi(2)
      do i=lo(1),hi(1)

      unp1 = cons(i,j,k+1,imz)/cons(i,j,k+1,irho)
      unp2 = cons(i,j,k+2,imz)/cons(i,j,k+2,irho)
      unp3 = cons(i,j,k+3,imz)/cons(i,j,k+3,irho)
      unp4 = cons(i,j,k+4,imz)/cons(i,j,k+4,irho)

      unm1 = cons(i,j,k-1,imz)/cons(i,j,k-1,irho)
      unm2 = cons(i,j,k-2,imz)/cons(i,j,k-2,irho)
      unm3 = cons(i,j,k-3,imz)/cons(i,j,k-3,irho)
      unm4 = cons(i,j,k-4,imz)/cons(i,j,k-4,irho)

      flux(i,j,k,irho)=flux(i,j,k,irho) -
     1   (alp*(cons(i,j,k+1,imz)-cons(i,j,k-1,imz))
     1  + bet*(cons(i,j,k+2,imz)-cons(i,j,k-2,imz))
     1  + gam*(cons(i,j,k+3,imz)-cons(i,j,k-3,imz))
     1  + del*(cons(i,j,k+4,imz)-cons(i,j,k-4,imz)))/dx(3)

      flux(i,j,k,imx)=flux(i,j,k,imx) -
     1   (alp*(cons(i,j,k+1,imx)*unp1-cons(i,j,k-1,imx)*unm1)
     1  + bet*(cons(i,j,k+2,imx)*unp2-cons(i,j,k-2,imx)*unm2)
     1  + gam*(cons(i,j,k+3,imx)*unp3-cons(i,j,k-3,imx)*unm3)
     1  + del*(cons(i,j,k+4,imx)*unp4-cons(i,j,k-4,imx)*unm4))/dx(3)

      flux(i,j,k,imy)=flux(i,j,k,imy) -
     1   (alp*(cons(i,j,k+1,imy)*unp1-cons(i,j,k-1,imy)*unm1)
     1  + bet*(cons(i,j,k+2,imy)*unp2-cons(i,j,k-2,imy)*unm2)
     1  + gam*(cons(i,j,k+3,imy)*unp3-cons(i,j,k-3,imy)*unm3)
     1  + del*(cons(i,j,k+4,imy)*unp4-cons(i,j,k-4,imy)*unm4))/dx(3)

      flux(i,j,k,imz)=flux(i,j,k,imz) -
     1   (alp*(cons(i,j,k+1,imz)*unp1-cons(i,j,k-1,imz)*unm1
     1      + (pres(i,j,k+1)-pres(i,j,k-1)))
     1  + bet*(cons(i,j,k+2,imz)*unp2-cons(i,j,k-2,imz)*unm2
     1      + (pres(i,j,k+2)-pres(i,j,k-2)))
     1  + gam*(cons(i,j,k+3,imz)*unp3-cons(i,j,k-3,imz)*unm3
     1      + (pres(i,j,k+3)-pres(i,j,k-3)))
     1  + del*(cons(i,j,k+4,imz)*unp4-cons(i,j,k-4,imz)*unm4
     1      + (pres(i,j,k+4)-pres(i,j,k-4)))) / dx(3)

      flux(i,j,k,iene)=flux(i,j,k,iene) -
     1   (alp*(cons(i,j,k+1,iene)*unp1-cons(i,j,k-1,iene)*unm1
     1      + (pres(i,j,k+1)*unp1-pres(i,j,k-1)*unm1))
     1  + bet*(cons(i,j,k+2,iene)*unp2-cons(i,j,k-2,iene)*unm2
     1      + (pres(i,j,k+2)*unp2-pres(i,j,k-2)*unm2))
     1  + gam*(cons(i,j,k+3,iene)*unp3-cons(i,j,k-3,iene)*unm3
     1      + (pres(i,j,k+3)*unp3-pres(i,j,k-3)*unm3))
     1  + del*(cons(i,j,k+4,iene)*unp4-cons(i,j,k-4,iene)*unm4
     1      + (pres(i,j,k+4)*unp4-pres(i,j,k-4)*unm4))) / dx(3)

      do nsp = 0,nspec-1

         flux(i,j,k,isp+nsp)=flux(i,j,k,isp+nsp) -
     1      (alp*(cons(i,j,k+1,isp+nsp)*unp1-cons(i,j,k-1,isp+nsp)*unm1)
     1     + bet*(cons(i,j,k+2,isp+nsp)*unp2-cons(i,j,k-2,isp+nsp)*unm2)
     1     + gam*(cons(i,j,k+3,isp+nsp)*unp3-cons(i,j,k-3,isp+nsp)*unm3)
     1     + del*(cons(i,j,k+4,isp+nsp)*unp4-cons(i,j,k-4,isp+nsp)*unm4)
     1       )/dx(3)

      enddo

      enddo
      enddo
      enddo


      return
      end

      subroutine init(lo,hi,ng,dx,nspec,cons,pres)

      use constants

      implicit none

c     inputs:  lo,hi,ng,nspec
c     outputs: cons,pres

      integer lo(3),hi(3),ng,nspec

      double precision dx(3)

      double precision pres(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,
     1             -ng+lo(3):hi(3)+ng)
      double precision cons(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,
     1             -ng+lo(3):hi(3)+ng,nspec+5)

      double precision xloc,yloc,zloc,rholoc,ploc
      double precision uvel,vvel,wvel
      double precision scale

      integer i,j,k

      scale = .02d0

      cons(:,:,:,isp:nspec+5) = 0d0

      do k=lo(3)-ng,hi(3)+ng
        zloc = dfloat(k)*dx(3)
        do j=lo(2)-ng,hi(2)+ng
          yloc = dfloat(j)*dx(2)
            do i=lo(1)-ng,hi(1)+ng
              xloc = dfloat(i)*dx(1)

              uvel = sin(xloc/scale)*sin(2.d0*yloc/scale)*
     1                 sin(3.d0*zloc/scale)
              vvel = sin(2.d0*xloc/scale)*sin(4.d0*yloc/scale)
     1                *sin(1.d0*zloc/scale)
              wvel = sin(3.d0*xloc/scale)*cos(2.d0*yloc/scale)
     1                *sin(2.d0*zloc/scale)
              rholoc = 1.d-3+1.d-5* sin(xloc/scale)*
     1                cos(2.d0*yloc/scale)*cos(3.d0*zloc/scale)
              ploc = 1.d6+ .001*sin(2.d0*xloc/scale)*
     1                cos(2.d0*yloc/scale)*sin(2.d0*zloc/scale)

              cons(i,j,k,irho) = rholoc
              cons(i,j,k,imx) = rholoc*uvel
              cons(i,j,k,imy) = rholoc*vvel
              cons(i,j,k,imz) = rholoc*wvel
              cons(i,j,k,iene) = ploc/0.4d0+
     1                rholoc*(uvel**2+vvel**2+wvel**2)/2.d0
              

              pres(i,j,k) = ploc

              cons(i,j,k,isp) = .2d0 + 0.1d0*uvel
              cons(i,j,k,isp+1) = .2d0 + 0.05d0*vvel
              cons(i,j,k,isp+2) = .2d0 + 0.03d0*wvel
              cons(i,j,k,isp+3) = 1.d0-cons(i,j,k,isp)
     1                -cons(i,j,k,isp+1)-cons(i,j,k,isp+2)

          enddo
        enddo
      enddo

      end

      program main

      implicit none

      integer, parameter :: ng    = 4
      integer, parameter :: nspec = 9

      integer, parameter :: lo(3) = (/  1,  1,  1 /)
      integer, parameter :: hi(3) = (/ 32, 32, 32 /)

      double precision, parameter :: dx(3) = (/ 1d-3, 1d-3, 1d-3 /)

      double precision, allocatable :: pres(:,:,:)
      double precision, allocatable :: cons(:,:,:,:)
      double precision, allocatable :: flux(:,:,:,:)

      double precision :: fluxmag(nspec+5)

      integer i,j,k,ns

      allocate(pres(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,
     1             -ng+lo(3):hi(3)+ng))

      allocate(cons(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,
     1             -ng+lo(3):hi(3)+ng,nspec+5))

      allocate(flux(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),nspec+5))

      call init(lo,hi,ng,dx,nspec,cons,pres)

      call hypterm(lo,hi,ng,nspec,dx,cons,pres,flux)

      fluxmag = 0.d0

      do ns = 1,nspec+5
         do k=lo(3),hi(3)
            do j=lo(2),hi(2)
               do i=lo(1),hi(1)
                  fluxmag(ns) = fluxmag(ns)+flux(i,j,k,ns)**2
               enddo
            enddo
         enddo
      enddo

      do ns = 1,nspec+5
        write(6,*)"component, fluxmag",ns,fluxmag(ns)
      enddo
 
      end program main
