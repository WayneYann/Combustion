      module diff_constants

      double precision, parameter :: center = -205.d0/72.d0
      double precision, parameter :: off1 = 8.d0 / 5.d0
      double precision, parameter :: off2 = -.2d0
      double precision, parameter :: off3 = 8.d0/315.d0
      double precision, parameter :: off4 = -1.d0/560.d0

      double precision, parameter :: alp = .8d0
      double precision, parameter :: bet = -.2d0
      double precision, parameter :: gam = 4.d0/105.d0
      double precision, parameter :: del = -1.d0/280.d0

      integer, parameter :: irho = 1
      integer, parameter :: imx  = 2
      integer, parameter :: imy  = 3
      integer, parameter :: imz  = 4
      integer, parameter :: iene = 5

      integer, parameter :: qrho = 1
      integer, parameter :: qu  = 2
      integer, parameter :: qv  = 3
      integer, parameter :: qw  = 4
      integer, parameter :: qpres = 5
      integer, parameter :: qT = 5

      end module constants

      end module diff_constants

      subroutine diffterm(lo,hi,ng,dx,cons,q,difflux,eta,alam)

      use diff_constants

      implicit none


      integer lo(3),hi(3),ng

      double precision dx(3)

      double precision eta, alam

      double precision cons(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,
     1             -ng+lo(3):hi(3)+ng,5)
      double precision q(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,
     1             -ng+lo(3):hi(3)+ng,6)
      double precision difflux(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),5)

      double precision, allocatable :: ux(:,:,:)
      double precision, allocatable :: uy(:,:,:)
      double precision, allocatable :: uz(:,:,:)
      double precision, allocatable :: vx(:,:,:)
      double precision, allocatable :: vy(:,:,:)
      double precision, allocatable :: vz(:,:,:)
      double precision, allocatable :: wx(:,:,:)
      double precision, allocatable :: wy(:,:,:)
      double precision, allocatable :: wz(:,:,:)
 
      double precision tauxx,tauyy,tauzz,tauxy,tauxz,tauyz
      double precision divu, uxx,uyy,uzz,vxx,vyy,vzz,wxx,wyy,wzz
      double precision mechwork
      double precision uxy,uxz,v

      integer i,j,k

      allocate(ux(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,
     1             -ng+lo(3):hi(3)+ng))
      allocate(uy(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,
     1             -ng+lo(3):hi(3)+ng))
      allocate(uz(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,
     1             -ng+lo(3):hi(3)+ng))
      allocate(vx(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,
     1             -ng+lo(3):hi(3)+ng))
      allocate(vy(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,
     1             -ng+lo(3):hi(3)+ng))
      allocate(vz(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,
     1             -ng+lo(3):hi(3)+ng))
      allocate(wx(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,
     1             -ng+lo(3):hi(3)+ng))
      allocate(wy(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,
     1             -ng+lo(3):hi(3)+ng))
      allocate(wz(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,
     1             -ng+lo(3):hi(3)+ng))


      do k=lo(3)-ng,hi(3)+ng
      do j=lo(2)-ng,hi(2)+ng
      do i=lo(1),hi(1)
 
      ux(i,j,k)= (alp*(q(i+1,j,k,qu)-q(i-1,j,k,qu))
     1  + bet*(q(i+2,j,k,qu)-q(i-2,j,k,qu))
     1  + gam*(q(i+3,j,k,qu)-q(i-3,j,k,qu))
     1  + del*(q(i+4,j,k,qu)-q(i-4,j,k,qu)))/dx(1)
 
      vx(i,j,k)= (alp*(q(i+1,j,k,qv)-q(i-1,j,k,qv))
     1  + bet*(q(i+2,j,k,qv)-q(i-2,j,k,qv))
     1  + gam*(q(i+3,j,k,qv)-q(i-3,j,k,qv))
     1  + del*(q(i+4,j,k,qv)-q(i-4,j,k,qv)))/dx(1)
 
      wx(i,j,k)= (alp*(q(i+1,j,k,qw)-q(i-1,j,k,qw))
     1  + bet*(q(i+2,j,k,qw)-q(i-2,j,k,qw))
     1  + gam*(q(i+3,j,k,qw)-q(i-3,j,k,qw))
     1  + del*(q(i+4,j,k,qw)-q(i-4,j,k,qw)))/dx(1)


      enddo
      enddo
      enddo

      do k=lo(3)-ng,hi(3)+ng
      do j=lo(2)   ,hi(2)   
      do i=lo(1)-ng,hi(1)+ng
 
      vy(i,j,k)= (alp*(q(i,j+1,k,qv)-q(i,j-1,k,qv))
     1  + bet*(q(i,j+2,k,qv)-q(i,j-2,k,qv))
     1  + gam*(q(i,j+3,k,qv)-q(i,j-3,k,qv))
     1  + del*(q(i,j+4,k,qv)-q(i,j-4,k,qv)))/dx(2)
 
      uy(i,j,k)= (alp*(q(i,j+1,k,qu)-q(i,j-1,k,qu))
     1  + bet*(q(i,j+2,k,qu)-q(i,j-2,k,qu))
     1  + gam*(q(i,j+3,k,qu)-q(i,j-3,k,qu))
     1  + del*(q(i,j+4,k,qu)-q(i,j-4,k,qu)))/dx(2)
 
      wy(i,j,k)= (alp*(q(i,j+1,k,qw)-q(i,j-1,k,qw))
     1  + bet*(q(i,j+2,k,qw)-q(i,j-2,k,qw))
     1  + gam*(q(i,j+3,k,qw)-q(i,j-3,k,qw))
     1  + del*(q(i,j+4,k,qw)-q(i,j-4,k,qw)))/dx(2)


      enddo
      enddo
      enddo

      do k=lo(3),hi(3)
      do j=lo(2)-ng,hi(2)+ng
      do i=lo(1)-ng,hi(1)+ng
 
      wz(i,j,k)= (alp*(q(i,j,k+1,qw)-q(i,j,k-1,qw))
     1  + bet*(q(i,j,k+2,qw)-q(i,j,k-2,qw))
     1  + gam*(q(i,j,k+3,qw)-q(i,j,k-3,qw))
     1  + del*(q(i,j,k+4,qw)-q(i,j,k-4,qw)))/dx(3)
 
      uz(i,j,k)= (alp*(q(i,j,k+1,qu)-q(i,j,k-1,qu))
     1  + bet*(q(i,j,k+2,qu)-q(i,j,k-2,qu))
     1  + gam*(q(i,j,k+3,qu)-q(i,j,k-3,qu))
     1  + del*(q(i,j,k+4,qu)-q(i,j,k-4,qu)))/dx(3)
 
      vz(i,j,k)= (alp*(q(i,j,k+1,qv)-q(i,j,k-1,qv))
     1  + bet*(q(i,j,k+2,qv)-q(i,j,k-2,qv))
     1  + gam*(q(i,j,k+3,qv)-q(i,j,k-3,qv))
     1  + del*(q(i,j,k+4,qv)-q(i,j,k-4,qv)))/dx(3)


      enddo
      enddo
      enddo


      do k=lo(3),hi(3)
      do j=lo(2),hi(2)
      do i=lo(1),hi(1)


      uxx = (center*q(i,j,k,qu)
     1  + off1*(q(i+1,j,k,qu)+q(i-1,j,k,qu))
     1  + off2*(q(i+2,j,k,qu)+q(i-2,j,k,qu))
     1  + off3*(q(i+3,j,k,qu)+q(i-3,j,k,qu))
     1  + off4*(q(i+4,j,k,qu)+q(i-4,j,k,qu)))/dx(1)**2

      uyy = (center*q(i,j,k,qu)
     1  + off1*(q(i,j+1,k,qu)+q(i,j-1,k,qu))
     1  + off2*(q(i,j+2,k,qu)+q(i,j-2,k,qu))
     1  + off3*(q(i,j+3,k,qu)+q(i,j-3,k,qu))
     1  + off4*(q(i,j+4,k,qu)+q(i,j-4,k,qu)))/dx(2)**2

      uzz = (center*q(i,j,k,qu)
     1  + off1*(q(i,j,k+1,qu)+q(i,j,k-1,qu))
     1  + off2*(q(i,j,k+2,qu)+q(i,j,k-2,qu))
     1  + off3*(q(i,j,k+3,qu)+q(i,j,k-3,qu))
     1  + off4*(q(i,j,k+4,qu)+q(i,j,k-4,qu)))/dx(3)**2


      vxx = (center*q(i,j,k,qv)
     1  + off1*(q(i+1,j,k,qv)+q(i-1,j,k,qv))
     1  + off2*(q(i+2,j,k,qv)+q(i-2,j,k,qv))
     1  + off3*(q(i+3,j,k,qv)+q(i-3,j,k,qv))
     1  + off4*(q(i+4,j,k,qv)+q(i-4,j,k,qv)))/dx(1)**2

      vyy = (center*q(i,j,k,qv)
     1  + off1*(q(i,j+1,k,qv)+q(i,j-1,k,qv))
     1  + off2*(q(i,j+2,k,qv)+q(i,j-2,k,qv))
     1  + off3*(q(i,j+3,k,qv)+q(i,j-3,k,qv))
     1  + off4*(q(i,j+4,k,qv)+q(i,j-4,k,qv)))/dx(2)**2

      vzz = (center*q(i,j,k,qv)
     1  + off1*(q(i,j,k+1,qv)+q(i,j,k-1,qv))
     1  + off2*(q(i,j,k+2,qv)+q(i,j,k-2,qv))
     1  + off3*(q(i,j,k+3,qv)+q(i,j,k-3,qv))
     1  + off4*(q(i,j,k+4,qv)+q(i,j,k-4,qv)))/dx(3)**2



      wxx = (center*q(i,j,k,qw)
     1  + off1*(q(i+1,j,k,qw)+q(i-1,j,k,qw))
     1  + off2*(q(i+2,j,k,qw)+q(i-2,j,k,qw))
     1  + off3*(q(i+3,j,k,qw)+q(i-3,j,k,qw))
     1  + off4*(q(i+4,j,k,qw)+q(i-4,j,k,qw)))/dx(1)**2

      wyy = (center*q(i,j,k,qw)
     1  + off1*(q(i,j+1,k,qw)+q(i,j-1,k,qw))
     1  + off2*(q(i,j+2,k,qw)+q(i,j-2,k,qw))
     1  + off3*(q(i,j+3,k,qw)+q(i,j-3,k,qw))
     1  + off4*(q(i,j+4,k,qw)+q(i,j-4,k,qw)))/dx(2)**2

      wzz = (center*q(i,j,k,qw)
     1  + off1*(q(i,j,k+1,qw)+q(i,j,k-1,qw))
     1  + off2*(q(i,j,k+2,qw)+q(i,j,k-2,qw))
     1  + off3*(q(i,j,k+3,qw)+q(i,j,k-3,qw))
     1  + off4*(q(i,j,k+4,qw)+q(i,j,k-4,qw)))/dx(3)**2


      vyx= (alp*(vy(i+1,j,k)-vy(i-1,j,k))
     1  + bet*(vy(i+2,j,k)-vy(i-2,j,k))
     1  + gam*(vy(i+3,j,k)-vy(i-3,j,k))
     1  + del*(vy(i+4,j,k)-vy(i-4,j,k)))/dx(1)

      wzx= (alp*(wz(i+1,j,k)-wz(i-1,j,k))
     1  + bet*(wz(i+2,j,k)-wz(i-2,j,k))
     1  + gam*(wz(i+3,j,k)-wz(i-3,j,k))
     1  + del*(wz(i+4,j,k)-wz(i-4,j,k)))/dx(1)


      uxy = (alp*(ux(i,j+1,k)-ux(i,j-1,k))
     1  + bet*(ux(i,j+2,k)-ux(i,j-2,k))
     1  + gam*(ux(i,j+3,k)-ux(i,j-3,k))
     1  + del*(ux(i,j+4,k)-ux(i,j-4,k)))/dx(2)

      wzy = (alp*(wz(i,j+1,k)-wz(i,j-1,k))
     1  + bet*(wz(i,j+2,k)-wz(i,j-2,k))
     1  + gam*(wz(i,j+3,k)-wz(i,j-3,k))
     1  + del*(wz(i,j+4,k)-wz(i,j-4,k)))/dx(2)

      uxz= (alp*(ux(i,j,k+1)-ux(i,j,k-1))
     1  + bet*(ux(i,j,k+2)-ux(i,j,k-2))
     1  + gam*(ux(i,j,k+3)-ux(i,j,k-3))
     1  + del*(ux(i,j,k+4)-ux(i,j,k-4)))/dx(3)

      vyz= (alp*(vy(i,j,k+1)-vy(i,j,k-1))
     1  + bet*(vy(i,j,k+2)-vy(i,j,k-2))
     1  + gam*(vy(i,j,k+3)-vy(i,j,k-3))
     1  + del*(vy(i,j,k+4)-vy(i,j,k-4)))/dx(3)

      difflux(i,j,k,imx) = eta * (4.d0*uxx/3.d0 + uyy + uzz
     1     (vyx+wzx)/3.d0)

      difflux(i,j,k,imy) = eta * (vxx + 4.d0*vyy/3.d0 + vzz
     1     (uxy+wzy)/3.d0)

      difflux(i,j,k,imz) = eta * (wxx + wyy + 4.d0*wzz/3.d0
     1     (uxz+vyz)/3.d0)

      txx = (center*Temp(i,j,k)
     1  + off1*(Temp(i+1,j,k)+Temp(i-1,j,k))
     1  + off2*(Temp(i+2,j,k)+Temp(i-2,j,k))
     1  + off3*(Temp(i+3,j,k)+Temp(i-3,j,k))
     1  + off4*(Temp(i+4,j,k)+Temp(i-4,j,k)))/dx(1)**2

      tyy = (center*Temp(i,j,k)
     1  + off1*(Temp(i,j+1,k)+Temp(i,j-1,k))
     1  + off2*(Temp(i,j+2,k)+Temp(i,j-2,k))
     1  + off3*(Temp(i,j+3,k)+Temp(i,j-3,k))
     1  + off4*(Temp(i,j+4,k)+Temp(i,j-4,k)))/dx(2)**2

      tzz = (center*Temp(i,j,k)
     1  + off1*(Temp(i,j,k+1)+Temp(i,j,k-1))
     1  + off2*(Temp(i,j,k+2)+Temp(i,j,k-2))
     1  + off3*(Temp(i,j,k+3)+Temp(i,j,k-3))
     1  + off4*(Temp(i,j,k+4)+Temp(i,j,k-4)))/dx(3)**2

      
      difflux(i,j,k,iene)= alam*(txx+tyy+tzz)

c    differentiate out \nabla \cdot \tau u
c    this is rather horrific

     divu = ux(i,j,k)+vy(i,j,k)+wz(i,j,k)
     tauxx = 2.d0*ux(i,j,k) - 2.d0*divu/3.d0
     tauyy = 2.d0*vy(i,j,k) - 2.d0*divu/3.d0
     tauzx = 2.d0*wz(i,j,k) - 2.d0*divu/3.d0
     tauxy = uy(i,j,k)+vx(i,j,k)
     tauxz = uz(i,j,k)+wx(i,j,k)
     tauyx = vz(i,j,k)+wy(i,j,k)

      mechwork = tauxx *ux(i,j,k) + 
     1   tauyy*vy(i,j,k)+tauzz*wz(i,j,k)+
     1   tauxy*2+tauxz**2+tauyz**2

      mechwork = eta*mechwork + difflux(i,j,k,imx)*q(i,j,k,qu)
     1   + difflux(i,j,k,imy)*q(i,j,k,qv)
     1   + difflux(i,j,k,imz)*q(i,j,k,qw)


      enddo
      enddo
      enddo

      call deallocate(ux)
      call deallocate(uy)
      call deallocate(uz)
      call deallocate(vx)
      call deallocate(vy)
      call deallocate(vz)
      call deallocate(wx)
      call deallocate(wy)
      call deallocate(wz)

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
