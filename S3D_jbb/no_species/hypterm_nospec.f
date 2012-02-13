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

      integer, parameter :: qrho = 1
      integer, parameter :: qu  = 2
      integer, parameter :: qv  = 3
      integer, parameter :: qw  = 4
      integer, parameter :: qpres = 5
      integer, parameter :: qT = 5

      end module constants

      subroutine hypterm(lo,hi,ng,dx,cons,q,flux)

      use constants

      implicit none

c     inputs:  lo,hi,ng,cons,q
c     cons -- rho, rho u , rho v, rho w, rho E
c     cons -- rho, u ,  v,  w, p, T
c     outputs: flux

      integer lo(3),hi(3),ng

      double precision dx(3)

      double precision cons(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,
     1             -ng+lo(3):hi(3)+ng,5)
      double precision q(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,
     1             -ng+lo(3):hi(3)+ng,6)
      double precision flux(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),5)

      integer i,j,k,nsp

      double precision unp1,unp2,unp3,unp4,unm1,unm2,unm3,unm4

      do k=lo(3),hi(3)
      do j=lo(2),hi(2)
      do i=lo(1),hi(1)

      unp1 = q(i+1,j,k,qu)
      unp2 = q(i+2,j,k,qu)
      unp3 = q(i+3,j,k,qu)
      unp4 = q(i+4,j,k,qu)

      unm1 = q(i-1,j,k,qu)
      unm2 = q(i-2,j,k,qu)
      unm3 = q(i-3,j,k,qu)
      unm4 = q(i-4,j,k,qu)

      flux(i,j,k,irho)= -(alp*(cons(i+1,j,k,imx)-cons(i-1,j,k,imx))
     1  + bet*(cons(i+2,j,k,imx)-cons(i-2,j,k,imx))
     1  + gam*(cons(i+3,j,k,imx)-cons(i-3,j,k,imx))
     1  + del*(cons(i+4,j,k,imx)-cons(i-4,j,k,imx)))/dx(1)

      flux(i,j,k,imx)= -(alp*(cons(i+1,j,k,imx)*unp1-cons(i-1,j,k,imx)*unm1
     1      + (q(i+1,j,k,qpres)-q(i-1,j,k,qpres)))
     1  + bet*(cons(i+2,j,k,imx)*unp2-cons(i-2,j,k,imx)*unm2
     1      + (q(i+2,j,k,qpres)-q(i-2,j,k,qpres)))
     1  + gam*(cons(i+3,j,k,imx)*unp3-cons(i-3,j,k,imx)*unm3
     1      + (q(i+3,j,k,qpres)-q(i-3,j,k,qpres)))
     1  + del*(cons(i+4,j,k,imx)*unp4-cons(i-4,j,k,imx)*unm4
     1      + (q(i+4,j,k,qpres)-q(i-4,j,k,qpres)))) / dx(1)

      flux(i,j,k,imy)= -
     1   (alp*(cons(i+1,j,k,imy)*unp1-cons(i-1,j,k,imy)*unm1)
     1  + bet*(cons(i+2,j,k,imy)*unp2-cons(i-2,j,k,imy)*unm2)
     1  + gam*(cons(i+3,j,k,imy)*unp3-cons(i-3,j,k,imy)*unm3)
     1  + del*(cons(i+4,j,k,imy)*unp4-cons(i-4,j,k,imy)*unm4))/dx(1)

      flux(i,j,k,imz)= -
     1   (alp*(cons(i+1,j,k,imz)*unp1-cons(i-1,j,k,imz)*unm1)
     1  + bet*(cons(i+2,j,k,imz)*unp2-cons(i-2,j,k,imz)*unm2)
     1  + gam*(cons(i+3,j,k,imz)*unp3-cons(i-3,j,k,imz)*unm3)
     1  + del*(cons(i+4,j,k,imz)*unp4-cons(i-4,j,k,imz)*unm4))/dx(1)

      flux(i,j,k,iene)= -
     1   (alp*(cons(i+1,j,k,iene)*unp1-cons(i-1,j,k,iene)*unm1
     1      + (q(i+1,j,k,qpres)*unp1-q(i-1,j,k,qpres)*unm1))
     1  + bet*(cons(i+2,j,k,iene)*unp2-cons(i-2,j,k,iene)*unm2
     1      + (q(i+2,j,k,qpres)*unp2-q(i-2,j,k,qpres)*unm2))
     1  + gam*(cons(i+3,j,k,iene)*unp3-cons(i-3,j,k,iene)*unm3
     1      + (q(i+3,j,k,qpres)*unp3-q(i-3,j,k,qpres)*unm3))
     1  + del*(cons(i+4,j,k,iene)*unp4-cons(i-4,j,k,iene)*unm4
     1      + (q(i+4,j,k,qpres)*unp4-q(i-4,j,k,qpres)*unm4))) / dx(1)


      enddo
      enddo
      enddo

      do k=lo(3),hi(3)
      do j=lo(2),hi(2)
      do i=lo(1),hi(1)

      unp1 = q(i,j+1,k,qv)
      unp2 = q(i,j+2,k,qv)
      unp3 = q(i,j+3,k,qv)
      unp4 = q(i,j+4,k,qv)

      unm1 = q(i,j-1,k,qv)
      unm2 = q(i,j-2,k,qv)
      unm3 = q(i,j-3,k,qv)
      unm4 = q(i,j-4,k,qv)

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
     1      + (q(i,j+1,k,qpres)-q(i,j-1,k,qpres)))
     1  + bet*(cons(i,j+2,k,imy)*unp2-cons(i,j-2,k,imy)*unm2
     1      + (q(i,j+2,k,qpres)-q(i,j-2,k,qpres)))
     1  + gam*(cons(i,j+3,k,imy)*unp3-cons(i,j-3,k,imy)*unm3
     1      + (q(i,j+3,k,qpres)-q(i,j-3,k,qpres)))
     1  + del*(cons(i,j+4,k,imy)*unp4-cons(i,j-4,k,imy)*unm4
     1      + (q(i,j+4,k,qpres)-q(i,j-4,k,qpres)))) / dx(2)

      flux(i,j,k,imz)=flux(i,j,k,imz) -
     1   (alp*(cons(i,j+1,k,imz)*unp1-cons(i,j-1,k,imz)*unm1)
     1  + bet*(cons(i,j+2,k,imz)*unp2-cons(i,j-2,k,imz)*unm2)
     1  + gam*(cons(i,j+3,k,imz)*unp3-cons(i,j-3,k,imz)*unm3)
     1  + del*(cons(i,j+4,k,imz)*unp4-cons(i,j-4,k,imz)*unm4))/dx(2)

      flux(i,j,k,iene)=flux(i,j,k,iene) -
     1   (alp*(cons(i,j+1,k,iene)*unp1-cons(i,j-1,k,iene)*unm1
     1      + (q(i,j+1,k,qpres)*unp1-q(i,j-1,k,qpres)*unm1))
     1  + bet*(cons(i,j+2,k,iene)*unp2-cons(i,j-2,k,iene)*unm2
     1      + (q(i,j+2,k,qpres)*unp2-q(i,j-2,k,qpres)*unm2))
     1  + gam*(cons(i,j+3,k,iene)*unp3-cons(i,j-3,k,iene)*unm3
     1      + (q(i,j+3,k,qpres)*unp3-q(i,j-3,k,qpres)*unm3))
     1  + del*(cons(i,j+4,k,iene)*unp4-cons(i,j-4,k,iene)*unm4
     1      + (q(i,j+4,k,qpres)*unp4-q(i,j-4,k,qpres)*unm4))) / dx(2)


      enddo
      enddo
      enddo


      do k=lo(3),hi(3)
      do j=lo(2),hi(2)
      do i=lo(1),hi(1)

      unp1 = q(i,j,k+1,qw)
      unp2 = q(i,j,k+2,qw)
      unp3 = q(i,j,k+3,qw)
      unp4 = q(i,j,k+4,qw)

      unm1 = q(i,j,k-1,qw)
      unm2 = q(i,j,k-2,qw)
      unm3 = q(i,j,k-3,qw)
      unm4 = q(i,j,k-4,qw)

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
     1      + (q(i,j,k+1,qpres)-q(i,j,k-1,qpres)))
     1  + bet*(cons(i,j,k+2,imz)*unp2-cons(i,j,k-2,imz)*unm2
     1      + (q(i,j,k+2,qpres)-q(i,j,k-2,qpres)))
     1  + gam*(cons(i,j,k+3,imz)*unp3-cons(i,j,k-3,imz)*unm3
     1      + (q(i,j,k+3,qpres)-q(i,j,k-3,qpres)))
     1  + del*(cons(i,j,k+4,imz)*unp4-cons(i,j,k-4,imz)*unm4
     1      + (q(i,j,k+4,qpres)-q(i,j,k-4,qpres)))) / dx(3)

      flux(i,j,k,iene)=flux(i,j,k,iene) -
     1   (alp*(cons(i,j,k+1,iene)*unp1-cons(i,j,k-1,iene)*unm1
     1      + (q(i,j,k+1,qpres)*unp1-q(i,j,k-1,qpres)*unm1))
     1  + bet*(cons(i,j,k+2,iene)*unp2-cons(i,j,k-2,iene)*unm2
     1      + (q(i,j,k+2,qpres)*unp2-q(i,j,k-2,qpres)*unm2))
     1  + gam*(cons(i,j,k+3,iene)*unp3-cons(i,j,k-3,iene)*unm3
     1      + (q(i,j,k+3,qpres)*unp3-q(i,j,k-3,qpres)*unm3))
     1  + del*(cons(i,j,k+4,iene)*unp4-cons(i,j,k-4,iene)*unm4
     1      + (q(i,j,k+4,qpres)*unp4-q(i,j,k-4,qpres)*unm4))) / dx(3)


      enddo
      enddo
      enddo


      return
      end
