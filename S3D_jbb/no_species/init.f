
      subroutine init(lo,hi,ng,dx,cons,pres)

      use constants

      implicit none


      integer lo(3),hi(3),ng

      double precision dx(3)

      double precision pres(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,
     1             -ng+lo(3):hi(3)+ng)
      double precision cons(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,
     1             -ng+lo(3):hi(3)+ng,5)

      double precision xloc,yloc,zloc,rholoc,ploc
      double precision uvel,vvel,wvel
      double precision scale

      integer i,j,k

      scale = .02d0

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

          enddo
        enddo
      enddo

      end

