#include "DDOp_F.H"
#include "ChemKinDriver_F.H"
#include "ArrayLim.H"
#include "CONSTANTS.H"

#define SDIM 2

      subroutine FORT_DIFFFORCE(lo,hi,d,DIMS(d),
     &                          X,DIMS(X),T,DIMS(T),dx,dir)
      implicit none
#include "ckdwrk.H"
      integer lo(SDIM)
      integer hi(SDIM)
      integer DIMDEC(d)
      integer DIMDEC(X)
      integer DIMDEC(T)
      REAL_T d(DIMV(d),1)
      REAL_T X(DIMV(X),1)
      REAL_T T(DIMV(X))
      REAL_T dx(SDIM)
      integer dir

      REAL_T dxinv
      integer i,j,n

      dxinv = one / dx(dir+1)

      if (dir .EQ. 0) then
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)+1
               do n=1,Nspec
                  d(i,j,n) = (X(i,j,n) - X(i-1,j,n))*dxinv
               end do
               d(i,j,Nspec+1) = (T(i,j) - T(i-1,j))*dxinv
            end do
         end do
      else
         do i=lo(1),hi(1)
            do j=lo(2),hi(2)+1
               do n=1,Nspec
                  d(i,j,n) = (X(i,j,n) - X(i,j-1,n))*dxinv
               end do
               d(i,j,Nspec+1) = (T(i,j) - T(i,j-1))*dxinv
            end do
         end do
      end if
      end

      subroutine FORT_DDC2E(lo, hi, 
     &                      DIMS(cfab), cfab,
     &                      DIMS(efab), efab, nc, dir)
      implicit none
      integer lo(SDIM), hi(SDIM), nc, dir
      integer DIMDEC(cfab)
      integer DIMDEC(efab)
      REAL_T  cfab(DIMV(cfab), nc)
      REAL_T  efab(DIMV(efab), nc)

      integer i,j,n

      if (dir .EQ. 0) then
         do n = 1,nc
            do j = lo(2)+1, hi(2)-1
               do i = lo(1)+1, hi(1)
                  efab(i,j,n) = half*(cfab(i,j,n) + cfab(i-1,j,n))
               end do
            end do
         end do
      else
         do n = 1,nc
            do i = lo(1)+1, hi(1)-1
               do j = lo(2)+1, hi(2)
                  efab(i,j,n) = half*(cfab(i,j,n) + cfab(i,j-1,n))
               end do
            end do
         end do
      end if
      end

      subroutine FORT_FLUX(lo, hi, fluxH, DIMS(fluxH),
     &                     fluxY, DIMS(fluxY), d, DIMS(d),
     &                     Ye, DIMS(Ye), Te, DIMS(Te), He, DIMS(He), dir)
      implicit none
#include "ckdwrk.H"
      integer lo(SDIM), hi(SDIM), dir
      integer DIMDEC(fluxH)
      integer DIMDEC(fluxY)
      integer DIMDEC(d)
      integer DIMDEC(Ye)
      integer DIMDEC(Te)
      integer DIMDEC(He)
      integer DIMDEC(Y)
      integer DIMDEC(T)
      REAL_T  fluxH(DIMV(fluxH))
      REAL_T  fluxY(DIMV(fluxY),1)
      REAL_T      d(DIMV(d)    ,1)
      REAL_T     Ye(DIMV(Ye)   ,1)
      REAL_T     Te(DIMV(Te))
      REAL_T     He(DIMV(He)   ,1)

      integer i,j,lo1,lo2,hi1,hi2,K,L
      REAL_T Yt(maxspec), WW, CPMS(maxspec), X(maxspec)
      REAL_T rhoD, rhoD1(maxspec*maxspec), rhoTD(maxspec), PTC
      REAL_T sum_rhoYV, sum_rhoHYV, sum_rhoThetaD, P1atm, RU, RUC

      rhoD(K,L) = rhoDt( (L-1)*Nspec + K )
      
      lo1 = lo(1)
      lo2 = lo(2)
      if (dir .EQ. 0) then
         hi1 = hi(1) + 1
         hi2 = hi(2)
      else 
         hi1 = hi(1)
         hi2 = hi(2) + 1
      end if
      
      CALL CKRP(IWRK(ckbi), RWRK(ckbr), RU, RUC, P1atm)
      
      do j = lo2, hi2
         do i = lo1, hi1

            do K = 1,Nspec
               Yt(K) = Ye(i,j,K)
            end do

            CALL CKCPMS(Te(i,j),IWRK(ckbi),RWRK(ckbr),CPMS)
            CALL CKYTX(Yt,IWRK(ckbi),RWRK(ckbr),X)
            CALL CKMMWY(Yt, IWRK(ckbr), RWRK(ckbr), WW)
            CALL EGSPAR(Te(i,j),X,Yt,CPMS,RWRK(egbr),IWRK(egbi))
            CALL EGSLTDR5(Te(i,j),Yt,WW,RWRK(egbr),IWRK(egbi),
     &                    PTC,rhoTD,rhoDt)
            
c     Note: scalings get cgs quantities from EGLib to MKS (assume He,d come in MKS)
            sum_rhoYV     = zero
            sum_rhoHYV    = zero
            sum_rhoThetaD = zero

            do K= 1,Nspec

               sum_rhoYV = - Ye(i,j,K)*rhoTD(K)*d(i,j,Nspec+1)/Te(i,j)
               do L= 1,Nspec
                  sum_rhoYV = sum_rhoYV - rhoD(K,L)*d(i,j,L)
               end do
               fluxY(i,j,K) = sum_rhoYV*1.d-1

               sum_rhoHYV = sum_rhoHYV + He(i,j,K)*fluxY(i,j,K)
               sum_rhoThetaD = sum_rhoThetaD + rhoTD(K)*d(i,j,K)

            end do

            fluxH(i,j) = sum_rhoHYV - (PTC*1.d-5)*d(i,j,Nspec+1)
     &           - (RU*Te(i,j)/WW)*(sum_rhoThetaD*1.d-1)

         end do
      end do
      
      end
      
      subroutine FORT_INCRDIV(lo, hi, div, DIMS(div), V, DIMS(V),
     &                        dx, dir, nc)
      implicit none
#include "ckdwrk.H"
      integer lo(SDIM), hi(SDIM), dir, nc
      integer DIMDEC(div)
      integer DIMDEC(V)
      REAL_T div(DIMV(div) ,1)
      REAL_T   V(DIMV(V)   ,1)
      REAL_T dx(SDIM)

      integer i,j,n

c     Note: We're incrementing with div(dx.Area.Flux), so the dx's
c     cancel, and we assume V is extensive (ie. Flux.Area)
      if (dir .EQ. 0) then
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               do n= 1,nc
                  div(i,j,n) = div(i,j,n)
     &                 + (V(i+1,j,n) - V(i,j,n))
               end do
            end do               
         end do
      else 
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               do n= 1,nc
                  div(i,j,n) = div(i,j,n)
     &                 + (V(i,j+1,n) - V(i,j,n))
               end do
            end do               
         end do
      end if
      end


      subroutine FORT_CPSCALE(lo,hi,lambda,DIMS(lambda),
     &                        Y,DIMS(Y),T,DIMS(T))
      implicit none
#include "ckdwrk.H"
      integer lo(SDIM)
      integer hi(SDIM)
      integer DIMDEC(lambda)
      integer DIMDEC(Y)
      integer DIMDEC(T)
      REAL_T lambda(DIMV(lambda))
      REAL_T Y(DIMV(Y),1)
      REAL_T T(DIMV(T))

      REAL_T YT(maxspec),CPMIX
      integer i,j,K
      
      do j = lo(2), hi(2)
         do i = lo(1), hi(1)

            do K = 1,Nspec
               Yt(K) = Y(i,j,K)
            end do

            CALL CKCPBS(T(i,j),Yt,IWRK(ckbi),RWRK(ckbr),CPMIX)
            lambda(i,j) = lambda(i,j)/(CPMIX*1.d-4)
            
         enddo
      enddo
      end


      subroutine FORT_THERM(lo, hi, lambda, DIMS(lambda),
     &                     Ye, DIMS(Ye), Te, DIMS(Te), dir)
      implicit none
#include "ckdwrk.H"
      integer lo(SDIM), hi(SDIM), dir
      integer DIMDEC(Ye)
      integer DIMDEC(Te)
      integer DIMDEC(Y)
      integer DIMDEC(T)
      integer DIMDEC(lambda)
      REAL_T     Ye(DIMV(Ye)   ,1)
      REAL_T     Te(DIMV(Te))
      REAL_T     lambda(DIMV(lambda)   )

      integer i,j,lo1,lo2,hi1,hi2,K,L
      REAL_T Yt(maxspec), WW, CPMS(maxspec), X(maxspec), TC

c     Note: scaling is to get D1 units of g/(s.cm^2) into MKS
      
      lo1 = lo(1)
      lo2 = lo(2)
      if (dir .EQ. 0) then
         hi1 = hi(1) + 1
         hi2 = hi(2)
         
         do j = lo2, hi2
            do i = lo1, hi1
               
               do K = 1,Nspec
                  Yt(K) = Ye(i,j,K)
               end do
               
               CALL CKCPMS(Te(i,j),IWRK(ckbi),RWRK(ckbr),CPMS)
               CALL CKYTX(Yt,IWRK(ckbi),RWRK(ckbr),X)
               CALL EGSPAR(Te(i,j),X,Yt,CPMS,RWRK(egbr),IWRK(egbi))
               CALL EGSL4(Te(i,j),Yt,RWRK(egbr),IWRK(egbi),TC)
               
               if(i-1.ge.lo(1))then
                  lambda(i-1,j) = lambda(i-1,j)+TC*1.d-5
               endif
               if(i.le.hi(1))then
                  lambda(i,j) = lambda(i,j)+TC*1.d-5
               endif

            end do
         end do

      else 

         hi1 = hi(1)
         hi2 = hi(2) + 1
         do j = lo2, hi2
            do i = lo1, hi1

               do K = 1,Nspec
                  Yt(K) = Ye(i,j,K)
               end do

               CALL CKCPMS(Te(i,j),IWRK(ckbi),RWRK(ckbr),CPMS)
               CALL CKYTX(Yt,IWRK(ckbi),RWRK(ckbr),X)
               CALL EGSPAR(Te(i,j),X,Yt,CPMS,RWRK(egbr),IWRK(egbi))
               CALL EGSL4(Te(i,j),Yt,RWRK(egbr),IWRK(egbi),TC)
            
               if(j-1.ge.lo(2))then
                  lambda(i,j-1) = lambda(i,j-1)+TC*1.d-5
               endif
               if(j.le.hi(2))then
                  lambda(i,j) = lambda(i,j)+TC*1.d-5
               endif

            end do
         end do
      end if
      
      end

      subroutine FORT_CRSNCCBND(lo, hi, fine, DIMS(fine),
     &                          crse, DIMS(crse), nc, face)
      integer lo(SDIM), hi(SDIM), nc
      integer DIMDEC(fine)
      integer DIMDEC(crse)
      REAL_T fine(DIMV(fine),nc)
      REAL_T crse(DIMV(crse),nc)
      integer face
      
      integer i,j,n, clo1, chi1, clo2, chi2, offset, dir
      
      if ((face.EQ.0).OR.(face.EQ.2)) then
         dir = 0
      else
         dir = 1
      end if

      if (face .LT. SDIM) then
         offset = 1
      else
         offset = 0
      end if
      
      if (dir .eq. 0) then
         clo1 = (lo(1)-offset)/2
         clo2 = lo(2)/2
         chi1 = (hi(1)-offset)/2
         chi2 = (hi(2)+1)/2 - 1
         do j=clo2,chi2
            do i=clo1,chi1
               do n=1,nc
                  crse(i,j,n) =
     &                 half*(fine(2*i+offset,2*j,n)+fine(2*i+offset,2*j+1,n))
               end do
            end do
         end do
      else
         clo1 = lo(1)/2
         clo2 = (lo(2)-offset)/2
         chi1 = (hi(1)+1)/2 - 1
         chi2 = (hi(2)-offset)/2
         do j=clo2,chi2
            do i=clo1,chi1
               do n=1,nc
                  crse(i,j,n) =
     &                 half*(fine(2*i,2*j+offset,n)+fine(2*i+1,2*j+offset,n))
               end do
            end do
         end do
      end if
      end

      subroutine FORT_DDCCAVG(c, DIMS(c), f, DIMS(f), lo, hi, nc)
      integer nc
      integer DIMDEC(f)
      integer DIMDEC(c)
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      REAL_T f(DIMV(f),nc)
      REAL_T c(DIMV(c),nc)
c
      integer i
      integer j
      integer n
      REAL_T denom
      parameter(denom=fourth)
c
      do n = 1, nc
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               c(i,j,n) = (
     $              f(2*i+1,2*j+1,n) + f(2*i  ,2*j+1,n)
     $              + f(2*i+1,2*j,n ) + f(2*i  ,2*j ,n))*denom
            end do
         end do
      end do
c
      end

      subroutine FORT_DDCCINT(f, DIMS(f), c, DIMS(c), lo, hi, nc)
      integer nc
      integer DIMDEC(f)
      integer DIMDEC(c)
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      REAL_T f(DIMV(f),nc)
      REAL_T c(DIMV(c),nc)
c
      integer i
      integer j
      integer n
c
      do n = 1, nc
         do j = lo(2), hi(2)
cdir$   ivdep
            do i = lo(1), hi(1)
               f(2*i+1,2*j+1,n) = c(i,j,n) + f(2*i+1,2*j+1,n)
               f(2*i  ,2*j+1,n) = c(i,j,n) + f(2*i  ,2*j+1,n)
               f(2*i+1,2*j  ,n) = c(i,j,n) + f(2*i+1,2*j  ,n)
               f(2*i  ,2*j  ,n) = c(i,j,n) + f(2*i  ,2*j  ,n)
            end do
         end do
      end do
c
      end
      
