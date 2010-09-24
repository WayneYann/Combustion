#include "DDOp_F.H"
#include "ChemDriver_F.H"
#include "ArrayLim.H"
#include "CONSTANTS.H"

#define SDIM 2

      subroutine FORT_DIFFFORCE(lo,hi,d,DIMS(d),
     &                          X,DIMS(X),T,DIMS(T),dx,dir)
      implicit none
#include "cdwrk.H"
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
               enddo
               d(i,j,Nspec+1) = (T(i,j) - T(i-1,j))*dxinv
            enddo
         enddo
      else
         do i=lo(1),hi(1)
            do j=lo(2),hi(2)+1
               do n=1,Nspec
                  d(i,j,n) = (X(i,j,n) - X(i,j-1,n))*dxinv
               enddo
               d(i,j,Nspec+1) = (T(i,j) - T(i,j-1))*dxinv
            enddo
         enddo
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

c     Here, the cc data is assumed to be computable over the box (lo,hi) 
c     surroundingNodes of grow((lo,hi),-1) in the dir direction.
      if (dir .EQ. 0) then
         do n = 1,nc
            do j = lo(2)+1, hi(2)-1
               do i = lo(1)+1, hi(1)
                  efab(i,j,n) = half*(cfab(i-1,j,n) + cfab(i,j,n))
               enddo
            enddo
         enddo
      else
         do n = 1,nc
            do i = lo(1)+1, hi(1)-1
               do j = lo(2)+1, hi(2)
                  efab(i,j,n) = half*(cfab(i,j-1,n) + cfab(i,j,n))
               enddo
            enddo
         enddo
      end if
      end

      subroutine FORT_FLUX(lo, hi, fluxH, DIMS(fluxH),
     &                     fluxY, DIMS(fluxY), d, DIMS(d),
     &                     Ye, DIMS(Ye), Te, DIMS(Te), He, DIMS(He),
     &                     dir, add_enth_flux)
      implicit none
#include "cdwrk.H"
#include "conp.H"
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
      integer add_enth_flux

      integer i,j,lo1,lo2,hi1,hi2,K,L,n
      REAL_T Yt(maxspec), WW, CPMS(maxspec), X(maxspec)
      REAL_T rhoD, rhoDt(maxspec*maxspec), rhoTD(maxspec), PTC
      REAL_T sum_rhoYV, sum_rhoHYV, sum_rhoThetaD, P1atm, RU, RUC, sum

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

c            do K = 1,Nspec
c               Yt(K) = Ye(i,j,K)
c            enddo
            sum = 0.d0
            do n=1,Nspec
               Yt(n) = MAX(Ye(i,j,n),0.d0)
               sum = sum+Yt(n)
            end do
            if (iN2 .gt. 0) then
               Yt(iN2) = Yt(iN2)+one-sum
            endif

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

            do K = 1,Nspec

               sum_rhoYV = - Ye(i,j,K)*rhoTD(K)*d(i,j,Nspec+1)/Te(i,j)
               do L= 1,Nspec
                  sum_rhoYV = sum_rhoYV - rhoD(K,L)*d(i,j,L)
               enddo
               fluxY(i,j,K) = sum_rhoYV*1.d-1

               sum_rhoHYV = sum_rhoHYV + He(i,j,K)*fluxY(i,j,K)
               sum_rhoThetaD = sum_rhoThetaD + rhoTD(K)*d(i,j,K)

            enddo

            if (add_enth_flux.eq.1) then
               fluxH(i,j) = sum_rhoHYV
            else
               fluxH(i,j) = 0.d0
            endif

c           Now add all the other terms
            fluxH(i,j) = - (PTC*1.d-5)*d(i,j,Nspec+1)
     &           - (RU*Te(i,j)/WW)*(sum_rhoThetaD*1.d-1)

         enddo
      enddo
      
      end
      
      subroutine FORT_INCRDIV(lo, hi, div, DIMS(div), V, DIMS(V),
     &                        dir, nc)
      implicit none
#include "cdwrk.H"
      integer lo(SDIM), hi(SDIM), dir, nc
      integer DIMDEC(div)
      integer DIMDEC(V)
      REAL_T div(DIMV(div) ,1)
      REAL_T   V(DIMV(V)   ,1)

      integer i,j,n

c     Note: We're incrementing with div(dx.Area.Flux), so the dx's
c     cancel, and we assume V is extensive (ie. Flux.Area)
      if (dir .EQ. 0) then
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               do n= 1,nc
                  div(i,j,n) = div(i,j,n)
     &                 + (V(i+1,j,n) - V(i,j,n))
               enddo
            enddo               
         enddo
      else 
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               do n= 1,nc
                  div(i,j,n) = div(i,j,n)
     &                 + (V(i,j+1,n) - V(i,j,n))
               enddo
            enddo               
         enddo
      end if
      end


      subroutine FORT_CPSCALE(lo,hi,lambda,DIMS(lambda),
     &                        Y,DIMS(Y),T,DIMS(T))
      implicit none
#include "cdwrk.H"
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
            enddo

            CALL CKCPBS(T(i,j),Yt,IWRK(ckbi),RWRK(ckbr),CPMIX)
            lambda(i,j) = lambda(i,j)/(CPMIX*1.d-4)
            
         enddo
      enddo
      end


      subroutine FORT_THERM(lo, hi, lambda, DIMS(lambda),
     &                     Ye, DIMS(Ye), Te, DIMS(Te), dir)
      implicit none
#include "cdwrk.H"
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
               enddo
               
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

            enddo
         enddo

      else 

         hi1 = hi(1)
         hi2 = hi(2) + 1
         do j = lo2, hi2
            do i = lo1, hi1

               do K = 1,Nspec
                  Yt(K) = Ye(i,j,K)
               enddo

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

            enddo
         enddo
      end if
      
      end

      subroutine FORT_CRSNCCBND(lo, hi, fine, DIMS(fine),
     &                          crse, DIMS(crse), nc, face, ratio)
      integer lo(SDIM), hi(SDIM), nc
      integer DIMDEC(fine)
      integer DIMDEC(crse)
      REAL_T fine(DIMV(fine),nc)
      REAL_T crse(DIMV(crse),nc)
      integer face, ratio(SDIM)
      
      integer icrse,jcrse,n,ii,jj,ifine,jfine, clo1, chi1, clo2, chi2, offset, dir
      
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
      
c     NOTE: need to generalize to rz
      if (dir .eq. 0) then
         clo1 = (lo(1)+offset)/ratio(1)-offset
         clo2 = lo(2)/ratio(2)
         chi1 = (hi(1)+offset)/ratio(1)-offset
         chi2 = (hi(2)+1)/ratio(2) - 1
         do jcrse=clo2,chi2
            do icrse=clo1,chi1
               do n=1,nc
                  crse(icrse,jcrse,n) = 0.d0
                  ifine = ratio(1)*(icrse + offset) - offset
                  do jj = 1,ratio(2)
                     jfine = ratio(2)*jcrse + jj - 1
                     crse(icrse,jcrse,n) = crse(icrse,jcrse,n) + fine(ifine,jfine,n)
                  enddo
                  crse(icrse,jcrse,n) = crse(icrse,jcrse,n)/ratio(2)
               enddo
            enddo
         enddo
      else
         clo1 = lo(1)/ratio(1)
         clo2 = (lo(2)+offset)/ratio(2)-offset
         chi1 = (hi(1)+1)/ratio(1) - 1
         chi2 = (hi(2)+offset)/ratio(2)-offset
         do jcrse=clo2,chi2
            do icrse=clo1,chi1
               do n=1,nc
                  crse(icrse,jcrse,n) = 0.d0
                  jfine = ratio(2)*(jcrse + offset) - offset
                  do ii = 1,ratio(1)
                     ifine = ratio(1)*icrse + ii - 1
                     crse(icrse,jcrse,n) = crse(icrse,jcrse,n) + fine(ifine,jfine,n) 
                  enddo
                  crse(icrse,jcrse,n) = crse(icrse,jcrse,n)/ratio(1)
               enddo
            enddo
         enddo
      end if
      end

      subroutine FORT_DDCCAVG(c, DIMS(c), f, DIMS(f), lo, hi, nc, ratio)
      integer nc
      integer DIMDEC(f)
      integer DIMDEC(c)
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      REAL_T f(DIMV(f),nc)
      REAL_T c(DIMV(c),nc)
      integer ratio(SDIM)

      integer i,j,n,ii,jj
      REAL_T fac

      fac = 1.d0 / (ratio(1)*ratio(2))

c     NOTE: generalize to rz
      do n = 1, nc
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               c(i,j,n) = 0.d0
               do jj = 1,ratio(2)
                  do ii = 1,ratio(1)
                     c(i,j,n) = c(i,j,n) + f(2*i+ii-1,2*j+jj-1,n)
                  enddo
               enddo
               c(i,j,n) = c(i,j,n)*fac
            enddo
         enddo
      enddo
      end

      subroutine FORT_DDCCINT(f, DIMS(f), c, DIMS(c), lo, hi, nc, ratio)
      integer nc
      integer DIMDEC(f)
      integer DIMDEC(c)
      integer lo(BL_SPACEDIM)
      integer hi(BL_SPACEDIM)
      REAL_T f(DIMV(f),nc)
      REAL_T c(DIMV(c),nc)
      integer ratio(SDIM)
c
      integer i,j,n,ii,jj

c     NOTE: generalize to rz
      do n = 1, nc
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               do jj = 1,ratio(2)
                  do ii = 1,ratio(1)
                     f(2*i+ii-1,2*j+jj-1,n) = c(i,j,n) + f(2*i+ii-1,2*j+jj-1,n)
                  enddo
               enddo
            enddo
         enddo
      enddo
c
      end
      
