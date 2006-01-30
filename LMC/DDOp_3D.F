#include "DDOp_F.H"
#include "ChemDriver_F.H"
#include "ArrayLim.H"
#include "CONSTANTS.H"

#define SDIM 3

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
      integer i,j,k,n

      dxinv = one / dx(dir+1)

      if (dir .EQ. 0) then
         do k=lo(3),hi(3)
            do j=lo(2),hi(2)
               do i=lo(1),hi(1)+1
                  do n=1,Nspec
                     d(i,j,k,n) = (X(i,j,k,n) - X(i-1,j,k,n))*dxinv
                  enddo
                  d(i,j,k,Nspec+1) = (T(i,j,k) - T(i-1,j,k))*dxinv
               enddo
            enddo
         enddo
      else if (dir .EQ. 1) then
         do k=lo(3),hi(3)
            do j=lo(2),hi(2)+1
               do i=lo(1),hi(1)
                  do n=1,Nspec
                     d(i,j,k,n) = (X(i,j,k,n) - X(i,j-1,k,n))*dxinv
                  enddo
                  d(i,j,k,Nspec+1) = (T(i,j,k) - T(i,j-1,k))*dxinv
               enddo
            enddo
         enddo
      else
         do k=lo(3),hi(3)+1
            do j=lo(2),hi(2)
               do i=lo(1),hi(1)
                  do n=1,Nspec
                     d(i,j,k,n) = (X(i,j,k,n) - X(i,j,k-1,n))*dxinv
                  enddo
                  d(i,j,k,Nspec+1) = (T(i,j,k) - T(i,j,k-1))*dxinv
               enddo
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

      integer i,j,k,n

c     Here, the cc data is assumed to be computable over the box (lo,hi) 
c     surroundingNodes of grow((lo,hi),-1) in the dir direction.
      if (dir .EQ. 0) then
         do n = 1,nc
            do k = lo(3)+1, hi(3)-1
               do j = lo(2)+1, hi(2)-1
                  do i = lo(1)+1, hi(1)
                     efab(i,j,k,n) = half*(cfab(i-1,j,k,n) + cfab(i,j,k,n))
                  enddo
               enddo
            enddo
         enddo
      else if (dir .EQ. 1) then
         do n = 1,nc
            do k = lo(3)+1, hi(3)-1
               do j = lo(2)+1, hi(2)
                  do i = lo(1)+1, hi(1)-1
                     efab(i,j,k,n) = half*(cfab(i,j-1,k,n) + cfab(i,j,k,n))
                  enddo
               enddo
            enddo
         enddo
      else
         do n = 1,nc
            do k = lo(3)+1, hi(3)
               do j = lo(2)+1, hi(2)-1
                  do i = lo(1)+1, hi(1)-1
                     efab(i,j,k,n) = half*(cfab(i,j,k-1,n) + cfab(i,j,k,n))
                  enddo
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

      integer i,j,k,lo1,lo2,lo3,hi1,hi2,hi3,L,M
      REAL_T Yt(maxspec), WW, CPMS(maxspec), X(maxspec)
      REAL_T rhoD, rhoDt(maxspec*maxspec), rhoTD(maxspec), PTC
      REAL_T sum_rhoYV, sum_rhoHYV, sum_rhoThetaD, P1atm, RU, RUC

      rhoD(L,M) = rhoDt( (M-1)*Nspec + L )
      
      lo1 = lo(1)
      lo2 = lo(2)
      lo3 = lo(3)
      if (dir .EQ. 0) then
         hi1 = hi(1) + 1
         hi2 = hi(2)
         hi3 = hi(3)
      elseif (dir .EQ. 1) then
         hi1 = hi(1)
         hi2 = hi(2) + 1
         hi3 = hi(3)
      else 
         hi1 = hi(1)
         hi2 = hi(2)
         hi3 = hi(3) + 1
      end if
      
      CALL CKRP(IWRK(ckbi), RWRK(ckbr), RU, RUC, P1atm)
      
      do k = lo3, hi3
         do j = lo2, hi2
            do i = lo1, hi1

               do L = 1,Nspec
                  Yt(L) = Ye(i,j,k,L)
               enddo
               
               CALL CKCPMS(Te(i,j,k),IWRK(ckbi),RWRK(ckbr),CPMS)
               CALL CKYTX(Yt,IWRK(ckbi),RWRK(ckbr),X)
               CALL CKMMWY(Yt, IWRK(ckbr), RWRK(ckbr), WW)
               CALL EGSPAR(Te(i,j,k),X,Yt,CPMS,RWRK(egbr),IWRK(egbi))
               CALL EGSLTDR5(Te(i,j,k),Yt,WW,RWRK(egbr),IWRK(egbi),
     &                       PTC,rhoTD,rhoDt)
            
c     Note: scalings get cgs quantities from EGLib to MKS (assume He,d come in MKS)
               sum_rhoYV     = zero
               sum_rhoHYV    = zero
               sum_rhoThetaD = zero
               
               do L= 1,Nspec
                  
                  sum_rhoYV = - Ye(i,j,k,L)*rhoTD(L)*d(i,j,k,Nspec+1)/Te(i,j,k)
                  do M = 1,Nspec
                     sum_rhoYV = sum_rhoYV - rhoD(L,M)*d(i,j,k,M)
                  enddo
                  fluxY(i,j,k,L) = sum_rhoYV*1.d-1
                  
                  sum_rhoHYV = sum_rhoHYV + He(i,j,k,L)*fluxY(i,j,k,L)
                  sum_rhoThetaD = sum_rhoThetaD + rhoTD(K)*d(i,j,k,L)
                  
               enddo
               
               if (add_enth_flux.eq.1) then
                  fluxH(i,j,k) = sum_rhoHYV
               else
                  fluxH(i,j,k) = 0.d0
               endif
                  
c           Now add all the other terms
               fluxH(i,j,k) = - (PTC*1.d-5)*d(i,j,k,Nspec+1)
     &              - (RU*Te(i,j,k)/WW)*(sum_rhoThetaD*1.d-1)

            enddo
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

      integer i,j,k,n

c     Note: We're incrementing with div(dx.Area.Flux), so the dx's
c     cancel, and we assume V is extensive (ie. Flux.Area)
      if (dir .EQ. 0) then
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  do n= 1,nc
                     div(i,j,k,n) = div(i,j,k,n)
     &                    + (V(i+1,j,k,n) - V(i,j,k,n))
                  enddo
               enddo
            enddo               
         enddo
      elseif (dir .EQ. 1) then
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  do n= 1,nc
                     div(i,j,k,n) = div(i,j,k,n)
     &                    + (V(i,j+1,k,n) - V(i,j,k,n))
                  enddo
               enddo
            enddo               
         enddo
      else 
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  do n= 1,nc
                     div(i,j,k,n) = div(i,j,k,n)
     &                    + (V(i,j,k+1,n) - V(i,j,k,n))
                  enddo
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
      integer i,j,k,n
      
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               
               do n = 1,Nspec
                  Yt(n) = Y(i,j,k,n)
               enddo
               
               CALL CKCPBS(T(i,j,k),Yt,IWRK(ckbi),RWRK(ckbr),CPMIX)
               lambda(i,j,k) = lambda(i,j,k)/(CPMIX*1.d-4)
               
            enddo
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

      integer i,j,k,n,lo1,lo2,lo3,hi1,hi2,hi3
      REAL_T Yt(maxspec), WW, CPMS(maxspec), X(maxspec), TC

c     Note: scaling is to get D1 units of g/(s.cm^2) into MKS
      
      lo1 = lo(1)
      lo2 = lo(2)
      lo3 = lo(3)
      if (dir .EQ. 0) then
         hi1 = hi(1) + 1
         hi2 = hi(2)
         hi3 = hi(3)
         
         do k = lo3, hi3
            do j = lo2, hi2
               do i = lo1, hi1
               
                  do n = 1,Nspec
                     Yt(n) = Ye(i,j,k,n)
                  enddo
                  
                  CALL CKCPMS(Te(i,j,k),IWRK(ckbi),RWRK(ckbr),CPMS)
                  CALL CKYTX(Yt,IWRK(ckbi),RWRK(ckbr),X)
                  CALL EGSPAR(Te(i,j,k),X,Yt,CPMS,RWRK(egbr),IWRK(egbi))
                  CALL EGSL4(Te(i,j,k),Yt,RWRK(egbr),IWRK(egbi),TC)
                  
                  if(i-1.ge.lo(1))then
                     lambda(i-1,j,k) = lambda(i-1,j,k)+TC*1.d-5
                  endif
                  if(i.le.hi(1))then
                     lambda(i,j,k) = lambda(i,j,k)+TC*1.d-5
                  endif

               enddo
            enddo
         enddo

      elseif (dir .EQ. 1) then
         hi1 = hi(1)
         hi2 = hi(2) + 1
         hi3 = hi(3)
         
         do k = lo3, hi3
            do j = lo2, hi2
               do i = lo1, hi1
               
                  do n = 1,Nspec
                     Yt(n) = Ye(i,j,k,n)
                  enddo
                  
                  CALL CKCPMS(Te(i,j,k),IWRK(ckbi),RWRK(ckbr),CPMS)
                  CALL CKYTX(Yt,IWRK(ckbi),RWRK(ckbr),X)
                  CALL EGSPAR(Te(i,j,k),X,Yt,CPMS,RWRK(egbr),IWRK(egbi))
                  CALL EGSL4(Te(i,j,k),Yt,RWRK(egbr),IWRK(egbi),TC)
                  
                  if(j-1.ge.lo(2))then
                     lambda(i,j-1,k) = lambda(i,j-1,k)+TC*1.d-5
                  endif
                  if(j.le.hi(2))then
                     lambda(i,j,k) = lambda(i,j,k)+TC*1.d-5
                  endif

               enddo
            enddo
         enddo

      else 

         hi1 = hi(1)
         hi2 = hi(2)
         hi3 = hi(3) + 1

         do k = lo3, hi3
            do j = lo2, hi2
               do i = lo1, hi1

                  do n = 1,Nspec
                     Yt(n) = Ye(i,j,k,n)
                  enddo
                  
                  CALL CKCPMS(Te(i,j,k),IWRK(ckbi),RWRK(ckbr),CPMS)
                  CALL CKYTX(Yt,IWRK(ckbi),RWRK(ckbr),X)
                  CALL EGSPAR(Te(i,j,k),X,Yt,CPMS,RWRK(egbr),IWRK(egbi))
                  CALL EGSL4(Te(i,j,k),Yt,RWRK(egbr),IWRK(egbi),TC)
                  
                  if(k-1.ge.lo(3))then
                     lambda(i,j,k-1) = lambda(i,j,k-1)+TC*1.d-5
                  endif
                  if(k.le.hi(3))then
                     lambda(i,j,k) = lambda(i,j,k)+TC*1.d-5
                  endif
                  
               enddo
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
      
      integer i,j,k,n,ii,jj,kk, clo1, chi1, clo2, chi2, clo3, chi3, offset, dir
      
      if ((face.EQ.0).OR.(face.EQ.3)) then
         dir = 0
      elseif ((face.EQ.1).OR.(face.EQ.4)) then
         dir = 1
      else
         dir = 2
      end if

      if (face .LT. SDIM) then
         offset = 1
      else
         offset = 0
      end if
      
      if (dir .eq. 0) then
         clo1 = (lo(1)-offset)/2
         clo2 = lo(2)/2
         clo3 = lo(3)/2
         chi1 = (hi(1)-offset)/2
         chi2 = (hi(2)+1)/2 - 1
         chi3 = (hi(3)+1)/2 - 1
         do k=clo3,chi3
            do j=clo2,chi2
               do i=clo1,chi1
                  do n=1,nc
                     crse(i,j,k,n) = 0.d0
                     do jj=1,ratio(2)
                        do kk=1,ratio(3)
                           crse(i,j,k,n) = crse(i,j,k,n) + fine(2*i+offset,2*j+jj-1,2*k+kk-1,n)
                        enddo
                     enddo
                     crse(i,j,k,n) = crse(i,j,k,n) / (ratio(2)*ratio(3))
                  enddo
               enddo
            enddo
         enddo
      elseif (dir .eq. 1) then
         clo1 = lo(1)/2
         clo2 = (lo(2)-offset)/2
         clo3 = lo(3)/2
         chi1 = (hi(1)+1)/2 - 1
         chi2 = (hi(2)-offset)/2
         chi3 = (hi(3)+1)/2 - 1
         do k=clo3,chi3
            do j=clo2,chi2
               do i=clo1,chi1
                  do n=1,nc
                     crse(i,j,k,n) = 0.d0
                     do ii=1,ratio(1)
                        do kk=1,ratio(3)
                           crse(i,j,k,n) = crse(i,j,k,n) + fine(2*i+ii-1,2*j+offset,2*k+kk-1,n)
                        enddo
                     enddo
                     crse(i,j,k,n) = crse(i,j,k,n) / (ratio(1)*ratio(3))
                  enddo
               enddo
            enddo
         enddo
      else
         clo1 = lo(1)/2
         clo2 = lo(2)/2
         clo3 = (lo(3)-offset)/2
         chi1 = (hi(1)+1)/2 - 1
         chi2 = (hi(2)+1)/2 - 1
         chi3 = (hi(3)-offset)/2
         do k=clo3,chi3
            do j=clo2,chi2
               do i=clo1,chi1
                  do n=1,nc
                     crse(i,j,k,n) = 0.d0
                     do ii=1,ratio(1)
                        do jj=1,ratio(2)
                           crse(i,j,k,n) = crse(i,j,k,n) + fine(2*i+ii-1,2*j+jj-1,2*k+offset,n)
                        enddo
                     enddo
                     crse(i,j,k,n) = crse(i,j,k,n) / (ratio(1)*ratio(2))
                  enddo
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
c
      integer i,j,k,n,ii,jj,kk
      REAL_T fac

      fac = 1.d0/(ratio(1)*ratio(2)*ratio(3))
c     NOTE: generalize to rz
      do n = 1, nc
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  c(i,j,k,n) = 0.d0
                  do kk = 1,ratio(3)
                     do jj = 1,ratio(2)
                         do ii = 1,ratio(1)
                            c(i,j,k,n) = c(i,j,k,n) + f(2*i+ii-1,2*j+jj-1,2*k+kk-1,n)
                         enddo
                     enddo
                  enddo
                  c(i,j,k,n) = c(i,j,k,n)*fac
               enddo
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
      integer i,j,k,n,ii,jj,kk

c     NOTE: generalize to rz
      do n = 1, nc
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  do kk = 1,ratio(3)
                     do jj = 1,ratio(2)
                        do ii = 1,ratio(1)
                           f(2*i+ii-1,2*j+jj-1,2*k+kk-1,n) = c(i,j,k,n) + f(2*i+ii-1,2*j+jj-1,2*k+kk-1,n)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      end
      
