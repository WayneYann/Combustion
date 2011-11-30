module chemsolv_module

   implicit none

   contains

      subroutine consteFY(N, TIME, Z, ZP, RPAR, IPAR)

      use  conp_module
      use cdwrk_module

      implicit none

      double precision :: TIME, Z(NEQ), ZP(NEQ), RPAR(*)
      integer          :: N, IPAR(*)
      
      double precision :: RHOcgs, CPB, SUM, H, WDOT, WT, EINT, TEMP
      integer          :: K, Niter

      double precision :: Conc(maxspec),Enthalpy(maxspec),wdots(maxspec)

      EINT = RPAR(NP)
      Niter = T_from_eY(Z(1),Z(2),EINT)

      if (Niter.lt.0) &
         call bl_error('consteFY: T_from_eY failed!!!!')

      CALL CKCPBS(Z(1),Z(2),IPAR(ckbi),RPAR(ckbr),CPB)

      RHOcgs = RWRK(NRHO)*1.d-3
      CALL CKYTCR(RHOcgs, Z(1), Z(2), IPAR(ckbi), RPAR(ckbr), Conc)
      TEMP = MAX(1.d0,Z(1))

      do K=1,Nspec
         Conc(K) = max(Conc(K),0d0)
      end do

      CALL CKWC(TEMP, Conc, IPAR(ckbi), RPAR(ckbr), wdots)
      CALL CKHMS(Z(1), IPAR(ckbi), RPAR(ckbr), Enthalpy)
      SUM = 0.d0

      do K = 1, Nspec
         H    = Enthalpy(K)
         WDOT = wdots(K)
         WT   = RPAR(NWT   + K - 1)
         ZP(K+1) = WDOT * WT / RHOcgs
         SUM = SUM + H * WDOT * WT
      end do

      ZP(1) = -SUM / (RHOcgs*CPB)

      end subroutine consteFY

!     
!     *********************************************************************************************     
!     
      integer function T_from_eY(T,Y,ein)

      use cdwrk_module

      implicit none

      double precision T,Y(1),e,ein
      double precision TMIN,TMAX,errMAX
      integer NiterMAX,Niter,n,NiterDAMP
!      parameter (TMIN=10, TMAX=8000, errMAX=1.e-8, NiterMAX=20)
      parameter (TMIN=300, TMAX=6500, errMAX=1.e-8, NiterMAX=20)
!      parameter (NiterDAMP = NiterMAX/2)
      parameter (NiterDAMP = NiterMAX)
      double precision  T0,h,cp,cv,de,temp,RoverWbar,Wbar,RU,RUC,P1ATM
      double precision res(0:NiterMAX-1),dT, etarg
      logical out_of_bounds, converged, soln_bad, stalled
      double precision e300,cv300,e6500,cv6500
      integer ihitlo,ihithi

      out_of_bounds(temp) = (temp.lt.TMIN-1.d0) .or. (temp.gt.TMAX)

      if ((T.GE.TMIN).and.(T.LE.TMAX)) then
         T0 = T
      else
         T0 = 0.5d0*(TMIN+TMAX)
         T = T0
      end if
      Niter = 0
      de = 0.d0
      soln_bad = .FALSE.
      etarg = ein * 1.d4
      ihitlo = 0
      ihithi = 0

      CALL CKUBMS(T,Y,IWRK(ckbi),RWRK(ckbr),e)
!     CALL CKHBMS(T,Y,IWRK(ckbi),RWRK(ckbr),h)
!     CALL CKRP(IWRK(ckbi),RWRK(ckbr),RU,RUC,P1ATM)
!     call CKMMWY(Y,IWRK(ckbi),RWRK(ckbr),Wbar)
!     RoverWbar = RU / Wbar
!     e  = h - T*RoverWbar
      de = 2.d0*ABS(e - etarg)/(1.d0 + ABS(e) + ABS(etarg))
      res(Niter) = de
      converged = de.le.errMAX

      do while ((.not.converged) .and. (.not.soln_bad))
!        CALL CKCPBS(T,Y,IWRK(ckbi),RWRK(ckbr),cp)
         CALL CKCVBS(T,Y,IWRK(ckbi),RWRK(ckbr),cv)
         dT = (etarg - e)/cv
         if ((Niter.le.NiterDAMP).and.(T+dT.ge.TMAX)) then
!           T = 0.5d0*(T + TMAX)
            T = TMAX
            ihithi = 1
         else if ((Niter.le.NiterDAMP).and.(T+dT.le.TMIN)) then
!           T = 0.5d0*(T + TMIN)
            T = TMIN
            ihitlo = 1
         else
            T = T + dT
         end if
         soln_bad = out_of_bounds(T)
         if (soln_bad) then
            T_from_eY = -1
            goto 100
         else
!           CALL CKHBMS(T,Y,IWRK(ckbi),RWRK(ckbr),h)
!           e  = h - T*RoverWbar
            CALL CKUBMS(T,Y,IWRK(ckbi),RWRK(ckbr),e)
            de = 2.d0*ABS(e - etarg)/(1.d0 + ABS(e) + ABS(etarg))
            res(Niter) = de
            Niter = Niter + 1
         end if
         if (Niter .ge. NiterMAX) then
            T_from_eY = -2
            goto 100
         endif
         converged = (de.le.errMAX) .or. (ABS(dT).le.errMAX)

         if((ihitlo.eq.1).and.(e.gt.etarg))then
            T = 300.d0
            CALL CKUBMS(T,Y,IWRK(ckbi),RWRK(ckbr),e300)
            CALL CKCVBS(T,Y,IWRK(ckbi),RWRK(ckbr),cv300)
            T=300.d0+(etarg-e300)/cv300
            converged = .true.
         endif
         if((ihithi.eq.1).and.(e.lt.etarg))then
            T = 6500.d0
            CALL CKUBMS(T,Y,IWRK(ckbi),RWRK(ckbr),e6500)
            CALL CKCVBS(T,Y,IWRK(ckbi),RWRK(ckbr),cv6500)
            T=6500.d0+(etarg-e6500)/cv6500
            converged = .true.
         endif

      end do

!     Set max iters taken during this solve, and exit
      T_from_eY = Niter
      return

!     Error condition....dump state and bail out
 100  continue

! #ifdef VERBOSE_FAILURE
!      write(6,997) 'T from (e,Y): failed'
!      write(6,997) 'iterations tried = ',Niter
!      write(6,998) 'initial T = ',T0
!      write(6,998) 'current T = ',T
!      write(6,998) 'species mass fracs:'
!      do n = 1,Nspec
!         write(6,998) '  ',Y(n)
!      end do
!      write(6,998)
!      write(6,998) 'residual = e - h + RT/Wbar [cgs]'
!      do n = 0,Niter-1
!         write(6,998) '  ',res(n)
!      end do
!
! 997  format(a,3(i4,a))
! 998  format(a,d21.12)
!#endif
      
      end function T_from_eY

end module chemsolv_module
