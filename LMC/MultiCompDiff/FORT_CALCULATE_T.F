      integer function FORT_CALCULATE_T (T,Hin,Y,errMax,NiterMAX,res)
c     JFG: FORT_CALCULATE_T is based on FORT_TfromHYpt.
c     if failure, then aborts rather than bails.
      implicit none
#include "cdwrk.H"
      REAL_T T,Y(*),H,Hin
      REAL_T TMIN,TMAX,errMAX
      integer NiterMAX,Niter,n,NiterDAMP
      parameter (TMIN=250, TMAX=5000)
      REAL_T  T0,cp,cv,dH,temp,RoverWbar,Wbar,RU,RUC,P1ATM
      REAL_T res(0:NiterMAX-1),dT, Htarg
      logical out_of_bounds, converged, soln_bad, stalled
      REAL_T h300,cp300,h6500,cp6500
      integer ihitlo,ihithi,j

      out_of_bounds(temp) = (temp.lt.TMIN-one) .or. (temp.gt.TMAX)

      NiterDAMP = NiterMAX
      if ((T.GE.TMIN).and.(T.LE.TMAX)) then
         T0 = T
      else
         T0 = half*(TMIN+TMAX)
         T = T0
      end if
      Niter = 0
      dH = zero
      soln_bad = .FALSE.
      Htarg = Hin * 1.d4
      ihitlo = 0
      ihithi = 0

      CALL CKHBMS(T,Y,IWRK(ckbi),RWRK(ckbr),H)
      dH = two*ABS(H - Htarg)/(one + ABS(H) + ABS(Htarg))
      res(Niter) = dH
      converged = dH.le.errMAX

      do while ((.not.converged) .and. (.not.soln_bad))

         CALL CKCPBS(T,Y,IWRK(ckbi),RWRK(ckbr),cp)
         dT = (Htarg - H)/cp
         if ((Niter.le.NiterDAMP).and.(T+dT.ge.TMAX)) then
            T = TMAX
            ihithi = 1
         else if ((Niter.le.NiterDAMP).and.(T+dT.le.TMIN)) then
            T = TMIN
            ihitlo = 1
         else
            T = T + dT
         end if
         soln_bad = out_of_bounds(T)
         if (soln_bad) then
            FORT_CALCULATE_T = -1
            goto 100
         else
            CALL CKHBMS(T,Y,IWRK(ckbi),RWRK(ckbr),H)
            dH = two*ABS(H - Htarg)/(one + ABS(H) + ABS(Htarg))
            res(Niter) = dH
            Niter = Niter + 1
         end if
         if (Niter .ge. NiterMAX) then
            FORT_CALCULATE_T = -2
            goto 100
         endif
         converged = (dH.le.errMAX) .or. (ABS(dT).le.errMAX)

         if ((ihitlo.eq.1).and.(H.gt.Htarg)) then
            T = TMIN
            CALL CKHBMS(T,Y,IWRK(ckbi),RWRK(ckbr),h300)
            CALL CKCPBS(T,Y,IWRK(ckbi),RWRK(ckbr),cp300)
            T=TMIN+(Htarg-h300)/cp300
            converged = .true.
         endif
         if ((ihithi.eq.1).and.(H.lt.Htarg)) then
            T = TMAX
            CALL CKHBMS(T,Y,IWRK(ckbi),RWRK(ckbr),h6500)
            CALL CKCPBS(T,Y,IWRK(ckbi),RWRK(ckbr),cp6500)
            T=TMAX+(Htarg-h6500)/cp6500
            converged = .true.
         endif
      end do
c
c     Set max iters taken during this solve and exit
c
      FORT_CALCULATE_T = Niter
      return
c
c     Error condition....dump state and abort
c
 100  continue

      write(6,997) 'FORT_CALCULATE_T: failed'
      write(6,997) 'iterations tried = ',Niter
      write(6,998) '  input h = ',Hin
      write(6,998) 'initial T = ',T0
      write(6,998) 'current T = ',T
      write(6,998) 'species mass fracs:'
      do n = 1,Nspec
         write(6,998) '  ',Y(n)
      end do
      write(6,998)
      write(6,998) 'residual:'
      do n = 0,Niter-1
         write(6,998) '  ',res(n)
      end do
      call bl_abort('FORT_CALCULATE_T: failed')

 997  format(a,3(i4,a))
 998  format(a,d21.12)
      end
