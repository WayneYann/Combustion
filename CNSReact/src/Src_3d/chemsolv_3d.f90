      subroutine chemsolv(lo, hi, u, u_l1, u_l2, u_l3, u_h1, u_h2, u_h3, dt, idbg)

      use meth_params_module, only : NTHERM, NVAR, nadv 

      use chemsolv_module
!     use    cdwrk_module
!     use     conp_module

      implicit none

      include "cdwrk.H"
      include "conp.H"

      integer          :: lo(:), hi(:)
      integer          :: u_l1,u_h1,u_l2,u_h2,u_l3,u_h3
      integer          :: idbg
      double precision :: u(u_l1:u_h1, u_l2:u_h2, u_l3:u_h3,NVAR)
      double precision :: dt

      ! Local variables

      external conpJY

      double precision   :: TT1, TT2
      integer            :: i, j, k, n, MF, ISTATE,  Niter
      integer            :: hardestCase, mostfs, lout, open_vode_failure_file
      character*(maxspnml) name

      integer, parameter :: ITOL  = 1
      integer, parameter :: IOPT  = 1
      integer, parameter :: ITASK = 1

      double precision, parameter :: RTOL    = 1.d-8
      double precision, parameter :: ATOLEPS = 1.d-8
      double precision, parameter :: spec_scalT = 2000.d0
      double precision, parameter :: NEWJ_TOL = 1.d100

      double precision :: deriv(maxspec+1)
      double precision :: ATOL(maxspec+1)

      double precision :: Ytemp(maxspec),Yres(maxspec),sum

      logical          :: newJ_triggered

      double precision YJ_SAVE(80)
      LOGICAL FIRST
      COMMON /VHACK/ YJ_SAVE, FIRST
      SAVE   /VHACK/
      double precision scale, ekin, u1, u2, u3
!     integer is_bad
      
!     ! NRHO is declared in conp_module, but set here
!     NRHO = 1
      
!     Set IOPT=1 parameter settings for VODE
      RWRK(dvbr+4) = 0
      RWRK(dvbr+5) = 0
      RWRK(dvbr+6) = 1.d-19
      IWRK(dvbi+4) = 0
      IWRK(dvbi+5) = max_vode_subcycles
      IWRK(dvbi+6) = 0
      
      
!     is_bad = 0
      TT2 = dt

      if (nstiff.eq.1) then
!   chemkin provided finite difference jacobian
         MF = 22
      else
         MF = 10
      endif

      ATOL(1) = ATOLEPS

!     find C2H2 and O2 in the list
!     iC2H2 = -1
!     iO2   = -1
!     do n = 1,Nspec
!        call get_spec_name(name,n)
!        if (name .eq. 'C2H2' ) iC2H2 = n
!        if (name .eq. 'O2'   ) iO2 = n
!     end do

      hardestCase = 0
      mostfs = 0
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)

!               if ((iC2H2.gt.0).and.(iO2.gt.0).and.
!    &              (u(i,j,k,NTHERM+nadv+iC2H2).gt. 0.d0) .and.
!    &              (u(i,j,k,NTHERM+nadv+iO2  ).gt. 0.d0) ) then

                  TT1    = 0.d0
                  ISTATE = 1

!     Note: pass eint in NP slot, and density in rho slot
                  u1 = u(i,j,k,2)/u(i,j,k,1)
                  u2 = u(i,j,k,3)/u(i,j,k,1)
                  u3 = u(i,j,k,4)/u(i,j,k,1)
                  ekin = 0.5d0*(u1**2 + u2**2 + u3**2)
                  RWRK(NP) = u(i,j,k,5)/u(i,j,k,1) - ekin
                  RWRK(NRHO) = u(i,j,k,1)
                     
                  sum = 0.d0
                  do n=1,Nspec
                     Ytemp(n) = MAX(u(i,j,k,NTHERM+NADV+n)/u(i,j,k,1),0.d0)
!                    Ytemp(n) = u(i,j,k,NTHERM+NADV+n)/u(i,j,k,1)
                     sum = sum+Ytemp(n)
                  end do
                  if (iN2 .gt. 0) then
                     Ytemp(iN2) = Ytemp(iN2) + 1.d0 - sum
                  endif
                
!     Note: ensure that T initialized to a good value, load Y,T into IC
                  do n = 1,Nspec
                     RWRK(NZ+n) = Ytemp(n)
                  end do

                  RWRK(NZ) = u(i,j,k,NTHERM+1)
!wqz                  Niter = T_from_eY(RWRK(NZ),RWRK(NZ+1),RWRK(NP))
                  if (Niter.lt.0) then
                     write(6,*) 'FORT_CHEMSOLV: T_from_eY failed!!!!'
                     write(6,*) u(i,j,k,NTHERM+1),RWRK(NP)
                     do n=1,Nspec
                     write(*,*)Ytemp(n)
                     enddo
                     stop
                  end if
                  
                  newJ_triggered = .FALSE.
                  sum = 0.d0
                  do n=1,NEQ
                     scale = spec_scalT
                     if (n.ne.1) scale = spec_scalY(n-1)
                     sum = sum + ABS(RWRK(NZ+n-1)-YJ_SAVE(n))/scale
                  end do
                  if (sum .gt. NEWJ_TOL) then
                     FIRST = .TRUE.
                     newJ_triggered = .TRUE.
                  end if

                  CALL DVODE(consteFY, NEQ, RWRK(NZ), TT1, TT2, ITOL, RTOL, ATOL, &
                             ITASK, ISTATE, IOPT, RWRK(dvbr), dvr, IWRK(dvbi), &
                             dvi, conpJY, MF, RWRK, IWRK)

                  hardestCase = MAX(hardestCase,IWRK(dvbi+10))
                  mostfs = MAX(mostfs,IWRK(dvbi+11))

!wqz                  Niter = T_from_eY(RWRK(NZ),RWRK(NZ+1),RWRK(NP))
                  if (Niter.lt.0) then
                     write(6,*) 'FORT_CHEMSOLV: T_from_eY failed after DVODE call!!!!'
                     stop
                  end if

                  if (verbose_vode .eq. 1) then
                     write(6,*) '......dvode done:'
                     write(6,*) ' last successful step size = ',RWRK(dvbr+10)
                     write(6,*) '          next step to try = ',RWRK(dvbr+11)
                     write(6,*) '   integrated time reached = ',RWRK(dvbr+12)
                     write(6,*) '      number of time steps = ',IWRK(dvbi+10)
                     write(6,*) '              number of fs = ',IWRK(dvbi+11)
                     write(6,*) '              number of Js = ',IWRK(dvbi+12)
                     write(6,*) '    method order last used = ',IWRK(dvbi+13)
                     write(6,*) '   method order to be used = ',IWRK(dvbi+14)
                     write(6,*) '            number of LUDs = ',IWRK(dvbi+18)
                     write(6,*) ' number of Newton iterations ',IWRK(dvbi+19)
                     write(6,*) ' number of Newton failures = ',IWRK(dvbi+20)
                     if (ISTATE.eq.-4 .or. ISTATE.eq.-5) then
                        call get_spec_name(name,IWRK(dvbi+15))
                        write(6,*) '   spec with largest error = ',name
                     end if
                  end if
                  
                  if (ISTATE .LE. -1) then
                     
!wqz                     call consteFY(Nspec, TT1, RWRK(NZ), deriv, RWRK, IWRK)
                     lout = open_vode_failure_file()
                     
                     write(lout,*)
                     write(lout,995) 'VODE Failed at (i,j,k) = (' , i, &
                          ',' ,j, ',' ,k, '),   Return code = ',ISTATE
                     write(lout,996) 'time(T2,Tl,dt)  ',dt, TT1, dt-TT1
                     write(lout,995) 'State ID, old, last, dY/dt, dY/dt*(dt)'
                     write(lout,996) 'T ', &             
                          u(i,j,k,NTHERM+1),u(i,j,k,NTHERM+1),RWRK(NZ),deriv(1)
                     do n=1,Nspec
                        call get_spec_name(name,n)
                        write(lout,996) name,u(i,j,k,NTHERM+NADV+n)/u(i,j,k,1), &
                             Ytemp(n), &
                             RWRK(NZ+n),deriv(n+1)
                     end do
                     write(6,*) 'VODE failed...see drop file...exiting...'
 995                 format(a,4(i4,a))
 996                 format(a16,1x,5e30.22)
                     
                     if (ISTATE .LE. -4) then
                        do n=1,Nspec
                           RWRK(NZ+n) = Ytemp(n)
                        end do
                     end if
                  end if
                  
                  ! Set result T and Y back into state
                  u(i,j,k,NTHERM+1) = RWRK(NZ)               
                  do n= 1,Nspec
                     Yres(n) = RWRK(NZ+n)
                  end do

                  do n=1,Nspec
                     u(i,j,k,NTHERM+NADV+n) = &
                          u(i,j,k,NTHERM+NADV+n)+(Yres(n)-Ytemp(n))*u(i,j,k,1)
                  end do
!               end if
            end do
         end do
      end do

!     write(6,*) 'FORT_CHEMSOLVE: hardest case took ',hardestCase,' time steps'
!     write(6,*) 'FORT_CHEMSOLVE: biggest number of f evals was ',mostfs

      end subroutine chemsolv

