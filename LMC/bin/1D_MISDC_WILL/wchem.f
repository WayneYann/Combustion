      !     Will's Backward Euler chemistry solve routine
      module wchem

        implicit none

        private
        
        include 'spec.h'

        !     Jacobian matrix and factorization
        double precision, allocatable, save :: Jac(:,:), A(:,:)
        !     Pivot vector returned from LINPACK
        integer, allocatable, save :: ipvt(:)
        
        integer :: num_rhs
        integer :: num_jac
        
        double precision :: prev_time
        
        integer, parameter :: j_skip = 1
        
        private :: wvode_j, wvode_rhs
        
        public :: bechem, wvode_chem, rk_chem

      contains

        
        subroutine rk_chem(rYh0, rYh, dt)
          double precision, intent(in ) :: rYh0(Nspec+1), dt
          double precision, intent(out) :: rYh(Nspec+1)
          
          integer, parameter :: M = 30
          integer            :: n, NEQ, ipar(1)
          double precision   :: dt_rk, rpar(1), t
          double precision, dimension(Nspec+1) :: k1,k2,k3,k4,z
          
          if (.not. allocated(A)) then
            allocate(Jac(Nspec+1,Nspec+1))
            allocate(A(Nspec+1,Nspec+1))
            allocate(ipvt(Nspec+1))
          end if
          
          NEQ = Nspec+1
          ipar(1) = 0
          rpar(1) = dt
          dt_rk = dt/M
          
          rYh = rYh0
          
          num_jac = 0
          
          do n=1,M
            t  = n*dt_rk
            z  = rYh
            call wvode_rhs(NEQ, t, z, k1, rpar, ipar)
            k1 = k1*dt_rk
            t  = t + 0.5*dt_rk
            z  = rYh + 0.5*k1
            call wvode_rhs(NEQ, t, z, k2, rpar, ipar)
            k2 = dt_rk*k2
            z  = rYh + 0.5*k2
            call wvode_rhs(NEQ, t, z, k3, rpar, ipar)
            k3 = dt_rk*k3
            t  = t+0.5*dt_rk
            z  = rYh + k3
            call wvode_rhs(NEQ, t, z, k4, rpar, ipar)
            k4 = dt_rk*k4
            
            rYh = rYh + (k1 + 2*k2 + 2*k3 + k4)/6.0
          end do
        end subroutine rk_chem
        
        subroutine WVODE_J(NEQ, TIME, Y, ML, MU, PD, NRPD, RPAR, IPAR)
          integer NEQ, NRPD, ML, MU, IPAR(*)
          double precision TIME, Y(NEQ), PD(NRPD,NEQ), RPAR(*)
          
          double precision f1(NEQ)
          double precision f0(NEQ)
          
          double precision perturb, Yj
          integer j
          
          IPAR(1) = 1
          
          CALL WVODE_RHS(NEQ, TIME, Y, f0, RPAR, IPAR)
          
          perturb = 1.0d-13
          
          do j=1,NEQ
            Yj = Y(j)
            Y(j) = Y(j)+perturb
            
            CALL WVODE_RHS(NEQ, TIME, Y, f1, RPAR, IPAR)
            
            PD(:,j) = (f1-f0)/perturb
            
            Y(j) = Yj
          end do
          
          IPAR(1) = 0
          
          !stop
        end subroutine WVODE_J

        ! WVODE_RHS supplies the right-hand side F for the ODE y' = F
        ! in our case, the ODE in question is
        !            y' = (I - t J)^-1 (Q + wdot)
        ! where J is the Jacobian supplied by DWDOT
        ! and wdot is the production rate
        ! the values of y are input in the argument Z
        ! and the right-hand side is output in the argument ZDOT 
        subroutine WVODE_RHS(NEQ, TIME, Z, ZDOT, RPAR, IPAR)
          double precision :: TIME, Z(Nspec+1), ZDOT(Nspec+1), RPAR(*)
          integer          :: NEQ, IPAR(*)
          
          double precision :: rho, rhoinv, T, cp, cpinv, cons_p
          double precision :: Y(Nspec+1), wdot(Nspec), C(Nspec)
          
          double precision :: rcond, rwrk(Nspec+1)
          
          integer n, i, j, IWRK
          
          integer NiterMAX, Niter
          parameter (NiterMAX = 40)
          double precision res(NiterMAX), errMAX
          
          errMax = hmix_TYP*1.e-21
          
          !if (IPAR(2) .eq. 135 .and. TIME .ne. prev_time) then
          !  print *,TIME,Z
          !end if
          
          rho = 0.d0
          do n = 1,Nspec
            rho = rho + Z(n)
          end do
          rhoinv = 1.d0/rho
          
          Y = Z*rhoinv
          T = T_INIT
          
          call FORT_TfromHYpt(T,Y(Nspec+1),Y,Nspec,errMax,NiterMAX,res,Niter)
          T_INIT = T
          
          if (Niter.lt.0) then
            print *,'WVODE_RHS: H to T solve failed'
            print *,'Niter =', Niter
            stop
          end if
          
          !     compute the production rates
          call CKWYR(rho, T, Y, iwrk, rwrk, wdot)
          
!          if ((prev_time .ne. TIME .or. num_rhs .eq. 0)
!     $         .and. MOD(num_rhs,50) .eq. 0) then
          if (IPAR(1) .eq. 1 .or. mod(num_jac, j_skip) .eq. 0) then
             num_jac = num_jac + 1
             !     working under constant pressure
             cons_p = 1
             !     compute the molar concentrations
             call CKYTCR(rho, T, Y, iwrk, rwrk, C)
             !     use the concentrarions to compute the reaction Jacobian
             call DWDOT(Jac, C, T, cons_p)
             !     compute C_p
             call CKCPBS(T, Y, iwrk, rwrk, cp)
             cpinv = 1.d0/cp

             !     convert to the appropriate units
             do j=1,Nspec
                do i=1,Nspec
                   Jac(i,j) = Jac(i,j) * mwt(i) * invmwt(j)
                end do
                !     last row is derivative of enthalpy
                i=Nspec+1
                Jac(i, j) = 0.d0
                !Jac(i,j) = Jac(i,j) * invmwt(j) * cp * rho
             end do
             Jac(Nspec+1,Nspec+1) = 0.d0
             
             !     last column is derivative wrt enhtalpy
             j = Nspec+1
             do i=1,Nspec
                Jac(i,j) = Jac(i,j) * mwt(i) * cpinv * rhoinv
             end do
             
             !     we need to compute the matrix given by (I - t J)
             A = -TIME*Jac
             do i=1,Nspec+1
                A(i,i) = 1.d0 + A(i,i)
             end do
            
             !     call LINPACK to get the LU factorization of the matrix A
             !     call dgefa(A, Nspec+1, Nspec+1, ipvt, info)
             call DGECO(A, Nspec+1, Nspec+1, ipvt, rcond, rwrk)
             
             prev_time = TIME
          end if
          
          !     form the right-hand side for the linear solve
          !     (I - tJ)x = Q + wdot
          do n=1,Nspec
            ZDOT(n) = c_0(n) + c_1(n)*TIME + wdot(n)*mwt(n)
          end do
          ZDOT(Nspec+1) = c_0(0) + c_1(0)*TIME
          
          !     call LINPACK to solve the linear system Ax = r
          !     using the output from the factorization given by dgefa
          call dgesl(A, Nspec+1, Nspec+1, ipvt, ZDOT, 0)
          
          num_rhs = num_rhs + 1
        end subroutine WVODE_RHS


        !     wvode_chem: use VODE to compute the solution 
        !                 to the Backward Euler step
        !     rYh0 : (in)  initial condition for the ODE
        !     rYh  : (out) solution to the ODE
        !     dt   : (in)  time step
        subroutine wvode_chem(rYh0, rYh, dt, i)
          double precision, intent(in ) :: rYh0(Nspec + 1)
          double precision, intent(out) ::  rYh(Nspec + 1)
          double precision, intent(in ) ::   dt
          integer,          intent(in ) ::    i
          
          ! VODE stuff
          integer NEQ, ITOL, IOPT, ITASK
          parameter (ITOL=1, IOPT=1, ITASK=1)
          double precision RTOL, ATOLEPS, TT1, TT2
          parameter (RTOL=1.0D-14, ATOLEPS=1.0D-14)
          integer MF, ISTATE
          
          ! need this
          double precision YJ_SAVE(80)
          logical FIRST
          common /VHACK/ YJ_SAVE, FIRST
          save   /VHACK/
          
          ! DVODE workspace requirements      
          integer dvr, dvi
          parameter (dvr = 22 + 9*(maxspec+1) + 2*(maxspec+1)**2)
          parameter (dvi = 30 + maxspec + 1)
          double precision DVRWRK(dvr)
          integer DVIWRK(dvi)

          double precision RPAR(1)
          integer IPAR(2)

          FIRST = .true.

          if (.not. allocated(A)) then
            allocate(Jac(Nspec+1,Nspec+1))
            allocate(A(Nspec+1,Nspec+1))
            allocate(ipvt(Nspec+1))
          end if
          
          DVRWRK(5)  = 0.d0
          DVRWRK(6)  = 0.d0
          DVRWRK(7)  = min_vode_timestep
          DVRWRK(8)  = 0.d0
          DVRWRK(9)  = 0.d0
          DVRWRK(10) = 0.d0
          DVIWRK(5)  = 2
          DVIWRK(6)  = max_vode_subcycles
          !DVIWRK(6)  = 5000
          !DVRWRK(7)  = dt/DVIWRK(6)
          DVIWRK(7)  = 0
          DVIWRK(8)  = 0
          DVIWRK(9)  = 0
          DVIWRK(10) = 0
          
          MF = 22
          ! user supplied Jacobian
          MF = 22
          
          ISTATE = 1
          NEQ = Nspec + 1
          
          rYh = rYh0
          
          TT1 = 0.d0
          TT2 = TT1 + dt
          
          num_rhs = 0
          prev_time = 0.d0
          num_jac = 0
          
          ! we use IPAR=0 to indicate a "normal" call to the right hand side
          ! when the rhs is called from the Jacobian compute routine, IPAR 
          ! is set to 1, which signals a full recompute
          IPAR(1) = 0
          IPAR(2) = i
          ! RPAR is set to dt
          RPAR = dt
          
          do while (TT1 .lt. dt)
             if (TT1 .ne. 0) then
               print *, 'trying again...'
               print *, 'i = ', i
             end if
             CALL DVODE(WVODE_RHS, NEQ, rYh, TT1, TT2, ITOL, RTOL, ATOLEPS,
     &                  ITASK, ISTATE, IOPT, DVRWRK, dvr, DVIWRK,
     &                  dvi, WVODE_J, MF, RPAR, IPAR)
             ISTATE = 1
          end do
          
          if (i .eq. 135) then
            print *,'calls to rhs = ', num_rhs
          end if
          
!          if (num_rhs .gt. 1000) then
!             print *,'calls to rhs: ', num_rhs, '  jacobian recomputes: ', num_jac
!          end if
        end subroutine wvode_chem

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         !     do a Backward Euler solve for the chemistry using Newton's method
         !     Y0 : initial guess
         !     rho: density
         !     YT : input: right-hand side
         !     output: solution
         !     dt
         !     ierr : (optional) error
         subroutine bechem(Y0, rho, YT, dt, ierr)
           double precision, intent(in   ) :: dt
           double precision, intent(in   ) :: rho
           double precision, intent(in   ) :: Y0(Nspec+1)
           double precision, intent(inout) :: YT(Nspec+1)
           integer, intent(out), optional :: ierr

           integer          :: iwrk, iter, n
           double precision :: rwrk, rmax, rho_inv, cp

           double precision, dimension(Nspec+1) :: rhs, r, dYTdt
           double precision, dimension(Nspec)   :: Y, enthalpies

           double precision :: rcond

           logical :: diag

           integer, parameter :: max_iter = 2000

           diag = .false.

           if (.not. allocated(A)) then
              allocate(Jac(Nspec+1,Nspec+1))
              allocate(A(Nspec+1,Nspec+1))
              allocate(ipvt(Nspec+1))
           end if

           rhs = YT
           YT = Y0

           rho_inv = 1.d0/rho

           !     maximum number of iterations is max_iter
           do iter = 0, max_iter
              !     Newton's method: iteratively solve J(x_{n+1} - x_n) = -F(x_n)
              !     F(x) is given by (I - wdot - rhs)
              
              do n = 1,Nspec
                 Y(n) = YT(n)*rho_inv
              enddo
              !     compute wdot
              call CKWYR(rho, YT(Nspec+1), Y, iwrk, rwrk, dYTdt)
              !     compute enthalpies
              call CKHMS(YT(Nspec+1), iwrk, rwrk, enthalpies)
              !     compute C_p
              call CKCPBS(YT(Nspec+1), Y, iwrk, rwrk, cp)

              dYTdt(Nspec+1) = 0.d0
              do n=1,Nspec
                 !     multiply by molecular weight to get the right units
                 dYTdt(n) = dYTdt(n) * mwt(n)
                 dYtdt(Nspec+1) = dYTdt(Nspec+1) - enthalpies(n)*dYTdt(n)
              end do
              dYTdt(Nspec+1) = dYTdt(Nspec+1)*rho_inv/cp

              r = -(YT - dt*dYTdt - rhs)

              !     rmax = maxval(abs(r(1:Nspec)))
              rmax = maxval(abs(r(1:Nspec)))
              if (isnan(rmax)) then
                 print *,'wchem: backward Euler solve returned NaN'
                 print *,'wchem: iteration: ', iter
                 print *,'wchem: reciprocal condition number of J = ', rcond
                 print *,'wchem: Cp = ', cp
                 print *,'wchem: density = ', rho
                 print *,'wchem: temperature = ', YT(Nspec+1)
                 stop
              endif
              !     if we have reached the desired tolerance then we are done
              if (rmax .le. 1.d-18) then
                 exit
              endif

              !     compute the Jacobian, and then compute its LU factorization
              call LUA

              !     call LINPACK to solve the linear system Ax = r
              !     using the output from the factorization given by dgefa

              call dgesl(A, Nspec+1, Nspec+1, ipvt, r, 0)

              !     solve for the difference x = x_{n+1} - x_n
              !     so the next iterate is given by x_{n+1} = x_n + x
              YT = YT + r
           end do

           if (iter .ge. max_iter) then
              print *,'uh oh! iter=',iter
              print *,'rmax=',rmax
              stop
           endif

           if (present(ierr)) then
              if (iter .gt. 100) then
                 ierr = iter
                 YT = rhs
              else
                 ierr = 0
              end if
           end if

         contains

           !     compute the LU factorization of the Jacobian
           subroutine LUA
             integer :: i, j, iwrk     !, info
             double precision :: rwrk
             double precision, dimension(Nspec) :: C
             integer, parameter :: consP = 1
             double precision :: z(Nspec+1)

             !     compute the molar concentrations
             call CKYTCR(rho, YT(Nspec+1), Y, iwrk, rwrk, C)
             !     use the concentrarions to compute the reaction Jacobian
             call DWDOT(Jac, C, YT(Nspec+1), consP)

             !     convert to the appropriate units
             do j=1,Nspec
                do i=1,Nspec
                   Jac(i,j) = Jac(i,j) * mwt(i) * invmwt(j)
                end do

                !     last row is derivative of temperatures
                i=Nspec+1
                Jac(i,j) = Jac(i,j) * invmwt(j) ! * rho! * cp
             end do

             !     last column is derivative wrt temperature
             j = Nspec+1
             do i=1,Nspec
                Jac(i,j) = Jac(i,j) * mwt(i) ! * rho_inv ! * cpinv
             end do

             !     we are computing the Jacobian of the function (I - dt w)
             !     so the Jacobian is given by I - dt Jac
             A = -dt*Jac

             do i=1,Nspec+1
                A(i,i) = 1.d0 + A(i,i)
             end do

             !     call LINPACK to get the LU factorization of the matrix A
             !     call dgefa(A, Nspec+1, Nspec+1, ipvt, info)
             call DGECO(A, Nspec+1, Nspec+1, ipvt, rcond, z)

           end subroutine LUA

         end subroutine bechem

      end module wchem
