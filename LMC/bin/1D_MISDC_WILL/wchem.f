      ! Will's Backward Euler chemistry solve routine
      module wchem
       
        implicit none
        
        private
        
        include 'spec.h'
        
        ! Jacobian matrix and factorization
        double precision, allocatable, save :: Jac(:,:), A(:,:)
        ! Pivot vector returned from LINPACK
        integer, allocatable, save :: ipvt(:)
        
        public :: bechem

      contains
      
        ! do a Backward Euler solve for the chemistry using Newton's method
        ! Y0 : initial guess
        ! Qrh: source term for enthalpy equation -- appears in temperature eqn
        ! Qry: source term for the species equation
        ! rho: density
        ! YT : input: right-hand side
        !      output: solution
        ! dt
        ! ierr : (optional) error
        subroutine bechem(Y0, Qrh, Qry, rho, YT, dt, ierr)
          double precision, intent(in   ) :: dt
          double precision, intent(in   ) :: Qrh
          double precision, intent(in   ) :: Qry(Nspec)
          double precision, intent(in   ) :: rho
          double precision, intent(in   ) :: Y0(Nspec+1)
          double precision, intent(inout) :: YT(Nspec+1)
          integer, intent(out), optional :: ierr
          
          integer          :: iwrk, iter, n
          double precision :: rwrk, rmax, rho_new, rho_inv, cp
          
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
          
          ! maximum number of iterations is max_iter
          do iter = 0, max_iter
             ! Newton's method: iteratively solve J(x_{n+1} - x_n) = -F(x_n)
             ! F(x) is given by (I - wdot - rhs)
             
             rho_new = 0.d0
             do n = 1,Nspec
                rho_new = rho_new + YT(n)
             enddo
             
             do n = 1,Nspec
                Y(n) = YT(n)*rho_inv
             enddo
             
             ! compute wdot
             call CKWYR(rho, YT(Nspec+1), Y, iwrk, rwrk, dYTdt)
             ! compute enthalpies
             call CKHMS(YT(Nspec+1), iwrk, rwrk, enthalpies)
             call CKCPBS(YT(Nspec+1), Y, iwrk, rwrk, cp)
             
             dYTdt(Nspec+1) = 0.d0
             do n=1,Nspec
                ! multiply by molecular weight to get the right units
                dYTdt(n) = dYTdt(n) * mwt(n)
                dYTdt(Nspec+1) = dYTdt(Nspec+1) + enthalpies(n)*(Qry(n) + dYTdt(n))
             end do
             dYTdt(Nspec+1) = (Qrh - dYTdt(Nspec+1))*rho_inv/cp
             
             r = -(YT - dt*dYTdt - rhs)
             
             !rmax = maxval(abs(r(1:Nspec)))
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
             ! if we have reached the desired tolerance then we are done
             if (rmax .le. 1.d-18) then
                exit
             endif
             
             ! compute the Jacobian, and then compute its LU factorization
             call LUA
             
             ! call LINPACK to solve the linear system Ax = r
             ! using the output from the factorization given by dgefa
             
             call dgesl(A, Nspec+1, Nspec+1, ipvt, r, 0)
             
             ! solve for the difference x = x_{n+1} - x_n
             ! so the next iterate is given by x_{n+1} = x_n + x
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
          
          ! compute the LU factorization of the Jacobian
          subroutine LUA
            integer :: i, j, iwrk!, info
            double precision :: rwrk
            double precision, dimension(Nspec) :: C
            integer, parameter :: consP = 1
            double precision :: z(Nspec+1)
            
            ! compute the molar concentrations
            call CKYTCR(rho, YT(Nspec+1), Y, iwrk, rwrk, C)
            ! use the concentrarions to compute the reaction Jacobian
            call DWDOT(Jac, C, YT(Nspec+1), consP)
            
            ! convert to the appropriate units
            do j=1,Nspec
               do i=1,Nspec
                  Jac(i,j) = Jac(i,j) * mwt(i) * invmwt(j)
               end do
               
               ! last row is derivative of temperatures
               i=Nspec+1
               Jac(i,j) = Jac(i,j) * invmwt(j)! * rho! * cp
            end do
            
            ! last column is derivative wrt temperature
            j = Nspec+1
            do i=1,Nspec
               Jac(i,j) = Jac(i,j) * mwt(i)! * rho_inv ! * cpinv
            end do
            
            ! we are computing the Jacobian of the function (I - dt w)
            ! so the Jacobian is given by I - dt Jac
            A = -dt*Jac
            
            do i=1,Nspec+1
               A(i,i) = 1.d0 + A(i,i)
            end do
            
            ! call LINPACK to get the LU factorization of the matrix A
            !call dgefa(A, Nspec+1, Nspec+1, ipvt, info)
            call DGECO(A, Nspec+1, Nspec+1, ipvt, rcond, z)
            
          end subroutine LUA

        end subroutine bechem

      end module wchem
