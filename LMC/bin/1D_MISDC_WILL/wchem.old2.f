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
        ! rho: density
        ! YT : input: right-hand side
        !      output: solution
        ! dt
        ! ierr : (optional) error
        subroutine bechem(Y0, T0, rho, YT, dt, ierr)
          double precision, intent(in   ) :: rho, dt
          double precision, intent(in   ) :: Y0(Nspec+1)
          double precision, intent(inout) :: YT(Nspec+1)
          double precision, intent(in   ) :: T0
          integer, intent(out), optional :: ierr
          
          integer :: iwrk, iter, n
          double precision :: rwrk, rhoinv, rmax, Tnew, rho_new, hmix
          double precision, dimension(Nspec+1) :: rhs, r, dYTdt
          
          double precision :: Y(Nspec)
          
          integer NiterMAX, Niter
          parameter (NiterMAX = 30)
          double precision res(NiterMAX), errMAX
          
          errMax = hmix_TYP*1.e-20
          
          if (.not. allocated(A)) then
             allocate(Jac(Nspec+1,Nspec+1))
             allocate(A(Nspec+1,Nspec+1))
             allocate(ipvt(Nspec+1))
          end if

          rhoinv = 1.d0/rho

          rhs = YT
          YT = Y0
          Tnew = T0
          
          ! maximum number of iterations is 100
          do iter = 0, 100
             ! Newton's method: iteratively solve J(x_{n+1} - x_n) = -F(x_n)
             ! F(x) is given by (I - wdot - rhs)
             
             rho_new = 0.d0
             do n = 1,Nspec
                rho_new = rho_new + YT(n)
             enddo
             
             do n = 1,Nspec
                Y(n) = YT(n)/rho_new
             enddo
             hmix = YT(Nspec + 1)/rho_new
             
             print *,'enthalpy hmix:'
             print *,'--------------'
             print *,hmix
             print *,''
             
             ! get the temperature
             call FORT_TfromHYpt(Tnew,hmix,Y,Nspec,errMax,NiterMAX,res,Niter)
             
             if (Niter.lt.0) then
               print *,'wchem: H to T solve failed'
               print *,'iter:', iter
               print *,'hmix:', hmix
               print *,'errMax:', errMax
               print *,'Niter=',Niter
               stop
            endif
             
             print *,'mass fractions Y:'
             print *,'-----------------'
             write(*,'(100G20.6E3)') Y(:)
             print *,''
             print *,'temperature T:'
             print *,'--------------'
             print *,Tnew
             print *,''
             print *,'density rho:'
             print *,'--------------'
             print *,rho
             print *,''
             
             ! compute wdot
             call CKWYR(rho, Tnew, Y, iwrk, rwrk, dYTdt)
             ! enthalpy term has no contribution
             dYTdt(Nspec+1) = 0.d0
             
             print *,'production rate wdot:'
             print *,'---------------------'
             write(*,'(100G20.6E3)') dYTdt(:)
             print *,''
             
             ! multiply by molecular weight to get the right units
             do n=1,Nspec
                dYTdt(n) = dYTdt(n) * mwt(n)! * rhoinv
c                dYTdt(Nspec+1) = dYTdt(Nspec+1) - dYTdt(n)*uk(n)
             end do

             !r = YT - dt * dYTdt - YT_init
             r = rhs - YT + dt*dYTdt
             
             !print *,'r = ', r
             
             rmax = maxval(abs(r(1:Nspec)))
             print *,'rmax = ', rmax
             ! if we have reached the desired tolerance (1e-14) then we are done
             !if (rmax .le. 1.d-14) then
             if (rmax .le. 1.d-12) then
                exit 
             endif
             
             ! compute the Jacobian, and then compute its LU factorization
             call LUA
             
             ! call LINPACK to solve the linear system Ax = r
             ! using the output from the factorization given by dgefa
             call dgesl(A, Nspec+1, Nspec+1, ipvt, r, 0)
             
             print *,'solution r:'
             print *,'---------'
             write(*,'(100G20.6E3)') r(:)
             print *,''
             
             ! solve for the difference x = x_{n+1} - x_n
             ! so the next iterate is given by x_{n+1} = x_n + x
             YT = YT + r
          end do

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
            integer :: i, j, iwrk, info
            double precision :: rwrk
            double precision :: cp, cpinv
            double precision, dimension(Nspec) :: C
            integer, parameter :: consP = 1
            
            ! call the chemistry package to get the analytical form of the derivatives
            call CKYTCR(rho, Tnew, Y, iwrk, rwrk, C)
            call CKCPBS(Tnew, Y, iwrk, rwrk, cp)
            call DWDOT(Jac, C, Tnew, consP)
            
            cpinv = 1.d0/cp
            
            ! convert to the appropriate units
            do j=1,Nspec
               do i=1,Nspec
                  Jac(i,j) = Jac(i,j) * mwt(i) * invmwt(j)
               end do
               
               ! last row is derivative of temperatures
               ! need to convert to derivative of enthalpy
               i=Nspec+1
               Jac(i,j) = Jac(i,j) * invmwt(j) * cp ! * rho
            end do
            
            ! last column is derivative wrt temperature
            ! need to convert to derivative wrt enthalpy
            j = Nspec+1
            do i=1,Nspec
               Jac(i,j) = Jac(i,j) * mwt(i) * cpinv ! * rhoinv
            end do
            
            ! we are computing the Jacobian of the function (I - dt w)
            ! so the Jacobian is given by I - dt Jac
            A = -dt*Jac
            
            do i=1,Nspec+1
               A(i,i) = 1.d0 + A(i,i)
            end do
            
            print *,'Jac matrix:'
            print *,'---------'
            do i=1,size(A,1)
               write(*,'(100G20.6E3)') Jac(i,:)
            end do 
            print *,''
            
            ! call LINPACK to get the LU factorization of the matrix A
            call dgefa(A, Nspec+1, Nspec+1, ipvt, info)
            
            print *,'A matrix:'
            print *,'---------'
            do i=1,size(A,1)
               write(*,'(100G20.6E3)') A(i,:)
            end do 
            print *,''
            
          end subroutine LUA

        end subroutine bechem

      end module wchem
