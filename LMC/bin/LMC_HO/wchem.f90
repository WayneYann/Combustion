module wchem_module
   
   implicit none
   private
   include 'spec.h'

   !     Jacobian matrix and factorization
   double precision, allocatable, save :: Jac(:,:), A(:,:)
   !     Pivot vector returned from LINPACK
   integer, allocatable, save :: ipvt(:)

   public :: bechem

contains

   !     do a Backward Euler solve for the chemistry using Newton's method
   !      Y : solution (output)
   !     Y0 : initial guess
   !    rho : density
   !   hmix : enthalpy
   !    rhs : input: right-hand side
   !     dt : timestep
   subroutine bechem(Y, Y0, rho, hmix, rhs, dt)
     double precision, intent(out) :: Y(Nspec)
     double precision, intent(in ) :: Y0(Nspec)
     double precision, intent(in ) :: rho
     double precision, intent(in ) :: hmix
     double precision, intent(in ) :: rhs(Nspec)
     double precision, intent(in ) :: dt
     
     integer          :: iwrk, iter, n
     double precision :: rwrk, rmax, rho_inv, cp, cp_inv
     double precision :: T
     
     double precision, dimension(Nspec) :: r, wdot
     
     double precision :: rcond
     
     integer, parameter :: max_iter = 100, NiterMAX = 40
     integer :: Niter
     double precision, parameter :: tol = 1.d-14
     double precision res(NiterMAX), errMAX
     
     errMax = hmix_TYP*1.e-21
   
     if (.not. allocated(A)) then
        allocate(Jac(Nspec+1,Nspec+1))
        allocate(A(Nspec,Nspec))
        allocate(ipvt(Nspec))
     end if
     
     rho_inv = 1.d0/rho
     
     ! start with the initial guess
     Y = Y0
     
     ! maximum number of iterations is max_iter
     do iter = 0, max_iter
        ! Newton's method: iteratively solve J(x_{n+1} - x_n) = -F(x_n)
        ! F(x) is given by (I - dt*wdot/rho - rhs/rho)
        
        ! get the temperature
        T = T_INIT
        call FORT_TfromHYpt(T,hmix,Y,Nspec,errMax,NiterMAX,res,Niter)
        T_INIT = T
        
        if (Niter.lt.0) then
          print *,'bechem: H to T solve failed'
          print *,'Niter =', Niter
          stop
        end if
        
        ! compute wdot
        call CKWYR(rho, T, Y, iwrk, rwrk, wdot)
        ! compute C_p and 1/C_p
        !call CKCPBS(T, Y, iwrk, rwrk, cp)
        call compute_cp(cp, T, Y)
        cp_inv = 1.d0/cp
        
        ! multiply by molecular weight to get the right units
        do n=1,Nspec
           wdot(n) = wdot(n) * mwt(n)
        end do
        r = -(Y - dt*wdot*rho_inv - rhs)

        rmax = maxval(abs(r))
        if (isnan(rmax)) then
           print *,'wchem: backward Euler solve returned NaN'
           print *,'wchem: iteration: ', iter
           print *,'wchem: reciprocal condition number of J = ', rcond
           print *,'wchem: Cp = ', cp
           print *,'wchem: density = ', rho
           print *,'wchem: temperature = ', T
           
           print *,'wchem: Y0:',Y0
           stop
        endif
        ! if we have reached the desired tolerance then we are done
        if (rmax .le. tol) then
           if (iter .gt. 10) then
           print *,'iters=',iter
           end if
           exit
        endif

        !     compute the Jacobian, and then compute its LU factorization
        call LUA
        !     call LINPACK to solve the linear system Ax = r
        !     using the output from the factorization given by dgefa
        call dgesl(A, Nspec, Nspec, ipvt, r, 0)
        
        !     solve for the difference x = x_{n+1} - x_n
        !     so the next iterate is given by x_{n+1} = x_n + x
        Y = Y + r
     end do

     if (iter .ge. max_iter) then
        print *,'wchem: Newton solve failed to converge'
        print *,'wchem: iter=',iter
        print *,'wchem: rmax=',rmax
        stop
     endif

   contains

     !     compute the LU factorization of the Jacobian
     subroutine LUA
       integer :: i, j, iwrk     !, info
       double precision :: rwrk
       double precision, dimension(Nspec) :: C
       integer, parameter :: consP = 1
       double precision :: z(Nspec)

       !    compute the molar concentrations
       call CKYTCR(rho, T, Y, iwrk, rwrk, C)
       !     use the concentrarions to compute the reaction Jacobian
       call DWDOT(Jac, C, T, consP)

       !     convert to the appropriate units
       do j=1,Nspec
          do i=1,Nspec
             Jac(i,j) = Jac(i,j) * mwt(i) * invmwt(j)
          end do
       end do

       !     we are computing the Jacobian of the function (I - dt w)
       !     so the Jacobian is given by I - dt Jac
       do j=1,Nspec
         do i=1,Nspec
            A(i,j) = -dt*Jac(i,j)
         end do
       end do

       do i=1,Nspec
          A(i,i) = 1.d0 + A(i,i)
       end do

       !     call LINPACK to get the LU factorization of the matrix A
       !     call dgefa(A, Nspec+1, Nspec+1, ipvt, info)
       call DGECO(A, Nspec, Nspec, ipvt, rcond, z)

     end subroutine LUA

   end subroutine bechem

end module wchem_module
