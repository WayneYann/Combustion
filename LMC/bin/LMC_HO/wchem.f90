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
   !     Y0 : initial guess
   !     rho: density
   !     YT : input: right-hand side
   !          output: solution
   !     dt : timestep
   !     ierr : (optional) error
   subroutine bechem(Y0, rho, rYh, dt, ierr)
     double precision, intent(in   ) :: dt
     double precision, intent(in   ) :: rho
     double precision, intent(in   ) :: Y0(Nspec+1)
     double precision, intent(inout) :: rYh(Nspec+1)
     integer, intent(out), optional  :: ierr

     integer          :: iwrk, iter, n
     double precision :: rwrk, rmax, rho_inv, cp, cp_inv, rho_new
     double precision :: T, hmix

     double precision, dimension(Nspec+1) :: rhs, r, wdot
     double precision, dimension(Nspec)   :: Y

     double precision :: rcond
     
     integer, parameter :: max_iter = 100, NiterMAX = 40
     integer :: Niter
     double precision, parameter :: tol = 1.d-18
     double precision res(NiterMAX), errMAX
   
     errMax = hmix_TYP*1.e-21
   
     if (.not. allocated(A)) then
        allocate(Jac(Nspec+1,Nspec+1))
        allocate(A(Nspec+1,Nspec+1))
        allocate(ipvt(Nspec+1))
     end if
     
     rhs = rYh
     rYH = Y0
     
     !     maximum number of iterations is max_iter
     do iter = 0, max_iter
        !     Newton's method: iteratively solve J(x_{n+1} - x_n) = -F(x_n)
        !     F(x) is given by (I - wdot - rhs)
        
        rho_new = 0.d0
        do n = 1,Nspec
          rho_new = rho_new + rYh(n)
        end do
        rho_inv = 1.d0/rho_new
        
        do n = 1,Nspec
           Y(n) = rYh(n)*rho_inv
        enddo
        hmix = rYh(Nspec+1)*rho_inv
        
        T = T_INIT
        call FORT_TfromHYpt(T,hmix,Y,Nspec,errMax,NiterMAX,res,Niter)
        T_INIT = T
        
        if (Niter.lt.0) then
          print *,'bechem: H to T solve failed'
          print *,'Niter =', Niter
          stop
        end if
        
        !     compute wdot
        call CKWYR(rho_new, T, Y, iwrk, rwrk, wdot)
        !     compute C_p
        call CKCPBS(T, Y, iwrk, rwrk, cp)
        
        cp_inv = 1.d0/cp
        
        wdot(Nspec+1) = 0
        do n=1,Nspec
           !     multiply by molecular weight to get the right units
           wdot(n) = wdot(n) * mwt(n)
        end do
        r = -(rYh - rhs - dt*wdot)

        !     rmax = maxval(abs(r(1:Nspec)))
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
        !     if we have reached the desired tolerance then we are done
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
        call dgesl(A, Nspec+1, Nspec+1, ipvt, r, 0)

        !     solve for the difference x = x_{n+1} - x_n
        !     so the next iterate is given by x_{n+1} = x_n + x
        rYh = rYh + r
     end do

     if (iter .ge. max_iter) then
        print *,'wchem: Newton solve failed to converge'
        print *,'wchem: iter=',iter
        print *,'wchem: rmax=',rmax
        stop
     endif

     if (present(ierr)) then
        if (iter .gt. 100) then
           ierr = iter
           rYh = rhs
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

       !    compute the molar concentrations
       call CKYTCR(rho_new, T, Y, iwrk, rwrk, C)
       !     use the concentrarions to compute the reaction Jacobian
       call DWDOT(Jac, C, T, consP)

       !     convert to the appropriate units
       do j=1,Nspec
          do i=1,Nspec
             Jac(i,j) = Jac(i,j) * mwt(i) * invmwt(j)
          end do

          !     last row is derivative of enthalpy
          i=Nspec+1
          Jac(i, j) = 0.d0
       end do
       Jac(Nspec+1, Nspec+1) = 0.d0
       
       !     last column is derivative wrt enthalpy
       j = Nspec+1
       do i=1,Nspec
          Jac(i,j) = Jac(i,j) * mwt(i) * rho_inv * cp_inv
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

end module wchem_module
