module chemsolv_module

  implicit none

  logical, save :: stiff = .true.

  private

  public :: chemsolv

contains

  subroutine chemsolv(loc, hic, u, u_l1, u_l2, u_l3, u_h1, u_h2, u_h3, dt)
    
    use meth_params_module, only : NVAR, URHO, UEDEN, UMX, UMY, UMZ, UTEMP, UFS
    use chemistry_module, only : nspec=>nspecies
    use vode_module, only : itol, rtol, atol, MF_NOSTIFF, MF_NUMERICAL_JAC, &
         voderwork, vodeiwork, lvoderwork, lvodeiwork, voderpar, vodeipar

    implicit none

    integer, intent(in) :: loc(3), hic(3)
    integer, intent(in) :: u_l1, u_l2, u_l3, u_h1, u_h2, u_h3
    double precision, intent(in) :: dt
    double precision, intent(inout) :: u(u_l1:u_h1,u_l2:u_h2,u_l3:u_h3,NVAR)

    external jac, f_rhs, dvode

    ! vode stuff
    integer, parameter :: itask=1, iopt=1
    integer :: MF, istate

    integer :: i, j, k, ckiwork
    double precision :: time, Y(nspec), rho, rhoinv, ei, ckrwork

    if (stiff) then
       MF = MF_NUMERICAL_JAC
    else
       MF = MF_NOSTIFF
    end if

    do k=loc(3),hic(3)
    do j=loc(2),hic(2)
    do i=loc(1),hic(1)

       rho = u(i,j,k,URHO)
       rhoinv = 1.d0/rho

       ei = rhoinv*( u(i,j,k,UEDEN) - 0.5d0*rhoinv* &
            (u(i,j,k,UMX)**2 + u(i,j,k,UMY)**2 + u(i,j,k,UMZ)**2) )

       voderpar(1) = rho
       voderpar(2) = ei
       
       Y = u(i,j,k,UFS:UFS+nspec-1) * rhoinv

       istate = 1
       time = 0.d0

       call dvode(f_rhs, nspec, Y, time, dt, itol, rtol, atol, itask, &
            istate, iopt, voderwork, lvoderwork, vodeiwork, lvodeiwork, &
            jac, MF, voderpar, vodeipar)

       if (istate < 0) then
          print *, 'chemsolv: VODE failed'
          print *, 'istate = ', istate, ' time =', time
          call bl_error("ERROR in chemsolv: VODE failed")
       end if

       call feeytt(ei, Y, ckiwork, ckrwork, u(i,j,k,UTEMP))

       u(i,j,k,UFS:UFS+nspec-1) = rho * Y
    end do
    end do
    end do

  end subroutine chemsolv

end module chemsolv_module


