module eigen_module
  implicit none

  private

  public :: get_eigen_matrices

contains

  subroutine get_eigen_matrices(rho, Y, T, vel, lv,rv)
    use meth_params_module, only : NSPEC, NCHARV, CFS
    use eos_module, only : eos_given_RTY, eos_get_eref

    double precision, intent(in) :: rho, Y(NSPEC), T, vel(3)
    double precision, intent(out) :: lv(NCHARV,NCHARV), rv(NCHARV,NCHARV)
    
    integer :: m, n
    double precision :: rhoInv, p, c, gamc, dpdr(NSPEC), dpde, e, ek, H
    double precision :: gt, b, d(NSPEC), gtinv, cinv
    double precision :: eref

    rhoInv = 1.d0/rho
    
    ek = 0.5d0*(vel(1)**2 + vel(2)**2 + vel(3)**2)

    call eos_given_RTY(e, p, c, gamc, dpdr, dpde, rho, T, Y)
    
    eref = eos_get_eref(Y)
    e = e - eref
    
    H = e + p*rhoInv + ek
    
    gt = dpde*rhoInv
    cinv = 1.d0/c
    b = gt*cinv*cinv
    gtinv = 1.d0/gt
    do n=1,nspec
       d(n) = b*(ek - e + dpdr(n)*gtinv)
    end do
    
    ! assemble left vectors
    lv(1,1) = -0.5d0*(cinv + b*vel(1))
    lv(2,1) = -0.5d0*b*vel(2)
    lv(3,1) = -0.5d0*b*vel(3)
    lv(4,1) =  0.5d0*b
    do n=1,nspec
       lv(CFS+n-1,1) = 0.5d0*(vel(1)*cinv + d(n))
    end do
    
    lv(1,2) =  0.5d0*(cinv - b*vel(1))
    lv(2,2) = -0.5d0*b*vel(2)
    lv(3,2) = -0.5d0*b*vel(3)
    lv(4,2) =  0.5d0*b
    do n=1,nspec
       lv(CFS+n-1,2) = 0.5d0*(-vel(1)*cinv + d(n))
    end do
    
    lv(1,3) = 0.d0
    lv(2,3) = 1.d0
    lv(3,3) = 0.d0
    lv(4,3) = 0.d0
    do n=1,nspec
       lv(CFS+n-1,3) = -vel(2)
    end do
    
    lv(1,4) = 0.d0
    lv(2,4) = 0.d0
    lv(3,4) = 1.d0
    lv(4,4) = 0.d0
    do n=1,nspec
       lv(CFS+n-1,4) = -vel(3)
    end do
    
    do m=1,nspec
       lv(1,CFS+m-1) =  Y(m)*b*vel(1)
       lv(2,CFS+m-1) =  Y(m)*b*vel(2)
       lv(3,CFS+m-1) =  Y(m)*b*vel(3)
       lv(4,CFS+m-1) = -Y(m)*b
       do n=1,nspec
          lv(CFS+n-1,CFS+m-1) = -Y(m)*d(n)
       end do
       lv(CFS+m-1,CFS+m-1) = lv(CFS+m-1,CFS+m-1) + 1.d0
    end do
    
    ! assemble right vectors
    rv(1,1) = vel(1) - c
    rv(2,1) = vel(2)
    rv(3,1) = vel(3)
    rv(4,1) = H - vel(1)*c
    do n=1,nspec
       rv(CFS+n-1,1) = Y(n)
    end do
    
    rv(1,2) = vel(1) + c
    rv(2,2) = vel(2)
    rv(3,2) = vel(3)
    rv(4,2) = H + vel(1)*c
    do n=1,nspec
       rv(CFS+n-1,2) = Y(n)
    end do
    
    rv(1,3) = 0.d0
    rv(2,3) = 1.d0
    rv(3,3) = 0.d0
    rv(4,3) = vel(2)
    do n=1,nspec
       rv(CFS+n-1,3) = 0.d0
    end do
    
    rv(1,4) = 0.d0
    rv(2,4) = 0.d0
    rv(3,4) = 1.d0
    rv(4,4) = vel(3)
    do n=1,nspec
       rv(CFS+n-1,4) = 0.d0
    end do
    
    do n=1,nspec
       rv(1,CFS+n-1) = vel(1)
       rv(2,CFS+n-1) = vel(2)
       rv(3,CFS+n-1) = vel(3)
       rv(4,CFS+n-1) = e + ek - dpdr(n)*gtinv
       do m=1,nspec
          rv(CFS+m-1,CFS+n-1) = 0.d0
       end do
       rv(CFS+n-1,CFS+n-1) = 1.d0
    end do

  end subroutine get_eigen_matrices

end module eigen_module
