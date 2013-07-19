  subroutine update_inlet_ylo_2d(lo,hi,qin,t,dx)

    use bl_constants_module, only : Pi=>M_PI
    use DME_jet_module
    use probin_module, only : prob_lo, prob_hi, inflow_period, inflow_vnmag, &
         xfrontw, splitx, vn_in, vn_co

    integer, intent(in) :: lo(2), hi(2)
    double precision, intent(in) :: t, dx(2)
    double precision, intent(inout) :: qin(nqin,lo(1):hi(1))

    integer :: i
    double precision :: x, vn0, facx, fact, sigma, eta
    
    facx = 2.d0*Pi/(prob_hi(1)-prob_lo(1))
    fact = sin(2.d0*Pi*t/inflow_period)

    sigma = 2.5d0*xfrontw*splitx

    do i=lo(1), hi(1)
       x = prob_lo(1) + dx(1)*(i + 0.5d0)

       eta = 0.5d0 * (tanh((x + splitx)/sigma)   &
            &       - tanh((x - splitx)/sigma))
       
       vn0 = eta *vn_in + (1.d0-eta) *vn_co

       ! sinusoidal variation of inflow
       qin(ivin,i) = vn0 + inflow_vnmag*eta*sin(x*facx)*fact
       
    end do

  end subroutine update_inlet_ylo_2d
