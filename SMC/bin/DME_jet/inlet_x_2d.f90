  subroutine update_inlet_xlo_2d(lo,hi,qin,t,dx)

    use bl_constants_module, only : Pi=>M_PI
    use DME_jet_module
    use probin_module, only : prob_lo, prob_hi, inflow_period, inflow_vnmag, &
         yfrontw, splity, vn_in, vn_co

    integer, intent(in) :: lo(2), hi(2)
    double precision, intent(in) :: t, dx(2)
    double precision, intent(inout) :: qin(nqin,lo(2):hi(2))

    integer :: j
    double precision :: y, vn0, facx, fact, sigma, eta
    
    facx = 2.d0*Pi/(prob_hi(2)-prob_lo(2))
    fact = sin(2.d0*Pi*t/inflow_period)

    sigma = 2.5d0*yfrontw*splity

    do j=lo(2), hi(2)
       y = prob_lo(2) + dx(2)*(j + 0.5d0)

       eta = 0.5d0 * (tanh((y + splity)/sigma)   &
            &       - tanh((y - splity)/sigma))
       
       vn0 = eta *vn_in + (1.d0-eta) *vn_co

       ! sinusoidal variation of inflow
       qin(iuin,j) = vn0 + inflow_vnmag*eta*sin(y*facx)*fact
       
    end do

  end subroutine update_inlet_xlo_2d

