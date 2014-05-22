module make_plotvar_1d_module

  implicit none

contains

  subroutine make_wbar_1d(lo, hi, w,  vlo, vhi, Q, qlo, qhi)
    use variables_module
    integer, intent(in) :: lo(1), hi(1), vlo(1), vhi(1), qlo(1), qhi(1)
    double precision, intent(inout) :: w(vlo(1):vhi(1))
    double precision, intent(in   ) :: Q(qlo(1):qhi(1),nprim)

    integer :: i,iwrk
    double precision :: rwrk, Yt(nspecies)

    do i=lo(1),hi(1)
       Yt = Q(i,qy1:qy1+nspecies-1)
       call ckmmwy(Yt, iwrk, rwrk, w(i))
    enddo

  end subroutine make_wbar_1d


  subroutine make_h_1d(lo, hi, h,  vlo, vhi, Q, qlo, qhi)
    use variables_module
    integer, intent(in) :: lo(1), hi(1), vlo(1), vhi(1), qlo(1), qhi(1)
    double precision, intent(inout) :: h(vlo(1):vhi(1))
    double precision, intent(in   ) :: Q(qlo(1):qhi(1),nprim)

    integer :: i, n, qyn, qhn

    do i=lo(1),hi(1)
       h(i) = 0.d0
       do n=1, nspecies
          qhn = qh1+n-1
          qyn = qy1+n-1
          h(i) = h(i) + q(i,qyn)*q(i,qhn)
       end do
    enddo

  end subroutine make_h_1d


  subroutine make_rhoh_1d(lo, hi, rh,  vlo, vhi, Q, qlo, qhi)
    use variables_module
    integer, intent(in) :: lo(1), hi(1), vlo(1), vhi(1), qlo(1), qhi(1)
    double precision, intent(inout) :: rh(vlo(1):vhi(1))
    double precision, intent(in   ) ::  Q(qlo(1):qhi(1),nprim)

    integer :: i, n, qyn, qhn

    do i=lo(1),hi(1)
       rh(i) = 0.d0
       do n=1, nspecies
          qhn = qh1+n-1
          qyn = qy1+n-1
          rh(i) = rh(i) + q(i,qyn)*q(i,qhn)
       end do
       rh(i) = rh(i) * q(i,qrho)
    enddo

  end subroutine make_rhoh_1d


  subroutine make_cs_1d(lo, hi, cs, vlo, vhi, Q, qlo, qhi)
    use variables_module
    use chemistry_module, only : Ru
    integer, intent(in) :: lo(1), hi(1), vlo(1), vhi(1), qlo(1), qhi(1)
    double precision, intent(inout) :: cs(vlo(1):vhi(1))
    double precision, intent(in   ) ::  Q(qlo(1):qhi(1),nprim)

    integer :: i, iwrk
    double precision :: Tt, Yt(nspecies), cv, cp, Wbar, gamma, rwrk

    do i=lo(1),hi(1)
       Tt = q(i,qtemp)
       Yt = q(i,qy1:qy1+nspecies-1)
       
       call ckcvbs(Tt, Yt, iwrk, rwrk, cv)
       call ckmmwy(Yt, iwrk, rwrk, Wbar)
       
       cp = cv + Ru/Wbar
       gamma = cp / cv
       cs(i) = sqrt(gamma*q(i,qpres)/q(i,qrho))
    enddo

  end subroutine make_cs_1d


  subroutine make_magvel_1d(lo, hi, v, vlo, vhi, Q, qlo, qhi)
    use variables_module
    integer, intent(in) :: lo(1), hi(1), vlo(1), vhi(1), qlo(1), qhi(1)
    double precision, intent(inout) ::  v(vlo(1):vhi(1))
    double precision, intent(in   ) ::  Q(qlo(1):qhi(1),nprim)

    integer :: i

    do i=lo(1),hi(1)
       v(i) = abs(q(i,qu))
    enddo

  end subroutine make_magvel_1d


  subroutine make_Mach_1d(lo, hi, Ma, vlo, vhi, Q, qlo, qhi)
    use variables_module
    use chemistry_module, only : Ru
    integer, intent(in) :: lo(1), hi(1), vlo(1), vhi(1), qlo(1), qhi(1)
    double precision, intent(inout) :: Ma(vlo(1):vhi(1))
    double precision, intent(in   ) ::  Q(qlo(1):qhi(1),nprim)

    integer :: i, iwrk
    double precision :: Tt, Yt(nspecies), cv, cp, Wbar, gamma, csinv, rwrk

    do i=lo(1),hi(1)
       Tt = q(i,qtemp)
       Yt = q(i,qy1:qy1+nspecies-1)
       
       call ckcvbs(Tt, Yt, iwrk, rwrk, cv)
       call ckmmwy(Yt, iwrk, rwrk, Wbar)
       
       cp = cv + Ru/Wbar
       gamma = cp / cv
       csinv = sqrt(q(i,qrho)/(gamma*q(i,qpres)))

       Ma(i) = q(i,qu) * csinv
    enddo

  end subroutine make_Mach_1d


  subroutine make_divu_1d(lo, hi, divu, vlo, vhi, Q, qlo, qhi, dx, dlo_g, dhi_g)
    use variables_module
    use derivative_stencil_module, only : stencil_ng, first_deriv_8, first_deriv_6, &
                first_deriv_4, first_deriv_l3, first_deriv_r3, first_deriv_rb, first_deriv_lb

    integer, intent(in) :: lo(1), hi(1), vlo(1), vhi(1), qlo(1), qhi(1), dlo_g(1), dhi_g(1)
    double precision, intent(inout) :: divu(vlo(1):vhi(1))
    double precision, intent(in   ) ::    Q(qlo(1):qhi(1),nprim)
    double precision, intent(in) :: dx(1)

    integer :: i
    double precision :: dxinv(1)
    integer :: slo(1), shi(1), dlo(1), dhi(1)

    ! Only the region bounded by [dlo,dhi] contains good data.
    ! [slo,shi] will be safe for 8th-order stencil
    dlo(1) = max(lo(1)-stencil_ng, dlo_g(1))
    dhi(1) = min(hi(1)+stencil_ng, dhi_g(1))
    slo(1) = dlo(1) + stencil_ng
    shi(1) = dhi(1) - stencil_ng

    dxinv(1) = 1.0d0 / dx(1)

    do i=slo(1),shi(1)
       divu(i) = dxinv(1)*first_deriv_8(q(i-4:i+4,qu))
    enddo
          
    ! lo-x boundary
    if (dlo(1) .eq. lo(1)) then
       i = lo(1)
       ! use completely right-biased stencil
       divu(i) = dxinv(1)*first_deriv_rb(q(i:i+3,qu))
       
       i = lo(1)+1
       ! use 3rd-order slightly right-biased stencil
       divu(i) = dxinv(1)*first_deriv_r3(q(i-1:i+2,qu))
       
       i = lo(1)+2
       ! use 4th-order stencil
       divu(i) = dxinv(1)*first_deriv_4(q(i-2:i+2,qu))
       
       i = lo(1)+3
       ! use 6th-order stencil
       divu(i) = dxinv(1)*first_deriv_6(q(i-3:i+3,qu))
    end if
    
    ! hi-x boundary
    if (dhi(1) .eq. hi(1)) then
       i = hi(1)-3
       ! use 6th-order stencil
       divu(i) = dxinv(1)*first_deriv_6(q(i-3:i+3,qu))
       
       i = hi(1)-2
       ! use 4th-order stencil
       divu(i) = dxinv(1)*first_deriv_4(q(i-2:i+2,qu))
       
       i = hi(1)-1
       ! use 3rd-order slightly left-biased stencil
       divu(i) = dxinv(1)*first_deriv_l3(q(i-2:i+1,qu))
       
       i = hi(1)
       ! use completely left-biased stencil
       divu(i) = dxinv(1)*first_deriv_lb(q(i-3:i,qu))
    end if

  end subroutine make_divu_1d


  subroutine make_burn_1d(lo, hi, burn, vlo, vhi, Q, qlo, qhi)
    use plotvar_index_module
    use variables_module
    use chemistry_module, only : molecular_weight, h0 => std_heat_formation

    integer, intent(in) :: lo(1), hi(1), vlo(1), vhi(1), qlo(1), qhi(1)
    double precision, intent(inout) :: burn(vlo(1):vhi(1),nburn)
    double precision, intent(in   ) ::    Q(qlo(1):qhi(1),nprim)

    integer :: i,n,np,iwrk
    double precision :: Yt(lo(1):hi(1),nspecies), wdot(lo(1):hi(1),nspecies), rwrk

    np = hi(1) - lo(1) + 1

    do n=1, nspecies
       do i=lo(1),hi(1)
          Yt(i,n) = q(i,qy1+n-1)
       end do
    end do
       
    call vckwyr(np, q(lo(1),qrho), q(lo(1),qtemp), Yt, iwrk, rwrk, wdot)
       
    if (ib_omegadot > 0) then
       do n=1, nspecies
          do i=lo(1),hi(1)
             burn(i,ib_omegadot+n-1) = wdot(i,n) * molecular_weight(n)
          end do
       end do
    end if
    
    if (ib_dYdt > 0) then
       do n=1, nspecies
          do i=lo(1),hi(1)
             burn(i,ib_dYdt+n-1) = wdot(i,n) * molecular_weight(n) / q(i,qrho)
          end do
       end do
    end if
    
    if (ib_heatRelease > 0) then
       do i=lo(1),hi(1)
          burn(i,ib_heatRelease) = 0.d0
       end do
       
       do n=1,nspecies
          do i=lo(1),hi(1)
             burn(i,ib_heatRelease) = burn(i,ib_heatRelease) &
                  - h0(n) * wdot(i,n) * molecular_weight(n)
          end do
       end do
    end if
    
    if (ib_fuelConsumption > 0) then
       do i=lo(1),hi(1)
          burn(i,ib_fuelConsumption) = -wdot(i,ifuel) * molecular_weight(ifuel)
       end do
    end if

  end subroutine make_burn_1d


  subroutine make_burn2_1d(lo, hi, burn, vlo, vhi, u0, u0lo, u0hi, u, ulo, uhi, dt)
    use plotvar_index_module
    use variables_module
    use chemistry_module, only : h0 => std_heat_formation

    integer, intent(in) :: lo(1),hi(1),vlo(1),vhi(1),u0lo(1),u0hi(1),ulo(1),uhi(1)
    double precision, intent(inout) :: burn( vlo(1): vhi(1),nburn)
    double precision, intent(in   ) ::   u0(u0lo(1):u0hi(1),ncons)
    double precision, intent(in   ) ::   u ( ulo(1): uhi(1),ncons)
    double precision, intent(in) :: dt

    integer :: i, n, iryn, ibn
    double precision :: dtinv

    dtinv = 1.d0/dt

    if (ib_omegadot > 0) then
       do n=1,nspecies
          iryn = iry1+n-1
          ibn  = ib_omegadot+n-1
          do i=lo(1),hi(1)
             burn(i,ibn) = (u(i,iryn)-u0(i,iryn))*dtinv
          end do
       end do
    end if

    if (ib_dYdt > 0) then
       do n=1,nspecies
          iryn = iry1+n-1
          ibn  = ib_dYdt+n-1
          do i=lo(1),hi(1)
             burn(i,ibn) = ( u(i,iryn)/ u(i,irho)  &
                  -           u0(i,iryn)/u0(i,irho) )*dtinv
          end do
       end do       
    end if

    if (ib_heatRelease > 0) then
       do i=lo(1),hi(1)
          burn(i,ib_heatRelease) = 0.d0
       end do
       
       do n=1,nspecies
          iryn = iry1+n-1
          do i=lo(1),hi(1)
             burn(i,ib_heatRelease) = burn(i,ib_heatRelease) &
                  - h0(n) * (u(i,iryn)-u0(i,iryn))*dtinv
          end do
       end do
    end if

    if (ib_fuelConsumption > 0) then
       iryn = iry1+ifuel-1
       do i=lo(1),hi(1)
          burn(i,ib_fuelConsumption) = -(u(i,iryn)-u0(i,iryn))*dtinv
       end do
    end if

  end subroutine make_burn2_1d

end module make_plotvar_1d_module
