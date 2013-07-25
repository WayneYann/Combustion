module make_plotvar_2d_module

  implicit none

contains

  subroutine make_h_2d(lo, hi, h,  vlo, vhi, Q, qlo, qhi)
    use variables_module
    integer, intent(in) :: lo(2), hi(2), vlo(2), vhi(2), qlo(2), qhi(2)
    double precision, intent(inout) :: h(vlo(1):vhi(1),vlo(2):vhi(2))
    double precision, intent(in   ) :: Q(qlo(1):qhi(1),qlo(2):qhi(2),nprim)

    integer :: i,j, n, qyn, qhn

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          h(i,j) = 0.d0
          do n=1, nspecies
             qhn = qh1+n-1
             qyn = qy1+n-1
             h(i,j) = h(i,j) + q(i,j,qyn)*q(i,j,qhn)
          end do
       enddo
    enddo

  end subroutine make_h_2d


  subroutine make_divu_2d(lo, hi, divu, vlo, vhi, Q, qlo, qhi, dx, dlo_g, dhi_g)
    use variables_module
    use derivative_stencil_module, only : stencil_ng, first_deriv_8, first_deriv_6, &
                first_deriv_4, first_deriv_l3, first_deriv_r3, first_deriv_rb, first_deriv_lb

    integer, intent(in) :: lo(2), hi(2), vlo(2), vhi(2), qlo(2), qhi(2), dlo_g(2), dhi_g(2)
    double precision, intent(inout) :: divu(vlo(1):vhi(1),vlo(2):vhi(2))
    double precision, intent(in   ) ::    Q(qlo(1):qhi(1),qlo(2):qhi(2),nprim)
    double precision, intent(in) :: dx(2)

    integer :: i,j
    double precision, allocatable :: ux(:,:), vy(:,:)
    double precision :: dxinv(2)
    integer :: slo(2), shi(2), dlo(2), dhi(2)

    ! Only the region bounded by [dlo,dhi] contains good data.
    ! [slo,shi] will be safe for 8th-order stencil
    do i=1,2
       dlo(i) = max(lo(i)-stencil_ng, dlo_g(i))
       dhi(i) = min(hi(i)+stencil_ng, dhi_g(i))
       slo(i) = dlo(i) + stencil_ng
       shi(i) = dhi(i) - stencil_ng
    end do

    do i=1,2
       dxinv(i) = 1.0d0 / dx(i)
    end do

    allocate(ux(lo(1):hi(1),lo(2):hi(2)))
    allocate(vy(lo(1):hi(1),lo(2):hi(2)))

    do j=lo(2),hi(2)

       do i=slo(1),shi(1)
          ux(i,j) = dxinv(1)*first_deriv_8(q(i-4:i+4,j,qu))
       enddo
          
       ! lo-x boundary
       if (dlo(1) .eq. lo(1)) then
          i = lo(1)
          ! use completely right-biased stencil
          ux(i,j) = dxinv(1)*first_deriv_rb(q(i:i+3,j,qu))
          
          i = lo(1)+1
          ! use 3rd-order slightly right-biased stencil
          ux(i,j) = dxinv(1)*first_deriv_r3(q(i-1:i+2,j,qu))
          
          i = lo(1)+2
          ! use 4th-order stencil
          ux(i,j) = dxinv(1)*first_deriv_4(q(i-2:i+2,j,qu))
          
          i = lo(1)+3
          ! use 6th-order stencil
          ux(i,j) = dxinv(1)*first_deriv_6(q(i-3:i+3,j,qu))
       end if
       
       ! hi-x boundary
       if (dhi(1) .eq. hi(1)) then
          i = hi(1)-3
          ! use 6th-order stencil
          ux(i,j) = dxinv(1)*first_deriv_6(q(i-3:i+3,j,qu))
          
          i = hi(1)-2
          ! use 4th-order stencil
          ux(i,j) = dxinv(1)*first_deriv_4(q(i-2:i+2,j,qu))
          
          i = hi(1)-1
          ! use 3rd-order slightly left-biased stencil
          ux(i,j) = dxinv(1)*first_deriv_l3(q(i-2:i+1,j,qu))
          
          i = hi(1)
          ! use completely left-biased stencil
          ux(i,j) = dxinv(1)*first_deriv_lb(q(i-3:i,j,qu))
       end if
       
    end do

    do j=slo(2),shi(2)
       do i=lo(1),hi(1)
          vy(i,j) = dxinv(2)*first_deriv_8(q(i,j-4:j+4,qv))
       enddo
    end do

    ! lo-y boundary
    if (dlo(2) .eq. lo(2)) then
       j = lo(2)
       ! use completely right-biased stencil
       do i=lo(1),hi(1)
          vy(i,j) = dxinv(2)*first_deriv_rb(q(i,j:j+3,qv))
       enddo
       
       j = lo(2)+1
       ! use 3rd-order slightly right-biased stencil
       do i=lo(1),hi(1)
          vy(i,j) = dxinv(2)*first_deriv_r3(q(i,j-1:j+2,qv))
       enddo
       
       j = lo(2)+2
       ! use 4th-order stencil
       do i=lo(1),hi(1)
          vy(i,j) = dxinv(2)*first_deriv_4(q(i,j-2:j+2,qv))
       enddo
       
       j = lo(2)+3
       ! use 6th-order stencil
       do i=lo(1),hi(1)
          vy(i,j) = dxinv(2)*first_deriv_6(q(i,j-3:j+3,qv))
       enddo
    end if
    
    ! hi-y boundary
    if (dhi(2) .eq. hi(2)) then
       j = hi(2)-3
       ! use 6th-order stencil
       do i=lo(1),hi(1)
          vy(i,j) = dxinv(2)*first_deriv_6(q(i,j-3:j+3,qv))
       enddo
       
       j = hi(2)-2
       ! use 4th-order stencil
       do i=lo(1),hi(1)
          vy(i,j) = dxinv(2)*first_deriv_4(q(i,j-2:j+2,qv))
       enddo
       
       j = hi(2)-1
       ! use 3rd-order slightly left-biased stencil
       do i=lo(1),hi(1)
          vy(i,j) = dxinv(2)*first_deriv_l3(q(i,j-2:j+1,qv))
       enddo
       
       j = hi(2)
       ! use completely left-biased stencil
       do i=lo(1),hi(1)
          vy(i,j) = dxinv(2)*first_deriv_lb(q(i,j-3:j,qv))
       enddo
    end if

    divu = ux + vy

    deallocate(ux,vy)

  end subroutine make_divu_2d


  subroutine make_magvort_2d(lo, hi, mvor, vlo, vhi, Q, qlo, qhi, dx, dlo_g, dhi_g)
    use variables_module
    use derivative_stencil_module, only : stencil_ng, first_deriv_8, first_deriv_6, &
                first_deriv_4, first_deriv_l3, first_deriv_r3, first_deriv_rb, first_deriv_lb

    integer, intent(in) :: lo(2), hi(2), vlo(2), vhi(2), qlo(2), qhi(2), dlo_g(2), dhi_g(2)
    double precision, intent(inout) :: mvor(vlo(1):vhi(1),vlo(2):vhi(2))
    double precision, intent(in   ) ::    Q(qlo(1):qhi(1),qlo(2):qhi(2),nprim)
    double precision, intent(in) :: dx(2)

    integer :: i,j
    double precision, dimension(:,:), allocatable :: uy,vx
    double precision :: dxinv(2)
    integer :: slo(2), shi(2), dlo(2), dhi(2)

    ! Only the region bounded by [dlo,dhi] contains good data.
    ! [slo,shi] will be safe for 8th-order stencil
    do i=1,2
       dlo(i) = max(lo(i)-stencil_ng, dlo_g(i))
       dhi(i) = min(hi(i)+stencil_ng, dhi_g(i))
       slo(i) = dlo(i) + stencil_ng
       shi(i) = dhi(i) - stencil_ng
    end do

    do i=1,2
       dxinv(i) = 1.0d0 / dx(i)
    end do

    allocate(uy(lo(1):hi(1),lo(2):hi(2)))
    allocate(vx(lo(1):hi(1),lo(2):hi(2)))

    do j=lo(2),hi(2)
       
       do i=slo(1),shi(1)
          vx(i,j) = dxinv(1)*first_deriv_8(q(i-4:i+4,j,qv))
       enddo
       
       ! lo-x boundary
       if (dlo(1) .eq. lo(1)) then
          i = lo(1)
          ! use completely right-biased stencil
          vx(i,j) = dxinv(1)*first_deriv_rb(q(i:i+3,j,qv))
          
          i = lo(1)+1
          ! use 3rd-order slightly right-biased stencil
          vx(i,j) = dxinv(1)*first_deriv_r3(q(i-1:i+2,j,qv))
          
          i = lo(1)+2
          ! use 4th-order stencil
          vx(i,j) = dxinv(1)*first_deriv_4(q(i-2:i+2,j,qv))
          
          i = lo(1)+3
          ! use 6th-order stencil
          vx(i,j) = dxinv(1)*first_deriv_6(q(i-3:i+3,j,qv))
       end if

       ! hi-x boundary
       if (dhi(1) .eq. hi(1)) then
          i = hi(1)-3
          ! use 6th-order stencil
          vx(i,j) = dxinv(1)*first_deriv_6(q(i-3:i+3,j,qv))

          i = hi(1)-2
          ! use 4th-order stencil
          vx(i,j) = dxinv(1)*first_deriv_4(q(i-2:i+2,j,qv))
          
          i = hi(1)-1
          ! use 3rd-order slightly left-biased stencil
          vx(i,j) = dxinv(1)*first_deriv_l3(q(i-2:i+1,j,qv))
          
          i = hi(1)
          ! use completely left-biased stencil
          vx(i,j) = dxinv(1)*first_deriv_lb(q(i-3:i,j,qv))
       end if

    end do

    do j=slo(2),shi(2)
       do i=lo(1),hi(1)
          uy(i,j) = dxinv(2)*first_deriv_8(q(i,j-4:j+4,qu))
       enddo
    end do

    ! lo-y boundary
    if (dlo(2) .eq. lo(2)) then
       j = lo(2)
       ! use completely right-biased stencil
       do i=lo(1),hi(1)
          uy(i,j) = dxinv(2)*first_deriv_rb(q(i,j:j+3,qu))
       enddo
       
       j = lo(2)+1
       ! use 3rd-order slightly right-biased stencil
       do i=lo(1),hi(1)
          uy(i,j) = dxinv(2)*first_deriv_r3(q(i,j-1:j+2,qu))
       enddo

       j = lo(2)+2
       ! use 4th-order stencil
       do i=lo(1),hi(1)
          uy(i,j) = dxinv(2)*first_deriv_4(q(i,j-2:j+2,qu))
       enddo
       
       j = lo(2)+3
       ! use 6th-order stencil
       do i=lo(1),hi(1)
          uy(i,j) = dxinv(2)*first_deriv_6(q(i,j-3:j+3,qu))
       enddo
    end if

    ! hi-y boundary
    if (dhi(2) .eq. hi(2)) then
       j = hi(2)-3
       ! use 6th-order stencil
       do i=lo(1),hi(1)
          uy(i,j) = dxinv(2)*first_deriv_6(q(i,j-3:j+3,qu))
       enddo

       j = hi(2)-2
       ! use 4th-order stencil
       do i=lo(1),hi(1)
          uy(i,j) = dxinv(2)*first_deriv_4(q(i,j-2:j+2,qu))
       enddo
       
       j = hi(2)-1
       ! use 3rd-order slightly left-biased stencil
       do i=lo(1),hi(1)
          uy(i,j) = dxinv(2)*first_deriv_l3(q(i,j-2:j+1,qu))
       enddo
       
       j = hi(2)
       ! use completely left-biased stencil
       do i=lo(1),hi(1)
          uy(i,j) = dxinv(2)*first_deriv_lb(q(i,j-3:j,qu))
       enddo
    end if

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          mvor(i,j) = abs(vx(i,j)-uy(i,j))
       end do
    end do

    deallocate(uy,vx)

  end subroutine make_magvort_2d


  subroutine make_omegadot_2d(lo, hi, odot, vlo, vhi, Q, qlo, qhi)
    use variables_module
    use chemistry_module, only : molecular_weight

    integer, intent(in) :: lo(2), hi(2), vlo(2), vhi(2), qlo(2), qhi(2)
    double precision, intent(inout) :: odot(vlo(1):vhi(1),vlo(2):vhi(2),nspecies)
    double precision, intent(in   ) ::    Q(qlo(1):qhi(1),qlo(2):qhi(2),nprim)

    integer :: i,j,n,np,iwrk
    double precision :: Yt(lo(1):hi(1),nspecies), wdot(lo(1):hi(1),nspecies), rwrk

    np = hi(1) - lo(1) + 1

    do j=lo(2),hi(2)
       
       do n=1, nspecies
          do i=lo(1),hi(1)
             Yt(i,n) = q(i,j,qy1+n-1)
          end do
       end do
       
       call vckwyr(np, q(lo(1),j,qrho), q(lo(1),j,qtemp), Yt, iwrk, rwrk, wdot)
       
       do n=1, nspecies
          do i=lo(1),hi(1)
             odot(i,j,n) = wdot(i,n) * molecular_weight(n)
          end do
       end do

    enddo

  end subroutine make_omegadot_2d


  subroutine make_dYdt_2d(lo, hi, Ydot, vlo, vhi, Q, qlo, qhi)
    use variables_module
    use chemistry_module, only : molecular_weight

    integer, intent(in) :: lo(2), hi(2), vlo(2), vhi(2), qlo(2), qhi(2)
    double precision, intent(inout) :: Ydot(vlo(1):vhi(1),vlo(2):vhi(2),nspecies)
    double precision, intent(in   ) ::    Q(qlo(1):qhi(1),qlo(2):qhi(2),nprim)

    integer :: i,j,n,np,iwrk
    double precision :: Yt(lo(1):hi(1),nspecies), wdot(lo(1):hi(1),nspecies), rwrk

    np = hi(1) - lo(1) + 1

    do j=lo(2),hi(2)
       
       do n=1, nspecies
          do i=lo(1),hi(1)
             Yt(i,n) = q(i,j,qy1+n-1)
          end do
       end do
       
       call vckwyr(np, q(lo(1),j,qrho), q(lo(1),j,qtemp), Yt, iwrk, rwrk, wdot)
       
       do n=1, nspecies
          do i=lo(1),hi(1)
             Ydot(i,j,n) = wdot(i,n) * molecular_weight(n) / q(i,j,qrho)
          end do
       end do

    enddo

  end subroutine make_dYdt_2d


  subroutine make_heatRelease_2d(lo, hi, heat, vlo, vhi, Q, qlo, qhi)
    use variables_module
    use chemistry_module, only : molecular_weight

    integer, intent(in) :: lo(2), hi(2), vlo(2), vhi(2), qlo(2), qhi(2)
    double precision, intent(inout) :: heat(vlo(1):vhi(1),vlo(2):vhi(2))
    double precision, intent(in   ) ::    Q(qlo(1):qhi(1),qlo(2):qhi(2),nprim)

    integer :: i,j,n,np,iwrk
    double precision :: Yt(lo(1):hi(1),nspecies), wdot(lo(1):hi(1),nspecies), rwrk

    heat(lo(1):hi(1),lo(2):hi(2)) = 0.d0

  end subroutine make_heatRelease_2d


  subroutine make_fuelCsmp_2d(lo, hi, fuel, vlo, vhi, Q, qlo, qhi, ifuel)
    use variables_module
    use chemistry_module, only : molecular_weight

    integer, intent(in) :: lo(2), hi(2), vlo(2), vhi(2), qlo(2), qhi(2), ifuel
    double precision, intent(inout) :: fuel(vlo(1):vhi(1),vlo(2):vhi(2))
    double precision, intent(in   ) ::    Q(qlo(1):qhi(1),qlo(2):qhi(2),nprim)

    integer :: i,j,n,np,iwrk
    double precision :: Yt(lo(1):hi(1),nspecies), wdot(lo(1):hi(1),nspecies), rwrk

    fuel(lo(1):hi(1),lo(2):hi(2)) = 0.d0

  end subroutine make_fuelCsmp_2d

end module make_plotvar_2d_module
