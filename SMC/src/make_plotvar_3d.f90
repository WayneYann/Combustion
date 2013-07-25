module make_plotvar_3d_module

  implicit none

contains

  subroutine make_h_3d(lo, hi, h,  vlo, vhi, Q, qlo, qhi)
    use variables_module
    integer, intent(in) :: lo(3), hi(3), vlo(3), vhi(3), qlo(3), qhi(3)
    double precision, intent(inout) :: h(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
    double precision, intent(in   ) :: Q(qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3),nprim)

    integer :: i,j,k, n, qyn, qhn

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             h(i,j,k) = 0.d0
             do n=1, nspecies
                qhn = qh1+n-1
                qyn = qy1+n-1
                h(i,j,k) = h(i,j,k) + q(i,j,k,qyn)*q(i,j,k,qhn)
             end do
          enddo
       enddo
    enddo

  end subroutine make_h_3d


  subroutine make_divu_3d(lo, hi, divu, vlo, vhi, Q, qlo, qhi, dx, dlo_g, dhi_g)
    use variables_module
    use derivative_stencil_module, only : stencil_ng, first_deriv_8, first_deriv_6, &
                first_deriv_4, first_deriv_l3, first_deriv_r3, first_deriv_rb, first_deriv_lb

    integer, intent(in) :: lo(3), hi(3), vlo(3), vhi(3), qlo(3), qhi(3), dlo_g(3), dhi_g(3)
    double precision, intent(inout) :: divu(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
    double precision, intent(in   ) ::    Q(qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3),nprim)
    double precision, intent(in) :: dx(3)

    integer :: i,j,k
    double precision, allocatable :: ux(:,:,:), vy(:,:,:), wz(:,:,:)
    double precision :: dxinv(3)
    integer :: slo(3), shi(3), dlo(3), dhi(3)

    ! Only the region bounded by [dlo,dhi] contains good data.
    ! [slo,shi] will be safe for 8th-order stencil
    do i=1,3
       dlo(i) = max(lo(i)-stencil_ng, dlo_g(i))
       dhi(i) = min(hi(i)+stencil_ng, dhi_g(i))
       slo(i) = dlo(i) + stencil_ng
       shi(i) = dhi(i) - stencil_ng
    end do

    do i=1,3
       dxinv(i) = 1.0d0 / dx(i)
    end do

    allocate(ux(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
    allocate(vy(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
    allocate(wz(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)

          do i=slo(1),shi(1)
             ux(i,j,k) = dxinv(1)*first_deriv_8(q(i-4:i+4,j,k,qu))
          enddo
          
          ! lo-x boundary
          if (dlo(1) .eq. lo(1)) then
             i = lo(1)
             ! use completely right-biased stencil
             ux(i,j,k) = dxinv(1)*first_deriv_rb(q(i:i+3,j,k,qu))

             i = lo(1)+1
             ! use 3rd-order slightly right-biased stencil
             ux(i,j,k) = dxinv(1)*first_deriv_r3(q(i-1:i+2,j,k,qu))

             i = lo(1)+2
             ! use 4th-order stencil
             ux(i,j,k) = dxinv(1)*first_deriv_4(q(i-2:i+2,j,k,qu))

             i = lo(1)+3
             ! use 6th-order stencil
             ux(i,j,k) = dxinv(1)*first_deriv_6(q(i-3:i+3,j,k,qu))
          end if

          ! hi-x boundary
          if (dhi(1) .eq. hi(1)) then
             i = hi(1)-3
             ! use 6th-order stencil
             ux(i,j,k) = dxinv(1)*first_deriv_6(q(i-3:i+3,j,k,qu))

             i = hi(1)-2
             ! use 4th-order stencil
             ux(i,j,k) = dxinv(1)*first_deriv_4(q(i-2:i+2,j,k,qu))

             i = hi(1)-1
             ! use 3rd-order slightly left-biased stencil
             ux(i,j,k) = dxinv(1)*first_deriv_l3(q(i-2:i+1,j,k,qu))

             i = hi(1)
             ! use completely left-biased stencil
             ux(i,j,k) = dxinv(1)*first_deriv_lb(q(i-3:i,j,k,qu))
          end if

       end do
    end do

    do k=lo(3),hi(3)

       do j=slo(2),shi(2)
          do i=lo(1),hi(1)
             vy(i,j,k) = dxinv(2)*first_deriv_8(q(i,j-4:j+4,k,qv))
          enddo
       end do

       ! lo-y boundary
       if (dlo(2) .eq. lo(2)) then
          j = lo(2)
          ! use completely right-biased stencil
          do i=lo(1),hi(1)
             vy(i,j,k) = dxinv(2)*first_deriv_rb(q(i,j:j+3,k,qv))
          enddo

          j = lo(2)+1
          ! use 3rd-order slightly right-biased stencil
          do i=lo(1),hi(1)
             vy(i,j,k) = dxinv(2)*first_deriv_r3(q(i,j-1:j+2,k,qv))
          enddo

          j = lo(2)+2
          ! use 4th-order stencil
          do i=lo(1),hi(1)
             vy(i,j,k) = dxinv(2)*first_deriv_4(q(i,j-2:j+2,k,qv))
          enddo

          j = lo(2)+3
          ! use 6th-order stencil
          do i=lo(1),hi(1)
             vy(i,j,k) = dxinv(2)*first_deriv_6(q(i,j-3:j+3,k,qv))
          enddo
       end if

       ! hi-y boundary
       if (dhi(2) .eq. hi(2)) then
          j = hi(2)-3
          ! use 6th-order stencil
          do i=lo(1),hi(1)
             vy(i,j,k) = dxinv(2)*first_deriv_6(q(i,j-3:j+3,k,qv))
          enddo

          j = hi(2)-2
          ! use 4th-order stencil
          do i=lo(1),hi(1)
             vy(i,j,k) = dxinv(2)*first_deriv_4(q(i,j-2:j+2,k,qv))
          enddo

          j = hi(2)-1
          ! use 3rd-order slightly left-biased stencil
          do i=lo(1),hi(1)
             vy(i,j,k) = dxinv(2)*first_deriv_l3(q(i,j-2:j+1,k,qv))
          enddo

          j = hi(2)
          ! use completely left-biased stencil
          do i=lo(1),hi(1)
             vy(i,j,k) = dxinv(2)*first_deriv_lb(q(i,j-3:j,k,qv))
          enddo
       end if

    end do

    do k=slo(3),shi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             wz(i,j,k) = dxinv(3)*first_deriv_8(q(i,j,k-4:k+4,qw))
          enddo
       end do
    end do

    ! lo-z boundary
    if (dlo(3) .eq. lo(3)) then
       k = lo(3)
       ! use completely right-biased stencil
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             wz(i,j,k) = dxinv(3)*first_deriv_rb(q(i,j,k:k+3,qw))
          enddo
       enddo

       k = lo(3)+1
       ! use 3rd-order slightly right-biased stencil
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             wz(i,j,k) = dxinv(3)*first_deriv_r3(q(i,j,k-1:k+2,qw))
          enddo
       enddo

       k = lo(3)+2
       ! use 4th-order stencil
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             wz(i,j,k) = dxinv(3)*first_deriv_4(q(i,j,k-2:k+2,qw))
          enddo
       enddo

       k = lo(3)+3
       ! use 6th-order stencil
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             wz(i,j,k) = dxinv(3)*first_deriv_6(q(i,j,k-3:k+3,qw))
          enddo
       enddo
    end if

    ! hi-z boundary
    if (dhi(3) .eq. hi(3)) then
       k = hi(3)-3
       ! use 6th-order stencil
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             wz(i,j,k) = dxinv(3)*first_deriv_6(q(i,j,k-3:k+3,qw))
          enddo
       enddo

       k = hi(3)-2
       ! use 4th-order stencil
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             wz(i,j,k) = dxinv(3)*first_deriv_4(q(i,j,k-2:k+2,qw))
          enddo
       enddo

       k = hi(3)-1
       ! use 3rd-order slightly left-biased stencil
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             wz(i,j,k) = dxinv(3)*first_deriv_l3(q(i,j,k-2:k+1,qw))
          enddo
       enddo

       k = hi(3)
       ! use completely left-biased stencil
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             wz(i,j,k) = dxinv(3)*first_deriv_lb(q(i,j,k-3:k,qw))
          enddo
       enddo
    end if

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             divu(i,j,k) = ux(i,j,k) + vy(i,j,k) + wz(i,j,k)
          end do
       enddo
    enddo

    deallocate(ux,vy,wz)

  end subroutine make_divu_3d


  subroutine make_magvort_3d(lo, hi, mvor, vlo, vhi, Q, qlo, qhi, dx, dlo_g, dhi_g)
    use variables_module
    use derivative_stencil_module, only : stencil_ng, first_deriv_8, first_deriv_6, &
                first_deriv_4, first_deriv_l3, first_deriv_r3, first_deriv_rb, first_deriv_lb

    integer, intent(in) :: lo(3), hi(3), vlo(3), vhi(3), qlo(3), qhi(3), dlo_g(3), dhi_g(3)
    double precision, intent(inout) :: mvor(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3))
    double precision, intent(in   ) ::    Q(qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3),nprim)
    double precision, intent(in) :: dx(3)

    integer :: i,j,k
    double precision, dimension(:,:,:), allocatable :: uy,uz,vx,vz,wx,wy
    double precision :: dxinv(3)
    integer :: slo(3), shi(3), dlo(3), dhi(3)

    ! Only the region bounded by [dlo,dhi] contains good data.
    ! [slo,shi] will be safe for 8th-order stencil
    do i=1,3
       dlo(i) = max(lo(i)-stencil_ng, dlo_g(i))
       dhi(i) = min(hi(i)+stencil_ng, dhi_g(i))
       slo(i) = dlo(i) + stencil_ng
       shi(i) = dhi(i) - stencil_ng
    end do

    do i=1,3
       dxinv(i) = 1.0d0 / dx(i)
    end do

    allocate(uy(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
    allocate(uz(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
    allocate(vx(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
    allocate(vz(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
    allocate(wx(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
    allocate(wy(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)

          do i=slo(1),shi(1)
             vx(i,j,k) = dxinv(1)*first_deriv_8(q(i-4:i+4,j,k,qv))
             wx(i,j,k) = dxinv(1)*first_deriv_8(q(i-4:i+4,j,k,qw))
          enddo
          
          ! lo-x boundary
          if (dlo(1) .eq. lo(1)) then
             i = lo(1)
             ! use completely right-biased stencil
             vx(i,j,k) = dxinv(1)*first_deriv_rb(q(i:i+3,j,k,qv))
             wx(i,j,k) = dxinv(1)*first_deriv_rb(q(i:i+3,j,k,qw))

             i = lo(1)+1
             ! use 3rd-order slightly right-biased stencil
             vx(i,j,k) = dxinv(1)*first_deriv_r3(q(i-1:i+2,j,k,qv))
             wx(i,j,k) = dxinv(1)*first_deriv_r3(q(i-1:i+2,j,k,qw))

             i = lo(1)+2
             ! use 4th-order stencil
             vx(i,j,k) = dxinv(1)*first_deriv_4(q(i-2:i+2,j,k,qv))
             wx(i,j,k) = dxinv(1)*first_deriv_4(q(i-2:i+2,j,k,qw))

             i = lo(1)+3
             ! use 6th-order stencil
             vx(i,j,k) = dxinv(1)*first_deriv_6(q(i-3:i+3,j,k,qv))
             wx(i,j,k) = dxinv(1)*first_deriv_6(q(i-3:i+3,j,k,qw))
          end if

          ! hi-x boundary
          if (dhi(1) .eq. hi(1)) then
             i = hi(1)-3
             ! use 6th-order stencil
             vx(i,j,k) = dxinv(1)*first_deriv_6(q(i-3:i+3,j,k,qv))
             wx(i,j,k) = dxinv(1)*first_deriv_6(q(i-3:i+3,j,k,qw))

             i = hi(1)-2
             ! use 4th-order stencil
             vx(i,j,k) = dxinv(1)*first_deriv_4(q(i-2:i+2,j,k,qv))
             wx(i,j,k) = dxinv(1)*first_deriv_4(q(i-2:i+2,j,k,qw))

             i = hi(1)-1
             ! use 3rd-order slightly left-biased stencil
             vx(i,j,k) = dxinv(1)*first_deriv_l3(q(i-2:i+1,j,k,qv))
             wx(i,j,k) = dxinv(1)*first_deriv_l3(q(i-2:i+1,j,k,qw))

             i = hi(1)
             ! use completely left-biased stencil
             vx(i,j,k) = dxinv(1)*first_deriv_lb(q(i-3:i,j,k,qv))
             wx(i,j,k) = dxinv(1)*first_deriv_lb(q(i-3:i,j,k,qw))
          end if

       end do
    end do

    do k=lo(3),hi(3)

       do j=slo(2),shi(2)
          do i=lo(1),hi(1)
             uy(i,j,k) = dxinv(2)*first_deriv_8(q(i,j-4:j+4,k,qu))
             wy(i,j,k) = dxinv(2)*first_deriv_8(q(i,j-4:j+4,k,qw))
          enddo
       end do

       ! lo-y boundary
       if (dlo(2) .eq. lo(2)) then
          j = lo(2)
          ! use completely right-biased stencil
          do i=lo(1),hi(1)
             uy(i,j,k) = dxinv(2)*first_deriv_rb(q(i,j:j+3,k,qu))
             wy(i,j,k) = dxinv(2)*first_deriv_rb(q(i,j:j+3,k,qw))
          enddo

          j = lo(2)+1
          ! use 3rd-order slightly right-biased stencil
          do i=lo(1),hi(1)
             uy(i,j,k) = dxinv(2)*first_deriv_r3(q(i,j-1:j+2,k,qu))
             wy(i,j,k) = dxinv(2)*first_deriv_r3(q(i,j-1:j+2,k,qw))
          enddo

          j = lo(2)+2
          ! use 4th-order stencil
          do i=lo(1),hi(1)
             uy(i,j,k) = dxinv(2)*first_deriv_4(q(i,j-2:j+2,k,qu))
             wy(i,j,k) = dxinv(2)*first_deriv_4(q(i,j-2:j+2,k,qw))
          enddo

          j = lo(2)+3
          ! use 6th-order stencil
          do i=lo(1),hi(1)
             uy(i,j,k) = dxinv(2)*first_deriv_6(q(i,j-3:j+3,k,qu))
             wy(i,j,k) = dxinv(2)*first_deriv_6(q(i,j-3:j+3,k,qw))
          enddo
       end if

       ! hi-y boundary
       if (dhi(2) .eq. hi(2)) then
          j = hi(2)-3
          ! use 6th-order stencil
          do i=lo(1),hi(1)
             uy(i,j,k) = dxinv(2)*first_deriv_6(q(i,j-3:j+3,k,qu))
             wy(i,j,k) = dxinv(2)*first_deriv_6(q(i,j-3:j+3,k,qw))
          enddo

          j = hi(2)-2
          ! use 4th-order stencil
          do i=lo(1),hi(1)
             uy(i,j,k) = dxinv(2)*first_deriv_4(q(i,j-2:j+2,k,qu))
             wy(i,j,k) = dxinv(2)*first_deriv_4(q(i,j-2:j+2,k,qw))
          enddo

          j = hi(2)-1
          ! use 3rd-order slightly left-biased stencil
          do i=lo(1),hi(1)
             uy(i,j,k) = dxinv(2)*first_deriv_l3(q(i,j-2:j+1,k,qu))
             wy(i,j,k) = dxinv(2)*first_deriv_l3(q(i,j-2:j+1,k,qw))
          enddo

          j = hi(2)
          ! use completely left-biased stencil
          do i=lo(1),hi(1)
             uy(i,j,k) = dxinv(2)*first_deriv_lb(q(i,j-3:j,k,qu))
             wy(i,j,k) = dxinv(2)*first_deriv_lb(q(i,j-3:j,k,qw))
          enddo
       end if

    end do

    do k=slo(3),shi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             uz(i,j,k) = dxinv(3)*first_deriv_8(q(i,j,k-4:k+4,qu))
             vz(i,j,k) = dxinv(3)*first_deriv_8(q(i,j,k-4:k+4,qv))
          enddo
       end do
    end do

    ! lo-z boundary
    if (dlo(3) .eq. lo(3)) then
       k = lo(3)
       ! use completely right-biased stencil
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             uz(i,j,k) = dxinv(3)*first_deriv_rb(q(i,j,k:k+3,qu))
             vz(i,j,k) = dxinv(3)*first_deriv_rb(q(i,j,k:k+3,qv))
          enddo
       enddo

       k = lo(3)+1
       ! use 3rd-order slightly right-biased stencil
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             uz(i,j,k) = dxinv(3)*first_deriv_r3(q(i,j,k-1:k+2,qu))
             vz(i,j,k) = dxinv(3)*first_deriv_r3(q(i,j,k-1:k+2,qv))
          enddo
       enddo

       k = lo(3)+2
       ! use 4th-order stencil
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             uz(i,j,k) = dxinv(3)*first_deriv_4(q(i,j,k-2:k+2,qu))
             vz(i,j,k) = dxinv(3)*first_deriv_4(q(i,j,k-2:k+2,qv))
          enddo
       enddo

       k = lo(3)+3
       ! use 6th-order stencil
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             uz(i,j,k) = dxinv(3)*first_deriv_6(q(i,j,k-3:k+3,qu))
             vz(i,j,k) = dxinv(3)*first_deriv_6(q(i,j,k-3:k+3,qv))
          enddo
       enddo
    end if

    ! hi-z boundary
    if (dhi(3) .eq. hi(3)) then
       k = hi(3)-3
       ! use 6th-order stencil
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             uz(i,j,k) = dxinv(3)*first_deriv_6(q(i,j,k-3:k+3,qu))
             vz(i,j,k) = dxinv(3)*first_deriv_6(q(i,j,k-3:k+3,qv))
          enddo
       enddo

       k = hi(3)-2
       ! use 4th-order stencil
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             uz(i,j,k) = dxinv(3)*first_deriv_4(q(i,j,k-2:k+2,qu))
             vz(i,j,k) = dxinv(3)*first_deriv_4(q(i,j,k-2:k+2,qv))
          enddo
       enddo

       k = hi(3)-1
       ! use 3rd-order slightly left-biased stencil
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             uz(i,j,k) = dxinv(3)*first_deriv_l3(q(i,j,k-2:k+1,qu))
             vz(i,j,k) = dxinv(3)*first_deriv_l3(q(i,j,k-2:k+1,qv))
          enddo
       enddo

       k = hi(3)
       ! use completely left-biased stencil
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             uz(i,j,k) = dxinv(3)*first_deriv_lb(q(i,j,k-3:k,qu))
             vz(i,j,k) = dxinv(3)*first_deriv_lb(q(i,j,k-3:k,qv))
          enddo
       enddo
    end if

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             mvor(i,j,k) = sqrt((wy(i,j,k)-vz(i,j,k))**2  &
                  +             (uz(i,j,k)-wx(i,j,k))**2  &
                  +             (vx(i,j,k)-uy(i,j,k))**2  )
          end do
       end do
    end do

    deallocate(uy,uz,vx,vz,wx,wy)

  end subroutine make_magvort_3d


  subroutine make_burn_3d(lo, hi, burn, vlo, vhi, Q, qlo, qhi)
    use plotvar_index_module
    use variables_module
    use chemistry_module, only : molecular_weight

    integer, intent(in) :: lo(3), hi(3), vlo(3), vhi(3), qlo(3), qhi(3)
    double precision, intent(inout) :: burn(vlo(1):vhi(1),vlo(2):vhi(2),vlo(3):vhi(3),nburn)
    double precision, intent(in   ) ::    Q(qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3),nprim)

    integer :: i,j,k,n,np,iwrk
    double precision :: Yt(lo(1):hi(1),nspecies), wdot(lo(1):hi(1),nspecies), rwrk

    np = hi(1) - lo(1) + 1

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)

          do n=1, nspecies
             do i=lo(1),hi(1)
                Yt(i,n) = q(i,j,k,qy1+n-1)
             end do
          end do
          
          call vckwyr(np, q(lo(1),j,k,qrho), q(lo(1),j,k,qtemp), Yt, iwrk, rwrk, wdot)

          if (ib_omegadot > 0) then
             do n=1, nspecies
                do i=lo(1),hi(1)
                   burn(i,j,k,ib_omegadot+n-1) = wdot(i,n) * molecular_weight(n)
                end do
             end do
          end if

          if (ib_dYdt > 0) then
             do n=1, nspecies
                do i=lo(1),hi(1)
                   burn(i,j,k,ib_dYdt+n-1) = wdot(i,n) * molecular_weight(n) / q(i,j,k,qrho)
                end do
             end do
          end if
  
          if (ib_heatRelease > 0) then
             do i=lo(1),hi(1)
                burn(i,j,k,ib_heatRelease) = 0.d0
             end do

             do n=1,nspecies
                burn(i,j,k,ib_heatRelease) = burn(i,j,k,ib_heatRelease) &
                     - q(i,j,k,qh1+n-1) * wdot(i,n) * molecular_weight(n)
             end do
          end if
          
          if (ib_fuelConsumption > 0) then
             do i=lo(1),hi(1)
                burn(i,j,k,ib_fuelConsumption) = wdot(i,ifuel) * molecular_weight(ifuel)
             end do
          end if
          
       enddo
    enddo

  end subroutine make_burn_3d

end module make_plotvar_3d_module
