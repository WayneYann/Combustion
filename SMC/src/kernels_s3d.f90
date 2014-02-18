module kernels_s3d_module

  use bc_module
  use chemistry_module, only : nspecies
  use derivative_stencil_module, only : stencil_ng, first_deriv_8, first_deriv_6, &
       first_deriv_4, first_deriv_l3, first_deriv_r3, first_deriv_rb, first_deriv_lb
  use variables_module

  implicit none

  private

  public :: s3d_diffterm_1_3d, s3d_diffterm_2_3d

contains

  subroutine s3d_diffterm_1_3d (lo,hi,dx,q,qlo,qhi,rhs,rlo,rhi, &
       qx,qxlo,qxhi,qy,qylo,qyhi,qz,qzlo,qzhi, &
       mu,xi,dlo_g,dhi_g,bclo,bchi)

    integer,         intent(in):: dlo_g(3),dhi_g(3),bclo(3),bchi(3)
    integer,         intent(in):: lo(3),hi(3),qlo(3),qhi(3),rlo(3),rhi(3)
    integer,         intent(in):: qxlo(3),qxhi(3),qylo(3),qyhi(3),qzlo(3),qzhi(3)
    double precision,intent(in):: dx(3)
    double precision,intent(in)   ::q  ( qlo(1): qhi(1), qlo(2): qhi(2), qlo(3): qhi(3),nprim)
    double precision,intent(in)   ::mu ( qlo(1): qhi(1), qlo(2): qhi(2), qlo(3): qhi(3))
    double precision,intent(in)   ::xi ( qlo(1): qhi(1), qlo(2): qhi(2), qlo(3): qhi(3))
    double precision,intent(inout)::rhs( rlo(1): rhi(1), rlo(2): rhi(2), rlo(3): rhi(3),ncons)
    double precision,intent(inout)::qx (qxlo(1):qxhi(1),qxlo(2):qxhi(2),qxlo(3):qxhi(3),ndq)
    double precision,intent(inout)::qy (qylo(1):qyhi(1),qylo(2):qyhi(2),qylo(3):qyhi(3),ndq)
    double precision,intent(inout)::qz (qzlo(1):qzhi(1),qzlo(2):qzhi(2),qzlo(3):qzhi(3),ndq)

    integer :: i,j,k,n, qxn, qdxn
    integer :: slo(3), shi(3), dlo(3), dhi(3)
    double precision :: dxinv(3)
    double precision, allocatable, dimension(:,:,:) :: vsm
    double precision, allocatable :: tmpx(:), tmpy(:,:),tmpz(:,:,:)
    double precision, dimension(lo(1):hi(1)) :: tauxx,tauyy,tauzz,divu

    ! used to turn off some terms
    double precision :: finlo(3), finhi(3)
    double precision :: foulo(3), fouhi(3)

    logical :: physbclo(3), physbchi(3)

    rhs(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:) = 0.d0

    ! Only the region bounded by [dlo_g,dhi_g] contains good data.
    ! [slo,shi] will be safe for 8th-order stencil
    do i=1,3
       dlo(i) = max(lo(i)-stencil_ng, dlo_g(i))
       dhi(i) = min(hi(i)+stencil_ng, dhi_g(i))
       slo(i) = dlo(i) + stencil_ng
       shi(i) = dhi(i) - stencil_ng

       if (dlo(i) .eq. lo(i)) then
          physbclo(i) = .true.
       else
          physbclo(i) = .false.
       end if

       if (dhi(i) .eq. hi(i)) then
          physbchi(i) = .true.
       else
          physbchi(i) = .false.
       end if
    end do

    finlo = 1.d0 
    finhi = 1.d0
    foulo = 1.d0 
    fouhi = 1.d0

    do i=1,3
       if (physbclo(i)) then 
          if (bclo(i) .eq. INLET) then
             finlo(i) = 0.d0
          else if (bclo(i) .eq. OUTLET) then
             foulo(i) = 0.d0
          end if
       end if

       if (physbchi(i)) then 
          if (bchi(i) .eq. INLET) then
             finhi(i) = 0.d0
          else if (bchi(i) .eq. OUTLET) then
             fouhi(i) = 0.d0
          end if
       end if
    end do

    do i = 1,3
       dxinv(i) = 1.0d0 / dx(i)
    end do

    allocate(vsm (dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3)))

    allocate(tmpx(dlo(1):dhi(1)))
    allocate(tmpy( lo(1): hi(1),dlo(2):dhi(2)))
    allocate(tmpz( lo(1): hi(1), lo(2): hi(2),dlo(3):dhi(3)))

    do k=dlo(3),dhi(3)
       do j=dlo(2),dhi(2)
          do i=dlo(1),dhi(1)
             vsm(i,j,k) = xi(i,j,k) -  TwoThirds*mu(i,j,k)
          enddo
       enddo
    enddo

    do k=dlo(3),dhi(3)
       do j=dlo(2),dhi(2)
          do i=slo(1),shi(1)
             qx(i,j,k,idu) = dxinv(1) * first_deriv_8( q(i-4:i+4,j,k,qu) )
             qx(i,j,k,idv) = dxinv(1) * first_deriv_8( q(i-4:i+4,j,k,qv) )
             qx(i,j,k,idw) = dxinv(1) * first_deriv_8( q(i-4:i+4,j,k,qw) )
          enddo
       enddo
    enddo

    do k=dlo(3),dhi(3)
       do j=slo(2),shi(2)   
          do i=dlo(1),dhi(1)
             qy(i,j,k,idu) = dxinv(2) * first_deriv_8( q(i,j-4:j+4,k,qu) )
             qy(i,j,k,idv) = dxinv(2) * first_deriv_8( q(i,j-4:j+4,k,qv) )
             qy(i,j,k,idw) = dxinv(2) * first_deriv_8( q(i,j-4:j+4,k,qw) )
          enddo
       enddo
    enddo

    do k=slo(3),shi(3)
       do j=dlo(2),dhi(2)
          do i=dlo(1),dhi(1)
             qz(i,j,k,idu) = dxinv(3) * first_deriv_8( q(i,j,k-4:k+4,qu) )
             qz(i,j,k,idv) = dxinv(3) * first_deriv_8( q(i,j,k-4:k+4,qv) )
             qz(i,j,k,idw) = dxinv(3) * first_deriv_8( q(i,j,k-4:k+4,qw) )
          enddo
       enddo
    enddo

    !
    ! lo-x boundary
    !
    if (physbclo(1)) then
       do k=dlo(3),dhi(3)
          do j=dlo(2),dhi(2)
             i = lo(1)
             ! use completely right-biased stencil
             qx(i,j,k,idu) = dxinv(1)*first_deriv_rb(q(i:i+3,j,k,qu))
             qx(i,j,k,idv) = dxinv(1)*first_deriv_rb(q(i:i+3,j,k,qv))
             qx(i,j,k,idw) = dxinv(1)*first_deriv_rb(q(i:i+3,j,k,qw))

             i = lo(1)+1
             ! use 3rd-order slightly right-biased stencil
             qx(i,j,k,idu) = dxinv(1)*first_deriv_r3(q(i-1:i+2,j,k,qu))
             qx(i,j,k,idv) = dxinv(1)*first_deriv_r3(q(i-1:i+2,j,k,qv))
             qx(i,j,k,idw) = dxinv(1)*first_deriv_r3(q(i-1:i+2,j,k,qw))

             i = lo(1)+2
             ! use 4th-order stencil
             qx(i,j,k,idu) = dxinv(1)*first_deriv_4(q(i-2:i+2,j,k,qu))
             qx(i,j,k,idv) = dxinv(1)*first_deriv_4(q(i-2:i+2,j,k,qv))
             qx(i,j,k,idw) = dxinv(1)*first_deriv_4(q(i-2:i+2,j,k,qw))

             i = lo(1)+3
             ! use 6th-order stencil
             qx(i,j,k,idu) = dxinv(1)*first_deriv_6(q(i-3:i+3,j,k,qu))
             qx(i,j,k,idv) = dxinv(1)*first_deriv_6(q(i-3:i+3,j,k,qv))
             qx(i,j,k,idw) = dxinv(1)*first_deriv_6(q(i-3:i+3,j,k,qw))
          end do
       end do
    end if

    !
    ! hi-x boundary
    !
    if (physbchi(1)) then
       do k=dlo(3),dhi(3)
          do j=dlo(2),dhi(2)
             i = hi(1)-3
             ! use 6th-order stencil
             qx(i,j,k,idu) = dxinv(1)*first_deriv_6(q(i-3:i+3,j,k,qu))
             qx(i,j,k,idv) = dxinv(1)*first_deriv_6(q(i-3:i+3,j,k,qv))
             qx(i,j,k,idw) = dxinv(1)*first_deriv_6(q(i-3:i+3,j,k,qw))

             i = hi(1)-2
             ! use 4th-order stencil
             qx(i,j,k,idu) = dxinv(1)*first_deriv_4(q(i-2:i+2,j,k,qu))
             qx(i,j,k,idv) = dxinv(1)*first_deriv_4(q(i-2:i+2,j,k,qv))
             qx(i,j,k,idw) = dxinv(1)*first_deriv_4(q(i-2:i+2,j,k,qw))

             i = hi(1)-1
             ! use 3rd-order slightly left-biased stencil
             qx(i,j,k,idu) = dxinv(1)*first_deriv_l3(q(i-2:i+1,j,k,qu))
             qx(i,j,k,idv) = dxinv(1)*first_deriv_l3(q(i-2:i+1,j,k,qv))
             qx(i,j,k,idw) = dxinv(1)*first_deriv_l3(q(i-2:i+1,j,k,qw))

             i = hi(1)
             ! use completely left-biased stencil
             qx(i,j,k,idu) = dxinv(1)*first_deriv_lb(q(i-3:i,j,k,qu))
             qx(i,j,k,idv) = dxinv(1)*first_deriv_lb(q(i-3:i,j,k,qv))
             qx(i,j,k,idw) = dxinv(1)*first_deriv_lb(q(i-3:i,j,k,qw))
          end do
       end do
    end if

    !
    ! lo-y boundary
    !
    if (physbclo(2)) then
       do k=dlo(3),dhi(3)
          j = lo(2)
          ! use completely right-biased stencil
          do i=dlo(1),dhi(1)
             qy(i,j,k,idu) = dxinv(2)*first_deriv_rb(q(i,j:j+3,k,qu))
             qy(i,j,k,idv) = dxinv(2)*first_deriv_rb(q(i,j:j+3,k,qv))
             qy(i,j,k,idw) = dxinv(2)*first_deriv_rb(q(i,j:j+3,k,qw))
          enddo

          j = lo(2)+1
          ! use 3rd-order slightly right-biased stencil
          do i=dlo(1),dhi(1)
             qy(i,j,k,idu) = dxinv(2)*first_deriv_r3(q(i,j-1:j+2,k,qu))
             qy(i,j,k,idv) = dxinv(2)*first_deriv_r3(q(i,j-1:j+2,k,qv))
             qy(i,j,k,idw) = dxinv(2)*first_deriv_r3(q(i,j-1:j+2,k,qw))
          enddo

          j = lo(2)+2
          ! use 4th-order stencil
          do i=dlo(1),dhi(1)
             qy(i,j,k,idu) = dxinv(2)*first_deriv_4(q(i,j-2:j+2,k,qu))
             qy(i,j,k,idv) = dxinv(2)*first_deriv_4(q(i,j-2:j+2,k,qv))
             qy(i,j,k,idw) = dxinv(2)*first_deriv_4(q(i,j-2:j+2,k,qw))
          enddo

          j = lo(2)+3
          ! use 6th-order stencil
          do i=dlo(1),dhi(1)
             qy(i,j,k,idu) = dxinv(2)*first_deriv_6(q(i,j-3:j+3,k,qu))
             qy(i,j,k,idv) = dxinv(2)*first_deriv_6(q(i,j-3:j+3,k,qv))
             qy(i,j,k,idw) = dxinv(2)*first_deriv_6(q(i,j-3:j+3,k,qw))
          enddo
       end do
    end if

    !
    ! hi-y boundary
    !
    if (physbchi(2)) then
       do k=dlo(3),dhi(3)
          j = hi(2)-3
          ! use 6th-order stencil
          do i=dlo(1),dhi(1)
             qy(i,j,k,idu) = dxinv(2)*first_deriv_6(q(i,j-3:j+3,k,qu))
             qy(i,j,k,idv) = dxinv(2)*first_deriv_6(q(i,j-3:j+3,k,qv))
             qy(i,j,k,idw) = dxinv(2)*first_deriv_6(q(i,j-3:j+3,k,qw))
          enddo

          j = hi(2)-2
          ! use 4th-order stencil
          do i=dlo(1),dhi(1)
             qy(i,j,k,idu) = dxinv(2)*first_deriv_4(q(i,j-2:j+2,k,qu))
             qy(i,j,k,idv) = dxinv(2)*first_deriv_4(q(i,j-2:j+2,k,qv))
             qy(i,j,k,idw) = dxinv(2)*first_deriv_4(q(i,j-2:j+2,k,qw))
          enddo

          j = hi(2)-1
          ! use 3rd-order slightly left-biased stencil
          do i=dlo(1),dhi(1)
             qy(i,j,k,idu) = dxinv(2)*first_deriv_l3(q(i,j-2:j+1,k,qu))
             qy(i,j,k,idv) = dxinv(2)*first_deriv_l3(q(i,j-2:j+1,k,qv))
             qy(i,j,k,idw) = dxinv(2)*first_deriv_l3(q(i,j-2:j+1,k,qw))
          enddo

          j = hi(2)
          ! use completely left-biased stencil
          do i=dlo(1),dhi(1)
             qy(i,j,k,idu) = dxinv(2)*first_deriv_lb(q(i,j-3:j,k,qu))
             qy(i,j,k,idv) = dxinv(2)*first_deriv_lb(q(i,j-3:j,k,qv))
             qy(i,j,k,idw) = dxinv(2)*first_deriv_lb(q(i,j-3:j,k,qw))
          enddo
       end do
    end if

    !
    ! lo-z boundary
    !
    if (physbclo(3)) then
       k = lo(3)
       ! use completely right-biased stencil
       do j=dlo(2),dhi(2)
          do i=dlo(1),dhi(1)
             qz(i,j,k,idu) = dxinv(3)*first_deriv_rb(q(i,j,k:k+3,qu))
             qz(i,j,k,idv) = dxinv(3)*first_deriv_rb(q(i,j,k:k+3,qv))
             qz(i,j,k,idw) = dxinv(3)*first_deriv_rb(q(i,j,k:k+3,qw))
          enddo
       enddo

       k = lo(3)+1
       ! use 3rd-order slightly right-biased stencil
       do j=dlo(2),dhi(2)
          do i=dlo(1),dhi(1)
             qz(i,j,k,idu) = dxinv(3)*first_deriv_r3(q(i,j,k-1:k+2,qu))
             qz(i,j,k,idv) = dxinv(3)*first_deriv_r3(q(i,j,k-1:k+2,qv))
             qz(i,j,k,idw) = dxinv(3)*first_deriv_r3(q(i,j,k-1:k+2,qw))
          enddo
       enddo

       k = lo(3)+2
       ! use 4th-order stencil
       do j=dlo(2),dhi(2)
          do i=dlo(1),dhi(1)
             qz(i,j,k,idu) = dxinv(3)*first_deriv_4(q(i,j,k-2:k+2,qu))
             qz(i,j,k,idv) = dxinv(3)*first_deriv_4(q(i,j,k-2:k+2,qv))
             qz(i,j,k,idw) = dxinv(3)*first_deriv_4(q(i,j,k-2:k+2,qw))
          enddo
       enddo

       k = lo(3)+3
       ! use 6th-order stencil
       do j=dlo(2),dhi(2)
          do i=dlo(1),dhi(1)
             qz(i,j,k,idu) = dxinv(3)*first_deriv_6(q(i,j,k-3:k+3,qu))
             qz(i,j,k,idv) = dxinv(3)*first_deriv_6(q(i,j,k-3:k+3,qv))
             qz(i,j,k,idw) = dxinv(3)*first_deriv_6(q(i,j,k-3:k+3,qw))
          enddo
       enddo
    end if

    !
    ! hi-z boundary
    !
    if (physbchi(3)) then
       k = hi(3)-3
       ! use 6th-order stencil
       do j=dlo(2),dhi(2)
          do i=dlo(1),dhi(1)
             qz(i,j,k,idu) = dxinv(3)*first_deriv_6(q(i,j,k-3:k+3,qu))
             qz(i,j,k,idv) = dxinv(3)*first_deriv_6(q(i,j,k-3:k+3,qv))
             qz(i,j,k,idw) = dxinv(3)*first_deriv_6(q(i,j,k-3:k+3,qw))
          enddo
       enddo

       k = hi(3)-2
       ! use 4th-order stencil
       do j=dlo(2),dhi(2)
          do i=dlo(1),dhi(1)
             qz(i,j,k,idu) = dxinv(3)*first_deriv_4(q(i,j,k-2:k+2,qu))
             qz(i,j,k,idv) = dxinv(3)*first_deriv_4(q(i,j,k-2:k+2,qv))
             qz(i,j,k,idw) = dxinv(3)*first_deriv_4(q(i,j,k-2:k+2,qw))
          enddo
       enddo

       k = hi(3)-1
       ! use 3rd-order slightly left-biased stencil
       do j=dlo(2),dhi(2)
          do i=dlo(1),dhi(1)
             qz(i,j,k,idu) = dxinv(3)*first_deriv_l3(q(i,j,k-2:k+1,qu))
             qz(i,j,k,idv) = dxinv(3)*first_deriv_l3(q(i,j,k-2:k+1,qv))
             qz(i,j,k,idw) = dxinv(3)*first_deriv_l3(q(i,j,k-2:k+1,qw))
          enddo
       enddo

       k = hi(3)
       ! use completely left-biased stencil
       do j=dlo(2),dhi(2)
          do i=dlo(1),dhi(1)
             qz(i,j,k,idu) = dxinv(3)*first_deriv_lb(q(i,j,k-3:k,qu))
             qz(i,j,k,idv) = dxinv(3)*first_deriv_lb(q(i,j,k-3:k,qv))
             qz(i,j,k,idw) = dxinv(3)*first_deriv_lb(q(i,j,k-3:k,qw))
          enddo
       enddo
    end if


    !----- mx -----

    !----- mx : d()/dx -----
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)

          do i=dlo(1),dhi(1)
             tmpx(i) = vsm(i,j,k)*(qy(i,j,k,idv)+qz(i,j,k,idw))
          end do

          !
          ! lo-x boundary
          !
          if (physbclo(1)) then
             i = lo(1)
             ! use completely right-biased stencil
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + finlo(1)*dxinv(1)*first_deriv_rb(tmpx(i:i+3))
             
             i = lo(1)+1
             ! use 3rd-order slightly right-biased stencil
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + dxinv(1)*first_deriv_r3(tmpx(i-1:i+2))

             i = lo(1)+2
             ! use 4th-order stencil
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + dxinv(1)*first_deriv_4(tmpx(i-2:i+2))

             i = lo(1)+3
             ! use 6th-order stencil
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + dxinv(1)*first_deriv_6(tmpx(i-3:i+3))
          end if

          do i=slo(1),shi(1)
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + dxinv(1)*first_deriv_8(tmpx(i-4:i+4)) 
          end do

          !
          ! hi-x boundary
          !
          if (physbchi(1)) then
             i = hi(1)-3
             ! use 6th-order stencil
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + dxinv(1)*first_deriv_6(tmpx(i-3:i+3))

             i = hi(1)-2
             ! use 4th-order stencil
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + dxinv(1)*first_deriv_4(tmpx(i-2:i+2))

             i = hi(1)-1
             ! use 3rd-order slightly left-biased stencil
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + dxinv(1)*first_deriv_l3(tmpx(i-2:i+1))

             i = hi(1)
             ! use completely left-biased stencil
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + finhi(1)*dxinv(1)*first_deriv_lb(tmpx(i-3:i))
          end if
       end do
    end do

    !----- mx : d()/dy -----
    do k=lo(3),hi(3)

       do j=dlo(2),dhi(2)
          do i=lo(1),hi(1)
             tmpy(i,j) = mu(i,j,k)*qx(i,j,k,idv)
          end do
       end do

       !
       ! lo-y boundary
       !
       if (physbclo(2)) then
          j = lo(2)
          ! use completely right-biased stencil
          do i=lo(1),hi(1)
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + foulo(2)*dxinv(2)*first_deriv_rb(tmpy(i,j:j+3))
          end do

          j = lo(2)+1
          ! use 3rd-order slightly right-biased stencil
          do i=lo(1),hi(1)
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + dxinv(2)*first_deriv_r3(tmpy(i,j-1:j+2))
          end do

          j = lo(2)+2
          ! use 4th-order stencil
          do i=lo(1),hi(1)
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + dxinv(2)*first_deriv_4(tmpy(i,j-2:j+2))
          end do

          j = lo(2)+3
          ! use 6th-order stencil
          do i=lo(1),hi(1)
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + dxinv(2)*first_deriv_6(tmpy(i,j-3:j+3))
          end do
       end if

       !
       ! interior
       !
       do j=slo(2),shi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + dxinv(2)*first_deriv_8(tmpy(i,j-4:j+4)) 
          end do
       end do

       !
       ! hi-y boundary
       !
       if (physbchi(2)) then
          j = hi(2)-3
          ! use 6th-order stencil
          do i=lo(1),hi(1)
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + dxinv(2)*first_deriv_6(tmpy(i,j-3:j+3))
          end do

          j = hi(2)-2
          ! use 4th-order stencil
          do i=lo(1),hi(1)
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + dxinv(2)*first_deriv_4(tmpy(i,j-2:j+2))
          end do

          j = hi(2)-1
          ! use 3rd-order slightly left-biased stencil
          do i=lo(1),hi(1)
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + dxinv(2)*first_deriv_l3(tmpy(i,j-2:j+1))
          end do

          j = hi(2)
          ! use completely left-biased stencil
          do i=lo(1),hi(1)
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + fouhi(2)*dxinv(2)*first_deriv_lb(tmpy(i,j-3:j))
          end do
       end if
    end do

    ! ----- mx : d()/dz -----

    do k=dlo(3),dhi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             tmpz(i,j,k) = mu(i,j,k)*qx(i,j,k,idw)
          end do
       end do
    end do

    ! lo-z boundary
    if (physbclo(3)) then
       k = lo(3)
       ! use completely right-biased stencil
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + foulo(3)*dxinv(3)*first_deriv_rb(tmpz(i,j,k:k+3))
          end do
       end do

       k = lo(3)+1
       ! use 3rd-order slightly right-biased stencil
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + dxinv(3)*first_deriv_r3(tmpz(i,j,k-1:k+2))
          end do
       end do

       k = lo(3)+2
       ! use 4th-order stencil
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + dxinv(3)*first_deriv_4(tmpz(i,j,k-2:k+2))
          end do
       end do

       k = lo(3)+3
       ! use 6th-order stencil
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + dxinv(3)*first_deriv_6(tmpz(i,j,k-3:k+3))
          end do
       end do
    end if

    !
    ! interior
    !
    do k=slo(3),shi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + dxinv(3)*first_deriv_8(tmpz(i,j,k-4:k+4)) 
          end do
       end do
    end do

    !
    ! hi-z boundary
    !
    if (physbchi(3)) then
       k = hi(3)-3
       ! use 6th-order stencil
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + dxinv(3)*first_deriv_6(tmpz(i,j,k-3:k+3))
          end do
       end do

       k = hi(3)-2
       ! use 4th-order stencil
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + dxinv(3)*first_deriv_4(tmpz(i,j,k-2:k+2))
          end do
       end do

       k = hi(3)-1
       ! use 3rd-order slightly left-biased stencil
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + dxinv(3)*first_deriv_l3(tmpz(i,j,k-2:k+1))
          end do
       end do

       k = hi(3)
       ! use completely left-biased stencil
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + fouhi(3)*dxinv(3)*first_deriv_lb(tmpz(i,j,k-3:k))
          end do
       end do
    end if

    !----- my -----

    !----- my : d()/dx -----
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)

          do i=dlo(1),dhi(1)
             tmpx(i) = mu(i,j,k)*qy(i,j,k,idu)
          end do

          !
          ! lo-x boundary
          !
          if (physbclo(1)) then
             i = lo(1)
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + foulo(1)*dxinv(1)*first_deriv_rb(tmpx(i:i+3))

             i = lo(1)+1
             ! use 3rd-order slightly right-biased stencil
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + dxinv(1)*first_deriv_r3(tmpx(i-1:i+2))

             i = lo(1)+2
             ! use 4th-order stencil
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + dxinv(1)*first_deriv_4(tmpx(i-2:i+2))

             i = lo(1)+3
             ! use 6th-order stencil
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + dxinv(1)*first_deriv_6(tmpx(i-3:i+3))
          end if

          !
          ! interior
          !
          do i=slo(1),shi(1)
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + dxinv(1)*first_deriv_8(tmpx(i-4:i+4)) 
          end do

          !
          ! hi-x boundary
          !
          if (physbchi(1)) then
             i = hi(1)-3
             ! use 6th-order stencil
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + dxinv(1)*first_deriv_6(tmpx(i-3:i+3))

             i = hi(1)-2
             ! use 4th-order stencil
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + dxinv(1)*first_deriv_4(tmpx(i-2:i+2))

             i = hi(1)-1
             ! use 3rd-order slightly left-biased stencil
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + dxinv(1)*first_deriv_l3(tmpx(i-2:i+1))

             i = hi(1)
             ! use completely left-biased stencil
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + fouhi(1)*dxinv(1)*first_deriv_lb(tmpx(i-3:i))
          end if
       end do
    end do

    !----- my : d()/dy -----
    do k=lo(3),hi(3)

       do j=dlo(2),dhi(2)
          do i=lo(1),hi(1)
             tmpy(i,j) = vsm(i,j,k)*(qx(i,j,k,idu)+qz(i,j,k,idw))
          end do
       end do

       !
       ! lo-y boundary
       !
       if (physbclo(2)) then
          j = lo(2)
          ! use completely right-biased stencil
          do i=lo(1),hi(1)
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + finlo(2)*dxinv(2)*first_deriv_rb(tmpy(i,j:j+3))
          end do

          j = lo(2)+1
          ! use 3rd-order slightly right-biased stencil
          do i=lo(1),hi(1)
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + dxinv(2)*first_deriv_r3(tmpy(i,j-1:j+2))
          end do

          j = lo(2)+2
          ! use 4th-order stencil
          do i=lo(1),hi(1)
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + dxinv(2)*first_deriv_4(tmpy(i,j-2:j+2))
          end do

          j = lo(2)+3
          ! use 6th-order stencil
          do i=lo(1),hi(1)
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + dxinv(2)*first_deriv_6(tmpy(i,j-3:j+3))
          end do
       end if

       !
       ! interior
       !
       do j=slo(2),shi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + dxinv(2)*first_deriv_8(tmpy(i,j-4:j+4)) 
          end do
       end do

       !
       ! hi-y boundary
       !
       if (physbchi(2)) then
          j = hi(2)-3
          ! use 6th-order stencil
          do i=lo(1),hi(1)
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + dxinv(2)*first_deriv_6(tmpy(i,j-3:j+3))
          end do

          j = hi(2)-2
          ! use 4th-order stencil
          do i=lo(1),hi(1)
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + dxinv(2)*first_deriv_4(tmpy(i,j-2:j+2))
          end do

          j = hi(2)-1
          ! use 3rd-order slightly left-biased stencil
          do i=lo(1),hi(1)
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + dxinv(2)*first_deriv_l3(tmpy(i,j-2:j+1))
          end do

          j = hi(2)
          ! use completely left-biased stencil
          do i=lo(1),hi(1)
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + finhi(2)*dxinv(2)*first_deriv_lb(tmpy(i,j-3:j))
          end do
       end if
    end do

    !----- my : d()/dz -----

    do k=dlo(3),dhi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             tmpz(i,j,k) = mu(i,j,k)*qy(i,j,k,idw)
          end do
       end do
    end do

    !
    ! lo-z boundary
    !
    if (physbclo(3)) then
       k = lo(3)
       ! use completely right-biased stencil
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + foulo(3)*dxinv(3)*first_deriv_rb(tmpz(i,j,k:k+3))
          end do
       end do

       k = lo(3)+1
       ! use 3rd-order slightly right-biased stencil
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + dxinv(3)*first_deriv_r3(tmpz(i,j,k-1:k+2))
          end do
       end do

       k = lo(3)+2
       ! use 4th-order stencil
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + dxinv(3)*first_deriv_4(tmpz(i,j,k-2:k+2))
          end do
       end do

       k = lo(3)+3
       ! use 6th-order stencil
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + dxinv(3)*first_deriv_6(tmpz(i,j,k-3:k+3))
          end do
       end do
    end if

    !
    ! interior
    !
    do k=slo(3),shi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + dxinv(3)*first_deriv_8(tmpz(i,j,k-4:k+4)) 
          end do
       end do
    end do

    !
    ! hi-z boundary
    !
    if (physbchi(3)) then
       k = hi(3)-3
       ! use 6th-order stencil
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + dxinv(3)*first_deriv_6(tmpz(i,j,k-3:k+3))
          end do
       end do

       k = hi(3)-2
       ! use 4th-order stencil
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + dxinv(3)*first_deriv_4(tmpz(i,j,k-2:k+2))
          end do
       end do

       k = hi(3)-1
       ! use 3rd-order slightly left-biased stencil
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + dxinv(3)*first_deriv_l3(tmpz(i,j,k-2:k+1))
          end do
       end do

       k = hi(3)
       ! use completely left-biased stencil
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + fouhi(3)*dxinv(3)*first_deriv_lb(tmpz(i,j,k-3:k))
          end do
       end do
    end if

    !----- mz -----

    !----- mz : d()/dx -------
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)

          do i=dlo(1),dhi(1)
             tmpx(i) = mu(i,j,k)*qz(i,j,k,idu)
          end do

          !
          ! lo-x boundary
          !
          if (physbclo(1)) then
             i = lo(1)
             ! use completely right-biased stencil
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + foulo(1)*dxinv(1)*first_deriv_rb(tmpx(i:i+3))

             i = lo(1)+1
             ! use 3rd-order slightly right-biased stencil
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + dxinv(1)*first_deriv_r3(tmpx(i-1:i+2))

             i = lo(1)+2
             ! use 4th-order stencil
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + dxinv(1)*first_deriv_4(tmpx(i-2:i+2))

             i = lo(1)+3
             ! use 6th-order stencil
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + dxinv(1)*first_deriv_6(tmpx(i-3:i+3))
          end if

          !
          ! interior
          !
          do i=slo(1),shi(1)
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + dxinv(1) * first_deriv_8(tmpx(i-4:i+4))
          end do

          !
          ! hi-x boundary
          !
          if (physbchi(1)) then
             i = hi(1)-3
             ! use 6th-order stencil
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + dxinv(1)*first_deriv_6(tmpx(i-3:i+3))

             i = hi(1)-2
             ! use 4th-order stencil
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + dxinv(1)*first_deriv_4(tmpx(i-2:i+2))

             i = hi(1)-1
             ! use 3rd-order slightly left-biased stencil
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + dxinv(1)*first_deriv_l3(tmpx(i-2:i+1))

             i = hi(1)
             ! use completely left-biased stencil
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + fouhi(1)*dxinv(1)*first_deriv_lb(tmpx(i-3:i))
          end if
       end do
    end do

    !----- mz : d()/dy -----
    do k=lo(3),hi(3)

       do j=dlo(2),dhi(2)
          do i=lo(1),hi(1)
             tmpy(i,j) = mu(i,j,k)*qz(i,j,k,idv)
          end do
       end do

       !
       ! lo-y boundary
       !
       if (physbclo(2)) then
          j = lo(2)
          ! use completely right-biased stencil
          do i=lo(1),hi(1)
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + foulo(2)*dxinv(2)*first_deriv_rb(tmpy(i,j:j+3))
          end do

          j = lo(2)+1
          ! use 3rd-order slightly right-biased stencil
          do i=lo(1),hi(1)
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + dxinv(2)*first_deriv_r3(tmpy(i,j-1:j+2))
          end do

          j = lo(2)+2
          ! use 4th-order stencil
          do i=lo(1),hi(1)
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + dxinv(2)*first_deriv_4(tmpy(i,j-2:j+2))
          end do

          j = lo(2)+3
          ! use 6th-order stencil
          do i=lo(1),hi(1)
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + dxinv(2)*first_deriv_6(tmpy(i,j-3:j+3))
          end do
       end if

       !
       ! interior
       !
       do j=slo(2),shi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + dxinv(2)*first_deriv_8(tmpy(i,j-4:j+4))
          end do
       end do

       !
       ! hi-y boundary
       !
       if (physbchi(2)) then
          j = hi(2)-3
          ! use 6th-order stencil
          do i=lo(1),hi(1)
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + dxinv(2)*first_deriv_6(tmpy(i,j-3:j+3))
          end do

          j = hi(2)-2
          ! use 4th-order stencil
          do i=lo(1),hi(1)
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + dxinv(2)*first_deriv_4(tmpy(i,j-2:j+2))
          end do

          j = hi(2)-1
          ! use 3rd-order slightly left-biased stencil
          do i=lo(1),hi(1)
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + dxinv(2)*first_deriv_l3(tmpy(i,j-2:j+1))
          end do

          j = hi(2)
          ! use completely left-biased stencil
          do i=lo(1),hi(1)
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + fouhi(2)*dxinv(2)*first_deriv_lb(tmpy(i,j-3:j))
          end do
       end if
    end do

    !----- mz : d()/dy -----

    do k=dlo(3),dhi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             tmpz(i,j,k) = vsm(i,j,k)*(qx(i,j,k,idu)+qy(i,j,k,idv))
          end do
       end do
    end do

    !
    ! lo-z boundary
    !
    if (physbclo(3)) then
       k = lo(3)
       ! use completely right-biased stencil
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + finlo(3)*dxinv(3)*first_deriv_rb(tmpz(i,j,k:k+3))
          end do
       end do

       k = lo(3)+1
       ! use 3rd-order slightly right-biased stencil
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + dxinv(3)*first_deriv_r3(tmpz(i,j,k-1:k+2))
          end do
       end do

       k = lo(3)+2
       ! use 4th-order stencil
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + dxinv(3)*first_deriv_4(tmpz(i,j,k-2:k+2))
          end do
       end do

       k = lo(3)+3
       ! use 6th-order stencil
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + dxinv(3)*first_deriv_6(tmpz(i,j,k-3:k+3))
          end do
       end do
    end if
    
    !
    ! interior
    !
    do k=slo(3),shi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + dxinv(3)*first_deriv_8(tmpz(i,j,k-4:k+4))
          end do
       end do
    end do

    !
    ! hi-z boundary
    !
    if (physbchi(3)) then
       k = hi(3)-3
       ! use 6th-order stencil
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + dxinv(3)*first_deriv_6(tmpz(i,j,k-3:k+3))
          end do
       end do

       k = hi(3)-2
       ! use 4th-order stencil
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + dxinv(3)*first_deriv_4(tmpz(i,j,k-2:k+2))
          end do
       end do

       k = hi(3)-1
       ! use 3rd-order slightly left-biased stencil
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + dxinv(3)*first_deriv_l3(tmpz(i,j,k-2:k+1))
          end do
       end do

       k = hi(3)
       ! use completely left-biased stencil
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + finhi(3)*dxinv(3)*first_deriv_lb(tmpz(i,j,k-3:k))
          end do
       end do
    end if

    !----- energy -----

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             divu(i) = (qx(i,j,k,idu)+qy(i,j,k,idv)+qz(i,j,k,idw))*vsm(i,j,k)
             tauxx(i) = 2.d0*mu(i,j,k)*qx(i,j,k,idu) + divu(i)
             tauyy(i) = 2.d0*mu(i,j,k)*qy(i,j,k,idv) + divu(i)
             tauzz(i) = 2.d0*mu(i,j,k)*qz(i,j,k,idw) + divu(i)
             
             ! change in internal energy
             rhs(i,j,k,iene) = tauxx(i)*qx(i,j,k,idu) &
                  +            tauyy(i)*qy(i,j,k,idv) &
                  +            tauzz(i)*qz(i,j,k,idw) &
                  + mu(i,j,k)*((qy(i,j,k,idu)+qx(i,j,k,idv))**2 &
                  &          + (qx(i,j,k,idw)+qz(i,j,k,idu))**2 &
                  &          + (qz(i,j,k,idv)+qy(i,j,k,idw))**2 )

          end do
       end do
    end do

    deallocate(tmpx,tmpy,tmpz)
    deallocate(vsm)

    ! grad T, grad p, and grad X

   do k=lo(3),hi(3)
       do j=lo(2),hi(2)

          !
          ! lo-x boundary
          !
          if (physbclo(1)) then
             i = lo(1)
             ! use completely right-biased stencil
             qx(i,j,k,idp) = dxinv(1) * first_deriv_rb( q(i:i+3,j,k,qpres) )
             qx(i,j,k,idT) = dxinv(1) * first_deriv_rb( q(i:i+3,j,k,qtemp) )
             
             i = lo(1)+1
             ! use 3rd-order slightly right-biased stencil
             qx(i,j,k,idp) = dxinv(1) * first_deriv_r3( q(i-1:i+2,j,k,qpres) )
             qx(i,j,k,idT) = dxinv(1) * first_deriv_r3( q(i-1:i+2,j,k,qtemp) )

             i = lo(1)+2
             ! use 4th-order stencil
             qx(i,j,k,idp) = dxinv(1) * first_deriv_4( q(i-2:i+2,j,k,qpres) )
             qx(i,j,k,idT) = dxinv(1) * first_deriv_4( q(i-2:i+2,j,k,qtemp) )

             i = lo(1)+3
             ! use 6th-order stencil
             qx(i,j,k,idp) = dxinv(1) * first_deriv_6( q(i-3:i+3,j,k,qpres) )
             qx(i,j,k,idT) = dxinv(1) * first_deriv_6( q(i-3:i+3,j,k,qtemp) )
          end if

          !
          ! interior
          !
          do i=slo(1),shi(1)
             qx(i,j,k,idp) = dxinv(1) * first_deriv_8( q(i-4:i+4,j,k,qpres) )
             qx(i,j,k,idT) = dxinv(1) * first_deriv_8( q(i-4:i+4,j,k,qtemp) )
          enddo

          !
          ! hi-x boundary
          !
          if (physbchi(1)) then
             i = hi(1)-3
             ! use 6th-order stencil
             qx(i,j,k,idp) = dxinv(1) * first_deriv_6( q(i-3:i+3,j,k,qpres) )
             qx(i,j,k,idT) = dxinv(1) * first_deriv_6( q(i-3:i+3,j,k,qtemp) )

             i = hi(1)-2
             ! use 4th-order stencil
             qx(i,j,k,idp) = dxinv(1) * first_deriv_4( q(i-2:i+2,j,k,qpres) )
             qx(i,j,k,idT) = dxinv(1) * first_deriv_4( q(i-2:i+2,j,k,qtemp) )

             i = hi(1)-1
             ! use 3rd-order slightly left-biased stencil
             qx(i,j,k,idp) = dxinv(1) * first_deriv_l3( q(i-2:i+1,j,k,qpres) )
             qx(i,j,k,idT) = dxinv(1) * first_deriv_l3( q(i-2:i+1,j,k,qtemp) )

             i = hi(1)
             ! use completely left-biased stencil
             qx(i,j,k,idp) = dxinv(1) * first_deriv_lb( q(i-3:i,j,k,qpres) )
             qx(i,j,k,idT) = dxinv(1) * first_deriv_lb( q(i-3:i,j,k,qtemp) )
          end if

       enddo
    enddo
     
    do k=lo(3),hi(3)

       !
       ! lo-y boundary
       !
       if (physbclo(2)) then
          j = lo(2)
          ! use completely right-biased stencil
          do i=lo(1),hi(1)
             qy(i,j,k,idp) = dxinv(2) * first_deriv_rb( q(i,j:j+3,k,qpres) )
             qy(i,j,k,idT) = dxinv(2) * first_deriv_rb( q(i,j:j+3,k,qtemp) )
          enddo

          j = lo(2)+1
          ! use 3rd-order slightly right-biased stencil
          do i=lo(1),hi(1)
             qy(i,j,k,idp) = dxinv(2) * first_deriv_r3( q(i,j-1:j+2,k,qpres) )
             qy(i,j,k,idT) = dxinv(2) * first_deriv_r3( q(i,j-1:j+2,k,qtemp) )
          enddo

          j = lo(2)+2
          ! use 4th-order stencil
          do i=lo(1),hi(1)
             qy(i,j,k,idp) = dxinv(2) * first_deriv_4( q(i,j-2:j+2,k,qpres) )
             qy(i,j,k,idT) = dxinv(2) * first_deriv_4( q(i,j-2:j+2,k,qtemp) )
          enddo

          j = lo(2)+3
          ! use 6th-order stencil
          do i=lo(1),hi(1)
             qy(i,j,k,idp) = dxinv(2) * first_deriv_6( q(i,j-3:j+3,k,qpres) )
             qy(i,j,k,idT) = dxinv(2) * first_deriv_6( q(i,j-3:j+3,k,qtemp) )
          enddo
       end if

       !
       ! interior
       !
       do j=slo(2),shi(2)
          do i=lo(1),hi(1)
             qy(i,j,k,idp) = dxinv(2) * first_deriv_8( q(i,j-4:j+4,k,qpres) )
             qy(i,j,k,idT) = dxinv(2) * first_deriv_8( q(i,j-4:j+4,k,qtemp) )
          enddo
       enddo

       !
       ! hi-y boundary
       !
       if (physbchi(2)) then
          j = hi(2)-3
          ! use 6th-order stencil
          do i=lo(1),hi(1)
             qy(i,j,k,idp) = dxinv(2) * first_deriv_6( q(i,j-3:j+3,k,qpres) )
             qy(i,j,k,idT) = dxinv(2) * first_deriv_6( q(i,j-3:j+3,k,qtemp) )
          enddo

          j = hi(2)-2
          ! use 4th-order stencil
          do i=lo(1),hi(1)
             qy(i,j,k,idp) = dxinv(2) * first_deriv_4( q(i,j-2:j+2,k,qpres) )
             qy(i,j,k,idT) = dxinv(2) * first_deriv_4( q(i,j-2:j+2,k,qtemp) )
          enddo

          j = hi(2)-1
          ! use 3rd-order slightly left-biased stencil
          do i=lo(1),hi(1)
             qy(i,j,k,idp) = dxinv(2) * first_deriv_l3( q(i,j-2:j+1,k,qpres) )
             qy(i,j,k,idT) = dxinv(2) * first_deriv_l3( q(i,j-2:j+1,k,qtemp) )
          enddo

          j = hi(2)
          ! use completely left-biased stencil
          do i=lo(1),hi(1)
             qy(i,j,k,idp) = dxinv(2) * first_deriv_lb( q(i,j-3:j,k,qpres) )
             qy(i,j,k,idT) = dxinv(2) * first_deriv_lb( q(i,j-3:j,k,qtemp) )
          enddo
       end if

    enddo

    ! lo-z boundary
    if (physbclo(3)) then
       k = lo(3)
       ! use completely right-biased stencil
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             qz(i,j,k,idp) = dxinv(3) * first_deriv_rb( q(i,j,k:k+3,qpres) )
             qz(i,j,k,idT) = dxinv(3) * first_deriv_rb( q(i,j,k:k+3,qtemp) )
          enddo
       enddo

       k = lo(3)+1
       ! use 3rd-order slightly right-biased stencil
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             qz(i,j,k,idp) = dxinv(3) * first_deriv_r3( q(i,j,k-1:k+2,qpres) )
             qz(i,j,k,idT) = dxinv(3) * first_deriv_r3( q(i,j,k-1:k+2,qtemp) )
          enddo
       enddo

       k = lo(3)+2
       ! use 4th-order stencil
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             qz(i,j,k,idp) = dxinv(3) * first_deriv_4( q(i,j,k-2:k+2,qpres) )
             qz(i,j,k,idT) = dxinv(3) * first_deriv_4( q(i,j,k-2:k+2,qtemp) )
          enddo
       enddo

       k = lo(3)+3
       ! use 6th-order stencil
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             qz(i,j,k,idp) = dxinv(3) * first_deriv_6( q(i,j,k-3:k+3,qpres) )
             qz(i,j,k,idT) = dxinv(3) * first_deriv_6( q(i,j,k-3:k+3,qtemp) )
          enddo
       enddo
    end if

    !
    ! interior
    !
    do k=slo(3),shi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             qz(i,j,k,idp) = dxinv(3) * first_deriv_8( q(i,j,k-4:k+4,qpres) )
             qz(i,j,k,idT) = dxinv(3) * first_deriv_8( q(i,j,k-4:k+4,qtemp) )
          enddo
       enddo
    enddo

    !
    ! hi-z boundary
    !
    if (physbchi(3)) then
       k = hi(3)-3
       ! use 6th-order stencil
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             qz(i,j,k,idp) = dxinv(3) * first_deriv_6( q(i,j,k-3:k+3,qpres) )
             qz(i,j,k,idT) = dxinv(3) * first_deriv_6( q(i,j,k-3:k+3,qtemp) )
          enddo
       enddo

       k = hi(3)-2
       ! use 4th-order stencil
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             qz(i,j,k,idp) = dxinv(3) * first_deriv_4( q(i,j,k-2:k+2,qpres) )
             qz(i,j,k,idT) = dxinv(3) * first_deriv_4( q(i,j,k-2:k+2,qtemp) )
          enddo
       enddo

       k = hi(3)-1
       ! use 3rd-order slightly left-biased stencil
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             qz(i,j,k,idp) = dxinv(3) * first_deriv_l3( q(i,j,k-2:k+1,qpres) )
             qz(i,j,k,idT) = dxinv(3) * first_deriv_l3( q(i,j,k-2:k+1,qtemp) )
          enddo
       enddo

       k = hi(3)
       ! use completely left-biased stencil
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             qz(i,j,k,idp) = dxinv(3) * first_deriv_lb( q(i,j,k-3:k,qpres) )
             qz(i,j,k,idT) = dxinv(3) * first_deriv_lb( q(i,j,k-3:k,qtemp) )
          enddo
       enddo
    end if

    do n=1,nspecies
       qxn = qx1 + n - 1
       qdxn = idX1 + n -1

       do k=lo(3),hi(3)
          do j=lo(2),hi(2)

             !
             ! lo-x boundary
             !
             if (physbclo(1)) then
                i = lo(1)
                ! use completely right-biased stencil
                qx(i,j,k,qdxn) = dxinv(1) * first_deriv_rb( q(i:i+3,j,k,qxn) )

                i = lo(1)+1
                ! use 3rd-order slightly right-biased stencil
                qx(i,j,k,qdxn) = dxinv(1) * first_deriv_r3( q(i-1:i+2,j,k,qxn) )

                i = lo(1)+2
                ! use 4th-order stencil
                qx(i,j,k,qdxn) = dxinv(1) * first_deriv_4( q(i-2:i+2,j,k,qxn) )
                
                i = lo(1)+3
                ! use 6th-order stencil
                qx(i,j,k,qdxn) = dxinv(1) * first_deriv_6( q(i-3:i+3,j,k,qxn) )
             end if

             !
             ! interior
             !
             do i=slo(1),shi(1)
                qx(i,j,k,qdxn) = dxinv(1) * first_deriv_8( q(i-4:i+4,j,k,qxn) )
             enddo

             !
             ! hi-x boundary
             !
             if (physbchi(1)) then
                i = hi(1)-3
                ! use 6th-order stencil
                qx(i,j,k,qdxn) = dxinv(1) * first_deriv_6( q(i-3:i+3,j,k,qxn) )

                i = hi(1)-2
                ! use 4th-order stencil
                qx(i,j,k,qdxn) = dxinv(1) * first_deriv_4( q(i-2:i+2,j,k,qxn) )
                
                i = hi(1)-1
                ! use 3rd-order slightly left-biased stencil
                qx(i,j,k,qdxn) = dxinv(1) * first_deriv_l3( q(i-2:i+1,j,k,qxn) )

                i = hi(1)
                ! use completely left-biased stencil
                qx(i,j,k,qdxn) = dxinv(1) * first_deriv_lb( q(i-3:i,j,k,qxn) )
             end if

          enddo
       enddo

       do k=lo(3),hi(3)

          !
          ! lo-y boundary
          !
          if (physbclo(2)) then
             j = lo(2)
             ! use completely right-biased stencil
             do i=lo(1),hi(1)
                qy(i,j,k,qdxn) = dxinv(2) * first_deriv_rb( q(i,j:j+3,k,qxn) )
             enddo
             
             j = lo(2)+1
             ! use 3rd-order slightly right-biased stencil
             do i=lo(1),hi(1)
                qy(i,j,k,qdxn) = dxinv(2) * first_deriv_r3( q(i,j-1:j+2,k,qxn) )
             enddo

             j = lo(2)+2
             ! use 4th-order stencil
             do i=lo(1),hi(1)
                qy(i,j,k,qdxn) = dxinv(2) * first_deriv_4( q(i,j-2:j+2,k,qxn) )
             enddo

             j = lo(2)+3
             ! use 6th-order stencil
             do i=lo(1),hi(1)
                qy(i,j,k,qdxn) = dxinv(2) * first_deriv_6( q(i,j-3:j+3,k,qxn) )
             enddo
          end if

          !
          ! interior
          !
          do j=slo(2),shi(2)
             do i=lo(1),hi(1)
                qy(i,j,k,qdxn) = dxinv(2) * first_deriv_8( q(i,j-4:j+4,k,qxn) )
             enddo
          enddo

          !
          ! hi-y boundary
          !
          if (physbchi(2)) then
             j = hi(2)-3
             ! use 6th-order stencil
             do i=lo(1),hi(1)
                qy(i,j,k,qdxn) = dxinv(2) * first_deriv_6( q(i,j-3:j+3,k,qxn) )
             enddo

             j = hi(2)-2
             ! use 4th-order stencil
             do i=lo(1),hi(1)
                qy(i,j,k,qdxn) = dxinv(2) * first_deriv_4( q(i,j-2:j+2,k,qxn) )
             enddo

             j = hi(2)-1
             ! use 3rd-order slightly left-biased stencil
             do i=lo(1),hi(1)
                qy(i,j,k,qdxn) = dxinv(2) * first_deriv_l3( q(i,j-2:j+1,k,qxn) )
             enddo

             j = hi(2)
             ! use completely left-biased stencil
             do i=lo(1),hi(1)
                qy(i,j,k,qdxn) = dxinv(2) * first_deriv_lb( q(i,j-3:j,k,qxn) )
             enddo
          end if

       enddo

       ! lo-z boundary
       if (physbclo(3)) then
          k = lo(3)
          ! use completely right-biased stencil
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                qz(i,j,k,qdxn) = dxinv(3) * first_deriv_rb( q(i,j,k:k+3,qxn) )
             enddo
          enddo

          k = lo(3)+1
          ! use 3rd-order slightly right-biased stencil
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                qz(i,j,k,qdxn) = dxinv(3) * first_deriv_r3( q(i,j,k-1:k+2,qxn) )
             enddo
          enddo

          k = lo(3)+2
          ! use 4th-order stencil
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                qz(i,j,k,qdxn) = dxinv(3) * first_deriv_4( q(i,j,k-2:k+2,qxn) )
             enddo
          enddo

          k = lo(3)+3
          ! use 6th-order stencil
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                qz(i,j,k,qdxn) = dxinv(3) * first_deriv_6( q(i,j,k-3:k+3,qxn) )
             enddo
          enddo
       end if

       !
       ! interior
       !
       do k=slo(3),shi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                qz(i,j,k,qdxn) = dxinv(3) * first_deriv_8( q(i,j,k-4:k+4,qxn) )
             enddo
          enddo
       enddo

       !
       ! hi-z boundary
       !
       if (physbchi(3)) then
          k = hi(3)-3
          ! use 6th-order stencil
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                qz(i,j,k,qdxn) = dxinv(3) * first_deriv_6( q(i,j,k-3:k+3,qxn) )
             enddo
          enddo

          k = hi(3)-2
          ! use 4th-order stencil
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                qz(i,j,k,qdxn) = dxinv(3) * first_deriv_4( q(i,j,k-2:k+2,qxn) )
             enddo
          enddo

          k = hi(3)-1
          ! use 3rd-order slightly left-biased stencil
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                qz(i,j,k,qdxn) = dxinv(3) * first_deriv_l3( q(i,j,k-2:k+1,qxn) )
             enddo
          enddo

          k = hi(3)
          ! use completely left-biased stencil
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                qz(i,j,k,qdxn) = dxinv(3) * first_deriv_lb( q(i,j,k-3:k,qxn) )
             enddo
          enddo          
       end if
    enddo
    
  end subroutine s3d_diffterm_1_3d


  subroutine s3d_diffterm_2_3d (lo,hi,dx,q,qlo,qhi,r_g,glo,ghi,rhs,rlo,rhi, &
       qx,qxlo,qxhi,qy,qylo,qyhi,qz,qzlo,qzhi, &
       mu,xi,lam,dxy,dlo_g,dhi_g,bclo,bchi)

    integer,         intent(in):: dlo_g(3),dhi_g(3),bclo(3),bchi(3)
    integer,         intent(in):: lo(3),hi(3),qlo(3),qhi(3),glo(3),ghi(3),rlo(3),rhi(3)
    integer,         intent(in):: qxlo(3),qxhi(3),qylo(3),qyhi(3),qzlo(3),qzhi(3)
    double precision,intent(in):: dx(3)
    double precision,intent(in)   ::q  ( qlo(1): qhi(1), qlo(2): qhi(2), qlo(3): qhi(3),nprim)
    double precision,intent(in)   ::qx (qxlo(1):qxhi(1),qxlo(2):qxhi(2),qxlo(3):qxhi(3),ndq)
    double precision,intent(in)   ::qy (qylo(1):qyhi(1),qylo(2):qyhi(2),qylo(3):qyhi(3),ndq)
    double precision,intent(in)   ::qz (qzlo(1):qzhi(1),qzlo(2):qzhi(2),qzlo(3):qzhi(3),ndq)
    double precision,intent(in)   ::mu ( qlo(1): qhi(1), qlo(2): qhi(2), qlo(3): qhi(3))
    double precision,intent(in)   ::xi ( qlo(1): qhi(1), qlo(2): qhi(2), qlo(3): qhi(3))
    double precision,intent(in)   ::lam( qlo(1): qhi(1), qlo(2): qhi(2), qlo(3): qhi(3))
    double precision,intent(in)   ::dxy( qlo(1): qhi(1), qlo(2): qhi(2), qlo(3): qhi(3),nspecies)
    double precision,intent(inout)::r_g( glo(1): ghi(1), glo(2): ghi(2), glo(3): ghi(3),ncons)
    double precision,intent(inout)::rhs( rlo(1): rhi(1), rlo(2): rhi(2), rlo(3): rhi(3),ncons)

    double precision, allocatable, dimension(:,:,:) :: vp, dpe, FE
    double precision, allocatable, dimension(:,:,:,:) :: dpy, FY
    ! dxy: diffusion coefficient of X in equation for Y
    ! dpy: diffusion coefficient of p in equation for Y
    ! NOT USING ! dxe: diffusion coefficient of X in equation for energy
    ! dpe: diffusion coefficient of p in equation for energy

    double precision :: dxinv(3)
    integer          :: i,j,k,n, qxn, qyn, qhn, idXn, iryn
    integer :: slo(3), shi(3), dlo(3), dhi(3)
    double precision, allocatable, dimension(:,:,:) :: rvc
    double precision, allocatable :: tmpx(:), tmpy(:,:),tmpz(:,:,:)

    ! used to turn off some terms
    double precision :: finlo(3), finhi(3)
    double precision :: foulo(3), fouhi(3)

    logical :: physbclo(3), physbchi(3)

    ! Only the region bounded by [dlo_g,dhi_g] contains good data.
    ! [slo,shi] will be safe for 8th-order stencil
    do i=1,3
       dlo(i) = max(lo(i)-stencil_ng, dlo_g(i))
       dhi(i) = min(hi(i)+stencil_ng, dhi_g(i))
       slo(i) = dlo(i) + stencil_ng
       shi(i) = dhi(i) - stencil_ng

       if (dlo(i) .eq. lo(i)) then
          physbclo(i) = .true.
       else
          physbclo(i) = .false.
       end if

       if (dhi(i) .eq. hi(i)) then
          physbchi(i) = .true.
       else
          physbchi(i) = .false.
       end if
    end do

    finlo = 1.d0 
    finhi = 1.d0
    foulo = 1.d0 
    fouhi = 1.d0

    do i=1,3
       if (physbclo(i)) then 
          if (bclo(i) .eq. INLET) then
             finlo(i) = 0.d0
          else if (bclo(i) .eq. OUTLET) then
             foulo(i) = 0.d0
          end if
       end if

       if (physbchi(i)) then 
          if (bchi(i) .eq. INLET) then
             finhi(i) = 0.d0
          else if (bchi(i) .eq. OUTLET) then
             fouhi(i) = 0.d0
          end if
       end if
    end do

    do i = 1,3
       dxinv(i) = 1.0d0 / dx(i)
    end do

    allocate(vp(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3)))

    allocate(dpy(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),nspecies))
    allocate(dpe(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3)))

    allocate(FY(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),nspecies))
    allocate(FE(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3)))

    allocate(rvc (dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3)))

    allocate(tmpx(dlo(1):dhi(1)))
    allocate(tmpy( lo(1): hi(1),dlo(2):dhi(2)))
    allocate(tmpz( lo(1): hi(1), lo(2): hi(2),dlo(3):dhi(3)))

    do k=dlo(3),dhi(3)
       do j=dlo(2),dhi(2)
          do i=dlo(1),dhi(1)
             vp(i,j,k) = xi(i,j,k) + FourThirds*mu(i,j,k)
          enddo
       enddo
    enddo

    dpe = 0.d0

    do n=1,nspecies
       qxn = qx1+n-1
       qyn = qy1+n-1
       qhn = qh1+n-1
       do k=dlo(3),dhi(3)
          do j=dlo(2),dhi(2)
             do i=dlo(1),dhi(1)
                dpy(i,j,k,n) = dxy(i,j,k,n)/q(i,j,k,qpres)*(q(i,j,k,qxn)-q(i,j,k,qyn))
                dpe(i,j,k) = dpe(i,j,k) + dpy(i,j,k,n)*q(i,j,k,qhn)
             end do
          end do
       end do
    end do

    ! x-direction

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)

          do i=lo(1)-4,hi(1)+4
             tmpx(i) = vp(i,j,k)*qx(i,j,k,idu) 
          end do
          do i=lo(1),hi(1)
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + dxinv(1) * first_deriv_8(tmpx(i-4:i+4))
          end do

          do i=lo(1)-4,hi(1)+4
             tmpx(i) = mu(i,j,k)*qx(i,j,k,idv)
          end do
          do i=lo(1),hi(1)
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + dxinv(1) * first_deriv_8(tmpx(i-4:i+4))
          end do

          do i=lo(1)-4,hi(1)+4
             tmpx(i) = mu(i,j,k)*qx(i,j,k,idw)
          end do
          do i=lo(1),hi(1)
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + dxinv(1) * first_deriv_8(tmpx(i-4:i+4))
          end do

          do i=lo(1)-4,hi(1)+4
             tmpx(i) = lam(i,j,k)*qx(i,j,k,idT)
          end do
          do i=lo(1),hi(1)
             rhs(i,j,k,iene) = rhs(i,j,k,iene) + dxinv(1) * first_deriv_8(tmpx(i-4:i+4))
          end do

       end do
    end do

    ! y-direction

    do k=lo(3),hi(3)

       do j=lo(2)-4,hi(2)+4
          do i=lo(1),hi(1)
             tmpy(i,j) = mu(i,j,k)*qy(i,j,k,idu)
          end do
       end do
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + dxinv(2) * first_deriv_8(tmpy(i,j-4:j+4))
          end do
       end do

       do j=lo(2)-4,hi(2)+4
          do i=lo(1),hi(1)
             tmpy(i,j) = vp(i,j,k)*qy(i,j,k,idv)
          end do
       end do
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + dxinv(2) * first_deriv_8(tmpy(i,j-4:j+4))
          end do
       end do

       do j=lo(2)-4,hi(2)+4
          do i=lo(1),hi(1)
             tmpy(i,j) = mu(i,j,k)*qy(i,j,k,idw)
          end do
       end do
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + dxinv(2) * first_deriv_8(tmpy(i,j-4:j+4))
          end do
       end do

       do j=lo(2)-4,hi(2)+4
          do i=lo(1),hi(1)
             tmpy(i,j) = lam(i,j,k)*qy(i,j,k,idT)
          end do
       end do
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,iene) = rhs(i,j,k,iene) + dxinv(2) * first_deriv_8(tmpy(i,j-4:j+4))
          end do
       end do

    end do

    do k=lo(3)-4,hi(3)+4
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             tmpz(i,j,k) = mu(i,j,k)*qz(i,j,k,idu)
          end do
       end do
    end do

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imx) = rhs(i,j,k,imx) + dxinv(3) * first_deriv_8(tmpz(i,j,k-4:k+4))
          end do
       end do
    end do

    do k=lo(3)-4,hi(3)+4
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             tmpz(i,j,k) = mu(i,j,k)*qz(i,j,k,idv)
          end do
       end do
    end do

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imy) = rhs(i,j,k,imy) + dxinv(3) * first_deriv_8(tmpz(i,j,k-4:k+4))
          end do
       end do
    end do

    do k=lo(3)-4,hi(3)+4
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             tmpz(i,j,k) = vp(i,j,k)*qz(i,j,k,idw)
          end do
       end do
    end do

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imz) = rhs(i,j,k,imz) + dxinv(3) * first_deriv_8(tmpz(i,j,k-4:k+4))
          end do
       end do
    end do

    do k=lo(3)-4,hi(3)+4
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             tmpz(i,j,k) = lam(i,j,k)*qz(i,j,k,idT)
          end do
       end do
    end do

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,iene) = rhs(i,j,k,iene) + dxinv(3) * first_deriv_8(tmpz(i,j,k-4:k+4))
          end do
       end do
    end do

    ! add kinetic energy
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,iene) = rhs(i,j,k,iene) &
                  + rhs(i,j,k,imx)*q(i,j,k,qu) &
                  + rhs(i,j,k,imy)*q(i,j,k,qv) &
                  + rhs(i,j,k,imz)*q(i,j,k,qw)
          end do
       end do
    end do

    ! x-direction
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=dlo(1),dhi(1)
             FE(i,j,k) = dpe(i,j,k) * qx(i,j,k,idp)
          end do
       end do
    end do

    rvc = 0.d0

    do n=1,nspecies
       idXn = idX1+n-1
       qhn = qh1+n-1
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=dlo(1),dhi(1)
                FE(i,j,k) = FE(i,j,k) + dxy(i,j,k,n)*qx(i,j,k,idXn)*q(i,j,k,qhn)
                FY(i,j,k,n) = dxy(i,j,k,n)*qx(i,j,k,idXn) + dpy(i,j,k,n)*qx(i,j,k,idp)
                rvc(i,j,k) = rvc(i,j,k) + FY(i,j,k,n)
             end do
          end do
       end do
    end do

    do n=1,nspecies
       qyn = qy1+n-1
       qhn = qh1+n-1
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=dlo(1),dhi(1)
                FY(i,j,k,n) = FY(i,j,k,n) - rvc(i,j,k)*q(i,j,k,qyn)
                FE(i,j,k) = FE(i,j,k) - rvc(i,j,k)*q(i,j,k,qyn)*q(i,j,k,qhn)
             end do
          end do
       end do
    end do

    do n=1,nspecies    
       iryn = iry1+n-1
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                rhs(i,j,k,iryn) = rhs(i,j,k,iryn) + &
                     dxinv(1) * first_deriv_8( FY(i-4:i+4,j,k,n) )
             end do
          end do
       end do
    end do
    
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,iene) = rhs(i,j,k,iene) + &
                  dxinv(1) * first_deriv_8( FE(i-4:i+4,j,k) )
          end do
       end do
    end do

    ! y-direction
    do k=lo(3),hi(3)
       do j=dlo(2),dhi(2)
          do i=lo(1),hi(1)
             FE(i,j,k) = dpe(i,j,k) * qy(i,j,k,idp)
          end do
       end do
    end do

    rvc = 0.d0

    do n=1,nspecies
       idXn = idX1+n-1
       qhn = qh1+n-1
       do k=lo(3),hi(3)
          do j=dlo(2),dhi(2)
             do i=lo(1),hi(1)
                FE(i,j,k) = FE(i,j,k) + dxy(i,j,k,n)*qy(i,j,k,idXn)*q(i,j,k,qhn)
                FY(i,j,k,n) = dxy(i,j,k,n)*qy(i,j,k,idXn) + dpy(i,j,k,n)*qy(i,j,k,idp)
                rvc(i,j,k) = rvc(i,j,k) + FY(i,j,k,n)                
             end do
          end do
       end do
    end do

    do n=1,nspecies
       qyn = qy1+n-1
       qhn = qh1+n-1
       do k=lo(3),hi(3)
          do j=dlo(2),dhi(2)
             do i=lo(1),hi(1)
                FY(i,j,k,n) = FY(i,j,k,n) - rvc(i,j,k)*q(i,j,k,qyn)
                FE(i,j,k) = FE(i,j,k) - rvc(i,j,k)*q(i,j,k,qyn)*q(i,j,k,qhn)                
             end do
          end do
       end do
    end do
    
    do n=1,nspecies    
       iryn = iry1+n-1
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                rhs(i,j,k,iryn) = rhs(i,j,k,iryn) + &
                     dxinv(2) * first_deriv_8( FY(i,j-4:j+4,k,n) )
             end do
          end do
       end do
    end do
    
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,iene) = rhs(i,j,k,iene) + &
                  dxinv(2) * first_deriv_8( FE(i,j-4:j+4,k) )
          end do
       end do
    end do

    ! z-direction
    do k=dlo(3),dhi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             FE(i,j,k) = dpe(i,j,k) * qz(i,j,k,idp)
          end do
       end do
    end do

    rvc = 0.d0

    do n=1,nspecies
       idXn = idX1+n-1
       qhn = qh1+n-1
       do k=dlo(3),dhi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                FE(i,j,k) = FE(i,j,k) + dxy(i,j,k,n)*qz(i,j,k,idXn)*q(i,j,k,qhn)
                FY(i,j,k,n) = dxy(i,j,k,n)*qz(i,j,k,idXn) + dpy(i,j,k,n)*qz(i,j,k,idp)
                rvc(i,j,k) = rvc(i,j,k) + FY(i,j,k,n)                
             end do
          end do
       end do
    end do

    do n=1,nspecies
       qyn = qy1+n-1
       qhn = qh1+n-1
       do k=dlo(3),dhi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                FY(i,j,k,n) = FY(i,j,k,n) - rvc(i,j,k)*q(i,j,k,qyn)
                FE(i,j,k) = FE(i,j,k) - rvc(i,j,k)*q(i,j,k,qyn)*q(i,j,k,qhn)
             end do
          end do
       end do
    end do

    do n=1,nspecies    
       iryn = iry1+n-1
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
            do i=lo(1),hi(1)
                rhs(i,j,k,iryn) = rhs(i,j,k,iryn) + &
                     dxinv(3) * first_deriv_8( FY(i,j,k-4:k+4,n) )
             end do
          end do
       end do
    end do

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,iene) = rhs(i,j,k,iene) + &
                  dxinv(3) * first_deriv_8( FE(i,j,k-4:k+4) )
          end do
       end do
    end do

    deallocate(vp,dpy,dpe,FY,FE,rvc,tmpx,tmpy,tmpz)

    r_g       (lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:) = &
         r_g  (lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:) &
         + rhs(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:)

  end subroutine s3d_diffterm_2_3d

end module kernels_s3d_module
