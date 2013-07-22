module kernels_module
  use bc_module
  use chemistry_module, only : nspecies, molecular_weight, Ru
  use derivative_stencil_module, only : stencil_ng, first_deriv_8, first_deriv_6, &
       first_deriv_4, first_deriv_l3, first_deriv_r3, first_deriv_rb, first_deriv_lb, &
       M8, M8T, M6, M6T, M4, M4T, M2, BRB, BLB, D8, D6, D4
  use variables_module
  implicit none

  private

  public :: hypterm_3d, narrow_diffterm_3d, chemterm_3d, comp_courno_3d

contains

  subroutine hypterm_3d (lo,hi,dx,cons,clo,chi,q,qlo,qhi,rhs_g,rlo,rhi,dlo_g,dhi_g,bclo,bchi)

    integer,         intent(in):: dlo_g(3),dhi_g(3),bclo(3),bchi(3)
    integer,         intent(in):: lo(3),hi(3),clo(3),chi(3),qlo(3),qhi(3),rlo(3),rhi(3)
    double precision,intent(in):: dx(3)
    double precision,intent(in):: cons(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3),ncons)
    double precision,intent(in)::    q(qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3),nprim)
    double precision           ::rhs_g(rlo(1):rhi(1),rlo(2):rhi(2),rlo(3):rhi(3),ncons)

    integer          :: i,j,k,n
    double precision :: dxinv(3)
    double precision :: un(-4:4)
    integer :: slo(3), shi(3), dlo(3), dhi(3)
    
    double precision, allocatable :: tmpx(:), tmpy(:,:),tmpz(:,:,:)
    double precision, allocatable :: rhs(:,:,:,:)

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

    do i=1,3
       dxinv(i) = 1.0d0 / dx(i)
    end do

    allocate(tmpx(dlo(1):dhi(1)))
    allocate(tmpy(lo(1) : hi(1),dlo(2):dhi(2)))
    allocate(tmpz(lo(1) : hi(1), lo(2): hi(2),dlo(3):dhi(3)))

    allocate(rhs(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),ncons))
    rhs = 0.d0

    ! ------- BEGIN x-direction -------

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)

          do i=slo(1),shi(1)
             rhs(i,j,k,irho) = rhs(i,j,k,irho) - dxinv(1) * &
                  first_deriv_8( cons(i-4:i+4,j,k,imx) ) 
          end do

          do i=dlo(1),dhi(1)
             tmpx(i) = cons(i,j,k,imx)*q(i,j,k,qu)+q(i,j,k,qpres)
          end do
          do i=slo(1),shi(1)
             rhs(i,j,k,imx) = rhs(i,j,k,imx) - dxinv(1) * first_deriv_8(tmpx(i-4:i+4))
          end do

          do i=dlo(1),dhi(1)
             tmpx(i) = cons(i,j,k,imy)*q(i,j,k,qu)
          end do
          do i=slo(1),shi(1)
             rhs(i,j,k,imy) = rhs(i,j,k,imy) - dxinv(1) * first_deriv_8(tmpx(i-4:i+4))
          end do

          do i=dlo(1),dhi(1)
             tmpx(i) = cons(i,j,k,imz)*q(i,j,k,qu)
          end do
          do i=slo(1),shi(1)
             rhs(i,j,k,imz) = rhs(i,j,k,imz) - dxinv(1) * first_deriv_8(tmpx(i-4:i+4))
          end do

          do i=dlo(1),dhi(1)
             tmpx(i) = (cons(i,j,k,iene)+q(i,j,k,qpres))*q(i,j,k,qu)
          end do
          do i=slo(1),shi(1)
             rhs(i,j,k,iene) = rhs(i,j,k,iene) - dxinv(1) * first_deriv_8(tmpx(i-4:i+4))
          end do

       enddo
    enddo

    do n = iry1, iry1+nspecies-1
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
    
             do i=dlo(1),dhi(1)
                tmpx(i) = cons(i,j,k,n)*q(i,j,k,qu)
             end do
             do i=slo(1),shi(1)
                rhs(i,j,k,n) = rhs(i,j,k,n) - dxinv(1) * first_deriv_8(tmpx(i-4:i+4))
             end do

          enddo
       enddo
    enddo

    ! ------- END x-direction -------

    ! ------- BEGIN y-direction -------

    do k=lo(3),hi(3)
       do j=slo(2),shi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,irho)=rhs(i,j,k,irho) - dxinv(2) * &
                  first_deriv_8( cons(i,j-4:j+4,k,imy) )
          end do
       end do

       do j=dlo(2),dhi(2)
          do i=lo(1),hi(1)
             tmpy(i,j) = cons(i,j,k,imx)*q(i,j,k,qv)
          end do
       end do
       do j=slo(2),shi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imx)=rhs(i,j,k,imx) - dxinv(2) * first_deriv_8(tmpy(i,j-4:j+4))
          enddo
       enddo

       do j=dlo(2),dhi(2)
          do i=lo(1),hi(1)
             tmpy(i,j) = cons(i,j,k,imy)*q(i,j,k,qv)+q(i,j,k,qpres)
          end do
       end do
       do j=slo(2),shi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imy)=rhs(i,j,k,imy) - dxinv(2) * first_deriv_8(tmpy(i,j-4:j+4))
          enddo
       enddo

       do j=dlo(2),dhi(2)
          do i=lo(1),hi(1)
             tmpy(i,j) = cons(i,j,k,imz)*q(i,j,k,qv)
          end do
       end do
       do j=slo(2),shi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imz)=rhs(i,j,k,imz) - dxinv(2) * first_deriv_8(tmpy(i,j-4:j+4))
          enddo
       enddo

       do j=dlo(2),dhi(2)
          do i=lo(1),hi(1)
             tmpy(i,j) = (cons(i,j,k,iene)+q(i,j,k,qpres))*q(i,j,k,qv)
          end do
       end do
       do j=slo(2),shi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,iene)=rhs(i,j,k,iene) - dxinv(2) * first_deriv_8(tmpy(i,j-4:j+4))
          enddo
       enddo
    enddo

    do n = iry1, iry1+nspecies-1
       do k=lo(3),hi(3)
          do j=dlo(2),dhi(2)
             do i=lo(1),hi(1)
                tmpy(i,j) = cons(i,j,k,n)*q(i,j,k,qv)
             end do
          end do
          do j=slo(2),shi(2)
             do i=lo(1),hi(1)
                rhs(i,j,k,n) = rhs(i,j,k,n) - dxinv(2) * first_deriv_8(tmpy(i,j-4:j+4))
             end do
          enddo
       enddo
    enddo

    ! ------- END y-direction -------

    ! ------- BEGIN z-direction -------

    do k=slo(3),shi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,irho)=rhs(i,j,k,irho) - dxinv(3) * &
                  first_deriv_8( cons(i,j,k-4:k+4,imz) )
          end do
       end do
    end do

    do k=dlo(3),dhi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             tmpz(i,j,k) = cons(i,j,k,imx) * q(i,j,k,qw)
          end do
       end do
    end do
    do k=slo(3),shi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imx)=rhs(i,j,k,imx) - dxinv(3) * first_deriv_8(tmpz(i,j,k-4:k+4))
          end do
       end do
    end do

    do k=dlo(3),dhi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             tmpz(i,j,k) = cons(i,j,k,imy) * q(i,j,k,qw)
          end do
       end do
    end do
    do k=slo(3),shi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imy)=rhs(i,j,k,imy) - dxinv(3) * first_deriv_8(tmpz(i,j,k-4:k+4))
          end do
       end do
    end do

    do k=dlo(3),dhi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             tmpz(i,j,k) = cons(i,j,k,imz)*q(i,j,k,qw) + q(i,j,k,qpres)
          end do
       end do
    end do
    do k=slo(3),shi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,imz)=rhs(i,j,k,imz) - dxinv(3) * first_deriv_8(tmpz(i,j,k-4:k+4))
          end do
       end do
    end do

    do k=dlo(3),dhi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             tmpz(i,j,k) = (cons(i,j,k,iene)+q(i,j,k,qpres))*q(i,j,k,qw)
          end do
       end do
    end do
    do k=slo(3),shi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhs(i,j,k,iene)=rhs(i,j,k,iene) - dxinv(3) * first_deriv_8(tmpz(i,j,k-4:k+4))
          end do
       end do
    end do

    do n = iry1, iry1+nspecies-1
       do k=dlo(3),dhi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                tmpz(i,j,k) = cons(i,j,k,n)*q(i,j,k,qw)
             end do
          end do
       end do
       do k=slo(3),shi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                rhs(i,j,k,n) = rhs(i,j,k,n) - dxinv(3) * first_deriv_8(tmpz(i,j,k-4:k+4))
             end do
          end do
       end do
    enddo

    ! ----------------- boundary -----------------------

    !
    ! ----- lo-x boundary -----
    !
    if (physbclo(1)) then 
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)

             ! if (bclo(1) .eq. WALL???) then
             !    i = lo(1)
             !    ! use completely right-biased stencil
             ! end if

             i = lo(1)+1
             ! use 3rd-order slightly right-biased stencil

             un(-1:2) = q(i-1:i+2,j,k,qu)

             rhs(i,j,k,irho) = rhs(i,j,k,irho) - dxinv(1) * &
                  first_deriv_r3( cons(i-1:i+2,j,k,imx) ) 

             rhs(i,j,k,imx) = rhs(i,j,k,imx) - dxinv(1) * &
                  first_deriv_r3( cons(i-1:i+2,j,k,imx)*un(-1:2)+q(i-1:i+2,j,k,qpres) )

             rhs(i,j,k,imy) = rhs(i,j,k,imy) - dxinv(1) * &
                  first_deriv_r3( cons(i-1:i+2,j,k,imy)*un(-1:2) ) 

             rhs(i,j,k,imz) = rhs(i,j,k,imz) - dxinv(1) * &
                  first_deriv_r3( cons(i-1:i+2,j,k,imz)*un(-1:2) ) 

             rhs(i,j,k,iene) = rhs(i,j,k,iene) - dxinv(1) * &
                  first_deriv_r3( (cons(i-1:i+2,j,k,iene)+q(i-1:i+2,j,k,qpres))*un(-1:2) )

             i = lo(1)+2
             ! use 4th-order stencil

             un(-2:2) = q(i-2:i+2,j,k,qu)

             rhs(i,j,k,irho) = rhs(i,j,k,irho) - dxinv(1) * &
                  first_deriv_4( cons(i-2:i+2,j,k,imx) ) 

             rhs(i,j,k,imx) = rhs(i,j,k,imx) - dxinv(1) * &
                  first_deriv_4( cons(i-2:i+2,j,k,imx)*un(-2:2)+q(i-2:i+2,j,k,qpres) )

             rhs(i,j,k,imy) = rhs(i,j,k,imy) - dxinv(1) * &
                  first_deriv_4( cons(i-2:i+2,j,k,imy)*un(-2:2) ) 

             rhs(i,j,k,imz) = rhs(i,j,k,imz) - dxinv(1) * &
                  first_deriv_4( cons(i-2:i+2,j,k,imz)*un(-2:2) ) 

             rhs(i,j,k,iene) = rhs(i,j,k,iene) - dxinv(1) * &
                  first_deriv_4( (cons(i-2:i+2,j,k,iene)+q(i-2:i+2,j,k,qpres))*un(-2:2) )


             i = lo(1)+3
             ! use 6th-order stencil

             un(-3:3) = q(i-3:i+3,j,k,qu)

             rhs(i,j,k,irho) = rhs(i,j,k,irho) - dxinv(1) * &
                  first_deriv_6( cons(i-3:i+3,j,k,imx) ) 

             rhs(i,j,k,imx) = rhs(i,j,k,imx) - dxinv(1) * &
                  first_deriv_6( cons(i-3:i+3,j,k,imx)*un(-3:3)+q(i-3:i+3,j,k,qpres) )

             rhs(i,j,k,imy) = rhs(i,j,k,imy) - dxinv(1) * &
                  first_deriv_6( cons(i-3:i+3,j,k,imy)*un(-3:3) ) 

             rhs(i,j,k,imz) = rhs(i,j,k,imz) - dxinv(1) * &
                  first_deriv_6( cons(i-3:i+3,j,k,imz)*un(-3:3) ) 

             rhs(i,j,k,iene) = rhs(i,j,k,iene) - dxinv(1) * &
                  first_deriv_6( (cons(i-3:i+3,j,k,iene)+q(i-3:i+3,j,k,qpres))*un(-3:3) )

          enddo
       enddo

       do n = iry1, iry1+nspecies-1
          do k=lo(3),hi(3)
             do j=lo(2),hi(2)

                ! if (bclo(1) .eq. WALL???) then
                !    i = lo(1)
                !    ! use completely right-biased stencil
                ! end if
                
                i = lo(1)+1
                ! use 3rd-order slightly right-biased stencil
                rhs(i,j,k,n) = rhs(i,j,k,n) - dxinv(1) * &
                     first_deriv_r3( cons(i-1:i+2,j,k,n)*q(i-1:i+2,j,k,qu) )

                i = lo(1)+2
                ! use 4th-order stencil
                rhs(i,j,k,n) = rhs(i,j,k,n) - dxinv(1) * &
                     first_deriv_4( cons(i-2:i+2,j,k,n)*q(i-2:i+2,j,k,qu) )

                i = lo(1)+3
                ! use 6th-order stencil
                rhs(i,j,k,n) = rhs(i,j,k,n) - dxinv(1) * &
                     first_deriv_6( cons(i-3:i+3,j,k,n)*q(i-3:i+3,j,k,qu) )
             end do
          end do
       end do
    end if

    !
    ! ----- hi-x boundary -----
    !
    if (physbchi(1)) then 
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             
             i = hi(1)-3
             ! use 6th-order stencil

             un(-3:3) = q(i-3:i+3,j,k,qu)

             rhs(i,j,k,irho) = rhs(i,j,k,irho) - dxinv(1) * &
                  first_deriv_6( cons(i-3:i+3,j,k,imx) ) 

             rhs(i,j,k,imx) = rhs(i,j,k,imx) - dxinv(1) * &
                  first_deriv_6( cons(i-3:i+3,j,k,imx)*un(-3:3)+q(i-3:i+3,j,k,qpres) )

             rhs(i,j,k,imy) = rhs(i,j,k,imy) - dxinv(1) * &
                  first_deriv_6( cons(i-3:i+3,j,k,imy)*un(-3:3) ) 

             rhs(i,j,k,imz) = rhs(i,j,k,imz) - dxinv(1) * &
                  first_deriv_6( cons(i-3:i+3,j,k,imz)*un(-3:3) ) 

             rhs(i,j,k,iene) = rhs(i,j,k,iene) - dxinv(1) * &
                  first_deriv_6( (cons(i-3:i+3,j,k,iene)+q(i-3:i+3,j,k,qpres))*un(-3:3) )

             i = hi(1)-2
             ! use 4th-order stencil

             un(-2:2) = q(i-2:i+2,j,k,qu)

             rhs(i,j,k,irho) = rhs(i,j,k,irho) - dxinv(1) * &
                  first_deriv_4( cons(i-2:i+2,j,k,imx) ) 

             rhs(i,j,k,imx) = rhs(i,j,k,imx) - dxinv(1) * &
                  first_deriv_4( cons(i-2:i+2,j,k,imx)*un(-2:2)+q(i-2:i+2,j,k,qpres) )

             rhs(i,j,k,imy) = rhs(i,j,k,imy) - dxinv(1) * &
                  first_deriv_4( cons(i-2:i+2,j,k,imy)*un(-2:2) ) 

             rhs(i,j,k,imz) = rhs(i,j,k,imz) - dxinv(1) * &
                  first_deriv_4( cons(i-2:i+2,j,k,imz)*un(-2:2) ) 

             rhs(i,j,k,iene) = rhs(i,j,k,iene) - dxinv(1) * &
                  first_deriv_4( (cons(i-2:i+2,j,k,iene)+q(i-2:i+2,j,k,qpres))*un(-2:2) )

             i = hi(1)-1
             ! use 3rd-order slightly left-biased stencil

             un(-2:1) = q(i-2:i+1,j,k,qu)

             rhs(i,j,k,irho) = rhs(i,j,k,irho) - dxinv(1) * &
                  first_deriv_l3( cons(i-2:i+1,j,k,imx) ) 

             rhs(i,j,k,imx) = rhs(i,j,k,imx) - dxinv(1) * &
                  first_deriv_l3( cons(i-2:i+1,j,k,imx)*un(-2:1)+q(i-2:i+1,j,k,qpres) )

             rhs(i,j,k,imy) = rhs(i,j,k,imy) - dxinv(1) * &
                  first_deriv_l3( cons(i-2:i+1,j,k,imy)*un(-2:1) ) 

             rhs(i,j,k,imz) = rhs(i,j,k,imz) - dxinv(1) * &
                  first_deriv_l3( cons(i-2:i+1,j,k,imz)*un(-2:1) ) 

             rhs(i,j,k,iene) = rhs(i,j,k,iene) - dxinv(1) * &
                  first_deriv_l3( (cons(i-2:i+1,j,k,iene)+q(i-2:i+1,j,k,qpres))*un(-2:1) )


             ! if (bchi(1) .eq. WALL???) then
             !    i = hi(1)
             !    ! use completely left-biased stencil
             ! end if

          enddo
       enddo
       
       do n = iry1, iry1+nspecies-1
          do k=lo(3),hi(3)
             do j=lo(2),hi(2)
                
                i = hi(1)-3
                ! use 6th-order stencil
                rhs(i,j,k,n) = rhs(i,j,k,n) - dxinv(1) * &
                     first_deriv_6( cons(i-3:i+3,j,k,n)*q(i-3:i+3,j,k,qu) )
                
                i = hi(1)-2
                ! use 4th-order stencil
                rhs(i,j,k,n) = rhs(i,j,k,n) - dxinv(1) * &
                     first_deriv_4( cons(i-2:i+2,j,k,n)*q(i-2:i+2,j,k,qu) )
                
                i = hi(1)-1
                ! use 3rd-order slightly left-biased stencil
                rhs(i,j,k,n) = rhs(i,j,k,n) - dxinv(1) * &
                     first_deriv_l3( cons(i-2:i+1,j,k,n)*q(i-2:i+1,j,k,qu) )
                
                ! if (bchi(1) .eq. WALL???) then
                !    i = hi(1)
                !    ! use completely left-biased stencil
                ! end if
                
             enddo
          enddo
       end do
    end if

    !
    ! ----- lo-y boundary -----
    !
    if (physbclo(2)) then 
       do k=lo(3),hi(3)

          ! if (bclo(2) .eq. WALL???) then
          !    j = lo(2)
          !    ! use completely right-biased stencil
          !    enddo
          ! end if

          j = lo(2)+1
          ! use 3rd-order slightly right-biased stencil

          do i=lo(1),hi(1)

             un(-1:2) = q(i,j-1:j+2,k,qv)

             rhs(i,j,k,irho)=rhs(i,j,k,irho) - dxinv(2) * &
                  first_deriv_r3( cons(i,j-1:j+2,k,imy) )

             rhs(i,j,k,imx)=rhs(i,j,k,imx) - dxinv(2) * &
                  first_deriv_r3( cons(i,j-1:j+2,k,imx)*un(-1:2) )

             rhs(i,j,k,imy)=rhs(i,j,k,imy) - dxinv(2) * &
                  first_deriv_r3( cons(i,j-1:j+2,k,imy)*un(-1:2)+q(i,j-1:j+2,k,qpres) )

             rhs(i,j,k,imz)=rhs(i,j,k,imz) - dxinv(2) * &
                  first_deriv_r3( cons(i,j-1:j+2,k,imz)*un(-1:2) )

             rhs(i,j,k,iene)=rhs(i,j,k,iene) - dxinv(2) * &
                  first_deriv_r3( (cons(i,j-1:j+2,k,iene)+q(i,j-1:j+2,k,qpres))*un(-1:2) )

          enddo

          j = lo(2)+2
          ! use 4th-order stencil

          do i=lo(1),hi(1)

             un(-2:2) = q(i,j-2:j+2,k,qv)

             rhs(i,j,k,irho)=rhs(i,j,k,irho) - dxinv(2) * &
                  first_deriv_4( cons(i,j-2:j+2,k,imy) )

             rhs(i,j,k,imx)=rhs(i,j,k,imx) - dxinv(2) * &
                  first_deriv_4( cons(i,j-2:j+2,k,imx)*un(-2:2) )

             rhs(i,j,k,imy)=rhs(i,j,k,imy) - dxinv(2) * &
                  first_deriv_4( cons(i,j-2:j+2,k,imy)*un(-2:2)+q(i,j-2:j+2,k,qpres) )

             rhs(i,j,k,imz)=rhs(i,j,k,imz) - dxinv(2) * &
                  first_deriv_4( cons(i,j-2:j+2,k,imz)*un(-2:2) )

             rhs(i,j,k,iene)=rhs(i,j,k,iene) - dxinv(2) * &
                  first_deriv_4( (cons(i,j-2:j+2,k,iene)+q(i,j-2:j+2,k,qpres))*un(-2:2) )

          enddo

          j = lo(2)+3
          ! use 6th-order stencil

          do i=lo(1),hi(1)

             un(-3:3) = q(i,j-3:j+3,k,qv)

             rhs(i,j,k,irho)=rhs(i,j,k,irho) - dxinv(2) * &
                  first_deriv_6( cons(i,j-3:j+3,k,imy) )

             rhs(i,j,k,imx)=rhs(i,j,k,imx) - dxinv(2) * &
                  first_deriv_6( cons(i,j-3:j+3,k,imx)*un(-3:3) )

             rhs(i,j,k,imy)=rhs(i,j,k,imy) - dxinv(2) * &
                  first_deriv_6( cons(i,j-3:j+3,k,imy)*un(-3:3)+q(i,j-3:j+3,k,qpres) )

             rhs(i,j,k,imz)=rhs(i,j,k,imz) - dxinv(2) * &
                  first_deriv_6( cons(i,j-3:j+3,k,imz)*un(-3:3) )

             rhs(i,j,k,iene)=rhs(i,j,k,iene) - dxinv(2) * &
                  first_deriv_6( (cons(i,j-3:j+3,k,iene)+q(i,j-3:j+3,k,qpres))*un(-3:3) )

          enddo

       enddo

       do n = iry1, iry1+nspecies-1
          do k=lo(3),hi(3)

             ! if (bclo(2) .eq. WALL???) then
             !    j = lo(2)
             !    ! use completely right-biased stencil
             !    enddo
             ! end if
             
             j = lo(2)+1
             ! use 3rd-order slightly right-biased stencil
             do i=lo(1),hi(1)
                rhs(i,j,k,n) = rhs(i,j,k,n) - dxinv(2) * &
                     first_deriv_r3( cons(i,j-1:j+2,k,n)*q(i,j-1:j+2,k,qv) )
             end do

             j = lo(2)+2
             ! use 4th-order stencil
             do i=lo(1),hi(1)
                rhs(i,j,k,n) = rhs(i,j,k,n) - dxinv(2) * &
                     first_deriv_4( cons(i,j-2:j+2,k,n)*q(i,j-2:j+2,k,qv) )
             end do

             j = lo(2)+3
             ! use 6th-order stencil
             do i=lo(1),hi(1)
                rhs(i,j,k,n) = rhs(i,j,k,n) - dxinv(2) * &
                     first_deriv_6( cons(i,j-3:j+3,k,n)*q(i,j-3:j+3,k,qv) )
             end do

          end do
       end do
    end if

    !
    ! ----- hi-y boundary -----
    !
    if (physbchi(2)) then 
       do k=lo(3),hi(3)

          j = hi(2)-3
          ! use 6th-order stencil

          do i=lo(1),hi(1)

             un(-3:3) = q(i,j-3:j+3,k,qv)

             rhs(i,j,k,irho)=rhs(i,j,k,irho) - dxinv(2) * &
                  first_deriv_6( cons(i,j-3:j+3,k,imy) )

             rhs(i,j,k,imx)=rhs(i,j,k,imx) - dxinv(2) * &
                  first_deriv_6( cons(i,j-3:j+3,k,imx)*un(-3:3) )

             rhs(i,j,k,imy)=rhs(i,j,k,imy) - dxinv(2) * &
                  first_deriv_6( cons(i,j-3:j+3,k,imy)*un(-3:3)+q(i,j-3:j+3,k,qpres) )

             rhs(i,j,k,imz)=rhs(i,j,k,imz) - dxinv(2) * &
                  first_deriv_6( cons(i,j-3:j+3,k,imz)*un(-3:3) )

             rhs(i,j,k,iene)=rhs(i,j,k,iene) - dxinv(2) * &
                  first_deriv_6( (cons(i,j-3:j+3,k,iene)+q(i,j-3:j+3,k,qpres))*un(-3:3) )

          enddo

          j = hi(2)-2
          ! use 4th-order stencil

          do i=lo(1),hi(1)

             un(-2:2) = q(i,j-2:j+2,k,qv)

             rhs(i,j,k,irho)=rhs(i,j,k,irho) - dxinv(2) * &
                  first_deriv_4( cons(i,j-2:j+2,k,imy) )

             rhs(i,j,k,imx)=rhs(i,j,k,imx) - dxinv(2) * &
                  first_deriv_4( cons(i,j-2:j+2,k,imx)*un(-2:2) )

             rhs(i,j,k,imy)=rhs(i,j,k,imy) - dxinv(2) * &
                  first_deriv_4( cons(i,j-2:j+2,k,imy)*un(-2:2)+q(i,j-2:j+2,k,qpres) )

             rhs(i,j,k,imz)=rhs(i,j,k,imz) - dxinv(2) * &
                  first_deriv_4( cons(i,j-2:j+2,k,imz)*un(-2:2) )

             rhs(i,j,k,iene)=rhs(i,j,k,iene) - dxinv(2) * &
                  first_deriv_4( (cons(i,j-2:j+2,k,iene)+q(i,j-2:j+2,k,qpres))*un(-2:2) )

          enddo

          j = hi(2)-1
          ! use 3rd-order slightly left-biased stencil

          do i=lo(1),hi(1)

             un(-2:1) = q(i,j-2:j+1,k,qv)

             rhs(i,j,k,irho)=rhs(i,j,k,irho) - dxinv(2) * &
                  first_deriv_l3( cons(i,j-2:j+1,k,imy) )

             rhs(i,j,k,imx)=rhs(i,j,k,imx) - dxinv(2) * &
                  first_deriv_l3( cons(i,j-2:j+1,k,imx)*un(-2:1) )

             rhs(i,j,k,imy)=rhs(i,j,k,imy) - dxinv(2) * &
                  first_deriv_l3( cons(i,j-2:j+1,k,imy)*un(-2:1)+q(i,j-2:j+1,k,qpres) )

             rhs(i,j,k,imz)=rhs(i,j,k,imz) - dxinv(2) * &
                  first_deriv_l3( cons(i,j-2:j+1,k,imz)*un(-2:1) )

             rhs(i,j,k,iene)=rhs(i,j,k,iene) - dxinv(2) * &
                  first_deriv_l3( (cons(i,j-2:j+1,k,iene)+q(i,j-2:j+1,k,qpres))*un(-2:1) )

          enddo

          ! if (bchi(2) .eq. WALL???) then
          !    j = hi(2)
          !    ! use completely left-biased stencil
          ! end if

       enddo

       do n = iry1, iry1+nspecies-1
          do k=lo(3),hi(3)

             j = hi(2)-3
             ! use 6th-order stencil
             do i=lo(1),hi(1)
                rhs(i,j,k,n) = rhs(i,j,k,n) - dxinv(2) * &
                     first_deriv_6( cons(i,j-3:j+3,k,n)*q(i,j-3:j+3,k,qv) )
             end do

             j = hi(2)-2
             ! use 4th-order stencil
             do i=lo(1),hi(1)
                rhs(i,j,k,n) = rhs(i,j,k,n) - dxinv(2) * &
                     first_deriv_4( cons(i,j-2:j+2,k,n)*q(i,j-2:j+2,k,qv) )
             end do

             j = hi(2)-1
             ! use 3rd-order slightly left-biased stencil
             do i=lo(1),hi(1)
                rhs(i,j,k,n) = rhs(i,j,k,n) - dxinv(2) * &
                     first_deriv_l3( cons(i,j-2:j+1,k,n)*q(i,j-2:j+1,k,qv) )
             end do

             ! if (bchi(2) .eq. WALL???) then
             !    j = hi(2)
             !    ! use completely left-biased stencil
             ! end if

          end do
       end do
    end if

    !
    ! ----- lo-z boundary -----
    !
    if (physbclo(3)) then

       ! if (bclo(3) .eq. WALL???) then
       !    k = lo(3)
       !    ! use completely right-biased stencil
       ! end if

       k = lo(3)+1
       ! use 3rd-order slightly right-biased stencil
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             un(-1:2) = q(i,j,k-1:k+2,qw)

             rhs(i,j,k,irho)=rhs(i,j,k,irho) - dxinv(3) * &
                  first_deriv_r3( cons(i,j,k-1:k+2,imz) )

             rhs(i,j,k,imx)=rhs(i,j,k,imx) - dxinv(3) * &
                  first_deriv_r3( cons(i,j,k-1:k+2,imx)*un(-1:2) )

             rhs(i,j,k,imy)=rhs(i,j,k,imy) - dxinv(3) * &
                  first_deriv_r3( cons(i,j,k-1:k+2,imy)*un(-1:2) )

             rhs(i,j,k,imz)=rhs(i,j,k,imz) - dxinv(3) * &
                  first_deriv_r3( cons(i,j,k-1:k+2,imz)*un(-1:2)+q(i,j,k-1:k+2,qpres) )

             rhs(i,j,k,iene)=rhs(i,j,k,iene) - dxinv(3) * &
                  first_deriv_r3( (cons(i,j,k-1:k+2,iene)+q(i,j,k-1:k+2,qpres))*un(-1:2) )

          enddo
       enddo

       do n = iry1, iry1+nspecies-1
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                rhs(i,j,k,n) = rhs(i,j,k,n) - dxinv(3) * &
                     first_deriv_r3( cons(i,j,k-1:k+2,n)*q(i,j,k-1:k+2,qw) )
             end do
          end do
       end do

       k = lo(3)+2
       ! use 4th-order stencil
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             un(-2:2) = q(i,j,k-2:k+2,qw)

             rhs(i,j,k,irho)=rhs(i,j,k,irho) - dxinv(3) * &
                  first_deriv_4( cons(i,j,k-2:k+2,imz) )

             rhs(i,j,k,imx)=rhs(i,j,k,imx) - dxinv(3) * &
                  first_deriv_4( cons(i,j,k-2:k+2,imx)*un(-2:2) )

             rhs(i,j,k,imy)=rhs(i,j,k,imy) - dxinv(3) * &
                  first_deriv_4( cons(i,j,k-2:k+2,imy)*un(-2:2) )

             rhs(i,j,k,imz)=rhs(i,j,k,imz) - dxinv(3) * &
                  first_deriv_4( cons(i,j,k-2:k+2,imz)*un(-2:2)+q(i,j,k-2:k+2,qpres) )

             rhs(i,j,k,iene)=rhs(i,j,k,iene) - dxinv(3) * &
                  first_deriv_4( (cons(i,j,k-2:k+2,iene)+q(i,j,k-2:k+2,qpres))*un(-2:2) )

          enddo
       enddo

       do n = iry1, iry1+nspecies-1
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                rhs(i,j,k,n) = rhs(i,j,k,n) - dxinv(3) * &
                     first_deriv_4( cons(i,j,k-2:k+2,n)*q(i,j,k-2:k+2,qw) )
             end do
          end do
       end do

       k = lo(3)+3
       ! use 6th-order stencil
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             un(-3:3) = q(i,j,k-3:k+3,qw)

             rhs(i,j,k,irho)=rhs(i,j,k,irho) - dxinv(3) * &
                  first_deriv_6( cons(i,j,k-3:k+3,imz) )

             rhs(i,j,k,imx)=rhs(i,j,k,imx) - dxinv(3) * &
                  first_deriv_6( cons(i,j,k-3:k+3,imx)*un(-3:3) )

             rhs(i,j,k,imy)=rhs(i,j,k,imy) - dxinv(3) * &
                  first_deriv_6( cons(i,j,k-3:k+3,imy)*un(-3:3) )

             rhs(i,j,k,imz)=rhs(i,j,k,imz) - dxinv(3) * &
                  first_deriv_6( cons(i,j,k-3:k+3,imz)*un(-3:3)+q(i,j,k-3:k+3,qpres) )

             rhs(i,j,k,iene)=rhs(i,j,k,iene) - dxinv(3) * &
                  first_deriv_6( (cons(i,j,k-3:k+3,iene)+q(i,j,k-3:k+3,qpres))*un(-3:3) )

          enddo
       enddo

       do n = iry1, iry1+nspecies-1
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                rhs(i,j,k,n) = rhs(i,j,k,n) - dxinv(3) * &
                     first_deriv_6( cons(i,j,k-3:k+3,n)*q(i,j,k-3:k+3,qw) )
             end do
          end do
       end do
    end if

    !
    ! ----- hi-z boundary -----
    !
    if (physbchi(3)) then

       k = hi(3)-3
       ! use 6th-order stencil
       
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             un(-3:3) = q(i,j,k-3:k+3,qw)

             rhs(i,j,k,irho)=rhs(i,j,k,irho) - dxinv(3) * &
                  first_deriv_6( cons(i,j,k-3:k+3,imz) )

             rhs(i,j,k,imx)=rhs(i,j,k,imx) - dxinv(3) * &
                  first_deriv_6( cons(i,j,k-3:k+3,imx)*un(-3:3) )

             rhs(i,j,k,imy)=rhs(i,j,k,imy) - dxinv(3) * &
                  first_deriv_6( cons(i,j,k-3:k+3,imy)*un(-3:3) )

             rhs(i,j,k,imz)=rhs(i,j,k,imz) - dxinv(3) * &
                  first_deriv_6( cons(i,j,k-3:k+3,imz)*un(-3:3)+q(i,j,k-3:k+3,qpres) )

             rhs(i,j,k,iene)=rhs(i,j,k,iene) - dxinv(3) * &
                  first_deriv_6( (cons(i,j,k-3:k+3,iene)+q(i,j,k-3:k+3,qpres))*un(-3:3) )

          enddo
       enddo

       do n = iry1, iry1+nspecies-1
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                rhs(i,j,k,n) = rhs(i,j,k,n) - dxinv(3) * &
                     first_deriv_6(cons(i,j,k-3:k+3,n)*q(i,j,k-3:k+3,qw))
             end do
          enddo
       end do

       k = hi(3)-2
       ! use 4th-order stencil
       
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             un(-2:2) = q(i,j,k-2:k+2,qw)

             rhs(i,j,k,irho)=rhs(i,j,k,irho) - dxinv(3) * &
                  first_deriv_4( cons(i,j,k-2:k+2,imz) )

             rhs(i,j,k,imx)=rhs(i,j,k,imx) - dxinv(3) * &
                  first_deriv_4( cons(i,j,k-2:k+2,imx)*un(-2:2) )

             rhs(i,j,k,imy)=rhs(i,j,k,imy) - dxinv(3) * &
                  first_deriv_4( cons(i,j,k-2:k+2,imy)*un(-2:2) )

             rhs(i,j,k,imz)=rhs(i,j,k,imz) - dxinv(3) * &
                  first_deriv_4( cons(i,j,k-2:k+2,imz)*un(-2:2)+q(i,j,k-2:k+2,qpres) )

             rhs(i,j,k,iene)=rhs(i,j,k,iene) - dxinv(3) * &
                  first_deriv_4( (cons(i,j,k-2:k+2,iene)+q(i,j,k-2:k+2,qpres))*un(-2:2) )

          enddo
       enddo

       do n = iry1, iry1+nspecies-1
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                rhs(i,j,k,n) = rhs(i,j,k,n) - dxinv(3) * &
                     first_deriv_4(cons(i,j,k-2:k+2,n)*q(i,j,k-2:k+2,qw))
             end do
          enddo
       enddo

       k = hi(3)-1
       ! use 3rd-order slightly left-biased stencil
       
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             un(-2:1) = q(i,j,k-2:k+1,qw)

             rhs(i,j,k,irho)=rhs(i,j,k,irho) - dxinv(3) * &
                  first_deriv_l3( cons(i,j,k-2:k+1,imz) )

             rhs(i,j,k,imx)=rhs(i,j,k,imx) - dxinv(3) * &
                  first_deriv_l3( cons(i,j,k-2:k+1,imx)*un(-2:1) )

             rhs(i,j,k,imy)=rhs(i,j,k,imy) - dxinv(3) * &
                  first_deriv_l3( cons(i,j,k-2:k+1,imy)*un(-2:1) )

             rhs(i,j,k,imz)=rhs(i,j,k,imz) - dxinv(3) * &
                  first_deriv_l3( cons(i,j,k-2:k+1,imz)*un(-2:1)+q(i,j,k-2:k+1,qpres) )

             rhs(i,j,k,iene)=rhs(i,j,k,iene) - dxinv(3) * &
                  first_deriv_l3( (cons(i,j,k-2:k+1,iene)+q(i,j,k-2:k+1,qpres))*un(-2:1) )

          enddo
       enddo

       do n = iry1, iry1+nspecies-1
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                rhs(i,j,k,n) = rhs(i,j,k,n) - dxinv(3) * &
                     first_deriv_l3(cons(i,j,k-2:k+1,n)*q(i,j,k-2:k+1,qw))
             end do
          enddo
       enddo

       ! if (bchi(3) .eq. WALL???) then
       !    k = hi(3)
       !    ! use completely left-biased stencil
       ! end if

    end if

    deallocate(tmpx,tmpy,tmpz)

    rhs_g(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:) = &
         rhs_g(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:) &
         + rhs(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:)
    deallocate(rhs)

  end subroutine hypterm_3d


  subroutine narrow_diffterm_3d (lo,hi,dx,q,qlo,qhi,rhs_g,glo,ghi,rhs,rlo,rhi, &
       mu,xi,lam,dxy,dlo_g,dhi_g,bclo,bchi)

    integer,         intent(in):: dlo_g(3),dhi_g(3),bclo(3),bchi(3)
    integer,         intent(in):: lo(3),hi(3),qlo(3),qhi(3),rlo(3),rhi(3),glo(3),ghi(3)
    double precision,intent(in):: dx(3)
    double precision,intent(in)   ::  q  (qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3),nprim)
    double precision,intent(in)   ::  mu (qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3))
    double precision,intent(in)   ::  xi (qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3))
    double precision,intent(in)   ::  lam(qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3))
    double precision,intent(in)   ::  dxy(qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3),nspecies)
    double precision,intent(inout)::rhs_g(glo(1):ghi(1),glo(2):ghi(2),glo(3):ghi(3),ncons)
    double precision,intent(inout)::rhs  (rlo(1):rhi(1),rlo(2):rhi(2),rlo(3):rhi(3),ncons)

    integer :: i
    integer :: slo(3), shi(3), dlo(3), dhi(3)
    double precision :: dxinv(3), dx2inv(3)

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
       dx2inv(i) = dxinv(i)**2
    end do

    rhs(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:) = 0.d0

    call diffterm_1(q,qlo,qhi,rhs,rlo,rhi,mu,xi, &
         lo,hi,slo,shi,dlo,dhi,finlo,finhi,foulo,fouhi,physbclo,physbchi,dxinv)

    call diffterm_2(q,qlo,qhi,rhs,rlo,rhi, mu,xi,lam,dxy, &
         lo,hi,slo,shi,dlo,dhi,finlo,finhi,foulo,fouhi,physbclo,physbchi,dxinv,dx2inv)

    rhs_g     (lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:) = &
         rhs_g(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:) &
         + rhs(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:)

  end subroutine narrow_diffterm_3d

  
  subroutine diffterm_1(q,qlo,qhi,rhs,rlo,rhi,mu,xi, &
       lo,hi,slo,shi,dlo,dhi,finlo,finhi,foulo,fouhi,physbclo,physbchi,dxinv)
    integer,         intent(in):: lo(3),hi(3),slo(3),shi(3),dlo(3),dhi(3)
    integer,         intent(in):: qlo(3),qhi(3),rlo(3),rhi(3)
    logical,         intent(in):: physbclo(3),physbchi(3)
    double precision,intent(in):: finlo(3),finhi(3),foulo(3),fouhi(3)
    double precision,intent(in):: dxinv(3)
    double precision,intent(in)   ::  q(qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3),nprim)
    double precision,intent(in)   :: mu(qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3))
    double precision,intent(in)   :: xi(qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3))
    double precision,intent(inout)::rhs(rlo(1):rhi(1),rlo(2):rhi(2),rlo(3):rhi(3),ncons)

    double precision, allocatable, dimension(:,:,:) :: ux,uy,uz,vx,vy,vz,wx,wy,wz
    double precision, allocatable :: tmpx(:), tmpy(:,:),tmpz(:,:,:)
    double precision, allocatable, dimension(:,:,:) :: vsm
    double precision, dimension(lo(1):hi(1)) :: tauxx,tauyy,tauzz,divu
    integer          :: i,j,k

    allocate(ux( lo(1): hi(1),dlo(2):dhi(2),dlo(3):dhi(3)))
    allocate(vx( lo(1): hi(1),dlo(2):dhi(2),dlo(3):dhi(3)))
    allocate(wx( lo(1): hi(1),dlo(2):dhi(2),dlo(3):dhi(3)))

    allocate(uy(dlo(1):dhi(1), lo(2): hi(2),dlo(3):dhi(3)))
    allocate(vy(dlo(1):dhi(1), lo(2): hi(2),dlo(3):dhi(3)))
    allocate(wy(dlo(1):dhi(1), lo(2): hi(2),dlo(3):dhi(3)))

    allocate(uz(dlo(1):dhi(1),dlo(2):dhi(2), lo(3): hi(3)))
    allocate(vz(dlo(1):dhi(1),dlo(2):dhi(2), lo(3): hi(3)))
    allocate(wz(dlo(1):dhi(1),dlo(2):dhi(2), lo(3): hi(3)))

    allocate(vsm(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3)))

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
             ux(i,j,k) = dxinv(1)*first_deriv_8(q(i-4:i+4,j,k,qu))
             vx(i,j,k) = dxinv(1)*first_deriv_8(q(i-4:i+4,j,k,qv))
             wx(i,j,k) = dxinv(1)*first_deriv_8(q(i-4:i+4,j,k,qw))
          enddo
       enddo
    enddo

    do k=dlo(3),dhi(3)
       do j=slo(2),shi(2)   
          do i=dlo(1),dhi(1)
             uy(i,j,k) = dxinv(2)*first_deriv_8(q(i,j-4:j+4,k,qu))
             vy(i,j,k) = dxinv(2)*first_deriv_8(q(i,j-4:j+4,k,qv))
             wy(i,j,k) = dxinv(2)*first_deriv_8(q(i,j-4:j+4,k,qw))
          enddo
       enddo
    enddo

    do k=slo(3),shi(3)
       do j=dlo(2),dhi(2)
          do i=dlo(1),dhi(1)
             uz(i,j,k) = dxinv(3)*first_deriv_8(q(i,j,k-4:k+4,qu))
             vz(i,j,k) = dxinv(3)*first_deriv_8(q(i,j,k-4:k+4,qv))
             wz(i,j,k) = dxinv(3)*first_deriv_8(q(i,j,k-4:k+4,qw))
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
             ux(i,j,k) = dxinv(1)*first_deriv_rb(q(i:i+3,j,k,qu))
             vx(i,j,k) = dxinv(1)*first_deriv_rb(q(i:i+3,j,k,qv))
             wx(i,j,k) = dxinv(1)*first_deriv_rb(q(i:i+3,j,k,qw))

             i = lo(1)+1
             ! use 3rd-order slightly right-biased stencil
             ux(i,j,k) = dxinv(1)*first_deriv_r3(q(i-1:i+2,j,k,qu))
             vx(i,j,k) = dxinv(1)*first_deriv_r3(q(i-1:i+2,j,k,qv))
             wx(i,j,k) = dxinv(1)*first_deriv_r3(q(i-1:i+2,j,k,qw))

             i = lo(1)+2
             ! use 4th-order stencil
             ux(i,j,k) = dxinv(1)*first_deriv_4(q(i-2:i+2,j,k,qu))
             vx(i,j,k) = dxinv(1)*first_deriv_4(q(i-2:i+2,j,k,qv))
             wx(i,j,k) = dxinv(1)*first_deriv_4(q(i-2:i+2,j,k,qw))

             i = lo(1)+3
             ! use 6th-order stencil
             ux(i,j,k) = dxinv(1)*first_deriv_6(q(i-3:i+3,j,k,qu))
             vx(i,j,k) = dxinv(1)*first_deriv_6(q(i-3:i+3,j,k,qv))
             wx(i,j,k) = dxinv(1)*first_deriv_6(q(i-3:i+3,j,k,qw))
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
             ux(i,j,k) = dxinv(1)*first_deriv_6(q(i-3:i+3,j,k,qu))
             vx(i,j,k) = dxinv(1)*first_deriv_6(q(i-3:i+3,j,k,qv))
             wx(i,j,k) = dxinv(1)*first_deriv_6(q(i-3:i+3,j,k,qw))

             i = hi(1)-2
             ! use 4th-order stencil
             ux(i,j,k) = dxinv(1)*first_deriv_4(q(i-2:i+2,j,k,qu))
             vx(i,j,k) = dxinv(1)*first_deriv_4(q(i-2:i+2,j,k,qv))
             wx(i,j,k) = dxinv(1)*first_deriv_4(q(i-2:i+2,j,k,qw))

             i = hi(1)-1
             ! use 3rd-order slightly left-biased stencil
             ux(i,j,k) = dxinv(1)*first_deriv_l3(q(i-2:i+1,j,k,qu))
             vx(i,j,k) = dxinv(1)*first_deriv_l3(q(i-2:i+1,j,k,qv))
             wx(i,j,k) = dxinv(1)*first_deriv_l3(q(i-2:i+1,j,k,qw))

             i = hi(1)
             ! use completely left-biased stencil
             ux(i,j,k) = dxinv(1)*first_deriv_lb(q(i-3:i,j,k,qu))
             vx(i,j,k) = dxinv(1)*first_deriv_lb(q(i-3:i,j,k,qv))
             wx(i,j,k) = dxinv(1)*first_deriv_lb(q(i-3:i,j,k,qw))
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
             uy(i,j,k) = dxinv(2)*first_deriv_rb(q(i,j:j+3,k,qu))
             vy(i,j,k) = dxinv(2)*first_deriv_rb(q(i,j:j+3,k,qv))
             wy(i,j,k) = dxinv(2)*first_deriv_rb(q(i,j:j+3,k,qw))
          enddo

          j = lo(2)+1
          ! use 3rd-order slightly right-biased stencil
          do i=dlo(1),dhi(1)
             uy(i,j,k) = dxinv(2)*first_deriv_r3(q(i,j-1:j+2,k,qu))
             vy(i,j,k) = dxinv(2)*first_deriv_r3(q(i,j-1:j+2,k,qv))
             wy(i,j,k) = dxinv(2)*first_deriv_r3(q(i,j-1:j+2,k,qw))
          enddo

          j = lo(2)+2
          ! use 4th-order stencil
          do i=dlo(1),dhi(1)
             uy(i,j,k) = dxinv(2)*first_deriv_4(q(i,j-2:j+2,k,qu))
             vy(i,j,k) = dxinv(2)*first_deriv_4(q(i,j-2:j+2,k,qv))
             wy(i,j,k) = dxinv(2)*first_deriv_4(q(i,j-2:j+2,k,qw))
          enddo

          j = lo(2)+3
          ! use 6th-order stencil
          do i=dlo(1),dhi(1)
             uy(i,j,k) = dxinv(2)*first_deriv_6(q(i,j-3:j+3,k,qu))
             vy(i,j,k) = dxinv(2)*first_deriv_6(q(i,j-3:j+3,k,qv))
             wy(i,j,k) = dxinv(2)*first_deriv_6(q(i,j-3:j+3,k,qw))
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
             uy(i,j,k) = dxinv(2)*first_deriv_6(q(i,j-3:j+3,k,qu))
             vy(i,j,k) = dxinv(2)*first_deriv_6(q(i,j-3:j+3,k,qv))
             wy(i,j,k) = dxinv(2)*first_deriv_6(q(i,j-3:j+3,k,qw))
          enddo

          j = hi(2)-2
          ! use 4th-order stencil
          do i=dlo(1),dhi(1)
             uy(i,j,k) = dxinv(2)*first_deriv_4(q(i,j-2:j+2,k,qu))
             vy(i,j,k) = dxinv(2)*first_deriv_4(q(i,j-2:j+2,k,qv))
             wy(i,j,k) = dxinv(2)*first_deriv_4(q(i,j-2:j+2,k,qw))
          enddo

          j = hi(2)-1
          ! use 3rd-order slightly left-biased stencil
          do i=dlo(1),dhi(1)
             uy(i,j,k) = dxinv(2)*first_deriv_l3(q(i,j-2:j+1,k,qu))
             vy(i,j,k) = dxinv(2)*first_deriv_l3(q(i,j-2:j+1,k,qv))
             wy(i,j,k) = dxinv(2)*first_deriv_l3(q(i,j-2:j+1,k,qw))
          enddo

          j = hi(2)
          ! use completely left-biased stencil
          do i=dlo(1),dhi(1)
             uy(i,j,k) = dxinv(2)*first_deriv_lb(q(i,j-3:j,k,qu))
             vy(i,j,k) = dxinv(2)*first_deriv_lb(q(i,j-3:j,k,qv))
             wy(i,j,k) = dxinv(2)*first_deriv_lb(q(i,j-3:j,k,qw))
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
             uz(i,j,k) = dxinv(3)*first_deriv_rb(q(i,j,k:k+3,qu))
             vz(i,j,k) = dxinv(3)*first_deriv_rb(q(i,j,k:k+3,qv))
             wz(i,j,k) = dxinv(3)*first_deriv_rb(q(i,j,k:k+3,qw))
          enddo
       enddo

       k = lo(3)+1
       ! use 3rd-order slightly right-biased stencil
       do j=dlo(2),dhi(2)
          do i=dlo(1),dhi(1)
             uz(i,j,k) = dxinv(3)*first_deriv_r3(q(i,j,k-1:k+2,qu))
             vz(i,j,k) = dxinv(3)*first_deriv_r3(q(i,j,k-1:k+2,qv))
             wz(i,j,k) = dxinv(3)*first_deriv_r3(q(i,j,k-1:k+2,qw))
          enddo
       enddo

       k = lo(3)+2
       ! use 4th-order stencil
       do j=dlo(2),dhi(2)
          do i=dlo(1),dhi(1)
             uz(i,j,k) = dxinv(3)*first_deriv_4(q(i,j,k-2:k+2,qu))
             vz(i,j,k) = dxinv(3)*first_deriv_4(q(i,j,k-2:k+2,qv))
             wz(i,j,k) = dxinv(3)*first_deriv_4(q(i,j,k-2:k+2,qw))
          enddo
       enddo

       k = lo(3)+3
       ! use 6th-order stencil
       do j=dlo(2),dhi(2)
          do i=dlo(1),dhi(1)
             uz(i,j,k) = dxinv(3)*first_deriv_6(q(i,j,k-3:k+3,qu))
             vz(i,j,k) = dxinv(3)*first_deriv_6(q(i,j,k-3:k+3,qv))
             wz(i,j,k) = dxinv(3)*first_deriv_6(q(i,j,k-3:k+3,qw))
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
             uz(i,j,k) = dxinv(3)*first_deriv_6(q(i,j,k-3:k+3,qu))
             vz(i,j,k) = dxinv(3)*first_deriv_6(q(i,j,k-3:k+3,qv))
             wz(i,j,k) = dxinv(3)*first_deriv_6(q(i,j,k-3:k+3,qw))
          enddo
       enddo

       k = hi(3)-2
       ! use 4th-order stencil
       do j=dlo(2),dhi(2)
          do i=dlo(1),dhi(1)
             uz(i,j,k) = dxinv(3)*first_deriv_4(q(i,j,k-2:k+2,qu))
             vz(i,j,k) = dxinv(3)*first_deriv_4(q(i,j,k-2:k+2,qv))
             wz(i,j,k) = dxinv(3)*first_deriv_4(q(i,j,k-2:k+2,qw))
          enddo
       enddo

       k = hi(3)-1
       ! use 3rd-order slightly left-biased stencil
       do j=dlo(2),dhi(2)
          do i=dlo(1),dhi(1)
             uz(i,j,k) = dxinv(3)*first_deriv_l3(q(i,j,k-2:k+1,qu))
             vz(i,j,k) = dxinv(3)*first_deriv_l3(q(i,j,k-2:k+1,qv))
             wz(i,j,k) = dxinv(3)*first_deriv_l3(q(i,j,k-2:k+1,qw))
          enddo
       enddo

       k = hi(3)
       ! use completely left-biased stencil
       do j=dlo(2),dhi(2)
          do i=dlo(1),dhi(1)
             uz(i,j,k) = dxinv(3)*first_deriv_lb(q(i,j,k-3:k,qu))
             vz(i,j,k) = dxinv(3)*first_deriv_lb(q(i,j,k-3:k,qv))
             wz(i,j,k) = dxinv(3)*first_deriv_lb(q(i,j,k-3:k,qw))
          enddo
       enddo
    end if


    !----- mx -----

    !----- mx : d()/dx -----
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)

          do i=dlo(1),dhi(1)
             tmpx(i) = vsm(i,j,k)*(vy(i,j,k)+wz(i,j,k))
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
             tmpy(i,j) = mu(i,j,k)*vx(i,j,k)
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
             tmpz(i,j,k) = mu(i,j,k)*wx(i,j,k)
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
             tmpx(i) = mu(i,j,k)*uy(i,j,k)
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
             tmpy(i,j) = vsm(i,j,k)*(ux(i,j,k)+wz(i,j,k))
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
             tmpz(i,j,k) = mu(i,j,k)*wy(i,j,k)
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
             tmpx(i) = mu(i,j,k)*uz(i,j,k)
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
             tmpy(i,j) = mu(i,j,k)*vz(i,j,k)
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
             tmpz(i,j,k) = vsm(i,j,k)*(ux(i,j,k)+vy(i,j,k))
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

             divu(i) = (ux(i,j,k)+vy(i,j,k)+wz(i,j,k))*vsm(i,j,k)
             tauxx(i) = 2.d0*mu(i,j,k)*ux(i,j,k) + divu(i)
             tauyy(i) = 2.d0*mu(i,j,k)*vy(i,j,k) + divu(i)
             tauzz(i) = 2.d0*mu(i,j,k)*wz(i,j,k) + divu(i)
             
             ! change in internal energy
             rhs(i,j,k,iene) = rhs(i,j,k,iene) + &
                  tauxx(i)*ux(i,j,k) + tauyy(i)*vy(i,j,k) + tauzz(i)*wz(i,j,k) &
                  + mu(i,j,k)*((uy(i,j,k)+vx(i,j,k))**2 &
                  &          + (wx(i,j,k)+uz(i,j,k))**2 &
                  &          + (vz(i,j,k)+wy(i,j,k))**2 )

          end do
       end do
    end do

    deallocate(tmpx,tmpy,tmpz)
    deallocate(vsm)
    deallocate(ux,uy,uz,vx,vy,vz,wx,wy,wz)

  end subroutine diffterm_1

  
  subroutine diffterm_2(q,qlo,qhi,rhs,rlo,rhi,mu,xi,lam,dxy, &
       lo,hi,slo,shi,dlo,dhi,finlo,finhi,foulo,fouhi,physbclo,physbchi,dxinv,dx2inv)
    use probin_module, only : reset_inactive_species
    integer,         intent(in):: lo(3),hi(3),slo(3),shi(3),dlo(3),dhi(3)
    integer,         intent(in):: qlo(3),qhi(3),rlo(3),rhi(3)
    logical,         intent(in):: physbclo(3),physbchi(3)
    double precision,intent(in):: finlo(3),finhi(3),foulo(3),fouhi(3)
    double precision,intent(in):: dxinv(3),dx2inv(3)
    double precision,intent(in)   :: q (qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3),nprim)
    double precision,intent(in)   :: mu(qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3))
    double precision,intent(in)   :: xi(qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3))
    double precision,intent(in)   ::lam(qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3))
    double precision,intent(in)   ::dxy(qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3),nspecies)
    double precision,intent(inout)::rhs(rlo(1):rhi(1),rlo(2):rhi(2),rlo(3):rhi(3),ncons)

    double precision, allocatable, dimension(:,:,:) :: vsp, dpe
    double precision, allocatable, dimension(:,:,:,:) :: Hg, dpy, dxe
    ! dxy: diffusion coefficient of X in equation for Y
    ! dpy: diffusion coefficient of p in equation for Y
    ! dxe: diffusion coefficient of X in equation for energy
    ! dpe: diffusion coefficient of p in equation for energy

    integer          :: i,j,k,n, qxn, qyn, qhn, iryn

    double precision :: mmtmp8(8,lo(1):hi(1)+1)
    double precision, allocatable, dimension(:,:,:,:) :: M8p
    double precision, allocatable, dimension(:,:,:) :: sumdrY, sumrYv, gradp
    double precision :: ry_c, ene_c

    double precision :: hhalf, sumdrytmp, sumryvtmp, gradptmp
    double precision :: Htot, Htmp(nspecies), Ytmp(nspecies)
    double precision :: M6p(6), M6X(6), mmtmp6(6)
    double precision :: M4p(4), M4X(4), mmtmp4(4)
    double precision :: M2p(2), M2X(2), mmtmp2(2)
    double precision :: BBp(4), BBX(4), mmtmpB(4)
    double precision :: rhstmp(nspecies), rhstot, rhsene
    double precision :: Hcell(0:1,2:ncons)
    integer :: iface

    logical :: add_v_correction
    add_v_correction = .not. reset_inactive_species

    allocate(vsp(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3)))

    allocate(dpy(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),nspecies))
    allocate(dxe(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),nspecies))
    allocate(dpe(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3)))

    allocate(Hg(lo(1):hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)+1,2:ncons))

    allocate(M8p(8,lo(1):hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)+1))

    allocate(sumdrY(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
    allocate(sumrYv(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
    allocate(gradp (lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))

    do k=dlo(3),dhi(3)
       do j=dlo(2),dhi(2)
          do i=dlo(1),dhi(1)
             vsp(i,j,k) = xi(i,j,k) + FourThirds*mu(i,j,k)
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
                dxe(i,j,k,n) = dxy(i,j,k,n)*q(i,j,k,qhn)
                dpe(i,j,k) = dpe(i,j,k) + dpy(i,j,k,n)*q(i,j,k,qhn)
             end do
          end do
       end do
    end do

    ! ------- BEGIN x-direction -------

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=slo(1),shi(1)+1
             mmtmp8(1:8,i) = matmul(vsp(i-4:i+3,j,k), M8)
             Hg(i,j,k,imx) = dot_product(mmtmp8(1:8,i), q(i-4:i+3,j,k,qu))
          end do
       end do
    end do

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=slo(1),shi(1)+1
             mmtmp8(1:8,i) = matmul(mu(i-4:i+3,j,k), M8)
             Hg(i,j,k,imy) = dot_product(mmtmp8(1:8,i), q(i-4:i+3,j,k,qv))
             Hg(i,j,k,imz) = dot_product(mmtmp8(1:8,i), q(i-4:i+3,j,k,qw))
          end do
       end do
    end do

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=slo(1),shi(1)+1
             mmtmp8(1:8,i) = matmul(lam(i-4:i+3,j,k), M8)
             Hg(i,j,k,iene) = dot_product(mmtmp8(1:8,i), q(i-4:i+3,j,k,qtemp))
          end do
       end do
    end do

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=slo(1),shi(1)+1
             mmtmp8(1:8,i) = matmul(M8T,  q(i-4:i+3,j,k,qpres))
             Hg(i,j,k,iene) = Hg(i,j,k,iene) + dot_product(dpe(i-4:i+3,j,k), mmtmp8(1:8,i))
          end do
          do i=slo(1),shi(1)+1
             M8p(:,i,j,k) = mmtmp8(1:8,i)
          end do
       end do
    end do

    do n = 1, nspecies

       if (n .eq. iias) cycle  ! inactive speices

       qxn = qx1+n-1
       iryn = iry1+n-1

       do k=lo(3),hi(3)
          do j=lo(2),hi(2)   
             do i=slo(1),shi(1)+1
                Hg(i,j,k,iryn) = dot_product(dpy(i-4:i+3,j,k,n), M8p(:,i,j,k))
             end do
          end do
       end do

       do k=lo(3),hi(3)
          do j=lo(2),hi(2)   
             do i=slo(1),shi(1)+1
                mmtmp8(1:8,i) = matmul(M8T, q(i-4:i+3,j,k,qxn))
                Hg(i,j,k,iene) = Hg(i,j,k,iene) + dot_product(dxe(i-4:i+3,j,k,n), mmtmp8(1:8,i))
                Hg(i,j,k,iryn) = Hg(i,j,k,iryn) &
                     + dot_product(dxy(i-4:i+3,j,k,n), mmtmp8(1:8,i))
             end do
          end do
       end do

    end do

    ! add x-direction rhs

    do n=2,iene
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=slo(1),shi(1)
                rhs(i,j,k,n) = rhs(i,j,k,n) + (Hg(i+1,j,k,n) - Hg(i,j,k,n)) * dx2inv(1)
             end do
          end do
       end do
    end do
       
    sumdrY = 0.d0
    do n=iry1,ncons

       if (n.eq.iry_ias) cycle

       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=slo(1),shi(1)
                sumdry(i,j,k) = sumdry(i,j,k) + (Hg(i+1,j,k,n) - Hg(i,j,k,n)) * dx2inv(1)
                rhs(i,j,k,n)  =  rhs(i,j,k,n) + (Hg(i+1,j,k,n) - Hg(i,j,k,n)) * dx2inv(1)
             end do
          end do
       end do

    end do

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=slo(1),shi(1)
             gradp(i,j,k) = dxinv(1) * first_deriv_8(q(i-4:i+4,j,k,qpres))
          end do
       end do
    end do
       
    sumryv = 0.d0
    do n = 1, nspecies

       if (n.eq.iias) cycle

       qxn = qx1+n-1
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=slo(1),shi(1)
                sumryv(i,j,k) = sumryv(i,j,k) + dpy(i,j,k,n)*gradp(i,j,k)  &
                     + dxy(i,j,k,n)*dxinv(1)*first_deriv_8(q(i-4:i+4,j,k,qxn))
             end do
          end do
       end do

    end do

    if (add_v_correction) then

       do n=1,nspecies
          qyn = qy1+n-1
          qhn = qh1+n-1
          iryn = iry1+n-1
          
          do k=lo(3),hi(3)
             do j=lo(2),hi(2)
                do i=slo(1),shi(1)
                   ry_c = q(i,j,k,qyn)*sumdry(i,j,k) + sumryv(i,j,k)*dxinv(1) * &
                        first_deriv_8(q(i-4:i+4,j,k,qyn))
                   ene_c = ry_c*q(i,j,k,qhn) + q(i,j,k,qyn)*sumryv(i,j,k)*dxinv(1)* &
                        first_deriv_8(q(i-4:i+4,j,k,qhn))
                   rhs(i,j,k,iene) = rhs(i,j,k,iene) - ene_c
                   rhs(i,j,k,iryn) = rhs(i,j,k,iryn) - ry_c
                end do
             end do
          end do
       end do

    else
    
       n = iias
       qhn = qh1+n-1
       iryn = iry1+n-1

       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=slo(1),shi(1)
                ene_c = sumdry(i,j,k)*q(i,j,k,qhn) + sumryv(i,j,k)*dxinv(1)* &
                     first_deriv_8(q(i-4:i+4,j,k,qhn))
                rhs(i,j,k,iene) = rhs(i,j,k,iene) - ene_c
                rhs(i,j,k,iryn) = rhs(i,j,k,iryn) - sumdry(i,j,k)
             end do
          end do
       end do
       
    end if

    ! ------- END x-direction -------

    ! ------- BEGIN y-direction -------

    do k=lo(3),hi(3)
       do j=slo(2),shi(2)+1
          do i=lo(1),hi(1)             
             mmtmp8(1:8,i) = matmul(mu(i,j-4:j+3,k), M8)
             Hg(i,j,k,imx) = dot_product(mmtmp8(1:8,i), q(i,j-4:j+3,k,qu))
             Hg(i,j,k,imz) = dot_product(mmtmp8(1:8,i), q(i,j-4:j+3,k,qw))
          end do
       end do
    end do

    do k=lo(3),hi(3)
       do j=slo(2),shi(2)+1
          do i=lo(1),hi(1)             
             mmtmp8(1:8,i) = matmul(vsp(i,j-4:j+3,k), M8)
             Hg(i,j,k,imy) = dot_product(mmtmp8(1:8,i), q(i,j-4:j+3,k,qv))
          end do
       end do
    end do

    do k=lo(3),hi(3)
       do j=slo(2),shi(2)+1
          do i=lo(1),hi(1)
             mmtmp8(1:8,i) = matmul(lam(i,j-4:j+3,k), M8)
             Hg(i,j,k,iene) = dot_product(mmtmp8(1:8,i), q(i,j-4:j+3,k,qtemp))
          end do
       end do
    end do

    do k=lo(3),hi(3)
       do j=slo(2),shi(2)+1
          do i=lo(1),hi(1)
             mmtmp8(1:8,i) = matmul(M8T, q(i,j-4:j+3,k,qpres))
             Hg(i,j,k,iene) = Hg(i,j,k,iene) + dot_product(dpe(i,j-4:j+3,k), mmtmp8(1:8,i))
          end do
          do i=lo(1),hi(1)
             M8p(:,i,j,k) = mmtmp8(1:8,i)
          end do
       end do
    end do

    do n = 1, nspecies

       if (n .eq. iias) cycle  ! inactive speices

       qxn = qx1+n-1
       iryn = iry1+n-1

       do k=lo(3),hi(3)
          do j=slo(2),shi(2)+1
             do i=lo(1),hi(1)
                Hg(i,j,k,iryn) = dot_product(dpy(i,j-4:j+3,k,n), M8p(:,i,j,k))
             end do
          end do
       end do

       do k=lo(3),hi(3)
          do j=slo(2),shi(2)+1
             do i=lo(1),hi(1)
                mmtmp8(1:8,i) = matmul(M8T, q(i,j-4:j+3,k,qxn))
                Hg(i,j,k,iene) = Hg(i,j,k,iene) + dot_product(dxe(i,j-4:j+3,k,n), mmtmp8(1:8,i))
                Hg(i,j,k,iryn) = Hg(i,j,k,iryn) &
                     + dot_product(dxy(i,j-4:j+3,k,n), mmtmp8(1:8,i))
             end do
          end do
       end do

    end do
       
    ! add y-direction rhs

    do n=2,iene
       do k=lo(3),hi(3)
          do j=slo(2),shi(2)
             do i=lo(1),hi(1)
                rhs(i,j,k,n) = rhs(i,j,k,n) + (Hg(i,j+1,k,n) - Hg(i,j,k,n)) * dx2inv(2)
             end do
          end do
       end do
    end do

    sumdrY = 0.d0
    do n=iry1,ncons

       if (n.eq.iry_ias) cycle

       do k=lo(3),hi(3)
          do j=slo(2),shi(2)
             do i=lo(1),hi(1)
                sumdry(i,j,k) = sumdry(i,j,k) + (Hg(i,j+1,k,n) - Hg(i,j,k,n)) * dx2inv(2)
                rhs(i,j,k,n)  =  rhs(i,j,k,n) + (Hg(i,j+1,k,n) - Hg(i,j,k,n)) * dx2inv(2)
             end do
          end do
       end do

    end do

    do k=lo(3),hi(3)
       do j=slo(2),shi(2)
          do i=lo(1),hi(1)
             gradp(i,j,k) = dxinv(2) * first_deriv_8(q(i,j-4:j+4,k,qpres))
          end do
       end do
    end do
    
    sumryv = 0.d0
    do n = 1, nspecies

       if (n.eq.iias) cycle

       qxn = qx1+n-1
       do k=lo(3),hi(3)
          do j=slo(2),shi(2)
             do i=lo(1),hi(1)
                sumryv(i,j,k) = sumryv(i,j,k) + dpy(i,j,k,n)*gradp(i,j,k)  &
                     + dxy(i,j,k,n)*dxinv(2)*first_deriv_8(q(i,j-4:j+4,k,qxn))
             end do
          end do
       end do

    end do

    if (add_v_correction) then

       do n=1,nspecies
          qyn = qy1+n-1
          qhn = qh1+n-1
          iryn = iry1+n-1

          do k=lo(3),hi(3)
             do j=slo(2),shi(2)
                do i=lo(1),hi(1)
                   ry_c = q(i,j,k,qyn)*sumdry(i,j,k) + sumryv(i,j,k)*dxinv(2) * &
                        first_deriv_8(q(i,j-4:j+4,k,qyn))
                   ene_c = ry_c*q(i,j,k,qhn) + q(i,j,k,qyn)*sumryv(i,j,k)*dxinv(2)* &
                        first_deriv_8(q(i,j-4:j+4,k,qhn))
                   rhs(i,j,k,iene) = rhs(i,j,k,iene) - ene_c
                   rhs(i,j,k,iryn) = rhs(i,j,k,iryn) - ry_c
                end do
             end do
          end do
       end do

    else

       n = iias
       qhn = qh1+n-1
       iryn = iry1+n-1

       do k=lo(3),hi(3)
          do j=slo(2),shi(2)
             do i=lo(1),hi(1)
                ene_c = sumdry(i,j,k)*q(i,j,k,qhn) + sumryv(i,j,k)*dxinv(2)* &
                     first_deriv_8(q(i,j-4:j+4,k,qhn))
                rhs(i,j,k,iene) = rhs(i,j,k,iene) - ene_c
                rhs(i,j,k,iryn) = rhs(i,j,k,iryn) - sumdry(i,j,k)
             end do
          end do
       end do

    end if

    ! ------- END y-direction -------

    ! ------- BEGIN z-direction -------

    do k=slo(3),shi(3)+1
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             mmtmp8(1:8,i) = matmul(mu(i,j,k-4:k+3), M8)
             Hg(i,j,k,imx) = dot_product(mmtmp8(1:8,i), q(i,j,k-4:k+3,qu))
             Hg(i,j,k,imy) = dot_product(mmtmp8(1:8,i), q(i,j,k-4:k+3,qv))
          end do
       end do
    end do

    do k=slo(3),shi(3)+1
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             mmtmp8(1:8,i) = matmul(vsp(i,j,k-4:k+3), M8)
             Hg(i,j,k,imz) = dot_product(mmtmp8(1:8,i), q(i,j,k-4:k+3,qw))
          end do
       end do
    end do

    do k=slo(3),shi(3)+1
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             mmtmp8(1:8,i) = matmul(lam(i,j,k-4:k+3), M8)
             Hg(i,j,k,iene) = dot_product(mmtmp8(1:8,i), q(i,j,k-4:k+3,qtemp))
          end do
       end do
    end do

    do k=slo(3),shi(3)+1
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             mmtmp8(1:8,i) = matmul(M8T, q(i,j,k-4:k+3,qpres))
             Hg(i,j,k,iene) = Hg(i,j,k,iene) + dot_product(dpe(i,j,k-4:k+3), mmtmp8(1:8,i))
          end do
          do i=lo(1),hi(1)
             M8p(:,i,j,k) = mmtmp8(1:8,i)
          end do
       end do
    end do

    do n = 1, nspecies

       if (n .eq. iias) cycle  ! inactive speices

       qxn = qx1+n-1
       iryn = iry1+n-1

       do k=slo(3),shi(3)+1
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                Hg(i,j,k,iryn) = dot_product(dpy(i,j,k-4:k+3,n), M8p(:,i,j,k))
             end do
          end do
       end do

       do k=slo(3),shi(3)+1
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                mmtmp8(1:8,i) = matmul(M8T, q(i,j,k-4:k+3,qxn))
                Hg(i,j,k,iene) = Hg(i,j,k,iene) + dot_product(dxe(i,j,k-4:k+3,n), mmtmp8(1:8,i))
                Hg(i,j,k,iryn) = Hg(i,j,k,iryn) &
                     + dot_product(dxy(i,j,k-4:k+3,n), mmtmp8(1:8,i))
             end do
          end do
       end do

    end do
    
    ! add z-direction rhs

    do n=2,iene
       do k=slo(3),shi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                rhs(i,j,k,n) = rhs(i,j,k,n) + (Hg(i,j,k+1,n) - Hg(i,j,k,n)) * dx2inv(3)
             end do
          end do
       end do
    end do
    
    sumdrY = 0.d0
    do n=iry1,ncons

       if (n.eq.iry_ias) cycle

       do k=slo(3),shi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                sumdry(i,j,k) = sumdry(i,j,k) + (Hg(i,j,k+1,n) - Hg(i,j,k,n)) * dx2inv(3)
                rhs(i,j,k,n)  =  rhs(i,j,k,n) + (Hg(i,j,k+1,n) - Hg(i,j,k,n)) * dx2inv(3)
             end do
          end do
       end do

    end do
    
    do k=slo(3),shi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             gradp(i,j,k) = dxinv(3) * first_deriv_8(q(i,j,k-4:k+4,qpres))
          end do
       end do
    end do
    
    sumryv = 0.d0
    do n = 1, nspecies

       if (n.eq.iias) cycle

       qxn = qx1+n-1
       do k=slo(3),shi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                sumryv(i,j,k) = sumryv(i,j,k) + dpy(i,j,k,n)*gradp(i,j,k)  &
                     + dxy(i,j,k,n)*dxinv(3)*first_deriv_8(q(i,j,k-4:k+4,qxn))
             end do
          end do
       end do

    end do

    if (add_v_correction) then

       do n=1,nspecies
          qyn = qy1+n-1
          qhn = qh1+n-1
          iryn = iry1+n-1

          do k=slo(3),shi(3)
             do j=lo(2),hi(2)
                do i=lo(1),hi(1)
                   ry_c = q(i,j,k,qyn)*sumdry(i,j,k) + sumryv(i,j,k)*dxinv(3) * &
                        first_deriv_8(q(i,j,k-4:k+4,qyn))
                   ene_c = ry_c*q(i,j,k,qhn) + q(i,j,k,qyn)*sumryv(i,j,k)*dxinv(3)* &
                        first_deriv_8(q(i,j,k-4:k+4,qhn))
                   rhs(i,j,k,iene) = rhs(i,j,k,iene) - ene_c
                   rhs(i,j,k,iryn) = rhs(i,j,k,iryn) - ry_c
                end do
             end do
          end do
       end do

    else

       n = iias
       qhn = qh1+n-1
       iryn = iry1+n-1

       do k=slo(3),shi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                ene_c = sumdry(i,j,k)*q(i,j,k,qhn) + sumryv(i,j,k)*dxinv(3)* &
                     first_deriv_8(q(i,j,k-4:k+4,qhn))
                rhs(i,j,k,iene) = rhs(i,j,k,iene) - ene_c
                rhs(i,j,k,iryn) = rhs(i,j,k,iryn) - sumdry(i,j,k)
             end do
          end do
       end do

    end if

    ! ------- END z-direction -------

    !
    ! lo-x boundary
    !
    if (physbclo(1)) then
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             i = lo(1)
             ! use completely right-biased stencil
             mmtmpB = matmul(vsp(i:i+3,j,k), BRB)
             rhs(i,j,k,imx) = rhs(i,j,k,imx)+finlo(1)*dx2inv(1)*dot_product(mmtmpB,q(i:i+3,j,k,qu))

             mmtmpB = matmul(mu(i:i+3,j,k), BRB)
             rhs(i,j,k,imy) = rhs(i,j,k,imy)+foulo(1)*dx2inv(1)*dot_product(mmtmpB,q(i:i+3,j,k,qv))
             rhs(i,j,k,imz) = rhs(i,j,k,imz)+foulo(1)*dx2inv(1)*dot_product(mmtmpB,q(i:i+3,j,k,qw))

             mmtmpB = matmul(lam(i:i+3,j,k), BRB)
             BBp = matmul(BRB, q(i:i+3,j,k,qpres))
             rhs(i,j,k,iene) = rhs(i,j,k,iene) + foulo(1)*dx2inv(1) * &
                  ( dot_product(mmtmpB, q(i:i+3,j,k,qtemp)) &
                  + dot_product(      dpe(i:i+3,j,k), BBp) )
             
             rhstot = 0.d0
             rhsene = 0.d0
             do n = 1, nspecies
                qxn = qx1+n-1
                qyn = qy1+n-1
                
                BBX = matmul(BRB, q(i:i+3,j,k,qxn))
                
                rhstmp(n) = dot_product(dpy(i:i+3,j,k,n), BBp) &
                     +      dot_product(dxy(i:i+3,j,k,n), BBX)
                
                rhsene = rhsene &
                     +      dot_product(dxe(i:i+3,j,k,n), BBX)
                
                rhstot = rhstot + rhstmp(n)
                Ytmp(n) = q(i,j,k,qyn)
             end do
             
             do n = 1, nspecies
                rhs(i,j,k,iry1+n-1) =  rhs(i,j,k,iry1+n-1) + &
                     foulo(1)*dx2inv(1) * (rhstmp(n) - Ytmp(n)*rhstot)
             end do
             
             do n = 1, nspecies
                qhn = qh1+n-1
                rhsene = rhsene - Ytmp(n) * q(i,j,k,qhn) * rhstot
             end do
             rhs(i,j,k,iene) = rhs(i,j,k,iene) + foulo(1)*dx2inv(1) * rhsene


             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! use 2nd-order stencil for cell lo(1)+1,j,k
             do iface=0,1 
                i = lo(1)+1 + iface

                mmtmp2 = matmul(vsp(i-1:i,j,k), M2)
                Hcell(iface,imx) = dot_product(mmtmp2, q(i-1:i,j,k,qu))

                mmtmp2 = matmul(mu(i-1:i,j,k), M2)
                Hcell(iface,imy) = dot_product(mmtmp2, q(i-1:i,j,k,qv))
                Hcell(iface,imz) = dot_product(mmtmp2, q(i-1:i,j,k,qw))

                mmtmp2 = matmul(lam(i-1:i,j,k), M2)
                M2p = matmul(M2,  q(i-1:i,j,k,qpres))
                Hcell(iface,iene) = dot_product(mmtmp2, q(i-1:i,j,k,qtemp)) &
                     &            + dot_product(      dpe(i-1:i,j,k), M2p)

                Htot = 0.d0
                do n = 1, nspecies
                   qxn = qx1+n-1
                   qyn = qy1+n-1

                   M2X = matmul(M2, q(i-1:i,j,k,qxn))
                
                   Htmp(n) = dot_product(dpy(i-1:i,j,k,n), M2p) &
                        +    dot_product(dxy(i-1:i,j,k,n), M2X)

                   Hcell(iface,iene) = Hcell(iface,iene) &
                        +    dot_product(dxe(i-1:i,j,k,n), M2X)

                   Htot = Htot + Htmp(n)
                   Ytmp(n) = 0.5d0*(q(i-1,j,k,qyn) + q(i,j,k,qyn))
                end do

                do n = 1, nspecies
                   Hcell(iface,iry1+n-1) = Htmp(n) - Ytmp(n)*Htot
                end do

                do n = 1, nspecies
                   qhn = qh1+n-1
                   hhalf = 0.5d0*(q(i-1,j,k,qhn) + q(i,j,k,qhn))
                   Hcell(iface,iene) =  Hcell(iface,iene) - Ytmp(n)*Htot*hhalf 
                end do
             end do

             i = lo(1)+1
             do n=2,ncons
                rhs(i,j,k,n) = rhs(i,j,k,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(1)
             end do

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! use 4th-order stencil for cell lo(1)+2,j,k
             do iface=0,1 
                i = lo(1)+2 + iface

                mmtmp4 = matmul(vsp(i-2:i+1,j,k), M4)
                Hcell(iface,imx) = dot_product(mmtmp4, q(i-2:i+1,j,k,qu))

                mmtmp4 = matmul(mu(i-2:i+1,j,k), M4)
                Hcell(iface,imy) = dot_product(mmtmp4, q(i-2:i+1,j,k,qv))
                Hcell(iface,imz) = dot_product(mmtmp4, q(i-2:i+1,j,k,qw))

                mmtmp4 = matmul(lam(i-2:i+1,j,k), M4)
                M4p = matmul(M4T,  q(i-2:i+1,j,k,qpres))
                Hcell(iface,iene) = dot_product(mmtmp4, q(i-2:i+1,j,k,qtemp))      &
                     &            + dot_product(      dpe(i-2:i+1,j,k), M4p)

                do n = 1, nspecies
                   if (n .eq. iias) cycle  ! inactive speices
                   
                   qxn = qx1+n-1
                   iryn = iry1+n-1
                   M4X = matmul(M4T, q(i-2:i+1,j,k,qxn))
                   Hcell(iface,iene) = Hcell(iface,iene) &
                        +    dot_product(dxe(i-2:i+1,j,k,n), M4X)
                   Hcell(iface,iryn) = dot_product(dpy(i-2:i+1,j,k,n), M4p) &
                        + dot_product(dxy(i-2:i+1,j,k,n), M4X)
                end do

             end do

             i = lo(1)+2

             do n=2,iene
                rhs(i,j,k,n) = rhs(i,j,k,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(1)
             end do

             sumdrytmp = 0.d0
             do n=iry1,ncons
                if (n.eq.iry_ias) cycle
                sumdrytmp = sumdrytmp + (Hcell(1,n) - Hcell(0,n)) * dx2inv(1)
                rhs(i,j,k,n)  =  rhs(i,j,k,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(1)
             end do

             gradptmp = dxinv(1) * first_deriv_4(q(i-2:i+2,j,k,qpres))

             sumryvtmp = 0.d0
             do n = 1, nspecies
                if (n.eq.iias) cycle
                qxn = qx1+n-1
                sumryvtmp = sumryvtmp + dpy(i,j,k,n)*gradptmp  &
                     + dxy(i,j,k,n)*dxinv(1)*first_deriv_4(q(i-2:i+2,j,k,qxn))
             end do

             if (add_v_correction) then

                do n=1,nspecies
                   qyn = qy1+n-1
                   qhn = qh1+n-1
                   iryn = iry1+n-1
          
                   ry_c = q(i,j,k,qyn)*sumdrytmp + sumryvtmp*dxinv(1) * &
                        first_deriv_4(q(i-2:i+2,j,k,qyn))
                   ene_c = ry_c*q(i,j,k,qhn) + q(i,j,k,qyn)*sumryvtmp*dxinv(1)* &
                        first_deriv_4(q(i-2:i+2,j,k,qhn))
                   rhs(i,j,k,iene) = rhs(i,j,k,iene) - ene_c
                   rhs(i,j,k,iryn) = rhs(i,j,k,iryn) - ry_c
                end do

             else
    
                n = iias
                qhn = qh1+n-1
                iryn = iry1+n-1

                ene_c = sumdrytmp*q(i,j,k,qhn) + sumryvtmp*dxinv(1)* &
                     first_deriv_4(q(i-2:i+2,j,k,qhn))
                rhs(i,j,k,iene) = rhs(i,j,k,iene) - ene_c
                rhs(i,j,k,iryn) = rhs(i,j,k,iryn) - sumdrytmp
       
             end if

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! use 6th-order stencil for cell lo(1)+3,j,k
             do iface=0,1 
                i = lo(1)+3 + iface

                mmtmp6 = matmul(vsp(i-3:i+2,j,k), M6)
                Hcell(iface,imx) = dot_product(mmtmp6, q(i-3:i+2,j,k,qu))

                mmtmp6 = matmul(mu(i-3:i+2,j,k), M6)
                Hcell(iface,imy) = dot_product(mmtmp6, q(i-3:i+2,j,k,qv))
                Hcell(iface,imz) = dot_product(mmtmp6, q(i-3:i+2,j,k,qw))

                mmtmp6 = matmul(lam(i-3:i+2,j,k), M6)
                M6p = matmul(M6T,  q(i-3:i+2,j,k,qpres))
                Hcell(iface,iene) = dot_product(mmtmp6, q(i-3:i+2,j,k,qtemp)) &
                     &            + dot_product(      dpe(i-3:i+2,j,k), M6p)

                do n = 1, nspecies
                   if (n .eq. iias) cycle  ! inactive speices
                   qxn = qx1+n-1
                   iryn = iry1+n-1
                   M6X = matmul(M6T, q(i-3:i+2,j,k,qxn))
                   Hcell(iface,iene) = Hcell(iface,iene) &
                        +    dot_product(dxe(i-3:i+2,j,k,n), M6X)
                   Hcell(iface,iryn) = dot_product(dpy(i-3:i+2,j,k,n), M6p) &
                        +    dot_product(dxy(i-3:i+2,j,k,n), M6X)
                end do

             end do

             i = lo(1)+3

             do n=2,iene
                rhs(i,j,k,n) = rhs(i,j,k,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(1)
             end do

             sumdrytmp = 0.d0
             do n=iry1,ncons
                if (n.eq.iry_ias) cycle
                sumdrytmp = sumdrytmp + (Hcell(1,n) - Hcell(0,n)) * dx2inv(1)
                rhs(i,j,k,n)  =  rhs(i,j,k,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(1)
             end do

             gradptmp = dxinv(1) * first_deriv_6(q(i-3:i+3,j,k,qpres))

             sumryvtmp = 0.d0
             do n = 1, nspecies
                if (n.eq.iias) cycle
                qxn = qx1+n-1
                sumryvtmp = sumryvtmp + dpy(i,j,k,n)*gradptmp  &
                     + dxy(i,j,k,n)*dxinv(1)*first_deriv_6(q(i-3:i+3,j,k,qxn))
             end do

             if (add_v_correction) then

                do n=1,nspecies
                   qyn = qy1+n-1
                   qhn = qh1+n-1
                   iryn = iry1+n-1
          
                   ry_c = q(i,j,k,qyn)*sumdrytmp + sumryvtmp*dxinv(1) * &
                        first_deriv_6(q(i-3:i+3,j,k,qyn))
                   ene_c = ry_c*q(i,j,k,qhn) + q(i,j,k,qyn)*sumryvtmp*dxinv(1)* &
                        first_deriv_6(q(i-3:i+3,j,k,qhn))
                   rhs(i,j,k,iene) = rhs(i,j,k,iene) - ene_c
                   rhs(i,j,k,iryn) = rhs(i,j,k,iryn) - ry_c
                end do

             else
    
                n = iias
                qhn = qh1+n-1
                iryn = iry1+n-1

                ene_c = sumdrytmp*q(i,j,k,qhn) + sumryvtmp*dxinv(1)* &
                     first_deriv_6(q(i-3:i+3,j,k,qhn))
                rhs(i,j,k,iene) = rhs(i,j,k,iene) - ene_c
                rhs(i,j,k,iryn) = rhs(i,j,k,iryn) - sumdrytmp
       
             end if

          end do
       end do
    end if

    !
    ! hi-x boundary
    !
    if (physbchi(1)) then
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! use 6th-order stencil for cell hi(1)-3,j,k
             do iface=0,1  ! two faces of 
                i = hi(1)-3 + iface

                mmtmp6 = matmul(vsp(i-3:i+2,j,k), M6)
                Hcell(iface,imx) = dot_product(mmtmp6, q(i-3:i+2,j,k,qu))

                mmtmp6 = matmul(mu(i-3:i+2,j,k), M6)
                Hcell(iface,imy) = dot_product(mmtmp6, q(i-3:i+2,j,k,qv))
                Hcell(iface,imz) = dot_product(mmtmp6, q(i-3:i+2,j,k,qw))

                mmtmp6 = matmul(lam(i-3:i+2,j,k), M6)
                M6p = matmul(M6T,  q(i-3:i+2,j,k,qpres))
                Hcell(iface,iene) = dot_product(mmtmp6, q(i-3:i+2,j,k,qtemp)) &
                     &            + dot_product(      dpe(i-3:i+2,j,k), M6p)

                do n = 1, nspecies
                   if (n .eq. iias) cycle  ! inactive speices
                   qxn = qx1+n-1
                   iryn = iry1+n-1
                   M6X = matmul(M6T, q(i-3:i+2,j,k,qxn))
                   Hcell(iface,iene) = Hcell(iface,iene) &
                        +    dot_product(dxe(i-3:i+2,j,k,n), M6X)                
                   Hcell(iface,iryn) = dot_product(dpy(i-3:i+2,j,k,n), M6p) &
                        +    dot_product(dxy(i-3:i+2,j,k,n), M6X)
                end do

             end do

             i = hi(1)-3

             do n=2,iene
                rhs(i,j,k,n) = rhs(i,j,k,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(1)
             end do

             sumdrytmp = 0.d0
             do n=iry1,ncons
                if (n.eq.iry_ias) cycle
                sumdrytmp = sumdrytmp + (Hcell(1,n) - Hcell(0,n)) * dx2inv(1)
                rhs(i,j,k,n)  =  rhs(i,j,k,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(1)
             end do

             gradptmp = dxinv(1) * first_deriv_6(q(i-3:i+3,j,k,qpres))

             sumryvtmp = 0.d0
             do n = 1, nspecies
                if (n.eq.iias) cycle
                qxn = qx1+n-1
                sumryvtmp = sumryvtmp + dpy(i,j,k,n)*gradptmp  &
                     + dxy(i,j,k,n)*dxinv(1)*first_deriv_6(q(i-3:i+3,j,k,qxn))
             end do

             if (add_v_correction) then

                do n=1,nspecies
                   qyn = qy1+n-1
                   qhn = qh1+n-1
                   iryn = iry1+n-1
          
                   ry_c = q(i,j,k,qyn)*sumdrytmp + sumryvtmp*dxinv(1) * &
                        first_deriv_6(q(i-3:i+3,j,k,qyn))
                   ene_c = ry_c*q(i,j,k,qhn) + q(i,j,k,qyn)*sumryvtmp*dxinv(1)* &
                        first_deriv_6(q(i-3:i+3,j,k,qhn))
                   rhs(i,j,k,iene) = rhs(i,j,k,iene) - ene_c
                   rhs(i,j,k,iryn) = rhs(i,j,k,iryn) - ry_c
                end do

             else
    
                n = iias
                qhn = qh1+n-1
                iryn = iry1+n-1

                ene_c = sumdrytmp*q(i,j,k,qhn) + sumryvtmp*dxinv(1)* &
                     first_deriv_6(q(i-3:i+3,j,k,qhn))
                rhs(i,j,k,iene) = rhs(i,j,k,iene) - ene_c
                rhs(i,j,k,iryn) = rhs(i,j,k,iryn) - sumdrytmp
       
             end if

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! use 4th-order stencil for cell hi(1)-2,j,k
             do iface=0,1 
                i = hi(1)-2 + iface

                mmtmp4 = matmul(vsp(i-2:i+1,j,k), M4)
                Hcell(iface,imx) = dot_product(mmtmp4, q(i-2:i+1,j,k,qu))

                mmtmp4 = matmul(mu(i-2:i+1,j,k), M4)
                Hcell(iface,imy) = dot_product(mmtmp4, q(i-2:i+1,j,k,qv))
                Hcell(iface,imz) = dot_product(mmtmp4, q(i-2:i+1,j,k,qw))

                mmtmp4 = matmul(lam(i-2:i+1,j,k), M4)
                M4p = matmul(M4T,  q(i-2:i+1,j,k,qpres))
                Hcell(iface,iene) = dot_product(mmtmp4, q(i-2:i+1,j,k,qtemp)) &
                     &            + dot_product(      dpe(i-2:i+1,j,k), M4p)

                do n = 1, nspecies
                   if (n .eq. iias) cycle  ! inactive speices
                   qxn = qx1+n-1
                   iryn = iry1+n-1

                   M4X = matmul(M4T, q(i-2:i+1,j,k,qxn))
                   Hcell(iface,iene) = Hcell(iface,iene) &
                        +    dot_product(dxe(i-2:i+1,j,k,n), M4X)
                   Hcell(iface,iryn) = dot_product(dpy(i-2:i+1,j,k,n), M4p) &
                        +    dot_product(dxy(i-2:i+1,j,k,n), M4X)
                end do

             end do

             i = hi(1)-2

             do n=2,iene
                rhs(i,j,k,n) = rhs(i,j,k,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(1)
             end do

             sumdrytmp = 0.d0
             do n=iry1,ncons
                if (n.eq.iry_ias) cycle
                sumdrytmp = sumdrytmp + (Hcell(1,n) - Hcell(0,n)) * dx2inv(1)
                rhs(i,j,k,n)  =  rhs(i,j,k,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(1)
             end do

             gradptmp = dxinv(1) * first_deriv_4(q(i-2:i+2,j,k,qpres))

             sumryvtmp = 0.d0
             do n = 1, nspecies
                if (n.eq.iias) cycle
                qxn = qx1+n-1
                sumryvtmp = sumryvtmp + dpy(i,j,k,n)*gradptmp  &
                     + dxy(i,j,k,n)*dxinv(1)*first_deriv_4(q(i-2:i+2,j,k,qxn))
             end do

             if (add_v_correction) then

                do n=1,nspecies
                   qyn = qy1+n-1
                   qhn = qh1+n-1
                   iryn = iry1+n-1
          
                   ry_c = q(i,j,k,qyn)*sumdrytmp + sumryvtmp*dxinv(1) * &
                        first_deriv_4(q(i-2:i+2,j,k,qyn))
                   ene_c = ry_c*q(i,j,k,qhn) + q(i,j,k,qyn)*sumryvtmp*dxinv(1)* &
                        first_deriv_4(q(i-2:i+2,j,k,qhn))
                   rhs(i,j,k,iene) = rhs(i,j,k,iene) - ene_c
                   rhs(i,j,k,iryn) = rhs(i,j,k,iryn) - ry_c
                end do

             else
    
                n = iias
                qhn = qh1+n-1
                iryn = iry1+n-1

                ene_c = sumdrytmp*q(i,j,k,qhn) + sumryvtmp*dxinv(1)* &
                     first_deriv_4(q(i-2:i+2,j,k,qhn))
                rhs(i,j,k,iene) = rhs(i,j,k,iene) - ene_c
                rhs(i,j,k,iryn) = rhs(i,j,k,iryn) - sumdrytmp
       
             end if

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! use 2nd-order stencil for cell hi(1)-1,j,k
             do iface=0,1 
                i = hi(1)-1 + iface

                mmtmp2 = matmul(vsp(i-1:i,j,k), M2)
                Hcell(iface,imx) = dot_product(mmtmp2, q(i-1:i,j,k,qu))

                mmtmp2 = matmul(mu(i-1:i,j,k), M2)
                Hcell(iface,imy) = dot_product(mmtmp2, q(i-1:i,j,k,qv))
                Hcell(iface,imz) = dot_product(mmtmp2, q(i-1:i,j,k,qw))

                mmtmp2 = matmul(lam(i-1:i,j,k), M2)
                M2p = matmul(M2,  q(i-1:i,j,k,qpres))
                Hcell(iface,iene) = dot_product(mmtmp2, q(i-1:i,j,k,qtemp)) &
                     &            + dot_product(      dpe(i-1:i,j,k), M2p)

                Htot = 0.d0
                do n = 1, nspecies
                   qxn = qx1+n-1
                   qyn = qy1+n-1

                   M2X = matmul(M2, q(i-1:i,j,k,qxn))
                
                   Htmp(n) = dot_product(dpy(i-1:i,j,k,n), M2p) &
                        +    dot_product(dxy(i-1:i,j,k,n), M2X)

                   Hcell(iface,iene) = Hcell(iface,iene) &
                        +    dot_product(dxe(i-1:i,j,k,n), M2X)

                   Htot = Htot + Htmp(n)
                   Ytmp(n) = 0.5d0*(q(i-1,j,k,qyn) + q(i,j,k,qyn))
                end do

                do n = 1, nspecies
                   Hcell(iface,iry1+n-1) = Htmp(n) - Ytmp(n)*Htot
                end do

                do n = 1, nspecies
                   qhn = qh1+n-1
                   hhalf = 0.5d0*(q(i-1,j,k,qhn) + q(i,j,k,qhn))
                   Hcell(iface,iene) =  Hcell(iface,iene) - Ytmp(n)*Htot*hhalf
                end do
             end do

             i = hi(1)-1
             do n=2,ncons
                rhs(i,j,k,n) = rhs(i,j,k,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(1)
             end do

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             i = hi(1)
             ! use completely left-biased stencil
             mmtmpB = matmul(vsp(i-3:i,j,k), BLB)
             rhs(i,j,k,imx) = rhs(i,j,k,imx)+finhi(1)*dx2inv(1)*dot_product(mmtmpB,q(i-3:i,j,k,qu))

             mmtmpB= matmul(mu(i-3:i,j,k), BLB)
             rhs(i,j,k,imy) = rhs(i,j,k,imy)+fouhi(1)*dx2inv(1)*dot_product(mmtmpB,q(i-3:i,j,k,qv))
             rhs(i,j,k,imz) = rhs(i,j,k,imz)+fouhi(1)*dx2inv(1)*dot_product(mmtmpB,q(i-3:i,j,k,qw))

             mmtmpB = matmul(lam(i-3:i,j,k), BLB)
             BBp = matmul(BLB, q(i-3:i,j,k,qpres))
             rhs(i,j,k,iene) = rhs(i,j,k,iene) + fouhi(1)*dx2inv(1) * &
                  ( dot_product(mmtmpB, q(i-3:i,j,k,qtemp)) &
                  + dot_product(      dpe(i-3:i,j,k), BBp) )
             
             rhstot = 0.d0
             rhsene = 0.d0
             do n = 1, nspecies
                qxn = qx1+n-1
                qyn = qy1+n-1

                BBX = matmul(BLB, q(i-3:i,j,k,qxn))
                
                rhstmp(n) = dot_product(dpy(i-3:i,j,k,n), BBp) &
                     +      dot_product(dxy(i-3:i,j,k,n), BBX)
                
                rhsene = rhsene &
                     +      dot_product(dxe(i-3:i,j,k,n), BBX)

                rhstot = rhstot + rhstmp(n)
                Ytmp(n) = q(i,j,k,qyn)
             end do
             
             do n = 1, nspecies
                rhs(i,j,k,iry1+n-1) =  rhs(i,j,k,iry1+n-1) + &
                     fouhi(1)*dx2inv(1) * (rhstmp(n) - Ytmp(n)*rhstot)
             end do

             do n = 1, nspecies
                qhn = qh1+n-1
                rhsene = rhsene - Ytmp(n) * q(i,j,k,qhn) * rhstot
             end do
             rhs(i,j,k,iene) = rhs(i,j,k,iene) + fouhi(1)*dx2inv(1) * rhsene
          end do
       end do
    end if

    !
    ! lo-y boundary
    !
    if (physbclo(2)) then
       do k=lo(3),hi(3)
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          j = lo(2)
          ! use completely right-biased stencil
          do i=lo(1),hi(1)

             mmtmpB = matmul(mu(i,j:j+3,k), BRB)
             rhs(i,j,k,imx) = rhs(i,j,k,imx)+foulo(2)*dx2inv(2)*dot_product(mmtmpB,q(i,j:j+3,k,qu))
             rhs(i,j,k,imz) = rhs(i,j,k,imz)+foulo(2)*dx2inv(2)*dot_product(mmtmpB,q(i,j:j+3,k,qw))

             mmtmpB = matmul(vsp(i,j:j+3,k), BRB)
             rhs(i,j,k,imy) = rhs(i,j,k,imy)+finlo(2)*dx2inv(2)*dot_product(mmtmpB,q(i,j:j+3,k,qv))

             mmtmpB = matmul(lam(i,j:j+3,k), BRB)
             BBp = matmul(BRB, q(i,j:j+3,k,qpres))
             rhs(i,j,k,iene) = rhs(i,j,k,iene) + foulo(2)*dx2inv(2) * &
                  ( dot_product(mmtmpB, q(i,j:j+3,k,qtemp)) &
                  + dot_product(      dpe(i,j:j+3,k), BBp) )
             
             rhstot = 0.d0
             rhsene = 0.d0
             do n = 1, nspecies
                qxn = qx1+n-1
                qyn = qy1+n-1

                BBX = matmul(BRB, q(i,j:j+3,k,qxn))
                
                rhstmp(n) = dot_product(dpy(i,j:j+3,k,n), BBp) &
                     +      dot_product(dxy(i,j:j+3,k,n), BBX)
                
                rhsene = rhsene &
                     +      dot_product(dxe(i,j:j+3,k,n), BBX)

                rhstot = rhstot + rhstmp(n)
                Ytmp(n) = q(i,j,k,qyn)
             end do
             
             do n = 1, nspecies
                rhs(i,j,k,iry1+n-1) =  rhs(i,j,k,iry1+n-1) + &
                     foulo(2)*dx2inv(2) * (rhstmp(n) - Ytmp(n)*rhstot)
             end do

             do n = 1, nspecies
                qhn = qh1+n-1
                rhsene = rhsene - Ytmp(n) * q(i,j,k,qhn) * rhstot
             end do
             rhs(i,j,k,iene) = rhs(i,j,k,iene) + foulo(2)*dx2inv(2) * rhsene
          end do

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! use 2nd-order stencil for cell i,lo(2)+1,k
          do i=lo(1),hi(1)
             do iface=0,1 
                j = lo(2)+1 + iface

                mmtmp2 = matmul(mu(i,j-1:j,k), M2)
                Hcell(iface,imx) = dot_product(mmtmp2, q(i,j-1:j,k,qu))
                Hcell(iface,imz) = dot_product(mmtmp2, q(i,j-1:j,k,qw))

                mmtmp2 = matmul(vsp(i,j-1:j,k), M2)
                Hcell(iface,imy) = dot_product(mmtmp2, q(i,j-1:j,k,qv))

                mmtmp2 = matmul(lam(i,j-1:j,k), M2)
                M2p = matmul(M2,  q(i,j-1:j,k,qpres))
                Hcell(iface,iene) = dot_product(mmtmp2, q(i,j-1:j,k,qtemp)) &
                     &            + dot_product(      dpe(i,j-1:j,k), M2p)

                Htot = 0.d0
                do n = 1, nspecies
                   qxn = qx1+n-1
                   qyn = qy1+n-1

                   M2X = matmul(M2, q(i,j-1:j,k,qxn))
                
                   Htmp(n) = dot_product(dpy(i,j-1:j,k,n), M2p) &
                        +    dot_product(dxy(i,j-1:j,k,n), M2X)

                   Hcell(iface,iene) = Hcell(iface,iene) &
                        +    dot_product(dxe(i,j-1:j,k,n), M2X)

                   Htot = Htot + Htmp(n)
                   Ytmp(n) = 0.5d0*(q(i,j-1,k,qyn) + q(i,j,k,qyn))
                end do

                do n = 1, nspecies
                   Hcell(iface,iry1+n-1) = Htmp(n) - Ytmp(n)*Htot
                end do

                do n = 1, nspecies
                   qhn = qh1+n-1
                   hhalf = 0.5d0*(q(i,j-1,k,qhn) + q(i,j,k,qhn))
                   Hcell(iface,iene) =  Hcell(iface,iene) - Ytmp(n)*Htot*hhalf
                end do
             end do

             j = lo(2)+1
             do n=2,ncons
                rhs(i,j,k,n) = rhs(i,j,k,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(2)
             end do
          end do

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! use 4th-order stencil for cell i,lo(2)+2,k
          do i=lo(1),hi(1)
             do iface=0,1 
                j = lo(2)+2 + iface

                mmtmp4 = matmul(mu(i,j-2:j+1,k), M4)
                Hcell(iface,imx) = dot_product(mmtmp4, q(i,j-2:j+1,k,qu))
                Hcell(iface,imz) = dot_product(mmtmp4, q(i,j-2:j+1,k,qw))

                mmtmp4 = matmul(vsp(i,j-2:j+1,k), M4)
                Hcell(iface,imy) = dot_product(mmtmp4, q(i,j-2:j+1,k,qv))

                mmtmp4 = matmul(lam(i,j-2:j+1,k), M4)
                M4p = matmul(M4T,  q(i,j-2:j+1,k,qpres))
                Hcell(iface,iene) = dot_product(mmtmp4, q(i,j-2:j+1,k,qtemp)) &
                     &            + dot_product(      dpe(i,j-2:j+1,k), M4p)

                do n = 1, nspecies
                   if (n .eq. iias) cycle  ! inactive speices

                   qxn = qx1+n-1
                   iryn = iry1+n-1
                   M4X = matmul(M4T, q(i,j-2:j+1,k,qxn))
                   Hcell(iface,iene) = Hcell(iface,iene) &
                        +    dot_product(dxe(i,j-2:j+1,k,n), M4X)
                   Hcell(iface,iryn) = dot_product(dpy(i,j-2:j+1,k,n), M4p) &
                        +    dot_product(dxy(i,j-2:j+1,k,n), M4X)
                end do

             end do

             j = lo(2)+2

             do n=2,iene
                rhs(i,j,k,n) = rhs(i,j,k,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(2)
             end do

             sumdrytmp = 0.d0
             do n=iry1,ncons
                if (n.eq.iry_ias) cycle
                sumdrytmp = sumdrytmp + (Hcell(1,n) - Hcell(0,n)) * dx2inv(2)
                rhs(i,j,k,n)  =  rhs(i,j,k,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(2)
             end do

             gradptmp = dxinv(2) * first_deriv_4(q(i,j-2:j+2,k,qpres))

             sumryvtmp = 0.d0
             do n = 1, nspecies
                if (n.eq.iias) cycle
                qxn = qx1+n-1
                sumryvtmp = sumryvtmp + dpy(i,j,k,n)*gradptmp  &
                     + dxy(i,j,k,n)*dxinv(2)*first_deriv_4(q(i,j-2:j+2,k,qxn))
             end do

             if (add_v_correction) then

                do n=1,nspecies
                   qyn = qy1+n-1
                   qhn = qh1+n-1
                   iryn = iry1+n-1
          
                   ry_c = q(i,j,k,qyn)*sumdrytmp + sumryvtmp*dxinv(2) * &
                        first_deriv_4(q(i,j-2:j+2,k,qyn))
                   ene_c = ry_c*q(i,j,k,qhn) + q(i,j,k,qyn)*sumryvtmp*dxinv(2)* &
                        first_deriv_4(q(i,j-2:j+2,k,qhn))
                   rhs(i,j,k,iene) = rhs(i,j,k,iene) - ene_c
                   rhs(i,j,k,iryn) = rhs(i,j,k,iryn) - ry_c
                end do

             else
    
                n = iias
                qhn = qh1+n-1
                iryn = iry1+n-1

                ene_c = sumdrytmp*q(i,j,k,qhn) + sumryvtmp*dxinv(2)* &
                     first_deriv_4(q(i,j-2:j+2,k,qhn))
                rhs(i,j,k,iene) = rhs(i,j,k,iene) - ene_c
                rhs(i,j,k,iryn) = rhs(i,j,k,iryn) - sumdrytmp
       
             end if
          end do

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! use 6th-order stencil for cell i,lo(2)+3,k
          do i=lo(1),hi(1)
             do iface=0,1 
                j = lo(2)+3 + iface

                mmtmp6 = matmul(mu(i,j-3:j+2,k), M6)
                Hcell(iface,imx) = dot_product(mmtmp6, q(i,j-3:j+2,k,qu))
                Hcell(iface,imz) = dot_product(mmtmp6, q(i,j-3:j+2,k,qw))

                mmtmp6 = matmul(vsp(i,j-3:j+2,k), M6)
                Hcell(iface,imy) = dot_product(mmtmp6, q(i,j-3:j+2,k,qv))

                mmtmp6 = matmul(lam(i,j-3:j+2,k), M6)
                M6p = matmul(M6T,  q(i,j-3:j+2,k,qpres))
                Hcell(iface,iene) = dot_product(mmtmp6, q(i,j-3:j+2,k,qtemp))      &
                     &            + dot_product(      dpe(i,j-3:j+2,k), M6p)

                do n = 1, nspecies
                   if (n .eq. iias) cycle  ! inactive speices
                   qxn = qx1+n-1
                   iryn = iry1+n-1
                   M6X = matmul(M6T, q(i,j-3:j+2,k,qxn))
                   Hcell(iface,iene) = Hcell(iface,iene) &
                        +    dot_product(dxe(i,j-3:j+2,k,n), M6X)
                   Hcell(iface,iryn) = dot_product(dpy(i,j-3:j+2,k,n), M6p) &
                        +    dot_product(dxy(i,j-3:j+2,k,n), M6X)
                end do

             end do

             j = lo(2)+3

             do n=2,iene
                rhs(i,j,k,n) = rhs(i,j,k,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(2)
             end do

             sumdrytmp = 0.d0
             do n=iry1,ncons
                if (n.eq.iry_ias) cycle
                sumdrytmp = sumdrytmp + (Hcell(1,n) - Hcell(0,n)) * dx2inv(2)
                rhs(i,j,k,n)  =  rhs(i,j,k,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(2)
             end do

             gradptmp = dxinv(2) * first_deriv_6(q(i,j-3:j+3,k,qpres))

             sumryvtmp = 0.d0
             do n = 1, nspecies
                if (n.eq.iias) cycle
                qxn = qx1+n-1
                sumryvtmp = sumryvtmp + dpy(i,j,k,n)*gradptmp  &
                     + dxy(i,j,k,n)*dxinv(2)*first_deriv_6(q(i,j-3:j+3,k,qxn))
             end do

             if (add_v_correction) then

                do n=1,nspecies
                   qyn = qy1+n-1
                   qhn = qh1+n-1
                   iryn = iry1+n-1
          
                   ry_c = q(i,j,k,qyn)*sumdrytmp + sumryvtmp*dxinv(2) * &
                        first_deriv_6(q(i,j-3:j+3,k,qyn))
                   ene_c = ry_c*q(i,j,k,qhn) + q(i,j,k,qyn)*sumryvtmp*dxinv(2)* &
                        first_deriv_6(q(i,j-3:j+3,k,qhn))
                   rhs(i,j,k,iene) = rhs(i,j,k,iene) - ene_c
                   rhs(i,j,k,iryn) = rhs(i,j,k,iryn) - ry_c
                end do

             else
    
                n = iias
                qhn = qh1+n-1
                iryn = iry1+n-1

                ene_c = sumdrytmp*q(i,j,k,qhn) + sumryvtmp*dxinv(2)* &
                     first_deriv_6(q(i,j-3:j+3,k,qhn))
                rhs(i,j,k,iene) = rhs(i,j,k,iene) - ene_c
                rhs(i,j,k,iryn) = rhs(i,j,k,iryn) - sumdrytmp
       
             end if
          end do
       end do
    end if

    !
    ! hi-y boundary
    !
    if (physbchi(2)) then
       do k=lo(3),hi(3)
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! use 6th-order stencil for cell i,hi(2)-3,k
          do i=lo(1),hi(1)
             do iface=0,1 
                j = hi(2)-3 + iface

                mmtmp6 = matmul(mu(i,j-3:j+2,k), M6)
                Hcell(iface,imx) = dot_product(mmtmp6, q(i,j-3:j+2,k,qu))
                Hcell(iface,imz) = dot_product(mmtmp6, q(i,j-3:j+2,k,qw))

                mmtmp6 = matmul(vsp(i,j-3:j+2,k), M6)
                Hcell(iface,imy) = dot_product(mmtmp6, q(i,j-3:j+2,k,qv))

                mmtmp6 = matmul(lam(i,j-3:j+2,k), M6)
                M6p = matmul(M6T,  q(i,j-3:j+2,k,qpres))
                Hcell(iface,iene) = dot_product(mmtmp6, q(i,j-3:j+2,k,qtemp))      &
                     &            + dot_product(      dpe(i,j-3:j+2,k), M6p)

                do n = 1, nspecies
                   if (n .eq. iias) cycle  ! inactive speices
                   qxn = qx1+n-1
                   iryn = iry1+n-1
                   M6X = matmul(M6T, q(i,j-3:j+2,k,qxn))
                   Hcell(iface,iene) = Hcell(iface,iene) &
                        +    dot_product(dxe(i,j-3:j+2,k,n), M6X)
                   Hcell(iface,iryn) = dot_product(dpy(i,j-3:j+2,k,n), M6p) &
                        +    dot_product(dxy(i,j-3:j+2,k,n), M6X)
                end do

             end do

             j = hi(2)-3

             do n=2,iene
                rhs(i,j,k,n) = rhs(i,j,k,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(2)
             end do

             sumdrytmp = 0.d0
             do n=iry1,ncons
                if (n.eq.iry_ias) cycle
                sumdrytmp = sumdrytmp + (Hcell(1,n) - Hcell(0,n)) * dx2inv(2)
                rhs(i,j,k,n)  =  rhs(i,j,k,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(2)
             end do

             gradptmp = dxinv(2) * first_deriv_6(q(i,j-3:j+3,k,qpres))

             sumryvtmp = 0.d0
             do n = 1, nspecies
                if (n.eq.iias) cycle
                qxn = qx1+n-1
                sumryvtmp = sumryvtmp + dpy(i,j,k,n)*gradptmp  &
                     + dxy(i,j,k,n)*dxinv(2)*first_deriv_6(q(i,j-3:j+3,k,qxn))
             end do

             if (add_v_correction) then

                do n=1,nspecies
                   qyn = qy1+n-1
                   qhn = qh1+n-1
                   iryn = iry1+n-1
          
                   ry_c = q(i,j,k,qyn)*sumdrytmp + sumryvtmp*dxinv(2) * &
                        first_deriv_6(q(i,j-3:j+3,k,qyn))
                   ene_c = ry_c*q(i,j,k,qhn) + q(i,j,k,qyn)*sumryvtmp*dxinv(2)* &
                        first_deriv_6(q(i,j-3:j+3,k,qhn))
                   rhs(i,j,k,iene) = rhs(i,j,k,iene) - ene_c
                   rhs(i,j,k,iryn) = rhs(i,j,k,iryn) - ry_c
                end do

             else
    
                n = iias
                qhn = qh1+n-1
                iryn = iry1+n-1

                ene_c = sumdrytmp*q(i,j,k,qhn) + sumryvtmp*dxinv(2)* &
                     first_deriv_6(q(i,j-3:j+3,k,qhn))
                rhs(i,j,k,iene) = rhs(i,j,k,iene) - ene_c
                rhs(i,j,k,iryn) = rhs(i,j,k,iryn) - sumdrytmp
       
             end if
          end do

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! use 4th-order stencil for cell i,hi(2)-2,k
          do i=lo(1),hi(1)
             do iface=0,1 
                j = hi(2)-2 + iface

                mmtmp4 = matmul(mu(i,j-2:j+1,k), M4)
                Hcell(iface,imx) = dot_product(mmtmp4, q(i,j-2:j+1,k,qu))
                Hcell(iface,imz) = dot_product(mmtmp4, q(i,j-2:j+1,k,qw))

                mmtmp4 = matmul(vsp(i,j-2:j+1,k), M4)
                Hcell(iface,imy) = dot_product(mmtmp4, q(i,j-2:j+1,k,qv))

                mmtmp4 = matmul(lam(i,j-2:j+1,k), M4)
                M4p = matmul(M4T,  q(i,j-2:j+1,k,qpres))
                Hcell(iface,iene) = dot_product(mmtmp4, q(i,j-2:j+1,k,qtemp)) &
                     &            + dot_product(      dpe(i,j-2:j+1,k), M4p)

                do n = 1, nspecies
                   if (n .eq. iias) cycle  ! inactive speices

                   qxn = qx1+n-1
                   iryn = iry1+n-1
                   M4X = matmul(M4T, q(i,j-2:j+1,k,qxn))
                   Hcell(iface,iene) = Hcell(iface,iene) &
                        +    dot_product(dxe(i,j-2:j+1,k,n), M4X)
                   Hcell(iface,iryn) = dot_product(dpy(i,j-2:j+1,k,n), M4p) &
                        +    dot_product(dxy(i,j-2:j+1,k,n), M4X)
                end do

             end do

             j = hi(2)-2

             do n=2,iene
                rhs(i,j,k,n) = rhs(i,j,k,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(2)
             end do

             sumdrytmp = 0.d0
             do n=iry1,ncons
                if (n.eq.iry_ias) cycle
                sumdrytmp = sumdrytmp + (Hcell(1,n) - Hcell(0,n)) * dx2inv(2)
                rhs(i,j,k,n)  =  rhs(i,j,k,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(2)
             end do

             gradptmp = dxinv(2) * first_deriv_4(q(i,j-2:j+2,k,qpres))

             sumryvtmp = 0.d0
             do n = 1, nspecies
                if (n.eq.iias) cycle
                qxn = qx1+n-1
                sumryvtmp = sumryvtmp + dpy(i,j,k,n)*gradptmp  &
                     + dxy(i,j,k,n)*dxinv(2)*first_deriv_4(q(i,j-2:j+2,k,qxn))
             end do

             if (add_v_correction) then

                do n=1,nspecies
                   qyn = qy1+n-1
                   qhn = qh1+n-1
                   iryn = iry1+n-1
          
                   ry_c = q(i,j,k,qyn)*sumdrytmp + sumryvtmp*dxinv(2) * &
                        first_deriv_4(q(i,j-2:j+2,k,qyn))
                   ene_c = ry_c*q(i,j,k,qhn) + q(i,j,k,qyn)*sumryvtmp*dxinv(2)* &
                        first_deriv_4(q(i,j-2:j+2,k,qhn))
                   rhs(i,j,k,iene) = rhs(i,j,k,iene) - ene_c
                   rhs(i,j,k,iryn) = rhs(i,j,k,iryn) - ry_c
                end do

             else
    
                n = iias
                qhn = qh1+n-1
                iryn = iry1+n-1

                ene_c = sumdrytmp*q(i,j,k,qhn) + sumryvtmp*dxinv(2)* &
                     first_deriv_4(q(i,j-2:j+2,k,qhn))
                rhs(i,j,k,iene) = rhs(i,j,k,iene) - ene_c
                rhs(i,j,k,iryn) = rhs(i,j,k,iryn) - sumdrytmp
       
             end if
          end do

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! use 2nd-order stencil for cell i,hi(2)-1,k
          do i=lo(1),hi(1)
             do iface=0,1 
                j = hi(2)-1 + iface

                mmtmp2 = matmul(mu(i,j-1:j,k), M2)
                Hcell(iface,imx) = dot_product(mmtmp2, q(i,j-1:j,k,qu))
                Hcell(iface,imz) = dot_product(mmtmp2, q(i,j-1:j,k,qw))

                mmtmp2 = matmul(vsp(i,j-1:j,k), M2)
                Hcell(iface,imy) = dot_product(mmtmp2, q(i,j-1:j,k,qv))

                mmtmp2 = matmul(lam(i,j-1:j,k), M2)
                M2p = matmul(M2,  q(i,j-1:j,k,qpres))
                Hcell(iface,iene) = dot_product(mmtmp2, q(i,j-1:j,k,qtemp))      &
                     &            + dot_product(      dpe(i,j-1:j,k), M2p)

                Htot = 0.d0
                do n = 1, nspecies
                   qxn = qx1+n-1
                   qyn = qy1+n-1

                   M2X = matmul(M2, q(i,j-1:j,k,qxn))
                
                   Htmp(n) = dot_product(dpy(i,j-1:j,k,n), M2p) &
                        +    dot_product(dxy(i,j-1:j,k,n), M2X)

                   Hcell(iface,iene) = Hcell(iface,iene) &
                        +    dot_product(dxe(i,j-1:j,k,n), M2X)

                   Htot = Htot + Htmp(n)
                   Ytmp(n) = 0.5d0*(q(i,j-1,k,qyn) + q(i,j,k,qyn))
                end do

                do n = 1, nspecies
                   Hcell(iface,iry1+n-1) = Htmp(n) - Ytmp(n)*Htot
                end do

                do n = 1, nspecies
                   qhn = qh1+n-1
                   hhalf = 0.5d0*(q(i,j-1,k,qhn) + q(i,j,k,qhn))
                   Hcell(iface,iene) =  Hcell(iface,iene) - Ytmp(n)*Htot*hhalf
                end do
             end do

             j = hi(2)-1
             do n=2,ncons
                rhs(i,j,k,n) = rhs(i,j,k,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(2)
             end do
          end do

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          j = hi(2)
          ! use completely right-biased stencil
          do i=lo(1),hi(1)

             mmtmpB = matmul(mu(i,j-3:j,k), BLB)
             rhs(i,j,k,imx) = rhs(i,j,k,imx)+fouhi(2)*dx2inv(2)*dot_product(mmtmpB,q(i,j-3:j,k,qu))
             rhs(i,j,k,imz) = rhs(i,j,k,imz)+fouhi(2)*dx2inv(2)*dot_product(mmtmpB,q(i,j-3:j,k,qw))

             mmtmpB = matmul(vsp(i,j-3:j,k), BLB)
             rhs(i,j,k,imy) = rhs(i,j,k,imy)+finhi(2)*dx2inv(2)*dot_product(mmtmpB,q(i,j-3:j,k,qv))

             mmtmpB = matmul(lam(i,j-3:j,k), BLB)
             BBp = matmul(BLB, q(i,j-3:j,k,qpres))
             rhs(i,j,k,iene) = rhs(i,j,k,iene) + fouhi(2)*dx2inv(2) * &
                  ( dot_product(mmtmpB, q(i,j-3:j,k,qtemp)) &
                  + dot_product(      dpe(i,j-3:j,k), BBp) )
             
             rhstot = 0.d0
             rhsene = 0.d0
             do n = 1, nspecies
                qxn = qx1+n-1
                qyn = qy1+n-1

                BBX = matmul(BLB, q(i,j-3:j,k,qxn))
                
                rhstmp(n) = dot_product(dpy(i,j-3:j,k,n), BBp) &
                     +      dot_product(dxy(i,j-3:j,k,n), BBX)
                
                rhsene = rhsene &
                     +      dot_product(dxe(i,j-3:j,k,n), BBX)

                rhstot = rhstot + rhstmp(n)
                Ytmp(n) = q(i,j,k,qyn)
             end do
             
             do n = 1, nspecies
                rhs(i,j,k,iry1+n-1) =  rhs(i,j,k,iry1+n-1) + &
                     fouhi(2)*dx2inv(2) * (rhstmp(n) - Ytmp(n)*rhstot)
             end do

             do n = 1, nspecies
                qhn = qh1+n-1
                rhsene = rhsene - Ytmp(n) * q(i,j,k,qhn) * rhstot
             end do
             rhs(i,j,k,iene) = rhs(i,j,k,iene) + fouhi(2)*dx2inv(2) * rhsene
          end do
       end do
    end if

    !
    ! lo-z boundary
    !
    if (physbclo(3)) then
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       k = lo(3)
       ! use completely right-biased stencil
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             mmtmpB = matmul(mu(i,j,k:k+3), BRB)
             rhs(i,j,k,imx) = rhs(i,j,k,imx)+foulo(3)*dx2inv(3)*dot_product(mmtmpB,q(i,j,k:k+3,qu))
             rhs(i,j,k,imy) = rhs(i,j,k,imy)+foulo(3)*dx2inv(3)*dot_product(mmtmpB,q(i,j,k:k+3,qv))

             mmtmpB = matmul(vsp(i,j,k:k+3), BRB)
             rhs(i,j,k,imz) = rhs(i,j,k,imz)+finlo(3)*dx2inv(3)*dot_product(mmtmpB,q(i,j,k:k+3,qw) )

             mmtmpB = matmul(lam(i,j,k:k+3), BRB)
             BBp = matmul(BRB, q(i,j,k:k+3,qpres))
             rhs(i,j,k,iene) = rhs(i,j,k,iene) + foulo(3)*dx2inv(3) * &
                  ( dot_product(mmtmpB, q(i,j,k:k+3,qtemp)) &
                  + dot_product(      dpe(i,j,k:k+3), BBp) )
             
             rhstot = 0.d0
             rhsene = 0.d0
             do n = 1, nspecies
                qxn = qx1+n-1
                qyn = qy1+n-1

                BBX = matmul(BRB, q(i,j,k:k+3,qxn))
                
                rhstmp(n) = dot_product(dpy(i,j,k:k+3,n), BBp) &
                     +      dot_product(dxy(i,j,k:k+3,n), BBX)
                
                rhsene = rhsene &
                     +      dot_product(dxe(i,j,k:k+3,n), BBX)

                rhstot = rhstot + rhstmp(n)
                Ytmp(n) = q(i,j,k,qyn)
             end do
             
             do n = 1, nspecies
                rhs(i,j,k,iry1+n-1) =  rhs(i,j,k,iry1+n-1) + &
                     foulo(3)*dx2inv(3) * (rhstmp(n) - Ytmp(n)*rhstot)
             end do

             do n = 1, nspecies
                qhn = qh1+n-1
                rhsene = rhsene - Ytmp(n) * q(i,j,k,qhn) * rhstot
             end do
             rhs(i,j,k,iene) = rhs(i,j,k,iene) + foulo(3)*dx2inv(3) * rhsene
             
          end do
       end do

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! use 2nd-order stencil for cell i,j,lo(3)+1
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             do iface=0,1 
                k = lo(3)+1 + iface

                mmtmp2 = matmul(mu(i,j,k-1:k), M2)
                Hcell(iface,imx) = dot_product(mmtmp2, q(i,j,k-1:k,qu))
                Hcell(iface,imy) = dot_product(mmtmp2, q(i,j,k-1:k,qv))

                mmtmp2 = matmul(vsp(i,j,k-1:k), M2)
                Hcell(iface,imz) = dot_product(mmtmp2, q(i,j,k-1:k,qw))

                mmtmp2 = matmul(lam(i,j,k-1:k), M2)
                M2p = matmul(M2,  q(i,j,k-1:k,qpres))
                Hcell(iface,iene) = dot_product(mmtmp2, q(i,j,k-1:k,qtemp)) &
                     &            + dot_product(       dpe(i,j,k-1:k), M2p)

                Htot = 0.d0
                do n = 1, nspecies
                   qxn = qx1+n-1
                   qyn = qy1+n-1

                   M2X = matmul(M2, q(i,j,k-1:k,qxn))
                
                   Htmp(n) = dot_product(dpy(i,j,k-1:k,n), M2p) &
                        +    dot_product(dxy(i,j,k-1:k,n), M2X)

                   Hcell(iface,iene) = Hcell(iface,iene) &
                        +    dot_product(dxe(i,j,k-1:k,n), M2X)

                   Htot = Htot + Htmp(n)
                   Ytmp(n) = 0.5d0*(q(i,j,k-1,qyn) + q(i,j,k,qyn))
                end do

                do n = 1, nspecies
                   Hcell(iface,iry1+n-1) = Htmp(n) - Ytmp(n)*Htot
                end do

                do n = 1, nspecies
                   qhn = qh1+n-1
                   hhalf = 0.5d0*(q(i,j,k-1,qhn) + q(i,j,k,qhn))
                   Hcell(iface,iene) =  Hcell(iface,iene) - Ytmp(n)*Htot*hhalf
                end do
             end do

             k = lo(3)+1
             do n=2,ncons
                rhs(i,j,k,n) = rhs(i,j,k,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(3)
             end do
          end do
       end do

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! use 4th-order stencil for cell i,j,lo(3)+2
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             do iface=0,1 
                k = lo(3)+2 + iface

                mmtmp4 = matmul(mu(i,j,k-2:k+1), M4)
                Hcell(iface,imx) = dot_product(mmtmp4, q(i,j,k-2:k+1,qu))
                Hcell(iface,imy) = dot_product(mmtmp4, q(i,j,k-2:k+1,qv))

                mmtmp4 = matmul(vsp(i,j,k-2:k+1), M4)
                Hcell(iface,imz) = dot_product(mmtmp4, q(i,j,k-2:k+1,qw))

                mmtmp4 = matmul(lam(i,j,k-2:k+1), M4)
                M4p = matmul(M4T,  q(i,j,k-2:k+1,qpres))
                Hcell(iface,iene) = dot_product(mmtmp4, q(i,j,k-2:k+1,qtemp)) &
                     &            + dot_product(      dpe(i,j,k-2:k+1), M4p)

                do n = 1, nspecies
                   if (n .eq. iias) cycle  ! inactive speices

                   qxn = qx1+n-1
                   iryn = iry1+n-1
                   M4X = matmul(M4T, q(i,j,k-2:k+1,qxn))
                   Hcell(iface,iene) = Hcell(iface,iene) &
                        +    dot_product(dxe(i,j,k-2:k+1,n), M4X)
                   Hcell(iface,iryn) = dot_product(dpy(i,j,k-2:k+1,n), M4p) &
                        +    dot_product(dxy(i,j,k-2:k+1,n), M4X)
                end do

             end do

             k = lo(3)+2

             do n=2,iene
                rhs(i,j,k,n) = rhs(i,j,k,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(3)
             end do

             sumdrytmp = 0.d0
             do n=iry1,ncons
                if (n.eq.iry_ias) cycle
                sumdrytmp = sumdrytmp + (Hcell(1,n) - Hcell(0,n)) * dx2inv(3)
                rhs(i,j,k,n)  =  rhs(i,j,k,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(3)
             end do

             gradptmp = dxinv(3) * first_deriv_4(q(i,j,k-2:k+2,qpres))

             sumryvtmp = 0.d0
             do n = 1, nspecies
                if (n.eq.iias) cycle
                qxn = qx1+n-1
                sumryvtmp = sumryvtmp + dpy(i,j,k,n)*gradptmp  &
                     + dxy(i,j,k,n)*dxinv(3)*first_deriv_4(q(i,j,k-2:k+2,qxn))
             end do

             if (add_v_correction) then

                do n=1,nspecies
                   qyn = qy1+n-1
                   qhn = qh1+n-1
                   iryn = iry1+n-1
          
                   ry_c = q(i,j,k,qyn)*sumdrytmp + sumryvtmp*dxinv(3) * &
                        first_deriv_4(q(i,j,k-2:k+2,qyn))
                   ene_c = ry_c*q(i,j,k,qhn) + q(i,j,k,qyn)*sumryvtmp*dxinv(3)* &
                        first_deriv_4(q(i,j,k-2:k+2,qhn))
                   rhs(i,j,k,iene) = rhs(i,j,k,iene) - ene_c
                   rhs(i,j,k,iryn) = rhs(i,j,k,iryn) - ry_c
                end do

             else
    
                n = iias
                qhn = qh1+n-1
                iryn = iry1+n-1

                ene_c = sumdrytmp*q(i,j,k,qhn) + sumryvtmp*dxinv(3)* &
                     first_deriv_4(q(i,j,k-2:k+2,qhn))
                rhs(i,j,k,iene) = rhs(i,j,k,iene) - ene_c
                rhs(i,j,k,iryn) = rhs(i,j,k,iryn) - sumdrytmp
       
             end if
          end do
       end do

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! use 6th-order stencil for cell i,j,lo(3)+3
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             do iface=0,1 
                k = lo(3)+3 + iface

                mmtmp6 = matmul(mu(i,j,k-3:k+2), M6)
                Hcell(iface,imx) = dot_product(mmtmp6, q(i,j,k-3:k+2,qu))
                Hcell(iface,imy) = dot_product(mmtmp6, q(i,j,k-3:k+2,qv))

                mmtmp6 = matmul(vsp(i,j,k-3:k+2), M6)
                Hcell(iface,imz) = dot_product(mmtmp6, q(i,j,k-3:k+2,qw))

                mmtmp6 = matmul(lam(i,j,k-3:k+2), M6)
                M6p = matmul(M6T,  q(i,j,k-3:k+2,qpres))
                Hcell(iface,iene) = dot_product(mmtmp6, q(i,j,k-3:k+2,qtemp)) &
                     &            + dot_product(      dpe(i,j,k-3:k+2), M6p)

                do n = 1, nspecies
                   if (n .eq. iias) cycle  ! inactive speices
                   qxn = qx1+n-1
                   iryn = iry1+n-1
                   M6X = matmul(M6T, q(i,j,k-3:k+2,qxn))
                   Hcell(iface,iene) = Hcell(iface,iene) &
                        +    dot_product(dxe(i,j,k-3:k+2,n), M6X)
                   Hcell(iface,iryn) = dot_product(dpy(i,j,k-3:k+2,n), M6p) &
                        +    dot_product(dxy(i,j,k-3:k+2,n), M6X)
                end do

             end do

             k = lo(3)+3

             do n=2,iene
                rhs(i,j,k,n) = rhs(i,j,k,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(3)
             end do

             sumdrytmp = 0.d0
             do n=iry1,ncons
                if (n.eq.iry_ias) cycle
                sumdrytmp = sumdrytmp + (Hcell(1,n) - Hcell(0,n)) * dx2inv(3)
                rhs(i,j,k,n)  =  rhs(i,j,k,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(3)
             end do

             gradptmp = dxinv(3) * first_deriv_6(q(i,j,k-3:k+3,qpres))

             sumryvtmp = 0.d0
             do n = 1, nspecies
                if (n.eq.iias) cycle
                qxn = qx1+n-1
                sumryvtmp = sumryvtmp + dpy(i,j,k,n)*gradptmp  &
                     + dxy(i,j,k,n)*dxinv(3)*first_deriv_6(q(i,j,k-3:k+3,qxn))
             end do

             if (add_v_correction) then

                do n=1,nspecies
                   qyn = qy1+n-1
                   qhn = qh1+n-1
                   iryn = iry1+n-1
          
                   ry_c = q(i,j,k,qyn)*sumdrytmp + sumryvtmp*dxinv(3) * &
                        first_deriv_6(q(i,j,k-3:k+3,qyn))
                   ene_c = ry_c*q(i,j,k,qhn) + q(i,j,k,qyn)*sumryvtmp*dxinv(3)* &
                        first_deriv_6(q(i,j,k-3:k+3,qhn))
                   rhs(i,j,k,iene) = rhs(i,j,k,iene) - ene_c
                   rhs(i,j,k,iryn) = rhs(i,j,k,iryn) - ry_c
                end do

             else
    
                n = iias
                qhn = qh1+n-1
                iryn = iry1+n-1

                ene_c = sumdrytmp*q(i,j,k,qhn) + sumryvtmp*dxinv(3)* &
                     first_deriv_6(q(i,j,k-3:k+3,qhn))
                rhs(i,j,k,iene) = rhs(i,j,k,iene) - ene_c
                rhs(i,j,k,iryn) = rhs(i,j,k,iryn) - sumdrytmp
       
             end if
          end do
       end do
    end if

    !
    ! hi-z boundary
    !
    if (physbchi(3)) then
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! use 6th-order stencil for cell i,j,hi(3)-3
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             do iface=0,1 
                k = hi(3)-3 + iface

                mmtmp6 = matmul(mu(i,j,k-3:k+2), M6)
                Hcell(iface,imx) = dot_product(mmtmp6, q(i,j,k-3:k+2,qu))
                Hcell(iface,imy) = dot_product(mmtmp6, q(i,j,k-3:k+2,qv))

                mmtmp6 = matmul(vsp(i,j,k-3:k+2), M6)
                Hcell(iface,imz) = dot_product(mmtmp6, q(i,j,k-3:k+2,qw))

                mmtmp6 = matmul(lam(i,j,k-3:k+2), M6)
                M6p = matmul(M6T,  q(i,j,k-3:k+2,qpres))
                Hcell(iface,iene) = dot_product(mmtmp6, q(i,j,k-3:k+2,qtemp)) &
                     &            + dot_product(      dpe(i,j,k-3:k+2), M6p)

                do n = 1, nspecies
                   if (n .eq. iias) cycle  ! inactive speices
                   qxn = qx1+n-1
                   iryn = iry1+n-1
                   M6X = matmul(M6T, q(i,j,k-3:k+2,qxn))
                   Hcell(iface,iene) = Hcell(iface,iene) &
                        +    dot_product(dxe(i,j,k-3:k+2,n), M6X)
                   Hcell(iface,iryn) = dot_product(dpy(i,j,k-3:k+2,n), M6p) &
                        +    dot_product(dxy(i,j,k-3:k+2,n), M6X)
                end do

             end do

             k = hi(3)-3

             do n=2,iene
                rhs(i,j,k,n) = rhs(i,j,k,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(3)
             end do

             sumdrytmp = 0.d0
             do n=iry1,ncons
                if (n.eq.iry_ias) cycle
                sumdrytmp = sumdrytmp + (Hcell(1,n) - Hcell(0,n)) * dx2inv(3)
                rhs(i,j,k,n)  =  rhs(i,j,k,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(3)
             end do

             gradptmp = dxinv(3) * first_deriv_6(q(i,j,k-3:k+3,qpres))

             sumryvtmp = 0.d0
             do n = 1, nspecies
                if (n.eq.iias) cycle
                qxn = qx1+n-1
                sumryvtmp = sumryvtmp + dpy(i,j,k,n)*gradptmp  &
                     + dxy(i,j,k,n)*dxinv(3)*first_deriv_6(q(i,j,k-3:k+3,qxn))
             end do

             if (add_v_correction) then

                do n=1,nspecies
                   qyn = qy1+n-1
                   qhn = qh1+n-1
                   iryn = iry1+n-1
          
                   ry_c = q(i,j,k,qyn)*sumdrytmp + sumryvtmp*dxinv(3) * &
                        first_deriv_6(q(i,j,k-3:k+3,qyn))
                   ene_c = ry_c*q(i,j,k,qhn) + q(i,j,k,qyn)*sumryvtmp*dxinv(3)* &
                        first_deriv_6(q(i,j,k-3:k+3,qhn))
                   rhs(i,j,k,iene) = rhs(i,j,k,iene) - ene_c
                   rhs(i,j,k,iryn) = rhs(i,j,k,iryn) - ry_c
                end do

             else
    
                n = iias
                qhn = qh1+n-1
                iryn = iry1+n-1

                ene_c = sumdrytmp*q(i,j,k,qhn) + sumryvtmp*dxinv(3)* &
                     first_deriv_6(q(i,j,k-3:k+3,qhn))
                rhs(i,j,k,iene) = rhs(i,j,k,iene) - ene_c
                rhs(i,j,k,iryn) = rhs(i,j,k,iryn) - sumdrytmp
       
             end if
          end do
       end do

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! use 4th-order stencil for cell i,j,hi(3)-2
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             do iface=0,1 
                k = hi(3)-2 + iface

                mmtmp4 = matmul(mu(i,j,k-2:k+1), M4)
                Hcell(iface,imx) = dot_product(mmtmp4, q(i,j,k-2:k+1,qu))
                Hcell(iface,imy) = dot_product(mmtmp4, q(i,j,k-2:k+1,qv))

                mmtmp4 = matmul(vsp(i,j,k-2:k+1), M4)
                Hcell(iface,imz) = dot_product(mmtmp4, q(i,j,k-2:k+1,qw))

                mmtmp4 = matmul(lam(i,j,k-2:k+1), M4)
                M4p = matmul(M4T,  q(i,j,k-2:k+1,qpres))
                Hcell(iface,iene) = dot_product(mmtmp4, q(i,j,k-2:k+1,qtemp)) &
                     &            + dot_product(      dpe(i,j,k-2:k+1), M4p)

                do n = 1, nspecies
                   if (n .eq. iias) cycle  ! inactive speices

                   qxn = qx1+n-1
                   iryn = iry1+n-1
                   M4X = matmul(M4T, q(i,j,k-2:k+1,qxn))
                   Hcell(iface,iene) = Hcell(iface,iene) &
                        +    dot_product(dxe(i,j,k-2:k+1,n), M4X)
                   Hcell(iface,iryn) = dot_product(dpy(i,j,k-2:k+1,n), M4p) &
                        +    dot_product(dxy(i,j,k-2:k+1,n), M4X)
                end do

             end do

             k = hi(3)-2

             do n=2,iene
                rhs(i,j,k,n) = rhs(i,j,k,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(3)
             end do

             sumdrytmp = 0.d0
             do n=iry1,ncons
                if (n.eq.iry_ias) cycle
                sumdrytmp = sumdrytmp + (Hcell(1,n) - Hcell(0,n)) * dx2inv(3)
                rhs(i,j,k,n)  =  rhs(i,j,k,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(3)
             end do

             gradptmp = dxinv(3) * first_deriv_4(q(i,j,k-2:k+2,qpres))

             sumryvtmp = 0.d0
             do n = 1, nspecies
                if (n.eq.iias) cycle
                qxn = qx1+n-1
                sumryvtmp = sumryvtmp + dpy(i,j,k,n)*gradptmp  &
                     + dxy(i,j,k,n)*dxinv(3)*first_deriv_4(q(i,j,k-2:k+2,qxn))
             end do

             if (add_v_correction) then

                do n=1,nspecies
                   qyn = qy1+n-1
                   qhn = qh1+n-1
                   iryn = iry1+n-1
          
                   ry_c = q(i,j,k,qyn)*sumdrytmp + sumryvtmp*dxinv(3) * &
                        first_deriv_4(q(i,j,k-2:k+2,qyn))
                   ene_c = ry_c*q(i,j,k,qhn) + q(i,j,k,qyn)*sumryvtmp*dxinv(3)* &
                        first_deriv_4(q(i,j,k-2:k+2,qhn))
                   rhs(i,j,k,iene) = rhs(i,j,k,iene) - ene_c
                   rhs(i,j,k,iryn) = rhs(i,j,k,iryn) - ry_c
                end do

             else
    
                n = iias
                qhn = qh1+n-1
                iryn = iry1+n-1

                ene_c = sumdrytmp*q(i,j,k,qhn) + sumryvtmp*dxinv(3)* &
                     first_deriv_4(q(i,j,k-2:k+2,qhn))
                rhs(i,j,k,iene) = rhs(i,j,k,iene) - ene_c
                rhs(i,j,k,iryn) = rhs(i,j,k,iryn) - sumdrytmp
       
             end if
          end do
       end do

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! use 2nd-order stencil for cell i,j,hi(3)-1
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             do iface=0,1 
                k = hi(3)-1 + iface

                mmtmp2 = matmul(mu(i,j,k-1:k), M2)
                Hcell(iface,imx) = dot_product(mmtmp2, q(i,j,k-1:k,qu))
                Hcell(iface,imy) = dot_product(mmtmp2, q(i,j,k-1:k,qv))

                mmtmp2 = matmul(vsp(i,j,k-1:k), M2)
                Hcell(iface,imz) = dot_product(mmtmp2, q(i,j,k-1:k,qw))

                mmtmp2 = matmul(lam(i,j,k-1:k), M2)
                M2p = matmul(M2,  q(i,j,k-1:k,qpres))
                Hcell(iface,iene) = dot_product(mmtmp2, q(i,j,k-1:k,qtemp)) &
                     &            + dot_product(      dpe(i,j,k-1:k), M2p)

                Htot = 0.d0
                do n = 1, nspecies
                   qxn = qx1+n-1
                   qyn = qy1+n-1

                   M2X = matmul(M2, q(i,j,k-1:k,qxn))
                
                   Htmp(n) = dot_product(dpy(i,j,k-1:k,n), M2p) &
                        +    dot_product(dxy(i,j,k-1:k,n), M2X)

                   Hcell(iface,iene) = Hcell(iface,iene) &
                        +    dot_product(dxe(i,j,k-1:k,n), M2X)

                   Htot = Htot + Htmp(n)
                   Ytmp(n) = 0.5d0*(q(i,j,k-1,qyn) + q(i,j,k,qyn))
                end do

                do n = 1, nspecies
                   Hcell(iface,iry1+n-1) = Htmp(n) - Ytmp(n)*Htot
                end do

                do n = 1, nspecies
                   qhn = qh1+n-1
                   hhalf = 0.5d0*(q(i,j,k-1,qhn) + q(i,j,k,qhn))
                   Hcell(iface,iene) =  Hcell(iface,iene) - Ytmp(n)*Htot*hhalf
                end do
             end do

             k = hi(3)-1
             do n=2,ncons
                rhs(i,j,k,n) = rhs(i,j,k,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(3)
             end do
          end do
       end do

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       k = hi(3)
       ! use completely right-biased stencil
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             mmtmpB = matmul(mu(i,j,k-3:k), BLB)
             rhs(i,j,k,imx) = rhs(i,j,k,imx)+fouhi(3)*dx2inv(3)*dot_product(mmtmpB,q(i,j,k-3:k,qu))
             rhs(i,j,k,imy) = rhs(i,j,k,imy)+fouhi(3)*dx2inv(3)*dot_product(mmtmpB,q(i,j,k-3:k,qv))

             mmtmpB = matmul(vsp(i,j,k-3:k), BLB)
             rhs(i,j,k,imz) = rhs(i,j,k,imz)+finhi(3)*dx2inv(3)*dot_product(mmtmpB,q(i,j,k-3:k,qw))

             mmtmpB = matmul(lam(i,j,k-3:k), BLB)
             BBp = matmul(BLB, q(i,j,k-3:k,qpres))
             rhs(i,j,k,iene) = rhs(i,j,k,iene) + fouhi(3)*dx2inv(3) * &
                  ( dot_product(mmtmpB, q(i,j,k-3:k,qtemp)) &
                  + dot_product(      dpe(i,j,k-3:k), BBp) )
             
             rhstot = 0.d0
             rhsene = 0.d0
             do n = 1, nspecies
                qxn = qx1+n-1
                qyn = qy1+n-1

                BBX = matmul(BLB, q(i,j,k-3:k,qxn))
                
                rhstmp(n) = dot_product(dpy(i,j,k-3:k,n), BBp) &
                     +      dot_product(dxy(i,j,k-3:k,n), BBX)
                
                rhsene = rhsene &
                     +      dot_product(dxe(i,j,k-3:k,n), BBX)

                rhstot = rhstot + rhstmp(n)
                Ytmp(n) = q(i,j,k,qyn)
             end do
             
             do n = 1, nspecies
                rhs(i,j,k,iry1+n-1) =  rhs(i,j,k,iry1+n-1) + &
                     fouhi(3)*dx2inv(3) * (rhstmp(n) - Ytmp(n)*rhstot)
             end do

             do n = 1, nspecies
                qhn = qh1+n-1
                rhsene = rhsene - Ytmp(n) * q(i,j,k,qhn) * rhstot
             end do
             rhs(i,j,k,iene) = rhs(i,j,k,iene) + fouhi(3)*dx2inv(3) * rhsene
             
          end do
       end do
    end if

    !
    ! add kinetic energy
    !
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

    deallocate(Hg,dpy,dxe,dpe,vsp,M8p)
    deallocate(sumdrY,sumryv,gradp)

  end subroutine diffterm_2


  subroutine chemterm_3d(lo,hi,q,qlo,qhi,up,uplo,uphi,upc,upclo,upchi,dt)
    use probin_module, only : use_vode
    use vode_module, only : verbose, itol, rtol, atol, vode_MF=>MF, &
         voderwork, vodeiwork, lvoderwork, lvodeiwork, voderpar, vodeipar

    double precision, intent(in) :: dt
    integer,         intent(in):: lo(3),hi(3),qlo(3),qhi(3),uplo(3),uphi(3),upclo(3),upchi(3)
    double precision,intent(in):: q  (  qlo(1):  qhi(1),  qlo(2):  qhi(2),  qlo(3):  qhi(3),nprim)
    double precision           :: up ( uplo(1): uphi(1), uplo(2): uphi(2), uplo(3): uphi(3),ncons)
    double precision           :: upc(upclo(1):upchi(1),upclo(2):upchi(2),upclo(3):upchi(3),nspecies)

    integer :: iwrk, i,j,k,n,np
    double precision :: Yt(lo(1):hi(1),nspecies), wdot(lo(1):hi(1),nspecies), rwrk
    double precision :: YTvode(nspecies+1), time, dtinv

    external f_jac, f_rhs, dvode

    ! vode stuff
    integer, parameter :: itask=1, iopt=1
    integer :: MF, istate, ifail

    if (use_vode) then

       dtinv = 1.d0/dt

       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)

                voderpar(1) = q(i,j,k,qrho)

                YTvode(1:nspecies) = q(i,j,k,qy1:qy1+nspecies-1)
                YTvode(nspecies+1) = q(i,j,k,qtemp)
                
                istate = 1
                
                time = 0.d0
                
                MF = vode_MF  ! vode might change its sign!
                call dvode(f_rhs, nspecies+1, YTvode, time, dt, itol, rtol, atol, itask, &
                     istate, iopt, voderwork, lvoderwork, vodeiwork, lvodeiwork, &
                     f_jac, MF, voderpar, vodeipar)
                
                if (verbose .ge. 1) then
                   write(6,*) '......dvode done:'
                   write(6,*) ' last successful step size = ',voderwork(11)
                   write(6,*) '          next step to try = ',voderwork(12)
                   write(6,*) '   integrated time reached = ',voderwork(13)
                   write(6,*) '      number of time steps = ',vodeiwork(11)
                   write(6,*) '              number of fs = ',vodeiwork(12)
                   write(6,*) '              number of Js = ',vodeiwork(13)
                   write(6,*) '    method order last used = ',vodeiwork(14)
                   write(6,*) '   method order to be used = ',vodeiwork(15)
                   write(6,*) '            number of LUDs = ',vodeiwork(19)
                   write(6,*) ' number of Newton iterations ',vodeiwork(20)
                   write(6,*) ' number of Newton failures = ',vodeiwork(21)
                   if (ISTATE.eq.-4 .or. ISTATE.eq.-5) then
                      ifail = vodeiwork(16)
                      if (ifail .eq. nspecies+1) then
                         write(6,*) '   T has the largest error'
                      else
                         write(6,*) '   spec with largest error is No. ', ifail
                      end if
                      call flush(6)
                   end if
                end if
                
                if (istate < 0) then
                   print *, 'chemsolv: VODE failed'
                   print *, 'istate = ', istate, ' time =', time
                   call bl_error("ERROR in burn: VODE failed")
                end if
                
                do n=1, nspecies
                   upc(i,j,k,n) = dtinv*q(i,j,k,qrho)*(YTvode(n)-q(i,j,k,qy1+n-1))
                   up (i,j,k,iry1+n-1) = upc(i,j,k,n)
                end do
                
             end do
          end do
       end do

    else

       np = hi(1) - lo(1) + 1
       
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             
             do n=1, nspecies
                do i=lo(1),hi(1)
                   Yt(i,n) = q(i,j,k,qy1+n-1)
                end do
             end do
             
             call vckwyr(np, q(lo(1),j,k,qrho), q(lo(1),j,k,qtemp), Yt, iwrk, rwrk, wdot)
             
             do n=1, nspecies
                do i=lo(1),hi(1)
                   upc(i,j,k,n) = wdot(i,n) * molecular_weight(n)
                   up (i,j,k,iry1+n-1) = upc(i,j,k,n)
                end do
             end do
             
          end do
       end do

    end if

  end subroutine chemterm_3d


  subroutine comp_courno_3d(lo,hi,dx,Q,qlo,qhi,courno)
    integer, intent(in) :: lo(3), hi(3), qlo(3), qhi(3)
    double precision, intent(in) :: dx(3)
    double precision, intent(in) :: q(qlo(1):qhi(1),qlo(2):qhi(2),qlo(3):qhi(3),nprim)
    double precision, intent(inout) :: courno

    integer :: i,j,k, iwrk
    double precision :: dxinv(3), c, rwrk, Cv, Cp
    double precision :: Tt, X(nspecies), gamma
    double precision :: courx, coury, courz

    do i=1,3
       dxinv(i) = 1.0d0 / dx(i)
    end do

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             Tt = q(i,j,k,qtemp)
             X  = q(i,j,k,qx1:qx1+nspecies-1)
             call ckcvbl(Tt, X, iwrk, rwrk, Cv)
             Cp = Cv + Ru
             gamma = Cp / Cv
             c = sqrt(gamma*q(i,j,k,qpres)/q(i,j,k,qrho))

             courx = (c+abs(q(i,j,k,qu))) * dxinv(1)
             coury = (c+abs(q(i,j,k,qv))) * dxinv(2)
             courz = (c+abs(q(i,j,k,qw))) * dxinv(3)
             
             courno = max( courx, coury, courz , courno )

          end do
       end do
    end do

  end subroutine comp_courno_3d

end module kernels_module
