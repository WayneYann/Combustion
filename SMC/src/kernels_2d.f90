module kernels_2d_module
  use bc_module
  use chemistry_module, only : nspecies, molecular_weight, Ru
  use derivative_stencil_module, only : stencil_ng, first_deriv_8, first_deriv_6, &
       first_deriv_4, first_deriv_l3, first_deriv_r3, first_deriv_rb, first_deriv_lb, &
       M8, M8T, M6, M6T, M4, M4T, M2, BRB, BLB, D8, D6, D4
  use variables_module
  implicit none

  private

  public :: hypterm_2d, narrow_diffterm_2d, chemterm_2d, comp_courno_2d

contains

  subroutine hypterm_2d (lo,hi,dx,cons,clo,chi,q,qlo,qhi,rhs_g,rlo,rhi,dlo_g,dhi_g,bclo,bchi)

    integer,         intent(in):: dlo_g(2),dhi_g(2),bclo(2),bchi(2)
    integer,         intent(in):: lo(2),hi(2),clo(2),chi(2),qlo(2),qhi(2),rlo(2),rhi(2)
    double precision,intent(in):: dx(2)
    double precision,intent(in):: cons(clo(1):chi(1),clo(2):chi(2),ncons)
    double precision,intent(in)::    q(qlo(1):qhi(1),qlo(2):qhi(2),nprim)
    double precision           ::rhs_g(rlo(1):rhi(1),rlo(2):rhi(2),ncons)

    integer          :: i,j,n
    double precision :: dxinv(2)
    double precision :: un(-4:4)
    integer :: slo(2), shi(2), dlo(2), dhi(2)
    
    double precision, allocatable :: tmpx(:), tmpy(:,:)
    double precision, allocatable :: rhs(:,:,:)

    logical :: physbclo(2), physbchi(2)

    ! Only the region bounded by [dlo_g,dhi_g] contains good data.
    ! [slo,shi] will be safe for 8th-order stencil
    do i=1,2
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

    do i=1,2
       dxinv(i) = 1.0d0 / dx(i)
    end do

    allocate(tmpx(dlo(1):dhi(1)))
    allocate(tmpy(lo(1) : hi(1),dlo(2):dhi(2)))

    allocate(rhs(lo(1):hi(1),lo(2):hi(2),ncons))
    rhs = 0.d0

    ! ------- BEGIN x-direction -------

    do j=lo(2),hi(2)
       
       do i=slo(1),shi(1)
          rhs(i,j,irho) = rhs(i,j,irho) - dxinv(1) * &
               first_deriv_8( cons(i-4:i+4,j,imx) ) 
       end do

       do i=dlo(1),dhi(1)
          tmpx(i) = cons(i,j,imx)*q(i,j,qu)+q(i,j,qpres)
       end do
       do i=slo(1),shi(1)
          rhs(i,j,imx) = rhs(i,j,imx) - dxinv(1) * first_deriv_8(tmpx(i-4:i+4))
       end do
       
       do i=dlo(1),dhi(1)
          tmpx(i) = cons(i,j,imy)*q(i,j,qu)
       end do
       do i=slo(1),shi(1)
          rhs(i,j,imy) = rhs(i,j,imy) - dxinv(1) * first_deriv_8(tmpx(i-4:i+4))
       end do
       
       do i=dlo(1),dhi(1)
          tmpx(i) = (cons(i,j,iene)+q(i,j,qpres))*q(i,j,qu)
       end do
       do i=slo(1),shi(1)
          rhs(i,j,iene) = rhs(i,j,iene) - dxinv(1) * first_deriv_8(tmpx(i-4:i+4))
       end do
       
    enddo

    do n = iry1, iry1+nspecies-1
       do j=lo(2),hi(2)
    
          do i=dlo(1),dhi(1)
             tmpx(i) = cons(i,j,n)*q(i,j,qu)
          end do
          do i=slo(1),shi(1)
             rhs(i,j,n) = rhs(i,j,n) - dxinv(1) * first_deriv_8(tmpx(i-4:i+4))
          end do
          
       enddo
    enddo

    ! ------- END x-direction -------

    ! ------- BEGIN y-direction -------

    do j=slo(2),shi(2)
       do i=lo(1),hi(1)
          rhs(i,j,irho)=rhs(i,j,irho) - dxinv(2) * &
               first_deriv_8( cons(i,j-4:j+4,imy) )
       end do
    end do

    do j=dlo(2),dhi(2)
       do i=lo(1),hi(1)
          tmpy(i,j) = cons(i,j,imx)*q(i,j,qv)
       end do
    end do
    do j=slo(2),shi(2)
       do i=lo(1),hi(1)
          rhs(i,j,imx)=rhs(i,j,imx) - dxinv(2) * first_deriv_8(tmpy(i,j-4:j+4))
       enddo
    enddo
    
    do j=dlo(2),dhi(2)
       do i=lo(1),hi(1)
          tmpy(i,j) = cons(i,j,imy)*q(i,j,qv)+q(i,j,qpres)
       end do
    end do
    do j=slo(2),shi(2)
       do i=lo(1),hi(1)
          rhs(i,j,imy)=rhs(i,j,imy) - dxinv(2) * first_deriv_8(tmpy(i,j-4:j+4))
       enddo
    enddo
    
    do j=dlo(2),dhi(2)
       do i=lo(1),hi(1)
          tmpy(i,j) = (cons(i,j,iene)+q(i,j,qpres))*q(i,j,qv)
       end do
    end do
    do j=slo(2),shi(2)
       do i=lo(1),hi(1)
          rhs(i,j,iene)=rhs(i,j,iene) - dxinv(2) * first_deriv_8(tmpy(i,j-4:j+4))
       enddo
    enddo
    
    do n = iry1, iry1+nspecies-1
       do j=dlo(2),dhi(2)
          do i=lo(1),hi(1)
             tmpy(i,j) = cons(i,j,n)*q(i,j,qv)
          end do
       end do
       do j=slo(2),shi(2)
          do i=lo(1),hi(1)
             rhs(i,j,n) = rhs(i,j,n) - dxinv(2) * first_deriv_8(tmpy(i,j-4:j+4))
          end do
       enddo
    enddo
    
    ! ------- END y-direction -------

    ! ----------------- boundary -----------------------

    !
    ! ----- lo-x boundary -----
    !
    if (physbclo(1)) then 
       do j=lo(2),hi(2)
          
          ! if (bclo(1) .eq. WALL???) then
          !    i = lo(1)
          !    ! use completely right-biased stencil
          ! end if
          
          i = lo(1)+1
          ! use 3rd-order slightly right-biased stencil
          
          un(-1:2) = q(i-1:i+2,j,qu)
          
          rhs(i,j,irho) = rhs(i,j,irho) - dxinv(1) * &
               first_deriv_r3( cons(i-1:i+2,j,imx) ) 
          
          rhs(i,j,imx) = rhs(i,j,imx) - dxinv(1) * &
               first_deriv_r3( cons(i-1:i+2,j,imx)*un(-1:2)+q(i-1:i+2,j,qpres) )
          
          rhs(i,j,imy) = rhs(i,j,imy) - dxinv(1) * &
               first_deriv_r3( cons(i-1:i+2,j,imy)*un(-1:2) ) 
          
          rhs(i,j,iene) = rhs(i,j,iene) - dxinv(1) * &
               first_deriv_r3( (cons(i-1:i+2,j,iene)+q(i-1:i+2,j,qpres))*un(-1:2) )

          i = lo(1)+2
          ! use 4th-order stencil

          un(-2:2) = q(i-2:i+2,j,qu)
          
          rhs(i,j,irho) = rhs(i,j,irho) - dxinv(1) * &
               first_deriv_4( cons(i-2:i+2,j,imx) ) 
          
          rhs(i,j,imx) = rhs(i,j,imx) - dxinv(1) * &
               first_deriv_4( cons(i-2:i+2,j,imx)*un(-2:2)+q(i-2:i+2,j,qpres) )
          
          rhs(i,j,imy) = rhs(i,j,imy) - dxinv(1) * &
               first_deriv_4( cons(i-2:i+2,j,imy)*un(-2:2) ) 
          
          rhs(i,j,iene) = rhs(i,j,iene) - dxinv(1) * &
               first_deriv_4( (cons(i-2:i+2,j,iene)+q(i-2:i+2,j,qpres))*un(-2:2) )


          i = lo(1)+3
          ! use 6th-order stencil
          
          un(-3:3) = q(i-3:i+3,j,qu)
          
          rhs(i,j,irho) = rhs(i,j,irho) - dxinv(1) * &
               first_deriv_6( cons(i-3:i+3,j,imx) ) 

          rhs(i,j,imx) = rhs(i,j,imx) - dxinv(1) * &
               first_deriv_6( cons(i-3:i+3,j,imx)*un(-3:3)+q(i-3:i+3,j,qpres) )

          rhs(i,j,imy) = rhs(i,j,imy) - dxinv(1) * &
               first_deriv_6( cons(i-3:i+3,j,imy)*un(-3:3) ) 

          rhs(i,j,iene) = rhs(i,j,iene) - dxinv(1) * &
               first_deriv_6( (cons(i-3:i+3,j,iene)+q(i-3:i+3,j,qpres))*un(-3:3) )
          
       enddo
       
       do n = iry1, iry1+nspecies-1
          do j=lo(2),hi(2)
             
             ! if (bclo(1) .eq. WALL???) then
             !    i = lo(1)
             !    ! use completely right-biased stencil
             ! end if
             
             i = lo(1)+1
             ! use 3rd-order slightly right-biased stencil
             rhs(i,j,n) = rhs(i,j,n) - dxinv(1) * &
                  first_deriv_r3( cons(i-1:i+2,j,n)*q(i-1:i+2,j,qu) )
             
             i = lo(1)+2
             ! use 4th-order stencil
             rhs(i,j,n) = rhs(i,j,n) - dxinv(1) * &
                  first_deriv_4( cons(i-2:i+2,j,n)*q(i-2:i+2,j,qu) )
             
             i = lo(1)+3
             ! use 6th-order stencil
             rhs(i,j,n) = rhs(i,j,n) - dxinv(1) * &
                  first_deriv_6( cons(i-3:i+3,j,n)*q(i-3:i+3,j,qu) )
          end do
       end do
    end if
    
    !
    ! ----- hi-x boundary -----
    !
    if (physbchi(1)) then 
       do j=lo(2),hi(2)
          
          i = hi(1)-3
          ! use 6th-order stencil
          
          un(-3:3) = q(i-3:i+3,j,qu)
          
          rhs(i,j,irho) = rhs(i,j,irho) - dxinv(1) * &
               first_deriv_6( cons(i-3:i+3,j,imx) ) 
          
          rhs(i,j,imx) = rhs(i,j,imx) - dxinv(1) * &
               first_deriv_6( cons(i-3:i+3,j,imx)*un(-3:3)+q(i-3:i+3,j,qpres) )
          
          rhs(i,j,imy) = rhs(i,j,imy) - dxinv(1) * &
               first_deriv_6( cons(i-3:i+3,j,imy)*un(-3:3) ) 
          
          rhs(i,j,iene) = rhs(i,j,iene) - dxinv(1) * &
               first_deriv_6( (cons(i-3:i+3,j,iene)+q(i-3:i+3,j,qpres))*un(-3:3) )
          
          i = hi(1)-2
          ! use 4th-order stencil
          
          un(-2:2) = q(i-2:i+2,j,qu)
          
          rhs(i,j,irho) = rhs(i,j,irho) - dxinv(1) * &
               first_deriv_4( cons(i-2:i+2,j,imx) ) 
          
          rhs(i,j,imx) = rhs(i,j,imx) - dxinv(1) * &
               first_deriv_4( cons(i-2:i+2,j,imx)*un(-2:2)+q(i-2:i+2,j,qpres) )
          
          rhs(i,j,imy) = rhs(i,j,imy) - dxinv(1) * &
               first_deriv_4( cons(i-2:i+2,j,imy)*un(-2:2) ) 
          
          rhs(i,j,iene) = rhs(i,j,iene) - dxinv(1) * &
               first_deriv_4( (cons(i-2:i+2,j,iene)+q(i-2:i+2,j,qpres))*un(-2:2) )
          
          i = hi(1)-1
          ! use 3rd-order slightly left-biased stencil
          
          un(-2:1) = q(i-2:i+1,j,qu)
          
          rhs(i,j,irho) = rhs(i,j,irho) - dxinv(1) * &
               first_deriv_l3( cons(i-2:i+1,j,imx) ) 
          
          rhs(i,j,imx) = rhs(i,j,imx) - dxinv(1) * &
               first_deriv_l3( cons(i-2:i+1,j,imx)*un(-2:1)+q(i-2:i+1,j,qpres) )
          
          rhs(i,j,imy) = rhs(i,j,imy) - dxinv(1) * &
               first_deriv_l3( cons(i-2:i+1,j,imy)*un(-2:1) ) 
          
          rhs(i,j,iene) = rhs(i,j,iene) - dxinv(1) * &
               first_deriv_l3( (cons(i-2:i+1,j,iene)+q(i-2:i+1,j,qpres))*un(-2:1) )
          
          
          ! if (bchi(1) .eq. WALL???) then
          !    i = hi(1)
          !    ! use completely left-biased stencil
          ! end if
          
       enddo
       
       do n = iry1, iry1+nspecies-1
          do j=lo(2),hi(2)
             
             i = hi(1)-3
             ! use 6th-order stencil
             rhs(i,j,n) = rhs(i,j,n) - dxinv(1) * &
                  first_deriv_6( cons(i-3:i+3,j,n)*q(i-3:i+3,j,qu) )
             
             i = hi(1)-2
             ! use 4th-order stencil
             rhs(i,j,n) = rhs(i,j,n) - dxinv(1) * &
                  first_deriv_4( cons(i-2:i+2,j,n)*q(i-2:i+2,j,qu) )
             
             i = hi(1)-1
             ! use 3rd-order slightly left-biased stencil
             rhs(i,j,n) = rhs(i,j,n) - dxinv(1) * &
                  first_deriv_l3( cons(i-2:i+1,j,n)*q(i-2:i+1,j,qu) )
             
             ! if (bchi(1) .eq. WALL???) then
             !    i = hi(1)
             !    ! use completely left-biased stencil
             ! end if
             
          enddo
       end do
    end if
    
    !
    ! ----- lo-y boundary -----
    !
    if (physbclo(2)) then 

       ! if (bclo(2) .eq. WALL???) then
       !    j = lo(2)
       !    ! use completely right-biased stencil
       !    enddo
       ! end if
       
       j = lo(2)+1
       ! use 3rd-order slightly right-biased stencil
       
       do i=lo(1),hi(1)
          
          un(-1:2) = q(i,j-1:j+2,qv)
          
          rhs(i,j,irho)=rhs(i,j,irho) - dxinv(2) * &
               first_deriv_r3( cons(i,j-1:j+2,imy) )
          
          rhs(i,j,imx)=rhs(i,j,imx) - dxinv(2) * &
               first_deriv_r3( cons(i,j-1:j+2,imx)*un(-1:2) )
          
          rhs(i,j,imy)=rhs(i,j,imy) - dxinv(2) * &
               first_deriv_r3( cons(i,j-1:j+2,imy)*un(-1:2)+q(i,j-1:j+2,qpres) )
          
          rhs(i,j,iene)=rhs(i,j,iene) - dxinv(2) * &
               first_deriv_r3( (cons(i,j-1:j+2,iene)+q(i,j-1:j+2,qpres))*un(-1:2) )
          
       enddo
       
       j = lo(2)+2
       ! use 4th-order stencil
       
       do i=lo(1),hi(1)
          
          un(-2:2) = q(i,j-2:j+2,qv)
          
          rhs(i,j,irho)=rhs(i,j,irho) - dxinv(2) * &
               first_deriv_4( cons(i,j-2:j+2,imy) )
          
          rhs(i,j,imx)=rhs(i,j,imx) - dxinv(2) * &
               first_deriv_4( cons(i,j-2:j+2,imx)*un(-2:2) )
          
          rhs(i,j,imy)=rhs(i,j,imy) - dxinv(2) * &
               first_deriv_4( cons(i,j-2:j+2,imy)*un(-2:2)+q(i,j-2:j+2,qpres) )
          
          rhs(i,j,iene)=rhs(i,j,iene) - dxinv(2) * &
               first_deriv_4( (cons(i,j-2:j+2,iene)+q(i,j-2:j+2,qpres))*un(-2:2) )
          
       enddo
       
       j = lo(2)+3
       ! use 6th-order stencil
       
       do i=lo(1),hi(1)
          
          un(-3:3) = q(i,j-3:j+3,qv)
          
          rhs(i,j,irho)=rhs(i,j,irho) - dxinv(2) * &
               first_deriv_6( cons(i,j-3:j+3,imy) )
          
          rhs(i,j,imx)=rhs(i,j,imx) - dxinv(2) * &
               first_deriv_6( cons(i,j-3:j+3,imx)*un(-3:3) )
          
          rhs(i,j,imy)=rhs(i,j,imy) - dxinv(2) * &
               first_deriv_6( cons(i,j-3:j+3,imy)*un(-3:3)+q(i,j-3:j+3,qpres) )
          
          rhs(i,j,iene)=rhs(i,j,iene) - dxinv(2) * &
               first_deriv_6( (cons(i,j-3:j+3,iene)+q(i,j-3:j+3,qpres))*un(-3:3) )
          
       enddo


       do n = iry1, iry1+nspecies-1
          
          ! if (bclo(2) .eq. WALL???) then
          !    j = lo(2)
          !    ! use completely right-biased stencil
          !    enddo
          ! end if
          
          j = lo(2)+1
          ! use 3rd-order slightly right-biased stencil
          do i=lo(1),hi(1)
             rhs(i,j,n) = rhs(i,j,n) - dxinv(2) * &
                  first_deriv_r3( cons(i,j-1:j+2,n)*q(i,j-1:j+2,qv) )
          end do
          
          j = lo(2)+2
          ! use 4th-order stencil
          do i=lo(1),hi(1)
             rhs(i,j,n) = rhs(i,j,n) - dxinv(2) * &
                  first_deriv_4( cons(i,j-2:j+2,n)*q(i,j-2:j+2,qv) )
          end do
          
          j = lo(2)+3
          ! use 6th-order stencil
          do i=lo(1),hi(1)
             rhs(i,j,n) = rhs(i,j,n) - dxinv(2) * &
                  first_deriv_6( cons(i,j-3:j+3,n)*q(i,j-3:j+3,qv) )
          end do
          
       end do
    end if

    !
    ! ----- hi-y boundary -----
    !
    if (physbchi(2)) then 
       
       j = hi(2)-3
       ! use 6th-order stencil
       
       do i=lo(1),hi(1)
          
          un(-3:3) = q(i,j-3:j+3,qv)
          
          rhs(i,j,irho)=rhs(i,j,irho) - dxinv(2) * &
               first_deriv_6( cons(i,j-3:j+3,imy) )
          
          rhs(i,j,imx)=rhs(i,j,imx) - dxinv(2) * &
               first_deriv_6( cons(i,j-3:j+3,imx)*un(-3:3) )
          
          rhs(i,j,imy)=rhs(i,j,imy) - dxinv(2) * &
               first_deriv_6( cons(i,j-3:j+3,imy)*un(-3:3)+q(i,j-3:j+3,qpres) )
          
          rhs(i,j,iene)=rhs(i,j,iene) - dxinv(2) * &
               first_deriv_6( (cons(i,j-3:j+3,iene)+q(i,j-3:j+3,qpres))*un(-3:3) )
          
       enddo
       
       j = hi(2)-2
       ! use 4th-order stencil
       
       do i=lo(1),hi(1)
          
          un(-2:2) = q(i,j-2:j+2,qv)
          
          rhs(i,j,irho)=rhs(i,j,irho) - dxinv(2) * &
               first_deriv_4( cons(i,j-2:j+2,imy) )
          
          rhs(i,j,imx)=rhs(i,j,imx) - dxinv(2) * &
               first_deriv_4( cons(i,j-2:j+2,imx)*un(-2:2) )
          
          rhs(i,j,imy)=rhs(i,j,imy) - dxinv(2) * &
               first_deriv_4( cons(i,j-2:j+2,imy)*un(-2:2)+q(i,j-2:j+2,qpres) )
          
          rhs(i,j,iene)=rhs(i,j,iene) - dxinv(2) * &
               first_deriv_4( (cons(i,j-2:j+2,iene)+q(i,j-2:j+2,qpres))*un(-2:2) )
          
       enddo
       
       j = hi(2)-1
       ! use 3rd-order slightly left-biased stencil
       
       do i=lo(1),hi(1)
          
          un(-2:1) = q(i,j-2:j+1,qv)
          
          rhs(i,j,irho)=rhs(i,j,irho) - dxinv(2) * &
               first_deriv_l3( cons(i,j-2:j+1,imy) )
          
          rhs(i,j,imx)=rhs(i,j,imx) - dxinv(2) * &
               first_deriv_l3( cons(i,j-2:j+1,imx)*un(-2:1) )
          
          rhs(i,j,imy)=rhs(i,j,imy) - dxinv(2) * &
               first_deriv_l3( cons(i,j-2:j+1,imy)*un(-2:1)+q(i,j-2:j+1,qpres) )
          
          rhs(i,j,iene)=rhs(i,j,iene) - dxinv(2) * &
               first_deriv_l3( (cons(i,j-2:j+1,iene)+q(i,j-2:j+1,qpres))*un(-2:1) )
          
       enddo
       
       ! if (bchi(2) .eq. WALL???) then
       !    j = hi(2)
       !    ! use completely left-biased stencil
       ! end if
       
       do n = iry1, iry1+nspecies-1
          
          j = hi(2)-3
          ! use 6th-order stencil
          do i=lo(1),hi(1)
             rhs(i,j,n) = rhs(i,j,n) - dxinv(2) * &
                  first_deriv_6( cons(i,j-3:j+3,n)*q(i,j-3:j+3,qv) )
          end do
          
          j = hi(2)-2
          ! use 4th-order stencil
          do i=lo(1),hi(1)
             rhs(i,j,n) = rhs(i,j,n) - dxinv(2) * &
                  first_deriv_4( cons(i,j-2:j+2,n)*q(i,j-2:j+2,qv) )
          end do
          
          j = hi(2)-1
          ! use 3rd-order slightly left-biased stencil
          do i=lo(1),hi(1)
             rhs(i,j,n) = rhs(i,j,n) - dxinv(2) * &
                  first_deriv_l3( cons(i,j-2:j+1,n)*q(i,j-2:j+1,qv) )
          end do
          
          ! if (bchi(2) .eq. WALL???) then
          !    j = hi(2)
          !    ! use completely left-biased stencil
          ! end if
          
       end do
    end if
    
    deallocate(tmpx,tmpy)
    
    rhs_g(lo(1):hi(1),lo(2):hi(2),:) = &
         rhs_g(lo(1):hi(1),lo(2):hi(2),:) &
         + rhs(lo(1):hi(1),lo(2):hi(2),:)
    deallocate(rhs)

  end subroutine hypterm_2d


  subroutine narrow_diffterm_2d (lo,hi,dx,q,qlo,qhi,rhs_g,glo,ghi,rhs,rlo,rhi, &
       mu,xi,lam,dxy,dlo_g,dhi_g,bclo,bchi)

    integer,         intent(in):: dlo_g(2),dhi_g(2),bclo(2),bchi(2)
    integer,         intent(in):: lo(2),hi(2),qlo(2),qhi(2),rlo(2),rhi(2),glo(2),ghi(2)
    double precision,intent(in):: dx(2)
    double precision,intent(in)   ::  q  (qlo(1):qhi(1),qlo(2):qhi(2),nprim)
    double precision,intent(in)   ::  mu (qlo(1):qhi(1),qlo(2):qhi(2))
    double precision,intent(in)   ::  xi (qlo(1):qhi(1),qlo(2):qhi(2))
    double precision,intent(in)   ::  lam(qlo(1):qhi(1),qlo(2):qhi(2))
    double precision,intent(in)   ::  dxy(qlo(1):qhi(1),qlo(2):qhi(2),nspecies)
    double precision,intent(inout)::rhs_g(glo(1):ghi(1),glo(2):ghi(2),ncons)
    double precision,intent(inout)::rhs  (rlo(1):rhi(1),rlo(2):rhi(2),ncons)

    integer :: i
    integer :: slo(2), shi(2), dlo(2), dhi(2)
    double precision :: dxinv(2), dx2inv(2)

    ! used to turn off some terms
    double precision :: finlo(2), finhi(2)
    double precision :: foulo(2), fouhi(2)

    logical :: physbclo(2), physbchi(2)

    ! Only the region bounded by [dlo_g,dhi_g] contains good data.
    ! [slo,shi] will be safe for 8th-order stencil
    do i=1,2
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

    do i=1,2
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

    do i = 1,2
       dxinv(i) = 1.0d0 / dx(i)
       dx2inv(i) = dxinv(i)**2
    end do

    rhs(lo(1):hi(1),lo(2):hi(2),:) = 0.d0

    call diffterm_1(q,qlo,qhi,rhs,rlo,rhi,mu,xi, &
         lo,hi,slo,shi,dlo,dhi,finlo,finhi,foulo,fouhi,physbclo,physbchi,dxinv)

    call diffterm_2(q,qlo,qhi,rhs,rlo,rhi, mu,xi,lam,dxy, &
         lo,hi,slo,shi,dlo,dhi,finlo,finhi,foulo,fouhi,physbclo,physbchi,dxinv,dx2inv)

    rhs_g     (lo(1):hi(1),lo(2):hi(2),:) = &
         rhs_g(lo(1):hi(1),lo(2):hi(2),:) &
         + rhs(lo(1):hi(1),lo(2):hi(2),:)

  end subroutine narrow_diffterm_2d

  
  subroutine diffterm_1(q,qlo,qhi,rhs,rlo,rhi,mu,xi, &
       lo,hi,slo,shi,dlo,dhi,finlo,finhi,foulo,fouhi,physbclo,physbchi,dxinv)
    integer,         intent(in):: lo(2),hi(2),slo(2),shi(2),dlo(2),dhi(2)
    integer,         intent(in):: qlo(2),qhi(2),rlo(2),rhi(2)
    logical,         intent(in):: physbclo(2),physbchi(2)
    double precision,intent(in):: finlo(2),finhi(2),foulo(2),fouhi(2)
    double precision,intent(in):: dxinv(2)
    double precision,intent(in)   ::  q(qlo(1):qhi(1),qlo(2):qhi(2),nprim)
    double precision,intent(in)   :: mu(qlo(1):qhi(1),qlo(2):qhi(2))
    double precision,intent(in)   :: xi(qlo(1):qhi(1),qlo(2):qhi(2))
    double precision,intent(inout)::rhs(rlo(1):rhi(1),rlo(2):rhi(2),ncons)

    double precision, allocatable, dimension(:,:) :: ux,uy,vx,vy
    double precision, allocatable :: tmpx(:), tmpy(:,:)
    double precision, allocatable, dimension(:,:) :: vsm
    double precision, dimension(lo(1):hi(1)) :: tauxx,tauyy,divu
    integer          :: i,j

    allocate(ux( lo(1): hi(1),dlo(2):dhi(2)))
    allocate(vx( lo(1): hi(1),dlo(2):dhi(2)))

    allocate(uy(dlo(1):dhi(1), lo(2): hi(2)))
    allocate(vy(dlo(1):dhi(1), lo(2): hi(2)))

    allocate(vsm(dlo(1):dhi(1),dlo(2):dhi(2)))

    allocate(tmpx(dlo(1):dhi(1)))
    allocate(tmpy( lo(1): hi(1),dlo(2):dhi(2)))

    do j=dlo(2),dhi(2)
       do i=dlo(1),dhi(1)
          vsm(i,j) = xi(i,j) -  TwoThirds*mu(i,j)
       enddo
    enddo

    do j=dlo(2),dhi(2)
       do i=slo(1),shi(1)
          ux(i,j) = dxinv(1)*first_deriv_8(q(i-4:i+4,j,qu))
          vx(i,j) = dxinv(1)*first_deriv_8(q(i-4:i+4,j,qv))
       enddo
    enddo

    do j=slo(2),shi(2)   
       do i=dlo(1),dhi(1)
          uy(i,j) = dxinv(2)*first_deriv_8(q(i,j-4:j+4,qu))
          vy(i,j) = dxinv(2)*first_deriv_8(q(i,j-4:j+4,qv))
       enddo
    enddo

    !
    ! lo-x boundary
    !
    if (physbclo(1)) then
       do j=dlo(2),dhi(2)
          i = lo(1)
          ! use completely right-biased stencil
          ux(i,j) = dxinv(1)*first_deriv_rb(q(i:i+3,j,qu))
          vx(i,j) = dxinv(1)*first_deriv_rb(q(i:i+3,j,qv))

          i = lo(1)+1
          ! use 3rd-order slightly right-biased stencil
          ux(i,j) = dxinv(1)*first_deriv_r3(q(i-1:i+2,j,qu))
          vx(i,j) = dxinv(1)*first_deriv_r3(q(i-1:i+2,j,qv))

          i = lo(1)+2
          ! use 4th-order stencil
          ux(i,j) = dxinv(1)*first_deriv_4(q(i-2:i+2,j,qu))
          vx(i,j) = dxinv(1)*first_deriv_4(q(i-2:i+2,j,qv))

          i = lo(1)+3
          ! use 6th-order stencil
          ux(i,j) = dxinv(1)*first_deriv_6(q(i-3:i+3,j,qu))
          vx(i,j) = dxinv(1)*first_deriv_6(q(i-3:i+3,j,qv))
       end do
    end if

    !
    ! hi-x boundary
    !
    if (physbchi(1)) then
       do j=dlo(2),dhi(2)
          i = hi(1)-3
          ! use 6th-order stencil
          ux(i,j) = dxinv(1)*first_deriv_6(q(i-3:i+3,j,qu))
          vx(i,j) = dxinv(1)*first_deriv_6(q(i-3:i+3,j,qv))

          i = hi(1)-2
          ! use 4th-order stencil
          ux(i,j) = dxinv(1)*first_deriv_4(q(i-2:i+2,j,qu))
          vx(i,j) = dxinv(1)*first_deriv_4(q(i-2:i+2,j,qv))

          i = hi(1)-1
          ! use 3rd-order slightly left-biased stencil
          ux(i,j) = dxinv(1)*first_deriv_l3(q(i-2:i+1,j,qu))
          vx(i,j) = dxinv(1)*first_deriv_l3(q(i-2:i+1,j,qv))

          i = hi(1)
          ! use completely left-biased stencil
          ux(i,j) = dxinv(1)*first_deriv_lb(q(i-3:i,j,qu))
          vx(i,j) = dxinv(1)*first_deriv_lb(q(i-3:i,j,qv))
       end do
    end if

    !
    ! lo-y boundary
    !
    if (physbclo(2)) then
       j = lo(2)
       ! use completely right-biased stencil
       do i=dlo(1),dhi(1)
          uy(i,j) = dxinv(2)*first_deriv_rb(q(i,j:j+3,qu))
          vy(i,j) = dxinv(2)*first_deriv_rb(q(i,j:j+3,qv))
       enddo
       
       j = lo(2)+1
       ! use 3rd-order slightly right-biased stencil
       do i=dlo(1),dhi(1)
          uy(i,j) = dxinv(2)*first_deriv_r3(q(i,j-1:j+2,qu))
          vy(i,j) = dxinv(2)*first_deriv_r3(q(i,j-1:j+2,qv))
       enddo

       j = lo(2)+2
       ! use 4th-order stencil
       do i=dlo(1),dhi(1)
          uy(i,j) = dxinv(2)*first_deriv_4(q(i,j-2:j+2,qu))
          vy(i,j) = dxinv(2)*first_deriv_4(q(i,j-2:j+2,qv))
       enddo
       
       j = lo(2)+3
       ! use 6th-order stencil
       do i=dlo(1),dhi(1)
          uy(i,j) = dxinv(2)*first_deriv_6(q(i,j-3:j+3,qu))
          vy(i,j) = dxinv(2)*first_deriv_6(q(i,j-3:j+3,qv))
       enddo
    end if

    !
    ! hi-y boundary
    !
    if (physbchi(2)) then
       j = hi(2)-3
       ! use 6th-order stencil
       do i=dlo(1),dhi(1)
          uy(i,j) = dxinv(2)*first_deriv_6(q(i,j-3:j+3,qu))
          vy(i,j) = dxinv(2)*first_deriv_6(q(i,j-3:j+3,qv))
       enddo
       
       j = hi(2)-2
       ! use 4th-order stencil
       do i=dlo(1),dhi(1)
          uy(i,j) = dxinv(2)*first_deriv_4(q(i,j-2:j+2,qu))
          vy(i,j) = dxinv(2)*first_deriv_4(q(i,j-2:j+2,qv))
       enddo
       
       j = hi(2)-1
       ! use 3rd-order slightly left-biased stencil
       do i=dlo(1),dhi(1)
          uy(i,j) = dxinv(2)*first_deriv_l3(q(i,j-2:j+1,qu))
          vy(i,j) = dxinv(2)*first_deriv_l3(q(i,j-2:j+1,qv))
       enddo
       
       j = hi(2)
       ! use completely left-biased stencil
       do i=dlo(1),dhi(1)
          uy(i,j) = dxinv(2)*first_deriv_lb(q(i,j-3:j,qu))
          vy(i,j) = dxinv(2)*first_deriv_lb(q(i,j-3:j,qv))
       enddo
    end if

    !----- mx -----

    !----- mx : d()/dx -----
    do j=lo(2),hi(2)

       do i=dlo(1),dhi(1)
          tmpx(i) = vsm(i,j)*vy(i,j)
       end do

       !
       ! lo-x boundary
       !
       if (physbclo(1)) then
          i = lo(1)
          ! use completely right-biased stencil
          rhs(i,j,imx) = rhs(i,j,imx) + finlo(1)*dxinv(1)*first_deriv_rb(tmpx(i:i+3))
          
          i = lo(1)+1
          ! use 3rd-order slightly right-biased stencil
          rhs(i,j,imx) = rhs(i,j,imx) + dxinv(1)*first_deriv_r3(tmpx(i-1:i+2))
          
          i = lo(1)+2
          ! use 4th-order stencil
          rhs(i,j,imx) = rhs(i,j,imx) + dxinv(1)*first_deriv_4(tmpx(i-2:i+2))
          
          i = lo(1)+3
          ! use 6th-order stencil
          rhs(i,j,imx) = rhs(i,j,imx) + dxinv(1)*first_deriv_6(tmpx(i-3:i+3))
       end if
       
       do i=slo(1),shi(1)
          rhs(i,j,imx) = rhs(i,j,imx) + dxinv(1)*first_deriv_8(tmpx(i-4:i+4)) 
       end do

       !
       ! hi-x boundary
       !
       if (physbchi(1)) then
          i = hi(1)-3
          ! use 6th-order stencil
          rhs(i,j,imx) = rhs(i,j,imx) + dxinv(1)*first_deriv_6(tmpx(i-3:i+3))
          
          i = hi(1)-2
          ! use 4th-order stencil
          rhs(i,j,imx) = rhs(i,j,imx) + dxinv(1)*first_deriv_4(tmpx(i-2:i+2))
          
          i = hi(1)-1
          ! use 3rd-order slightly left-biased stencil
          rhs(i,j,imx) = rhs(i,j,imx) + dxinv(1)*first_deriv_l3(tmpx(i-2:i+1))
          
          i = hi(1)
          ! use completely left-biased stencil
          rhs(i,j,imx) = rhs(i,j,imx) + finhi(1)*dxinv(1)*first_deriv_lb(tmpx(i-3:i))
       end if
    end do
    
    !----- mx : d()/dy -----
    do j=dlo(2),dhi(2)
       do i=lo(1),hi(1)
          tmpy(i,j) = mu(i,j)*vx(i,j)
       end do
    end do

    !
    ! lo-y boundary
    !
    if (physbclo(2)) then
       j = lo(2)
       ! use completely right-biased stencil
       do i=lo(1),hi(1)
          rhs(i,j,imx) = rhs(i,j,imx) + foulo(2)*dxinv(2)*first_deriv_rb(tmpy(i,j:j+3))
       end do
       
       j = lo(2)+1
       ! use 3rd-order slightly right-biased stencil
       do i=lo(1),hi(1)
          rhs(i,j,imx) = rhs(i,j,imx) + dxinv(2)*first_deriv_r3(tmpy(i,j-1:j+2))
       end do
       
       j = lo(2)+2
       ! use 4th-order stencil
       do i=lo(1),hi(1)
          rhs(i,j,imx) = rhs(i,j,imx) + dxinv(2)*first_deriv_4(tmpy(i,j-2:j+2))
       end do
       
       j = lo(2)+3
       ! use 6th-order stencil
       do i=lo(1),hi(1)
          rhs(i,j,imx) = rhs(i,j,imx) + dxinv(2)*first_deriv_6(tmpy(i,j-3:j+3))
       end do
    end if

    !
    ! interior
    !
    do j=slo(2),shi(2)
       do i=lo(1),hi(1)
          rhs(i,j,imx) = rhs(i,j,imx) + dxinv(2)*first_deriv_8(tmpy(i,j-4:j+4)) 
       end do
    end do
    
    !
    ! hi-y boundary
    !
    if (physbchi(2)) then
       j = hi(2)-3
       ! use 6th-order stencil
       do i=lo(1),hi(1)
          rhs(i,j,imx) = rhs(i,j,imx) + dxinv(2)*first_deriv_6(tmpy(i,j-3:j+3))
       end do
       
       j = hi(2)-2
       ! use 4th-order stencil
       do i=lo(1),hi(1)
          rhs(i,j,imx) = rhs(i,j,imx) + dxinv(2)*first_deriv_4(tmpy(i,j-2:j+2))
       end do
       
       j = hi(2)-1
       ! use 3rd-order slightly left-biased stencil
       do i=lo(1),hi(1)
          rhs(i,j,imx) = rhs(i,j,imx) + dxinv(2)*first_deriv_l3(tmpy(i,j-2:j+1))
       end do
       
       j = hi(2)
       ! use completely left-biased stencil
       do i=lo(1),hi(1)
          rhs(i,j,imx) = rhs(i,j,imx) + fouhi(2)*dxinv(2)*first_deriv_lb(tmpy(i,j-3:j))
       end do
    end if

    !----- my -----

    !----- my : d()/dx -----
    do j=lo(2),hi(2)
       
       do i=dlo(1),dhi(1)
          tmpx(i) = mu(i,j)*uy(i,j)
       end do
       
       !
       ! lo-x boundary
       !
       if (physbclo(1)) then
          i = lo(1)
          rhs(i,j,imy) = rhs(i,j,imy) + foulo(1)*dxinv(1)*first_deriv_rb(tmpx(i:i+3))
          
          i = lo(1)+1
          ! use 3rd-order slightly right-biased stencil
          rhs(i,j,imy) = rhs(i,j,imy) + dxinv(1)*first_deriv_r3(tmpx(i-1:i+2))
          
          i = lo(1)+2
          ! use 4th-order stencil
          rhs(i,j,imy) = rhs(i,j,imy) + dxinv(1)*first_deriv_4(tmpx(i-2:i+2))
          
          i = lo(1)+3
          ! use 6th-order stencil
          rhs(i,j,imy) = rhs(i,j,imy) + dxinv(1)*first_deriv_6(tmpx(i-3:i+3))
       end if
       
       !
       ! interior
       !
       do i=slo(1),shi(1)
          rhs(i,j,imy) = rhs(i,j,imy) + dxinv(1)*first_deriv_8(tmpx(i-4:i+4)) 
       end do
       
       !
       ! hi-x boundary
       !
       if (physbchi(1)) then
          i = hi(1)-3
          ! use 6th-order stencil
          rhs(i,j,imy) = rhs(i,j,imy) + dxinv(1)*first_deriv_6(tmpx(i-3:i+3))
          
          i = hi(1)-2
          ! use 4th-order stencil
          rhs(i,j,imy) = rhs(i,j,imy) + dxinv(1)*first_deriv_4(tmpx(i-2:i+2))
          
          i = hi(1)-1
          ! use 3rd-order slightly left-biased stencil
          rhs(i,j,imy) = rhs(i,j,imy) + dxinv(1)*first_deriv_l3(tmpx(i-2:i+1))
          
          i = hi(1)
          ! use completely left-biased stencil
          rhs(i,j,imy) = rhs(i,j,imy) + fouhi(1)*dxinv(1)*first_deriv_lb(tmpx(i-3:i))
       end if
    end do

    !----- my : d()/dy -----

    do j=dlo(2),dhi(2)
       do i=lo(1),hi(1)
          tmpy(i,j) = vsm(i,j)*ux(i,j)
       end do
    end do
    
    !
    ! lo-y boundary
    !
    if (physbclo(2)) then
       j = lo(2)
       ! use completely right-biased stencil
       do i=lo(1),hi(1)
          rhs(i,j,imy) = rhs(i,j,imy) + finlo(2)*dxinv(2)*first_deriv_rb(tmpy(i,j:j+3))
       end do
       
       j = lo(2)+1
       ! use 3rd-order slightly right-biased stencil
       do i=lo(1),hi(1)
          rhs(i,j,imy) = rhs(i,j,imy) + dxinv(2)*first_deriv_r3(tmpy(i,j-1:j+2))
       end do
       
       j = lo(2)+2
       ! use 4th-order stencil
       do i=lo(1),hi(1)
          rhs(i,j,imy) = rhs(i,j,imy) + dxinv(2)*first_deriv_4(tmpy(i,j-2:j+2))
       end do
       
       j = lo(2)+3
       ! use 6th-order stencil
       do i=lo(1),hi(1)
          rhs(i,j,imy) = rhs(i,j,imy) + dxinv(2)*first_deriv_6(tmpy(i,j-3:j+3))
       end do
    end if
    
    !
    ! interior
    !
    do j=slo(2),shi(2)
       do i=lo(1),hi(1)
          rhs(i,j,imy) = rhs(i,j,imy) + dxinv(2)*first_deriv_8(tmpy(i,j-4:j+4)) 
       end do
    end do
    
    !
    ! hi-y boundary
    !
    if (physbchi(2)) then
       j = hi(2)-3
       ! use 6th-order stencil
       do i=lo(1),hi(1)
          rhs(i,j,imy) = rhs(i,j,imy) + dxinv(2)*first_deriv_6(tmpy(i,j-3:j+3))
       end do
       
       j = hi(2)-2
       ! use 4th-order stencil
       do i=lo(1),hi(1)
          rhs(i,j,imy) = rhs(i,j,imy) + dxinv(2)*first_deriv_4(tmpy(i,j-2:j+2))
       end do
       
       j = hi(2)-1
       ! use 3rd-order slightly left-biased stencil
       do i=lo(1),hi(1)
          rhs(i,j,imy) = rhs(i,j,imy) + dxinv(2)*first_deriv_l3(tmpy(i,j-2:j+1))
       end do
       
       j = hi(2)
       ! use completely left-biased stencil
       do i=lo(1),hi(1)
          rhs(i,j,imy) = rhs(i,j,imy) + finhi(2)*dxinv(2)*first_deriv_lb(tmpy(i,j-3:j))
       end do
    end if

    !----- energy -----

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          
          divu(i) = (ux(i,j)+vy(i,j))*vsm(i,j)
          tauxx(i) = 2.d0*mu(i,j)*ux(i,j) + divu(i)
          tauyy(i) = 2.d0*mu(i,j)*vy(i,j) + divu(i)
          
          ! change in internal energy
          rhs(i,j,iene) = rhs(i,j,iene) + &
               tauxx(i)*ux(i,j) + tauyy(i)*vy(i,j) &
               + mu(i,j)*(uy(i,j)+vx(i,j))**2 

       end do
    end do

    deallocate(tmpx,tmpy)
    deallocate(vsm)
    deallocate(ux,uy,vx,vy)

  end subroutine diffterm_1

  
  subroutine diffterm_2(q,qlo,qhi,rhs,rlo,rhi,mu,xi,lam,dxy, &
       lo,hi,slo,shi,dlo,dhi,finlo,finhi,foulo,fouhi,physbclo,physbchi,dxinv,dx2inv)
    use probin_module, only : reset_inactive_species
    integer,         intent(in):: lo(2),hi(2),slo(2),shi(2),dlo(2),dhi(2)
    integer,         intent(in):: qlo(2),qhi(2),rlo(2),rhi(2)
    logical,         intent(in):: physbclo(2),physbchi(2)
    double precision,intent(in):: finlo(2),finhi(2),foulo(2),fouhi(2)
    double precision,intent(in):: dxinv(2),dx2inv(2)
    double precision,intent(in)   :: q (qlo(1):qhi(1),qlo(2):qhi(2),nprim)
    double precision,intent(in)   :: mu(qlo(1):qhi(1),qlo(2):qhi(2))
    double precision,intent(in)   :: xi(qlo(1):qhi(1),qlo(2):qhi(2))
    double precision,intent(in)   ::lam(qlo(1):qhi(1),qlo(2):qhi(2))
    double precision,intent(in)   ::dxy(qlo(1):qhi(1),qlo(2):qhi(2),nspecies)
    double precision,intent(inout)::rhs(rlo(1):rhi(1),rlo(2):rhi(2),ncons)

    double precision, allocatable, dimension(:,:) :: vsp, dpe
    double precision, allocatable, dimension(:,:,:) :: Hg, dpy, dxe
    ! dxy: diffusion coefficient of X in equation for Y
    ! dpy: diffusion coefficient of p in equation for Y
    ! dxe: diffusion coefficient of X in equation for energy
    ! dpe: diffusion coefficient of p in equation for energy

    integer          :: i,j,n, qxn, qyn, qhn, iryn

    double precision :: mmtmp8(8,lo(1):hi(1)+1)
    double precision, allocatable, dimension(:,:,:) :: M8p
    double precision, allocatable, dimension(:,:) :: sumdrY, sumrYv, gradp
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

    allocate(vsp(dlo(1):dhi(1),dlo(2):dhi(2)))

    allocate(dpy(dlo(1):dhi(1),dlo(2):dhi(2),nspecies))
    allocate(dxe(dlo(1):dhi(1),dlo(2):dhi(2),nspecies))
    allocate(dpe(dlo(1):dhi(1),dlo(2):dhi(2)))

    allocate(Hg(lo(1):hi(1)+1,lo(2):hi(2)+1,2:ncons))

    allocate(M8p(8,lo(1):hi(1)+1,lo(2):hi(2)+1))

    allocate(sumdrY(lo(1):hi(1),lo(2):hi(2)))
    allocate(sumrYv(lo(1):hi(1),lo(2):hi(2)))
    allocate(gradp (lo(1):hi(1),lo(2):hi(2)))

    do j=dlo(2),dhi(2)
       do i=dlo(1),dhi(1)
          vsp(i,j) = xi(i,j) + FourThirds*mu(i,j)
       enddo
    enddo

    dpe = 0.d0

    do n=1,nspecies
       qxn = qx1+n-1
       qyn = qy1+n-1
       qhn = qh1+n-1
       do j=dlo(2),dhi(2)
          do i=dlo(1),dhi(1)
             dpy(i,j,n) = dxy(i,j,n)/q(i,j,qpres)*(q(i,j,qxn)-q(i,j,qyn))
             dxe(i,j,n) = dxy(i,j,n)*q(i,j,qhn)
             dpe(i,j) = dpe(i,j) + dpy(i,j,n)*q(i,j,qhn)
          end do
       end do
    end do

    ! ------- BEGIN x-direction -------

    do j=lo(2),hi(2)
       do i=slo(1),shi(1)+1
          mmtmp8(1:8,i) = matmul(vsp(i-4:i+3,j), M8)
          Hg(i,j,imx) = dot_product(mmtmp8(1:8,i), q(i-4:i+3,j,qu))
       end do
    end do

    do j=lo(2),hi(2)
       do i=slo(1),shi(1)+1
          mmtmp8(1:8,i) = matmul(mu(i-4:i+3,j), M8)
          Hg(i,j,imy) = dot_product(mmtmp8(1:8,i), q(i-4:i+3,j,qv))
       end do
    end do

    do j=lo(2),hi(2)
       do i=slo(1),shi(1)+1
          mmtmp8(1:8,i) = matmul(lam(i-4:i+3,j), M8)
          Hg(i,j,iene) = dot_product(mmtmp8(1:8,i), q(i-4:i+3,j,qtemp))
       end do
    end do

    do j=lo(2),hi(2)
       do i=slo(1),shi(1)+1
          mmtmp8(1:8,i) = matmul(M8T,  q(i-4:i+3,j,qpres))
          Hg(i,j,iene) = Hg(i,j,iene) + dot_product(dpe(i-4:i+3,j), mmtmp8(1:8,i))
       end do
       do i=slo(1),shi(1)+1
          M8p(:,i,j) = mmtmp8(1:8,i)
       end do
    end do

    do n = 1, nspecies

       if (n .eq. iias) cycle  ! inactive speices

       qxn = qx1+n-1
       iryn = iry1+n-1

       do j=lo(2),hi(2)   
          do i=slo(1),shi(1)+1
             Hg(i,j,iryn) = dot_product(dpy(i-4:i+3,j,n), M8p(:,i,j))
          end do
       end do

       do j=lo(2),hi(2)   
          do i=slo(1),shi(1)+1
             mmtmp8(1:8,i) = matmul(M8T, q(i-4:i+3,j,qxn))
             Hg(i,j,iene) = Hg(i,j,iene) + dot_product(dxe(i-4:i+3,j,n), mmtmp8(1:8,i))
             Hg(i,j,iryn) = Hg(i,j,iryn) &
                  + dot_product(dxy(i-4:i+3,j,n), mmtmp8(1:8,i))
          end do
       end do

    end do

    ! add x-direction rhs

    do n=2,iene
       do j=lo(2),hi(2)
          do i=slo(1),shi(1)
             rhs(i,j,n) = rhs(i,j,n) + (Hg(i+1,j,n) - Hg(i,j,n)) * dx2inv(1)
          end do
       end do
    end do
       
    sumdrY = 0.d0
    do n=iry1,ncons

       if (n.eq.iry_ias) cycle

       do j=lo(2),hi(2)
          do i=slo(1),shi(1)
             sumdry(i,j) = sumdry(i,j) + (Hg(i+1,j,n) - Hg(i,j,n)) * dx2inv(1)
             rhs(i,j,n)  =  rhs(i,j,n) + (Hg(i+1,j,n) - Hg(i,j,n)) * dx2inv(1)
          end do
       end do

    end do

    do j=lo(2),hi(2)
       do i=slo(1),shi(1)
          gradp(i,j) = dxinv(1) * first_deriv_8(q(i-4:i+4,j,qpres))
       end do
    end do
       
    sumryv = 0.d0
    do n = 1, nspecies

       if (n.eq.iias) cycle

       qxn = qx1+n-1
       do j=lo(2),hi(2)
          do i=slo(1),shi(1)
             sumryv(i,j) = sumryv(i,j) + dpy(i,j,n)*gradp(i,j)  &
                  + dxy(i,j,n)*dxinv(1)*first_deriv_8(q(i-4:i+4,j,qxn))
          end do
       end do

    end do

    if (add_v_correction) then

       do n=1,nspecies
          qyn = qy1+n-1
          qhn = qh1+n-1
          iryn = iry1+n-1
          
          do j=lo(2),hi(2)
             do i=slo(1),shi(1)
                ry_c = q(i,j,qyn)*sumdry(i,j) + sumryv(i,j)*dxinv(1) * &
                     first_deriv_8(q(i-4:i+4,j,qyn))
                ene_c = ry_c*q(i,j,qhn) + q(i,j,qyn)*sumryv(i,j)*dxinv(1)* &
                     first_deriv_8(q(i-4:i+4,j,qhn))
                rhs(i,j,iene) = rhs(i,j,iene) - ene_c
                rhs(i,j,iryn) = rhs(i,j,iryn) - ry_c
             end do
          end do
       end do

    else
    
       n = iias
       qhn = qh1+n-1
       iryn = iry1+n-1

       do j=lo(2),hi(2)
          do i=slo(1),shi(1)
             ene_c = sumdry(i,j)*q(i,j,qhn) + sumryv(i,j)*dxinv(1)* &
                  first_deriv_8(q(i-4:i+4,j,qhn))
             rhs(i,j,iene) = rhs(i,j,iene) - ene_c
             rhs(i,j,iryn) = rhs(i,j,iryn) - sumdry(i,j)
          end do
       end do
       
    end if

    ! ------- END x-direction -------

    ! ------- BEGIN y-direction -------

    do j=slo(2),shi(2)+1
       do i=lo(1),hi(1)             
          mmtmp8(1:8,i) = matmul(mu(i,j-4:j+3), M8)
          Hg(i,j,imx) = dot_product(mmtmp8(1:8,i), q(i,j-4:j+3,qu))
       end do
    end do
    
    do j=slo(2),shi(2)+1
       do i=lo(1),hi(1)             
          mmtmp8(1:8,i) = matmul(vsp(i,j-4:j+3), M8)
          Hg(i,j,imy) = dot_product(mmtmp8(1:8,i), q(i,j-4:j+3,qv))
       end do
    end do

    do j=slo(2),shi(2)+1
       do i=lo(1),hi(1)
          mmtmp8(1:8,i) = matmul(lam(i,j-4:j+3), M8)
          Hg(i,j,iene) = dot_product(mmtmp8(1:8,i), q(i,j-4:j+3,qtemp))
       end do
    end do
    
    do j=slo(2),shi(2)+1
       do i=lo(1),hi(1)
          mmtmp8(1:8,i) = matmul(M8T, q(i,j-4:j+3,qpres))
          Hg(i,j,iene) = Hg(i,j,iene) + dot_product(dpe(i,j-4:j+3), mmtmp8(1:8,i))
       end do
       do i=lo(1),hi(1)
          M8p(:,i,j) = mmtmp8(1:8,i)
       end do
    end do
    
    do n = 1, nspecies

       if (n .eq. iias) cycle  ! inactive speices

       qxn = qx1+n-1
       iryn = iry1+n-1

       do j=slo(2),shi(2)+1
          do i=lo(1),hi(1)
             Hg(i,j,iryn) = dot_product(dpy(i,j-4:j+3,n), M8p(:,i,j))
          end do
       end do

       do j=slo(2),shi(2)+1
          do i=lo(1),hi(1)
             mmtmp8(1:8,i) = matmul(M8T, q(i,j-4:j+3,qxn))
             Hg(i,j,iene) = Hg(i,j,iene) + dot_product(dxe(i,j-4:j+3,n), mmtmp8(1:8,i))
             Hg(i,j,iryn) = Hg(i,j,iryn) &
                  + dot_product(dxy(i,j-4:j+3,n), mmtmp8(1:8,i))
          end do
       end do

    end do
       
    ! add y-direction rhs

    do n=2,iene
       do j=slo(2),shi(2)
          do i=lo(1),hi(1)
             rhs(i,j,n) = rhs(i,j,n) + (Hg(i,j+1,n) - Hg(i,j,n)) * dx2inv(2)
          end do
       end do
    end do

    sumdrY = 0.d0
    do n=iry1,ncons

       if (n.eq.iry_ias) cycle

       do j=slo(2),shi(2)
          do i=lo(1),hi(1)
             sumdry(i,j) = sumdry(i,j) + (Hg(i,j+1,n) - Hg(i,j,n)) * dx2inv(2)
             rhs(i,j,n)  =  rhs(i,j,n) + (Hg(i,j+1,n) - Hg(i,j,n)) * dx2inv(2)
          end do
       end do

    end do

    do j=slo(2),shi(2)
       do i=lo(1),hi(1)
          gradp(i,j) = dxinv(2) * first_deriv_8(q(i,j-4:j+4,qpres))
       end do
    end do
    
    sumryv = 0.d0
    do n = 1, nspecies

       if (n.eq.iias) cycle

       qxn = qx1+n-1
       do j=slo(2),shi(2)
          do i=lo(1),hi(1)
             sumryv(i,j) = sumryv(i,j) + dpy(i,j,n)*gradp(i,j)  &
                  + dxy(i,j,n)*dxinv(2)*first_deriv_8(q(i,j-4:j+4,qxn))
          end do
       end do

    end do

    if (add_v_correction) then

       do n=1,nspecies
          qyn = qy1+n-1
          qhn = qh1+n-1
          iryn = iry1+n-1

          do j=slo(2),shi(2)
             do i=lo(1),hi(1)
                ry_c = q(i,j,qyn)*sumdry(i,j) + sumryv(i,j)*dxinv(2) * &
                     first_deriv_8(q(i,j-4:j+4,qyn))
                ene_c = ry_c*q(i,j,qhn) + q(i,j,qyn)*sumryv(i,j)*dxinv(2)* &
                     first_deriv_8(q(i,j-4:j+4,qhn))
                rhs(i,j,iene) = rhs(i,j,iene) - ene_c
                rhs(i,j,iryn) = rhs(i,j,iryn) - ry_c
             end do
          end do
       end do

    else

       n = iias
       qhn = qh1+n-1
       iryn = iry1+n-1

       do j=slo(2),shi(2)
          do i=lo(1),hi(1)
             ene_c = sumdry(i,j)*q(i,j,qhn) + sumryv(i,j)*dxinv(2)* &
                  first_deriv_8(q(i,j-4:j+4,qhn))
             rhs(i,j,iene) = rhs(i,j,iene) - ene_c
             rhs(i,j,iryn) = rhs(i,j,iryn) - sumdry(i,j)
          end do
       end do

    end if

    ! ------- END y-direction -------
    
    !
    ! lo-x boundary
    !
    if (physbclo(1)) then
       do j=lo(2),hi(2)
          i = lo(1)
          ! use completely right-biased stencil
          mmtmpB = matmul(vsp(i:i+3,j), BRB)
          rhs(i,j,imx) = rhs(i,j,imx)+finlo(1)*dx2inv(1)*dot_product(mmtmpB,q(i:i+3,j,qu))
          
          mmtmpB = matmul(mu(i:i+3,j), BRB)
          rhs(i,j,imy) = rhs(i,j,imy)+foulo(1)*dx2inv(1)*dot_product(mmtmpB,q(i:i+3,j,qv))
          
          mmtmpB = matmul(lam(i:i+3,j), BRB)
          BBp = matmul(BRB, q(i:i+3,j,qpres))
          rhs(i,j,iene) = rhs(i,j,iene) + foulo(1)*dx2inv(1) * &
               ( dot_product(mmtmpB, q(i:i+3,j,qtemp)) &
               + dot_product(      dpe(i:i+3,j), BBp) )
          
          rhstot = 0.d0
          rhsene = 0.d0
          do n = 1, nspecies
             qxn = qx1+n-1
             qyn = qy1+n-1
             
             BBX = matmul(BRB, q(i:i+3,j,qxn))
             
             rhstmp(n) = dot_product(dpy(i:i+3,j,n), BBp) &
                  +      dot_product(dxy(i:i+3,j,n), BBX)
             
             rhsene = rhsene &
                  +      dot_product(dxe(i:i+3,j,n), BBX)
             
             rhstot = rhstot + rhstmp(n)
             Ytmp(n) = q(i,j,qyn)
          end do
          
          do n = 1, nspecies
             rhs(i,j,iry1+n-1) =  rhs(i,j,iry1+n-1) + &
                  foulo(1)*dx2inv(1) * (rhstmp(n) - Ytmp(n)*rhstot)
          end do
          
          do n = 1, nspecies
             qhn = qh1+n-1
             rhsene = rhsene - Ytmp(n) * q(i,j,qhn) * rhstot
          end do
          rhs(i,j,iene) = rhs(i,j,iene) + foulo(1)*dx2inv(1) * rhsene
          
          
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! use 2nd-order stencil for cell lo(1)+1,j,
          do iface=0,1 
             i = lo(1)+1 + iface
             
             mmtmp2 = matmul(vsp(i-1:i,j), M2)
             Hcell(iface,imx) = dot_product(mmtmp2, q(i-1:i,j,qu))
             
             mmtmp2 = matmul(mu(i-1:i,j), M2)
             Hcell(iface,imy) = dot_product(mmtmp2, q(i-1:i,j,qv))
             
             mmtmp2 = matmul(lam(i-1:i,j), M2)
             M2p = matmul(M2,  q(i-1:i,j,qpres))
             Hcell(iface,iene) = dot_product(mmtmp2, q(i-1:i,j,qtemp)) &
                  &            + dot_product(      dpe(i-1:i,j), M2p)
             
             Htot = 0.d0
             do n = 1, nspecies
                qxn = qx1+n-1
                qyn = qy1+n-1
                
                M2X = matmul(M2, q(i-1:i,j,qxn))
                
                Htmp(n) = dot_product(dpy(i-1:i,j,n), M2p) &
                     +    dot_product(dxy(i-1:i,j,n), M2X)
                
                Hcell(iface,iene) = Hcell(iface,iene) &
                     +    dot_product(dxe(i-1:i,j,n), M2X)
                
                Htot = Htot + Htmp(n)
                Ytmp(n) = 0.5d0*(q(i-1,j,qyn) + q(i,j,qyn))
             end do
             
             do n = 1, nspecies
                Hcell(iface,iry1+n-1) = Htmp(n) - Ytmp(n)*Htot
             end do
             
             do n = 1, nspecies
                qhn = qh1+n-1
                hhalf = 0.5d0*(q(i-1,j,qhn) + q(i,j,qhn))
                Hcell(iface,iene) =  Hcell(iface,iene) - Ytmp(n)*Htot*hhalf 
             end do
          end do

          i = lo(1)+1
          do n=2,ncons
             rhs(i,j,n) = rhs(i,j,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(1)
          end do
          
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! use 4th-order stencil for cell lo(1)+2,j,
          do iface=0,1 
             i = lo(1)+2 + iface

             mmtmp4 = matmul(vsp(i-2:i+1,j), M4)
             Hcell(iface,imx) = dot_product(mmtmp4, q(i-2:i+1,j,qu))
             
             mmtmp4 = matmul(mu(i-2:i+1,j), M4)
             Hcell(iface,imy) = dot_product(mmtmp4, q(i-2:i+1,j,qv))
             
             mmtmp4 = matmul(lam(i-2:i+1,j), M4)
             M4p = matmul(M4T,  q(i-2:i+1,j,qpres))
             Hcell(iface,iene) = dot_product(mmtmp4, q(i-2:i+1,j,qtemp))      &
                  &            + dot_product(      dpe(i-2:i+1,j), M4p)
             
             do n = 1, nspecies
                if (n .eq. iias) cycle  ! inactive speices
                
                qxn = qx1+n-1
                iryn = iry1+n-1
                M4X = matmul(M4T, q(i-2:i+1,j,qxn))
                Hcell(iface,iene) = Hcell(iface,iene) &
                     +    dot_product(dxe(i-2:i+1,j,n), M4X)
                Hcell(iface,iryn) = dot_product(dpy(i-2:i+1,j,n), M4p) &
                     + dot_product(dxy(i-2:i+1,j,n), M4X)
             end do
             
          end do
          
          i = lo(1)+2
          
          do n=2,iene
             rhs(i,j,n) = rhs(i,j,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(1)
          end do
          
          sumdrytmp = 0.d0
          do n=iry1,ncons
             if (n.eq.iry_ias) cycle
             sumdrytmp = sumdrytmp + (Hcell(1,n) - Hcell(0,n)) * dx2inv(1)
             rhs(i,j,n)  =  rhs(i,j,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(1)
          end do
          
          gradptmp = dxinv(1) * first_deriv_4(q(i-2:i+2,j,qpres))
          
          sumryvtmp = 0.d0
          do n = 1, nspecies
             if (n.eq.iias) cycle
             qxn = qx1+n-1
             sumryvtmp = sumryvtmp + dpy(i,j,n)*gradptmp  &
                  + dxy(i,j,n)*dxinv(1)*first_deriv_4(q(i-2:i+2,j,qxn))
          end do
          
          if (add_v_correction) then
             
             do n=1,nspecies
                qyn = qy1+n-1
                qhn = qh1+n-1
                iryn = iry1+n-1
                
                ry_c = q(i,j,qyn)*sumdrytmp + sumryvtmp*dxinv(1) * &
                     first_deriv_4(q(i-2:i+2,j,qyn))
                ene_c = ry_c*q(i,j,qhn) + q(i,j,qyn)*sumryvtmp*dxinv(1)* &
                     first_deriv_4(q(i-2:i+2,j,qhn))
                rhs(i,j,iene) = rhs(i,j,iene) - ene_c
                rhs(i,j,iryn) = rhs(i,j,iryn) - ry_c
             end do
             
          else
             
             n = iias
             qhn = qh1+n-1
             iryn = iry1+n-1
             
             ene_c = sumdrytmp*q(i,j,qhn) + sumryvtmp*dxinv(1)* &
                  first_deriv_4(q(i-2:i+2,j,qhn))
             rhs(i,j,iene) = rhs(i,j,iene) - ene_c
             rhs(i,j,iryn) = rhs(i,j,iryn) - sumdrytmp
             
          end if
          
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! use 6th-order stencil for cell lo(1)+3,j,
          do iface=0,1 
             i = lo(1)+3 + iface
             
             mmtmp6 = matmul(vsp(i-3:i+2,j), M6)
             Hcell(iface,imx) = dot_product(mmtmp6, q(i-3:i+2,j,qu))
             
             mmtmp6 = matmul(mu(i-3:i+2,j), M6)
             Hcell(iface,imy) = dot_product(mmtmp6, q(i-3:i+2,j,qv))
             
             mmtmp6 = matmul(lam(i-3:i+2,j), M6)
             M6p = matmul(M6T,  q(i-3:i+2,j,qpres))
             Hcell(iface,iene) = dot_product(mmtmp6, q(i-3:i+2,j,qtemp)) &
                  &            + dot_product(      dpe(i-3:i+2,j), M6p)
             
             do n = 1, nspecies
                if (n .eq. iias) cycle  ! inactive speices
                qxn = qx1+n-1
                iryn = iry1+n-1
                M6X = matmul(M6T, q(i-3:i+2,j,qxn))
                Hcell(iface,iene) = Hcell(iface,iene) &
                     +    dot_product(dxe(i-3:i+2,j,n), M6X)
                Hcell(iface,iryn) = dot_product(dpy(i-3:i+2,j,n), M6p) &
                     +    dot_product(dxy(i-3:i+2,j,n), M6X)
             end do
             
          end do
          
          i = lo(1)+3
          
          do n=2,iene
             rhs(i,j,n) = rhs(i,j,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(1)
          end do
          
          sumdrytmp = 0.d0
          do n=iry1,ncons
             if (n.eq.iry_ias) cycle
             sumdrytmp = sumdrytmp + (Hcell(1,n) - Hcell(0,n)) * dx2inv(1)
             rhs(i,j,n)  =  rhs(i,j,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(1)
          end do
          
          gradptmp = dxinv(1) * first_deriv_6(q(i-3:i+3,j,qpres))
          
          sumryvtmp = 0.d0
          do n = 1, nspecies
             if (n.eq.iias) cycle
             qxn = qx1+n-1
             sumryvtmp = sumryvtmp + dpy(i,j,n)*gradptmp  &
                  + dxy(i,j,n)*dxinv(1)*first_deriv_6(q(i-3:i+3,j,qxn))
          end do
          
          if (add_v_correction) then
             
             do n=1,nspecies
                qyn = qy1+n-1
                qhn = qh1+n-1
                iryn = iry1+n-1
                
                ry_c = q(i,j,qyn)*sumdrytmp + sumryvtmp*dxinv(1) * &
                     first_deriv_6(q(i-3:i+3,j,qyn))
                ene_c = ry_c*q(i,j,qhn) + q(i,j,qyn)*sumryvtmp*dxinv(1)* &
                     first_deriv_6(q(i-3:i+3,j,qhn))
                rhs(i,j,iene) = rhs(i,j,iene) - ene_c
                rhs(i,j,iryn) = rhs(i,j,iryn) - ry_c
             end do
             
          else
             
             n = iias
             qhn = qh1+n-1
             iryn = iry1+n-1
             
             ene_c = sumdrytmp*q(i,j,qhn) + sumryvtmp*dxinv(1)* &
                  first_deriv_6(q(i-3:i+3,j,qhn))
             rhs(i,j,iene) = rhs(i,j,iene) - ene_c
             rhs(i,j,iryn) = rhs(i,j,iryn) - sumdrytmp
             
          end if
          
       end do
    end if

    !
    ! hi-x boundary
    !
    if (physbchi(1)) then
       do j=lo(2),hi(2)
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! use 6th-order stencil for cell hi(1)-3,j,
          do iface=0,1  ! two faces of 
             i = hi(1)-3 + iface
             
             mmtmp6 = matmul(vsp(i-3:i+2,j), M6)
             Hcell(iface,imx) = dot_product(mmtmp6, q(i-3:i+2,j,qu))
             
             mmtmp6 = matmul(mu(i-3:i+2,j), M6)
             Hcell(iface,imy) = dot_product(mmtmp6, q(i-3:i+2,j,qv))
             
             mmtmp6 = matmul(lam(i-3:i+2,j), M6)
             M6p = matmul(M6T,  q(i-3:i+2,j,qpres))
             Hcell(iface,iene) = dot_product(mmtmp6, q(i-3:i+2,j,qtemp)) &
                  &            + dot_product(      dpe(i-3:i+2,j), M6p)
             
             do n = 1, nspecies
                if (n .eq. iias) cycle  ! inactive speices
                qxn = qx1+n-1
                iryn = iry1+n-1
                M6X = matmul(M6T, q(i-3:i+2,j,qxn))
                Hcell(iface,iene) = Hcell(iface,iene) &
                     +    dot_product(dxe(i-3:i+2,j,n), M6X)                
                Hcell(iface,iryn) = dot_product(dpy(i-3:i+2,j,n), M6p) &
                     +    dot_product(dxy(i-3:i+2,j,n), M6X)
             end do
             
          end do
          
          i = hi(1)-3
          
          do n=2,iene
             rhs(i,j,n) = rhs(i,j,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(1)
          end do
          
          sumdrytmp = 0.d0
          do n=iry1,ncons
             if (n.eq.iry_ias) cycle
             sumdrytmp = sumdrytmp + (Hcell(1,n) - Hcell(0,n)) * dx2inv(1)
             rhs(i,j,n)  =  rhs(i,j,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(1)
          end do
          
          gradptmp = dxinv(1) * first_deriv_6(q(i-3:i+3,j,qpres))
          
          sumryvtmp = 0.d0
          do n = 1, nspecies
             if (n.eq.iias) cycle
             qxn = qx1+n-1
             sumryvtmp = sumryvtmp + dpy(i,j,n)*gradptmp  &
                  + dxy(i,j,n)*dxinv(1)*first_deriv_6(q(i-3:i+3,j,qxn))
          end do
          
          if (add_v_correction) then
             
             do n=1,nspecies
                qyn = qy1+n-1
                qhn = qh1+n-1
                iryn = iry1+n-1
                
                ry_c = q(i,j,qyn)*sumdrytmp + sumryvtmp*dxinv(1) * &
                     first_deriv_6(q(i-3:i+3,j,qyn))
                ene_c = ry_c*q(i,j,qhn) + q(i,j,qyn)*sumryvtmp*dxinv(1)* &
                     first_deriv_6(q(i-3:i+3,j,qhn))
                rhs(i,j,iene) = rhs(i,j,iene) - ene_c
                rhs(i,j,iryn) = rhs(i,j,iryn) - ry_c
             end do
             
          else
             
             n = iias
             qhn = qh1+n-1
             iryn = iry1+n-1
             
             ene_c = sumdrytmp*q(i,j,qhn) + sumryvtmp*dxinv(1)* &
                  first_deriv_6(q(i-3:i+3,j,qhn))
             rhs(i,j,iene) = rhs(i,j,iene) - ene_c
             rhs(i,j,iryn) = rhs(i,j,iryn) - sumdrytmp
             
          end if
          
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! use 4th-order stencil for cell hi(1)-2,j,
          do iface=0,1 
             i = hi(1)-2 + iface
             
             mmtmp4 = matmul(vsp(i-2:i+1,j), M4)
             Hcell(iface,imx) = dot_product(mmtmp4, q(i-2:i+1,j,qu))
             
             mmtmp4 = matmul(mu(i-2:i+1,j), M4)
             Hcell(iface,imy) = dot_product(mmtmp4, q(i-2:i+1,j,qv))
             
             mmtmp4 = matmul(lam(i-2:i+1,j), M4)
             M4p = matmul(M4T,  q(i-2:i+1,j,qpres))
             Hcell(iface,iene) = dot_product(mmtmp4, q(i-2:i+1,j,qtemp)) &
                  &            + dot_product(      dpe(i-2:i+1,j), M4p)
             
             do n = 1, nspecies
                if (n .eq. iias) cycle  ! inactive speices
                qxn = qx1+n-1
                iryn = iry1+n-1
                
                M4X = matmul(M4T, q(i-2:i+1,j,qxn))
                Hcell(iface,iene) = Hcell(iface,iene) &
                     +    dot_product(dxe(i-2:i+1,j,n), M4X)
                Hcell(iface,iryn) = dot_product(dpy(i-2:i+1,j,n), M4p) &
                     +    dot_product(dxy(i-2:i+1,j,n), M4X)
             end do
             
          end do
          
          i = hi(1)-2
          
          do n=2,iene
             rhs(i,j,n) = rhs(i,j,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(1)
          end do
          
          sumdrytmp = 0.d0
          do n=iry1,ncons
             if (n.eq.iry_ias) cycle
             sumdrytmp = sumdrytmp + (Hcell(1,n) - Hcell(0,n)) * dx2inv(1)
             rhs(i,j,n)  =  rhs(i,j,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(1)
          end do
          
          gradptmp = dxinv(1) * first_deriv_4(q(i-2:i+2,j,qpres))
          
          sumryvtmp = 0.d0
          do n = 1, nspecies
             if (n.eq.iias) cycle
             qxn = qx1+n-1
             sumryvtmp = sumryvtmp + dpy(i,j,n)*gradptmp  &
                  + dxy(i,j,n)*dxinv(1)*first_deriv_4(q(i-2:i+2,j,qxn))
          end do
          
          if (add_v_correction) then
             
             do n=1,nspecies
                qyn = qy1+n-1
                qhn = qh1+n-1
                iryn = iry1+n-1
                
                ry_c = q(i,j,qyn)*sumdrytmp + sumryvtmp*dxinv(1) * &
                     first_deriv_4(q(i-2:i+2,j,qyn))
                ene_c = ry_c*q(i,j,qhn) + q(i,j,qyn)*sumryvtmp*dxinv(1)* &
                     first_deriv_4(q(i-2:i+2,j,qhn))
                rhs(i,j,iene) = rhs(i,j,iene) - ene_c
                rhs(i,j,iryn) = rhs(i,j,iryn) - ry_c
             end do
             
          else
             
             n = iias
             qhn = qh1+n-1
             iryn = iry1+n-1
             
             ene_c = sumdrytmp*q(i,j,qhn) + sumryvtmp*dxinv(1)* &
                  first_deriv_4(q(i-2:i+2,j,qhn))
             rhs(i,j,iene) = rhs(i,j,iene) - ene_c
             rhs(i,j,iryn) = rhs(i,j,iryn) - sumdrytmp
             
          end if
          
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! use 2nd-order stencil for cell hi(1)-1,j,
          do iface=0,1 
             i = hi(1)-1 + iface
             
             mmtmp2 = matmul(vsp(i-1:i,j), M2)
             Hcell(iface,imx) = dot_product(mmtmp2, q(i-1:i,j,qu))
             
             mmtmp2 = matmul(mu(i-1:i,j), M2)
             Hcell(iface,imy) = dot_product(mmtmp2, q(i-1:i,j,qv))
             
             mmtmp2 = matmul(lam(i-1:i,j), M2)
             M2p = matmul(M2,  q(i-1:i,j,qpres))
             Hcell(iface,iene) = dot_product(mmtmp2, q(i-1:i,j,qtemp)) &
                  &            + dot_product(      dpe(i-1:i,j), M2p)
             
             Htot = 0.d0
             do n = 1, nspecies
                qxn = qx1+n-1
                qyn = qy1+n-1
                
                M2X = matmul(M2, q(i-1:i,j,qxn))
                
                Htmp(n) = dot_product(dpy(i-1:i,j,n), M2p) &
                     +    dot_product(dxy(i-1:i,j,n), M2X)
                
                Hcell(iface,iene) = Hcell(iface,iene) &
                     +    dot_product(dxe(i-1:i,j,n), M2X)
                
                Htot = Htot + Htmp(n)
                Ytmp(n) = 0.5d0*(q(i-1,j,qyn) + q(i,j,qyn))
             end do
             
             do n = 1, nspecies
                Hcell(iface,iry1+n-1) = Htmp(n) - Ytmp(n)*Htot
             end do
             
             do n = 1, nspecies
                qhn = qh1+n-1
                hhalf = 0.5d0*(q(i-1,j,qhn) + q(i,j,qhn))
                Hcell(iface,iene) =  Hcell(iface,iene) - Ytmp(n)*Htot*hhalf
             end do
          end do
          
          i = hi(1)-1
          do n=2,ncons
             rhs(i,j,n) = rhs(i,j,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(1)
          end do
          
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          i = hi(1)
          ! use completely left-biased stencil
          mmtmpB = matmul(vsp(i-3:i,j), BLB)
          rhs(i,j,imx) = rhs(i,j,imx)+finhi(1)*dx2inv(1)*dot_product(mmtmpB,q(i-3:i,j,qu))
          
          mmtmpB= matmul(mu(i-3:i,j), BLB)
          rhs(i,j,imy) = rhs(i,j,imy)+fouhi(1)*dx2inv(1)*dot_product(mmtmpB,q(i-3:i,j,qv))
          
          mmtmpB = matmul(lam(i-3:i,j), BLB)
          BBp = matmul(BLB, q(i-3:i,j,qpres))
          rhs(i,j,iene) = rhs(i,j,iene) + fouhi(1)*dx2inv(1) * &
               ( dot_product(mmtmpB, q(i-3:i,j,qtemp)) &
               + dot_product(      dpe(i-3:i,j), BBp) )
          
          rhstot = 0.d0
          rhsene = 0.d0
          do n = 1, nspecies
             qxn = qx1+n-1
             qyn = qy1+n-1
             
             BBX = matmul(BLB, q(i-3:i,j,qxn))
             
             rhstmp(n) = dot_product(dpy(i-3:i,j,n), BBp) &
                  +      dot_product(dxy(i-3:i,j,n), BBX)
             
             rhsene = rhsene &
                  +      dot_product(dxe(i-3:i,j,n), BBX)
             
             rhstot = rhstot + rhstmp(n)
             Ytmp(n) = q(i,j,qyn)
          end do
          
          do n = 1, nspecies
             rhs(i,j,iry1+n-1) =  rhs(i,j,iry1+n-1) + &
                  fouhi(1)*dx2inv(1) * (rhstmp(n) - Ytmp(n)*rhstot)
          end do
          
          do n = 1, nspecies
             qhn = qh1+n-1
             rhsene = rhsene - Ytmp(n) * q(i,j,qhn) * rhstot
          end do
          rhs(i,j,iene) = rhs(i,j,iene) + fouhi(1)*dx2inv(1) * rhsene
       end do
    end if

    !
    ! lo-y boundary
    !
    if (physbclo(2)) then
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       j = lo(2)
       ! use completely right-biased stencil
       do i=lo(1),hi(1)
          
          mmtmpB = matmul(mu(i,j:j+3), BRB)
          rhs(i,j,imx) = rhs(i,j,imx)+foulo(2)*dx2inv(2)*dot_product(mmtmpB,q(i,j:j+3,qu))
          
          mmtmpB = matmul(vsp(i,j:j+3), BRB)
          rhs(i,j,imy) = rhs(i,j,imy)+finlo(2)*dx2inv(2)*dot_product(mmtmpB,q(i,j:j+3,qv))
          
          mmtmpB = matmul(lam(i,j:j+3), BRB)
          BBp = matmul(BRB, q(i,j:j+3,qpres))
          rhs(i,j,iene) = rhs(i,j,iene) + foulo(2)*dx2inv(2) * &
               ( dot_product(mmtmpB, q(i,j:j+3,qtemp)) &
               + dot_product(      dpe(i,j:j+3), BBp) )
          
          rhstot = 0.d0
          rhsene = 0.d0
          do n = 1, nspecies
             qxn = qx1+n-1
             qyn = qy1+n-1
             
             BBX = matmul(BRB, q(i,j:j+3,qxn))
             
             rhstmp(n) = dot_product(dpy(i,j:j+3,n), BBp) &
                  +      dot_product(dxy(i,j:j+3,n), BBX)
             
             rhsene = rhsene &
                  +      dot_product(dxe(i,j:j+3,n), BBX)
             
             rhstot = rhstot + rhstmp(n)
             Ytmp(n) = q(i,j,qyn)
          end do
          
          do n = 1, nspecies
             rhs(i,j,iry1+n-1) =  rhs(i,j,iry1+n-1) + &
                  foulo(2)*dx2inv(2) * (rhstmp(n) - Ytmp(n)*rhstot)
          end do
          
          do n = 1, nspecies
             qhn = qh1+n-1
             rhsene = rhsene - Ytmp(n) * q(i,j,qhn) * rhstot
          end do
          rhs(i,j,iene) = rhs(i,j,iene) + foulo(2)*dx2inv(2) * rhsene
       end do
       
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! use 2nd-order stencil for cell i,lo(2)+1
       do i=lo(1),hi(1)
          do iface=0,1 
             j = lo(2)+1 + iface
             
             mmtmp2 = matmul(mu(i,j-1:j), M2)
             Hcell(iface,imx) = dot_product(mmtmp2, q(i,j-1:j,qu))
             
             mmtmp2 = matmul(vsp(i,j-1:j), M2)
             Hcell(iface,imy) = dot_product(mmtmp2, q(i,j-1:j,qv))
             
             mmtmp2 = matmul(lam(i,j-1:j), M2)
             M2p = matmul(M2,  q(i,j-1:j,qpres))
             Hcell(iface,iene) = dot_product(mmtmp2, q(i,j-1:j,qtemp)) &
                  &            + dot_product(      dpe(i,j-1:j), M2p)
             
             Htot = 0.d0
             do n = 1, nspecies
                qxn = qx1+n-1
                qyn = qy1+n-1
                
                M2X = matmul(M2, q(i,j-1:j,qxn))
                
                Htmp(n) = dot_product(dpy(i,j-1:j,n), M2p) &
                     +    dot_product(dxy(i,j-1:j,n), M2X)
                
                Hcell(iface,iene) = Hcell(iface,iene) &
                     +    dot_product(dxe(i,j-1:j,n), M2X)
                
                Htot = Htot + Htmp(n)
                Ytmp(n) = 0.5d0*(q(i,j-1,qyn) + q(i,j,qyn))
             end do
             
             do n = 1, nspecies
                Hcell(iface,iry1+n-1) = Htmp(n) - Ytmp(n)*Htot
             end do
             
             do n = 1, nspecies
                qhn = qh1+n-1
                hhalf = 0.5d0*(q(i,j-1,qhn) + q(i,j,qhn))
                Hcell(iface,iene) =  Hcell(iface,iene) - Ytmp(n)*Htot*hhalf
             end do
          end do
          
          j = lo(2)+1
          do n=2,ncons
             rhs(i,j,n) = rhs(i,j,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(2)
          end do
       end do
       
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! use 4th-order stencil for cell i,lo(2)+2
       do i=lo(1),hi(1)
          do iface=0,1 
             j = lo(2)+2 + iface
             
             mmtmp4 = matmul(mu(i,j-2:j+1), M4)
             Hcell(iface,imx) = dot_product(mmtmp4, q(i,j-2:j+1,qu))
             
             mmtmp4 = matmul(vsp(i,j-2:j+1), M4)
             Hcell(iface,imy) = dot_product(mmtmp4, q(i,j-2:j+1,qv))
             
             mmtmp4 = matmul(lam(i,j-2:j+1), M4)
             M4p = matmul(M4T,  q(i,j-2:j+1,qpres))
             Hcell(iface,iene) = dot_product(mmtmp4, q(i,j-2:j+1,qtemp)) &
                  &            + dot_product(      dpe(i,j-2:j+1), M4p)
             
             do n = 1, nspecies
                if (n .eq. iias) cycle  ! inactive speices
                
                qxn = qx1+n-1
                iryn = iry1+n-1
                M4X = matmul(M4T, q(i,j-2:j+1,qxn))
                Hcell(iface,iene) = Hcell(iface,iene) &
                     +    dot_product(dxe(i,j-2:j+1,n), M4X)
                Hcell(iface,iryn) = dot_product(dpy(i,j-2:j+1,n), M4p) &
                     +    dot_product(dxy(i,j-2:j+1,n), M4X)
             end do
             
          end do
          
          j = lo(2)+2
          
          do n=2,iene
             rhs(i,j,n) = rhs(i,j,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(2)
          end do
          
          sumdrytmp = 0.d0
          do n=iry1,ncons
             if (n.eq.iry_ias) cycle
             sumdrytmp = sumdrytmp + (Hcell(1,n) - Hcell(0,n)) * dx2inv(2)
             rhs(i,j,n)  =  rhs(i,j,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(2)
          end do

          gradptmp = dxinv(2) * first_deriv_4(q(i,j-2:j+2,qpres))

          sumryvtmp = 0.d0
          do n = 1, nspecies
             if (n.eq.iias) cycle
             qxn = qx1+n-1
             sumryvtmp = sumryvtmp + dpy(i,j,n)*gradptmp  &
                  + dxy(i,j,n)*dxinv(2)*first_deriv_4(q(i,j-2:j+2,qxn))
          end do
          
          if (add_v_correction) then
             
             do n=1,nspecies
                qyn = qy1+n-1
                qhn = qh1+n-1
                iryn = iry1+n-1
                
                ry_c = q(i,j,qyn)*sumdrytmp + sumryvtmp*dxinv(2) * &
                     first_deriv_4(q(i,j-2:j+2,qyn))
                ene_c = ry_c*q(i,j,qhn) + q(i,j,qyn)*sumryvtmp*dxinv(2)* &
                     first_deriv_4(q(i,j-2:j+2,qhn))
                rhs(i,j,iene) = rhs(i,j,iene) - ene_c
                rhs(i,j,iryn) = rhs(i,j,iryn) - ry_c
             end do
             
          else
             
             n = iias
             qhn = qh1+n-1
             iryn = iry1+n-1
             
             ene_c = sumdrytmp*q(i,j,qhn) + sumryvtmp*dxinv(2)* &
                  first_deriv_4(q(i,j-2:j+2,qhn))
             rhs(i,j,iene) = rhs(i,j,iene) - ene_c
             rhs(i,j,iryn) = rhs(i,j,iryn) - sumdrytmp
             
          end if
       end do
       
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! use 6th-order stencil for cell i,lo(2)+3
       do i=lo(1),hi(1)
          do iface=0,1 
             j = lo(2)+3 + iface
             
             mmtmp6 = matmul(mu(i,j-3:j+2), M6)
             Hcell(iface,imx) = dot_product(mmtmp6, q(i,j-3:j+2,qu))
             
             mmtmp6 = matmul(vsp(i,j-3:j+2), M6)
             Hcell(iface,imy) = dot_product(mmtmp6, q(i,j-3:j+2,qv))
             
             mmtmp6 = matmul(lam(i,j-3:j+2), M6)
             M6p = matmul(M6T,  q(i,j-3:j+2,qpres))
             Hcell(iface,iene) = dot_product(mmtmp6, q(i,j-3:j+2,qtemp))      &
                  &            + dot_product(      dpe(i,j-3:j+2), M6p)
             
             do n = 1, nspecies
                if (n .eq. iias) cycle  ! inactive speices
                qxn = qx1+n-1
                iryn = iry1+n-1
                M6X = matmul(M6T, q(i,j-3:j+2,qxn))
                Hcell(iface,iene) = Hcell(iface,iene) &
                     +    dot_product(dxe(i,j-3:j+2,n), M6X)
                Hcell(iface,iryn) = dot_product(dpy(i,j-3:j+2,n), M6p) &
                     +    dot_product(dxy(i,j-3:j+2,n), M6X)
             end do
             
          end do
          
          j = lo(2)+3
          
          do n=2,iene
             rhs(i,j,n) = rhs(i,j,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(2)
          end do
          
          sumdrytmp = 0.d0
          do n=iry1,ncons
             if (n.eq.iry_ias) cycle
             sumdrytmp = sumdrytmp + (Hcell(1,n) - Hcell(0,n)) * dx2inv(2)
             rhs(i,j,n)  =  rhs(i,j,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(2)
          end do
          
          gradptmp = dxinv(2) * first_deriv_6(q(i,j-3:j+3,qpres))
          
          sumryvtmp = 0.d0
          do n = 1, nspecies
             if (n.eq.iias) cycle
             qxn = qx1+n-1
             sumryvtmp = sumryvtmp + dpy(i,j,n)*gradptmp  &
                  + dxy(i,j,n)*dxinv(2)*first_deriv_6(q(i,j-3:j+3,qxn))
          end do
          
          if (add_v_correction) then
             
             do n=1,nspecies
                qyn = qy1+n-1
                qhn = qh1+n-1
                iryn = iry1+n-1
                
                ry_c = q(i,j,qyn)*sumdrytmp + sumryvtmp*dxinv(2) * &
                     first_deriv_6(q(i,j-3:j+3,qyn))
                ene_c = ry_c*q(i,j,qhn) + q(i,j,qyn)*sumryvtmp*dxinv(2)* &
                     first_deriv_6(q(i,j-3:j+3,qhn))
                rhs(i,j,iene) = rhs(i,j,iene) - ene_c
                rhs(i,j,iryn) = rhs(i,j,iryn) - ry_c
             end do
             
          else
             
             n = iias
             qhn = qh1+n-1
             iryn = iry1+n-1
             
             ene_c = sumdrytmp*q(i,j,qhn) + sumryvtmp*dxinv(2)* &
                  first_deriv_6(q(i,j-3:j+3,qhn))
             rhs(i,j,iene) = rhs(i,j,iene) - ene_c
             rhs(i,j,iryn) = rhs(i,j,iryn) - sumdrytmp
             
          end if
       end do
    end if

    !
    ! hi-y boundary
    !
    if (physbchi(2)) then
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! use 6th-order stencil for cell i,hi(2)-3
       do i=lo(1),hi(1)
          do iface=0,1 
             j = hi(2)-3 + iface
             
             mmtmp6 = matmul(mu(i,j-3:j+2), M6)
             Hcell(iface,imx) = dot_product(mmtmp6, q(i,j-3:j+2,qu))
             
             mmtmp6 = matmul(vsp(i,j-3:j+2), M6)
             Hcell(iface,imy) = dot_product(mmtmp6, q(i,j-3:j+2,qv))
             
             mmtmp6 = matmul(lam(i,j-3:j+2), M6)
             M6p = matmul(M6T,  q(i,j-3:j+2,qpres))
             Hcell(iface,iene) = dot_product(mmtmp6, q(i,j-3:j+2,qtemp))      &
                  &            + dot_product(      dpe(i,j-3:j+2), M6p)
             
             do n = 1, nspecies
                if (n .eq. iias) cycle  ! inactive speices
                qxn = qx1+n-1
                iryn = iry1+n-1
                M6X = matmul(M6T, q(i,j-3:j+2,qxn))
                Hcell(iface,iene) = Hcell(iface,iene) &
                     +    dot_product(dxe(i,j-3:j+2,n), M6X)
                Hcell(iface,iryn) = dot_product(dpy(i,j-3:j+2,n), M6p) &
                     +    dot_product(dxy(i,j-3:j+2,n), M6X)
             end do
             
          end do
          
          j = hi(2)-3
          
          do n=2,iene
             rhs(i,j,n) = rhs(i,j,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(2)
          end do
          
          sumdrytmp = 0.d0
          do n=iry1,ncons
             if (n.eq.iry_ias) cycle
             sumdrytmp = sumdrytmp + (Hcell(1,n) - Hcell(0,n)) * dx2inv(2)
             rhs(i,j,n)  =  rhs(i,j,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(2)
          end do
          
          gradptmp = dxinv(2) * first_deriv_6(q(i,j-3:j+3,qpres))
          
          sumryvtmp = 0.d0
          do n = 1, nspecies
             if (n.eq.iias) cycle
             qxn = qx1+n-1
             sumryvtmp = sumryvtmp + dpy(i,j,n)*gradptmp  &
                  + dxy(i,j,n)*dxinv(2)*first_deriv_6(q(i,j-3:j+3,qxn))
          end do
          
          if (add_v_correction) then
             
             do n=1,nspecies
                qyn = qy1+n-1
                qhn = qh1+n-1
                iryn = iry1+n-1
                
                ry_c = q(i,j,qyn)*sumdrytmp + sumryvtmp*dxinv(2) * &
                     first_deriv_6(q(i,j-3:j+3,qyn))
                ene_c = ry_c*q(i,j,qhn) + q(i,j,qyn)*sumryvtmp*dxinv(2)* &
                     first_deriv_6(q(i,j-3:j+3,qhn))
                rhs(i,j,iene) = rhs(i,j,iene) - ene_c
                rhs(i,j,iryn) = rhs(i,j,iryn) - ry_c
             end do
             
          else
             
             n = iias
             qhn = qh1+n-1
             iryn = iry1+n-1
             
             ene_c = sumdrytmp*q(i,j,qhn) + sumryvtmp*dxinv(2)* &
                  first_deriv_6(q(i,j-3:j+3,qhn))
             rhs(i,j,iene) = rhs(i,j,iene) - ene_c
             rhs(i,j,iryn) = rhs(i,j,iryn) - sumdrytmp
             
          end if
       end do
       
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! use 4th-order stencil for cell i,hi(2)-2
       do i=lo(1),hi(1)
          do iface=0,1 
             j = hi(2)-2 + iface
             
             mmtmp4 = matmul(mu(i,j-2:j+1), M4)
             Hcell(iface,imx) = dot_product(mmtmp4, q(i,j-2:j+1,qu))
             
             mmtmp4 = matmul(vsp(i,j-2:j+1), M4)
             Hcell(iface,imy) = dot_product(mmtmp4, q(i,j-2:j+1,qv))
             
             mmtmp4 = matmul(lam(i,j-2:j+1), M4)
             M4p = matmul(M4T,  q(i,j-2:j+1,qpres))
             Hcell(iface,iene) = dot_product(mmtmp4, q(i,j-2:j+1,qtemp)) &
                  &            + dot_product(      dpe(i,j-2:j+1), M4p)
             
             do n = 1, nspecies
                if (n .eq. iias) cycle  ! inactive speices
                
                qxn = qx1+n-1
                iryn = iry1+n-1
                M4X = matmul(M4T, q(i,j-2:j+1,qxn))
                Hcell(iface,iene) = Hcell(iface,iene) &
                     +    dot_product(dxe(i,j-2:j+1,n), M4X)
                Hcell(iface,iryn) = dot_product(dpy(i,j-2:j+1,n), M4p) &
                     +    dot_product(dxy(i,j-2:j+1,n), M4X)
             end do
             
          end do
          
          j = hi(2)-2
          
          do n=2,iene
             rhs(i,j,n) = rhs(i,j,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(2)
          end do
          
          sumdrytmp = 0.d0
          do n=iry1,ncons
             if (n.eq.iry_ias) cycle
             sumdrytmp = sumdrytmp + (Hcell(1,n) - Hcell(0,n)) * dx2inv(2)
             rhs(i,j,n)  =  rhs(i,j,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(2)
          end do
          
          gradptmp = dxinv(2) * first_deriv_4(q(i,j-2:j+2,qpres))
          
          sumryvtmp = 0.d0
          do n = 1, nspecies
             if (n.eq.iias) cycle
             qxn = qx1+n-1
             sumryvtmp = sumryvtmp + dpy(i,j,n)*gradptmp  &
                  + dxy(i,j,n)*dxinv(2)*first_deriv_4(q(i,j-2:j+2,qxn))
          end do
          
          if (add_v_correction) then
             
             do n=1,nspecies
                qyn = qy1+n-1
                qhn = qh1+n-1
                iryn = iry1+n-1
                
                ry_c = q(i,j,qyn)*sumdrytmp + sumryvtmp*dxinv(2) * &
                     first_deriv_4(q(i,j-2:j+2,qyn))
                ene_c = ry_c*q(i,j,qhn) + q(i,j,qyn)*sumryvtmp*dxinv(2)* &
                     first_deriv_4(q(i,j-2:j+2,qhn))
                rhs(i,j,iene) = rhs(i,j,iene) - ene_c
                rhs(i,j,iryn) = rhs(i,j,iryn) - ry_c
             end do
             
          else
             
             n = iias
             qhn = qh1+n-1
             iryn = iry1+n-1
             
             ene_c = sumdrytmp*q(i,j,qhn) + sumryvtmp*dxinv(2)* &
                  first_deriv_4(q(i,j-2:j+2,qhn))
             rhs(i,j,iene) = rhs(i,j,iene) - ene_c
             rhs(i,j,iryn) = rhs(i,j,iryn) - sumdrytmp
             
          end if
       end do

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! use 2nd-order stencil for cell i,hi(2)-1
       do i=lo(1),hi(1)
          do iface=0,1 
             j = hi(2)-1 + iface
             
             mmtmp2 = matmul(mu(i,j-1:j), M2)
             Hcell(iface,imx) = dot_product(mmtmp2, q(i,j-1:j,qu))
             
             mmtmp2 = matmul(vsp(i,j-1:j), M2)
             Hcell(iface,imy) = dot_product(mmtmp2, q(i,j-1:j,qv))
             
             mmtmp2 = matmul(lam(i,j-1:j), M2)
             M2p = matmul(M2,  q(i,j-1:j,qpres))
             Hcell(iface,iene) = dot_product(mmtmp2, q(i,j-1:j,qtemp))      &
                  &            + dot_product(      dpe(i,j-1:j), M2p)
             
             Htot = 0.d0
             do n = 1, nspecies
                qxn = qx1+n-1
                qyn = qy1+n-1
                
                M2X = matmul(M2, q(i,j-1:j,qxn))
                
                Htmp(n) = dot_product(dpy(i,j-1:j,n), M2p) &
                     +    dot_product(dxy(i,j-1:j,n), M2X)
                
                Hcell(iface,iene) = Hcell(iface,iene) &
                     +    dot_product(dxe(i,j-1:j,n), M2X)
                
                Htot = Htot + Htmp(n)
                Ytmp(n) = 0.5d0*(q(i,j-1,qyn) + q(i,j,qyn))
             end do
             
             do n = 1, nspecies
                Hcell(iface,iry1+n-1) = Htmp(n) - Ytmp(n)*Htot
             end do
             
             do n = 1, nspecies
                qhn = qh1+n-1
                hhalf = 0.5d0*(q(i,j-1,qhn) + q(i,j,qhn))
                Hcell(iface,iene) =  Hcell(iface,iene) - Ytmp(n)*Htot*hhalf
             end do
          end do
          
          j = hi(2)-1
          do n=2,ncons
             rhs(i,j,n) = rhs(i,j,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(2)
          end do
       end do
       
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       j = hi(2)
       ! use completely right-biased stencil
       do i=lo(1),hi(1)

          mmtmpB = matmul(mu(i,j-3:j), BLB)
          rhs(i,j,imx) = rhs(i,j,imx)+fouhi(2)*dx2inv(2)*dot_product(mmtmpB,q(i,j-3:j,qu))
          
          mmtmpB = matmul(vsp(i,j-3:j), BLB)
          rhs(i,j,imy) = rhs(i,j,imy)+finhi(2)*dx2inv(2)*dot_product(mmtmpB,q(i,j-3:j,qv))
          
          mmtmpB = matmul(lam(i,j-3:j), BLB)
          BBp = matmul(BLB, q(i,j-3:j,qpres))
          rhs(i,j,iene) = rhs(i,j,iene) + fouhi(2)*dx2inv(2) * &
               ( dot_product(mmtmpB, q(i,j-3:j,qtemp)) &
               + dot_product(      dpe(i,j-3:j), BBp) )
          
          rhstot = 0.d0
          rhsene = 0.d0
          do n = 1, nspecies
             qxn = qx1+n-1
             qyn = qy1+n-1
             
             BBX = matmul(BLB, q(i,j-3:j,qxn))
             
             rhstmp(n) = dot_product(dpy(i,j-3:j,n), BBp) &
                  +      dot_product(dxy(i,j-3:j,n), BBX)
             
             rhsene = rhsene &
                  +      dot_product(dxe(i,j-3:j,n), BBX)
             
             rhstot = rhstot + rhstmp(n)
             Ytmp(n) = q(i,j,qyn)
          end do
          
          do n = 1, nspecies
             rhs(i,j,iry1+n-1) =  rhs(i,j,iry1+n-1) + &
                  fouhi(2)*dx2inv(2) * (rhstmp(n) - Ytmp(n)*rhstot)
          end do
          
          do n = 1, nspecies
             qhn = qh1+n-1
             rhsene = rhsene - Ytmp(n) * q(i,j,qhn) * rhstot
          end do
          rhs(i,j,iene) = rhs(i,j,iene) + fouhi(2)*dx2inv(2) * rhsene
       end do
    end if

    !
    ! add kinetic energy
    !
    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          rhs(i,j,iene) = rhs(i,j,iene) &
               + rhs(i,j,imx)*q(i,j,qu) &
               + rhs(i,j,imy)*q(i,j,qv)
       end do
    end do

    deallocate(Hg,dpy,dxe,dpe,vsp,M8p)
    deallocate(sumdrY,sumryv,gradp)

  end subroutine diffterm_2


  subroutine chemterm_2d(lo,hi,q,qlo,qhi,up,uplo,uphi,dt)
    use probin_module, only : use_vode
    use vode_module, only : verbose, itol, rtol, atol, vode_MF=>MF, &
         voderwork, vodeiwork, lvoderwork, lvodeiwork, voderpar, vodeipar

    double precision, intent(in) :: dt
    integer,         intent(in):: lo(2),hi(2),qlo(2),qhi(2),uplo(2),uphi(2)
    double precision,intent(in):: q ( qlo(1): qhi(1), qlo(2): qhi(2),nprim)
    double precision           :: up(uplo(1):uphi(1),uplo(2):uphi(2),ncons)

    integer :: iwrk, i,j,n,np
    double precision :: Yt(lo(1):hi(1),nspecies), wdot(lo(1):hi(1),nspecies), rwrk
    double precision :: YTvode(nspecies+1), time, dtinv

    external f_jac, f_rhs, dvode

    ! vode stuff
    integer, parameter :: itask=1, iopt=1
    integer :: MF, istate, ifail

    if (use_vode) then

       dtinv = 1.d0/dt

       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             voderpar(1) = q(i,j,qrho)

             YTvode(1:nspecies) = q(i,j,qy1:qy1+nspecies-1)
             YTvode(nspecies+1) = q(i,j,qtemp)
       
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
                up(i,j,iry1+n-1) = dtinv*q(i,j,qrho)*(YTvode(n)-q(i,j,qy1+n-1))
             end do

          end do
       end do

    else

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
                up(i,j,iry1+n-1) = wdot(i,n) * molecular_weight(n)
             end do
          end do
          
       end do
       
    end if

  end subroutine chemterm_2d


  subroutine comp_courno_2d(lo,hi,dx,Q,qlo,qhi,courno)
    integer, intent(in) :: lo(2), hi(2), qlo(2), qhi(2)
    double precision, intent(in) :: dx(2)
    double precision, intent(in) :: q(qlo(1):qhi(1),qlo(2):qhi(2),nprim)
    double precision, intent(inout) :: courno

    integer :: i,j, iwrk
    double precision :: dxinv(2), c, rwrk, Cv, Cp
    double precision :: Tt, X(nspecies), gamma
    double precision :: courx, coury

    do i=1,2
       dxinv(i) = 1.0d0 / dx(i)
    end do

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          
          Tt = q(i,j,qtemp)
          X  = q(i,j,qx1:qx1+nspecies-1)
          call ckcvbl(Tt, X, iwrk, rwrk, Cv)
          Cp = Cv + Ru
          gamma = Cp / Cv
          c = sqrt(gamma*q(i,j,qpres)/q(i,j,qrho))
          
          courx = (c+abs(q(i,j,qu))) * dxinv(1)
          coury = (c+abs(q(i,j,qv))) * dxinv(2)
             
          courno = max( courx, coury , courno )

       end do
    end do

  end subroutine comp_courno_2d

end module kernels_2d_module
