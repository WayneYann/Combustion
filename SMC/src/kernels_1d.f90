module kernels_1d_module
  use bc_module
  use chemistry_module, only : nspecies, molecular_weight, inv_mwt, Ru
  use derivative_stencil_module, only : stencil_ng, first_deriv_8, first_deriv_6, &
       first_deriv_4, first_deriv_l3, first_deriv_r3, first_deriv_rb, first_deriv_lb, &
       M8, M8T, M6, M6T, M4, M4T, M2, BRB, BLB, D8, D6, D4
  use variables_module
  implicit none

  private

  public :: hypterm_1d, narrow_diffterm_1d, chemterm_1d, comp_courno_1d

contains

  subroutine hypterm_1d (lo,hi,dx,cons,clo,chi,q,qlo,qhi,rhs_g,rlo,rhi,dlo_g,dhi_g,bclo,bchi)

    integer,         intent(in):: dlo_g(1),dhi_g(1),bclo(1),bchi(1)
    integer,         intent(in):: lo(1),hi(1),clo(1),chi(1),qlo(1),qhi(1),rlo(1),rhi(1)
    double precision,intent(in):: dx(1)
    double precision,intent(in):: cons(clo(1):chi(1),ncons)
    double precision,intent(in)::    q(qlo(1):qhi(1),nprim)
    double precision           ::rhs_g(rlo(1):rhi(1),ncons)

    integer          :: i,n
    double precision :: dxinv(1)
    double precision :: un(-4:4)
    integer :: slo(1), shi(1), dlo(1), dhi(1)
    
    double precision, allocatable :: tmpx(:)
    double precision, allocatable :: rhs(:,:)

    logical :: physbclo(1), physbchi(1)

    ! Only the region bounded by [dlo_g,dhi_g] contains good data.
    ! [slo,shi] will be safe for 8th-order stencil
    dlo(1) = max(lo(1)-stencil_ng, dlo_g(1))
    dhi(1) = min(hi(1)+stencil_ng, dhi_g(1))
    slo(1) = dlo(1) + stencil_ng
    shi(1) = dhi(1) - stencil_ng
    
    if (dlo(1) .eq. lo(1)) then
       physbclo(1) = .true.
    else
       physbclo(1) = .false.
    end if
    
    if (dhi(1) .eq. hi(1)) then
       physbchi(1) = .true.
    else
       physbchi(1) = .false.
    end if

    dxinv(1) = 1.0d0 / dx(1)

    allocate(tmpx(dlo(1):dhi(1)))
    allocate(rhs(lo(1):hi(1),ncons))
    rhs = 0.d0

    do i=slo(1),shi(1)
!expand       rhs(i,irho) = rhs(i,irho) - dxinv(1) * &
!expand            first_deriv_8( cons(i-4:i+4,imx) ) 
       rhs(i,irho) = rhs(i,irho) - dxinv(1) * &
            ( D8(1)*(cons(i+1,imx)-cons(i-1,imx)) &
            + D8(2)*(cons(i+2,imx)-cons(i-2,imx)) &
            + D8(3)*(cons(i+3,imx)-cons(i-3,imx)) &
            + D8(4)*(cons(i+4,imx)-cons(i-4,imx)) )
    end do

    do i=dlo(1),dhi(1)
       tmpx(i) = cons(i,imx)*q(i,qu)+q(i,qpres)
    end do
    do i=slo(1),shi(1)
!expand       rhs(i,imx) = rhs(i,imx) - dxinv(1) * first_deriv_8(tmpx(i-4:i+4))
       rhs(i,imx) = rhs(i,imx) - dxinv(1) * &
            ( D8(1)*(tmpx(i+1)-tmpx(i-1)) &
            + D8(2)*(tmpx(i+2)-tmpx(i-2)) &
            + D8(3)*(tmpx(i+3)-tmpx(i-3)) &
            + D8(4)*(tmpx(i+4)-tmpx(i-4)) )
    end do
       
    do i=dlo(1),dhi(1)
       tmpx(i) = (cons(i,iene)+q(i,qpres))*q(i,qu)
    end do
    do i=slo(1),shi(1)
!expand       rhs(i,iene) = rhs(i,iene) - dxinv(1) * first_deriv_8(tmpx(i-4:i+4))
       rhs(i,iene) = rhs(i,iene) - dxinv(1) * &
            ( D8(1)*(tmpx(i+1)-tmpx(i-1)) &
            + D8(2)*(tmpx(i+2)-tmpx(i-2)) &
            + D8(3)*(tmpx(i+3)-tmpx(i-3)) &
            + D8(4)*(tmpx(i+4)-tmpx(i-4)) )
    end do

    do n = iry1, iry1+nspecies-1
       do i=dlo(1),dhi(1)
          tmpx(i) = cons(i,n)*q(i,qu)
       end do
       do i=slo(1),shi(1)
!expand          rhs(i,n) = rhs(i,n) - dxinv(1) * first_deriv_8(tmpx(i-4:i+4))
          rhs(i,n) = rhs(i,n) - dxinv(1) * &
               ( D8(1)*(tmpx(i+1)-tmpx(i-1)) &
               + D8(2)*(tmpx(i+2)-tmpx(i-2)) &
               + D8(3)*(tmpx(i+3)-tmpx(i-3)) &
               + D8(4)*(tmpx(i+4)-tmpx(i-4)) )
       end do
    enddo

    !
    ! ----- lo-x boundary -----
    !
    if (physbclo(1)) then 
          
       ! if (bclo(1) .eq. WALL???) then
       !    i = lo(1)
       !    ! use completely right-biased stencil
       ! end if
          
       i = lo(1)+1
       ! use 3rd-order slightly right-biased stencil
       
       un(-1:2) = q(i-1:i+2,qu)
          
       rhs(i,irho) = rhs(i,irho) - dxinv(1) * &
            first_deriv_r3( cons(i-1:i+2,imx) ) 
          
       rhs(i,imx) = rhs(i,imx) - dxinv(1) * &
            first_deriv_r3( cons(i-1:i+2,imx)*un(-1:2)+q(i-1:i+2,qpres) )
          
       rhs(i,iene) = rhs(i,iene) - dxinv(1) * &
            first_deriv_r3( (cons(i-1:i+2,iene)+q(i-1:i+2,qpres))*un(-1:2) )

       i = lo(1)+2
       ! use 4th-order stencil

       un(-2:2) = q(i-2:i+2,qu)
          
!expand       rhs(i,irho) = rhs(i,irho) - dxinv(1) * &
!expand            first_deriv_4( cons(i-2:i+2,imx) ) 
       rhs(i,irho) = rhs(i,irho) - dxinv(1) * &
            ( D4(1)*(cons(i+1,imx)-cons(i-1,imx)) &
            + D4(2)*(cons(i+2,imx)-cons(i-2,imx)) )
          
!expand       rhs(i,imx) = rhs(i,imx) - dxinv(1) * &
!expand            first_deriv_4( cons(i-2:i+2,imx)*un(-2:2)+q(i-2:i+2,qpres) )
       rhs(i,imx) = rhs(i,imx) - dxinv(1) * &
            ( D4(1)*((cons(i+1,imx)*un(1)+q(i+1,qpres))-(cons(i-1,imx)*un(-1)+q(i-1,qpres))) &
            + D4(2)*((cons(i+2,imx)*un(2)+q(i+2,qpres))-(cons(i-2,imx)*un(-2)+q(i-2,qpres))) )
          
!expand       rhs(i,iene) = rhs(i,iene) - dxinv(1) * &
!expand            first_deriv_4( (cons(i-2:i+2,iene)+q(i-2:i+2,qpres))*un(-2:2) )
       rhs(i,iene) = rhs(i,iene) - dxinv(1) * &
            ( D4(1)*((cons(i+1,iene)+q(i+1,qpres))*un(1)-(cons(i-1,iene)+q(i-1,qpres))*un(-1)) &
            + D4(2)*((cons(i+2,iene)+q(i+2,qpres))*un(2)-(cons(i-2,iene)+q(i-2,qpres))*un(-2)) )
       
       i = lo(1)+3
       ! use 6th-order stencil
       
       un(-3:3) = q(i-3:i+3,qu)
          
!expand       rhs(i,irho) = rhs(i,irho) - dxinv(1) * &
!expand            first_deriv_6( cons(i-3:i+3,imx) ) 
       rhs(i,irho) = rhs(i,irho) - dxinv(1) * &
            ( D6(1)*(cons(i+1,imx)-cons(i-1,imx)) &
            + D6(2)*(cons(i+2,imx)-cons(i-2,imx)) &
            + D6(3)*(cons(i+3,imx)-cons(i-3,imx)) )
       
!expand       rhs(i,imx) = rhs(i,imx) - dxinv(1) * &
!expand               first_deriv_6( cons(i-3:i+3,imx)*un(-3:3)+q(i-3:i+3,qpres) )
       rhs(i,imx) = rhs(i,imx) - dxinv(1) * &
            ( D6(1)*((cons(i+1,imx)*un(1)+q(i+1,qpres))-(cons(i-1,imx)*un(-1)+q(i-1,qpres))) &
            + D6(2)*((cons(i+2,imx)*un(2)+q(i+2,qpres))-(cons(i-2,imx)*un(-2)+q(i-2,qpres))) &
            + D6(3)*((cons(i+3,imx)*un(3)+q(i+3,qpres))-(cons(i-3,imx)*un(-3)+q(i-3,qpres))) )

!expand       rhs(i,iene) = rhs(i,iene) - dxinv(1) * &
!expand            first_deriv_6( (cons(i-3:i+3,iene)+q(i-3:i+3,qpres))*un(-3:3) )
       rhs(i,iene) = rhs(i,iene) - dxinv(1) * &
            ( D6(1)*((cons(i+1,iene)+q(i+1,qpres))*un(1)-(cons(i-1,iene)+q(i-1,qpres))*un(-1)) &
            + D6(2)*((cons(i+2,iene)+q(i+2,qpres))*un(2)-(cons(i-2,iene)+q(i-2,qpres))*un(-2)) &
            + D6(3)*((cons(i+3,iene)+q(i+3,qpres))*un(3)-(cons(i-3,iene)+q(i-3,qpres))*un(-3)) )
       
       do n = iry1, iry1+nspecies-1
          ! if (bclo(1) .eq. WALL???) then
          !    i = lo(1)
          !    ! use completely right-biased stencil
          ! end if
          
          i = lo(1)+1
          ! use 3rd-order slightly right-biased stencil
          rhs(i,n) = rhs(i,n) - dxinv(1) * &
               first_deriv_r3( cons(i-1:i+2,n)*q(i-1:i+2,qu) )
             
          i = lo(1)+2
          ! use 4th-order stencil
!expand          rhs(i,n) = rhs(i,n) - dxinv(1) * &
!expand               first_deriv_4( cons(i-2:i+2,n)*q(i-2:i+2,qu) )
          rhs(i,n) = rhs(i,n) - dxinv(1) * &
               ( D4(1)*(cons(i+1,n)*q(i+1,qu)-cons(i-1,n)*q(i-1,qu)) &
               + D4(2)*(cons(i+2,n)*q(i+2,qu)-cons(i-2,n)*q(i-2,qu)) )
             
          i = lo(1)+3
          ! use 6th-order stencil
!expand          rhs(i,n) = rhs(i,n) - dxinv(1) * &
!expand               first_deriv_6( cons(i-3:i+3,n)*q(i-3:i+3,qu) )
          rhs(i,n) = rhs(i,n) - dxinv(1) * &
               ( D6(1)*(cons(i+1,n)*q(i+1,qu)-cons(i-1,n)*q(i-1,qu)) &
               + D6(2)*(cons(i+2,n)*q(i+2,qu)-cons(i-2,n)*q(i-2,qu)) &
               + D6(3)*(cons(i+3,n)*q(i+3,qu)-cons(i-3,n)*q(i-3,qu)) )
       end do
    end if
    
    !
    ! ----- hi-x boundary -----
    !
    if (physbchi(1)) then 
          
       i = hi(1)-3
       ! use 6th-order stencil
       
       un(-3:3) = q(i-3:i+3,qu)
       
!expand       rhs(i,irho) = rhs(i,irho) - dxinv(1) * &
!expand            first_deriv_6( cons(i-3:i+3,imx) ) 
       rhs(i,irho) = rhs(i,irho) - dxinv(1) * &
            ( D6(1)*(cons(i+1,imx)-cons(i-1,imx)) &
            + D6(2)*(cons(i+2,imx)-cons(i-2,imx)) &
            + D6(3)*(cons(i+3,imx)-cons(i-3,imx)) )
       
!expand       rhs(i,imx) = rhs(i,imx) - dxinv(1) * &
!expand            first_deriv_6( cons(i-3:i+3,imx)*un(-3:3)+q(i-3:i+3,qpres) )
       rhs(i,imx) = rhs(i,imx) - dxinv(1) * &
            ( D6(1)*((cons(i+1,imx)*un(1)+q(i+1,qpres))-(cons(i-1,imx)*un(-1)+q(i-1,qpres))) &
            + D6(2)*((cons(i+2,imx)*un(2)+q(i+2,qpres))-(cons(i-2,imx)*un(-2)+q(i-2,qpres))) &
            + D6(3)*((cons(i+3,imx)*un(3)+q(i+3,qpres))-(cons(i-3,imx)*un(-3)+q(i-3,qpres))) )

!expand       rhs(i,iene) = rhs(i,iene) - dxinv(1) * &
!expand            first_deriv_6( (cons(i-3:i+3,iene)+q(i-3:i+3,qpres))*un(-3:3) )
       rhs(i,iene) = rhs(i,iene) - dxinv(1) * &
            ( D6(1)*((cons(i+1,iene)+q(i+1,qpres))*un(1)-(cons(i-1,iene)+q(i-1,qpres))*un(-1)) &
            + D6(2)*((cons(i+2,iene)+q(i+2,qpres))*un(2)-(cons(i-2,iene)+q(i-2,qpres))*un(-2)) &
            + D6(3)*((cons(i+3,iene)+q(i+3,qpres))*un(3)-(cons(i-3,iene)+q(i-3,qpres))*un(-3)) )
       
       i = hi(1)-2
       ! use 4th-order stencil
       
       un(-2:2) = q(i-2:i+2,qu)
       
!expand       rhs(i,irho) = rhs(i,irho) - dxinv(1) * &
!expand            first_deriv_4( cons(i-2:i+2,imx) ) 
       rhs(i,irho) = rhs(i,irho) - dxinv(1) * &
            ( D4(1)*(cons(i+1,imx)-cons(i-1,imx)) &
            + D4(2)*(cons(i+2,imx)-cons(i-2,imx)) )
       
!expand       rhs(i,imx) = rhs(i,imx) - dxinv(1) * &
!expand            first_deriv_4( cons(i-2:i+2,imx)*un(-2:2)+q(i-2:i+2,qpres) )
       rhs(i,imx) = rhs(i,imx) - dxinv(1) * &
            ( D4(1)*((cons(i+1,imx)*un(1)+q(i+1,qpres))-(cons(i-1,imx)*un(-1)+q(i-1,qpres))) &
            + D4(2)*((cons(i+2,imx)*un(2)+q(i+2,qpres))-(cons(i-2,imx)*un(-2)+q(i-2,qpres))) )
          
!expand       rhs(i,iene) = rhs(i,iene) - dxinv(1) * &
!expand            first_deriv_4( (cons(i-2:i+2,iene)+q(i-2:i+2,qpres))*un(-2:2) )
       rhs(i,iene) = rhs(i,iene) - dxinv(1) * &
            ( D4(1)*((cons(i+1,iene)+q(i+1,qpres))*un(1)-(cons(i-1,iene)+q(i-1,qpres))*un(-1)) &
            + D4(2)*((cons(i+2,iene)+q(i+2,qpres))*un(2)-(cons(i-2,iene)+q(i-2,qpres))*un(-2)) )
          
       i = hi(1)-1
       ! use 3rd-order slightly left-biased stencil
          
       un(-2:1) = q(i-2:i+1,qu)
          
       rhs(i,irho) = rhs(i,irho) - dxinv(1) * &
            first_deriv_l3( cons(i-2:i+1,imx) ) 
          
       rhs(i,imx) = rhs(i,imx) - dxinv(1) * &
            first_deriv_l3( cons(i-2:i+1,imx)*un(-2:1)+q(i-2:i+1,qpres) )
          
       rhs(i,iene) = rhs(i,iene) - dxinv(1) * &
            first_deriv_l3( (cons(i-2:i+1,iene)+q(i-2:i+1,qpres))*un(-2:1) )
          
       ! if (bchi(1) .eq. WALL???) then
       !    i = hi(1)
       !    ! use completely left-biased stencil
       ! end if
       
       do n = iry1, iry1+nspecies-1
             
          i = hi(1)-3
          ! use 6th-order stencil
!expand          rhs(i,n) = rhs(i,n) - dxinv(1) * &
!expand               first_deriv_6( cons(i-3:i+3,n)*q(i-3:i+3,qu) )
          rhs(i,n) = rhs(i,n) - dxinv(1) * &
               ( D6(1)*(cons(i+1,n)*q(i+1,qu)-cons(i-1,n)*q(i-1,qu)) &
               + D6(2)*(cons(i+2,n)*q(i+2,qu)-cons(i-2,n)*q(i-2,qu)) &
               + D6(3)*(cons(i+3,n)*q(i+3,qu)-cons(i-3,n)*q(i-3,qu)) )
          
          i = hi(1)-2
          ! use 4th-order stencil
!expand          rhs(i,n) = rhs(i,n) - dxinv(1) * &
!expand               first_deriv_4( cons(i-2:i+2,n)*q(i-2:i+2,qu) )
          rhs(i,n) = rhs(i,n) - dxinv(1) * &
               ( D4(1)*(cons(i+1,n)*q(i+1,qu)-cons(i-1,n)*q(i-1,qu)) &
               + D4(2)*(cons(i+2,n)*q(i+2,qu)-cons(i-2,n)*q(i-2,qu)) )
          
          i = hi(1)-1
          ! use 3rd-order slightly left-biased stencil
          rhs(i,n) = rhs(i,n) - dxinv(1) * &
               first_deriv_l3( cons(i-2:i+1,n)*q(i-2:i+1,qu) )
          
          ! if (bchi(1) .eq. WALL???) then
          !    i = hi(1)
          !    ! use completely left-biased stencil
          ! end if

       end do
    end if
    
    deallocate(tmpx)
    
    rhs_g(lo(1):hi(1),:) = rhs_g(lo(1):hi(1),:) + rhs(lo(1):hi(1),:)
    deallocate(rhs)

  end subroutine hypterm_1d


  subroutine narrow_diffterm_1d (lo,hi,dx,q,qlo,qhi,rhs_g,glo,ghi,rhs,rlo,rhi, &
       mu,xi,lam,dxy,dlo_g,dhi_g,bclo,bchi)

    integer,         intent(in):: dlo_g(1),dhi_g(1),bclo(1),bchi(1)
    integer,         intent(in):: lo(1),hi(1),qlo(1),qhi(1),rlo(1),rhi(1),glo(1),ghi(1)
    double precision,intent(in):: dx(1)
    double precision,intent(in)   ::  q  (qlo(1):qhi(1),nprim)
    double precision,intent(in)   ::  mu (qlo(1):qhi(1))
    double precision,intent(in)   ::  xi (qlo(1):qhi(1))
    double precision,intent(in)   ::  lam(qlo(1):qhi(1))
    double precision,intent(in)   ::  dxy(qlo(1):qhi(1),nspecies)
    double precision,intent(inout)::rhs_g(glo(1):ghi(1),ncons)
    double precision,intent(inout)::rhs  (rlo(1):rhi(1),ncons)

    integer :: slo(1), shi(1), dlo(1), dhi(1)
    double precision :: dxinv(1), dx2inv(1)

    ! used to turn off some terms
    double precision :: finlo(1), finhi(1)
    double precision :: foulo(1), fouhi(1)

    logical :: physbclo(1), physbchi(1)

    ! Only the region bounded by [dlo_g,dhi_g] contains good data.
    ! [slo,shi] will be safe for 8th-order stencil
    dlo(1) = max(lo(1)-stencil_ng, dlo_g(1))
    dhi(1) = min(hi(1)+stencil_ng, dhi_g(1))
    slo(1) = dlo(1) + stencil_ng
    shi(1) = dhi(1) - stencil_ng
    
    if (dlo(1) .eq. lo(1)) then
       physbclo(1) = .true.
    else
       physbclo(1) = .false.
    end if
    
    if (dhi(1) .eq. hi(1)) then
       physbchi(1) = .true.
    else
       physbchi(1) = .false.
    end if

    finlo = 1.d0 
    finhi = 1.d0
    foulo = 1.d0 
    fouhi = 1.d0

    if (physbclo(1)) then 
       if (bclo(1) .eq. INLET) then
          finlo(1) = 0.d0
       else if (bclo(1) .eq. OUTLET) then
          foulo(1) = 0.d0
       end if
    end if
    
    if (physbchi(1)) then 
       if (bchi(1) .eq. INLET) then
          finhi(1) = 0.d0
       else if (bchi(1) .eq. OUTLET) then
          fouhi(1) = 0.d0
       end if
    end if

    dxinv(1) = 1.0d0 / dx(1)
    dx2inv(1) = dxinv(1)**2

    rhs(lo(1):hi(1),:) = 0.d0

    call diffterm_2(q,qlo,qhi,rhs,rlo,rhi, mu,xi,lam,dxy, &
         lo,hi,slo,shi,dlo,dhi,finlo,finhi,foulo,fouhi,physbclo,physbchi,dxinv,dx2inv)

    rhs_g(lo(1):hi(1),:) = rhs_g(lo(1):hi(1),:) + rhs(lo(1):hi(1),:)

  end subroutine narrow_diffterm_1d

    
  subroutine diffterm_2(q,qlo,qhi,rhs,rlo,rhi,mu,xi,lam,dxy0, &
       lo,hi,slo,shi,dlo,dhi,finlo,finhi,foulo,fouhi,physbclo,physbchi,dxinv,dx2inv)
    use probin_module, only : reset_inactive_species, diff_gradY
    integer,         intent(in):: lo(1),hi(1),slo(1),shi(1),dlo(1),dhi(1)
    integer,         intent(in):: qlo(1),qhi(1),rlo(1),rhi(1)
    logical,         intent(in):: physbclo(1),physbchi(1)
    double precision,intent(in):: finlo(1),finhi(1),foulo(1),fouhi(1)
    double precision,intent(in):: dxinv(1),dx2inv(1)
    double precision,intent(in)   :: q (qlo(1):qhi(1),nprim)
    double precision,intent(in)   :: mu(qlo(1):qhi(1))
    double precision,intent(in)   :: xi(qlo(1):qhi(1))
    double precision,intent(in)   ::lam(qlo(1):qhi(1))
    double precision,target,intent(in)::dxy0(qlo(1):qhi(1),nspecies)
    double precision,intent(inout)::rhs(rlo(1):rhi(1),ncons)

    double precision, allocatable, dimension(:) :: vsp, dpe
    double precision, allocatable, dimension(:,:) :: Hg, dpy, dxe
    double precision, pointer :: dxy(:,:)
    ! dxy: diffusion coefficient of X in equation for Y
    ! dpy: diffusion coefficient of p in equation for Y
    ! dxe: diffusion coefficient of X in equation for energy
    ! dpe: diffusion coefficient of p in equation for energy

    integer          :: i,n, qxn, qyn, qhn, iryn, qxy1, qhias

    double precision :: mmtmp8(8,lo(1):hi(1)+1)
    double precision, allocatable, dimension(:,:) :: M8p
    double precision, allocatable, dimension(:) :: sumdrY, sumrYv, gradp
    double precision :: ry_c, ene_c

    integer :: iwrk
    double precision :: hhalf, sumdrytmp, sumryvtmp, gradptmp, Wbar, rwrk
    double precision :: Htot, Htmp(nspecies), Ytmp(nspecies), Xtmp(nspecies)
    double precision :: M6p(6), M6X(6), mmtmp6(6)
    double precision :: M4p(4), M4X(4), mmtmp4(4)
    double precision :: M2p(2), M2X(2), mmtmp2(2)
    double precision :: BBp(4), BBX(4), mmtmpB(4)
    double precision :: rhstmp(nspecies), rhstot, rhsene
    double precision :: Hcell(0:1,2:ncons)
    integer :: iface

    if (diff_gradY) then  ! original LMC formulation
       qxy1 = qy1
       allocate(dxy(dlo(1):dhi(1),nspecies))
       dxy = dxy0(dlo(1):dhi(1),:)
    else                  ! Original SMC formulation (i.e., LMC w/ Wbar)
       qxy1 = qx1
       dxy => dxy0
    end if

    allocate(vsp(dlo(1):dhi(1)))

    allocate(dpy(dlo(1):dhi(1),nspecies))
    allocate(dxe(dlo(1):dhi(1),nspecies))
    allocate(dpe(dlo(1):dhi(1)))

    allocate(Hg(lo(1):hi(1)+1,2:ncons))

    allocate(M8p(8,lo(1):hi(1)+1))

    allocate(sumdrY(lo(1):hi(1)))
    allocate(sumrYv(lo(1):hi(1)))
    allocate(gradp (lo(1):hi(1)))

    do i=dlo(1),dhi(1)
       vsp(i) = xi(i) + FourThirds*mu(i)
    enddo

    dpe = 0.d0

    if (reset_inactive_species) then
       
       qhias = qh1 + iias - 1

       do n=1,nspecies
          if (n .eq. iias) then  ! inactive speices
             do i=dlo(1),dhi(1)
                dpy(i,n) = 0.d0
                dxe(i,n) = 0.d0
             end do
          else
             qxn = qx1+n-1
             qyn = qy1+n-1
             qhn = qh1+n-1
             do i=dlo(1),dhi(1)
                dpy(i,n) = dxy(i,n)/q(i,qpres)*(q(i,qxn)-q(i,qyn))
                dxe(i,n) = dxy(i,n)*(q(i,qhn)-q(i,qhias))
                dpe(i) = dpe(i) + dpy(i,n)*(q(i,qhn)-q(i,qhias))
             end do
          end if
       end do

    else

       do n=1,nspecies
          qxn = qx1+n-1
          qyn = qy1+n-1
          qhn = qh1+n-1
          do i=dlo(1),dhi(1)
             dpy(i,n) = dxy(i,n)/q(i,qpres)*(q(i,qxn)-q(i,qyn))
             dxe(i,n) = dxy(i,n)*q(i,qhn)
             dpe(i) = dpe(i) + dpy(i,n)*q(i,qhn)
          end do
       end do

    end if

    if (diff_gradY) then
       do i=dlo(1),dhi(1)
          Xtmp = q(i,qx1:qx1+nspecies-1)
          call ckmmwx(Xtmp, iwrk, rwrk, Wbar)
          do n=1,nspecies
             dxy(i,n) = dxy(i,n) * Wbar * inv_mwt(n)
             dxe(i,n) = dxe(i,n) * Wbar * inv_mwt(n)
          end do
       end do
    end if
    
    ! ------- BEGIN x-direction -------

    do i=slo(1),shi(1)+1
!expand       mmtmp8(1:8,i) = matmul(vsp(i-4:i+3), M8)
       mmtmp8(1,i) = vsp(i-4) * M8(1,1) &
                   + vsp(i-3) * M8(2,1) &
                   + vsp(i-2) * M8(3,1) &
                   + vsp(i-1) * M8(4,1) &
                   + vsp(i  ) * M8(5,1)
       mmtmp8(2,i) = vsp(i-4) * M8(1,2) &
                   + vsp(i-3) * M8(2,2) &
                   + vsp(i-2) * M8(3,2) &
                   + vsp(i-1) * M8(4,2) &
                   + vsp(i  ) * M8(5,2) &
                   + vsp(i+1) * M8(6,2)
       mmtmp8(3,i) = vsp(i-4) * M8(1,3) &
                   + vsp(i-3) * M8(2,3) &
                   + vsp(i-2) * M8(3,3) &
                   + vsp(i-1) * M8(4,3) &
                   + vsp(i  ) * M8(5,3) &
                   + vsp(i+1) * M8(6,3) &
                   + vsp(i+2) * M8(7,3)
       mmtmp8(4,i) = vsp(i-4) * M8(1,4) &
                   + vsp(i-3) * M8(2,4) &
                   + vsp(i-2) * M8(3,4) &
                   + vsp(i-1) * M8(4,4) &
                   + vsp(i  ) * M8(5,4) &
                   + vsp(i+1) * M8(6,4) &
                   + vsp(i+2) * M8(7,4) &
                   + vsp(i+3) * M8(8,4)
       mmtmp8(5,i) = vsp(i-4) * M8(1,5) &
                   + vsp(i-3) * M8(2,5) &
                   + vsp(i-2) * M8(3,5) &
                   + vsp(i-1) * M8(4,5) &
                   + vsp(i  ) * M8(5,5) &
                   + vsp(i+1) * M8(6,5) &
                   + vsp(i+2) * M8(7,5) &
                   + vsp(i+3) * M8(8,5)
       mmtmp8(6,i) = vsp(i-3) * M8(2,6) &
                   + vsp(i-2) * M8(3,6) &
                   + vsp(i-1) * M8(4,6) &
                   + vsp(i  ) * M8(5,6) &
                   + vsp(i+1) * M8(6,6) &
                   + vsp(i+2) * M8(7,6) &
                   + vsp(i+3) * M8(8,6)
       mmtmp8(7,i) = vsp(i-2) * M8(3,7) &
                   + vsp(i-1) * M8(4,7) &
                   + vsp(i  ) * M8(5,7) &
                   + vsp(i+1) * M8(6,7) &
                   + vsp(i+2) * M8(7,7) &
                   + vsp(i+3) * M8(8,7)
       mmtmp8(8,i) = vsp(i-1) * M8(4,8) &
                   + vsp(i  ) * M8(5,8) &
                   + vsp(i+1) * M8(6,8) &
                   + vsp(i+2) * M8(7,8) &
                   + vsp(i+3) * M8(8,8)
!expand       Hg(i,imx) = dot_product(mmtmp8(1:8,i), q(i-4:i+3,qu))
       Hg(i,imx) =  &
            ( mmtmp8(1,i)*q(i-4,qu) + mmtmp8(2,i)*q(i-3,qu) &
            + mmtmp8(3,i)*q(i-2,qu) + mmtmp8(4,i)*q(i-1,qu) &
            + mmtmp8(5,i)*q(i  ,qu) + mmtmp8(6,i)*q(i+1,qu) &
            + mmtmp8(7,i)*q(i+2,qu) + mmtmp8(8,i)*q(i+3,qu) )
    end do

    do i=slo(1),shi(1)+1
!expand       mmtmp8(1:8,i) = matmul(lam(i-4:i+3), M8)
       mmtmp8(1,i) = lam(i-4) * M8(1,1) &
                   + lam(i-3) * M8(2,1) &
                   + lam(i-2) * M8(3,1) &
                   + lam(i-1) * M8(4,1) &
                   + lam(i  ) * M8(5,1)
       mmtmp8(2,i) = lam(i-4) * M8(1,2) &
                   + lam(i-3) * M8(2,2) &
                   + lam(i-2) * M8(3,2) &
                   + lam(i-1) * M8(4,2) &
                   + lam(i  ) * M8(5,2) &
                   + lam(i+1) * M8(6,2)
       mmtmp8(3,i) = lam(i-4) * M8(1,3) &
                   + lam(i-3) * M8(2,3) &
                   + lam(i-2) * M8(3,3) &
                   + lam(i-1) * M8(4,3) &
                   + lam(i  ) * M8(5,3) &
                   + lam(i+1) * M8(6,3) &
                   + lam(i+2) * M8(7,3)
       mmtmp8(4,i) = lam(i-4) * M8(1,4) &
                   + lam(i-3) * M8(2,4) &
                   + lam(i-2) * M8(3,4) &
                   + lam(i-1) * M8(4,4) &
                   + lam(i  ) * M8(5,4) &
                   + lam(i+1) * M8(6,4) &
                   + lam(i+2) * M8(7,4) &
                   + lam(i+3) * M8(8,4)
       mmtmp8(5,i) = lam(i-4) * M8(1,5) &
                   + lam(i-3) * M8(2,5) &
                   + lam(i-2) * M8(3,5) &
                   + lam(i-1) * M8(4,5) &
                   + lam(i  ) * M8(5,5) &
                   + lam(i+1) * M8(6,5) &
                   + lam(i+2) * M8(7,5) &
                   + lam(i+3) * M8(8,5)
       mmtmp8(6,i) = lam(i-3) * M8(2,6) &
                   + lam(i-2) * M8(3,6) &
                   + lam(i-1) * M8(4,6) &
                   + lam(i  ) * M8(5,6) &
                   + lam(i+1) * M8(6,6) &
                   + lam(i+2) * M8(7,6) &
                   + lam(i+3) * M8(8,6)
       mmtmp8(7,i) = lam(i-2) * M8(3,7) &
                   + lam(i-1) * M8(4,7) &
                   + lam(i  ) * M8(5,7) &
                   + lam(i+1) * M8(6,7) &
                   + lam(i+2) * M8(7,7) &
                   + lam(i+3) * M8(8,7)
       mmtmp8(8,i) = lam(i-1) * M8(4,8) &
                   + lam(i  ) * M8(5,8) &
                   + lam(i+1) * M8(6,8) &
                   + lam(i+2) * M8(7,8) &
                   + lam(i+3) * M8(8,8)
!expand       Hg(i,iene) = dot_product(mmtmp8(1:8,i), q(i-4:i+3,qtemp))
       Hg(i,iene) =  &
            ( mmtmp8(1,i)*q(i-4,qtemp) + mmtmp8(2,i)*q(i-3,qtemp) &
            + mmtmp8(3,i)*q(i-2,qtemp) + mmtmp8(4,i)*q(i-1,qtemp) &
            + mmtmp8(5,i)*q(i  ,qtemp) + mmtmp8(6,i)*q(i+1,qtemp) &
            + mmtmp8(7,i)*q(i+2,qtemp) + mmtmp8(8,i)*q(i+3,qtemp) )
    end do

    do i=slo(1),shi(1)+1
!expand       mmtmp8(1:8,i) = matmul(M8, q(i-4:i+3,qpres))
       mmtmp8(1,i) = M8T(1,1) * q(i-4,qpres) &
                   + M8T(2,1) * q(i-3,qpres) &
                   + M8T(3,1) * q(i-2,qpres) &
                   + M8T(4,1) * q(i-1,qpres) &
                   + M8T(5,1) * q(i  ,qpres)
       mmtmp8(2,i) = M8T(1,2) * q(i-4,qpres) &
                   + M8T(2,2) * q(i-3,qpres) &
                   + M8T(3,2) * q(i-2,qpres) &
                   + M8T(4,2) * q(i-1,qpres) &
                   + M8T(5,2) * q(i  ,qpres) &
                   + M8T(6,2) * q(i+1,qpres)
       mmtmp8(3,i) = M8T(1,3) * q(i-4,qpres) &
                   + M8T(2,3) * q(i-3,qpres) &
                   + M8T(3,3) * q(i-2,qpres) &
                   + M8T(4,3) * q(i-1,qpres) &
                   + M8T(5,3) * q(i  ,qpres) &
                   + M8T(6,3) * q(i+1,qpres) &
                   + M8T(7,3) * q(i+2,qpres)
       mmtmp8(4,i) = M8T(1,4) * q(i-4,qpres) &
                   + M8T(2,4) * q(i-3,qpres) &
                   + M8T(3,4) * q(i-2,qpres) &
                   + M8T(4,4) * q(i-1,qpres) &
                   + M8T(5,4) * q(i  ,qpres) &
                   + M8T(6,4) * q(i+1,qpres) &
                   + M8T(7,4) * q(i+2,qpres) &
                   + M8T(8,4) * q(i+3,qpres)
       mmtmp8(5,i) = M8T(1,5) * q(i-4,qpres) &
                   + M8T(2,5) * q(i-3,qpres) &
                   + M8T(3,5) * q(i-2,qpres) &
                   + M8T(4,5) * q(i-1,qpres) &
                   + M8T(5,5) * q(i  ,qpres) &
                   + M8T(6,5) * q(i+1,qpres) &
                   + M8T(7,5) * q(i+2,qpres) &
                   + M8T(8,5) * q(i+3,qpres)
       mmtmp8(6,i) = M8T(2,6) * q(i-3,qpres) &
                   + M8T(3,6) * q(i-2,qpres) &
                   + M8T(4,6) * q(i-1,qpres) &
                   + M8T(5,6) * q(i  ,qpres) &
                   + M8T(6,6) * q(i+1,qpres) &
                   + M8T(7,6) * q(i+2,qpres) &
                   + M8T(8,6) * q(i+3,qpres)
       mmtmp8(7,i) = M8T(3,7) * q(i-2,qpres) &
                   + M8T(4,7) * q(i-1,qpres) &
                   + M8T(5,7) * q(i  ,qpres) &
                   + M8T(6,7) * q(i+1,qpres) &
                   + M8T(7,7) * q(i+2,qpres) &
                   + M8T(8,7) * q(i+3,qpres)
       mmtmp8(8,i) = M8T(4,8) * q(i-1,qpres) &
                   + M8T(5,8) * q(i  ,qpres) &
                   + M8T(6,8) * q(i+1,qpres) &
                   + M8T(7,8) * q(i+2,qpres) &
                   + M8T(8,8) * q(i+3,qpres)
!expand       Hg(i,iene) = Hg(i,iene) + dot_product(dpe(i-4:i+3), mmtmp8(1:8,i))
       Hg(i,iene) = Hg(i,iene)+ &
            ( dpe(i-4)*mmtmp8(1,i) + dpe(i-3)*mmtmp8(2,i) &
            + dpe(i-2)*mmtmp8(3,i) + dpe(i-1)*mmtmp8(4,i) &
            + dpe(i  )*mmtmp8(5,i) + dpe(i+1)*mmtmp8(6,i) &
            + dpe(i+2)*mmtmp8(7,i) + dpe(i+3)*mmtmp8(8,i) )
    end do
    do i=slo(1),shi(1)+1
       M8p(:,i) = mmtmp8(1:8,i)
    end do

    do n = 1, nspecies

       if (n .eq. iias) cycle  ! inactive speices

       qxn = qxy1+n-1    ! qxn might point to Y!
       iryn = iry1+n-1

       do i=slo(1),shi(1)+1
!expand          Hg(i,iryn) = dot_product(dpy(i-4:i+3,n), M8p(:,i))
          Hg(i,iryn) =  &
               ( dpy(i-4,n)*M8p(1,i) + dpy(i-3,n)*M8p(2,i) &
               + dpy(i-2,n)*M8p(3,i) + dpy(i-1,n)*M8p(4,i) &
               + dpy(i  ,n)*M8p(5,i) + dpy(i+1,n)*M8p(6,i) &
               + dpy(i+2,n)*M8p(7,i) + dpy(i+3,n)*M8p(8,i) )
       end do

       do i=slo(1),shi(1)+1
!expand          mmtmp8(1:8,i) = matmul(M8, q(i-4:i+3,qxn))
          mmtmp8(1,i) = M8T(1,1) * q(i-4,qxn) &
                      + M8T(2,1) * q(i-3,qxn) &
                      + M8T(3,1) * q(i-2,qxn) &
                      + M8T(4,1) * q(i-1,qxn) &
                      + M8T(5,1) * q(i  ,qxn)
          mmtmp8(2,i) = M8T(1,2) * q(i-4,qxn) &
                      + M8T(2,2) * q(i-3,qxn) &
                      + M8T(3,2) * q(i-2,qxn) &
                      + M8T(4,2) * q(i-1,qxn) &
                      + M8T(5,2) * q(i  ,qxn) &
                      + M8T(6,2) * q(i+1,qxn)
          mmtmp8(3,i) = M8T(1,3) * q(i-4,qxn) &
                      + M8T(2,3) * q(i-3,qxn) &
                      + M8T(3,3) * q(i-2,qxn) &
                      + M8T(4,3) * q(i-1,qxn) &
                      + M8T(5,3) * q(i  ,qxn) &
                      + M8T(6,3) * q(i+1,qxn) &
                      + M8T(7,3) * q(i+2,qxn)
          mmtmp8(4,i) = M8T(1,4) * q(i-4,qxn) &
                      + M8T(2,4) * q(i-3,qxn) &
                      + M8T(3,4) * q(i-2,qxn) &
                      + M8T(4,4) * q(i-1,qxn) &
                      + M8T(5,4) * q(i  ,qxn) &
                      + M8T(6,4) * q(i+1,qxn) &
                      + M8T(7,4) * q(i+2,qxn) &
                      + M8T(8,4) * q(i+3,qxn)
          mmtmp8(5,i) = M8T(1,5) * q(i-4,qxn) &
                      + M8T(2,5) * q(i-3,qxn) &
                      + M8T(3,5) * q(i-2,qxn) &
                      + M8T(4,5) * q(i-1,qxn) &
                      + M8T(5,5) * q(i  ,qxn) &
                      + M8T(6,5) * q(i+1,qxn) &
                      + M8T(7,5) * q(i+2,qxn) &
                      + M8T(8,5) * q(i+3,qxn)
          mmtmp8(6,i) = M8T(2,6) * q(i-3,qxn) &
                      + M8T(3,6) * q(i-2,qxn) &
                      + M8T(4,6) * q(i-1,qxn) &
                      + M8T(5,6) * q(i  ,qxn) &
                      + M8T(6,6) * q(i+1,qxn) &
                      + M8T(7,6) * q(i+2,qxn) &
                      + M8T(8,6) * q(i+3,qxn)
          mmtmp8(7,i) = M8T(3,7) * q(i-2,qxn) &
                      + M8T(4,7) * q(i-1,qxn) &
                      + M8T(5,7) * q(i  ,qxn) &
                      + M8T(6,7) * q(i+1,qxn) &
                      + M8T(7,7) * q(i+2,qxn) &
                      + M8T(8,7) * q(i+3,qxn)
          mmtmp8(8,i) = M8T(4,8) * q(i-1,qxn) &
                      + M8T(5,8) * q(i  ,qxn) &
                      + M8T(6,8) * q(i+1,qxn) &
                      + M8T(7,8) * q(i+2,qxn) &
                      + M8T(8,8) * q(i+3,qxn)
!expand          Hg(i,iene) = Hg(i,iene) + dot_product(dxe(i-4:i+3,n), mmtmp8(1:8,i))
          Hg(i,iene) = Hg(i,iene)+ &
               ( dxe(i-4,n)*mmtmp8(1,i) + dxe(i-3,n)*mmtmp8(2,i) &
               + dxe(i-2,n)*mmtmp8(3,i) + dxe(i-1,n)*mmtmp8(4,i) &
               + dxe(i  ,n)*mmtmp8(5,i) + dxe(i+1,n)*mmtmp8(6,i) &
               + dxe(i+2,n)*mmtmp8(7,i) + dxe(i+3,n)*mmtmp8(8,i) )
!expand          Hg(i,iryn) = Hg(i,iryn) &
!expand               + dot_product(dxy(i-4:i+3,n), mmtmp8(1:8,i))
          Hg(i,iryn) = Hg(i,iryn)+ &
               ( dxy(i-4,n)*mmtmp8(1,i) + dxy(i-3,n)*mmtmp8(2,i) &
               + dxy(i-2,n)*mmtmp8(3,i) + dxy(i-1,n)*mmtmp8(4,i) &
               + dxy(i  ,n)*mmtmp8(5,i) + dxy(i+1,n)*mmtmp8(6,i) &
               + dxy(i+2,n)*mmtmp8(7,i) + dxy(i+3,n)*mmtmp8(8,i) )
       end do

    end do

    ! add x-direction rhs

    do n=2,iene
       do i=slo(1),shi(1)
          rhs(i,n) = rhs(i,n) + (Hg(i+1,n) - Hg(i,n)) * dx2inv(1)
       end do
    end do
       
    sumdrY = 0.d0
    do n=iry1,ncons

       if (n.eq.iry_ias) cycle

       do i=slo(1),shi(1)
          sumdry(i) = sumdry(i) + (Hg(i+1,n) - Hg(i,n)) * dx2inv(1)
          rhs(i,n)  =  rhs(i,n) + (Hg(i+1,n) - Hg(i,n)) * dx2inv(1)
       end do

    end do

    if (reset_inactive_species) then

       do i=slo(1),shi(1)
          rhs(i,iry_ias) = rhs(i,iry_ias) - sumdry(i)
       end do

    else

       do i=slo(1),shi(1)
!expand          gradp(i) = dxinv(1) * first_deriv_8(q(i-4:i+4,qpres))
          gradp(i) = dxinv(1) * &
               ( D8(1)*(q(i+1,qpres)-q(i-1,qpres)) &
               + D8(2)*(q(i+2,qpres)-q(i-2,qpres)) &
               + D8(3)*(q(i+3,qpres)-q(i-3,qpres)) &
               + D8(4)*(q(i+4,qpres)-q(i-4,qpres)) )
       end do
       
       sumryv = 0.d0
       do n = 1, nspecies
          qxn = qxy1+n-1    ! qxn might point to Y!
          do i=slo(1),shi(1)
!expand             sumryv(i) = sumryv(i) + dpy(i,n)*gradp(i)  &
!expand                  + dxy(i,n)*dxinv(1)*first_deriv_8(q(i-4:i+4,qxn))
             sumryv(i) = sumryv(i)+dpy(i,n)*gradp(i)+dxy(i,n) * dxinv(1) * &
                  ( D8(1)*(q(i+1,qxn)-q(i-1,qxn)) &
                  + D8(2)*(q(i+2,qxn)-q(i-2,qxn)) &
                  + D8(3)*(q(i+3,qxn)-q(i-3,qxn)) &
                  + D8(4)*(q(i+4,qxn)-q(i-4,qxn)) )
          end do
       end do

       do n=1,nspecies
          qyn = qy1+n-1
          qhn = qh1+n-1
          iryn = iry1+n-1
          
          do i=slo(1),shi(1)
!expand             ry_c = q(i,qyn)*sumdry(i) + sumryv(i)*dxinv(1) * &
!expand                  first_deriv_8(q(i-4:i+4,qyn))
             ry_c = q(i,qyn)*sumdry(i)+sumryv(i) * dxinv(1) * &
                  ( D8(1)*(q(i+1,qyn)-q(i-1,qyn)) &
                  + D8(2)*(q(i+2,qyn)-q(i-2,qyn)) &
                  + D8(3)*(q(i+3,qyn)-q(i-3,qyn)) &
                  + D8(4)*(q(i+4,qyn)-q(i-4,qyn)) )
!expand             ene_c = ry_c*q(i,qhn) + q(i,qyn)*sumryv(i)*dxinv(1)* &
!expand                  first_deriv_8(q(i-4:i+4,qhn))
             ene_c = ry_c*q(i,qhn)+q(i,qyn)*sumryv(i) * dxinv(1) * &
                  ( D8(1)*(q(i+1,qhn)-q(i-1,qhn)) &
                  + D8(2)*(q(i+2,qhn)-q(i-2,qhn)) &
                  + D8(3)*(q(i+3,qhn)-q(i-3,qhn)) &
                  + D8(4)*(q(i+4,qhn)-q(i-4,qhn)) )
             rhs(i,iene) = rhs(i,iene) - ene_c
             rhs(i,iryn) = rhs(i,iryn) - ry_c
          end do
       end do

    end if

    ! ------- END x-direction -------

    !
    ! lo-x boundary
    !
    if (physbclo(1)) then
       i = lo(1)
       ! use completely right-biased stencil
       mmtmpB = matmul(vsp(i:i+3), BRB)
!expand       rhs(i,imx) = rhs(i,imx)+finlo(1)*dx2inv(1)*dot_product(mmtmpB,q(i:i+3,qu))
       rhs(i,imx) = rhs(i,imx)+finlo(1)*dx2inv(1)* &
            ( q(i  ,qu)*mmtmpB(1) + q(i+1,qu)*mmtmpB(2) &
            + q(i+2,qu)*mmtmpB(3) + q(i+3,qu)*mmtmpB(4) )
          
       mmtmpB = matmul(lam(i:i+3), BRB)
       BBp = matmul(BRB, q(i:i+3,qpres))
!expand       rhs(i,iene) = rhs(i,iene) + foulo(1)*dx2inv(1) * &
!expand            ( dot_product(mmtmpB, q(i:i+3,qtemp)) &
!expand            + dot_product(      dpe(i:i+3), BBp) )
       rhs(i,iene) = rhs(i,iene)+foulo(1)*dx2inv(1)* &
            ( ( q(i  ,qtemp)*mmtmpB(1) + q(i+1,qtemp)*mmtmpB(2) &
              + q(i+2,qtemp)*mmtmpB(3) + q(i+3,qtemp)*mmtmpB(4) ) &
            + ( dpe(i  )*BBp(1) + dpe(i+1)*BBp(2) &
              + dpe(i+2)*BBp(3) + dpe(i+3)*BBp(4) ) )
          
       rhstot = 0.d0
       rhsene = 0.d0
       do n = 1, nspecies
          
          if (n .eq. iias) cycle  ! inactive speices

          qxn = qxy1+n-1
          qyn = qy1+n-1
          
          BBX = matmul(BRB, q(i:i+3,qxn))
          
!expand          rhstmp(n) = dot_product(dpy(i:i+3,n), BBp) &
!expand               +      dot_product(dxy(i:i+3,n), BBX)
          rhstmp(n) =  &
               ( ( dpy(i  ,n)*BBp(1) + dpy(i+1,n)*BBp(2) &
                 + dpy(i+2,n)*BBp(3) + dpy(i+3,n)*BBp(4) ) &
               + ( dxy(i  ,n)*BBX(1) + dxy(i+1,n)*BBX(2) &
                 + dxy(i+2,n)*BBX(3) + dxy(i+3,n)*BBX(4) ) )
             
!expand          rhsene = rhsene &
!expand               +      dot_product(dxe(i:i+3,n), BBX)
          rhsene = rhsene+ &
               ( dxe(i  ,n)*BBX(1) + dxe(i+1,n)*BBX(2) &
               + dxe(i+2,n)*BBX(3) + dxe(i+3,n)*BBX(4) )
          
          rhstot = rhstot + rhstmp(n)
          Ytmp(n) = q(i,qyn)
          
          rhs(i,iry1+n-1) = rhs(i,iry1+n-1) + foulo(1)*dx2inv(1)*rhstmp(n)
       end do

       rhs(i,iene) = rhs(i,iene) + foulo(1)*dx2inv(1) * rhsene
       
       if (reset_inactive_species) then
          
          rhs(i,iry_ias) = rhs(i,iry_ias) - foulo(1)*dx2inv(1)*rhstot
          
       else
          
          do n = 1, nspecies
             rhs(i,iry1+n-1) =  rhs(i,iry1+n-1) - &
                  foulo(1)*dx2inv(1) * Ytmp(n) * rhstot
          end do
          
          rhsene = 0.d0
          do n = 1, nspecies
             qhn = qh1+n-1
             rhsene = rhsene - Ytmp(n) * q(i,qhn) * rhstot
          end do
          rhs(i,iene) = rhs(i,iene) + foulo(1)*dx2inv(1) * rhsene
          
       end if
          
          
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! use 2nd-order stencil for cell lo(1)+1,
       do iface=0,1 
          i = lo(1)+1 + iface
             
          mmtmp2 = matmul(vsp(i-1:i), M2)
!expand          Hcell(iface,imx) = dot_product(mmtmp2, q(i-1:i,qu))
          Hcell(iface,imx) =  &
               ( q(i-1,qu)*mmtmp2(1) + q(i  ,qu)*mmtmp2(2) )
             
          mmtmp2 = matmul(lam(i-1:i), M2)
          M2p = matmul(M2,  q(i-1:i,qpres))
!expand          Hcell(iface,iene) = dot_product(mmtmp2, q(i-1:i,qtemp)) &
!expand               &            + dot_product(      dpe(i-1:i), M2p)
          Hcell(iface,iene) =  &
               ( ( q(i-1,qtemp)*mmtmp2(1) + q(i  ,qtemp)*mmtmp2(2) ) &
               + ( dpe(i-1)*M2p(1) + dpe(i  )*M2p(2) ) )
             
          Htot = 0.d0
          do n = 1, nspecies
             
             if (n .eq. iias) cycle  ! inactive speices
             
             qxn = qxy1+n-1
             qyn = qy1+n-1
             
             M2X = matmul(M2, q(i-1:i,qxn))
             
!expand             Htmp(n) = dot_product(dpy(i-1:i,n), M2p) &
!expand                  +    dot_product(dxy(i-1:i,n), M2X)
             Htmp(n) =  &
                  ( ( dpy(i-1,n)*M2p(1) + dpy(i  ,n)*M2p(2) ) &
                  + ( dxy(i-1,n)*M2X(1) + dxy(i  ,n)*M2X(2) ) )
             
!expand             Hcell(iface,iene) = Hcell(iface,iene) &
!expand                  +    dot_product(dxe(i-1:i,n), M2X)
             Hcell(iface,iene) = Hcell(iface,iene)+ &
                  ( dxe(i-1,n)*M2X(1) + dxe(i  ,n)*M2X(2) )
             
             Htot = Htot + Htmp(n)
             Ytmp(n) = 0.5d0*(q(i-1,qyn) + q(i,qyn))
             
             Hcell(iface,iry1+n-1) = Htmp(n)
          end do
          
          if (.not. reset_inactive_species) then
             do n = 1, nspecies
                Hcell(iface,iry1+n-1) = Hcell(iface,iry1+n-1) - Ytmp(n)*Htot
             end do
             
             do n = 1, nspecies
                qhn = qh1+n-1
                hhalf = 0.5d0*(q(i-1,qhn) + q(i,qhn))
                Hcell(iface,iene) =  Hcell(iface,iene) - Ytmp(n)*Htot*hhalf 
             end do
          end if
       end do

       i = lo(1)+1
       if (reset_inactive_species) then
          do n=2,iene
             rhs(i,n) = rhs(i,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(1)
          end do
          sumdrytmp = 0.d0
          do n=iry1,ncons
             if (n.eq.iry_ias) cycle
             sumdrytmp  = sumdrytmp  + (Hcell(1,n) - Hcell(0,n)) * dx2inv(1)
             rhs(i,n) = rhs(i,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(1)
          end do
          rhs(i,iry_ias) = rhs(i,iry_ias) - sumdrytmp
       else
          do n=2,ncons
             rhs(i,n) = rhs(i,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(1)
          end do
       end if
          
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! use 4th-order stencil for cell lo(1)+2,
       do iface=0,1 
          i = lo(1)+2 + iface
          
!expand          mmtmp4 = matmul(vsp(i-2:i+1), M4)
          mmtmp4(1) = vsp(i-2) * M4(1,1) &
                    + vsp(i-1) * M4(2,1) &
                    + vsp(i  ) * M4(3,1)
          mmtmp4(2) = vsp(i-2) * M4(1,2) &
                    + vsp(i-1) * M4(2,2) &
                    + vsp(i  ) * M4(3,2) &
                    + vsp(i+1) * M4(4,2)
          mmtmp4(3) = vsp(i-2) * M4(1,3) &
                    + vsp(i-1) * M4(2,3) &
                    + vsp(i  ) * M4(3,3) &
                    + vsp(i+1) * M4(4,3)
          mmtmp4(4) = vsp(i-1) * M4(2,4) &
                    + vsp(i  ) * M4(3,4) &
                    + vsp(i+1) * M4(4,4)
!expand          Hcell(iface,imx) = dot_product(mmtmp4, q(i-2:i+1,qu))
          Hcell(iface,imx) =  &
               ( q(i-2,qu)*mmtmp4(1) + q(i-1,qu)*mmtmp4(2) &
               + q(i  ,qu)*mmtmp4(3) + q(i+1,qu)*mmtmp4(4) )
          
!expand          mmtmp4 = matmul(lam(i-2:i+1), M4)
          mmtmp4(1) = lam(i-2) * M4(1,1) &
                    + lam(i-1) * M4(2,1) &
                    + lam(i  ) * M4(3,1)
          mmtmp4(2) = lam(i-2) * M4(1,2) &
                    + lam(i-1) * M4(2,2) &
                    + lam(i  ) * M4(3,2) &
                    + lam(i+1) * M4(4,2)
          mmtmp4(3) = lam(i-2) * M4(1,3) &
                    + lam(i-1) * M4(2,3) &
                    + lam(i  ) * M4(3,3) &
                    + lam(i+1) * M4(4,3)
          mmtmp4(4) = lam(i-1) * M4(2,4) &
                    + lam(i  ) * M4(3,4) &
                    + lam(i+1) * M4(4,4)
!expand          M4p = matmul(M4,  q(i-2:i+1,qpres))
          M4p(1) = M4T(1,1) * q(i-2,qpres) &
                 + M4T(2,1) * q(i-1,qpres) &
                 + M4T(3,1) * q(i  ,qpres)
          M4p(2) = M4T(1,2) * q(i-2,qpres) &
                 + M4T(2,2) * q(i-1,qpres) &
                 + M4T(3,2) * q(i  ,qpres) &
                 + M4T(4,2) * q(i+1,qpres)
          M4p(3) = M4T(1,3) * q(i-2,qpres) &
                 + M4T(2,3) * q(i-1,qpres) &
                 + M4T(3,3) * q(i  ,qpres) &
                 + M4T(4,3) * q(i+1,qpres)
          M4p(4) = M4T(2,4) * q(i-1,qpres) &
                 + M4T(3,4) * q(i  ,qpres) &
                 + M4T(4,4) * q(i+1,qpres)
!expand          Hcell(iface,iene) = dot_product(mmtmp4, q(i-2:i+1,qtemp))      &
!expand               &            + dot_product(      dpe(i-2:i+1), M4p)
          Hcell(iface,iene) =  &
               ( ( q(i-2,qtemp)*mmtmp4(1) + q(i-1,qtemp)*mmtmp4(2) &
                 + q(i  ,qtemp)*mmtmp4(3) + q(i+1,qtemp)*mmtmp4(4) ) &
               + ( dpe(i-2)*M4p(1) + dpe(i-1)*M4p(2) &
                 + dpe(i  )*M4p(3) + dpe(i+1)*M4p(4) ) )
          
          do n = 1, nspecies
             if (n .eq. iias) cycle  ! inactive speices
             
             qxn = qxy1+n-1
             iryn = iry1+n-1
!expand             M4X = matmul(M4, q(i-2:i+1,qxn))
             M4X(1) = M4T(1,1) * q(i-2,qxn) &
                    + M4T(2,1) * q(i-1,qxn) &
                    + M4T(3,1) * q(i  ,qxn)
             M4X(2) = M4T(1,2) * q(i-2,qxn) &
                    + M4T(2,2) * q(i-1,qxn) &
                    + M4T(3,2) * q(i  ,qxn) &
                    + M4T(4,2) * q(i+1,qxn)
             M4X(3) = M4T(1,3) * q(i-2,qxn) &
                    + M4T(2,3) * q(i-1,qxn) &
                    + M4T(3,3) * q(i  ,qxn) &
                    + M4T(4,3) * q(i+1,qxn)
             M4X(4) = M4T(2,4) * q(i-1,qxn) &
                    + M4T(3,4) * q(i  ,qxn) &
                    + M4T(4,4) * q(i+1,qxn)
!expand             Hcell(iface,iene) = Hcell(iface,iene) &
!expand                  +    dot_product(dxe(i-2:i+1,n), M4X)
             Hcell(iface,iene) = Hcell(iface,iene)+ &
                  ( dxe(i-2,n)*M4X(1) + dxe(i-1,n)*M4X(2) &
                  + dxe(i  ,n)*M4X(3) + dxe(i+1,n)*M4X(4) )
!expand             Hcell(iface,iryn) = dot_product(dpy(i-2:i+1,n), M4p) &
!expand                  + dot_product(dxy(i-2:i+1,n), M4X)
             Hcell(iface,iryn) =  &
                  ( ( dpy(i-2,n)*M4p(1) + dpy(i-1,n)*M4p(2) &
                    + dpy(i  ,n)*M4p(3) + dpy(i+1,n)*M4p(4) ) &
                  + ( dxy(i-2,n)*M4X(1) + dxy(i-1,n)*M4X(2) &
                    + dxy(i  ,n)*M4X(3) + dxy(i+1,n)*M4X(4) ) )
          end do
          
       end do
          
       i = lo(1)+2
       
       do n=2,iene
          rhs(i,n) = rhs(i,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(1)
       end do
       
       sumdrytmp = 0.d0
       do n=iry1,ncons
          if (n.eq.iry_ias) cycle
          sumdrytmp  =  sumdrytmp + (Hcell(1,n) - Hcell(0,n)) * dx2inv(1)
          rhs(i,n) = rhs(i,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(1)
       end do
       
       if (reset_inactive_species) then
          
          rhs(i,iry_ias) = rhs(i,iry_ias) - sumdrytmp
          
       else
          
!expand          gradptmp = dxinv(1) * first_deriv_4(q(i-2:i+2,qpres))
          gradptmp = dxinv(1) * &
               ( D4(1)*(q(i+1,qpres)-q(i-1,qpres)) &
               + D4(2)*(q(i+2,qpres)-q(i-2,qpres)) )
          
          sumryvtmp = 0.d0
          do n = 1, nspecies
             qxn = qxy1+n-1
!expand             sumryvtmp = sumryvtmp + dpy(i,n)*gradptmp  &
!expand                  + dxy(i,n)*dxinv(1)*first_deriv_4(q(i-2:i+2,qxn))
             sumryvtmp = sumryvtmp+dpy(i,n)*gradptmp+dxy(i,n) * dxinv(1) * &
                  ( D4(1)*(q(i+1,qxn)-q(i-1,qxn)) &
                  + D4(2)*(q(i+2,qxn)-q(i-2,qxn)) )
          end do
          
          do n=1,nspecies
             qyn = qy1+n-1
             qhn = qh1+n-1
             iryn = iry1+n-1
             
!expand             ry_c = q(i,qyn)*sumdrytmp + sumryvtmp*dxinv(1) * &
!expand                  first_deriv_4(q(i-2:i+2,qyn))
             ry_c = q(i,qyn)*sumdrytmp+sumryvtmp * dxinv(1) * &
                  ( D4(1)*(q(i+1,qyn)-q(i-1,qyn)) &
                  + D4(2)*(q(i+2,qyn)-q(i-2,qyn)) )
!expand             ene_c = ry_c*q(i,qhn) + q(i,qyn)*sumryvtmp*dxinv(1)* &
!expand                  first_deriv_4(q(i-2:i+2,qhn))
             ene_c = ry_c*q(i,qhn)+q(i,qyn)*sumryvtmp * dxinv(1) * &
                  ( D4(1)*(q(i+1,qhn)-q(i-1,qhn)) &
                  + D4(2)*(q(i+2,qhn)-q(i-2,qhn)) )
             rhs(i,iene) = rhs(i,iene) - ene_c
             rhs(i,iryn) = rhs(i,iryn) - ry_c
          end do
          
       end if
          
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! use 6th-order stencil for cell lo(1)+3,
       do iface=0,1 
          i = lo(1)+3 + iface
          
!expand          mmtmp6 = matmul(vsp(i-3:i+2), M6)
          mmtmp6(1) = vsp(i-3) * M6(1,1) &
                    + vsp(i-2) * M6(2,1) &
                    + vsp(i-1) * M6(3,1) &
                    + vsp(i  ) * M6(4,1)
          mmtmp6(2) = vsp(i-3) * M6(1,2) &
                    + vsp(i-2) * M6(2,2) &
                    + vsp(i-1) * M6(3,2) &
                    + vsp(i  ) * M6(4,2) &
                    + vsp(i+1) * M6(5,2)
          mmtmp6(3) = vsp(i-3) * M6(1,3) &
                    + vsp(i-2) * M6(2,3) &
                    + vsp(i-1) * M6(3,3) &
                    + vsp(i  ) * M6(4,3) &
                    + vsp(i+1) * M6(5,3) &
                    + vsp(i+2) * M6(6,3)
          mmtmp6(4) = vsp(i-3) * M6(1,4) &
                    + vsp(i-2) * M6(2,4) &
                    + vsp(i-1) * M6(3,4) &
                    + vsp(i  ) * M6(4,4) &
                    + vsp(i+1) * M6(5,4) &
                    + vsp(i+2) * M6(6,4)
          mmtmp6(5) = vsp(i-2) * M6(2,5) &
                    + vsp(i-1) * M6(3,5) &
                    + vsp(i  ) * M6(4,5) &
                    + vsp(i+1) * M6(5,5) &
                    + vsp(i+2) * M6(6,5)
          mmtmp6(6) = vsp(i-1) * M6(3,6) &
                    + vsp(i  ) * M6(4,6) &
                    + vsp(i+1) * M6(5,6) &
                    + vsp(i+2) * M6(6,6)
!expand          Hcell(iface,imx) = dot_product(mmtmp6, q(i-3:i+2,qu))
          Hcell(iface,imx) =  &
               ( q(i-3,qu)*mmtmp6(1) + q(i-2,qu)*mmtmp6(2) &
               + q(i-1,qu)*mmtmp6(3) + q(i  ,qu)*mmtmp6(4) &
               + q(i+1,qu)*mmtmp6(5) + q(i+2,qu)*mmtmp6(6) )
             
!expand          mmtmp6 = matmul(lam(i-3:i+2), M6)
          mmtmp6(1) = lam(i-3) * M6(1,1) &
                    + lam(i-2) * M6(2,1) &
                    + lam(i-1) * M6(3,1) &
                    + lam(i  ) * M6(4,1)
          mmtmp6(2) = lam(i-3) * M6(1,2) &
                    + lam(i-2) * M6(2,2) &
                    + lam(i-1) * M6(3,2) &
                    + lam(i  ) * M6(4,2) &
                    + lam(i+1) * M6(5,2)
          mmtmp6(3) = lam(i-3) * M6(1,3) &
                    + lam(i-2) * M6(2,3) &
                    + lam(i-1) * M6(3,3) &
                    + lam(i  ) * M6(4,3) &
                    + lam(i+1) * M6(5,3) &
                    + lam(i+2) * M6(6,3)
          mmtmp6(4) = lam(i-3) * M6(1,4) &
                    + lam(i-2) * M6(2,4) &
                    + lam(i-1) * M6(3,4) &
                    + lam(i  ) * M6(4,4) &
                    + lam(i+1) * M6(5,4) &
                    + lam(i+2) * M6(6,4)
          mmtmp6(5) = lam(i-2) * M6(2,5) &
                    + lam(i-1) * M6(3,5) &
                    + lam(i  ) * M6(4,5) &
                    + lam(i+1) * M6(5,5) &
                    + lam(i+2) * M6(6,5)
          mmtmp6(6) = lam(i-1) * M6(3,6) &
                    + lam(i  ) * M6(4,6) &
                    + lam(i+1) * M6(5,6) &
                    + lam(i+2) * M6(6,6)
!expand          M6p = matmul(M6,  q(i-3:i+2,qpres))
          M6p(1) = M6T(1,1) * q(i-3,qpres) &
                 + M6T(2,1) * q(i-2,qpres) &
                 + M6T(3,1) * q(i-1,qpres) &
                 + M6T(4,1) * q(i  ,qpres)
          M6p(2) = M6T(1,2) * q(i-3,qpres) &
                 + M6T(2,2) * q(i-2,qpres) &
                 + M6T(3,2) * q(i-1,qpres) &
                 + M6T(4,2) * q(i  ,qpres) &
                 + M6T(5,2) * q(i+1,qpres)
          M6p(3) = M6T(1,3) * q(i-3,qpres) &
                 + M6T(2,3) * q(i-2,qpres) &
                 + M6T(3,3) * q(i-1,qpres) &
                 + M6T(4,3) * q(i  ,qpres) &
                 + M6T(5,3) * q(i+1,qpres) &
                 + M6T(6,3) * q(i+2,qpres)
          M6p(4) = M6T(1,4) * q(i-3,qpres) &
                 + M6T(2,4) * q(i-2,qpres) &
                 + M6T(3,4) * q(i-1,qpres) &
                 + M6T(4,4) * q(i  ,qpres) &
                 + M6T(5,4) * q(i+1,qpres) &
                 + M6T(6,4) * q(i+2,qpres)
          M6p(5) = M6T(2,5) * q(i-2,qpres) &
                 + M6T(3,5) * q(i-1,qpres) &
                 + M6T(4,5) * q(i  ,qpres) &
                 + M6T(5,5) * q(i+1,qpres) &
                 + M6T(6,5) * q(i+2,qpres)
          M6p(6) = M6T(3,6) * q(i-1,qpres) &
                 + M6T(4,6) * q(i  ,qpres) &
                 + M6T(5,6) * q(i+1,qpres) &
                 + M6T(6,6) * q(i+2,qpres)
!expand          Hcell(iface,iene) = dot_product(mmtmp6, q(i-3:i+2,qtemp)) &
!expand               &            + dot_product(      dpe(i-3:i+2), M6p)
          Hcell(iface,iene) =  &
               ( ( q(i-3,qtemp)*mmtmp6(1) + q(i-2,qtemp)*mmtmp6(2) &
                 + q(i-1,qtemp)*mmtmp6(3) + q(i  ,qtemp)*mmtmp6(4) &
                 + q(i+1,qtemp)*mmtmp6(5) + q(i+2,qtemp)*mmtmp6(6) ) &
               + ( dpe(i-3)*M6p(1) + dpe(i-2)*M6p(2) &
                 + dpe(i-1)*M6p(3) + dpe(i  )*M6p(4) &
                 + dpe(i+1)*M6p(5) + dpe(i+2)*M6p(6) ) )
          
          do n = 1, nspecies
             if (n .eq. iias) cycle  ! inactive speices
             qxn = qxy1+n-1
             iryn = iry1+n-1
!expand             M6X = matmul(M6, q(i-3:i+2,qxn))
             M6X(1) = M6T(1,1) * q(i-3,qxn) &
                    + M6T(2,1) * q(i-2,qxn) &
                    + M6T(3,1) * q(i-1,qxn) &
                    + M6T(4,1) * q(i  ,qxn)
             M6X(2) = M6T(1,2) * q(i-3,qxn) &
                    + M6T(2,2) * q(i-2,qxn) &
                    + M6T(3,2) * q(i-1,qxn) &
                    + M6T(4,2) * q(i  ,qxn) &
                    + M6T(5,2) * q(i+1,qxn)
             M6X(3) = M6T(1,3) * q(i-3,qxn) &
                    + M6T(2,3) * q(i-2,qxn) &
                    + M6T(3,3) * q(i-1,qxn) &
                    + M6T(4,3) * q(i  ,qxn) &
                    + M6T(5,3) * q(i+1,qxn) &
                    + M6T(6,3) * q(i+2,qxn)
             M6X(4) = M6T(1,4) * q(i-3,qxn) &
                    + M6T(2,4) * q(i-2,qxn) &
                    + M6T(3,4) * q(i-1,qxn) &
                    + M6T(4,4) * q(i  ,qxn) &
                    + M6T(5,4) * q(i+1,qxn) &
                    + M6T(6,4) * q(i+2,qxn)
             M6X(5) = M6T(2,5) * q(i-2,qxn) &
                    + M6T(3,5) * q(i-1,qxn) &
                    + M6T(4,5) * q(i  ,qxn) &
                    + M6T(5,5) * q(i+1,qxn) &
                    + M6T(6,5) * q(i+2,qxn)
             M6X(6) = M6T(3,6) * q(i-1,qxn) &
                    + M6T(4,6) * q(i  ,qxn) &
                    + M6T(5,6) * q(i+1,qxn) &
                    + M6T(6,6) * q(i+2,qxn)
!expand             Hcell(iface,iene) = Hcell(iface,iene) &
!expand                  +    dot_product(dxe(i-3:i+2,n), M6X)
             Hcell(iface,iene) = Hcell(iface,iene)+ &
                  ( dxe(i-3,n)*M6X(1) + dxe(i-2,n)*M6X(2) &
                  + dxe(i-1,n)*M6X(3) + dxe(i  ,n)*M6X(4) &
                  + dxe(i+1,n)*M6X(5) + dxe(i+2,n)*M6X(6) )
!expand             Hcell(iface,iryn) = dot_product(dpy(i-3:i+2,n), M6p) &
!expand                  +    dot_product(dxy(i-3:i+2,n), M6X)
             Hcell(iface,iryn) =  &
                  ( ( dpy(i-3,n)*M6p(1) + dpy(i-2,n)*M6p(2) &
                    + dpy(i-1,n)*M6p(3) + dpy(i  ,n)*M6p(4) &
                    + dpy(i+1,n)*M6p(5) + dpy(i+2,n)*M6p(6) ) &
                  + ( dxy(i-3,n)*M6X(1) + dxy(i-2,n)*M6X(2) &
                    + dxy(i-1,n)*M6X(3) + dxy(i  ,n)*M6X(4) &
                    + dxy(i+1,n)*M6X(5) + dxy(i+2,n)*M6X(6) ) )
          end do
          
       end do
          
       i = lo(1)+3
       
       do n=2,iene
          rhs(i,n) = rhs(i,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(1)
       end do
       
       sumdrytmp = 0.d0
       do n=iry1,ncons
          if (n.eq.iry_ias) cycle
          sumdrytmp  = sumdrytmp  + (Hcell(1,n) - Hcell(0,n)) * dx2inv(1)
          rhs(i,n) = rhs(i,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(1)
       end do
       
       if (reset_inactive_species) then
          
          rhs(i,iry_ias) = rhs(i,iry_ias) - sumdrytmp
          
       else
          
!expand          gradptmp = dxinv(1) * first_deriv_6(q(i-3:i+3,qpres))
          gradptmp = dxinv(1) * &
               ( D6(1)*(q(i+1,qpres)-q(i-1,qpres)) &
               + D6(2)*(q(i+2,qpres)-q(i-2,qpres)) &
               + D6(3)*(q(i+3,qpres)-q(i-3,qpres)) )
          
          sumryvtmp = 0.d0
          do n = 1, nspecies
             qxn = qxy1+n-1
!expand             sumryvtmp = sumryvtmp + dpy(i,n)*gradptmp  &
!expand                  + dxy(i,n)*dxinv(1)*first_deriv_6(q(i-3:i+3,qxn))
             sumryvtmp = sumryvtmp+dpy(i,n)*gradptmp+dxy(i,n) * dxinv(1) * &
                  ( D6(1)*(q(i+1,qxn)-q(i-1,qxn)) &
                  + D6(2)*(q(i+2,qxn)-q(i-2,qxn)) &
                  + D6(3)*(q(i+3,qxn)-q(i-3,qxn)) )
          end do
             
          do n=1,nspecies
             qyn = qy1+n-1
             qhn = qh1+n-1
             iryn = iry1+n-1
             
!expand             ry_c = q(i,qyn)*sumdrytmp + sumryvtmp*dxinv(1) * &
!expand                  first_deriv_6(q(i-3:i+3,qyn))
             ry_c = q(i,qyn)*sumdrytmp+sumryvtmp * dxinv(1) * &
                  ( D6(1)*(q(i+1,qyn)-q(i-1,qyn)) &
                  + D6(2)*(q(i+2,qyn)-q(i-2,qyn)) &
                  + D6(3)*(q(i+3,qyn)-q(i-3,qyn)) )
!expand             ene_c = ry_c*q(i,qhn) + q(i,qyn)*sumryvtmp*dxinv(1)* &
!expand                  first_deriv_6(q(i-3:i+3,qhn))
             ene_c = ry_c*q(i,qhn)+q(i,qyn)*sumryvtmp * dxinv(1) * &
                  ( D6(1)*(q(i+1,qhn)-q(i-1,qhn)) &
                  + D6(2)*(q(i+2,qhn)-q(i-2,qhn)) &
                  + D6(3)*(q(i+3,qhn)-q(i-3,qhn)) )
             rhs(i,iene) = rhs(i,iene) - ene_c
             rhs(i,iryn) = rhs(i,iryn) - ry_c
          end do
             
       end if

    end if

    !
    ! hi-x boundary
    !
    if (physbchi(1)) then
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! use 6th-order stencil for cell hi(1)-3,
       do iface=0,1  ! two faces of 
          i = hi(1)-3 + iface
          
!expand          mmtmp6 = matmul(vsp(i-3:i+2), M6)
          mmtmp6(1) = vsp(i-3) * M6(1,1) &
                    + vsp(i-2) * M6(2,1) &
                    + vsp(i-1) * M6(3,1) &
                    + vsp(i  ) * M6(4,1)
          mmtmp6(2) = vsp(i-3) * M6(1,2) &
                    + vsp(i-2) * M6(2,2) &
                    + vsp(i-1) * M6(3,2) &
                    + vsp(i  ) * M6(4,2) &
                    + vsp(i+1) * M6(5,2)
          mmtmp6(3) = vsp(i-3) * M6(1,3) &
                    + vsp(i-2) * M6(2,3) &
                    + vsp(i-1) * M6(3,3) &
                    + vsp(i  ) * M6(4,3) &
                    + vsp(i+1) * M6(5,3) &
                    + vsp(i+2) * M6(6,3)
          mmtmp6(4) = vsp(i-3) * M6(1,4) &
                    + vsp(i-2) * M6(2,4) &
                    + vsp(i-1) * M6(3,4) &
                    + vsp(i  ) * M6(4,4) &
                    + vsp(i+1) * M6(5,4) &
                    + vsp(i+2) * M6(6,4)
          mmtmp6(5) = vsp(i-2) * M6(2,5) &
                    + vsp(i-1) * M6(3,5) &
                    + vsp(i  ) * M6(4,5) &
                    + vsp(i+1) * M6(5,5) &
                    + vsp(i+2) * M6(6,5)
          mmtmp6(6) = vsp(i-1) * M6(3,6) &
                    + vsp(i  ) * M6(4,6) &
                    + vsp(i+1) * M6(5,6) &
                    + vsp(i+2) * M6(6,6)
!expand          Hcell(iface,imx) = dot_product(mmtmp6, q(i-3:i+2,qu))
          Hcell(iface,imx) =  &
               ( q(i-3,qu)*mmtmp6(1) + q(i-2,qu)*mmtmp6(2) &
               + q(i-1,qu)*mmtmp6(3) + q(i  ,qu)*mmtmp6(4) &
               + q(i+1,qu)*mmtmp6(5) + q(i+2,qu)*mmtmp6(6) )
             
!expand          mmtmp6 = matmul(lam(i-3:i+2), M6)
          mmtmp6(1) = lam(i-3) * M6(1,1) &
                    + lam(i-2) * M6(2,1) &
                    + lam(i-1) * M6(3,1) &
                    + lam(i  ) * M6(4,1)
          mmtmp6(2) = lam(i-3) * M6(1,2) &
                    + lam(i-2) * M6(2,2) &
                    + lam(i-1) * M6(3,2) &
                    + lam(i  ) * M6(4,2) &
                    + lam(i+1) * M6(5,2)
          mmtmp6(3) = lam(i-3) * M6(1,3) &
                    + lam(i-2) * M6(2,3) &
                    + lam(i-1) * M6(3,3) &
                    + lam(i  ) * M6(4,3) &
                    + lam(i+1) * M6(5,3) &
                    + lam(i+2) * M6(6,3)
          mmtmp6(4) = lam(i-3) * M6(1,4) &
                    + lam(i-2) * M6(2,4) &
                    + lam(i-1) * M6(3,4) &
                    + lam(i  ) * M6(4,4) &
                    + lam(i+1) * M6(5,4) &
                    + lam(i+2) * M6(6,4)
          mmtmp6(5) = lam(i-2) * M6(2,5) &
                    + lam(i-1) * M6(3,5) &
                    + lam(i  ) * M6(4,5) &
                    + lam(i+1) * M6(5,5) &
                    + lam(i+2) * M6(6,5)
          mmtmp6(6) = lam(i-1) * M6(3,6) &
                    + lam(i  ) * M6(4,6) &
                    + lam(i+1) * M6(5,6) &
                    + lam(i+2) * M6(6,6)
!expand          M6p = matmul(M6,  q(i-3:i+2,qpres))
          M6p(1) = M6T(1,1) * q(i-3,qpres) &
                 + M6T(2,1) * q(i-2,qpres) &
                 + M6T(3,1) * q(i-1,qpres) &
                 + M6T(4,1) * q(i  ,qpres)
          M6p(2) = M6T(1,2) * q(i-3,qpres) &
                 + M6T(2,2) * q(i-2,qpres) &
                 + M6T(3,2) * q(i-1,qpres) &
                 + M6T(4,2) * q(i  ,qpres) &
                 + M6T(5,2) * q(i+1,qpres)
          M6p(3) = M6T(1,3) * q(i-3,qpres) &
                 + M6T(2,3) * q(i-2,qpres) &
                 + M6T(3,3) * q(i-1,qpres) &
                 + M6T(4,3) * q(i  ,qpres) &
                 + M6T(5,3) * q(i+1,qpres) &
                 + M6T(6,3) * q(i+2,qpres)
          M6p(4) = M6T(1,4) * q(i-3,qpres) &
                 + M6T(2,4) * q(i-2,qpres) &
                 + M6T(3,4) * q(i-1,qpres) &
                 + M6T(4,4) * q(i  ,qpres) &
                 + M6T(5,4) * q(i+1,qpres) &
                 + M6T(6,4) * q(i+2,qpres)
          M6p(5) = M6T(2,5) * q(i-2,qpres) &
                 + M6T(3,5) * q(i-1,qpres) &
                 + M6T(4,5) * q(i  ,qpres) &
                 + M6T(5,5) * q(i+1,qpres) &
                 + M6T(6,5) * q(i+2,qpres)
          M6p(6) = M6T(3,6) * q(i-1,qpres) &
                 + M6T(4,6) * q(i  ,qpres) &
                 + M6T(5,6) * q(i+1,qpres) &
                 + M6T(6,6) * q(i+2,qpres)
!expand          Hcell(iface,iene) = dot_product(mmtmp6, q(i-3:i+2,qtemp)) &
!expand               &            + dot_product(      dpe(i-3:i+2), M6p)
          Hcell(iface,iene) =  &
               ( ( q(i-3,qtemp)*mmtmp6(1) + q(i-2,qtemp)*mmtmp6(2) &
                 + q(i-1,qtemp)*mmtmp6(3) + q(i  ,qtemp)*mmtmp6(4) &
                 + q(i+1,qtemp)*mmtmp6(5) + q(i+2,qtemp)*mmtmp6(6) ) &
               + ( dpe(i-3)*M6p(1) + dpe(i-2)*M6p(2) &
                 + dpe(i-1)*M6p(3) + dpe(i  )*M6p(4) &
                 + dpe(i+1)*M6p(5) + dpe(i+2)*M6p(6) ) )
             
          do n = 1, nspecies
             if (n .eq. iias) cycle  ! inactive speices
             qxn = qxy1+n-1
             iryn = iry1+n-1
!expand             M6X = matmul(M6, q(i-3:i+2,qxn))
             M6X(1) = M6T(1,1) * q(i-3,qxn) &
                    + M6T(2,1) * q(i-2,qxn) &
                    + M6T(3,1) * q(i-1,qxn) &
                    + M6T(4,1) * q(i  ,qxn)
             M6X(2) = M6T(1,2) * q(i-3,qxn) &
                    + M6T(2,2) * q(i-2,qxn) &
                    + M6T(3,2) * q(i-1,qxn) &
                    + M6T(4,2) * q(i  ,qxn) &
                    + M6T(5,2) * q(i+1,qxn)
             M6X(3) = M6T(1,3) * q(i-3,qxn) &
                    + M6T(2,3) * q(i-2,qxn) &
                    + M6T(3,3) * q(i-1,qxn) &
                    + M6T(4,3) * q(i  ,qxn) &
                    + M6T(5,3) * q(i+1,qxn) &
                    + M6T(6,3) * q(i+2,qxn)
             M6X(4) = M6T(1,4) * q(i-3,qxn) &
                    + M6T(2,4) * q(i-2,qxn) &
                    + M6T(3,4) * q(i-1,qxn) &
                    + M6T(4,4) * q(i  ,qxn) &
                    + M6T(5,4) * q(i+1,qxn) &
                    + M6T(6,4) * q(i+2,qxn)
             M6X(5) = M6T(2,5) * q(i-2,qxn) &
                    + M6T(3,5) * q(i-1,qxn) &
                    + M6T(4,5) * q(i  ,qxn) &
                    + M6T(5,5) * q(i+1,qxn) &
                    + M6T(6,5) * q(i+2,qxn)
             M6X(6) = M6T(3,6) * q(i-1,qxn) &
                    + M6T(4,6) * q(i  ,qxn) &
                    + M6T(5,6) * q(i+1,qxn) &
                    + M6T(6,6) * q(i+2,qxn)
!expand             Hcell(iface,iene) = Hcell(iface,iene) &
!expand                  +    dot_product(dxe(i-3:i+2,n), M6X)                
             Hcell(iface,iene) = Hcell(iface,iene)+ &
                  ( dxe(i-3,n)*M6X(1) + dxe(i-2,n)*M6X(2) &
                  + dxe(i-1,n)*M6X(3) + dxe(i  ,n)*M6X(4) &
                  + dxe(i+1,n)*M6X(5) + dxe(i+2,n)*M6X(6) )
!expand             Hcell(iface,iryn) = dot_product(dpy(i-3:i+2,n), M6p) &
!expand                  +    dot_product(dxy(i-3:i+2,n), M6X)
             Hcell(iface,iryn) =  &
                  ( ( dpy(i-3,n)*M6p(1) + dpy(i-2,n)*M6p(2) &
                    + dpy(i-1,n)*M6p(3) + dpy(i  ,n)*M6p(4) &
                    + dpy(i+1,n)*M6p(5) + dpy(i+2,n)*M6p(6) ) &
                  + ( dxy(i-3,n)*M6X(1) + dxy(i-2,n)*M6X(2) &
                    + dxy(i-1,n)*M6X(3) + dxy(i  ,n)*M6X(4) &
                    + dxy(i+1,n)*M6X(5) + dxy(i+2,n)*M6X(6) ) )
          end do
             
       end do
          
       i = hi(1)-3
          
       do n=2,iene
          rhs(i,n) = rhs(i,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(1)
       end do
       
       sumdrytmp = 0.d0
       do n=iry1,ncons
          if (n.eq.iry_ias) cycle
          sumdrytmp  = sumdrytmp  + (Hcell(1,n) - Hcell(0,n)) * dx2inv(1)
          rhs(i,n) = rhs(i,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(1)
       end do
       
       if (reset_inactive_species) then
          
          rhs(i,iry_ias) = rhs(i,iry_ias) - sumdrytmp
          
       else
          
!expand          gradptmp = dxinv(1) * first_deriv_6(q(i-3:i+3,qpres))
          gradptmp = dxinv(1) * &
               ( D6(1)*(q(i+1,qpres)-q(i-1,qpres)) &
               + D6(2)*(q(i+2,qpres)-q(i-2,qpres)) &
               + D6(3)*(q(i+3,qpres)-q(i-3,qpres)) )
          
          sumryvtmp = 0.d0
          do n = 1, nspecies
             qxn = qxy1+n-1
!expand             sumryvtmp = sumryvtmp + dpy(i,n)*gradptmp  &
!expand                  + dxy(i,n)*dxinv(1)*first_deriv_6(q(i-3:i+3,qxn))
             sumryvtmp = sumryvtmp+dpy(i,n)*gradptmp+dxy(i,n) * dxinv(1) * &
                  ( D6(1)*(q(i+1,qxn)-q(i-1,qxn)) &
                  + D6(2)*(q(i+2,qxn)-q(i-2,qxn)) &
                  + D6(3)*(q(i+3,qxn)-q(i-3,qxn)) )
          end do
          
          do n=1,nspecies
             qyn = qy1+n-1
             qhn = qh1+n-1
             iryn = iry1+n-1
             
!expand             ry_c = q(i,qyn)*sumdrytmp + sumryvtmp*dxinv(1) * &
!expand                  first_deriv_6(q(i-3:i+3,qyn))
             ry_c = q(i,qyn)*sumdrytmp+sumryvtmp * dxinv(1) * &
                  ( D6(1)*(q(i+1,qyn)-q(i-1,qyn)) &
                  + D6(2)*(q(i+2,qyn)-q(i-2,qyn)) &
                  + D6(3)*(q(i+3,qyn)-q(i-3,qyn)) )
!expand             ene_c = ry_c*q(i,qhn) + q(i,qyn)*sumryvtmp*dxinv(1)* &
!expand                  first_deriv_6(q(i-3:i+3,qhn))
             ene_c = ry_c*q(i,qhn)+q(i,qyn)*sumryvtmp * dxinv(1) * &
                  ( D6(1)*(q(i+1,qhn)-q(i-1,qhn)) &
                  + D6(2)*(q(i+2,qhn)-q(i-2,qhn)) &
                  + D6(3)*(q(i+3,qhn)-q(i-3,qhn)) )
             rhs(i,iene) = rhs(i,iene) - ene_c
             rhs(i,iryn) = rhs(i,iryn) - ry_c
          end do
          
       end if

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! use 4th-order stencil for cell hi(1)-2,
       do iface=0,1 
          i = hi(1)-2 + iface
          
!expand          mmtmp4 = matmul(vsp(i-2:i+1), M4)
          mmtmp4(1) = vsp(i-2) * M4(1,1) &
                    + vsp(i-1) * M4(2,1) &
                    + vsp(i  ) * M4(3,1)
          mmtmp4(2) = vsp(i-2) * M4(1,2) &
                    + vsp(i-1) * M4(2,2) &
                    + vsp(i  ) * M4(3,2) &
                    + vsp(i+1) * M4(4,2)
          mmtmp4(3) = vsp(i-2) * M4(1,3) &
                    + vsp(i-1) * M4(2,3) &
                    + vsp(i  ) * M4(3,3) &
                    + vsp(i+1) * M4(4,3)
          mmtmp4(4) = vsp(i-1) * M4(2,4) &
                    + vsp(i  ) * M4(3,4) &
                    + vsp(i+1) * M4(4,4)
!expand          Hcell(iface,imx) = dot_product(mmtmp4, q(i-2:i+1,qu))
          Hcell(iface,imx) =  &
               ( q(i-2,qu)*mmtmp4(1) + q(i-1,qu)*mmtmp4(2) &
               + q(i  ,qu)*mmtmp4(3) + q(i+1,qu)*mmtmp4(4) )
             
!expand          mmtmp4 = matmul(lam(i-2:i+1), M4)
          mmtmp4(1) = lam(i-2) * M4(1,1) &
                    + lam(i-1) * M4(2,1) &
                    + lam(i  ) * M4(3,1)
          mmtmp4(2) = lam(i-2) * M4(1,2) &
                    + lam(i-1) * M4(2,2) &
                    + lam(i  ) * M4(3,2) &
                    + lam(i+1) * M4(4,2)
          mmtmp4(3) = lam(i-2) * M4(1,3) &
                    + lam(i-1) * M4(2,3) &
                    + lam(i  ) * M4(3,3) &
                    + lam(i+1) * M4(4,3)
          mmtmp4(4) = lam(i-1) * M4(2,4) &
                    + lam(i  ) * M4(3,4) &
                    + lam(i+1) * M4(4,4)
!expand          M4p = matmul(M4,  q(i-2:i+1,qpres))
          M4p(1) = M4T(1,1) * q(i-2,qpres) &
                 + M4T(2,1) * q(i-1,qpres) &
                 + M4T(3,1) * q(i  ,qpres)
          M4p(2) = M4T(1,2) * q(i-2,qpres) &
                 + M4T(2,2) * q(i-1,qpres) &
                 + M4T(3,2) * q(i  ,qpres) &
                 + M4T(4,2) * q(i+1,qpres)
          M4p(3) = M4T(1,3) * q(i-2,qpres) &
                 + M4T(2,3) * q(i-1,qpres) &
                 + M4T(3,3) * q(i  ,qpres) &
                 + M4T(4,3) * q(i+1,qpres)
          M4p(4) = M4T(2,4) * q(i-1,qpres) &
                 + M4T(3,4) * q(i  ,qpres) &
                 + M4T(4,4) * q(i+1,qpres)
!expand          Hcell(iface,iene) = dot_product(mmtmp4, q(i-2:i+1,qtemp)) &
!expand               &            + dot_product(      dpe(i-2:i+1), M4p)
          Hcell(iface,iene) =  &
               ( ( q(i-2,qtemp)*mmtmp4(1) + q(i-1,qtemp)*mmtmp4(2) &
                 + q(i  ,qtemp)*mmtmp4(3) + q(i+1,qtemp)*mmtmp4(4) ) &
               + ( dpe(i-2)*M4p(1) + dpe(i-1)*M4p(2) &
                 + dpe(i  )*M4p(3) + dpe(i+1)*M4p(4) ) )
             
          do n = 1, nspecies
             if (n .eq. iias) cycle  ! inactive speices
             qxn = qxy1+n-1
             iryn = iry1+n-1
             
!expand             M4X = matmul(M4, q(i-2:i+1,qxn))
             M4X(1) = M4T(1,1) * q(i-2,qxn) &
                    + M4T(2,1) * q(i-1,qxn) &
                    + M4T(3,1) * q(i  ,qxn)
             M4X(2) = M4T(1,2) * q(i-2,qxn) &
                    + M4T(2,2) * q(i-1,qxn) &
                    + M4T(3,2) * q(i  ,qxn) &
                    + M4T(4,2) * q(i+1,qxn)
             M4X(3) = M4T(1,3) * q(i-2,qxn) &
                    + M4T(2,3) * q(i-1,qxn) &
                    + M4T(3,3) * q(i  ,qxn) &
                    + M4T(4,3) * q(i+1,qxn)
             M4X(4) = M4T(2,4) * q(i-1,qxn) &
                    + M4T(3,4) * q(i  ,qxn) &
                    + M4T(4,4) * q(i+1,qxn)
!expand             Hcell(iface,iene) = Hcell(iface,iene) &
!expand                  +    dot_product(dxe(i-2:i+1,n), M4X)
             Hcell(iface,iene) = Hcell(iface,iene)+ &
                  ( dxe(i-2,n)*M4X(1) + dxe(i-1,n)*M4X(2) &
                  + dxe(i  ,n)*M4X(3) + dxe(i+1,n)*M4X(4) )
!expand             Hcell(iface,iryn) = dot_product(dpy(i-2:i+1,n), M4p) &
!expand                  +    dot_product(dxy(i-2:i+1,n), M4X)
             Hcell(iface,iryn) =  &
                  ( ( dpy(i-2,n)*M4p(1) + dpy(i-1,n)*M4p(2) &
                    + dpy(i  ,n)*M4p(3) + dpy(i+1,n)*M4p(4) ) &
                  + ( dxy(i-2,n)*M4X(1) + dxy(i-1,n)*M4X(2) &
                    + dxy(i  ,n)*M4X(3) + dxy(i+1,n)*M4X(4) ) )
          end do
          
       end do
       
       i = hi(1)-2
       
       do n=2,iene
          rhs(i,n) = rhs(i,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(1)
       end do
       
       sumdrytmp = 0.d0
       do n=iry1,ncons
          if (n.eq.iry_ias) cycle
          sumdrytmp  = sumdrytmp  + (Hcell(1,n) - Hcell(0,n)) * dx2inv(1)
          rhs(i,n) = rhs(i,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(1)
       end do
       
       if (reset_inactive_species) then
          
          rhs(i,iry_ias) = rhs(i,iry_ias) - sumdrytmp
          
       else
          
!expand          gradptmp = dxinv(1) * first_deriv_4(q(i-2:i+2,qpres))
          gradptmp = dxinv(1) * &
               ( D4(1)*(q(i+1,qpres)-q(i-1,qpres)) &
               + D4(2)*(q(i+2,qpres)-q(i-2,qpres)) )
          
          sumryvtmp = 0.d0
          do n = 1, nspecies
             if (n.eq.iias) cycle
             qxn = qxy1+n-1
!expand             sumryvtmp = sumryvtmp + dpy(i,n)*gradptmp  &
!expand                  + dxy(i,n)*dxinv(1)*first_deriv_4(q(i-2:i+2,qxn))
             sumryvtmp = sumryvtmp+dpy(i,n)*gradptmp+dxy(i,n) * dxinv(1) * &
                  ( D4(1)*(q(i+1,qxn)-q(i-1,qxn)) &
                  + D4(2)*(q(i+2,qxn)-q(i-2,qxn)) )
          end do
          
          do n=1,nspecies
             qyn = qy1+n-1
             qhn = qh1+n-1
             iryn = iry1+n-1
             
!expand             ry_c = q(i,qyn)*sumdrytmp + sumryvtmp*dxinv(1) * &
!expand                  first_deriv_4(q(i-2:i+2,qyn))
             ry_c = q(i,qyn)*sumdrytmp+sumryvtmp * dxinv(1) * &
                  ( D4(1)*(q(i+1,qyn)-q(i-1,qyn)) &
                  + D4(2)*(q(i+2,qyn)-q(i-2,qyn)) )
!expand             ene_c = ry_c*q(i,qhn) + q(i,qyn)*sumryvtmp*dxinv(1)* &
!expand                  first_deriv_4(q(i-2:i+2,qhn))
             ene_c = ry_c*q(i,qhn)+q(i,qyn)*sumryvtmp * dxinv(1) * &
                  ( D4(1)*(q(i+1,qhn)-q(i-1,qhn)) &
                  + D4(2)*(q(i+2,qhn)-q(i-2,qhn)) )
             rhs(i,iene) = rhs(i,iene) - ene_c
             rhs(i,iryn) = rhs(i,iryn) - ry_c
          end do
          
       end if
       
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! use 2nd-order stencil for cell hi(1)-1,
       do iface=0,1 
          i = hi(1)-1 + iface
          
          mmtmp2 = matmul(vsp(i-1:i), M2)
!expand          Hcell(iface,imx) = dot_product(mmtmp2, q(i-1:i,qu))
          Hcell(iface,imx) =  &
               ( q(i-1,qu)*mmtmp2(1) + q(i  ,qu)*mmtmp2(2) )
             
          mmtmp2 = matmul(lam(i-1:i), M2)
          M2p = matmul(M2,  q(i-1:i,qpres))
!expand          Hcell(iface,iene) = dot_product(mmtmp2, q(i-1:i,qtemp)) &
!expand               &            + dot_product(      dpe(i-1:i), M2p)
          Hcell(iface,iene) =  &
               ( ( q(i-1,qtemp)*mmtmp2(1) + q(i  ,qtemp)*mmtmp2(2) ) &
               + ( dpe(i-1)*M2p(1) + dpe(i  )*M2p(2) ) )
             
          Htot = 0.d0
          do n = 1, nspecies
 
             if (n .eq. iias) cycle  ! inactive speices
             
             qxn = qxy1+n-1
             qyn = qy1+n-1
             
             M2X = matmul(M2, q(i-1:i,qxn))
             
!expand             Htmp(n) = dot_product(dpy(i-1:i,n), M2p) &
!expand                  +    dot_product(dxy(i-1:i,n), M2X)
             Htmp(n) =  &
                  ( ( dpy(i-1,n)*M2p(1) + dpy(i  ,n)*M2p(2) ) &
                  + ( dxy(i-1,n)*M2X(1) + dxy(i  ,n)*M2X(2) ) )
             
!expand             Hcell(iface,iene) = Hcell(iface,iene) &
!expand                  +    dot_product(dxe(i-1:i,n), M2X)
             Hcell(iface,iene) = Hcell(iface,iene)+ &
                  ( dxe(i-1,n)*M2X(1) + dxe(i  ,n)*M2X(2) )
             
             Htot = Htot + Htmp(n)
             Ytmp(n) = 0.5d0*(q(i-1,qyn) + q(i,qyn))
             
             Hcell(iface,iry1+n-1) = Htmp(n)
          end do
          
          if (.not. reset_inactive_species) then
             do n = 1, nspecies
                Hcell(iface,iry1+n-1) = Hcell(iface,iry1+n-1) - Ytmp(n)*Htot
             end do
             
             do n = 1, nspecies
                qhn = qh1+n-1
                hhalf = 0.5d0*(q(i-1,qhn) + q(i,qhn))
                Hcell(iface,iene) =  Hcell(iface,iene) - Ytmp(n)*Htot*hhalf
             end do
          end if
       end do
       
       i = hi(1)-1
       if (reset_inactive_species) then
          do n=2,iene
             rhs(i,n) = rhs(i,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(1)
          end do
          sumdrytmp = 0.d0
          do n=iry1,ncons
             if (n.eq.iry_ias) cycle
             sumdrytmp  = sumdrytmp  + (Hcell(1,n) - Hcell(0,n)) * dx2inv(1)
             rhs(i,n) = rhs(i,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(1)
          end do
          rhs(i,iry_ias) = rhs(i,iry_ias) - sumdrytmp
       else
          do n=2,ncons
             rhs(i,n) = rhs(i,n) + (Hcell(1,n) - Hcell(0,n)) * dx2inv(1)
          end do
       end if
       
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       i = hi(1)
       ! use completely left-biased stencil
       mmtmpB = matmul(vsp(i-3:i), BLB)
!expand       rhs(i,imx) = rhs(i,imx)+finhi(1)*dx2inv(1)*dot_product(mmtmpB,q(i-3:i,qu))
       rhs(i,imx) = rhs(i,imx)+finhi(1)*dx2inv(1)* &
            ( q(i-3,qu)*mmtmpB(1) + q(i-2,qu)*mmtmpB(2) &
            + q(i-1,qu)*mmtmpB(3) + q(i  ,qu)*mmtmpB(4) )
       
       mmtmpB = matmul(lam(i-3:i), BLB)
       BBp = matmul(BLB, q(i-3:i,qpres))
!expand       rhs(i,iene) = rhs(i,iene) + fouhi(1)*dx2inv(1) * &
!expand            ( dot_product(mmtmpB, q(i-3:i,qtemp)) &
!expand            + dot_product(      dpe(i-3:i), BBp) )
       rhs(i,iene) = rhs(i,iene)+fouhi(1)*dx2inv(1)* &
            ( ( q(i-3,qtemp)*mmtmpB(1) + q(i-2,qtemp)*mmtmpB(2) &
              + q(i-1,qtemp)*mmtmpB(3) + q(i  ,qtemp)*mmtmpB(4) ) &
            + ( dpe(i-3)*BBp(1) + dpe(i-2)*BBp(2) &
              + dpe(i-1)*BBp(3) + dpe(i  )*BBp(4) ) )
       
       rhstot = 0.d0
       rhsene = 0.d0
       do n = 1, nspecies
          
          if (n .eq. iias) cycle  ! inactive speices
          
          qxn = qxy1+n-1
          qyn = qy1+n-1
          
          BBX = matmul(BLB, q(i-3:i,qxn))
          
!expand          rhstmp(n) = dot_product(dpy(i-3:i,n), BBp) &
!expand               +      dot_product(dxy(i-3:i,n), BBX)
          rhstmp(n) =  &
               ( ( dpy(i-3,n)*BBp(1) + dpy(i-2,n)*BBp(2) &
                 + dpy(i-1,n)*BBp(3) + dpy(i  ,n)*BBp(4) ) &
               + ( dxy(i-3,n)*BBX(1) + dxy(i-2,n)*BBX(2) &
                 + dxy(i-1,n)*BBX(3) + dxy(i  ,n)*BBX(4) ) )
          
!expand          rhsene = rhsene &
!expand               +      dot_product(dxe(i-3:i,n), BBX)
          rhsene = rhsene+ &
               ( dxe(i-3,n)*BBX(1) + dxe(i-2,n)*BBX(2) &
               + dxe(i-1,n)*BBX(3) + dxe(i  ,n)*BBX(4) )
          
          rhstot = rhstot + rhstmp(n)
          Ytmp(n) = q(i,qyn)
          
          rhs(i,iry1+n-1) =  rhs(i,iry1+n-1) + fouhi(1)*dx2inv(1)*rhstmp(n)
       end do

       rhs(i,iene) = rhs(i,iene) + fouhi(1)*dx2inv(1) * rhsene
       
       if (reset_inactive_species) then
          
          rhs(i,iry_ias) = rhs(i,iry_ias) - fouhi(1)*dx2inv(1)*rhstot
          
       else
          
          do n = 1, nspecies
             rhs(i,iry1+n-1) =  rhs(i,iry1+n-1) - &
                  fouhi(1)*dx2inv(1) * Ytmp(n)*rhstot
          end do
          
          rhsene = 0.d0
          do n = 1, nspecies
             qhn = qh1+n-1
             rhsene = rhsene - Ytmp(n) * q(i,qhn) * rhstot
          end do
          rhs(i,iene) = rhs(i,iene) + fouhi(1)*dx2inv(1) * rhsene
       end if

    end if

    !
    ! add kinetic energy
    !
    do i=lo(1),hi(1)
       rhs(i,iene) = rhs(i,iene) + rhs(i,imx)*q(i,qu)
    end do

    deallocate(Hg,dpy,dxe,dpe,vsp,M8p)
    deallocate(sumdrY,sumryv,gradp)

    if (diff_gradY) then
       deallocate(dxy)
    else
       nullify(dxy)
    end if

  end subroutine diffterm_2


  subroutine chemterm_1d(lo,hi,q,qlo,qhi,up,uplo,uphi,upc,upclo,upchi,dt)
    use probin_module, only : use_vode, vode_always_new_j
    use vode_module, only : verbose, itol, rtol, atol, vode_MF=>MF, &
         voderwork, vodeiwork, lvoderwork, lvodeiwork, voderpar, vodeipar

    double precision, intent(in) :: dt
    integer,         intent(in):: lo(1),hi(1),qlo(1),qhi(1),uplo(1),uphi(1),upclo(1),upchi(1)
    double precision,intent(in):: q  (  qlo(1):  qhi(1),nprim)
    double precision           :: up ( uplo(1): uphi(1),ncons)
    double precision           :: upc(upclo(1):upchi(1),ncons)

    integer :: iwrk, i,n,np, iryn
    double precision :: Yt(lo(1):hi(1),nspecies), wdot(lo(1):hi(1),nspecies), rwrk
    double precision :: YTvode(nspecies+1), time, dtinv

    external f_jac, f_rhs, dvode

    ! vode stuff
    integer, parameter :: itask=1, iopt=1
    integer :: MF, istate, ifail

    if (use_vode) then

       dtinv = 1.d0/dt

       call setfirst(.true.)

       do i=lo(1),hi(1)

          voderpar(1) = q(i,qrho)

          YTvode(1:nspecies) = q(i,qy1:qy1+nspecies-1)
          YTvode(nspecies+1) = q(i,qtemp)
       
          istate = 1
          
          time = 0.d0
          
          if (vode_always_new_j) call setfirst(.true.)
          
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
             iryn = iry1+n-1
             upc(i,iryn) = dtinv*q(i,qrho)*(YTvode(n)-q(i,qy1+n-1))
          end do
          
       end do

       do n=1, nspecies
          iryn = iry1+n-1
          do i=lo(1),hi(1)
             up(i,iryn) = upc(i,iryn)                
          end do
       end do

    else

       np = hi(1) - lo(1) + 1
       
       do n=1, nspecies
          do i=lo(1),hi(1)
             Yt(i,n) = q(i,qy1+n-1)
          end do
       end do
          
       call vckwyr(np, q(lo(1),qrho), q(lo(1),qtemp), Yt, iwrk, rwrk, wdot)
          
       do n=1, nspecies
          iryn = iry1+n-1
          do i=lo(1),hi(1)
             upc(i,iryn) = wdot(i,n) * molecular_weight(n)
             up (i,iryn) = upc(i,iryn)
          end do
       end do
       
    end if

  end subroutine chemterm_1d


  subroutine comp_courno_1d(lo,hi,dx,Q,qlo,qhi,courno)
    integer, intent(in) :: lo(1), hi(1), qlo(1), qhi(1)
    double precision, intent(in) :: dx(1)
    double precision, intent(in) :: q(qlo(1):qhi(1),nprim)
    double precision, intent(inout) :: courno

    integer :: i, iwrk
    double precision :: dxinv(1), c, rwrk, Cv, Cp
    double precision :: Tt, X(nspecies), gamma
    double precision :: courx

    dxinv(1) = 1.0d0 / dx(1)

    do i=lo(1),hi(1)
          
       Tt = q(i,qtemp)
       X  = q(i,qx1:qx1+nspecies-1)
       call ckcvbl(Tt, X, iwrk, rwrk, Cv)
       Cp = Cv + Ru
       gamma = Cp / Cv
       c = sqrt(gamma*q(i,qpres)/q(i,qrho))
       
       courx = (c+abs(q(i,qu))) * dxinv(1)
             
       courno = max( courx, courno )

    end do

  end subroutine comp_courno_1d

end module kernels_1d_module
