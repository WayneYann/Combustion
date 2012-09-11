module init_data_module

  use multifab_module
  use bl_constants_module, only : Pi=>M_PI

  implicit none

  private

  public :: init_data

contains

  subroutine init_data(data,dx,plo,phi)

    type(multifab),   intent(inout) :: data
    double precision, intent(in   ) :: dx(data%dim)
    double precision, intent(in   ) :: plo(data%dim), phi(data%dim)

    integer                   :: lo(data%dim), hi(data%dim), ng, i
    double precision, pointer :: dp(:,:,:,:)

    ng = data%ng

    do i=1,nfabs(data)
       dp => dataptr(data,i)
       lo = lwb(get_box(data,i))
       hi = upb(get_box(data,i))

       select case(data%dim)
       case (2)
          call bl_error('We only support 3-D')
       case (3)
          call init_data_3d(lo,hi,ng,dx,dp,plo,phi)
       end select
    end do

  end subroutine init_data

  subroutine init_data_3d(lo,hi,ng,dx,cons,phlo,phhi)

    use variables_module, only : irho, imx,imy,imz,iene,iry1,ncons
    use chemistry_module, only : nspecies
    use probin_module,    only : prob_type, pertmag, rfire, Tinit, uinit, vinit, winit

    integer,          intent(in   ) :: lo(3),hi(3),ng
    double precision, intent(in   ) :: dx(3),phlo(3),phhi(3)
    double precision, intent(inout) :: cons(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,ncons)

    integer          :: i,j,k,n
    double precision :: x, y, z, r

    double precision pmf_vals(nspecies+3)
    double precision Xt(nspecies), Yt(nspecies)
    double precision rhot,u1t,u2t,u3t,Tt,et
    integer :: iwrk
    double precision :: rwrk

    double precision, parameter :: patmos = 1.01325d6

    integer, parameter :: lmodemin=14, lmodemax=20
    double precision :: alphalm(lmodemin:lmodemax,0:lmodemax)
    double precision ::  betalm(lmodemin:lmodemax,0:lmodemax)
    double precision :: gammalm(lmodemin:lmodemax,0:lmodemax)
    double precision :: rfront, phi, theta, xtemp, xloc, yloc, zloc
    integer :: l, m, ctr, seed_size
    integer, allocatable :: seed(:)

    call random_seed(size=seed_size)
    allocate(seed(seed_size))
    seed = 134527
    call random_seed(put=seed)
    do i = lmodemin,lmodemax
       do j = 0,i
          call random_number(alphalm(i,j))
          call random_number(betalm (i,j))
          call random_number(gammalm(i,j))
       end do
    end do
    deallocate(seed)

    !$omp parallel do &
    !$omp private(i,j,k,n,x,y,z,r,pmf_vals) &
    !$omp private(Xt,Yt,rhot,u1t,u2t,u3t,Tt,et,iwrk,rwrk) &
    !$omp private(rfront,phi,theta,xtemp,xloc,yloc,zloc,l,m,ctr)
    do k=lo(3),hi(3)
       z = phlo(3) + dx(3)*(k + 0.5d0)
       do j=lo(2),hi(2)
          y = phlo(2) + dx(2)*(j + 0.5d0)
          do i=lo(1),hi(1)
             x = phlo(1) + dx(1)*(i + 0.5d0)

             if (prob_type .eq. 1) then
                r = sqrt(x**2+y**2+z**2)

                if (pertmag .gt. 0.d0) then
                   rfront = 0.d0
                   ctr = 0
                   do l = lmodemin,lmodemax
                      do m = 0,l
                         ctr = ctr + 1
                         
                         ! set local coordinates
                         if (r.LE.dx(1)) then
                            theta = 0.d0
                            phi = 0.d0
                         else
                            xtemp = cos(2.d0*Pi*alphalm(l,m))*x + sin(2.d0*Pi*alphalm(l,m))*y
                            yloc = -sin(2.d0*Pi*alphalm(l,m))*x + cos(2.d0*Pi*alphalm(l,m))*y
                            
                            xloc = cos(2.d0*Pi*gammalm(l,m))*xtemp + sin(2.d0*Pi*gammalm(l,m))*z
                            zloc = -sin(2.d0*Pi*gammalm(l,m))*xtemp + cos(2.d0*Pi*gammalm(l,m))*z
                            
                            theta = ATAN2(SQRT(xloc**2+yloc**2)/r,zloc/r)-Pi
                            phi = ATAN2(yloc/r,xloc/r)
                         endif
                         
                         rfront = rfront + betalm(l,m) * Ylm(l,m,phi,COS(theta))
                         
                      end do
                   end do

                   if (ctr.gt.0) rfront = rfront/ctr

                else
                   rfront = 0.d0
                end if
                
                rfront = rfire - pertmag*rfront*rfire - r + 3.011d0  
                ! 3.011d0 is roughly the sufrace of fire for pmf.

             else if (prob_type .eq. 2) then
                rfront = 0.d0
             else
                call bl_error("Unknown prob_type")
             end if

             call pmf(rfront,rfront,pmf_vals,n)

             if (n.ne.nspecies+3) then
                write(6,*)"n,nspecies",n,nspecies
                call bl_error('INITDATA: n .ne. nspecies+3')
             endif

             if (prob_type .eq. 1) then
                Tt = pmf_vals(1)
                u1t = 0.d0 ! pmf_vals(2) * x/r
                u2t = 0.d0 ! pmf_vals(2) * y/r
                u3t = 0.d0 ! pmf_vals(2) * z/r
             else if (prob_type .eq. 2) then
                Tt = Tinit
                u1t = uinit
                u2t = vinit
                u3t = winit
             else
                call bl_error("Unknown prob_type")                
             end if

             do n = 1,nspecies
                Xt(n) = pmf_vals(3+n)
             end do

             CALL CKXTY (Xt, IWRK, RWRK, Yt)
             CALL CKRHOY(patmos,Tt,Yt,IWRK,RWRK,rhot)
             call CKUBMS(Tt,Yt,IWRK,RWRK,et)
          
             cons(i,j,k,irho) = rhot
             cons(i,j,k,imx)  = rhot*u1t
             cons(i,j,k,imy)  = rhot*u2t
             cons(i,j,k,imz)  = rhot*u3t
             cons(i,j,k,iene) = rhot*(et + 0.5d0*(u1t**2 + u2t**2 + u3t**2))

             do n=1,nspecies
                cons(i,j,k,iry1-1+n) = Yt(n)*rhot
             end do

          enddo
       enddo
    enddo
    !$omp end parallel do

  end subroutine init_data_3d


  double precision function Ylm(l,m,phi,costheta)
    double precision phi, costheta, facN, facD
    integer m, l, n
    facN = l - m
    if (facN .EQ. 0) facN=1
    do n = l-m-1,2,-1
       facN = facN*n
    end do
    facD = l+m
    if (facD.EQ.0) facD=1
    do n = l+m-1,2,-1
       facD = facD*n
    end do
    Ylm = SQRT((2*l+1)*facN/(4*Pi*facD))*plgndr(l,m,costheta)*COS(m*phi)
  end function Ylm
 
  double precision FUNCTION plgndr(l,m,x)
    INTEGER l,m
    double precision x
    INTEGER i,ll
    double precision fact,pll,pmm,pmmp1,somx2
    if(m.lt.0.or.m.gt.l.or.abs(x).gt.1.) then
       call bl_error('bad arguments in plgndr')
    endif
    pmm=1.d0
    if(m.gt.0) then
       somx2=sqrt((1.d0-x)*(1.d0+x))
       fact=1.d0
       do i=1,m
          pmm=-pmm*fact*somx2
          fact=fact+2.d0
       end do
    endif
    if(l.eq.m) then
       plgndr=pmm
    else
       pmmp1=x*(2*m+1)*pmm
       if(l.eq.m+1) then
          plgndr=pmmp1
       else
          do ll=m+2,l
             pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
             pmm=pmmp1
             pmmp1=pll
          end do
          plgndr=pll
       endif
    endif
  end FUNCTION plgndr
  
end module init_data_module
