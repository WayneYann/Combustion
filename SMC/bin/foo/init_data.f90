module init_data_module

  use multifab_module

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

    do i=1,nboxes(data)
       if ( multifab_remote(data,i) ) cycle

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

  subroutine init_data_3d(lo,hi,ng,dx,cons,plo,phi)

    use variables, only : irho, imx,imy,imz,iene,iry1,ncons
    use chemistry_module, only : nspecies

    integer,          intent(in   ) :: lo(3),hi(3),ng
    double precision, intent(in   ) :: dx(3),plo(3),phi(3)
    double precision, intent(inout) :: cons(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,ncons)

    integer          :: i,j,k,n
    double precision :: x, y, z, r, rmin, rmax

    double precision pmf_vals(nspecies+3)
    double precision Xt(nspecies), Yt(nspecies)
    double precision rhot,u1t,u2t,u3t,Tt,et
    integer :: iwrk
    double precision :: rwrk
    double precision :: pert, Lx, Ly

    double precision, parameter :: patmos = 1.01325d6
    double precision, parameter :: pertmag = 1.d-2
    double precision, parameter :: PI = 4.d0*atan(1.d0)

    do k=lo(3),hi(3)
       z = plo(3) + dx(3)*(k + 0.5d0)
       do j=lo(2),hi(2)
          y = plo(2) + dx(2)*(j + 0.5d0)
          do i=lo(1),hi(1)
             x = plo(1) + dx(1)*(i + 0.5d0)

             Lx = phi(1) - plo(1)
             Ly = phi(2) - plo(2)
             pert = pertmag*(1.000 * sin(2*Pi*4*x/Lx)             * sin(2*Pi*5*y/Ly) &
                  &        + 1.023 * sin(2*Pi*2*(x-.004598)/Lx)   * sin(2*Pi*4*(y-.0053765)/Ly) &
                  &        + 0.945 * sin(2*Pi*3*(x-.00712435)/Lx) * sin(2*Pi*3*(y-.02137)/Ly) &
                  &        + 1.017 * sin(2*Pi*5*(x-.0033)/Lx)     * sin(2*Pi*6*(y-.018)/Ly) &
                  &         + .982 * sin(2*Pi*5*(x-.014234)/Lx) )
             
             r = abs(z) + 2.5d0 + pert
             rmin = (r - 0.5d0*dx(3))
             rmax = (r + 0.5d0*dx(3))

             call pmf(rmin,rmax,pmf_vals,n)

             if (n.ne.nspecies+3) then
                write(6,*)"n,nspecies",n,nspecies
                call bl_error('INITDATA: n .ne. nspecies+3')
             endif

             Tt = pmf_vals(1)
             do n = 1,nspecies
                Xt(n) = pmf_vals(3+n)
             end do
             u1t = 0.d0 ! pmf_vals(2) * x/r
             u2t = 0.d0 ! pmf_vals(2) * y/r
             u3t = 0.d0 ! pmf_vals(2) * z/r
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

  end subroutine init_data_3d
  
end module init_data_module
