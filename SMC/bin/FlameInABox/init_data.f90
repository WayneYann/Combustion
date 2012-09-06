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

  subroutine init_data_3d(lo,hi,ng,dx,cons,phlo,phhi)

    use variables_module, only : irho, imx,imy,imz,iene,iry1,ncons
    use chemistry_module, only : nspecies
    use probin_module,    only : pertmag, xfire

    integer,          intent(in   ) :: lo(3),hi(3),ng
    double precision, intent(in   ) :: dx(3),phlo(3),phhi(3)
    double precision, intent(inout) :: cons(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,ncons)

    integer          :: i,j,k,n
    double precision :: x, y, z

    double precision pmf_vals(nspecies+3)
    double precision Xt(nspecies), Yt(nspecies)
    double precision rhot,u1t,u2t,u3t,Tt,et
    integer :: iwrk
    double precision :: rwrk
    double precision :: pert, Ly, Lz, x1, x2

    double precision, parameter :: patmos = 1.01325d6

    !$omp parallel do &
    !$omp private(i,j,k,n,x,y,z,pmf_vals) &
    !$omp private(Xt,Yt,rhot,u1t,u2t,u3t,Tt,et,iwrk,rwrk) &
    !$omp private(pert,Ly,Lz,x1,x2)
    do k=lo(3),hi(3)
       z = phlo(3) + dx(3)*(k + 0.5d0)
       do j=lo(2),hi(2)
          y = phlo(2) + dx(2)*(j + 0.5d0)
          do i=lo(1),hi(1)
             x = phlo(1) + dx(1)*(i + 0.5d0)

             pert = 0.d0
             if (pertmag .gt. 0.d0) then
                Ly = phhi(2) - phlo(2)
                Lz = phhi(3) - phlo(3)

                pert = pertmag*(1.000d0 * sin(2*Pi*4*y/Ly)             * sin(2*Pi*5*z/Lz) &
                     &        + 1.023d0 * sin(2*Pi*2*(y-.004598)/Ly)   * sin(2*Pi*4*(z-.0053765)/Lz) &
                     &        + 0.945d0 * sin(2*Pi*3*(y-.00712435)/Ly) * sin(2*Pi*3*(z-.02137)/Lz) &
                     &        + 1.017d0 * sin(2*Pi*5*(y-.0033)/Ly)     * sin(2*Pi*6*(z-.018)/Lz) & 
                     &        + 0.982d0 * sin(2*Pi*5*(y-.014234)/Ly) )
             endif

             x1 = (x + 3.011d0 + xfire - 0.5d0*dx(1) + pert)
             x2 = (x + 3.011d0 + xfire + 0.5d0*dx(1) + pert)

             call pmf(x1,x2,pmf_vals,n)

             if (n.ne.nspecies+3) then
                write(6,*)"n,nspecies",n,nspecies
                call bl_error('INITDATA: n .ne. nspecies+3')
             endif

             Tt = pmf_vals(1)
             u1t = pmf_vals(2) 
             u2t = 0.d0 
             u3t = 0.d0 

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
  
end module init_data_module
