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
    use probin_module,    only : Minfty, Rvortex, Cvortex

    integer,          intent(in   ) :: lo(3),hi(3),ng
    double precision, intent(in   ) :: dx(3),phlo(3),phhi(3)
    double precision, intent(inout) :: cons(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,-ng+lo(3):hi(3)+ng,ncons)

    integer          :: i,j,k,n
    double precision :: x, y, z

    double precision pmf_vals(nspecies+3)
    double precision Xt(nspecies), Yt(nspecies)
    double precision rhot,u1t,u2t,u3t,Tt,et,Pt
    integer :: iwrk
    double precision :: rwrk, Lx, exptmp, Wbar, Rc, Cc, cs, Cv, Cp, gamma, uinfty

    double precision, parameter :: patmos = 1.01325d6, Ru = 8.31451d7

    !$omp parallel do &
    !$omp private(i,j,k,n,x,y,z,pmf_vals) &
    !$omp private(Xt,Yt,rhot,u1t,u2t,u3t,Tt,et,iwrk,rwrk) &
    !$omp private(Lx,exptmp, Wbar, Rc, Cc, cs, Cv, Cp, gamma, uinfty)
    do k=lo(3),hi(3)
       z = phlo(3) + dx(3)*(k + 0.5d0)
       do j=lo(2),hi(2)
          y = phlo(2) + dx(2)*(j + 0.5d0)
          do i=lo(1),hi(1)
             x = phlo(1) + dx(1)*(i + 0.5d0)

             Lx = phhi(1) - phlo(1)

             call pmf(-1.d0,-1.d0,pmf_vals,n)

             if (n.ne.nspecies+3) then
                write(6,*)"n,nspecies",n,nspecies
                call bl_error('INITDATA: n .ne. nspecies+3')
             endif

             Tt = pmf_vals(1)

             do n = 1,nspecies
                Xt(n) = pmf_vals(3+n)
             end do
             CALL CKXTY (Xt, IWRK, RWRK, Yt)
             CALL CKRHOY(patmos,Tt,Yt,IWRK,RWRK,rhot)
             ! we now have rho

             call CKCVBL(Tt, Xt, iwrk, rwrk, Cv)
             Cp = Cv + Ru
             gamma = Cp / Cv
             cs = sqrt(gamma*patmos/rhot)

             uinfty = Minfty*cs
             Rc = Rvortex*Lx
             Cc = Cvortex*cs*Lx

             exptmp = exp(-(x**2+y**2)/(2.d0*Rc**2))
             Pt = patmos - rhot*(Cc/Rc)**2*exptmp

             u1t = uinfty - Cc*exptmp*y/Rc**2
             u2t = Cc*exptmp*x/Rc**2
             u3t = 0.d0

             call CKMMWX(Xt, iwrk, rwrk, Wbar)
             Tt = Pt*Wbar / (rhot*Ru)

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
