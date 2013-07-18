module init_data_module

  use multifab_module
  use bl_constants_module, only : Pi=>M_PI

  use DME_jet_module

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

    if (.not. initialized) then
       call init_DME_jet()
    end if

    ng = data%ng

    do i=1,nfabs(data)
       dp => dataptr(data,i)
       lo = lwb(get_box(data,i))
       hi = upb(get_box(data,i))

       select case(data%dim)
       case (2)
          call init_data_2d(lo,hi,ng,dx,dp,plo,phi)
       case (3)
          call bl_error('We only support 2-D for this test')
       end select
    end do

  end subroutine init_data

  subroutine init_data_2d(lo,hi,ng,dx,cons,phlo,phhi)

    use variables_module, only : irho, imx,imy,iene,iry1,ncons
    use chemistry_module, only : nspecies
    use probin_module,    only : pamb, T_in, vn_in, T_co, vn_co, splity, yfrontw

    integer,          intent(in   ) :: lo(2),hi(2),ng
    double precision, intent(in   ) :: dx(2),phlo(2),phhi(2)
    double precision, intent(inout) :: cons(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,ncons)

    integer          :: i,j,n
    double precision :: y

    double precision Yt(nspecies)
    double precision rhot,u1t,u2t,Tt,et, sigma, eta
    integer :: iwrk
    double precision :: rwrk

    sigma = 2.5d0*yfrontw*splity

    do j=lo(2),hi(2)
       y = phlo(2) + dx(2)*(j + 0.5d0)

       eta = 0.5d0 * (tanh((y + splity)/sigma)   &
            &       - tanh((y - splity)/sigma))
       
       do n=1,nspecies
          Yt(n) = eta*fuel_Y(n) + (1.d0-eta)*air_Y(n)
       end do
       Tt  = eta * T_in + (1.d0-eta) * T_co
       u1t = eta *vn_in + (1.d0-eta) *vn_co
       u2t = 0.d0
       
       CALL CKRHOY(pamb,Tt,Yt,IWRK,RWRK,rhot)
       call CKUBMS(Tt,Yt,IWRK,RWRK,et)
       
       do i=lo(1),hi(1)
          
          cons(i,j,irho) = rhot
          cons(i,j,imx)  = rhot*u1t
          cons(i,j,imy)  = rhot*u2t
          cons(i,j,iene) = rhot*(et + 0.5d0*(u1t**2 + u2t**2))
          
          do n=1,nspecies
             cons(i,j,iry1-1+n) = Yt(n)*rhot
          end do
          
       enddo
    enddo

  end subroutine init_data_2d
  
end module init_data_module
