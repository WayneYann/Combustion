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
    use probin_module,    only : prob_type, pamb, T_in, vn_in, T_co, vn_co, &
         splitx, xfrontw, Tfrontw, vt_in, vt_co, blobr, blobx, bloby, blobT

    integer,          intent(in   ) :: lo(2),hi(2),ng
    double precision, intent(in   ) :: dx(2),phlo(2),phhi(2)
    double precision, intent(inout) :: cons(-ng+lo(1):hi(1)+ng,-ng+lo(2):hi(2)+ng,ncons)

    integer          :: i,j,n
    double precision :: x, y

    double precision Yt(nspecies)
    double precision rhot,u1t,u2t,Tt,et, sigma, eta, eta1, r
    integer :: iwrk
    double precision :: rwrk

    if (prob_type .ne. 0 .and. prob_type .ne. 1) then
       call bl_error('unsupported prob_type')
    end if

    sigma = 2.5d0*xfrontw*splitx

    do j=lo(2),hi(2)

       y = phlo(2) + dx(2)*(j + 0.5d0)

       do i=lo(1),hi(1)

          x = phlo(1) + dx(1)*(i + 0.5d0)

          if (prob_type .eq. 0) then
             eta = 0.5d0 * (tanh((x + splitx)/sigma)   &
                  &       - tanh((x - splitx)/sigma))
          else if (prob_type .eq. 1) then
             eta = 0.5d0 * (tanh((x + splitx)/Tfrontw)  &
                  &       - tanh((x - splitx)/Tfrontw))
             eta1 = 0.5d0 * (tanh((x + blobr)/xfrontw)  &
                  &        - tanh((x - blobr)/xfrontw))
          end if
       
          do n=1,nspecies
             Yt(n) = eta*fuel_Y(n) + (1.d0-eta)*air_Y(n)
          end do
          Tt  = eta  *  T_in + (1.d0-eta ) * T_co
          if (prob_type .eq. 0) then 
             u1t = eta  * vt_in + (1.d0-eta ) * vt_co
             u2t = eta  * vn_in + (1.d0-eta ) * vn_co
          else 
             u1t = eta1 * vt_in + (1.d0-eta1) * vt_co
             u2t = eta1 * vn_in + (1.d0-eta1) * vn_co
          end if
       
          if (blobr .gt. 0.d0) then
             eta = 0.5d0*(1.d0 - TANH(-2.d0*(y-bloby)/Tfrontw))
             Tt  = eta * T_co + (1.d0-eta) * Tt
             do n=1,nspecies
                Yt(n) = eta*air_Y(n) + (1.d0-eta)*Yt(n)
             end do

             ! Superimpose blob of hot air
             r = SQRT((x-blobx)**2 + (y-bloby)**2)
             eta = 0.5d0*(1.d0 - TANH(2.d0*(r-blobr)/Tfrontw))
             do n=1,nspecies
                Yt(n) = eta*air_Y(n) + (1.d0-eta)*Yt(n)
             enddo
             Tt  = eta * blobT + (1.d0-eta) * Tt
          end if

          CALL CKRHOY(pamb,Tt,Yt,IWRK,RWRK,rhot)
          call CKUBMS(Tt,Yt,IWRK,RWRK,et)
                 
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
