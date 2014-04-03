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

    integer                   :: lo(1), hi(1), ng, i, dm
    double precision, pointer :: dp(:,:,:,:)

    ng = data%ng
    dm = data%dim
    if (dm.ne.1) then
       call bl_error("This problem is 1D only.")
    end if

    do i=1,nfabs(data)
       dp => dataptr(data,i)
       lo = lwb(get_box(data,i))
       hi = upb(get_box(data,i))

       call init_data_1d(lo,hi,ng,dx,dp,plo,phi)
    end do

  end subroutine init_data

  subroutine init_data_1d(lo,hi,ng,dx,cons,phlo,phhi)

    use variables_module, only : irho, imx,iene,iry1,ncons, iH2, iO2, iN2
    use chemistry_module, only : nspecies, patm
    use probin_module

    integer,          intent(in   ) :: lo(1),hi(1),ng
    double precision, intent(in   ) :: dx(1),phlo(1),phhi(1)
    double precision, intent(inout) :: cons(-ng+lo(1):hi(1)+ng,ncons)

    integer          :: i,n
    double precision :: x, r

    double precision pmf_vals(nspecies+3)
    double precision Xt(nspecies), Yt(nspecies)
    double precision rhot,u1t,Tt,et,Pt
    integer :: iwrk
    double precision :: rwrk

    cons = 0.0d0

    !$omp parallel do &
    !$omp private(i,n,x,r,pmf_vals,Pt,Xt,Yt,rhot,u1t,Tt,et,iwrk,rwrk)
    do i=lo(1),hi(1)
       x = phlo(1) + dx(1)*(i+0.5d0)
       r = x + 3.0d0
       ! the magic number is roughly the sufrace of fire for pmf.

       call pmf(r,r,pmf_vals,n)

       if (n .ne. nspecies+3) then
          write(6,*)"n,nspecies",n,nspecies
          call bl_error('INITDATA: n .ne. nspecies+3')
       end if

       Pt = patm
       Tt = pmf_vals(1)
       u1t = pmf_vals(2)

       do n = 1,nspecies
          Xt(n) = pmf_vals(3+n)
       end do

       CALL CKXTY (Xt, IWRK, RWRK, Yt)
       CALL CKRHOY(Pt,Tt,Yt,IWRK,RWRK,rhot)
       call CKUBMS(Tt,Yt,IWRK,RWRK,et)

       cons(i,irho) = rhot
       cons(i,imx)  = rhot*u1t
       cons(i,iene) = rhot*(et + 0.5d0*u1t**2)
       do n=1,nspecies
          cons(i,iry1-1+n) = Yt(n)*rhot
       end do
    end do
    !$omp end parallel do

  end subroutine init_data_1d

end module init_data_module
