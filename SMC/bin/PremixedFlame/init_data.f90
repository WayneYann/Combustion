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
    use chemistry_module, only : nspecies, patm, get_species_index
    use probin_module

    integer,          intent(in   ) :: lo(1),hi(1),ng
    double precision, intent(in   ) :: dx(1),phlo(1),phhi(1)
    double precision, intent(inout) :: cons(-ng+lo(1):hi(1)+ng,ncons)

    integer          :: i,n
    double precision :: x

    double precision Xt(nspecies), Yt(nspecies), Yti(nspecies)
    double precision rhot,Tt,et,Pt, eta, eta01, eta00
    integer :: iwrk
    double precision :: rwrk

    double precision, parameter :: flrt = 0.07d0
    double precision, parameter :: pres = 1.d0
    !
    double precision, parameter :: reac_H2  = 0.296d0
    double precision, parameter :: reac_O2  = 0.148d0
    double precision, parameter :: reac_N2  = 0.556d0
    !
    double precision, parameter :: prod_H2O = 0.296d0
    double precision, parameter :: prod_O2  = 0.000d0
    double precision, parameter :: prod_N2  = 0.556d0
    !
    double precision, parameter :: intm_HO2  = 0.0001d0
    double precision, parameter :: intm_O    = 0.0001d0
    double precision, parameter :: intm_H2O2 = 0.0001d0
    double precision, parameter :: intm_H    = 0.0100d0
    double precision, parameter :: intm_OH   = 0.0100d0
    !
    double precision, parameter :: temp_r = 298.d0
    double precision, parameter :: temp_p = 2300.d0
    !
    double precision, parameter :: flame_width = 0.01d0
    double precision, parameter :: flame_pos = 0.d0

    integer, save :: iH2O=-100, iHO2, iO, iH2O2, iH, iOH
    double precision, dimension(:), allocatable, save :: Xr, Xp, Xi

    if (iH2O .eq. -100) then
       iH2O = get_species_index("H2O")
       iHO2 = get_species_index("HO2")
       iO = get_species_index("O")
       iH2O2 = get_species_index("H2O2")
       iH = get_species_index("H")
       iOH = get_species_index("OH")

       allocate(Xr(nspecies))
       allocate(Xp(nspecies))
       allocate(Xi(nspecies))

       Xr = 0.d0
       Xr(iH2) = reac_H2
       Xr(iO2) = reac_O2
       Xr(iN2) = reac_N2
       
       Xp = 0.d0
       Xp(iH2O) = prod_H2O
       Xp(iO2)  = prod_O2
       Xp(iN2)  = prod_N2

       Xi = 0.d0
       Xi(iHO2)  = intm_HO2
       Xi(iO)    = intm_O
       Xi(iH2O2) = intm_H2O2
       Xi(iH)    = intm_H
       Xi(iOH)   = intm_OH
    end if

    cons = 0.0d0

    do i=lo(1),hi(1)
       x = phlo(1) + dx(1)*(i+0.5d0) - flame_pos

       Pt = pres*patm

       if (x .lt. -10.d0*flame_width) then
          Tt = temp_r
          Xt = Xr
          Xt = Xt/sum(Xt)
          CALL CKXTY (Xt, IWRK, RWRK, Yt)
       else if (x .gt. 10.d0*flame_width) then
          Tt = temp_p
          Xt = Xp
          Xt = Xt/sum(Xt)
          CALL CKXTY (Xt, IWRK, RWRK, Yt)
       else
          eta = tanh(x/flame_width)
          eta01 = (eta+1.d0)*0.5d0
          Tt = (1.d0-eta01)*temp_r + eta01*temp_p
          Xt = (1.d0-eta01)*Xr + eta01*Xp
          Xt = Xt/sum(Xt)
          CALL CKXTY (Xt, IWRK, RWRK, Yti)
          ! add flame
          eta00 = 1.d0 - abs(eta)
          Xt = Xt + eta00*Xi
          Xt = Xt/sum(Xt)
          CALL CKXTY (Xt, IWRK, RWRK, Yt)
          Yt = Yt*(1.d0-Yti(iN2))/(1.d0-Yt(iN2))
          Yt(iN2) = Yti(iN2)          
       end if       

       CALL CKRHOY(Pt,Tt,Yt,IWRK,RWRK,rhot)
       call CKUBMS(Tt,Yt,IWRK,RWRK,et)

       cons(i,irho) = rhot
       cons(i,imx)  = flrt
       cons(i,iene) = rhot*et + 0.5d0*cons(i,imx)**2/cons(i,irho)
       do n=1,nspecies
          cons(i,iry1-1+n) = Yt(n)*rhot
       end do
    end do

  end subroutine init_data_1d

end module init_data_module
