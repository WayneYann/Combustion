module inflow_module

  implicit none
  
  double precision, allocatable, save :: inflow_state(:)

contains

  subroutine init_inflow()

    use probdata_module
    use meth_params_module, only : NVAR, URHO, UMX, UEDEN, UTEMP, UFS, NSPEC

    double precision :: Pt, rhot, Tt, et, Xt(nspec), Yt(nspec), rwrk
    integer :: n,iwrk

    if (.not. allocated(inflow_state)) then
       allocate(inflow_state(NVAR))
    end if

    Pt = pres
    Tt = T_reac
    Xt = X_reac
    Xt = Xt/sum(Xt)
    CALL CKXTY (Xt, IWRK, RWRK, Yt)

    CALL CKRHOY(Pt,Tt,Yt,IWRK,RWRK,rhot)
    call CKUBMS(Tt,Yt,IWRK,RWRK,et)

    inflow_state(URHO) = rhot
    inflow_state(UMX)  = massFlux
    inflow_state(UEDEN) = rhot*et + 0.5d0*inflow_state(UMX)**2/inflow_state(URHO)
    inflow_state(UTEMP) = Tt
    do n=1,nspec
       inflow_state(UFS-1+n) = Yt(n)*rhot
    end do
    
  end subroutine init_inflow
  
end module inflow_module
