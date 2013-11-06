module lmc
  use probin
  implicit none

  integer, parameter :: maxspec = 53
  integer, parameter :: maxreac = 325
  integer, parameter :: unlim = 0

  double precision, save :: pcgs, p1atm, invmwt(maxspec), mwt(maxspec)
  double precision, save :: T_bc(0:1), rho_bc(0:1), Y_bc(maxspec,0:1), h_bc(0:1), u_bc(0:1)

  integer, save :: nscal, nspec, nelt, nreac, nfit
  integer, save :: density, temp, rhort, rhoh, firstspec, lastspec

contains

  subroutine init_chem()

    double precision :: ru, ruc, rwrk, rdummy
    integer          :: iwrk, n, idummy

    ! this sets nspec
    call ckindx(idummy,rdummy,nelt,nspec,nreac,nfit)

    ! lout = 6
    ! call EGINICD(eg_nodes, lout, eg_IFLAG, eg_ITLS,  EGRWRK, egr, EGIWRK, egi)

    ! other useful things
    call ckrp(iwrk, rwrk, ru, ruc, p1atm)
    call ckwt(iwrk, rwrk, mwt)

    do n=1,Nspec
       invmwt(n) = 1.d0 / mwt(n)
    end do

    ! call cksyms(kname,maxspnml)
    ! do n=1,Nspec
    !    offset = (n-1)*maxspnml+1
    !    call convStr(kname(offset),maxspnml,specNames(n),len)
    !    if (specNames(n).eq.'H2')  iH2=n
    !    if (specNames(n).eq.'O2')  iO2=n
    !    if (specNames(n).eq.'N2')  iN2=n
    !    if (specNames(n).eq.'CH4') iCH4=n
    !    if (specNames(n).eq.'CH3OCH3') iCH3OCH3=n
    !    if (specNames(n).eq.'CO2') iCO2=n
    !    if (specNames(n).eq.'H2O') iH2O=n
    ! enddo
  end subroutine init_chem

  subroutine init_prob(domnlo)
    double precision, intent(in) :: domnlo

    integer          :: npmf, iwrk
    double precision :: xpmf, valspmf(maxspec+3), rwrk

    Density   = 1
    Temp      = 2
    RhoH      = 3
    RhoRT     = 4
    FirstSpec = 5
    LastSpec  = FirstSpec + (nspec-1)
    nscal     = LastSpec

    ! set boundary data
    xpmf = domnlo - flame_offset
    call pmf(xpmf,xpmf,valspmf,npmf)
    call ckxty(valspmf(4),iwrk,rwrk,y_bc(1,0))
    t_bc(0) = valspmf(1)
    if (abs(v_in) .lt. 1.d20) then
       u_bc(0) = v_in
    else
       u_bc(0) = valspmf(2)
    endif

    call ckhbms(t_bc(0),y_bc(1,0),iwrk,rwrk,h_bc(0))
    call ckrhoy(pcgs,t_bc(0),y_bc(1,0),iwrk,rwrk,rho_bc(0))
  end subroutine init_prob

end module lmc
