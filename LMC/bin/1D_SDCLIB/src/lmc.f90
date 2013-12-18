module lmc
  use probin
  implicit none

  integer, parameter :: maxspec = 53
  integer, parameter :: maxreac = 325
  integer, parameter :: maxspnml = 16
  integer, parameter :: unlim = 0
  integer, parameter :: on_lo = 0, on_hi = 1

  double precision, save :: Pcgs, P1ATM, invmwt(maxspec), mwt(maxspec)
  double precision, save :: T_bc(0:1), rho_bc(0:1), Y_bc(maxspec,0:1), h_bc(0:1), u_bc(0:1)
  double precision, save :: vel_TYP, hmix_TYP, c_0(0:maxspec), c_1(0:maxspec)

  integer, save :: nscal, nspec, nelt, nreac, nfit
  integer, save :: Density, Temp, RhoRT, RhoH, firstspec, lastspec
  integer, save :: iH2, iO2, iN2, iCH4, iCH3OCH3, iCO2, iH2O

  integer, parameter :: eg_nodes = 1
  integer, parameter :: egi = maxspec
  integer, parameter :: egr = 23 + 14*maxspec + 32*maxspec**2 + 13*eg_nodes &
       + 30*eg_nodes*maxspec + 5*eg_nodes*maxspec**2
  integer          :: egiwrk(egi)
  double precision :: egrwrk(egr)

  integer, parameter :: eg_ITLS = 1, eg_IFLAG = 3
contains

  subroutine init_chem()

    double precision    :: ru, ruc, rwrk, rdummy
    integer             :: iwrk, n, idummy, kname(maxspnml*maxspec), offset, len
    character(maxspnml) :: specNames(maxspec)

    ! this sets nspec
    call ckindx(idummy,rdummy,nelt,nspec,nreac,nfit)

    call EGINICD(eg_nodes, 6, eg_IFLAG, eg_ITLS,  EGRWRK, egr, EGIWRK, egi)

    ! other useful things
    call ckrp(iwrk, rwrk, ru, ruc, p1atm)
    call ckwt(iwrk, rwrk, mwt)

    Pcgs = P1ATM

    do n=1,Nspec
       invmwt(n) = 1.d0 / mwt(n)
    end do

    call cksyms(kname,maxspnml)
    do n=1,Nspec
       offset = (n-1)*maxspnml+1
       call convStr(kname(offset),maxspnml,specNames(n),len)
       if (specNames(n).eq.'H2')  iH2=n
       if (specNames(n).eq.'O2')  iO2=n
       if (specNames(n).eq.'N2')  iN2=n
       if (specNames(n).eq.'CH4') iCH4=n
       if (specNames(n).eq.'CH3OCH3') iCH3OCH3=n
       if (specNames(n).eq.'CO2') iCO2=n
       if (specNames(n).eq.'H2O') iH2O=n
    enddo
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
    end if

    call ckhbms(t_bc(0),y_bc(1,0),iwrk,rwrk,h_bc(0))
    call ckrhoy(pcgs,t_bc(0),y_bc(1,0),iwrk,rwrk,rho_bc(0))
  end subroutine init_prob

  subroutine init_data(vel,scal,dx,lo,hi,bc)
    integer,          intent(in   ) :: lo, hi, bc(2)
    double precision, intent(  out) :: vel(lo-2:hi+2)
    double precision, intent(  out) :: scal(lo-2:hi+2,nscal)
    double precision, intent(in   ) :: dx

    double precision :: x, rho, Y(Nspec), T, h
    double precision :: xPMFlo, xPMFhi
    double precision :: valsPMF(Nspec+3), RWRK, sum, sigma, cent
    integer          :: i, n, l, nPMF, IWRK

    do i=lo,hi
       x = (dble(i)+0.5d0)*dx

       if (probtype.eq.1) then
          xPMFlo = x - flame_offset - 0.5*dx
          xPMFhi = x - flame_offset + 0.5*dx
          call pmf(xPMFlo,xPMFhi,valsPMF,nPMF)
       else if (probtype.eq.2) then
          sigma = 0.0451d0
          cent = 0.d0
          valsPMF(1) = 5.d2 + 1.d3*dexp(-(x-cent)**2/(2.d0*sigma**2))
          valsPMF(2) = 0.d0
          do n=1,Nspec
             valsPMF(3+n) = 0.d0
             if (n.eq.iCH3OCH3) valsPMF(3+n) = .0374d0
             if (n.eq.iCO2)     valsPMF(3+n) = .0522d0
             if (n.eq.iO2)      valsPMF(3+n) = .1128d0
             if (n.eq.iN2)      valsPMF(3+n) = .7192d0
             if (n.eq.iH2O)     valsPMF(3+n) = .0784d0
          end do
       end if

       call CKXTY(valsPMF(4),IWRK,RWRK,Y)
       T = valsPMF(1)

       if (iN2.gt.0  .and.  iN2.le.Nspec) then
          sum = 0.d0
          do n=1,Nspec
             sum = sum + Y(n)
          end do
          Y(iN2) = Y(iN2) - 1.d0 + sum
       end if
       call CKHBMS(T,Y,IWRK,RWRK,h)
       call CKRHOY(Pcgs,T,Y,IWRK,RWRK,rho)
       do n=1,Nspec
          scal(i,FirstSpec+n-1) = rho * Y(n)
       end do
       scal(i,Density) = rho
       scal(i,Temp) = T
       scal(i,RhoH) = rho * h
       vel(i) = valsPMF(2)
       if (l .eq. 0 .and. i .eq. lo) then
          hmix_TYP = ABS(h)
       else
          hmix_TYP = MAX(hmix_TYP,ABS(h))
       end if

    end do

    c_0 = 0.d0
    c_1 = 0.d0

    !     fill coarse level ghost cells
    call set_bc_s(scal,lo,hi,bc)
    call set_bc_v(vel,lo,hi,bc)

    vel_TYP = ABS(vel(0))
    do i=lo,hi
       vel_TYP = MAX(vel_TYP,ABS(vel(i)))
    end do

  end subroutine init_data

  subroutine set_bc_s(scal,lo,hi,bc)
    integer,          intent(in   ) :: lo,hi,bc(2)
    double precision, intent(inout) :: scal(lo-2:hi+2,nscal)

    integer          :: n, is, HorL
    double precision :: u, rho, Y(Nspec), T, hmix

    if (probtype.eq.1) then
       if (bc(1) .eq. 1) then
          ! lo:  Dirichlet values for rho, Y, T, h
          HorL = on_lo
          call bcfunction(HorL, u, rho, Y, T, hmix)

          scal(lo-2:lo-1,Density) = rho
          scal(lo-2:lo-1,Temp) = T
          do n=1,Nspec
             is = FirstSpec + n - 1
             scal(lo-2:lo-1,is) = Y(n) * rho
          enddo
          scal(lo-2:lo-1,RhoH) = hmix * rho
       end if
    else if (probtype.eq.2) then
       ! lo:  Neumann for all
       scal(lo-2:lo-1,Density) = scal(lo,Density)
       scal(lo-2:lo-1,Temp) = scal(lo,Temp)
       do n=1,Nspec
          is = FirstSpec + n - 1
          scal(lo-2:lo-1,is) = scal(lo,is)
       enddo
       scal(lo-2:lo-1,RhoH) = scal(lo,RhoH)
    else
       print *,'Unknown probtype',probtype
       stop
    endif

    if (bc(2) .eq. 2) then
       ! hi:  Neumann for all
       scal(hi+1:hi+2,Density) = scal(hi,Density)
       scal(hi+1:hi+2,Temp) = scal(hi,Temp)
       do n=1,Nspec
          is = FirstSpec + n - 1
          scal(hi+1:hi+2,is) = scal(hi,is)
       enddo
       scal(hi+1:hi+2,RhoH) = scal(hi,RhoH)
    end if

  end subroutine set_bc_s

  subroutine set_bc_v(vel,lo,hi,bc)
    integer,          intent(in   ) :: lo,hi,bc(2)
    double precision, intent(inout) :: vel(lo-2:hi+2)

    integer          :: HorL
    double precision :: u, rho, Y(Nspec), T, hmix

    if (bc(1) .eq. 1) then
       ! lo:  Dirichlet
       HorL = on_lo
       call bcfunction(HorL, u, rho, Y, T, hmix)
       vel(lo-2:lo-1) = u
    end if

    if (bc(2) .eq. 2) then
       ! hi:  Neumann for all
       vel(hi+1:hi+2) = vel(hi)
    end if

  end subroutine set_bc_v

  subroutine bcfunction(LorR, u, rho, Y, T, h)
    integer,          intent(in   ) :: LorR
    double precision, intent(  out) :: u, rho, T, h, Y(:)

    integer :: n

    do n = 1,Nspec
       Y(n) = Y_bc(n,LorR)
    end do
    T = T_bc(LorR)
    rho = rho_bc(LorR)
    h = h_bc(LorR)
    u = u_bc(LorR)

  end subroutine bcfunction

  subroutine convStr(codedString,maxLen,string,strLen)
    implicit none
    integer codedString(*), maxLen, k, strLen
    character*(*) string
    do k=1,maxLen
       string(k:k) = ' '
       if (codedString(k) .ne. ICHAR(' ')) then
          string(k:k) = CHAR(codedString(k))
          strLen = k
       endif
    enddo
  end subroutine convStr

end module lmc
