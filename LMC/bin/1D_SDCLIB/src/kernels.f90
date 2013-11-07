module kernels
  use probin, only: nx
  use lmc
  implicit none
  double precision, parameter :: eps = 1.d-6
contains

  !
  ! Calculate cell-centered grad pi from nodal pi
  !
  subroutine calc_grad_pi(gp, press, lo, hi, dx)
    integer,          intent(in   ) :: lo, hi
    double precision, intent(  out) :: gp(lo-1:hi+1)
    double precision, intent(in   ) :: press(lo-1:hi+2), dx

    integer :: i

    do i = lo-1, hi+1
       gp(i) = (press(i+1) - press(i)) / dx
    end do
  end subroutine calc_grad_pi

  !
  ! Calculate edge velocities
  !
  subroutine calc_edge_vel(edgevel, vel, lo, hi, bc)
    integer,          intent(in   ) :: lo, hi, bc(2)
    double precision, intent(  out) :: edgevel(lo:hi+1)
    double precision, intent(in   ) :: vel(lo-2:hi+2)

    integer          :: i
    double precision :: slo, shi, slope(lo-1:hi+1)

    call mkslopes(vel, slope, lo, hi, bc)

    do i = lo, hi+1
       slo = vel(i-1) + 0.5d0 * slope(i-1)
       shi = vel(i  ) - 0.5d0 * slope(i  )
       if ( (slo+shi) .gt. eps) then
          edgevel(i) = slo
       else if ( (slo+shi) .lt. -eps) then
          edgevel(i) = shi
       else ! XXX if ( (abs(slo+shi) .le. eps) .or. (slo .le. 0.d0 .and. shi .ge. 0.d0)) then
          edgevel(i) = 0.d0
       end if
       if (i .eq. lo   .and. bc(1) .eq. 1) edgevel(i) = vel(i-1)
       if (i .eq. hi+1 .and. bc(2) .eq. 2) edgevel(i) = slo
    end do
  end subroutine calc_edge_vel

  !
  ! Project.
  !
  subroutine mac_project(macvel,rho,divu,dx,lo,hi,bc)
    integer,          intent(in   ) :: lo, hi, bc(2)
    double precision, intent(inout) :: macvel(lo:hi+1)
    double precision, intent(in   ) :: rho(lo-2:hi+2), divu(lo-1:hi+1), dx

    integer          :: i
    double precision :: a(nx), b(nx), c(nx), r(nx), u(nx), gam(nx), phi(-1:nx+1)

    do i=lo,hi
       u(i+1) = 0.d0
       r(i+1) = (macvel(i+1)-macvel(i))/dx - divu(i)
       a(i+1) =  (2.d0/dx**2)*(1.d0/(rho(i-1)+rho(i)))
       b(i+1) = -(2.d0/dx**2)*(1.d0/(rho(i-1)+rho(i)) + 1.d0/(rho(i)+rho(i+1)))
       c(i+1) =  (2.d0/dx**2)*(1.d0/(rho(i)+rho(i+1)))

       if (i .eq. lo .and. bc(1) .eq. 1) then
          a(i+1) = 0.d0
          b(i+1) = -(2.d0/dx**2)*(1.d0/(rho(i)+rho(i+1)))
          c(i+1) = -b(i+1)
       end if

       if (i .eq. hi .and. bc(2) .eq. 2) then
          a(i+1) = (1.d0/(3.d0*dx**2))*(1.d0/rho(i)) + (2.d0/dx**2)*(1.d0/(rho(i-1)+rho(i)))
          b(i+1) = -(3.d0/dx**2)*(1.d0/rho(i))       - (2.d0/dx**2)*(1.d0/(rho(i-1)+rho(i)))
          c(i+1) = 0.d0
       end if
    end do

    call tridiag(a,b,c,r,u,gam,hi-lo+1)

    phi(lo:hi) = u(lo+1:hi+1)
    if (bc(1) .eq. 1) phi(lo-1) = phi(lo)
    if (bc(2) .eq. 2) phi(hi+1) = -2.d0*phi(hi) + (1.d0/3.d0)*phi(hi-1)

    do i=lo,hi+1
       macvel(i) = macvel(i) - (2.d0/(rho(i-1)+rho(i)))*(phi(i)-phi(i-1))/dx
    end do
  end subroutine mac_project


  !
  ! Scalar advection.
  !
  subroutine scal_aofs(scal,macvel,aofs,divu,dx,lo,hi,bc)
    integer,          intent(in   ) :: lo, hi, bc(2)
    double precision, intent(in   ) :: dx, scal(lo-2:hi+2,nscal), macvel(lo:hi+1), divu(lo-1:hi+1)
    double precision, intent(  out) :: aofs(lo:hi,nscal)

    double precision :: sedge(lo:hi+1,nscal), slope(lo-1:hi+1), slo, shi, Y(nspec), rwrk, hmix
    integer          :: i, n, iwrk, ispec
    logical          :: compute_comp(nscal)

    compute_comp = .true.
    compute_comp(Density) = .false.
    compute_comp(RhoRT)   = .false. ! XXX: MWE: should this be RhoH?
    compute_comp(Temp)    = .false.

    do n = 1,nscal
       if (.not. compute_comp(n)) cycle

       call mkslopes(scal(:,n),slope,lo,hi,bc)

       do i=lo,hi+1
          slo = scal(i-1,n) + 0.5d0 * slope(i-1)
          shi = scal(i  ,n) - 0.5d0 * slope(i  )

          if ( macvel(i) .gt. eps) then
             sedge(i,n) = slo
          else if ( macvel(i) .lt. -eps) then
             sedge(i,n) = shi
          else if ( abs(macvel(i)) .le. eps) then
             sedge(i,n) = 0.5d0 * (slo + shi)
          endif

          if (i .eq. lo   .and. bc(1) .eq. 1) sedge(i,n) = scal(i-1,n)
          if (i .eq. hi+1 .and. bc(2) .eq. 2) sedge(i,n) = slo
       end do
    end do

    do i=lo,hi+1
       sedge(i,Density) = 0.d0
       ! compute Rho on edges as sum of (Rho Y_i) on edges,
       do n = 1,nspec
          ispec = FirstSpec-1+n
          sedge(i,Density) = sedge(i,Density) + sedge(i,ispec)
       end do
       ! compute rho.hmix as sum of (H_i.Rho.Y_i)
       if (.not. compute_comp(RhoH) ) then
          ! XXX: MWE: it seems like we want to do this no matter
          ! what..., otherwise sedge(:,Temp) is bogus down below...
          do n = 1,nspec
             ispec = FirstSpec-1+n
             Y(n) = sedge(i,ispec) / sedge(i,Density)
          end do
          call CKHBMS(sedge(i,Temp),Y,IWRK,RWRK,Hmix)
          sedge(i,RhoH) = Hmix * sedge(i,Density)
       end if
       ! XXX: MWE: what about sedge(:,RhoRT)?
    end do

    do n = 1,nscal
       if (n.eq.Temp) then
          do i=lo,hi
             aofs(i,n) = ( macvel(i+1)*sedge(i+1,n) - macvel(i  )*sedge(i  ,n)) / dx
             aofs(i,n) = aofs(i,n) - &
                  (macvel(i+1)  - macvel(i)) * 0.5d0 * ( sedge(i,n) + sedge(i+1,n)) / dx
          end do
       else
          do i=lo,hi
             aofs(i,n) = ( macvel(i+1)*sedge(i+1,n) - macvel(i  )*sedge(i  ,n)) / dx
          end do
       end if

       ! XXX
       ! make these negative here so we can add as source terms later.
       do i=lo,hi
          aofs(i,n) = -aofs(i,n)
       end do
    end do
  end subroutine scal_aofs

  subroutine calc_divu(scal,beta,I_R,divu,dx,lo,hi)
    integer,          intent(in   ) :: lo, hi
    double precision, intent(in   ) :: scal(lo-2:hi+2,nscal)
    double precision, intent(in   ) :: beta(lo-1:hi+1,nscal)
    double precision, intent(in   ) :: I_R(lo:hi,nspec)
    double precision, intent(  out) :: divu(lo-1:hi+1)
    double precision, intent(in   ) :: dx

    double precision :: Y(Nspec), HK(Nspec), cpmix, mwmix
    double precision :: diff(lo-1:hi+1,nscal), diffdiff(lo-1:hi+1)
    double precision :: RWRK,rho,T, gamma_lo(lo:hi,Nspec), gamma_hi(lo:hi,Nspec)
    integer          :: IWRK,i,n

    ! compute Gamma_m
    call get_spec_visc_terms(scal,beta,diff(:,FirstSpec:), gamma_lo,gamma_hi,dx,lo,hi)

    ! compute div lambda grad T
    diff(:,Temp) = 0.d0
    call addDivLambdaGradT(scal,beta,diff(:,Temp),dx,lo,hi)

    ! compute div h_m Gamma_m
    call get_diffdiff_terms(scal,gamma_lo,gamma_hi, diffdiff,dx,lo,hi)

    ! combine div lambda grad T + div h_m Gamma_m
    do i=lo,hi
       diff(i,Temp) = diff(i,Temp) + diffdiff(i)
    end do

    do i=lo,hi
       rho = scal(i,Density)
       do n = 1,Nspec
          Y(n) = scal(i,FirstSpec + n - 1) / rho
       end do
       T = scal(i,Temp)
       call CKMMWY(Y,IWRK,RWRK,mwmix)
       call CKCPBS(T,Y,IWRK,RWRK,cpmix)
       call CKHMS(T,IWRK,RWRK,HK)

       divu(i) = diff(i,Temp)/(rho*cpmix*T)

       do n=1,Nspec
          divu(i) = divu(i) &
                     + (diff(i,FirstSpec+n-1) + I_R(i,n)) &
                     *(invmwt(n)*mwmix/rho - HK(n)/(rho*cpmix*T))
       end do
    end do

  end subroutine calc_divu


  subroutine get_spec_visc_terms(scal,beta,visc,gamma_lo, gamma_hi,dx,lo,hi)
    integer,          intent(in   ) :: lo, hi
    double precision, intent(in   ) :: scal(lo-2:hi+2,nscal)
    double precision, intent(in   ) :: beta(lo-1:hi+1,nscal)
    double precision, intent(  out) :: visc(lo-1:hi+1,nscal)
    double precision, intent(  out) :: gamma_lo(lo:hi,nspec)
    double precision, intent(  out) :: gamma_hi(lo:hi,nspec)
    double precision, intent(in   ) :: dx

    integer          :: i,n,is,IWRK
    double precision ::  beta_lo,beta_hi,RWRK
    double precision ::  dxsqinv
    double precision ::  Y(lo-1:hi+1,Nspec), sum_gamma_lo, sum_gamma_hi, sumRhoY_lo, sumRhoY_hi
    double precision ::  RhoYe_lo, RhoYe_hi
    double precision ::  X(lo-1:hi+1,Nspec)
    double precision ::  scal_X(lo-2:hi+2,nscal)

    dxsqinv = 1.d0/(dx*dx)

    do i=lo-1,hi+1

       ! compute Y = rho*Y / rho
       do n=1,Nspec
          Y(i,n) = scal(i,FirstSpec+n-1)/scal(i,Density)
       enddo

       ! convert Y to X
       CALL CKYTX(Y(i,:),IWRK,RWRK,X(i,:))

       ! compute rho*X
       do n=1,Nspec
          scal_X(i,FirstSpec+n-1) = scal(i,Density)*X(i,n)
       end do

    enddo

    do i=lo,hi

       sum_gamma_lo = 0.d0
       sum_gamma_hi = 0.d0
       sumRhoY_lo = 0.d0
       sumRhoY_hi = 0.d0

       do n=1,Nspec
          is = FirstSpec + n - 1

          ! compute beta on edges
          if (coef_avg_harm.eq.1) then
             beta_lo = 2.d0 / (1.d0/beta(i,is)+1.d0/beta(i-1,is))
             beta_hi = 2.d0 / (1.d0/beta(i,is)+1.d0/beta(i+1,is))
          else
             beta_lo = 0.5d0*(beta(i,is) + beta(i-1,is))
             beta_hi = 0.5d0*(beta(i,is) + beta(i+1,is))
          endif

          ! compute gamma
          gamma_hi(i,n) = beta_hi*(X(i+1,n) - X(i  ,n))
          gamma_lo(i,n) = beta_lo*(X(i  ,n) - X(i-1,n))

          ! compute div(gamma).  If non-unity Le we overwrite this later
          visc(i,n) = (gamma_hi(i,n)-gamma_lo(i,n))*dxsqinv

          if (LeEQ1 .eq. 0) then

             ! need to correct fluxes so they add to zero on each face
             ! build up the sum of species fluxes on lo and hi faces
             ! this will be "rho * V_c"
             sum_gamma_lo = sum_gamma_lo + gamma_lo(i,n)
             sum_gamma_hi = sum_gamma_hi + gamma_hi(i,n)

             ! build up the sum of rho*Y_m
             ! this will be the density
             sumRhoY_lo = sumRhoY_lo+0.5d0*(scal(i-1,is)+scal(i,is))
             sumRhoY_hi = sumRhoY_hi+0.5d0*(scal(i,is)+scal(i+1,is))

          end if

       enddo

       if (LeEQ1 .eq. 0) then
          ! correct the fluxes so they add up to zero before computing visc
          do n=1,Nspec
             is = FirstSpec + n - 1

             ! compute rho*Y_m on each face
             RhoYe_lo = .5d0*(scal(i-1,is)+scal(i,is))
             RhoYe_hi = .5d0*(scal(i,is)+scal(i+1,is))

             ! set flux = flux - (rho*V_c)*(rho*Y_m)/rho
             gamma_lo(i,n) = gamma_lo(i,n) - sum_gamma_lo*RhoYe_lo/sumRhoY_lo
             gamma_hi(i,n) = gamma_hi(i,n) - sum_gamma_hi*RhoYe_hi/sumRhoY_hi

             ! compute div(gamma)
             visc(i,n) = (gamma_hi(i,n)-gamma_lo(i,n))*dxsqinv

          end do
       end if

    end do

  end subroutine get_spec_visc_terms

  subroutine addDivLambdaGradT(scal,beta,visc,dx,lo,hi)
    integer,          intent(in   ) :: lo, hi
    double precision, intent(in   ) :: scal(lo-2:hi+2,nscal)
    double precision, intent(in   ) :: beta(lo-1:hi+1,nscal)
    double precision, intent(  out) :: visc(lo-1:hi+1)
    double precision, intent(in   ) :: dx

    integer          :: i
    double precision ::  beta_lo,beta_hi
    double precision ::  flux_lo,flux_hi
    double precision ::  dxsqinv

    dxsqinv = 1.d0/(dx*dx)
    do i=lo,hi
       if (coef_avg_harm.eq.1) then
          beta_lo = 2.d0 / (1.d0/beta(i,Temp)+1.d0/beta(i-1,Temp))
          beta_hi = 2.d0 / (1.d0/beta(i,Temp)+1.d0/beta(i+1,Temp))
       else
          beta_lo = 0.5*(beta(i,Temp) + beta(i-1,Temp))
          beta_hi = 0.5*(beta(i,Temp) + beta(i+1,Temp))
       endif

       flux_hi = beta_hi*(scal(i+1,Temp) - scal(i  ,Temp))
       flux_lo = beta_lo*(scal(i  ,Temp) - scal(i-1,Temp))
       visc(i) = visc(i) + (flux_hi - flux_lo) * dxsqinv

    enddo

  end subroutine addDivLambdaGradT


  subroutine get_diffdiff_terms(scal,gamma_lo,gamma_hi, diffdiff,dx,lo,hi)
    integer,          intent(in   ) :: lo, hi
    double precision, intent(in   ) :: scal(lo-2:hi+2,nscal)
    double precision, intent(  out) :: diffdiff(lo-1:hi+1)
    double precision, intent(in   ) :: gamma_lo(lo:hi,nspec)
    double precision, intent(in   ) :: gamma_hi(lo:hi,nspec)
    double precision, intent(in   ) :: dx


    integer          :: i,is,n,IWRK
    double precision ::  dxsqinv,RWRK
    double precision ::  hm(Nspec,lo-1:hi+1)
    double precision ::  flux_lo(Nspec),flux_hi(Nspec)

    dxsqinv = 1.d0/(dx*dx)

    diffdiff = 0.d0

    do i=lo-1,hi+1
       ! compute cell-centered h_m
       call CKHMS(scal(i,Temp),IWRK,RWRK,hm(1,i))
    end do

    do i=lo,hi
       do n=1,Nspec
          is = FirstSpec + n - 1

          ! set face fluxes to h_m * gamma_m
          flux_lo(n) = gamma_lo(i,n)*(hm(n,i-1)+hm(n,i))/2.d0
          flux_hi(n) = gamma_hi(i,n)*(hm(n,i+1)+hm(n,i))/2.d0

          ! differential diffusion is divergence of face fluxes
          diffdiff(i) = diffdiff(i) + (flux_hi(n) - flux_lo(n))*dxsqinv

       end do
    end do
  end subroutine get_diffdiff_terms

  subroutine calc_omega_dot(scal, omegadot, lo, hi)
    integer,          intent(in   ) :: lo, hi
    double precision, intent(in   ) :: scal(lo-2:hi+2,nscal)
    double precision, intent(  out) :: omegadot(lo:hi,nspec)

    double precision :: wdotk(nspec), C(nspec), rwrk
    integer          :: i, n, iwrk

    do i=lo,hi
       do n=1,Nspec
          C(n) = scal(i,FirstSpec+n-1)*invmwt(n)
       end do
       call CKWC(scal(i,Temp),C,IWRK,RWRK,WDOTK)
       do n=1,Nspec
          omegadot(i,n) = WDOTK(n)*mwt(n)
       end do
    end do
  end subroutine calc_omega_dot


  subroutine calc_diffusivities(scal, beta, beta_for_Y, beta_for_Wbar, mu, lo, hi)
    integer,          intent(in   ) :: lo, hi
    double precision, intent(in   ) :: scal(lo-2:hi+2,nscal)
    double precision, intent(  out) :: beta(lo-1:hi+1,nscal)
    double precision, intent(  out) :: beta_for_Y(lo-1:hi+1,nscal)
    double precision, intent(  out) :: beta_for_Wbar(lo-1:hi+1,nscal)
    double precision, intent(  out) :: mu(lo-1:hi+1)

    double precision :: Dt(Nspec), CPMS(Nspec), Y(Nspec)
    double precision :: Tt, Wavg, rho
    double precision :: X(Nspec), alpha, l1, l2, cpmix, RWRK
    integer          ::n, i, IWRK

    double precision, parameter :: fourThirds = 4.d0/3.d0

    ! !  Ensure chem/tran initialized
    ! if (traninit.lt.0) call initchem()

    if (LeEQ1 .eq. 0) then

       do i=lo-1,hi+1
          Tt = MAX(scal(i,Temp),TMIN_TRANS)
          rho = 0.d0
          do n=1,Nspec
             rho = rho + scal(i,FirstSpec+n-1)
          enddo
          do n=1,Nspec
             Y(n) = scal(i,FirstSpec+n-1) / rho
          enddo

          ! given y[species]: maxx fractions
          ! returns mean molecular weight (gm/mole)
          CALL CKMMWY(Y,IWRK,RWRK,Wavg)

          ! returns the specific heats at constant pressure
          ! in mass units
          CALL CKCPMS(Tt,IWRK,RWRK,CPMS)

          ! convert y[species] (mass fracs) to x[species] (mole fracs)
          CALL CKYTX(Y,IWRK,RWRK,X)

          ! initialize the thermomolecular parameters that are needed in order
          ! to evaluate the transport linear systems
          CALL EGSPAR(Tt,X,Y,CPMS,EGRWRK,EGIWRK)

          ! compute flux diffusion coefficients
          CALL EGSV1(Pcgs,Tt,Y,Wavg,EGRWRK,Dt)

          do n=1,Nspec
             beta         (i,FirstSpec+n-1) = rho*Dt(n)
             beta_for_Y   (i,FirstSpec+n-1) = rho*Wavg*invmwt(n)*Dt(n)
             beta_for_Wbar(i,FirstSpec+n-1) = rho*Y(n)*invmwt(n)*Dt(n)
          end do

          alpha = 1.0D0
          ! compute thermal conductivity
          CALL EGSL1(alpha, Tt, X, EGRWRK, l1)
          alpha = -1.0D0
          ! compute thermal conductivity with a different averating parameters
          CALL EGSL1(alpha, Tt, X, EGRWRK, l2)
          beta(i,Temp) = .5d0 * (l1 + l2)
          ! Returns the mean specific heat at CP
          CALL CKCPBS(scal(i,Temp),Y,IWRK,RWRK,CPMIX)
          beta(i,RhoH) = beta(i,Temp) / CPMIX

          ! compute shear viscosity
          CALL EGSE3(Tt, Y, EGRWRK, mu(i))
          mu(i) = fourThirds*mu(i)
       enddo
    else
       do i=lo-1,hi+1
          ! Kanuary, Combustion Phenomena (Wiley, New York) 1982:  mu [g/(cm.s)] = 10 mu[kg/(m.s)]
          mu(i) = 10.d0 * 1.85d-5*(MAX(scal(i,Temp),1.d0)/298.d0)**.7d0
          ! For Le=1, rho.D = lambda/cp = mu/Pr  (in general, Le = Sc/Pr)
          rho = 0.d0
          do n=1,Nspec
             rho = rho + scal(i,FirstSpec+n-1)
          enddo

          do n=1,Nspec
             Y(n) = scal(i,FirstSpec+n-1) / rho
          enddo

          ! given y[species]: maxx fractions
          ! returns mean molecular weight (gm/mole)
          CALL CKMMWY(Y,IWRK,RWRK,Wavg)

          ! returns the specific heats at constant pressure
          ! in mass units
          CALL CKCPMS(Tt,IWRK,RWRK,CPMS)

          ! convert y[species] (mass fracs) to x[species] (mole fracs)
          CALL CKYTX(Y,IWRK,RWRK,X)

          ! initialize the thermomolecular parameters that are needed in order
          ! to evaluate the transport linear systems
          CALL EGSPAR(Tt,X,Y,CPMS,EGRWRK,EGIWRK)

          ! compute flux diffusion coefficients
          CALL EGSV1(Pcgs,Tt,Y,Wavg,EGRWRK,Dt)

          do n=1,Nspec
             beta(i,FirstSpec+n-1) = mu(i) / (Sc * Wavg * invmwt(n))
             beta_for_Y(i,FirstSpec+n-1) = mu(i) / Sc
             beta_for_Wbar(i,FirstSpec+n-1) = mu(i) * Y(n) / (Sc * Wavg)
          end do

          ! Returns the mean specific heat at CP
          CALL CKCPBS(scal(i,Temp),Y,IWRK,RWRK,CPMIX)
          beta(i,RhoH) = mu(i) / Pr
          beta(i,Temp) = beta(i,RhoH) * CPMIX

          mu(i) = fourThirds*mu(i)
       enddo
    endif

  end subroutine calc_diffusivities

  subroutine get_vel_visc_terms(vel,beta,visc,dx,lo,hi)
    integer,          intent(in   ) :: lo, hi
    double precision, intent(in   ) :: vel(lo-2:hi+2)
    double precision, intent(in   ) :: beta(lo-1:hi+1)
    double precision, intent(  out) :: visc(lo-1:hi+1)
    double precision, intent(in   ) :: dx

    integer          :: i
    double precision :: beta_lo,beta_hi
    double precision :: flux_lo,flux_hi
    double precision :: dxsqinv

    ! Compute D(tau) = d/dx ( a . du/dx ), a=4.mu/3

    dxsqinv = 1.d0/(dx*dx)
    do i=lo,hi
       if (coef_avg_harm.eq.1) then
          beta_lo = 2.d0 / (1.d0/beta(i)+1.d0/beta(i-1))
          beta_hi = 2.d0 / (1.d0/beta(i)+1.d0/beta(i+1))
       else
          beta_lo = 0.5*(beta(i) + beta(i-1))
          beta_hi = 0.5*(beta(i) + beta(i+1))
       end if

       flux_hi = beta_hi*(vel(i+1) - vel(i  ))
       flux_lo = beta_lo*(vel(i  ) - vel(i-1))
       visc(i) = (flux_hi - flux_lo) * dxsqinv
    end do
  end subroutine get_vel_visc_terms

  ! subroutine get_spec_visc_terms_Wbar(scal,beta_for_Wbar,visc,
  !   &                                    gamma_Wbar_lo,gamma_Wbar_hi,
  !   &                                    dx,lo,hi)

  !   implicit none
  !   include 'spec.h'
  !   real*8          scal(-2:nfine+1,nscal)
  !   real*8 beta_for_Wbar(-1:nfine  ,nscal)
  !   real*8          visc(-1:nfine  ,Nspec)
  !   real*8 gamma_Wbar_lo( 0:nfine-1,Nspec)
  !   real*8 gamma_Wbar_hi( 0:nfine-1,Nspec)
  !   real*8 dx
  !   integer lo,hi

  !   integer i,n,is,IWRK
  !   real*8 beta_for_Wbar_lo,beta_for_Wbar_hi,RWRK
  !   real*8 dxsqinv
  !   real*8 Y(-1:nfine,Nspec)
  !   real*8 Wbar(-1:nfine)
  !   real*8 sum_gamma_lo, sum_gamma_hi, sumRhoY_lo, sumRhoY_hi
  !   real*8 RhoYe_lo, RhoYe_hi

  !   dxsqinv = 1.d0/(dx*dx)

  !   do i=lo-1,hi+1

  !      c     compute Y = rho*Y / rho
  !      do n=1,Nspec
  !         Y(i,n) = scal(i,FirstSpec+n-1)/scal(i,Density)
  !      enddo

  !      c     convert Y to Wbar
  !      CALL CKMMWY(Y(i,:),IWRK,RWRK,Wbar(i))

  !   enddo

  !   do i=lo,hi

  !      sum_gamma_lo = 0.d0
  !      sum_gamma_hi = 0.d0
  !      sumRhoY_lo = 0.d0
  !      sumRhoY_hi = 0.d0

  !      do n=1,Nspec
  !         is = FirstSpec + n - 1

  !         c     compute beta on edges
  !         if (coef_avg_harm.eq.1) then
  !            beta_for_Wbar_lo = 2.d0 / (1.d0/beta_for_Wbar(i,is)+1.d0/beta_for_Wbar(i-1,is))
  !            beta_for_Wbar_hi = 2.d0 / (1.d0/beta_for_Wbar(i,is)+1.d0/beta_for_Wbar(i+1,is))
  !         else
  !            beta_for_Wbar_lo = 0.5d0*(beta_for_Wbar(i,is) + beta_for_Wbar(i-1,is))
  !            beta_for_Wbar_hi = 0.5d0*(beta_for_Wbar(i,is) + beta_for_Wbar(i+1,is))
  !         endif

  !         c     compute gamma
  !         gamma_Wbar_hi(i,n) = beta_for_Wbar_hi*(Wbar(i+1) - Wbar(i  ))
  !         gamma_Wbar_lo(i,n) = beta_for_Wbar_lo*(Wbar(i  ) - Wbar(i-1))

  !         c     compute div(gamma).
  !         c     no need to conservatively correct these
  !         c     we will correct beta grad X after the species diffusion solve
  !         c     in fact the algorithm is more stable without the correction here
  !         visc(i,n) = (gamma_Wbar_hi(i,n)-gamma_Wbar_lo(i,n))*dxsqinv

  !      enddo

  !   end do

  ! end subroutine get_spec_visc_terms_Wbar

  ! subroutine get_spec_visc_terms_Y_and_Wbar(scal,beta_for_Y,visc,
  !   &                                          gamma_Wbar_lo,
  !   &                                          gamma_Wbar_hi,
  !   &                                          gamma_lo,
  !   &                                          gamma_hi,
  !   &                                          dx,lo,hi)

  !   c     compute
  !   c     gamma_m = beta_for_y grad Y + gamma_Wbar
  !   c     conservatively correct this, then set visc = (1/dxsq)*div(gamma_m)

  !   implicit none
  !   include 'spec.h'
  !   real*8          scal(-2:nfine+1,nscal)
  !   real*8    beta_for_Y(-1:nfine  ,nscal)
  !   real*8          visc(-1:nfine  ,Nspec)
  !   real*8 gamma_Wbar_lo( 0:nfine-1,Nspec)
  !   real*8 gamma_Wbar_hi( 0:nfine-1,Nspec)
  !   real*8 gamma_lo( 0:nfine-1,Nspec)
  !   real*8 gamma_hi( 0:nfine-1,Nspec)
  !   real*8 dx
  !   integer lo,hi

  !   integer i,n,is,IWRK
  !   real*8 beta_for_Y_lo,beta_for_Y_hi,RWRK
  !   real*8 dxsqinv
  !   real*8 Y(-1:nfine,Nspec), sum_gamma_lo, sum_gamma_hi, sumRhoY_lo, sumRhoY_hi
  !   real*8 RhoYe_lo, RhoYe_hi
  !   real*8 X(-1:nfine,Nspec)
  !   real*8 scal_X(-2:nfine+1,nscal)

  !   dxsqinv = 1.d0/(dx*dx)

  !   do i=lo-1,hi+1

  !      c     compute Y = rho*Y / rho
  !      do n=1,Nspec
  !         Y(i,n) = scal(i,FirstSpec+n-1)/scal(i,Density)
  !      enddo

  !      c     convert Y to X
  !      CALL CKYTX(Y(i,:),IWRK,RWRK,X(i,:))

  !      c     compute rho*X
  !      do n=1,Nspec
  !         scal_X(i,FirstSpec+n-1) = scal(i,Density)*X(i,n)
  !      end do

  !   enddo

  !   do i=lo,hi

  !      sum_gamma_lo = 0.d0
  !      sum_gamma_hi = 0.d0
  !      sumRhoY_lo = 0.d0
  !      sumRhoY_hi = 0.d0

  !      do n=1,Nspec
  !         is = FirstSpec + n - 1

  !         c     compute beta on edges
  !         if (coef_avg_harm.eq.1) then
  !            beta_for_Y_lo = 2.d0 / (1.d0/beta_for_Y(i,is)+1.d0/beta_for_Y(i-1,is))
  !            beta_for_Y_hi = 2.d0 / (1.d0/beta_for_Y(i,is)+1.d0/beta_for_Y(i+1,is))
  !         else
  !            beta_for_Y_lo = 0.5d0*(beta_for_Y(i,is) + beta_for_Y(i-1,is))
  !            beta_for_Y_hi = 0.5d0*(beta_for_Y(i,is) + beta_for_Y(i+1,is))
  !         endif

  !         c     compute gamma
  !         gamma_hi(i,n) = beta_for_Y_hi*(Y(i+1,n) - Y(i  ,n))
  !         gamma_lo(i,n) = beta_for_Y_lo*(Y(i  ,n) - Y(i-1,n))

  !         gamma_hi(i,n) = gamma_hi(i,n) + gamma_Wbar_hi(i,n)
  !         gamma_lo(i,n) = gamma_lo(i,n) + gamma_Wbar_lo(i,n)

  !         c     compute div(gamma).  If non-unity Le we overwrite this later
  !         visc(i,n) = (gamma_hi(i,n)-gamma_lo(i,n))*dxsqinv

  !         if (LeEQ1 .eq. 0) then

  !            c              need to correct fluxes so they add to zero on each face
  !            c              build up the sum of species fluxes on lo and hi faces
  !            c              this will be "rho * V_c"
  !            sum_gamma_lo = sum_gamma_lo + gamma_lo(i,n)
  !            sum_gamma_hi = sum_gamma_hi + gamma_hi(i,n)

  !            c              build up the sum of rho*Y_m
  !            c              this will be the density
  !            sumRhoY_lo = sumRhoY_lo+0.5d0*(scal(i-1,is)+scal(i,is))
  !            sumRhoY_hi = sumRhoY_hi+0.5d0*(scal(i,is)+scal(i+1,is))

  !         end if

  !      enddo

  !      if (LeEQ1 .eq. 0) then
  !         c           correct the fluxes so they add up to zero before computing visc
  !         do n=1,Nspec
  !            is = FirstSpec + n - 1

  !            c              compute rho*Y_m on each face
  !            RhoYe_lo = .5d0*(scal(i-1,is)+scal(i,is))
  !            RhoYe_hi = .5d0*(scal(i,is)+scal(i+1,is))

  !            c              set flux = flux - (rho*V_c)*(rho*Y_m)/rho
  !            gamma_lo(i,n) = gamma_lo(i,n)
  !            $              - sum_gamma_lo*RhoYe_lo/sumRhoY_lo
  !            gamma_hi(i,n) = gamma_hi(i,n)
  !            $              - sum_gamma_hi*RhoYe_hi/sumRhoY_hi

  !            c              compute div(gamma)
  !            visc(i,n) = (gamma_hi(i,n)-gamma_lo(i,n))*dxsqinv

  !         end do
  !      end if

  !   end do

  ! end subroutine get_spec_visc_terms_Y_and_Wbar

  !
  ! Solve XXX.
  !
  subroutine tridiag(a,b,c,r,u,gam,n)
    integer,          intent(in   ) :: n
    double precision, intent(in   ) :: a(n),b(n),c(n),r(n)
    double precision, intent(inout) :: u(n), gam(n)

    integer          :: j
    double precision :: bet

    if (b(1) .eq. 0) stop "CANT HAVE B(1) = ZERO"

    bet  = b(1)
    u(1) = r(1)/bet

    do j = 2,n
       gam(j) = c(j-1)/bet
       bet    = b(j) - a(j)*gam(j)
       if (bet .eq. 0)  stop 'TRIDIAG FAILED'
       u(j)   = (r(j)-a(j)*u(j-1))/bet
    end do

    do j = n-1,1,-1
       u(j) = u(j) - gam(j+1)*u(j+1)
    end do
  end subroutine tridiag

  !
  ! XXX.
  !
  subroutine mkslopes(scal,slope,lo,hi,bc)
    use lmc, only: unlim
    integer,          intent(in   ) :: lo, hi, bc(2)
    double precision, intent(in   ) :: scal(lo-2:hi+2)
    double precision, intent(  out) :: slope(lo-1:hi+1)

    double precision   :: slo,shi,slim,smid
    integer            :: i

    integer, parameter :: cen=1, lim=2, flag=3, fromm=4

    if (unlim .eq. 0) then
       do i=lo-1,hi+1
          shi = 2.0d0*(scal(i+1) - scal(i  ))
          slo = 2.0d0*(scal(i  ) - scal(i-1))
          smid = 0.5d0*(scal(i+1) - scal(i-1))
          if (i .eq. lo .and. bc(1) .eq. 1) then
             ! inflow: value in ghost cell is value at inflow face
             smid = (scal(1)+3.d0*scal(0)-4.d0*scal(-1))/3.d0
          end if
          slim = min(abs(slo),abs(shi))
          slim = min(abs(smid),slim)
          if (slo*shi .lt. 0.d0) then
             slope(i) = 0.d0
          else
             slope(i) = sign(1.d0,smid)*slim
          end if
       end do
    else
       do i=lo-1,hi+1
          slope(i) = 0.5d0*(scal(i+1) - scal(i-1))
       end do
    end if

    if (bc(1) .eq. 1) slope(lo-1) = 0.d0
    if (bc(2) .eq. 2) slope(hi:)  = 0.d0
  end subroutine mkslopes

end module kernels
