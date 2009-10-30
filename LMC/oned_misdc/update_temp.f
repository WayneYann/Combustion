        subroutine update_temp(nx,scal_old,scal_new,
     $                         aofs,alpha,beta,Rhs,dx,dt,be_cn_theta)

        implicit none

        include 'nums.fi'
        include 'sndata.fi'

        integer nx

c       Quantities passed in
        real*8 scal_old(-1:nx  ,nscal)
        real*8 scal_new(-1:nx  ,nscal)
        real*8     aofs(0 :nx-1,nscal)
        real*8    alpha(0 :nx-1)
        real*8     beta(-1:nx  ,nscal)
        real*8      Rhs(0 :nx-1)
        real*8       cp(0 :nx-1)
        real*8 dx
        real*8 dt
        real*8 be_cn_theta

c       Local variables
        real*8  dth,dxsqinv
        real*8  scal_mid(-1:nx  ,nscal)
        real*8  visc_old(0 :nx-1)
        real*8  visc_new(0 :nx-1)
        real*8  delta_rhs(0 :nx-1)
        real*8 rhohalf
        real*8 beta_lo,beta_hi
        real*8 visc_term
        integer i,n,nn

        n        = Temp

        dth = 0.5d0 * dt
        dxsqinv = 1.d0/(dx*dx)

        be_cn_theta = 0.5d0

c*************************************************************************
c       Create Snew = Sold - dt*aofs
c*************************************************************************

        do i = 0,nx-1
          scal_new(i,n) = scal_old(i,n) + dt * aofs(i,n)
c         write(6,1000) i,aofs(i,n)
        enddo
1000    format('AOFS',i4,2x,f40.20)

c*************************************************************************
c       Create RHS = time n diffusive term.
c*************************************************************************
        do i = 0,nx-1

          beta_lo = 2.d0 / (1.d0/beta(i,Temp)+1.d0/beta(i-1,Temp))
          beta_hi = 2.d0 / (1.d0/beta(i,Temp)+1.d0/beta(i+1,Temp))

          visc_term =  
     $       beta_hi*(scal_old(i+1,Temp)-scal_old(i  ,Temp)) -
     $       beta_lo*(scal_old(i  ,Temp)-scal_old(i-1,Temp))
          Rhs(i) = visc_term * dxsqinv * dth

        enddo

c*************************************************************************
c       Construct alpha.
c*************************************************************************
        do i = -1,nx
          do nn = 1,nscal
            scal_mid(i,nn) = 0.5d0 * (scal_old(i,nn) + scal_new(i,nn))
          enddo
        enddo
        call compute_cp(nx,cp,scal_mid)
        do i = 0,nx-1
          rhohalf = 0.5d0 * (scal_old(i,Density) + scal_new(i,Density))
          alpha(i) = rhohalf * cp(i)
        enddo


c*************************************************************************
c       Update Rhs.
c*************************************************************************

c       do i = 0,nx-1
c         Rhs(i) = Rhs(i) + dt*delta_rhs(i)
c       enddo

        do i = 0,nx-1
          Rhs(i) = Rhs(i) + alpha(i) * scal_new(i,Temp)
        enddo

        end

c*************************************************************************
c*************************************************************************
c       Compute c_p
c*************************************************************************
c*************************************************************************

        subroutine compute_cp(nx,cp,scal)
        implicit none

        include 'nums.fi'
        include 'sndata.fi'

c       Quantities passed in
        integer nx
        real*8   cp(0 :nx-1)
        real*8 scal(-1:nx  ,nscal)

c       Local variables
        real*8 pres, enthalpy, eint, c_v, c_p, ne
        real*8 eta, pele, dpdt, dpdr, dedt, dedr
        real*8 zion(nspec),aion(nspec),xmass(nspec)
        real*8 rho,tmp
        integer i,n

        call nb_get_az(aion, zion, nspec)

        do i = 0,nx-1
          rho  = scal(i,Density)
          tmp  = scal(i,Temp)
          do n = 1, nspec
             xmass(n) = scal(i,FirstSpec-1+n)/rho
          end do
          call eos(1, rho, tmp, nspec, xmass,
     &         aion, zion,
     &         pres, enthalpy, eint, c_v, c_p, ne, eta, pele,
     &         dpdt,dpdr,dedt,dedr)
          cp(i) = c_p
        enddo

        end
