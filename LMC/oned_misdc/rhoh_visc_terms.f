        subroutine rhoh_visc_terms(nx,scal,beta,visc,dx)

        implicit none
 
        include 'nums.fi'
        include 'sndata.fi'
        include 'nkbrn.fi'

c       Quantities passed in.
        integer nx
        real*8            scal(-1:nx  ,nscal)
        real*8            beta(-1:nx  ,nscal)
        real*8 spec_grad_terms(0 :nx-1)
        real*8            visc(0 :nx-1)
        real*8           dx

c       Local variables
        integer i,n
        integer ispec
        real*8 beta_lo,beta_hi
        real*8 flux_lo,flux_hi
        real*8 dxsqinv

        dxsqinv = 1.d0/(dx*dx)

        do i = -1,nx
          scal(i,RhoH) = scal(i,RhoH) / scal(i,Density)
        enddo

        do i = 0,nx-1
            beta_lo = 2.d0 / (1.d0/beta(i,RhoH)+1.d0/beta(i-1,RhoH))
            beta_hi = 2.d0 / (1.d0/beta(i,RhoH)+1.d0/beta(i+1,RhoH))
            flux_hi = beta_hi*(scal(i+1,RhoH) - scal(i  ,RhoH)) 
            flux_lo = beta_lo*(scal(i  ,RhoH) - scal(i-1,RhoH)) 
            visc(i) =  flux_hi - flux_lo
            visc(i) = visc(i) * dxsqinv
        end do

        do i = -1,nx
          scal(i,RhoH) = scal(i,RhoH) * scal(i,Density)
        enddo

        call get_spec_grad_terms(nx,scal,beta,spec_grad_terms,dx)

        do i = 0,nx-1
            visc(i) = visc(i) - spec_grad_terms(i)
        end do

        end

        subroutine get_spec_grad_terms(nx,scal,beta,spec_grad_terms,dx)

        implicit none
 
        include 'nums.fi'
        include 'sndata.fi'
        include 'nkbrn.fi'

c       Quantities passed in/out.
        integer nx
        real*8            scal(-1:nx  ,nscal)
        real*8            beta(-1:nx  ,nscal)
        real*8 spec_grad_terms(0 :nx-1)
        real*8           dx

c       Local variables
        integer i,n
        integer ispec
        real*8 beta2(-1:nx,nspec)
        real*8 beta_lo
        real*8 beta_hi
        real*8 dxsqinv

        dxsqinv = 1.d0/(dx*dx)

        call sn_dhdxi(nx,scal,beta2)

        do i = 0,nx
          do n = 1,nspec
            beta2(i,n) = beta2(i,n) * beta(i,RhoH)
          end do
        end do

        do i = 0,nx-1
           spec_grad_terms(i) = 0.d0
        enddo

        do n = 1,nspec
          do i = 0,nx
             ispec = FirstSpec-1+n
             scal(i,ispec) = scal(i,ispec) / scal(i,Density)
          enddo
        enddo

        do n = 1,nspec
          do i = 0,nx-1
             ispec = FirstSpec-1+n
             beta_lo = 2.d0 / (1.d0/beta2(i,n)+1.d0/beta2(i-1,n))
             beta_hi = 2.d0 / (1.d0/beta2(i,n)+1.d0/beta2(i+1,n))
             spec_grad_terms(i) = spec_grad_terms(i) + 
     $        ( beta_hi*(scal(i+1,ispec) - scal(i  ,ispec))
     $        - beta_lo*(scal(i  ,ispec) - scal(i-1,ispec)) )*dxsqinv
          end do
        end do

        do n = 1,nspec
          do i = 0,nx
             ispec = FirstSpec-1+n
             scal(i,ispec) = scal(i,ispec) * scal(i,Density)
          enddo
        enddo

        end

        subroutine sn_dhdxi(nx,scal,beta2)

        implicit none
 
        include 'nums.fi'
        include 'sndata.fi'
        include 'nkbrn.fi'

c       Quantities passed in.
        integer nx
        real*8   scal(-1:nx,nscal)
        real*8  beta2(-1:nx,nspec)

c       Local variables
        integer i,n
        integer ins(MAXREAC), outs(MAXREAC)
        real*8 zion(nspec),aion(nspec),xmass(nspec)
        real*8 qs(MAXREAC)
        real*8 pres, enthalpy, eint, c_v, c_p, ne
        real*8 eta, pele, dpdt, dpdr, dedt, dedr
        real*8 rho
        real*8 t1,t3

        call nb_get_az(aion, zion, nspec)
        call nb_get_reac(ins,outs,qs,nreac)

        do i = 0,nx
           do n = 1, nspec
              xmass(n) = scal(i,FirstSpec-1+n) / scal(i,Density)
           end do
           rho = scal(i,Density)
           call eos(1, rho, scal(i,Temp), nspec, xmass,
     &              aion, zion,
     &              pres, enthalpy, eint, c_v, c_p, ne, eta, pele,
     &               dpdt,dpdr,dedt,dedr)


            t1 = kerg*scal(i,Temp)/amu
            t3 = t1*(1.5D0-(rho*dedr-pres/rho)/dpdr)
            do n = 1, nspec
               beta2(i,n) = t3/aion(n)
            end do
            do n = 1, nreac
               beta2(i,ins(n)) = beta2(i,ins(n)) + qs(n)
            end do
        end do
        end
