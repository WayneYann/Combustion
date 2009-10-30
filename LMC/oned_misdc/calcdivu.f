        subroutine calc_divu(nx,scal,beta,Ydot,divu,dx)

        implicit none
 
        include 'nums.fi'
        include 'sndata.fi'
        include 'nkbrn.fi'

c       Quantities passed in
        integer nx
        real*8 scal(-1:nx,nscal)
        real*8 beta(-1:nx,nscal)
        real*8 divu(0:nx-1)
        real*8 Ydot(0:nx-1,nspec)
        real*8 dx

c       Local variables
        real*8 visc(0:nx-1)
        real*8   cp(0:nx-1)
        real*8 zion(nspec),aion(nspec),xmass(nspec)
        real*8 divu_node
        real*8 rho,tmp
        real*8 pres, enthalpy, eint, c_v, c_p, ne
        real*8 eta, pele, dpdt, dpdr, dedt, dedr
        real*8 dpdx
        real*8 dhdx,qreac
        real*8 t1, t2, t3

        integer i,n
        integer ins(nreac),outs(nreac)
        real*8  qs(nreac)

c       First create the del dot lambda terms 
        call rhoh_visc_terms(nx,scal,beta,visc,dx)
c       NOTE: WE DONT MODIFY WITH CP OR RHO HERE BECAUSE
c             THAT IS HANDLED BELOW.

        do n = FirstSpec,LastSpec
        do i = -1,nx
          scal(i,n) = scal(i,n) / scal(i,Density)
        end do
        end do

        do i = 0,nx-1
          divu(i) = visc(i)
        end do

        call nb_get_az(aion, zion, nspec)
        call nb_get_reac(ins,outs,qs,nreac)

        do i = 0,nx-1
            do n = 1, nspec
               xmass(n) = scal(i,FirstSpec-1+n)
            end do 
            rho  = scal(i,Density)
            tmp = scal(i,Temp)
            call eos(1, rho, tmp, nspec, xmass,
     &           aion, zion, 
     &           pres, enthalpy, eint, c_v, c_p, ne, eta, pele,
     &           dpdt,dpdr,dedt,dedr)

c           New code that used to be in addl_temp_visc_terms.
            t1 = kerg*tmp/amu
            t2 = t1*(1.5d0-(rho*dedr-pres/rho)/dpdr)
            dhdx = 0.d0
            do n = 1, nspec
               dhdx = dhdx + Ydot(i,n)/aion(n)
            end do
            qreac = 0.d0
            do n = 1, nreac
               qreac = qreac + qs(n)*Ydot(i,ins(n))
            end do
            dhdx = dhdx*t2 + qreac
c           Note that commenting out this line balances out the not
c                dividing by rho above.
c           dhdx = rho*dhdx
c           Note that dividing divu as well as dhdx by c_p here takes
c                care of not dividing divu by c_p above. 
            divu(i) = (divu(i) - dhdx)/c_p
c           End of new code.

            t1 = kerg*tmp/amu
            t2 = t1*rho
            dpdx = 0
            do n = 1, nspec
               dpdx = dpdx + Ydot(i,n)/aion(n)
            end do
c           dpdx = rho*dpdx*t2
            dpdx = dpdx*t2
            divu(i) = divu(i)*dpdt + dpdx
            divu(i) = divu(i)/rho**2/dpdr
        end do

        do n = FirstSpec,LastSpec
        do i = -1,nx
          scal(i,n) = scal(i,n) * scal(i,Density)
        end do
        end do

        divu(nx-1) = 0.d0

        end
