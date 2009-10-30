        subroutine addl_temp_visc_terms(nx,scal,visc,dx)

        implicit none
 
        include 'nums.fi'
        include 'sndata.fi'
        include 'nkbrn.fi'

c       Quantities passed in
        integer nx
        real*8 scal(-1:nx  ,nscal)
        real*8 visc(0 :nx-1)
        real*8 dx

c       Local variables
        integer ins(nreac), outs(nreac)
        real*8  qs(nreac)
        real*8 pres, enthalpy, eint, c_v, c_p, ne
        real*8 eta, pele, dpdt, dpdr, dedt, dedr
        real*8 zion(nspec),aion(nspec),xmass(nspec)

        real*8 tmp,rho,t1,t2,dhdx,qreac
        integer i,n

        call nb_get_reac(ins,outs,qs,nreac)
        call nb_get_az(aion, zion, nspec)

        do i = 0,nx-1
            rho  = scal(i,Density)
            tmp  = scal(i,Temp)
            do n = 1, nspec
               xmass(n) = scal(i,FirstSpec-1+n)/rho
            end do 
            call eos(1, rho, tmp, nspec, xmass,
     &           aion, zion, 
     &           pres, enthalpy, eint, c_v, c_p, ne, eta, pele,
     &           dpdt,dpdr,dedt,dedr)

            visc(i) = visc(i) / c_p
        end do

        end
