        subroutine compute_rhort(nx,scal,rhort)

        implicit none
 
        include 'nums.fi'
        include 'sndata.fi'
        include 'nkbrn.fi'

c       Quantities passed in
        integer nx
        real*8  scal(-1:nx,nscal)
        real*8 rhort(-1:nx)

c       Local variables
        real*8 aion(nspec),zion(nspec),xmass(nspec)
        real*8 pres, enthalpy, eint, c_v, c_p, ne
        real*8 eta, pele, dpdt, dpdr, dedt, dedr
        real*8 rhort_min,rhort_max
        integer i,n
        integer ispec

        call nb_get_az(aion,zion,nspec)

c       Define p = rho*R*t from call to eos.
        do i = 0,nx-1
          do n = 1,nspec
            ispec = FirstSpec-1+n
            xmass(n) = scal(i,ispec) / scal(i,Density)
          end do
          call eos(1,scal(i,Density),scal(i,Temp),
     $             nspec,xmass,aion,zion,
     $             rhort(i),enthalpy,eint,c_v,c_p,ne,eta,pele,
     $             dpdt,dpdr,dedt,dedr)
          if (i.eq.0) then
            rhort_min = rhort(i)
            rhort_max = rhort(i)
          else
            rhort_min = min(rhort_min,rhort(i))
            rhort_max = max(rhort_max,rhort(i))
          endif
        end do
        rhort(-1) = rhort(0)
        rhort(nx) = rhort(nx-1)

        print *,'RHORT MIN/MAX ',rhort_min,rhort_max

        end
