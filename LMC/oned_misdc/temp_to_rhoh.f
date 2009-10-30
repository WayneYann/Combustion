        subroutine temp_to_rhoh(lo,hi,lo_lim,hi_lim,scal)

        implicit none
 
        include 'nums.fi'
        include 'sndata.fi'
        include 'nkbrn.fi'

c       Quantities passed in.
        integer lo
        integer hi
        integer lo_lim
        integer hi_lim
        real*8 scal(lo_lim:hi_lim,nscal)

c       Local variables
        integer ins(nspec),outs(nspec)
        real*8  qs(nspec)
        real*8 aion(nspec),zion(nspec),xmass(nspec)
        real*8 pres,enthalpy,eint,c_v,c_p,ne,eta,pele
        real*8 dpdt,dpdr,dedt,dedr
        real*8 qreac,ttt
        integer i,n
        integer ispec

        do i = lo,hi
          scal(i,RhoH) = scal(i,RhoH) / scal(i,Density)          
          do n = 1,nspec
            ispec = FirstSpec - 1 + n
            scal(i,ispec) = scal(i,ispec) / scal(i,Density)          
          enddo
        enddo

        call nb_get_az(aion,zion,nspec)
        call nb_get_reac(ins,outs,qs,nreac)

        do i = lo,hi
           do n = 1,nspec
             ispec = FirstSpec - 1 + n
             xmass(n) = scal(i,ispec)
           enddo
           call eos(1, scal(i,Density), scal(i,Temp),
     &              nspec, xmass, aion, zion,
     &              pres, enthalpy, eint, c_v, c_p, ne, eta, pele,
     &              dpdt,dpdr,dedt,dedr)
           qreac = 0.0d0
           do n = 1, nreac
              qreac = qreac + scal(i,FirstSpec-1+ins(n))*qs(n)
           end do
           scal(i,RhoH) = enthalpy + qreac
        enddo

        do i = lo,hi
          scal(i,RhoH) = scal(i,RhoH) * scal(i,Density)          
          do n = 1,nspec
            ispec = FirstSpec - 1 + n
            scal(i,ispec) = scal(i,ispec) * scal(i,Density)          
          enddo
        enddo
        end
