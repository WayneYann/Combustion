        subroutine rhoh_to_temp(nx,scal)

        implicit none
 
        include 'nums.fi'
        include 'sndata.fi'
        include 'nkbrn.fi'

c       Quantities passed in.
        integer nx
        real*8 scal(-1:nx,nscal)

c       Local variables
        integer ins(nspec),outs(nspec)
        real*8  qs(nspec)
        real*8 aion(nspec),zion(nspec),xmass(nspec)
        real*8 pres,enthalpy,eint,c_v,c_p,ne,eta,pele
        real*8 dpdt,dpdr,dedt,dedr
        real*8 qreac
        integer i,n
        integer ispec
        integer cnt,itemp

        do i = -1,nx
          scal(i,RhoH) = scal(i,RhoH) / scal(i,Density)          
          do n = 1,nspec
            ispec = FirstSpec - 1 + n
            scal(i,ispec) = scal(i,ispec) / scal(i,Density)          
          enddo
        enddo
 

        call nb_get_az(aion,zion,nspec)
        call nb_get_reac(ins,outs,qs,nreac)

        cnt = 0
        do i = -1,nx
           do n = 1,nspec
             ispec = FirstSpec - 1 + n
             xmass(n) = scal(i,ispec)
           enddo
           qreac = 0.0d0
           do n = 1, nreac
              ispec = FirstSpec-1+ins(n)
              qreac = qreac + scal(i,ispec)*qs(n)
           end do
           enthalpy = scal(i,RhoH) - qreac
           itemp = 2
           call eos(itemp, scal(i,Density), scal(i,Temp),
     &              nspec, xmass, aion, zion,
     &              pres, enthalpy, eint, c_v, c_p, ne, eta, pele,
     &              dpdt,dpdr,dedt,dedr)
            cnt = max(cnt,itemp)
        enddo

        if (cnt .lt. 0) then
          print *,'ERROR: RHOH_TO_TEMP '
          stop
        else
          print *,'RHOH_TO_TEMP: max_iters = ',cnt
        endif

        do i = -1,nx
          scal(i,RhoH) = scal(i,RhoH) * scal(i,Density)          
          do n = 1,nspec
            ispec = FirstSpec - 1 + n
            scal(i,ispec) = scal(i,ispec) * scal(i,Density)          
          enddo
        enddo

        scal(nx,Temp) = scal(nx-1,Temp)

        end
