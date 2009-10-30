        subroutine calc_diffusivities(nx,scal,beta)

        implicit none
 
        include 'nums.fi'
        include 'sndata.fi'
        include 'nkbrn.fi'

c       Quantities passed in.
        integer nx
        real*8  scal(-1:nx,nscal)
        real*8  beta(-1:nx,nscal)

c       Local variables
        integer i,n
        real*8 zion(nspec),aion(nspec),xmass(nspec)
        real*8 cond,diff,c_p

        call nb_get_az(aion, zion, nspec)

        do n = 1,nscal
            do i = -1,nx
               beta(i,n) = 0.d0
            enddo
        enddo

        do i = -1,nx
          do n = 1, nspec
            xmass(n) = scal(i,FirstSpec-1+n) / scal(i,Density)
          end do
          call conductivity(scal(i,Temp),scal(i,Density),
     &                      xmass,zion,aion,nspec,cond,diff,c_p)
          beta(i,RhoH) = cond/c_p
          beta(i,Temp) = cond
        end do

        end
