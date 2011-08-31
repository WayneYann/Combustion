        subroutine macproj(nx,macvel,divu,dx)

        implicit none
 
c       Quantities passed in
        integer nx
        real*8 macvel(0 :nx)
        real*8   divu(0 :nx-1)
        real*8 dx

c       Local variables
        integer i

        print *,'... mac_projection'

        do i = 1,nx
          macvel(i) = macvel(i-1) + divu(i-1)*dx
        end do

c       macvel(0) was set in predict_vel and remains unchanged

        end
