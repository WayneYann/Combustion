        subroutine writecheck(nsteps,nx,vel,scal,press,
     $                        Ydot,divu,dsdt,intra,
     $                        dx,time,dt_old,cfl_used)

        implicit none
 
        include 'nums.fi'
        include 'sndata.fi'

        integer nsteps
        integer nx
        real*8   vel(-1:nx)
        real*8  scal(-1:nx,nscal)
        real*8 press(0:nx)
        real*8  Ydot( 0:nx-1,nscal)
        real*8  divu( 0:nx-1)
        real*8  dsdt( 0:nx-1)
        real*8 intra( 0:nx-1,nscal)
        real*8 dx
        real*8 time
        real*8 dt_old
        real*8 cfl_used

        character chkfile*(8)
        character char_of_int*(5)
        integer i,n

        chkfile(1:3) = 'chk'

        write(char_of_int,1005) nsteps
        chkfile(4:8) = char_of_int
1005    format(i5.5)
        print *,'...writing checkpoint ',chkfile

        open(10,file=chkfile,form='unformatted')
        write(10) nsteps
        write(10) time
        write(10) dt_old
        write(10) cfl_used
        do i = -1,nx
          write(10) (i+.5)*dx,scal(i,FirstSpec),
     $                        scal(i,FirstSpec+1),
     $                        scal(i,FirstSpec+2),
     $                        scal(i,Density),
     $                        scal(i,Temp),
     $                        scal(i,RhoH),
     $                        vel(i)
        enddo
        do i = 0,nx-1
          write(10) (i+.5)*dx,Ydot(i,1),Ydot(i,2),Ydot(i,3),
     $                        divu(i),dsdt(i)
        enddo
        do i = 0,nx
          write(10) (i+.5)*dx,press(i)
        enddo
        do n = 1,nscal
        do i = 0,nx-1
          write(10) intra(i,n)
        enddo
        enddo

        end
