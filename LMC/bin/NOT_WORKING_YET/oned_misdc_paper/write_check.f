      subroutine write_check(nsteps,vel,scal,press,
     $                       I_R,divu,dsdt,dx,time,dt_old,cfl_used)
      implicit none
      include 'spec.h'
      integer nsteps
      real*8   vel(-1:nx)
      real*8  scal(-1:nx,nscal)
      real*8 press(0:nx)
      real*8   I_R(0:nx-1,0:maxspec)
      real*8  divu( 0:nx-1)
      real*8  dsdt( 0:nx-1)
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
 1005 format(i5.5)
      print *,'...writing checkpoint ',chkfile
      
      open(10,file=chkfile,form='unformatted')
      write(10) nsteps
      write(10) time
      write(10) dt_old
      write(10) cfl_used
      write(10) hmix_TYP
      do i = -1,nx
         write(10) (i+.5)*dx,(scal(i,FirstSpec+n),n=0,Nspec-1),
     $                        scal(i,Density),
     $                        scal(i,Temp),
     $                        scal(i,RhoH),
     $                        vel(i)
      enddo
      do i = 0,nx-1
         write(10) (i+.5)*dx,(I_R(i,n),n=0,Nspec),divu(i),dsdt(i)
      enddo
      do i = 0,nx
         write(10) i*dx,press(i)
      enddo

      close(10)

      end
