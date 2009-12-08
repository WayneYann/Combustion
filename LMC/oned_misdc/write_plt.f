      subroutine write_plt(vel,scal,press,dx,nsteps,time)
      implicit none
      include 'spec.h'

      integer nsteps
      real*8   vel(-1:nx)
      real*8  scal(-1:nx,nscal)
      real*8 press( 0:nx)
      real*8 ptherm(-1:nx)
      real*8 dx
      real*8 time
      
      character pltfile*(8)
      character char_of_int*(5)
      
      integer i,n,nvars
      pltfile(1:3) = 'plt'
      
      write(char_of_int,1005) nsteps
      pltfile(4:8) = char_of_int
 1005 format(i5.5)
 1006 FORMAT(16(E15.8,1X))      
      call compute_pthermo(scal,ptherm)
      
      open(10,file=pltfile,form='formatted')
      print *,'...writing data to ',pltfile
      write(10,*) nsteps
      write(10,*) nx
      write(10,*) time
C CEG:: changing this to output ghost cells too
      do i = -1,nx
         write(10,1006) (i+.5)*dx,(scal(i,FirstSpec+n),n=0,Nspec-1),
     $                         scal(i,Density),
     $                         scal(i,Temp),
     $                         scal(i,RhoH),
     $                         vel(i),
     $                         ptherm(i),
     $                         (ptherm(i)-Pcgs)
      enddo

      end
