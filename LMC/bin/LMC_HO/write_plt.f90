      subroutine write_plt(vel,scal,divu,dx,nsteps,time)
      implicit none
      include 'spec.h'

      integer nsteps
      real*8   vel( 0:nx)
      real*8  scal(-2:nx+1,nscal)
      real*8  divu( 0:nx-1)
      real*8    dx
      real*8 time
      
      real*8 Y(Nspec)
      character pltfile*(8)
      character char_of_int*(5)
      
      integer i,n
      pltfile(1:3) = 'rht'
      write(char_of_int,1005) nsteps
      pltfile(4:8) = char_of_int
 1005 format(i5.5)
 1006 FORMAT(200(E23.15E3,1X))
 
      call compute_pthermo(scal(:,:))

      open(10,file=pltfile,form='formatted')
      print *,'...writing data to ',pltfile
      write(10,*) nsteps
      write(10,*) nx
      write(10,*) time

      do i=0,nx-1
         do n=1,Nspec
            Y(n) = scal(i,FirstSpec+n-1)/scal(i,Density)
         enddo
         write(10,1006) (i+.5)*dx,(Y(n),n=1,Nspec),&
                           scal(i,Density),&
                           scal(i,Temp),&
                           scal(i,RhoH),&
                           vel(i),&
                           scal(i,RhoRT),&
                           divu(i)
      enddo

      close(10)

      end
