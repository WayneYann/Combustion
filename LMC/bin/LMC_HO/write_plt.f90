      subroutine write_plt(vel,scal_cc,divu,dx,nsteps,time)
      use cell_conversions_module
      implicit none
      include 'spec.h'

      integer nsteps
      real*8   vel( 0:nx)
      real*8  scal_cc(-2:nx+1,nscal)
      real*8  divu( 0:nx-1)
      real*8    dx
      real*8 time
      
      double precision :: scal_avg(-2:nx+1,nscal)
      
      real*8 Y(Nspec)
      character pltfile*(8)
      character char_of_int*(5)
      
      integer i,n
      pltfile(1:3) = 'rht'
      write(char_of_int,1005) nsteps
      pltfile(4:8) = char_of_int
 1005 format(i5.5)
 1006 FORMAT(200(E23.15E3,1X))
 
      call compute_pthermo(scal_cc(:,:))

      call scal_cc_to_avg(scal_avg, scal_cc)

      open(10,file=pltfile,form='formatted')
      print *,'...writing data to ',pltfile
      write(10,*) nsteps
      write(10,*) nx
      write(10,*) time

      do i=0,nx-1
         write(10,1006) (i+.5d0)*dx,(scal_avg(i,FirstSpec+n-1),n=1,Nspec),&
                           scal_avg(i,Density),&
                           scal_avg(i,Temp),&
                           scal_avg(i,RhoH),&
                           vel(i),&
                           scal_avg(i,RhoRT),&
                           divu(i)

      enddo

      close(10)

      end
