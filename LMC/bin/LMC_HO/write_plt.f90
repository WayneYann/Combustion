      subroutine write_plt(vel,scal_avg,divu,dx,nsteps,time)
      use cell_conversions_module
      implicit none
      include 'spec.h'

      integer nsteps
      real*8   vel( 0:nx)
      real*8  scal_avg(-2:nx+1,nscal)
      real*8  divu( 0:nx-1)
      real*8    dx
      real*8 time
      
      double precision :: Y(-2:nx+1, Nspec)
      
      !double precision :: scal_avg(-2:nx+1,nscal)
      double precision :: scal_cc(-2:nx+1,nscal)
      
      double precision :: h_avg(-2:nx+1)
      
      character pltfile*(8)
      character char_of_int*(5)
      
      integer i,n
      
      pltfile(1:3) = 'rht'
      write(char_of_int,1005) nsteps
      pltfile(4:8) = char_of_int
 1005 format(i5.5)
 1006 FORMAT(200(E23.15E3,1X))
 
      
      call scal_avg_to_cc(scal_cc, scal_avg)
      call compute_pthermo(scal_cc)
      call extrapolate_cc_to_avg(scal_avg(0:,RhoRT), scal_cc(0:,RhoRT))
      
      open(10,file=pltfile,form='formatted')
      print *,'...writing data to ',pltfile
      write(10,*) nsteps
      write(10,*) nx
      write(10,*) time
      
      do n=1,Nspec
         call divide_avgs(Y(:,n), scal_avg(:,FirstSpec+n-1), scal_avg(:,Density), Y_bc(n,on_lo))
      end do
      
      call divide_avgs(h_avg, scal_avg(:,RhoH), scal_avg(:,Density), h_bc(on_lo))
      
      do i=0,nx-1
         write(10,1006) (i+.5d0)*dx,(Y(i,n),n=1,Nspec),&
                           scal_avg(i,Density),&
                           scal_avg(i,Temp),&
                           scal_avg(i,RhoH),&
                           vel(i),&
                           scal_avg(i,RhoRT),&
                           divu(i),&
                           i*dx

      enddo
      
      close(10)

      end
