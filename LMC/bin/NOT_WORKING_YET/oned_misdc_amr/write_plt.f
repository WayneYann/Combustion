      subroutine write_plt(vel,scal,press,divu,I_R,dx,nsteps,time,lo,hi)
      implicit none
      include 'spec.h'

      integer nsteps
      real*8   vel(0:nlevs-1,-2:nx+1)
      real*8  scal(0:nlevs-1,-2:nx+1,nscal)
      real*8 press(0:nlevs-1,-1:nx+1)
      real*8  divu(0:nlevs-1, 0:nx-1)
      real*8   I_R(0:nlevs-1,-1:nx  ,0:Nspec)
      real*8    dx(0:nlevs-1)
      real*8 time
      integer lo(0:nlevs-1)
      integer hi(0:nlevs-1)
      
      real*8 Y(Nspec)
      character pltfile*(8)
      character char_of_int*(5)
      
      integer i,l,n
      pltfile(1:3) = 'rht'
      write(char_of_int,1005) nsteps
      pltfile(4:8) = char_of_int
 1005 format(i5.5)
 1006 FORMAT(200(E23.15E3,1X))
      do l=0,nlevs-1
         call compute_pthermo(scal(l,:,:),scal(l,:,RhoRT),lo(l),hi(l))
      end do

      open(10,file=pltfile,form='formatted')
      print *,'...writing data to ',pltfile
      write(10,*) nsteps
      write(10,*) nx
      write(10,*) time

      do l=0,nlevs-1
         do i=lo(l),hi(l)
            do n=1,Nspec
               Y(n) = scal(l,i,FirstSpec+n-1)/scal(l,i,Density)
            enddo
            write(10,1006) (i+.5)*dx(l),(Y(n),n=1,Nspec),
     $                     scal(l,i,Density),
     $                     scal(l,i,Temp),
     $                     scal(l,i,RhoH),
     $                     vel(l,i),
     $                     press(l,i),
     $                     scal(l,i,RhoRT),
     $                     divu(l,i),
     $                     (I_R(l,i,n),n=1,Nspec)
         enddo
      enddo

      close(10)

      end
