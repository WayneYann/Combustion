      subroutine write_check(nsteps,vel,scal,press,
     $                       I_R,divu,dsdt,dx,time,dt_old,lo,hi)
      implicit none
      include 'spec.h'
      integer nsteps
      real*8   vel(0:nlevs-1,-2:nfine+1)
      real*8  scal(0:nlevs-1,-2:nfine+1,nscal)
      real*8 press(0:nlevs-1,-1:nfine+1)
      real*8   I_R(0:nlevs-1,-1:nfine  ,0:Nspec)
      real*8  divu(0:nlevs-1, 0:nfine-1)
      real*8  dsdt(0:nlevs-1, 0:nfine-1)
      real*8 dx
      real*8 time
      real*8 dt_old
      integer lo(0:nlevs-1), hi(0:nlevs-1)
      
      character chkfile*(8)
      character char_of_int*(5)
      integer i,l,n

      chkfile(1:3) = 'chk'
      
      write(char_of_int,1005) nsteps
      chkfile(4:8) = char_of_int
 1005 format(i5.5)
      print *,'...writing checkpoint ',chkfile
      
      open(10,file=chkfile,form='unformatted')
      write(10) nsteps
      write(10) time
      write(10) dt_old
      write(10) hmix_TYP
      write(10) vel_TYP

      do l=0,nlevs-1

!     cell-centered, 2 ghost cells
         do i=lo(l)-2,hi(l)+2
            write(10) (i+.5)*dx,(scal(l,i,FirstSpec+n),n=0,Nspec-1),
     $                           scal(l,i,Density),
     $                           scal(l,i,Temp),
     $                           scal(l,i,RhoH),
     $                           vel(l,i)
         enddo

!     cell-centered, 1 ghost cell
         do i=lo(l)-1,hi(l)+1
            write(10) (i+.5)*dx,(I_R(l,i,n),n=0,Nspec)
         enddo

!     cell-centered, no ghost cells
         do i=lo(l),hi(l)
            write(10) (i+.5)*dx,divu(l,i),dsdt(l,i)
         enddo

!     nodal, 1 ghost cell
         do i=lo(l)-1,hi(l)+2
            write(10) i*dx,press(l,i)
         enddo

      enddo

      close(10)

      end
