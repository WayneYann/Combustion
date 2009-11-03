      subroutine read_check(chkfile,nx,vel,scal,press,
     $                      Ydot,divu,dsdt,intra,
     $                      time,at_nstep,dt_old,cfl_used)
      implicit none
      include 'spec.h'
      character chkfile*(16)
      integer at_nstep
      integer nx
      real*8   vel(-1:nx)
      real*8  scal(-1:nx  ,nscal)
      real*8  Ydot( 0:nx-1,maxspec)
      real*8  divu( 0:nx-1)
      real*8  dsdt( 0:nx-1)
      real*8 intra( 0:nx-1,nscal)
      real*8 press(0:nx)
      real*8 time
      real*8 dt_old
      real*8 cfl_used
      
      real*8 x
      integer i,n
      
      print *,'READING CHKFILE ',chkfile
      
      open(10,file=chkfile,form='unformatted',status='old')
      read(10) at_nstep
      read(10) time
      read(10) dt_old
      read(10) cfl_used
      do i = -1,nx
         read(10) x,(scal(i,FirstSpec+n),n=0,Nspec-1),
     $            scal(i,Density),
     $            scal(i,Temp),
     $            scal(i,RhoH),
     $            vel(i)
      enddo
      do i = 0,nx-1
         read(10) x,(Ydot(i,n),n=1,Nspec),divu(i),dsdt(i)
      enddo
      do i = 0,nx
         read(10) x,press(i)
      enddo
      do n = 1,nscal
         do i = 0,nx-1
            read(10) intra(i,n)
         enddo
      enddo
      end
      
