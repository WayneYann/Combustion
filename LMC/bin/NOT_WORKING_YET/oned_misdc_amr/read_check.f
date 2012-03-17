      subroutine read_check(chkfile,vel,scal,press,
     $                      I_R,divu,dsdt,
     $                      time,at_nstep,dt_old,lo,hi)
      implicit none
      include 'spec.h'
      character chkfile*(16)
      integer at_nstep
      real*8   vel(0:nlevs-1,-2:nfine+1)
      real*8  scal(0:nlevs-1,-2:nfine+1,nscal)
      real*8 press(0:nlevs-1,-1:nfine+1)
      real*8   I_R(0:nlevs-1,-1:nfine  ,0:Nspec)
      real*8  divu(0:nlevs-1, 0:nfine-1)
      real*8  dsdt(0:nlevs-1, 0:nfine-1)
      real*8 time
      real*8 dt_old
      integer lo(0:nlevs-1), hi(0:nlevs-1)
      
      real*8 x
      integer i,l,n
      
      print *,'READING CHKFILE ',chkfile
      
      open(10,file=chkfile,form='unformatted',status='old')
      read(10) at_nstep
      read(10) time
      read(10) dt_old
      read(10) hmix_TYP
      read(10) vel_TYP

      do l=0,nlevs-1

!     cell-centered, 2 ghost cells
         do i=lo(l)-2,hi(l)+2
            read(10) x,(scal(l,i,FirstSpec+n),n=0,Nspec-1),
     $                  scal(l,i,Density),
     $                  scal(l,i,Temp),
     $                  scal(l,i,RhoH),
     $                  vel(l,i)
         enddo

!     cell-centered, 1 ghost cell
         do i=lo(l)-1,hi(l)+1
            read(10) x,(I_R(l,i,n),n=0,Nspec)
         enddo

!     cell-centered, no ghost cells
         do i=lo(l),hi(l)
            read(10) x,divu(l,i),dsdt(l,i)
         enddo

!     nodal, 1 ghost cell
         do i=lo(l)-1,hi(l)+2
            read(10) x,press(l,i)
         enddo

      enddo

      end
      
