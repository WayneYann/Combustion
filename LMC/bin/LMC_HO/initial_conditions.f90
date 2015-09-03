module initial_conditions_module
   use ghost_cells_module
   implicit none
   private
   include 'spec.h'
   public :: load_initial_conditions
   
contains
   
   subroutine load_initial_conditions(scal_avg, dx)
      double precision, intent(out) :: scal_avg(-2:nx+1, nscal)
      double precision, intent(in ) :: dx
      
      integer :: nsteps_ic
      integer :: nx_ic
      double precision :: time_ic
      double precision, allocatable :: data_ic(:,:)
      
      integer :: f, i, j, n, size
      
      open(10,file='ic.dat',form='formatted')
      read(10,*) nsteps_ic
      read(10,*) nx_ic
      read(10,*) time_ic
      
      size = (2*Nspec) + 8
      
      allocate(data_ic(0:nx_ic-1, 1:size))
      
      f = nx_ic/nx
      
      do i=0,nx_ic-1
        read(10,*) data_ic(i,1:size)
      end do
      
      scal_avg = 0
      do i=0,nx-1
         do j=i*f,(i+1)*f - 1
            scal_avg(i,Density) = scal_avg(i,Density) + data_ic(j,Nspec+2)
         end do
         scal_avg(i,Density) = scal_avg(i,Density)/dble(f)
         
         do n=1,Nspec
            do j=i*f,(i+1)*f - 1
               scal_avg(i,FirstSpec+n-1) = scal_avg(i,FirstSpec+n-1) &
                  + data_ic(j,n+1)
            end do
            scal_avg(i,FirstSpec+n-1) = scal_avg(i,FirstSpec+n-1)/dble(f)
            
            ! todo: is this the right thing to do?
            if (scal_avg(i,FirstSpec+n-1) .lt. 0) then
               scal_avg(i,FirstSpec+n-1) = 0
            end if
            
         end do
         
         do j=i*f,(i+1)*f - 1
            scal_avg(i,RhoH) = scal_avg(i,RhoH) + data_ic(j,Nspec+4)
         end do
         scal_avg(i,RhoH) = scal_avg(i,RhoH)/dble(f)
         
         do j=i*f,(i+1)*f - 1
            scal_avg(i,Temp) = scal_avg(i,Temp) + data_ic(j,Nspec+3)
         end do
         scal_avg(i,Temp) = scal_avg(i,Temp)/dble(f)
         
         if (i .eq. 0) then
            hmix_TYP = ABS(scal_avg(i,RhoH)/scal_avg(i,Density))
         else
            hmix_TYP = MAX(hmix_TYP,ABS(scal_avg(i,RhoH)/scal_avg(i,Density)))
         endif
      end do
      
      c_0 = 0.d0
      c_1 = 0.d0
      
      deallocate(data_ic)
      close(10)
      
      scal_avg(:,RhoRT) = 0
      call fill_scal_avg_ghost_cells(scal_avg)
   end subroutine load_initial_conditions
   
end module initial_conditions_module
