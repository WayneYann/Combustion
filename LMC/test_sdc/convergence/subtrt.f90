!
! Computes error for convergence study
!

module difference

  use bl_types
  use multifab_module
  use ml_layout_module
  use average
  use fabio_module

  implicit none
 
contains
 
  subroutine conv(sdc_loc,strang_loc,files,cfl)
    integer, intent(in   ) :: sdc_loc(:)
    integer, intent(in   ) :: strang_loc(:)
    integer, intent(in   ) :: files(:)
    real(dp_t), intent(in   ) :: cfl

    !local variables   
    integer, parameter :: dm        = 2, &
                          levs      = 1, &
                          nscal     = 4, &
                          nspec     = 3, &
                          n_grids   = 4, &
                          ng_cell   = 3, &
                          ng_grow   = 1
! dont'r forget to check  dx
    integer                   :: nlevs,n,i
    type(ml_layout)           :: mla(n_grids*2)
    real(dp_t)                :: time,dt
    type(multifab)            :: sdc(n_grids),strang(n_grids)
    type(multifab)            :: u_sdc,u_strang
    type(ml_boxarray)         :: mba
    logical                   :: nodal(dm)
    logical                   :: pmask(dm)
    real(dp_t)                :: l2_error(n_grids,nscal+dm+1), l1_error(n_grids,nscal+dm+1)
    real(dp_t)                :: dx(n_grids,dm)
    character(len=20)         :: plot_names(nscal+dm)
    character(len=20)         :: sd_name
    type(multifab), allocatable  :: plotdata(:)


    nlevs = 1
    nodal(:) = .true.
    pmask(:) = .true.
    
    dx(1,:) = 1./32.

    do n = 2, n_grids
       dx(n,:) = dx(n-1,:)/2.
    end do

    plot_names(1) = 'density'
    plot_names(2) = 'species A'
    plot_names(3) = 'species B'
    plot_names(4) = 'species C'
    plot_names(5) = 'x-velocity'
    plot_names(6) = 'y-velocity'



    do n = 1, n_grids
       ! read  in data
       call initialize_from_restart(mla(n),sdc_loc(n),time,dt,pmask,u_sdc,&
                                    sdc(n))

       i = n+n_grids
       call initialize_from_restart(mla(i),strang_loc(n),time,dt,pmask,u_strang,&
                                    strang(n))

       call multifab_sub_sub(sdc(n),strang(n))
       call multifab_sub_sub(u_sdc,u_strang)
       call write_plotfile(n, mla(n),sdc(n),u_sdc)

       ! calculate l1 error
       l1_error(n,nscal+1) = 0.d0
       do i = 1,nscal
          l1_error(n,i) = multifab_norm_l1_c(sdc(n),i,1,all=.false.)
          l1_error(n,nscal+1) = l1_error(n,nscal+1) + l1_error(n,i)
       enddo
       do i = 1,dm
          l1_error(n,nscal+1+i) = multifab_norm_l1_c(u_sdc,i,1,all=.false.)
       enddo

       l1_error(n,:) = l1_error(n,:)*(dx(n,1)*dx(n,2))
       l1_error(n,nscal+1) = l1_error(n,nscal+1)/dble(nscal)

       ! calculate l2 error
       l2_error(n,nscal+1) = 0.d0
       do i = 1,nscal
          l2_error(n,i) = multifab_norm_l2_c(sdc(n),i,1,all=.false.)
          l2_error(n,nscal+1) = l2_error(n,nscal+1) + l2_error(n,i)
       enddo
       do i = 1,dm
          l2_error(n,nscal+1+i) = multifab_norm_l2_c(u_sdc,i,1,all=.false.)
       enddo

       l2_error(n,:) = l2_error(n,:)*sqrt(dx(n,1)*dx(n,2))
       l2_error(n,nscal+1) = l2_error(n,nscal+1)/dble(nscal)

       call destroy(u_sdc)
       call destroy(u_strang)
    enddo

    write(files(1),*)'#Convergence study data -- L1 norm'
    write(files(1),*)'#    dt     SpeciesA     SpeciesB    SpeciesC    Combined'
    write(files(2),*)'#Convergence study data -- L2 norm'
    write(files(2),*)'#    dt     SpeciesA     SpeciesB    SpeciesC    Combined'

    do n = 1, n_grids
       write(files(1),1000) cfl*dx(n,1),l1_error(n,1),l1_error(n,2),l1_error(n,3),&
                   l1_error(n,4), l1_error(n,5),l1_error(n,6),l1_error(n,7), l1_error(n,8)
       write(files(2),1000) cfl*dx(n,1),l2_error(n,1),l2_error(n,2),l2_error(n,3),&
                   l2_error(n,4), l2_error(n,5),l2_error(n,6),l2_error(n,7),l2_error(n,8)
    end do

    do n = 1,n_grids*2
       call destroy(mla(n))
    end do
    
1000 FORMAT(11(E15.8,1X)) 

contains

   subroutine initialize_from_restart(mla,restart,time,dt,pmask,uold,sold)
 
     type(ml_layout),intent(out)   :: mla
     integer       , intent(in   ) :: restart
     real(dp_t)    , intent(  out) :: time,dt
     logical       , intent(in   ) :: pmask(:)
     type(multifab)       :: sold,uold

     type(ml_boxarray)         :: mba
     type(multifab), pointer   :: chkdata(:)

     integer :: n


     call fill_restart_data(restart,mba,chkdata)

     call ml_layout_build(mla,mba,pmask)

     nlevs = mba%nlevel

     call multifab_build(   uold, mla%la(1),    dm, ng_cell)
     call multifab_build(   sold, mla%la(1), nscal, ng_cell)

     call multifab_copy_c(uold,1,chkdata(1),1         ,dm)
     call multifab_copy_c(sold,1,chkdata(1),1+dm      ,nscal)
     !
     ! The layouts for chkdata and chk_p are built standalone, level
     ! by level, and need to be destroy()d as such as well.
     !
     do n = 1,nlevs     
        call destroy(chkdata(n)%la)
        call multifab_destroy(chkdata(n))
     end do
     deallocate(chkdata)
     call destroy(mba)

  end subroutine initialize_from_restart

  subroutine fill_restart_data(restart_int,mba,chkdata)

    integer          , intent(in   ) :: restart_int
    type(ml_boxarray), intent(  out) :: mba
    type(multifab)   , pointer       :: chkdata(:)

    character(len=7)                  :: sd_name
    integer                           :: n


    write(unit=sd_name,fmt='("plt",i4.4)') restart_int
    print *,'Reading ',sd_name,' to get state data for restart'
    call data_read(chkdata, sd_name)

    call build(mba,levs,dm)
    mba%pd(1) =  bbox(get_boxarray(chkdata(1)))
!     do n = 2,levs
!       mba%pd(n) = refine(mba%pd(n-1),2)
!       mba%rr(n-1,:) = rrs(n-1)
!     end do
    do n = 1,levs
      call boxarray_build_copy(mba%bas(n), get_boxarray(chkdata(n))) 
    end do

  end subroutine fill_restart_data

  subroutine data_read(mfs, dirname)

    use bl_IO_module
    use fab_module
    use fabio_module, only: fabio_ml_multifab_read_d
    use parallel

    type(multifab  ),                pointer :: mfs(:)
    character(len=*), intent(in   )          :: dirname

!    character(len=128) :: sd_name


!   Read the state data into a multilevel multifab.
    write(unit=sd_name, fmt='(a)') trim(dirname)
    call fabio_ml_multifab_read_d(mfs, sd_name)
    
  end subroutine data_read

  subroutine write_plotfile(istep_to_write, mla, mf, u)

    integer,         intent(in   ) :: istep_to_write
    type(ml_layout), intent(in   ) :: mla
    type(multifab),  intent(in   ) :: mf,u
  
    integer                        :: n,n_plot_comps

    allocate(plotdata(levs))
    n_plot_comps = nscal+dm

    do n = 1,levs
       call multifab_build(plotdata(n), mla%la(n), n_plot_comps, 0)
       call multifab_copy_c(plotdata(n),1         ,mf,1,nscal)
       call multifab_copy_c(plotdata(n),nscal+1   ,u ,1,dm)
    end do

    write(unit=sd_name,fmt='("plt",i4.4)') istep_to_write
    call fabio_ml_multifab_write_d(plotdata, mla%mba%rr(:,1), sd_name, plot_names, &
                                   mla%mba%pd(1), time, dx(1,:))

    do n = 1,levs
      call multifab_destroy(plotdata(n))
    end do
    deallocate(plotdata)

  end subroutine write_plotfile

end subroutine conv

end module difference
