!
! Computes error for convergence study
!
! NOTE: this does not quite work for data from multilevel simulations
!  Would need to fix data_read() and fill_restart_data()
!  See boxlib/fabio.f90 fabio_ml_multifab_read_d() L583
!  can extract needed data and fill vars (ie. nlevs) from here
!

module convergence

  use bl_types
  use multifab_module
  use ml_layout_module
  use average
  use fabio_module

  implicit none
 
contains
 
  subroutine conv(loc,files)
    integer, intent(in   ) :: loc(:)
    integer, intent(in   ) :: files(:)

    !local variables   
    integer, parameter :: dm        = 2, &
                          levs      = 1, &
                          ! H chem vlaue
!                          nscal     = 19, &
                          ! CH4 2step value
                          nscal     = 16, &
                          ! H chem vlaue
!                          nspec     = 9, &
                          ! CH4 2step value
                          nspec     = 6, &
! indices of the data i care about
                          rho       = 1, &
                          rhoH      = 2, &
                          Temp      = 4, &
                          first_spec= 9, &
!
                          n_grids   = 3, &
                          n_err     = nspec+4, &
! not really concerned with ghost cells, throw away
                          ng_cell   = 0, &
                          exact_loc = 1024
! don't forget to check ir, dx

    integer                   :: nlevs,n,i,index
    integer                   :: ir(n_grids,dm)
    type(ml_layout)           :: mla(n_grids+1)
    real(dp_t)                :: time,dummy
    type(multifab)            :: s_exact,s(n_grids)
    type(multifab)            :: uold
    type(ml_boxarray)         :: mba
    logical                   :: nodal(dm)
    logical                   :: pmask(dm)
    type(multifab)            :: avg(n_grids)
    real(dp_t)                :: l2_error(n_grids,n_err), l1_error(n_grids,n_err)
    real(dp_t)                :: dx(n_grids,dm)
    real(dp_t)                :: dt(n_grids)
    character(len=20)         :: plot_names(nscal)
    character(len=20)         :: sd_name
    type(multifab), allocatable  :: plotdata(:)


    nlevs = 1
    nodal(:) = .true.
    pmask(:) = .true.
    
    ir(1,:) = 16
!    dx(1,:) = 1.d0/32.d0 !1.d0/dble(ir(1,:))
!    dt(1)   = 3.2d-5
    dx(1,:) = 1.d-3/16.d0 !1.d0/dble(ir(1,:))
    dt(1)   = 3.44d-7

    do n = 2, n_grids
       ir(n,:) = ir(n-1,:)/2.d0
       dx(n,:) = dx(n-1,:)/2.d0
       dt(n)   = dt(n-1)/2.d0
    end do

! could potentially get plot names from restart data
    plot_names(1) = 'density'
    plot_names(2) = 'rhoH'
    plot_names(3) = 'tracer'
    plot_names(4) = 'Temp'
    plot_names(5) = 'RhoRT'
    plot_names(6) = 'divu'
    plot_names(7) = 'dsdt'
    plot_names(8) = 'FuncCount'
!    plot_names(9) = 'X(H2)'
!    plot_names(10) = 'X(H)'
!    plot_names(11) = 'X(O)'
!    plot_names(12) = 'X(O2)'
!    plot_names(13) = 'X(OH)'
!    plot_names(14) = 'X(H2O)'
!    plot_names(15) = 'X(HO2)'
!    plot_names(16) = 'X(H2O2)'
!    plot_names(17) = 'X(N2)'
!    plot_names(18) = 'mag_vort'
!    plot_names(19) = 'HeatRelease'

    plot_names(9) = 'Y(02)'
    plot_names(10) = 'Y(H20)'
    plot_names(11) = 'Y(CH4)'
    plot_names(12) = 'Y(CO)'
    plot_names(13) = 'Y(CO2)'
    plot_names(14) = 'Y(N2)'
    plot_names(15) = 'mag_vort'
    plot_names(16) = 'HeatRelease'

    write(*,*) loc(1),loc(2),loc(3),loc(4)
    write(*,*) files(1),files(2)
! read  in data
    call initialize_from_restart(mla(n_grids+1),exact_loc,time,dummy,pmask,uold,&
                                 s_exact)
    call destroy(uold)

    do n = 1, n_grids
       call initialize_from_restart(mla(n),loc(n),time,dummy,pmask,uold,&
                                    s(n))
       call destroy(uold)
    enddo

! avg the exact data down onto course grids
    do n = 1,n_grids
       call multifab_build(avg(n),s(n)%la,nscal)
       call ml_cc_restriction(avg(n),s_exact,ir(n,:))
       call write_plotfile(n, mla(n),avg(n))
       call multifab_sub_sub(avg(n),s(n))
       call write_plotfile(n+n_grids, mla(n),avg(n))
    enddo

! calculate l1 error
    do i = 1,n_grids
       ! Density
       l1_error(i,1) = multifab_norm_l1_c(avg(i),rho,1,all=.false.)
       ! RhoH
       l1_error(i,2) = multifab_norm_l1_c(avg(i),rhoH,1,all=.false.)
       ! Temp
       l1_error(i,3) = multifab_norm_l1_c(avg(i),Temp,1,all=.false.)
       ! Species
       l1_error(i,n_err) = 0.d0
       do n = 0,nspec-1
          index = first_spec+n
          l1_error(i,4+n) = multifab_norm_l1_c(avg(i),index,1,all=.false.)
          l1_error(i,n_err) = l1_error(i,n_err) + l1_error(i,4+n)
       enddo
       l1_error(i,:) = l1_error(i,:)*(dx(i,1)*dx(i,2))
       l1_error(i,n_err) = l1_error(i,n_err)/dble(nspec)
    end do


! calculate l2 error
    do i = 1,n_grids
       ! Density
       l2_error(i,1) = multifab_norm_l2_c(avg(i),rho,1,all=.false.)
       ! RhoH
       l2_error(i,2) = multifab_norm_l2_c(avg(i),rhoH,1,all=.false.)
       ! Temp
       l2_error(i,3) = multifab_norm_l2_c(avg(i),Temp,1,all=.false.)
       l2_error(i,n_err) = 0.d0
       do n = 0,nspec-1
          index = first_spec+n
          l2_error(i,4+n) = multifab_norm_l2_c(avg(i),index,1,all=.false.)
          l2_error(i,n_err) = l2_error(i,n_err) + l2_error(i,4+n)
       enddo
       l2_error(i,:) = l2_error(i,:)*sqrt(dx(i,1)*dx(i,2))
       l2_error(i,n_err) = l2_error(i,n_err)/dble(nspec)
    end do

    write(files(1),*)'#Convergence study data -- L1 norm'
    write(files(1),*)'#   dt      Density      RhoH      Temp      Species      Combined'
    do n = 1, n_grids
       write(files(1),1000) dt(n),l1_error(n,1),l1_error(n,2),l1_error(n,3),&
                   l1_error(n,4),l1_error(n,5),l1_error(n,6),l1_error(n,7),&
                   l1_error(n,8),l1_error(n,9),l1_error(n,10)!,l1_error(n,11),&
!                   l1_error(n,12),l1_error(n,13)
    end do

    write(files(2),*)'#Convergence study data -- L2 norm'
    write(files(2),*)'#   dt      Density      RhoH      Temp      Species      Combined'
    do n = 1, n_grids
       write(files(2),1000) dt(n),l2_error(n,1),l2_error(n,2),l2_error(n,3),&
                   l2_error(n,4),l2_error(n,5),l2_error(n,6),l2_error(n,7),&
                   l2_error(n,8),l2_error(n,9),l2_error(n,10)!,l2_error(n,11),&
!                   l2_error(n,12),l2_error(n,13)
    end do

    do n = 1,n_grids+1
       call destroy(mla(n))
    end do
    
1000 FORMAT(14(E15.8,1X)) 

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
! could write own verion of mf_read to get more data out of the plotfile
    write(unit=sd_name, fmt='(a)') trim(dirname)
    call fabio_ml_multifab_read_d(mfs, sd_name)
    
  end subroutine data_read

  subroutine write_plotfile(istep_to_write, mla, mf)

    integer,         intent(in   ) :: istep_to_write
    type(ml_layout), intent(in   ) :: mla
    type(multifab),  intent(in   ) :: mf
  
    integer                        :: n,n_plot_comps

    allocate(plotdata(levs))
    n_plot_comps = nscal

    do n = 1,levs
       call multifab_build(plotdata(n), mla%la(n), n_plot_comps, 0)
       call multifab_copy_c(plotdata(n),1        ,mf,1,nscal)

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

end module convergence
