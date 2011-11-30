! ::: 
! ::: ----------------------------------------------------------------
! ::: 

      subroutine get_method_params(nGrowHyp)

        ! Passing data from f90 back to C++

        use meth_params_module

        implicit none 

        integer, intent(out) :: ngrowHyp

        nGrowHyp = NHYP

      end subroutine get_method_params

! ::: 
! ::: ----------------------------------------------------------------
! ::: 

      subroutine set_method_params(dm,Density,Xmom,Eden,Eint,Temp, &
                                   FirstAdv,FirstSpec,numadv, &
                                   small_dens_in, small_temp_in, small_pres_in, &
                                   allow_negative_energy_in,ppm_type_in, &
                                   normalize_species_in)

        ! Passing data from C++ into f90

        use meth_params_module
        use cdwrk_module, only : nspec
        use eos_module

        implicit none 
 
        integer, intent(in) :: dm
        integer, intent(in) :: Density, Xmom, Eden, Eint, Temp, FirstAdv, FirstSpec
        integer, intent(in) :: numadv
        integer, intent(in) :: allow_negative_energy_in, ppm_type_in
        double precision, intent(in) :: small_dens_in, small_temp_in, small_pres_in
        integer, intent(in) :: normalize_species_in

        integer             :: QLAST

        iorder = 2 
        difmag = 0.1d0

        ! NTHERM: number of thermodynamic variables
        ! NVAR  : number of total variables in initial system
        ! dm refers to mometum components, '4' refers to rho, rhoE, rhoe and T
        NTHERM = dm + 4
        NVAR = NTHERM + nspec +  numadv

        nadv = numadv

        ! We use these to index into the state "U"
        URHO  = Density   + 1
        UMX   = Xmom      + 1
        if (dm .ge. 2) UMY = UMX + 1
        if (dm .eq. 3) UMZ = UMY + 1
        UEDEN = Eden      + 1
        UEINT = Eint      + 1
        UTEMP = Temp      + 1

        if (numadv .ge. 1) then
          UFA   = FirstAdv  + 1
        else 
          UFA = 1
        end if

        UFS   = FirstSpec + 1

        ! QTHERM: number of primitive variables, which includes pressure (+1), but
        !         not little e (-1)
        ! QVAR  : number of total variables in primitive form

        QTHERM = NTHERM
        QVAR = QTHERM + nspec + numadv

        ! We use these to index into the state "Q"
        QRHO  = 1

        QU    = 2
        QLAST = 2

        if (dm .ge. 2) then
           QV    = 3
           QLAST = 3
        end if

        if (dm .eq. 3) then
           QW    = 4
           QLAST = 4
        end if

        QPRES   = QLAST + 1
        QREINT  = QLAST + 2
        QTEMP   = QTHERM 

        if (numadv .ge. 1) then
          QFA = QTHERM + 1
          QFS = QFA + numadv
        else 
          QFA = 1
          QFS = QTHERM + 1
        end if

        if (small_pres_in > 0.d0) then
          small_pres = small_pres_in
        else
          small_pres = 1.d-8
        end if

        call eos_init(small_dens=small_dens_in, small_temp=small_temp_in)

        call eos_get_small_dens(small_dens)
        call eos_get_small_temp(small_temp)

        allow_negative_energy = allow_negative_energy_in
        ppm_type               = ppm_type_in
        normalize_species     = normalize_species_in

      end subroutine set_method_params

! ::: 
! ::: ----------------------------------------------------------------
! ::: 

      subroutine set_problem_params(dm,physbc_lo_in,physbc_hi_in, phys_prob_lo_in,   &
         phys_prob_hi_in, Outflow_in,Symmetry_in,coord_type_in)

        ! Passing data from C++ into f90

        use prob_params_module

        implicit none 
 
        integer, intent(in) :: dm
        integer, intent(in) :: physbc_lo_in(dm),physbc_hi_in(dm)
        integer, intent(in) :: Outflow_in
        integer, intent(in) :: Symmetry_in
        integer, intent(in) :: coord_type_in
        double precision, intent(in) :: phys_prob_lo_in(dm),phys_prob_hi_in(dm)

        allocate(physbc_lo(dm))
        allocate(physbc_hi(dm))
        allocate(phys_prob_lo(dm))
        allocate(phys_prob_hi(dm))

        physbc_lo(:) = physbc_lo_in(:)
        physbc_hi(:) = physbc_hi_in(:)

        phys_prob_lo(:) = phys_prob_lo_in(:)
        phys_prob_hi(:) = phys_prob_hi_in(:)

        Outflow  = Outflow_in
        Symmetry = Symmetry_in

        coord_type = coord_type_in

      end subroutine set_problem_params

! ::: 
! ::: ----------------------------------------------------------------
! ::: 

      subroutine ca_set_special_tagging_flag(dummy,flag) 
      use probdata_module
      double precision :: dummy 
      integer          :: flag
      end subroutine ca_set_special_tagging_flag
