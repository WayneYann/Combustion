! :::
! ::: ----------------------------------------------------------------
! :::

subroutine cns_umdrv(lo,hi,&
     uin  ,  uin_l1,  uin_l2,  uin_l3,  uin_h1,  uin_h2,  uin_h3, &
     uout , uout_l1, uout_l2, uout_l3, uout_h1, uout_h2, uout_h3, &
     delta,dt, &
     flux1,flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3, &
     flux2,flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3, &
     flux3,flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3, &
     area1,area1_l1,area1_l2,area1_l3,area1_h1,area1_h2,area1_h3, &
     area2,area2_l1,area2_l2,area2_l3,area2_h1,area2_h2,area2_h3, &
     area3,area3_l1,area3_l2,area3_l3,area3_h1,area3_h2,area3_h3, &
     vol  ,  vol_l1,  vol_l2,  vol_l3,  vol_h1,  vol_h2,  vol_h3, &
     courno,verbose)

  use meth_params_module, only : URHO, QVAR, NVAR, NHYP, NDIF, &
       QPRES,QU,QV,QW,normalize_species
  use advection_module, only : umeth3d, ctoprim, ptoderiv, srctosrcQ, uflaten, divu, consup, &
       enforce_minimum_density, normalize_new_species
  use diff_flux_module, only : diffFlux, fluxtosrc, diffup
  use chemsolv_module, only : chemsolv

  implicit none

  integer,intent(in):: lo(3),hi(3),verbose
  integer,intent(in)::   uin_l1,   uin_l2,   uin_l3,   uin_h1,   uin_h2,   uin_h3
  integer,intent(in)::  uout_l1,  uout_l2,  uout_l3,  uout_h1,  uout_h2,  uout_h3
  integer,intent(in):: flux1_l1, flux1_l2, flux1_l3, flux1_h1, flux1_h2, flux1_h3
  integer,intent(in):: flux2_l1, flux2_l2, flux2_l3, flux2_h1, flux2_h2, flux2_h3
  integer,intent(in):: flux3_l1, flux3_l2, flux3_l3, flux3_h1, flux3_h2, flux3_h3
  integer,intent(in):: area1_l1, area1_l2, area1_l3, area1_h1, area1_h2, area1_h3
  integer,intent(in):: area2_l1, area2_l2, area2_l3, area2_h1, area2_h2, area2_h3
  integer,intent(in):: area3_l1, area3_l2, area3_l3, area3_h1, area3_h2, area3_h3
  integer,intent(in)::   vol_l1,   vol_l2,   vol_l3,   vol_h1,   vol_h2,   vol_h3

  double precision,intent(inout)::  uin(  uin_l1:  uin_h1,  uin_l2:  uin_h2,  uin_l3:  uin_h3,NVAR)
  double precision,intent(  out):: uout( uout_l1: uout_h1, uout_l2: uout_h2, uout_l3: uout_h3,NVAR)
  double precision,intent(  out)::flux1(flux1_l1:flux1_h1,flux1_l2:flux1_h2,flux1_l3:flux1_h3,NVAR)
  double precision,intent(  out)::flux2(flux2_l1:flux2_h1,flux2_l2:flux2_h2,flux2_l3:flux2_h3,NVAR)
  double precision,intent(  out)::flux3(flux3_l1:flux3_h1,flux3_l2:flux3_h2,flux3_l3:flux3_h3,NVAR)
  double precision,intent(in   )::area1(area1_l1:area1_h1,area1_l2:area1_h2,area1_l3:area1_h3)
  double precision,intent(in   )::area2(area2_l1:area2_h1,area2_l2:area2_h2,area2_l3:area2_h3)
  double precision,intent(in   )::area3(area3_l1:area3_h1,area3_l2:area3_h2,area3_l3:area3_h3)
  double precision,intent(in   )::  vol(  vol_l1:  vol_h1,  vol_l2:  vol_h2,  vol_l3:  vol_h3)

  double precision, intent(in) :: delta(3),dt
  double precision, intent(inout) :: courno

  integer :: lo_work(3),hi_work(3)
  integer :: lo_diff(3),hi_diff(3)
  integer :: lo_hyp(3), hi_hyp(3)

  double precision, allocatable:: q(:,:,:,:)
  double precision, allocatable:: dpdr(:,:,:)
  double precision, allocatable:: dpde(:,:,:)
  double precision, allocatable:: gamc(:,:,:)
  double precision, allocatable:: flatn(:,:,:)
  double precision, allocatable:: c(:,:,:)
  double precision, allocatable:: csml(:,:,:)
  double precision, allocatable:: div(:,:,:)
  double precision, allocatable:: pdivu(:,:,:)
  double precision, allocatable:: src (:,:,:,:)
  double precision, allocatable:: srcQ(:,:,:,:)

  double precision, allocatable:: dfluxx(:,:,:,:)
  double precision, allocatable:: dfluxy(:,:,:,:)
  double precision, allocatable:: dfluxz(:,:,:,:)

  double precision :: dx,dy,dz, dt_flux, dt_src
  integer :: iflaten, dflux_timer

  allocate(    q(uin_l1:uin_h1,uin_l2:uin_h2,uin_l3:uin_h3,QVAR))
  allocate( dpdr(uin_l1:uin_h1,uin_l2:uin_h2,uin_l3:uin_h3))
  allocate( dpde(uin_l1:uin_h1,uin_l2:uin_h2,uin_l3:uin_h3))
  allocate( gamc(uin_l1:uin_h1,uin_l2:uin_h2,uin_l3:uin_h3))
  allocate(flatn(uin_l1:uin_h1,uin_l2:uin_h2,uin_l3:uin_h3))
  allocate(    c(uin_l1:uin_h1,uin_l2:uin_h2,uin_l3:uin_h3))
  allocate( csml(uin_l1:uin_h1,uin_l2:uin_h2,uin_l3:uin_h3))

  lo_hyp(1:3)= lo(1:3)-NDIF   ! hyperbolic solve on this domain
  hi_hyp(1:3)= hi(1:3)+NDIF  

  allocate(  div(lo_hyp(1):hi_hyp(1)+1,lo_hyp(2):hi_hyp(2)+1,lo_hyp(3):hi_hyp(3)+1))
  allocate(pdivu(lo_hyp(1):hi_hyp(1)  ,lo_hyp(2):hi_hyp(2)  ,lo_hyp(3):hi_hyp(3)  ))

  lo_diff(1:3)= lo(1:3)-NDIF-1 ! domain for the first diffusion solve
  hi_diff(1:3)= hi(1:3)+NDIF+1 

  allocate(dfluxx(lo_diff(1):hi_diff(1)+1,lo_diff(2):hi_diff(2)  ,lo_diff(3):hi_diff(3)  ,NVAR))
  allocate(dfluxy(lo_diff(1):hi_diff(1)  ,lo_diff(2):hi_diff(2)+1,lo_diff(3):hi_diff(3)  ,NVAR))
  allocate(dfluxz(lo_diff(1):hi_diff(1)  ,lo_diff(2):hi_diff(2)  ,lo_diff(3):hi_diff(3)+1,NVAR))

  allocate(src (lo_diff(1):hi_diff(1),lo_diff(2):hi_diff(2),lo_diff(3):hi_diff(3),NVAR))
  allocate(srcQ(lo_diff(1):hi_diff(1),lo_diff(2):hi_diff(2),lo_diff(3):hi_diff(3),QVAR))

  dx = delta(1)
  dy = delta(2)
  dz = delta(3)

  ! Chemically react input state for half time step
  lo_work = lo-NHYP-NDIF
  hi_work = hi+NHYP+NDIF
  call chemsolv(lo_work, hi_work, &
       uin,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, 0.5d0*dt)

  lo_work = lo_hyp - NHYP
  hi_work = hi_hyp + NHYP
  call ctoprim(lo_work, hi_work, &
       uin,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
       q  ,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3)

  call ptoderiv(lo_work, hi_work, q, c, gamc, csml, dpdr, dpde, &
       uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
       courno,dx,dy,dz,dt,lo,hi)

  iflaten = 1
  if (iflaten .eq. 1) then
     lo_work = lo_hyp - 1
     hi_work = hi_hyp + 1
     call uflaten(lo_work, hi_work, &
          q(:,:,:,QPRES), &
          q(:,:,:,QU), &
          q(:,:,:,QV), &
          q(:,:,:,QW), &
          flatn,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3)
  else
     flatn = 1.d0
  end if

  dflux_timer = 1
  call diffFlux(lo_diff,hi_diff, &
       q,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
       dfluxx,lo_diff(1),lo_diff(2),lo_diff(3),hi_diff(1)+1,hi_diff(2),hi_diff(3),&
       dfluxy,lo_diff(1),lo_diff(2),lo_diff(3),hi_diff(1),hi_diff(2)+1,hi_diff(3),&
       dfluxz,lo_diff(1),lo_diff(2),lo_diff(3),hi_diff(1),hi_diff(2),hi_diff(3)+1,&
       dx,dy,dz,dflux_timer)

  ! Compute source terms due to diffusion for hyperbolic advection step 
  call fluxtosrc(lo_diff, hi_diff, &
       src,lo_diff(1),lo_diff(2),lo_diff(3),hi_diff(1),hi_diff(2),hi_diff(3),&
       dfluxx,lo_diff(1),lo_diff(2),lo_diff(3),hi_diff(1)+1,hi_diff(2),hi_diff(3),&
       dfluxy,lo_diff(1),lo_diff(2),lo_diff(3),hi_diff(1),hi_diff(2)+1,hi_diff(3),&
       dfluxz,lo_diff(1),lo_diff(2),lo_diff(3),hi_diff(1),hi_diff(2),hi_diff(3)+1,&
       area1,area1_l1,area1_l2,area1_l3,area1_h1,area1_h2,area1_h3, &
       area2,area2_l1,area2_l2,area2_l3,area2_h1,area2_h2,area2_h3, &
       area3,area3_l1,area3_l2,area3_l3,area3_h1,area3_h2,area3_h3, &
       vol  ,  vol_l1,  vol_l2,  vol_l3,  vol_h1,  vol_h2,  vol_h3)

  call srctosrcQ(src,srcQ,lo_diff(1),lo_diff(2),lo_diff(3),hi_diff(1),hi_diff(2),hi_diff(3),&
       q,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
       dpdr,dpde,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3)
  
  ! Compute hyperbolic fluxes using unsplit Godunov
  call umeth3d(q,c,gamc,csml,flatn,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
       srcQ,lo_diff(1),lo_diff(2),lo_diff(3),hi_diff(1),hi_diff(2),hi_diff(3),&
       lo_hyp(1),lo_hyp(2),lo_hyp(3),hi_hyp(1),hi_hyp(2),hi_hyp(3),dx,dy,dz,dt, &
       flux1,flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3, &
       flux2,flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3, &
       flux3,flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3, &
       pdivu)
  
  ! Compute divergence of velocity field (on surroundingNodes(lo,hi))
  call divu(lo_hyp,hi_hyp,q,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
       dx,dy,dz,div,lo_hyp(1),lo_hyp(2),lo_hyp(3),hi_hyp(1)+1,hi_hyp(2)+1,hi_hyp(3)+1)
  
  ! Conservative update using hyperbolic flux (w diffusion source terms)
  call consup(uin,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
       uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3, &
       src,lo_diff(1),lo_diff(2),lo_diff(3),hi_diff(1),hi_diff(2),hi_diff(3),&
       flux1,flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3, &
       flux2,flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3, &
       flux3,flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3, &
       area1,area1_l1,area1_l2,area1_l3,area1_h1,area1_h2,area1_h3, &
       area2,area2_l1,area2_l2,area2_l3,area2_h1,area2_h2,area2_h3, &
       area3,area3_l1,area3_l2,area3_l3,area3_h1,area3_h2,area3_h3, &
       vol,vol_l1,vol_l2,vol_l3,vol_h1,vol_h2,vol_h3, &
       div,pdivu,lo_hyp,hi_hyp,dx,dy,dz,dt)

  ! Get primitives of the "guess" state
  call ctoprim(lo_hyp, hi_hyp, &
       uin,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
       q  ,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3)

  dflux_timer = 2
  call diffFlux(lo, hi, &
       q,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
       dfluxx,lo_diff(1),lo_diff(2),lo_diff(3),hi_diff(1)+1,hi_diff(2),hi_diff(3),&
       dfluxy,lo_diff(1),lo_diff(2),lo_diff(3),hi_diff(1),hi_diff(2)+1,hi_diff(3),&
       dfluxz,lo_diff(1),lo_diff(2),lo_diff(3),hi_diff(1),hi_diff(2),hi_diff(3)+1,&
       dx,dy,dz,dflux_timer)

  ! Add diffusion flux with half dt because the flux stores the sum of old and new fluxes,
  ! and subtract src (i.e., the old diffusion flux)
  dt_flux = 0.5d0*dt
  dt_src = -1.0d0*dt
  call diffup(lo, hi, &
       uout,uout_l1,uout_l2,uout_l3, uout_h1,uout_h2,uout_h3, &
       src,lo_diff(1),lo_diff(2),lo_diff(3),hi_diff(1),hi_diff(2),hi_diff(3),&
       dfluxx,lo_diff(1),lo_diff(2),lo_diff(3),hi_diff(1)+1,hi_diff(2),hi_diff(3),&
       dfluxy,lo_diff(1),lo_diff(2),lo_diff(3),hi_diff(1),hi_diff(2)+1,hi_diff(3),&
       dfluxz,lo_diff(1),lo_diff(2),lo_diff(3),hi_diff(1),hi_diff(2),hi_diff(3)+1,&
       area1,area1_l1,area1_l2,area1_l3,area1_h1,area1_h2,area1_h3, &
       area2,area2_l1,area2_l2,area2_l3,area2_h1,area2_h2,area2_h3, &
       area3,area3_l1,area3_l2,area3_l3,area3_h1,area3_h2,area3_h3, &
       vol  ,  vol_l1,  vol_l2,  vol_l3,  vol_h1,  vol_h2,  vol_h3, &
       dx,dy,dz,dt_flux, dt_src)

  ! Chemically react input state for half time step
  call chemsolv(lo, hi, &
       uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3, 0.5d0*dt)

  ! Add diffusion fluxes to hyperbolic fluxes to pass back to AMR
  flux1 = flux1 + dfluxx(flux1_l1:flux1_h1,flux1_l2:flux1_h2,flux1_l3:flux1_h3,:)
  flux2 = flux2 + dfluxy(flux2_l1:flux2_h1,flux2_l2:flux2_h2,flux2_l3:flux2_h3,:)
  flux3 = flux3 + dfluxz(flux3_l1:flux3_h1,flux3_l2:flux3_h2,flux3_l3:flux3_h3,:)
  
  ! Enforce the density >= small_dens.
  call enforce_minimum_density(uin, uin_l1, uin_l2, uin_l3, uin_h1, uin_h2, uin_h3, &
       uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3, &
       lo,hi,verbose)
  
  ! Enforce species >= 0
  call cns_enforce_nonnegative_species(uout,uout_l1,uout_l2,uout_l3, &
       uout_h1,uout_h2,uout_h3,lo,hi)
  
  ! Re-normalize the species
  if (normalize_species .eq. 1) then
     call normalize_new_species(uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3, &
          lo,hi)
  end if
  
  deallocate(q,dpdr,dpde,gamc,flatn,c,csml,div,src,srcQ,pdivu)
  deallocate(dfluxx,dfluxy,dfluxz)
  
end subroutine cns_umdrv

! :: ----------------------------------------------------------
! :: Volume-weight average the fine grid data onto the coarse
! :: grid.  Overlap is given in coarse grid coordinates.
! ::
! :: INPUTS / OUTPUTS:
! ::  crse      <=  coarse grid data
! ::  clo,chi    => index limits of crse array interior
! ::  nvar	 => number of components in arrays
! ::  fine       => fine grid data
! ::  flo,fhi    => index limits of fine array interior
! ::  lo,hi      => index limits of overlap (crse grid)
! ::  lrat       => refinement ratio
! ::
! :: NOTE:
! ::  Assumes all data cell centered
! :: ----------------------------------------------------------
! ::
subroutine cns_avgdown(crse,c_l1,c_l2,c_l3,c_h1,c_h2,c_h3,nvar, &
     cv,cv_l1,cv_l2,cv_l3,cv_h1,cv_h2,cv_h3, &
     fine,f_l1,f_l2,f_l3,f_h1,f_h2,f_h3, &
     fv,fv_l1,fv_l2,fv_l3,fv_h1,fv_h2,fv_h3,lo,hi,lrat)
  
  implicit none
  integer c_l1,c_l2,c_l3,c_h1,c_h2,c_h3
  integer cv_l1,cv_l2,cv_l3,cv_h1,cv_h2,cv_h3
  integer f_l1,f_l2,f_l3,f_h1,f_h2,f_h3
  integer fv_l1,fv_l2,fv_l3,fv_h1,fv_h2,fv_h3
  integer lo(3), hi(3)
  integer nvar, lrat(3)
  double precision crse(c_l1:c_h1,c_l2:c_h2,c_l3:c_h3,nvar)
  double precision cv(cv_l1:cv_h1,cv_l2:cv_h2,cv_l3:cv_h3)
  double precision fine(f_l1:f_h1,f_l2:f_h2,f_l3:f_h3,nvar)
  double precision fv(fv_l1:fv_h1,fv_l2:fv_h2,fv_l3:fv_h3)
  
  integer i, j, k, n, ic, jc, kc, ioff, joff, koff
  integer lratx, lraty, lratz
  double precision   volfrac
  
  lratx   = lrat(1)
  lraty   = lrat(2)
  lratz   = lrat(3)
  volfrac = 1.d0/dble(lrat(1)*lrat(2)*lrat(3))
  
  do n = 1, nvar
     !
     ! Set coarse grid to zero on overlap.
     !
     do kc = lo(3), hi(3)
        do jc = lo(2), hi(2)
           do ic = lo(1), hi(1)
              crse(ic,jc,kc,n) = 0.d0
           enddo
        enddo
     enddo
     !
     ! Sum fine data.
     !
     do koff = 0, lratz-1
        !$OMP PARALLEL DO PRIVATE(i,j,k,ic,jc,kc,ioff,joff)
        do kc = lo(3),hi(3)
           k = kc*lratz + koff
           do joff = 0, lraty-1
              do jc = lo(2), hi(2)
                 j = jc*lraty + joff
                 do ioff = 0, lratx-1
                    do ic = lo(1), hi(1)
                       i = ic*lratx + ioff
                       crse(ic,jc,kc,n) = crse(ic,jc,kc,n) + fine(i,j,k,n)
                    enddo
                 enddo
              enddo
           enddo
        enddo
        !$OMP END PARALLEL DO
     enddo
     !
     ! Divide out by volume weight.
     !
     !$OMP PARALLEL DO PRIVATE(ic,jc,kc)
     do kc = lo(3), hi(3)
        do jc = lo(2), hi(2)
           do ic = lo(1), hi(1)
              crse(ic,jc,kc,n) = volfrac*crse(ic,jc,kc,n)
           enddo
        enddo
     enddo
     !$OMP END PARALLEL DO
     
  enddo
  
end subroutine cns_avgdown

! :::
! ::: ------------------------------------------------------------------
! :::

subroutine cns_compute_temp(lo,hi,state,state_l1,state_l2,state_l3, &
                                 state_h1,state_h2,state_h3)

  use chemistry_module, only : nspecies
  use eos_module, only : eos_get_T
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEINT, UTEMP, &
       UFS, small_temp

  implicit none

  integer, intent(in) :: lo(3),hi(3)
  integer, intent(in) :: state_l1,state_l2,state_l3
  integer, intent(in) :: state_h1,state_h2,state_h3
  double precision, intent(inout) :: state(state_l1:state_h1, &
       &                                   state_l2:state_h2, &
       &                                   state_l3:state_h3,NVAR)

  integer          :: i,j,k
  double precision :: rhoInv,eint,xn(nspecies)
  integer          :: pt_index(3)

  do k = lo(3),hi(3)
     do j = lo(2),hi(2)
        do i = lo(1),hi(1)
           if (state(i,j,k,URHO) <= 0.d0) then
              print *,'   '
              print *,'>>> Error: CNSReact_3d::cns_compute_temp ',i,j,k
              print *,'>>> ... negative density ',state(i,j,k,URHO)
              print *,'    '
              call bl_error("Error:: CNSReact_3d.f90 :: cns_compute_temp")
           end if
        enddo
     enddo
  enddo
  
  !$OMP PARALLEL DO PRIVATE(i,j,k,rhoInv,xn,eint,pt_index)
  do k = lo(3),hi(3)
  do j = lo(2),hi(2)
  do i = lo(1),hi(1)

     rhoInv = 1.d0 / state(i,j,k,URHO)

     xn(1:nspecies)  = state(i,j,k,UFS:UFS+nspecies-1) * rhoInv

     eint = state(i,j,k,UEINT) * rhoInv

     pt_index(1) = i
     pt_index(2) = j
     pt_index(3) = k
     call eos_get_T(state(i,j,k,UTEMP), eint, xn, pt_index)

  enddo
  enddo
  enddo
  !$OMP END PARALLEL DO

end subroutine cns_compute_temp








