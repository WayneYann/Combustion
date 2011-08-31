module ppm_module

  implicit none

  public

  private :: ppm_type1, ppm_type2

contains
  !
  ! characteristics based on u
  !
  subroutine ppm(s,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3,rho,u,cspd,Ip,Im, &
                 grav,gv_l1,gv_l2,gv_l3,gv_h1,gv_h2,gv_h3, &
                 ilo1,ilo2,ihi1,ihi2,dx,dy,dz,dt,k3d,kc,ivar)

    use meth_params_module, only : ppm_type

    implicit none

    integer          qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
    integer          gv_l1,gv_l2,gv_l3,gv_h1,gv_h2,gv_h3
    integer          ilo1,ilo2,ihi1,ihi2

    double precision    s(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
    double precision  rho(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
    double precision    u(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,3)
    double precision cspd(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)

    double precision Ip(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3)
    double precision Im(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3)

    double precision grav(gv_l1:gv_h1,gv_l2:gv_h2,gv_l3:gv_h3,3)

    double precision dx,dy,dz,dt
    integer          k3d,kc,ivar

    if (ppm_type .eq. 1) then 

       call ppm_type1(s,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3,rho,u,cspd,Ip,Im, &
                      grav,gv_l1,gv_l2,gv_l3,gv_h1,gv_h2,gv_h3, &
                      ilo1,ilo2,ihi1,ihi2,dx,dy,dz,dt,k3d,kc,ivar)

    else if (ppm_type .eq. 2) then 

       call ppm_type2(s,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3,u,cspd,Ip,Im, &
                      ilo1,ilo2,ihi1,ihi2,dx,dy,dz,dt,k3d,kc,ivar)

    end if

  end subroutine ppm

  ! :::
  ! ::: ----------------------------------------------------------------
  ! :::
  
  subroutine ppm_type1(s,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3,rho,u,cspd,Ip,Im, &
                       grav,gv_l1,gv_l2,gv_l3,gv_h1,gv_h2,gv_h3, &
                       ilo1,ilo2,ihi1,ihi2,dx,dy,dz,dt,k3d,kc,ivar)

    use meth_params_module, only : ppm_type, QPRES

    implicit none

    integer          qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
    integer          gv_l1,gv_l2,gv_l3,gv_h1,gv_h2,gv_h3
    integer          ilo1,ilo2,ihi1,ihi2

    double precision    s(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
    double precision  rho(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
    double precision    u(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,3)
    double precision cspd(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)

    double precision Ip(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3)
    double precision Im(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3)

    double precision grav(gv_l1:gv_h1,gv_l2:gv_h2,gv_l3:gv_h3,3)

    double precision dx,dy,dz,dt
    integer          k3d,kc,ivar

    ! local
    integer i,j,k

    double precision dsl, dsr, dsc
    double precision sgn, sigma, s6

    ! s_{\ib,+}, s_{\ib,-}
    double precision, allocatable :: sp(:,:)
    double precision, allocatable :: sm(:,:)

    ! \delta s_{\ib}^{vL}
    double precision, allocatable :: dsvl(:,:)
    double precision, allocatable :: dsvlm(:,:)
    double precision, allocatable :: dsvlp(:,:)

    ! s_{i+\half}^{H.O.}
    double precision, allocatable :: sedge(:,:)
    double precision, allocatable :: sedgez(:,:,:)

    ! cell-centered indexing
    allocate(sp(ilo1-1:ihi1+1,ilo2-1:ihi2+1))
    allocate(sm(ilo1-1:ihi1+1,ilo2-1:ihi2+1))

    if (ppm_type .ne. 1) &
         call bl_error("Should have ppm_type = 1 in ppm_type1")

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! x-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! cell-centered indexing w/extra x-ghost cell
    allocate(dsvl(ilo1-2:ihi1+2,ilo2-1:ihi2+1))

    ! edge-centered indexing for x-faces -- ppm_type = 1 only
    allocate(sedge(ilo1-1:ihi1+2,ilo2-1:ihi2+1))

    ! compute s at x-edges

    ! compute van Leer slopes in x-direction
    dsvl = 0.d0
    do j=ilo2-1,ihi2+1
       do i=ilo1-2,ihi1+2
          dsc = 0.5d0 * (s(i+1,j,k3d) - s(i-1,j,k3d))
          dsl = 2.d0  * (s(i  ,j,k3d) - s(i-1,j,k3d))
          dsr = 2.d0  * (s(i+1,j,k3d) - s(i  ,j,k3d))
          if (dsl*dsr .gt. 0.d0) &
               dsvl(i,j) = sign(1.d0,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
          !                if (ivar.eq.QPRES) then
          !                   if ( (grav(i  ,j,k3d,1)*grav(i-1,j,k3d,1) .lt. 0.d0 ) ) then
          !                       print *,'CHANGE OF SIGN AT ',i,j,k3d
          !                       print *,'GRAV:LEFT MID ',grav(i-1,j,k3d,1), grav(i,j,k3d,1)
          !                       print *,'PRES:LT MD RT ',s(i-1,j,k3d), s(i,j,k3d), s(i+1,j,k3d) 
          !                       print *,'DSL  ',dsl
          !                       print *,'DSR  ',dsr
          !                       print *,'DSC  ',dsc
          !                       print *,'DSVL ',dsvl(i,j)
          !                       print *,' '
          !                   end if
          !                 end if
       end do
    end do

    ! interpolate s to x-edges
    do j=ilo2-1,ihi2+1
       do i=ilo1-1,ihi1+2
          sedge(i,j) = 0.5d0*(s(i,j,k3d)+s(i-1,j,k3d)) &
               - (1.d0/6.d0)*(dsvl(i,j)-dsvl(i-1,j))
          ! make sure sedge lies in between adjacent cell-centered values
          sedge(i,j) = max(sedge(i,j),min(s(i,j,k3d),s(i-1,j,k3d)))
          sedge(i,j) = min(sedge(i,j),max(s(i,j,k3d),s(i-1,j,k3d)))
       end do
    end do

    !$OMP PARALLEL DO PRIVATE(i,j,s6,sigma)
    do j=ilo2-1,ihi2+1
       do i=ilo1-1,ihi1+1

          ! copy sedge into sp and sm
          sp(i,j) = sedge(i+1,j)
          sm(i,j) = sedge(i  ,j)

          ! modify using quadratic limiters
          if ((sp(i,j)-s(i,j,k3d))*(s(i,j,k3d)-sm(i,j)) .le. 0.d0) then
             sp(i,j) = s(i,j,k3d)
             sm(i,j) = s(i,j,k3d)
          else if (abs(sp(i,j)-s(i,j,k3d)) .ge. 2.d0*abs(sm(i,j)-s(i,j,k3d))) then
             sp(i,j) = 3.d0*s(i,j,k3d) - 2.d0*sm(i,j)
          else if (abs(sm(i,j)-s(i,j,k3d)) .ge. 2.d0*abs(sp(i,j)-s(i,j,k3d))) then
             sm(i,j) = 3.d0*s(i,j,k3d) - 2.d0*sp(i,j)
          end if

          ! compute x-component of Ip and Im
          s6 = 6.0d0*s(i,j,k3d) - 3.0d0*(sm(i,j)+sp(i,j))
          sigma = abs(u(i,j,k3d,1)-cspd(i,j,k3d))*dt/dx
          Ip(i,j,kc,1,1) = sp(i,j) - &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)-(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
          Im(i,j,kc,1,1) = sm(i,j) + &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)+(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
          sigma = abs(u(i,j,k3d,1))*dt/dx
          Ip(i,j,kc,1,2) = sp(i,j) - &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)-(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
          Im(i,j,kc,1,2) = sm(i,j) + &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)+(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
          sigma = abs(u(i,j,k3d,1)+cspd(i,j,k3d))*dt/dx
          Ip(i,j,kc,1,3) = sp(i,j) - &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)-(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
          Im(i,j,kc,1,3) = sm(i,j) + &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)+(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
       end do
    end do
    !$OMP END PARALLEL DO

    deallocate(sedge,dsvl)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! y-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! cell-centered indexing w/extra y-ghost cell
    allocate( dsvl(ilo1-1:ihi1+1,ilo2-2:ihi2+2))

    ! edge-centered indexing for y-faces
    allocate(sedge(ilo1-1:ihi1+1,ilo2-1:ihi2+2))

    ! compute s at y-edges

    ! compute van Leer slopes in y-direction
    dsvl = 0.d0
    do j=ilo2-2,ihi2+2
       do i=ilo1-1,ihi1+1
          dsc = 0.5d0 * (s(i,j+1,k3d) - s(i,j-1,k3d))
          dsl = 2.d0  * (s(i,j  ,k3d) - s(i,j-1,k3d))
          dsr = 2.d0  * (s(i,j+1,k3d) - s(i,j  ,k3d))
          if (dsl*dsr .gt. 0.d0) &
               dsvl(i,j) = sign(1.d0,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
       end do
    end do

    ! interpolate s to y-edges
    do j=ilo2-1,ihi2+2
       do i=ilo1-1,ihi1+1
          sedge(i,j) = 0.5d0*(s(i,j,k3d)+s(i,j-1,k3d)) &
               - (1.d0/6.d0)*(dsvl(i,j)-dsvl(i,j-1))
          ! make sure sedge lies in between adjacent cell-centered values
          sedge(i,j) = max(sedge(i,j),min(s(i,j,k3d),s(i,j-1,k3d)))
          sedge(i,j) = min(sedge(i,j),max(s(i,j,k3d),s(i,j-1,k3d)))
       end do
    end do

    !$OMP PARALLEL DO PRIVATE(i,j,s6,sigma)
    do j=ilo2-1,ihi2+1
       do i=ilo1-1,ihi1+1

          ! copy sedge into sp and sm
          sp(i,j) = sedge(i,j+1)
          sm(i,j) = sedge(i,j  )

          ! modify using quadratic limiters
          if ((sp(i,j)-s(i,j,k3d))*(s(i,j,k3d)-sm(i,j)) .le. 0.d0) then
             sp(i,j) = s(i,j,k3d)
             sm(i,j) = s(i,j,k3d)
          else if (abs(sp(i,j)-s(i,j,k3d)) .ge. 2.d0*abs(sm(i,j)-s(i,j,k3d))) then
             sp(i,j) = 3.d0*s(i,j,k3d) - 2.d0*sm(i,j)
          else if (abs(sm(i,j)-s(i,j,k3d)) .ge. 2.d0*abs(sp(i,j)-s(i,j,k3d))) then
             sm(i,j) = 3.d0*s(i,j,k3d) - 2.d0*sp(i,j)
          end if

          ! compute y-component of Ip and Im
          s6 = 6.0d0*s(i,j,k3d) - 3.0d0*(sm(i,j)+sp(i,j))
          sigma = abs(u(i,j,k3d,2)-cspd(i,j,k3d))*dt/dy
          Ip(i,j,kc,2,1) = sp(i,j) - &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)-(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
          Im(i,j,kc,2,1) = sm(i,j) + &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)+(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
          sigma = abs(u(i,j,k3d,2))*dt/dy
          Ip(i,j,kc,2,2) = sp(i,j) - &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)-(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
          Im(i,j,kc,2,2) = sm(i,j) + &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)+(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
          sigma = abs(u(i,j,k3d,2)+cspd(i,j,k3d))*dt/dy
          Ip(i,j,kc,2,3) = sp(i,j) - &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)-(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
          Im(i,j,kc,2,3) = sm(i,j) + &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)+(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
       end do
    end do
    !$OMP END PARALLEL DO

    deallocate(dsvl,sedge)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! z-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! cell-centered indexing
    allocate( dsvl(ilo1-1:ihi1+1,ilo2-1:ihi2+1))
    allocate(dsvlm(ilo1-1:ihi1+1,ilo2-1:ihi2+1))
    allocate(dsvlp(ilo1-1:ihi1+1,ilo2-1:ihi2+1))

    allocate(sedgez(ilo1-1:ihi1+1,ilo2-2:ihi2+3,k3d-1:k3d+2))

    ! compute s at z-edges

    ! compute van Leer slopes in z-direction
    dsvl  = 0.d0
    dsvlm = 0.d0
    dsvlp = 0.d0

    !$OMP PARALLEL DO PRIVATE(i,j,k,dsc,dsl,dsr,s6,sigma)
    do j=ilo2-1,ihi2+1
       do i=ilo1-1,ihi1+1

          ! compute on slab below
          k = k3d-1
          dsc = 0.5d0 * (s(i,j,k+1) - s(i,j,k-1))
          dsl = 2.0d0 * (s(i,j,k  ) - s(i,j,k-1))
          dsr = 2.0d0 * (s(i,j,k+1) - s(i,j,k  ))
          if (dsl*dsr .gt. 0.d0) &
               dsvlm(i,j) = sign(1.0d0,dsc)*min(abs(dsc),abs(dsl),abs(dsr))

          ! compute on slab above
          k = k3d+1
          dsc = 0.5d0 * (s(i,j,k+1) - s(i,j,k-1))
          dsl = 2.0d0 * (s(i,j,k  ) - s(i,j,k-1))
          dsr = 2.0d0 * (s(i,j,k+1) - s(i,j,k  ))
          if (dsl*dsr .gt. 0.d0) &
               dsvlp(i,j) = sign(1.0d0,dsc)*min(abs(dsc),abs(dsl),abs(dsr))

          ! compute on current slab
          k = k3d
          dsc = 0.5d0 * (s(i,j,k+1) - s(i,j,k-1))
          dsl = 2.0d0 * (s(i,j,k  ) - s(i,j,k-1))
          dsr = 2.0d0 * (s(i,j,k+1) - s(i,j,k  ))
          if (dsl*dsr .gt. 0.d0) &
               dsvl(i,j) = sign(1.0d0,dsc)*min(abs(dsc),abs(dsl),abs(dsr))

          ! interpolate to lo face
          k = k3d
          sm(i,j) = 0.5d0*(s(i,j,k)+s(i,j,k-1)) - (1.0d0/6.0d0)*(dsvl(i,j)-dsvlm(i,j))
          ! make sure sedge lies in between adjacent cell-centered values
          sm(i,j) = max(sm(i,j),min(s(i,j,k),s(i,j,k-1)))
          sm(i,j) = min(sm(i,j),max(s(i,j,k),s(i,j,k-1)))

          ! interpolate to hi face
          k = k3d+1
          sp(i,j) = 0.5d0*(s(i,j,k)+s(i,j,k-1)) - (1.0d0/6.0d0)*(dsvlp(i,j)-dsvl(i,j))
          ! make sure sedge lies in between adjacent cell-centered values
          sp(i,j) = max(sp(i,j),min(s(i,j,k),s(i,j,k-1)))
          sp(i,j) = min(sp(i,j),max(s(i,j,k),s(i,j,k-1)))

          ! modify using quadratic limiters
          if ((sp(i,j)-s(i,j,k3d))*(s(i,j,k3d)-sm(i,j)) .le. 0.0d0) then
             sp(i,j) = s(i,j,k3d)
             sm(i,j) = s(i,j,k3d)
          else if (abs(sp(i,j)-s(i,j,k3d)) .ge. 2.0d0*abs(sm(i,j)-s(i,j,k3d))) then
             sp(i,j) = 3.0d0*s(i,j,k3d) - 2.0d0*sm(i,j)
          else if (abs(sm(i,j)-s(i,j,k3d)) .ge. 2.0d0*abs(sp(i,j)-s(i,j,k3d))) then
             sm(i,j) = 3.0d0*s(i,j,k3d) - 2.0d0*sp(i,j)
          end if

          ! compute z-component of Ip and Im
          s6 = 6.0d0*s(i,j,k3d) - 3.0d0*(sm(i,j)+sp(i,j))
          sigma = abs(u(i,j,k3d,3)-cspd(i,j,k3d))*dt/dz
          Ip(i,j,kc,3,1) = sp(i,j) - &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)-(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
          Im(i,j,kc,3,1) = sm(i,j) + &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)+(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
          sigma = abs(u(i,j,k3d,3))*dt/dz
          Ip(i,j,kc,3,2) = sp(i,j) - &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)-(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
          Im(i,j,kc,3,2) = sm(i,j) + &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)+(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
          sigma = abs(u(i,j,k3d,3)+cspd(i,j,k3d))*dt/dz
          Ip(i,j,kc,3,3) = sp(i,j) - &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)-(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
          Im(i,j,kc,3,3) = sm(i,j) + &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)+(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
       end do
    end do
    !$OMP END PARALLEL DO

    deallocate(dsvl,dsvlm,dsvlp,sp,sm,sedgez)

  end subroutine ppm_type1

  ! :::
  ! ::: ----------------------------------------------------------------
  ! :::

  subroutine ppm_type2(s,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3,u,cspd,Ip,Im, &
                       ilo1,ilo2,ihi1,ihi2,dx,dy,dz,dt,k3d,kc,ivar)

    use meth_params_module, only : ppm_type, QPRES

    implicit none

    integer          qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
    integer          ilo1,ilo2,ihi1,ihi2
    double precision s(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
    double precision u(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,1:3)
    double precision cspd(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
    double precision Ip(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3)
    double precision Im(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3)
    double precision dx,dy,dz,dt
    integer          k3d,kc,ivar

    ! local
    integer i,j,k
    logical extremum, bigp, bigm

    double precision D2, D2C, D2L, D2R, D2LIM, alphap, alpham
    double precision sgn, sigma, s6
    double precision dafacem, dafacep, dabarm, dabarp, dafacemin, dabarmin
    double precision dachkm, dachkp
    double precision amax, delam, delap

    ! s_{\ib,+}, s_{\ib,-}
    double precision, allocatable :: sp(:,:)
    double precision, allocatable :: sm(:,:)

    ! \delta s_{\ib}^{vL}
    double precision, allocatable :: dsvl(:,:)
    double precision, allocatable :: dsvlm(:,:)
    double precision, allocatable :: dsvlp(:,:)

    ! s_{i+\half}^{H.O.}
    double precision, allocatable :: sedge(:,:)
    double precision, allocatable :: sedgez(:,:,:)

    ! constant used in Colella 2008
    double precision, parameter :: C = 1.25d0

    ! cell-centered indexing
    allocate(sp(ilo1-1:ihi1+1,ilo2-1:ihi2+1))
    allocate(sm(ilo1-1:ihi1+1,ilo2-1:ihi2+1))

    if (ppm_type .ne. 2) &
         call bl_error("Should have ppm_type = 2 in ppm_type2")

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! x-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! cell-centered indexing w/extra x-ghost cell
    allocate(dsvl(ilo1-2:ihi1+2,ilo2-1:ihi2+1))

    ! edge-centered indexing for x-faces
    allocate(sedge(ilo1-2:ihi1+3,ilo2-1:ihi2+1))

    ! compute s at x-edges

    ! interpolate s to x-edges
    do j=ilo2-1,ihi2+1
       do i=ilo1-2,ihi1+3
          sedge(i,j) = (7.d0/12.d0)*(s(i-1,j,k3d)+s(i  ,j,k3d)) &
               - (1.d0/12.d0)*(s(i-2,j,k3d)+s(i+1,j,k3d))
          !
          ! limit sedge
          !
          if ((sedge(i,j)-s(i-1,j,k3d))*(s(i,j,k3d)-sedge(i,j)) .lt. 0.d0) then
             D2  = 3.d0*(s(i-1,j,k3d)-2.d0*sedge(i,j)+s(i,j,k3d))
             D2L = s(i-2,j,k3d)-2.d0*s(i-1,j,k3d)+s(i,j,k3d)
             D2R = s(i-1,j,k3d)-2.d0*s(i,j,k3d)+s(i+1,j,k3d)
             sgn = sign(1.d0,D2)
             D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),0.d0)
             sedge(i,j) = 0.5d0*(s(i-1,j,k3d)+s(i,j,k3d)) - (1.d0/6.d0)*D2LIM
          end if
       end do
    end do
    !
    ! Use Colella 2008 limiters.
    !
    ! This is a new version of the algorithm to eliminate sensitivity to roundoff.
    !
    !$OMP PARALLEL DO PRIVATE(i,j,alphap,alpham,bigp,bigm,extremum,dafacem) &
    !$OMP PRIVATE(dafacep,dabarm,dabarp,dafacemin,dabarmin,dachkm,dachkp,D2,D2L) &
    !$OMP PRIVATE(D2R,D2C,sgn,D2LIM,amax,delam,delap,s6,sigma)
    do j=ilo2-1,ihi2+1
       do i=ilo1-1,ihi1+1

          alphap   = sedge(i+1,j)-s(i,j,k3d)
          alpham   = sedge(i  ,j)-s(i,j,k3d)
          bigp     = abs(alphap).gt.2.d0*abs(alpham)
          bigm     = abs(alpham).gt.2.d0*abs(alphap)
          extremum = .false.

          if (alpham*alphap .ge. 0.d0) then
             extremum = .true.
          else if (bigp .or. bigm) then
             !
             ! Possible extremum. We look at cell centered values and face
             ! centered values for a change in sign in the differences adjacent to
             ! the cell. We use the pair of differences whose minimum magnitude is the
             ! largest, and thus least susceptible to sensitivity to roundoff.
             !
             dafacem   = sedge(i,j) - sedge(i-1,j)
             dafacep   = sedge(i+2,j) - sedge(i+1,j)
             dabarm    = s(i,j,k3d) - s(i-1,j,k3d)
             dabarp    = s(i+1,j,k3d) - s(i,j,k3d)
             dafacemin = min(abs(dafacem),abs(dafacep))
             dabarmin  = min(abs(dabarm),abs(dabarp))
             if (dafacemin.ge.dabarmin) then
                dachkm = dafacem
                dachkp = dafacep
             else
                dachkm = dabarm
                dachkp = dabarp
             endif
             extremum = (dachkm*dachkp .le. 0.d0)
          end if

          if (extremum) then
             D2     = 6.d0*(alpham + alphap)
             D2L    = s(i-2,j,k3d)-2.d0*s(i-1,j,k3d)+s(i,j,k3d)
             D2R    = s(i,j,k3d)-2.d0*s(i+1,j,k3d)+s(i+2,j,k3d)
             D2C    = s(i-1,j,k3d)-2.d0*s(i,j,k3d)+s(i+1,j,k3d)
             sgn    = sign(1.d0,D2)
             D2LIM  = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),0.d0)
             alpham = alpham*D2LIM/max(abs(D2),1.d-10)
             alphap = alphap*D2LIM/max(abs(D2),1.d-10)
          else
             if (bigp) then
                sgn   = sign(1.d0,alpham)
                amax  = -alphap**2 / (4*(alpham + alphap))
                delam = s(i-1,j,k3d) - s(i,j,k3d)
                if (sgn*amax .ge. sgn*delam) then
                   if (sgn*(delam - alpham).ge.1.d-10) then
                      alphap = (-2.d0*delam - 2.d0*sgn*sqrt(delam**2 - delam*alpham))
                   else 
                      alphap = -2.d0*alpham
                   endif
                endif
             end if
             if (bigm) then
                sgn   = sign(1.d0,alphap)
                amax  = -alpham**2 / (4*(alpham + alphap))
                delap = s(i+1,j,k3d) - s(i,j,k3d)
                if (sgn*amax .ge. sgn*delap) then
                   if (sgn*(delap - alphap).ge.1.d-10) then
                      alpham = (-2.d0*delap - 2.d0*sgn*sqrt(delap**2 - delap*alphap))
                   else
                      alpham = -2.d0*alphap
                   endif
                endif
             end if
          end if

          sm(i,j) = s(i,j,k3d) + alpham
          sp(i,j) = s(i,j,k3d) + alphap
          !
          ! Compute x-component of Ip and Im.
          !
          s6    = 6.0d0*s(i,j,k3d) - 3.0d0*(sm(i,j)+sp(i,j))
          sigma = abs(u(i,j,k3d,1)-cspd(i,j,k3d))*dt/dx

          Ip(i,j,kc,1,1) = sp(i,j) - &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)-(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
          Im(i,j,kc,1,1) = sm(i,j) + &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)+(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
          sigma = abs(u(i,j,k3d,1))*dt/dx
          Ip(i,j,kc,1,2) = sp(i,j) - &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)-(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
          Im(i,j,kc,1,2) = sm(i,j) + &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)+(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
          sigma = abs(u(i,j,k3d,1)+cspd(i,j,k3d))*dt/dx
          Ip(i,j,kc,1,3) = sp(i,j) - &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)-(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
          Im(i,j,kc,1,3) = sm(i,j) + &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)+(1.0d0-(2.0d0/3.0d0)*sigma)*s6)

       end do
    end do
    !$OMP END PARALLEL DO

    deallocate(sedge,dsvl)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! y-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! cell-centered indexing w/extra y-ghost cell
    allocate( dsvl(ilo1-1:ihi1+1,ilo2-2:ihi2+2))

    ! edge-centered indexing for y-faces
    allocate(sedge(ilo1-1:ihi1+1,ilo2-2:ihi2+3))

    ! compute s at y-edges

    ! interpolate s to y-edges
    do j=ilo2-2,ihi2+3
       do i=ilo1-1,ihi1+1
          sedge(i,j) = (7.d0/12.d0)*(s(i,j-1,k3d)+s(i,j,k3d)) &
               - (1.d0/12.d0)*(s(i,j-2,k3d)+s(i,j+1,k3d))
          !
          ! limit sedge
          !
          if ((sedge(i,j)-s(i,j-1,k3d))*(s(i,j,k3d)-sedge(i,j)) .lt. 0.d0) then
             D2  = 3.d0*(s(i,j-1,k3d)-2.d0*sedge(i,j)+s(i,j,k3d))
             D2L = s(i,j-2,k3d)-2.d0*s(i,j-1,k3d)+s(i,j,k3d)
             D2R = s(i,j-1,k3d)-2.d0*s(i,j,k3d)+s(i,j+1,k3d)
             sgn = sign(1.d0,D2)
             D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),0.d0)
             sedge(i,j) = 0.5d0*(s(i,j-1,k3d)+s(i,j,k3d)) - (1.d0/6.d0)*D2LIM
          end if
       end do
    end do
    !
    ! Use Colella 2008 limiters.
    !
    ! This is a new version of the algorithm to eliminate sensitivity to roundoff.
    !
    !$OMP PARALLEL DO PRIVATE(i,j,alphap,alpham,bigp,bigm,extremum,dafacem,dafacep,dabarm,dabarp,dafacemin) &
    !$OMP PRIVATE(dabarmin,dachkm,dachkp,D2,D2L,D2R,D2C,sgn,D2LIM,amax,delam,delap,s6,sigma)
    do j=ilo2-1,ihi2+1
       do i=ilo1-1,ihi1+1

          alphap   = sedge(i,j+1)-s(i,j,k3d)
          alpham   = sedge(i,j  )-s(i,j,k3d)
          bigp     = abs(alphap).gt.2.d0*abs(alpham)
          bigm     = abs(alpham).gt.2.d0*abs(alphap)
          extremum = .false.

          if (alpham*alphap .ge. 0.d0) then
             extremum = .true.
          else if (bigp .or. bigm) then
             !
             ! Possible extremum. We look at cell centered values and face
             ! centered values for a change in sign in the differences adjacent to
             ! the cell. We use the pair of differences whose minimum magnitude is the
             ! largest, and thus least susceptible to sensitivity to roundoff.
             !
             dafacem   = sedge(i,j) - sedge(i,j-1)
             dafacep   = sedge(i,j+2) - sedge(i,j+1)
             dabarm    = s(i,j,k3d) - s(i,j-1,k3d)
             dabarp    = s(i,j+1,k3d) - s(i,j,k3d)
             dafacemin = min(abs(dafacem),abs(dafacep))
             dabarmin  = min(abs(dabarm),abs(dabarp))
             if (dafacemin.ge.dabarmin) then
                dachkm = dafacem
                dachkp = dafacep
             else
                dachkm = dabarm
                dachkp = dabarp
             endif
             extremum = (dachkm*dachkp .le. 0.d0)
          end if

          if (extremum) then
             D2     = 6.d0*(alpham + alphap)
             D2L    = s(i,j-2,k3d)-2.d0*s(i,j-1,k3d)+s(i,j,k3d)
             D2R    = s(i,j,k3d)-2.d0*s(i,j+1,k3d)+s(i,j+2,k3d)
             D2C    = s(i,j-1,k3d)-2.d0*s(i,j,k3d)+s(i,j+1,k3d)
             sgn    = sign(1.d0,D2)
             D2LIM  = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),0.d0)
             alpham = alpham*D2LIM/max(abs(D2),1.d-10)
             alphap = alphap*D2LIM/max(abs(D2),1.d-10)
          else
             if (bigp) then
                sgn   = sign(1.d0,alpham)
                amax  = -alphap**2 / (4*(alpham + alphap))
                delam = s(i,j-1,k3d) - s(i,j,k3d)
                if (sgn*amax .ge. sgn*delam) then
                   if (sgn*(delam - alpham).ge.1.d-10) then
                      alphap = (-2.d0*delam - 2.d0*sgn*sqrt(delam**2 - delam*alpham))
                   else 
                      alphap = -2.d0*alpham
                   endif
                endif
             end if
             if (bigm) then
                sgn   = sign(1.d0,alphap)
                amax  = -alpham**2 / (4*(alpham + alphap))
                delap = s(i,j+1,k3d) - s(i,j,k3d)
                if (sgn*amax .ge. sgn*delap) then
                   if (sgn*(delap - alphap).ge.1.d-10) then
                      alpham = (-2.d0*delap - 2.d0*sgn*sqrt(delap**2 - delap*alphap))
                   else
                      alpham = -2.d0*alphap
                   endif
                endif
             end if
          end if

          sm(i,j) = s(i,j,k3d) + alpham
          sp(i,j) = s(i,j,k3d) + alphap
          !
          ! Compute y-component of Ip and Im.
          !
          s6    = 6.0d0*s(i,j,k3d) - 3.0d0*(sm(i,j)+sp(i,j))
          sigma = abs(u(i,j,k3d,2)-cspd(i,j,k3d))*dt/dy

          Ip(i,j,kc,2,1) = sp(i,j) - &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)-(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
          Im(i,j,kc,2,1) = sm(i,j) + &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)+(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
          sigma = abs(u(i,j,k3d,2))*dt/dy
          Ip(i,j,kc,2,2) = sp(i,j) - &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)-(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
          Im(i,j,kc,2,2) = sm(i,j) + &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)+(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
          sigma = abs(u(i,j,k3d,2)+cspd(i,j,k3d))*dt/dy
          Ip(i,j,kc,2,3) = sp(i,j) - &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)-(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
          Im(i,j,kc,2,3) = sm(i,j) + &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)+(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
       end do
    end do
    !$OMP END PARALLEL DO

    deallocate(dsvl,sedge)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! z-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! cell-centered indexing
    allocate( dsvl(ilo1-1:ihi1+1,ilo2-1:ihi2+1))
    allocate(dsvlm(ilo1-1:ihi1+1,ilo2-1:ihi2+1))
    allocate(dsvlp(ilo1-1:ihi1+1,ilo2-1:ihi2+1))

    allocate(sedgez(ilo1-1:ihi1+1,ilo2-2:ihi2+3,k3d-1:k3d+2))

    ! compute s at z-edges

    ! interpolate s to z-edges
    !$OMP PARALLEL DO PRIVATE(i,j,k,D2,D2L,D2R,sgn,D2LIM)
    do k=k3d-1,k3d+2
       do j=ilo2-1,ihi2+1
          do i=ilo1-1,ihi1+1
             sedgez(i,j,k) = (7.d0/12.d0)*(s(i,j,k-1)+s(i,j,k)) &
                  - (1.d0/12.d0)*(s(i,j,k-2)+s(i,j,k+1))
             !
             ! limit sedgez
             !
             if ((sedgez(i,j,k)-s(i,j,k-1))*(s(i,j,k)-sedgez(i,j,k)) .lt. 0.d0) then
                D2  = 3.d0*(s(i,j,k-1)-2.d0*sedgez(i,j,k)+s(i,j,k))
                D2L = s(i,j,k-2)-2.d0*s(i,j,k-1)+s(i,j,k)
                D2R = s(i,j,k-1)-2.d0*s(i,j,k)+s(i,j,k+1)
                sgn = sign(1.d0,D2)
                D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),0.d0)
                sedgez(i,j,k) = 0.5d0*(s(i,j,k-1)+s(i,j,k)) - (1.d0/6.d0)*D2LIM
             end if
          end do
       end do
    end do
    !$OMP END PARALLEL DO
    !
    ! Use Colella 2008 limiters.
    !
    ! This is a new version of the algorithm to eliminate sensitivity to roundoff.
    !
    k = k3d
    !$OMP PARALLEL DO PRIVATE(i,j,alphap,alpham,bigp,bigm,extremum,dafacem,dafacep,dabarm,dabarp,dafacemin) &
    !$OMP PRIVATE(dabarmin,dachkm,dachkp,D2,D2L,D2R,D2C,sgn,D2LIM,amax,delam,delap,s6,sigma)
    do j=ilo2-1,ihi2+1
       do i=ilo1-1,ihi1+1

          alphap   = sedgez(i,j,k+1)-s(i,j,k)
          alpham   = sedgez(i,j,k  )-s(i,j,k)
          bigp     = abs(alphap).gt.2.d0*abs(alpham)
          bigm     = abs(alpham).gt.2.d0*abs(alphap)
          extremum = .false.

          if (alpham*alphap .ge. 0.d0) then
             extremum = .true.
          else if (bigp .or. bigm) then
             !
             ! Possible extremum. We look at cell centered values and face
             ! centered values for a change in sign in the differences adjacent to
             ! the cell. We use the pair of differences whose minimum magnitude is the
             ! largest, and thus least susceptible to sensitivity to roundoff.
             !
             dafacem   = sedgez(i,j,k) - sedgez(i,j,k-1)
             dafacep   = sedgez(i,j,k+2) - sedgez(i,j,k+1)
             dabarm    = s(i,j,k) - s(i,j,k-1)
             dabarp    = s(i,j,k+1) - s(i,j,k)
             dafacemin = min(abs(dafacem),abs(dafacep))
             dabarmin  = min(abs(dabarm),abs(dabarp))
             if (dafacemin.ge.dabarmin) then
                dachkm = dafacem
                dachkp = dafacep
             else
                dachkm = dabarm
                dachkp = dabarp
             endif
             extremum = (dachkm*dachkp .le. 0.d0)
          end if

          if (extremum) then
             D2     = 6.d0*(alpham + alphap)
             D2L    = s(i,j,k-2)-2.d0*s(i,j,k-1)+s(i,j,k)
             D2R    = s(i,j,k)-2.d0*s(i,j,k+1)+s(i,j,k+2)
             D2C    = s(i,j,k-1)-2.d0*s(i,j,k)+s(i,j,k+1)
             sgn    = sign(1.d0,D2)
             D2LIM  = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),0.d0)
             alpham = alpham*D2LIM/max(abs(D2),1.d-10)
             alphap = alphap*D2LIM/max(abs(D2),1.d-10)
          else
             if (bigp) then
                sgn   = sign(1.d0,alpham)
                amax  = -alphap**2 / (4*(alpham + alphap))
                delam = s(i,j,k-1) - s(i,j,k)
                if (sgn*amax .ge. sgn*delam) then
                   if (sgn*(delam - alpham).ge.1.d-10) then
                      alphap = (-2.d0*delam - 2.d0*sgn*sqrt(delam**2 - delam*alpham))
                   else 
                      alphap = -2.d0*alpham
                   endif
                endif
             end if
             if (bigm) then
                sgn   = sign(1.d0,alphap)
                amax  = -alpham**2 / (4*(alpham + alphap))
                delap = s(i,j,k+1) - s(i,j,k)
                if (sgn*amax .ge. sgn*delap) then
                   if (sgn*(delap - alphap).ge.1.d-10) then
                      alpham = (-2.d0*delap - 2.d0*sgn*sqrt(delap**2 - delap*alphap))
                   else
                      alpham = -2.d0*alphap
                   endif
                endif
             end if
          end if

          sm(i,j) = s(i,j,k) + alpham
          sp(i,j) = s(i,j,k) + alphap
          !
          ! Compute z-component of Ip and Im.
          !
          s6    = 6.0d0*s(i,j,k3d) - 3.0d0*(sm(i,j)+sp(i,j))
          sigma = abs(u(i,j,k3d,3)-cspd(i,j,k3d))*dt/dz

          Ip(i,j,kc,3,1) = sp(i,j) - &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)-(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
          Im(i,j,kc,3,1) = sm(i,j) + &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)+(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
          sigma = abs(u(i,j,k3d,3))*dt/dz
          Ip(i,j,kc,3,2) = sp(i,j) - &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)-(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
          Im(i,j,kc,3,2) = sm(i,j) + &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)+(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
          sigma = abs(u(i,j,k3d,3)+cspd(i,j,k3d))*dt/dz
          Ip(i,j,kc,3,3) = sp(i,j) - &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)-(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
          Im(i,j,kc,3,3) = sm(i,j) + &
               (sigma/2.0d0)*(sp(i,j)-sm(i,j)+(1.0d0-(2.0d0/3.0d0)*sigma)*s6)

       end do
    end do
    !$OMP END PARALLEL DO

    deallocate(dsvl,dsvlm,dsvlp,sp,sm,sedgez)

  end subroutine ppm_type2

  ! ::: 
  ! ::: ------------------------------------------------------------------
  ! ::: 

  subroutine tracexy_ppm(q,c,flatn,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                         Ip,Im, &
                         qxm,qxp,qym,qyp,qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
                         ilo1,ilo2,ihi1,ihi2,dx,dy,dt,kc,k3d)

    use network, only : nspec
    use meth_params_module, only : iorder, QVAR, QRHO, QU, QV, QW, &
         QREINT, QPRES, QFA, QFS, nadv, small_dens, &
         ppm_type

    implicit none

    integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
    integer qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3
    integer ilo1,ilo2,ihi1,ihi2
    integer kc,k3d

    double precision     q(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
    double precision     c(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
    double precision flatn(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)

    double precision   Ip(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3,QVAR)
    double precision   Im(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3,QVAR)

    double precision qxm(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
    double precision qxp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
    double precision qym(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
    double precision qyp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
    double precision dx, dy, dt

    ! Local variables
    integer i, j
    integer n, iadv
    integer ns, ispec, iaux

    double precision dtdx, dtdy
    double precision cc, csq, rho, u, v, w, p, rhoe
    double precision drho, du, dv, dw, dp, drhoe
    double precision dup, dvp, dwp, dpp
    double precision dum, dvm, dwm, dpm

    double precision enth, alpham, alphap, alpha0r, alpha0e
    double precision alpha0u, alpha0v, alpha0w
    double precision apright, amright, azrright, azeright
    double precision azu1rght, azv1rght, azw1rght
    double precision apleft, amleft, azrleft, azeleft
    double precision azu1left, azv1left, azw1left

    if (ppm_type .eq. 0) then
       print *,'Oops -- shouldnt be in tracexy_ppm with ppm_type = 0'
       call bl_error("Error:: ppm_3d.f90 :: tracexy_ppm")
    end if

    dtdx = dt/dx
    dtdy = dt/dy

    !!!!!!!!!!!!!!!
    ! PPM CODE
    !!!!!!!!!!!!!!!

    ! Trace to left and right edges using upwind PPM
    !$OMP PARALLEL DO PRIVATE(i,j,cc,csq,rho,u,v,w,p,rhoe,enth,dum,dvm,dwm,dpm) &
    !$OMP PRIVATE(drho,du,dv,dw,dp,drhoe,dup,dvp,dwp,dpp,alpham,alphap,alpha0r) &
    !$OMP PRIVATE(alpha0e,alpha0v,alpha0w,amright,apright,azrright,azeright,azv1rght,azw1rght) &
    !$OMP PRIVATE(amleft,apleft,azrleft,azeleft,azv1left,azw1left)
    do j = ilo2-1, ihi2+1
       do i = ilo1-1, ihi1+1

          cc = c(i,j,k3d)
          csq = cc**2
          rho = q(i,j,k3d,QRHO)
          u = q(i,j,k3d,QU)
          v = q(i,j,k3d,QV)
          w = q(i,j,k3d,QW)
          p = q(i,j,k3d,QPRES)
          rhoe = q(i,j,k3d,QREINT)
          enth = ( (rhoe+p)/rho )/csq

          ! plus state on face i
          dum    = flatn(i,j,k3d)*(u    - Im(i,j,kc,1,1,QU))
          dvm    = flatn(i,j,k3d)*(v    - Im(i,j,kc,1,1,QV))
          dwm    = flatn(i,j,k3d)*(w    - Im(i,j,kc,1,1,QW))
          dpm    = flatn(i,j,k3d)*(p    - Im(i,j,kc,1,1,QPRES))

          drho  = flatn(i,j,k3d)*(rho  - Im(i,j,kc,1,2,QRHO))
          du    = flatn(i,j,k3d)*(u    - Im(i,j,kc,1,2,QU))
          dv    = flatn(i,j,k3d)*(v    - Im(i,j,kc,1,2,QV))
          dw    = flatn(i,j,k3d)*(w    - Im(i,j,kc,1,2,QW))
          dp    = flatn(i,j,k3d)*(p    - Im(i,j,kc,1,2,QPRES))
          drhoe = flatn(i,j,k3d)*(rhoe - Im(i,j,kc,1,2,QREINT))

          dup    = flatn(i,j,k3d)*(u    - Im(i,j,kc,1,3,QU))
          dvp    = flatn(i,j,k3d)*(v    - Im(i,j,kc,1,3,QV))
          dwp    = flatn(i,j,k3d)*(w    - Im(i,j,kc,1,3,QW))
          dpp    = flatn(i,j,k3d)*(p    - Im(i,j,kc,1,3,QPRES))

          alpham = 0.5d0*(dpm/(rho*cc) - dum)*rho/cc
          alphap = 0.5d0*(dpp/(rho*cc) + dup)*rho/cc
          alpha0r = drho - dp/csq
          alpha0e = drhoe - dp*enth
          alpha0v = dv
          alpha0w = dw

          if (u-cc .gt. 0.d0) then
             amright = 0.d0
          else if (u-cc .lt. 0.d0) then
             amright = -alpham
          else
             amright = -0.5d0*alpham
          endif
          if (u+cc .gt. 0.d0) then
             apright = 0.d0
          else if (u+cc .lt. 0.d0) then
             apright = -alphap
          else
             apright = -0.5d0*alphap
          endif
          if (u .gt. 0.d0) then
             azrright = 0.d0
             azeright = 0.d0
             azv1rght = 0.d0
             azw1rght = 0.d0
          else if (u .lt. 0.d0) then
             azrright = -alpha0r
             azeright = -alpha0e
             azv1rght = -alpha0v
             azw1rght = -alpha0w
          else
             azrright = -0.5d0*alpha0r
             azeright = -0.5d0*alpha0e
             azv1rght = -0.5d0*alpha0v
             azw1rght = -0.5d0*alpha0w
          endif

          if (i .ge. ilo1) then
             qxp(i,j,kc,QRHO) = rho + apright + amright + azrright
             qxp(i,j,kc,QRHO) = max(small_dens,qxp(i,j,kc,QRHO))
             qxp(i,j,kc,QU) = u + (apright - amright)*cc/rho
             qxp(i,j,kc,QV) = v + azv1rght
             qxp(i,j,kc,QW) = w + azw1rght
             qxp(i,j,kc,QPRES) = p + (apright + amright)*csq
             qxp(i,j,kc,QREINT) = rhoe + (apright + amright)*enth*csq + azeright
          end if

          ! minus state on face i+1
          dum    = flatn(i,j,k3d)*(u    - Ip(i,j,kc,1,1,QU))
          dvm    = flatn(i,j,k3d)*(v    - Ip(i,j,kc,1,1,QV))
          dwm    = flatn(i,j,k3d)*(w    - Ip(i,j,kc,1,1,QW))
          dpm    = flatn(i,j,k3d)*(p    - Ip(i,j,kc,1,1,QPRES))

          drho  = flatn(i,j,k3d)*(rho  - Ip(i,j,kc,1,2,QRHO))
          du    = flatn(i,j,k3d)*(u    - Ip(i,j,kc,1,2,QU))
          dv    = flatn(i,j,k3d)*(v    - Ip(i,j,kc,1,2,QV))
          dw    = flatn(i,j,k3d)*(w    - Ip(i,j,kc,1,2,QW))
          dp    = flatn(i,j,k3d)*(p    - Ip(i,j,kc,1,2,QPRES))
          drhoe = flatn(i,j,k3d)*(rhoe - Ip(i,j,kc,1,2,QREINT))

          dup    = flatn(i,j,k3d)*(u    - Ip(i,j,kc,1,3,QU))
          dvp    = flatn(i,j,k3d)*(v    - Ip(i,j,kc,1,3,QV))
          dwp    = flatn(i,j,k3d)*(w    - Ip(i,j,kc,1,3,QW))
          dpp    = flatn(i,j,k3d)*(p    - Ip(i,j,kc,1,3,QPRES))

          alpham = 0.5d0*(dpm/(rho*cc) - dum)*rho/cc
          alphap = 0.5d0*(dpp/(rho*cc) + dup)*rho/cc
          alpha0r = drho - dp/csq
          alpha0e = drhoe - dp*enth
          alpha0v = dv
          alpha0w = dw

          if (u-cc .gt. 0.d0) then
             amleft = -alpham
          else if (u-cc .lt. 0.d0) then
             amleft = 0.d0
          else
             amleft = -0.5d0*alpham
          endif
          if (u+cc .gt. 0.d0) then
             apleft = -alphap
          else if (u+cc .lt. 0.d0) then
             apleft = 0.d0
          else
             apleft = -0.5d0*alphap
          endif
          if (u .gt. 0.d0) then
             azrleft = -alpha0r
             azeleft = -alpha0e
             azv1left = -alpha0v
             azw1left = -alpha0w
          else if (u .lt. 0.d0) then
             azrleft = 0.d0
             azeleft = 0.d0
             azv1left = 0.d0
             azw1left = 0.d0
          else
             azrleft = -0.5d0*alpha0r
             azeleft = -0.5d0*alpha0e
             azv1left = -0.5d0*alpha0v
             azw1left = -0.5d0*alpha0w
          endif

          if (i .le. ihi1) then
             qxm(i+1,j,kc,QRHO) = rho + apleft + amleft + azrleft
             qxm(i+1,j,kc,QRHO) = max(qxm(i+1,j,kc,QRHO),small_dens)
             qxm(i+1,j,kc,QU) = u + (apleft - amleft)*cc/rho
             qxm(i+1,j,kc,QV) = v + azv1left
             qxm(i+1,j,kc,QW) = w + azw1left
             qxm(i+1,j,kc,QPRES) = p + (apleft + amleft)*csq
             qxm(i+1,j,kc,QREINT) = rhoe + (apleft + amleft)*enth*csq + azeleft
          end if

       end do
    end do
    !$OMP END PARALLEL DO

    ! Now do the passively advected quantities
    !$OMP PARALLEL DO PRIVATE(iadv,n,i,j,u) IF(nadv.gt.1)
    do iadv = 1, nadv
       n = QFA + iadv - 1
       do j = ilo2-1, ihi2+1

          ! plus state on face i
          do i = ilo1, ihi1+1
             u = q(i,j,k3d,QU)
             if (u .gt. 0.d0) then
                qxp(i,j,kc,n) = q(i,j,k3d,n)
             else if (u .lt. 0.d0) then
                qxp(i,j,kc,n) = q(i,j,k3d,n) &
                     + flatn(i,j,k3d)*(Im(i,j,kc,1,2,n) - q(i,j,k3d,n))
             else
                qxp(i,j,kc,n) = q(i,j,k3d,n) &
                     + 0.5d0*flatn(i,j,k3d)*(Im(i,j,kc,1,2,n) - q(i,j,k3d,n))
             endif
          enddo

          ! minus state on face i+1
          do i = ilo1-1, ihi1
             u = q(i,j,k3d,QU)
             if (u .gt. 0.d0) then
                qxm(i+1,j,kc,n) = q(i,j,k3d,n) &
                     + flatn(i,j,k3d)*(Ip(i,j,kc,1,2,n) - q(i,j,k3d,n))
             else if (u .lt. 0.d0) then
                qxm(i+1,j,kc,n) = q(i,j,k3d,n)
             else
                qxm(i+1,j,kc,n) = q(i,j,k3d,n) &
                     + 0.5d0*flatn(i,j,k3d)*(Ip(i,j,kc,1,2,n) - q(i,j,k3d,n))
             endif
          enddo

       enddo
    enddo
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(ispec,ns,i,j,u) IF(nspec.gt.1)
    do ispec = 1, nspec
       ns = QFS + ispec - 1

       do j = ilo2-1, ihi2+1

          ! plus state on face i
          do i = ilo1, ihi1+1
             u = q(i,j,k3d,QU)
             if (u .gt. 0.d0) then
                qxp(i,j,kc,ns) = q(i,j,k3d,ns)
             else if (u .lt. 0.d0) then
                qxp(i,j,kc,ns) = q(i,j,k3d,ns) &
                     + flatn(i,j,k3d)*(Im(i,j,kc,1,2,ns) - q(i,j,k3d,ns))
             else
                qxp(i,j,kc,ns) = q(i,j,k3d,ns) &
                     + 0.5d0*flatn(i,j,k3d)*(Im(i,j,kc,1,2,ns) - q(i,j,k3d,ns))
             endif
          enddo

          ! minus state on face i+1
          do i = ilo1-1, ihi1
             u = q(i,j,k3d,QU)
             if (u .gt. 0.d0) then
                qxm(i+1,j,kc,ns) = q(i,j,k3d,ns) &
                     + flatn(i,j,k3d)*(Ip(i,j,kc,1,2,ns) - q(i,j,k3d,ns))
             else if (u .lt. 0.d0) then
                qxm(i+1,j,kc,ns) = q(i,j,k3d,ns)
             else
                qxm(i+1,j,kc,ns) = q(i,j,k3d,ns) &
                     + 0.5d0*flatn(i,j,k3d)*(Ip(i,j,kc,1,2,ns) - q(i,j,k3d,ns))
             endif
          enddo

       enddo
    enddo
    !$OMP END PARALLEL DO

    ! Trace to bottom and top edges using upwind PPM
    !$OMP PARALLEL DO PRIVATE(i,j,cc,csq,rho,u,v,w,p,rhoe,enth,dum,dvm,dwm,dpm) &
    !$OMP PRIVATE(drho,du,dv,dw,dp,drhoe,dup,dvp,dwp,dpp,alpham,alphap,alpha0r) &
    !$OMP PRIVATE(alpha0e,alpha0u,alpha0w,amright,apright,azrright,azeright,azu1rght,azw1rght,amleft) &
    !$OMP PRIVATE(apleft,azrleft,azeleft,azu1left,azw1left)
    do j = ilo2-1, ihi2+1
       do i = ilo1-1, ihi1+1

          cc = c(i,j,k3d)
          csq = cc**2
          rho = q(i,j,k3d,QRHO)
          u = q(i,j,k3d,QU)
          v = q(i,j,k3d,QV)
          w = q(i,j,k3d,QW)
          p = q(i,j,k3d,QPRES)
          rhoe = q(i,j,k3d,QREINT)
          enth = ( (rhoe+p)/rho )/csq

          ! plus state on face j
          dum    = flatn(i,j,k3d)*(u    - Im(i,j,kc,2,1,QU))
          dvm    = flatn(i,j,k3d)*(v    - Im(i,j,kc,2,1,QV))
          dwm    = flatn(i,j,k3d)*(w    - Im(i,j,kc,2,1,QW))
          dpm    = flatn(i,j,k3d)*(p    - Im(i,j,kc,2,1,QPRES))

          drho  = flatn(i,j,k3d)*(rho  - Im(i,j,kc,2,2,QRHO))
          du    = flatn(i,j,k3d)*(u    - Im(i,j,kc,2,2,QU))
          dv    = flatn(i,j,k3d)*(v    - Im(i,j,kc,2,2,QV))
          dw    = flatn(i,j,k3d)*(w    - Im(i,j,kc,2,2,QW))
          dp    = flatn(i,j,k3d)*(p    - Im(i,j,kc,2,2,QPRES))
          drhoe = flatn(i,j,k3d)*(rhoe - Im(i,j,kc,2,2,QREINT))

          dup    = flatn(i,j,k3d)*(u    - Im(i,j,kc,2,3,QU))
          dvp    = flatn(i,j,k3d)*(v    - Im(i,j,kc,2,3,QV))
          dwp    = flatn(i,j,k3d)*(w    - Im(i,j,kc,2,3,QW))
          dpp    = flatn(i,j,k3d)*(p    - Im(i,j,kc,2,3,QPRES))

          alpham = 0.5d0*(dpm/(rho*cc) - dvm)*rho/cc
          alphap = 0.5d0*(dpp/(rho*cc) + dvp)*rho/cc
          alpha0r = drho - dp/csq
          alpha0e = drhoe - dp*enth
          alpha0u = du
          alpha0w = dw

          if (v-cc .gt. 0.d0) then
             amright = 0.d0
          else if (v-cc .lt. 0.d0) then
             amright = -alpham
          else
             amright = -0.5d0*alpham
          endif
          if (v+cc .gt. 0.d0) then
             apright = 0.d0
          else if (v+cc .lt. 0.d0) then
             apright = -alphap
          else
             apright = -0.5d0*alphap
          endif
          if (v .gt. 0.d0) then
             azrright = 0.d0
             azeright = 0.d0
             azu1rght = 0.d0
             azw1rght = 0.d0
          else if (v .lt. 0.d0) then
             azrright = -alpha0r
             azeright = -alpha0e
             azu1rght = -alpha0u
             azw1rght = -alpha0w
          else
             azrright = -0.5d0*alpha0r
             azeright = -0.5d0*alpha0e
             azu1rght = -0.5d0*alpha0u
             azw1rght = -0.5d0*alpha0w
          endif

          if (j .ge. ilo2) then
             qyp(i,j,kc,QRHO) = rho + apright + amright + azrright
             qyp(i,j,kc,QRHO) = max(small_dens, qyp(i,j,kc,QRHO))
             qyp(i,j,kc,QV) = v + (apright - amright)*cc/rho
             qyp(i,j,kc,QU) = u + azu1rght
             qyp(i,j,kc,QW) = w + azw1rght
             qyp(i,j,kc,QPRES) = p + (apright + amright)*csq
             qyp(i,j,kc,QREINT) = rhoe + (apright + amright)*enth*csq + azeright
          end if

          ! minus state on face j+1
          dum    = flatn(i,j,k3d)*(u    - Ip(i,j,kc,2,1,QU))
          dvm    = flatn(i,j,k3d)*(v    - Ip(i,j,kc,2,1,QV))
          dwm    = flatn(i,j,k3d)*(w    - Ip(i,j,kc,2,1,QW))
          dpm    = flatn(i,j,k3d)*(p    - Ip(i,j,kc,2,1,QPRES))

          drho  = flatn(i,j,k3d)*(rho  - Ip(i,j,kc,2,2,QRHO))
          du    = flatn(i,j,k3d)*(u    - Ip(i,j,kc,2,2,QU))
          dv    = flatn(i,j,k3d)*(v    - Ip(i,j,kc,2,2,QV))
          dw    = flatn(i,j,k3d)*(w    - Ip(i,j,kc,2,2,QW))
          dp    = flatn(i,j,k3d)*(p    - Ip(i,j,kc,2,2,QPRES))
          drhoe = flatn(i,j,k3d)*(rhoe - Ip(i,j,kc,2,2,QREINT))

          dup    = flatn(i,j,k3d)*(u    - Ip(i,j,kc,2,3,QU))
          dvp    = flatn(i,j,k3d)*(v    - Ip(i,j,kc,2,3,QV))
          dwp    = flatn(i,j,k3d)*(w    - Ip(i,j,kc,2,3,QW))
          dpp    = flatn(i,j,k3d)*(p    - Ip(i,j,kc,2,3,QPRES))

          alpham = 0.5d0*(dpm/(rho*cc) - dvm)*rho/cc
          alphap = 0.5d0*(dpp/(rho*cc) + dvp)*rho/cc
          alpha0r = drho - dp/csq
          alpha0e = drhoe - dp*enth
          alpha0u = du
          alpha0w = dw

          if (v-cc .gt. 0.d0) then
             amleft = -alpham
          else if (v-cc .lt. 0.d0) then
             amleft = 0.d0
          else
             amleft = -0.5d0*alpham
          endif
          if (v+cc .gt. 0.d0) then
             apleft = -alphap
          else if (v+cc .lt. 0.d0) then
             apleft = 0.d0
          else
             apleft = -0.5d0*alphap
          endif
          if (v .gt. 0.d0) then
             azrleft = -alpha0r
             azeleft = -alpha0e
             azu1left = -alpha0u
             azw1left = -alpha0w
          else if (v .lt. 0.d0) then
             azrleft = 0.d0
             azeleft = 0.d0
             azu1left = 0.d0
             azw1left = 0.d0
          else
             azrleft = -0.5d0*alpha0r
             azeleft = -0.5d0*alpha0e
             azu1left = -0.5d0*alpha0u
             azw1left = -0.5d0*alpha0w
          endif

          if (j .le. ihi2) then
             qym(i,j+1,kc,QRHO) = rho + apleft + amleft + azrleft
             qym(i,j+1,kc,QRHO) = max(small_dens, qym(i,j+1,kc,QRHO))
             qym(i,j+1,kc,QV) = v + (apleft - amleft)*cc/rho
             qym(i,j+1,kc,QU) = u + azu1left
             qym(i,j+1,kc,QW) = w + azw1left
             qym(i,j+1,kc,QPRES) = p + (apleft + amleft)*csq
             qym(i,j+1,kc,QREINT) = rhoe + (apleft + amleft)*enth*csq + azeleft
          end if

       end do
    end do
    !$OMP END PARALLEL DO

    ! Now do the passively advected quantities
    !$OMP PARALLEL DO PRIVATE(iadv,n,i,j,v) IF(nadv.gt.1)
    do iadv = 1, nadv
       n = QFA + iadv - 1
       do i = ilo1-1, ihi1+1

          ! plus state on face j
          do j = ilo2, ihi2+1
             v = q(i,j,k3d,QV)
             if (v .gt. 0.d0) then
                qyp(i,j,kc,n) = q(i,j,k3d,n)
             else if (v .lt. 0.d0) then
                qyp(i,j,kc,n) = q(i,j,k3d,n) &
                     + flatn(i,j,k3d)*(Im(i,j,kc,2,2,n) - q(i,j,k3d,n))
             else
                qyp(i,j,kc,n) = q(i,j,k3d,n) &
                     + 0.5d0*flatn(i,j,k3d)*(Im(i,j,kc,2,2,n) - q(i,j,k3d,n))
             endif
          enddo

          ! minus state on face j+1
          do j = ilo2-1, ihi2
             v = q(i,j,k3d,QV)
             if (v .gt. 0.d0) then
                qym(i,j+1,kc,n) = q(i,j,k3d,n) &
                     + flatn(i,j,k3d)*(Ip(i,j,kc,2,2,n) - q(i,j,k3d,n))
             else if (v .lt. 0.d0) then
                qym(i,j+1,kc,n) = q(i,j,k3d,n)
             else
                qym(i,j+1,kc,n) = q(i,j,k3d,n) &
                     + 0.5d0*flatn(i,j,k3d)*(Ip(i,j,kc,2,2,n) - q(i,j,k3d,n))
             endif
          enddo

       enddo
    enddo
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(ispec,ns,i,j,v) IF(nspec.gt.1)
    do ispec = 1, nspec
       ns = QFS + ispec - 1
       do i = ilo1-1, ihi1+1

          ! plus state on face j
          do j = ilo2, ihi2+1
             v = q(i,j,k3d,QV)
             if (v .gt. 0.d0) then
                qyp(i,j,kc,ns) = q(i,j,k3d,ns)
             else if (v .lt. 0.d0) then
                qyp(i,j,kc,ns) = q(i,j,k3d,ns) &
                     + flatn(i,j,k3d)*(Im(i,j,kc,2,2,ns) - q(i,j,k3d,ns))
             else
                qyp(i,j,kc,ns) = q(i,j,k3d,ns) &
                     + 0.5d0*flatn(i,j,k3d)*(Im(i,j,kc,2,2,ns) - q(i,j,k3d,ns))
             endif
          enddo

          ! minus state on face j+1
          do j = ilo2-1, ihi2
             v = q(i,j,k3d,QV)
             if (v .gt. 0.d0) then
                qym(i,j+1,kc,ns) = q(i,j,k3d,ns) &
                     + flatn(i,j,k3d)*(Ip(i,j,kc,2,2,ns) - q(i,j,k3d,ns))
             else if (v .lt. 0.d0) then
                qym(i,j+1,kc,ns) = q(i,j,k3d,ns)
             else
                qym(i,j+1,kc,ns) = q(i,j,k3d,ns) &
                     + 0.5d0*flatn(i,j,k3d)*(Ip(i,j,kc,2,2,ns) - q(i,j,k3d,ns))
             endif
          enddo

       enddo
    enddo
    !$OMP END PARALLEL DO

  end subroutine tracexy_ppm

  ! ::: 
  ! ::: ------------------------------------------------------------------
  ! ::: 

  subroutine tracez_ppm(q,c,flatn,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                        Ip,Im, &
                        qzm,qzp,qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
                        ilo1,ilo2,ihi1,ihi2,dz,dt,km,kc,k3d)

    use network, only : nspec
    use meth_params_module, only : iorder, QVAR, QRHO, QU, QV, QW, &
         QREINT, QPRES, QFA, QFS, nadv, small_dens, &
         ppm_type

    implicit none

    integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
    integer qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3
    integer ilo1,ilo2,ihi1,ihi2
    integer km,kc,k3d

    double precision     q(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,QVAR)
    double precision     c(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
    double precision flatn(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)

    double precision   Ip(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3,QVAR)
    double precision   Im(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3,QVAR)
    double precision qzm(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
    double precision qzp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
    double precision dz, dt

    !     Local variables
    integer i, j
    integer n, iadv
    integer ns, ispec, iaux

    double precision dtdz
    double precision cc, csq, rho, u, v, w, p, rhoe
    double precision dup, dvp, dwp, dpp
    double precision dum, dvm, dwm, dpm

    double precision drho, du, dv, dw, dp, drhoe
    double precision enth, alpham, alphap, alpha0r, alpha0e
    double precision alpha0u, alpha0v
    double precision apright, amright, azrright, azeright
    double precision azu1rght, azv1rght
    double precision apleft, amleft, azrleft, azeleft
    double precision azu1left, azv1left

    if (ppm_type .eq. 0) then
       print *,'Oops -- shouldnt be in tracez_ppm with ppm_type = 0'
       call bl_error("Error:: ppm_3d.f90 :: tracez_ppm")
    end if

    dtdz = dt/dz

    !!!!!!!!!!!!!!!
    ! PPM CODE
    !!!!!!!!!!!!!!!

    ! Trace to left and right edges using upwind PPM
    !$OMP PARALLEL DO PRIVATE(i,j,cc,csq,rho,u,v,w,p,rhoe,enth,dum,dvm,dwm,dpm) &
    !$OMP PRIVATE(drho,du,dv,dw,dp,drhoe,dup,dvp,dwp,dpp,alpham,alphap,alpha0r,alpha0e) &
    !$OMP PRIVATE(alpha0u,alpha0v,amright,apright,azrright,azeright,azu1rght,azv1rght,amleft,apleft)&
    !$OMP PRIVATE(azrleft,azeleft,azu1left,azv1left)
    do j = ilo2-1, ihi2+1
       do i = ilo1-1, ihi1+1

          cc = c(i,j,k3d)
          csq = cc**2
          rho = q(i,j,k3d,QRHO)
          u = q(i,j,k3d,QU)
          v = q(i,j,k3d,QV)
          w = q(i,j,k3d,QW)
          p = q(i,j,k3d,QPRES)
          rhoe = q(i,j,k3d,QREINT)
          enth = ( (rhoe+p)/rho )/csq

          ! plus state on face kc
          dum    = flatn(i,j,k3d)*(u    - Im(i,j,kc,3,1,QU))
          dvm    = flatn(i,j,k3d)*(v    - Im(i,j,kc,3,1,QV))
          dwm    = flatn(i,j,k3d)*(w    - Im(i,j,kc,3,1,QW))
          dpm    = flatn(i,j,k3d)*(p    - Im(i,j,kc,3,1,QPRES))

          drho  = flatn(i,j,k3d)*(rho  - Im(i,j,kc,3,2,QRHO))
          du    = flatn(i,j,k3d)*(u    - Im(i,j,kc,3,2,QU))
          dv    = flatn(i,j,k3d)*(v    - Im(i,j,kc,3,2,QV))
          dw    = flatn(i,j,k3d)*(w    - Im(i,j,kc,3,2,QW))
          dp    = flatn(i,j,k3d)*(p    - Im(i,j,kc,3,2,QPRES))
          drhoe = flatn(i,j,k3d)*(rhoe - Im(i,j,kc,3,2,QREINT))

          dup    = flatn(i,j,k3d)*(u    - Im(i,j,kc,3,3,QU))
          dvp    = flatn(i,j,k3d)*(v    - Im(i,j,kc,3,3,QV))
          dwp    = flatn(i,j,k3d)*(w    - Im(i,j,kc,3,3,QW))
          dpp    = flatn(i,j,k3d)*(p    - Im(i,j,kc,3,3,QPRES))

          alpham = 0.5d0*(dpm/(rho*cc) - dwm)*rho/cc
          alphap = 0.5d0*(dpp/(rho*cc) + dwp)*rho/cc
          alpha0r = drho - dp/csq
          alpha0e = drhoe - dp*enth
          alpha0u = du
          alpha0v = dv

          if (w-cc .gt. 0.d0) then
             amright = 0.d0
          else if (w-cc .lt. 0.d0) then
             amright = -alpham
          else
             amright = -0.5d0*alpham
          endif
          if (w+cc .gt. 0.d0) then
             apright = 0.d0
          else if (w+cc .lt. 0.d0) then
             apright = -alphap
          else
             apright = -0.5d0*alphap
          endif
          if (w .gt. 0.d0) then
             azrright = 0.d0
             azeright = 0.d0
             azu1rght = 0.d0
             azv1rght = 0.d0
          else if (w .lt. 0.d0) then
             azrright = -alpha0r
             azeright = -alpha0e
             azu1rght = -alpha0u
             azv1rght = -alpha0v
          else
             azrright = -0.5d0*alpha0r
             azeright = -0.5d0*alpha0e
             azu1rght = -0.5d0*alpha0u
             azv1rght = -0.5d0*alpha0v
          endif

          qzp(i,j,kc,QRHO) = rho + apright + amright + azrright
          qzp(i,j,kc,QRHO) = max(small_dens, qzp(i,j,kc,QRHO))
          qzp(i,j,kc,QW) = w + (apright - amright)*cc/rho
          qzp(i,j,kc,QU) = u + azu1rght
          qzp(i,j,kc,QV) = v + azv1rght
          qzp(i,j,kc,QPRES) = p + (apright + amright)*csq
          qzp(i,j,kc,QREINT) = rhoe + (apright + amright)*enth*csq + azeright

          ! minus state on face kc
          ! note this is different from how we do 1D, 2D, and the
          ! x and y-faces in 3D, where the analogous thing would have
          ! been to find the minus state on face kc+1
          cc = c(i,j,k3d-1)
          csq = cc**2
          rho = q(i,j,k3d-1,QRHO)
          u = q(i,j,k3d-1,QU)
          v = q(i,j,k3d-1,QV)
          w = q(i,j,k3d-1,QW)
          p = q(i,j,k3d-1,QPRES)
          rhoe = q(i,j,k3d-1,QREINT)
          enth = ( (rhoe+p)/rho )/csq

          dum    = flatn(i,j,k3d-1)*(u    - Ip(i,j,km,3,1,QU))
          dvm    = flatn(i,j,k3d-1)*(v    - Ip(i,j,km,3,1,QV))
          dwm    = flatn(i,j,k3d-1)*(w    - Ip(i,j,km,3,1,QW))
          dpm    = flatn(i,j,k3d-1)*(p    - Ip(i,j,km,3,1,QPRES))

          drho  = flatn(i,j,k3d-1)*(rho  - Ip(i,j,km,3,2,QRHO))
          du    = flatn(i,j,k3d-1)*(u    - Ip(i,j,km,3,2,QU))
          dv    = flatn(i,j,k3d-1)*(v    - Ip(i,j,km,3,2,QV))
          dw    = flatn(i,j,k3d-1)*(w    - Ip(i,j,km,3,2,QW))
          dp    = flatn(i,j,k3d-1)*(p    - Ip(i,j,km,3,2,QPRES))
          drhoe = flatn(i,j,k3d-1)*(rhoe - Ip(i,j,km,3,2,QREINT))

          dup    = flatn(i,j,k3d-1)*(u    - Ip(i,j,km,3,3,QU))
          dvp    = flatn(i,j,k3d-1)*(v    - Ip(i,j,km,3,3,QV))
          dwp    = flatn(i,j,k3d-1)*(w    - Ip(i,j,km,3,3,QW))
          dpp    = flatn(i,j,k3d-1)*(p    - Ip(i,j,km,3,3,QPRES))

          alpham = 0.5d0*(dpm/(rho*cc) - dwm)*rho/cc
          alphap = 0.5d0*(dpp/(rho*cc) + dwp)*rho/cc
          alpha0r = drho - dp/csq
          alpha0e = drhoe - dp*enth
          alpha0u = du
          alpha0v = dv

          if (w-cc .gt. 0.d0) then
             amleft = -alpham
          else if (w-cc .lt. 0.d0) then
             amleft = 0.d0
          else
             amleft = -0.5d0*alpham
          endif
          if (w+cc .gt. 0.d0) then
             apleft = -alphap
          else if (w+cc .lt. 0.d0) then
             apleft = 0.d0
          else
             apleft = -0.5d0*alphap
          endif
          if (w .gt. 0.d0) then
             azrleft = -alpha0r
             azeleft = -alpha0e
             azu1left = -alpha0u
             azv1left = -alpha0v
          else if (w .lt. 0.d0) then
             azrleft = 0.d0
             azeleft = 0.d0
             azu1left = 0.d0
             azv1left = 0.d0
          else
             azrleft = -0.5d0*alpha0r
             azeleft = -0.5d0*alpha0e
             azu1left = -0.5d0*alpha0u
             azv1left = -0.5d0*alpha0v
          endif

          qzm(i,j,kc,QRHO) = rho + apleft + amleft + azrleft
          qzm(i,j,kc,QRHO) = max(small_dens, qzm(i,j,kc,QRHO))
          qzm(i,j,kc,QW) = w + (apleft - amleft)*cc/rho
          qzm(i,j,kc,QU) = u + azu1left
          qzm(i,j,kc,QV) = v + azv1left
          qzm(i,j,kc,QPRES) = p + (apleft + amleft)*csq
          qzm(i,j,kc,QREINT) = rhoe + (apleft + amleft)*enth*csq + azeleft

       end do
    end do
    !$OMP END PARALLEL DO

    ! Now do the passively advected quantities
    !$OMP PARALLEL DO PRIVATE(iadv,n,i,j,w) IF(nadv.gt.1)
    do iadv = 1, nadv
       n = QFA + iadv - 1
       do j = ilo2-1, ihi2+1
          do i = ilo1-1, ihi1+1

             ! plus state on face kc
             w = q(i,j,k3d,QW)
             if (w .gt. 0.d0) then
                qzp(i,j,kc,n) = q(i,j,k3d,n)
             else if (w .lt. 0.d0) then
                qzp(i,j,kc,n) = q(i,j,k3d,n) &
                     + flatn(i,j,k3d)*(Im(i,j,kc,3,2,n) - q(i,j,k3d,n))
             else
                qzp(i,j,kc,n) = q(i,j,k3d,n) &
                     + 0.5d0*flatn(i,j,k3d)*(Im(i,j,kc,3,2,n) - q(i,j,k3d,n))
             endif

             ! minus state on face k
             w = q(i,j,k3d-1,QW)
             if (w .gt. 0.d0) then
                qzm(i,j,kc,n) = q(i,j,k3d-1,n) &
                     + flatn(i,j,k3d-1)*(Ip(i,j,km,3,2,n) - q(i,j,k3d-1,n))
             else if (w .lt. 0.d0) then
                qzm(i,j,kc,n) = q(i,j,k3d-1,n)
             else
                qzm(i,j,kc,n) = q(i,j,k3d-1,n) &
                     + 0.5d0*flatn(i,j,k3d-1)*(Ip(i,j,km,3,2,n) - q(i,j,k3d-1,n))
             endif

          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(ispec,ns,i,j,w) IF(nspec.gt.1)
    do ispec = 1, nspec
       ns = QFS + ispec - 1
       do j = ilo2-1, ihi2+1
          do i = ilo1-1, ihi1+1

             ! plus state on face kc
             w = q(i,j,k3d,QW)
             if (w .gt. 0.d0) then
                qzp(i,j,kc,ns) = q(i,j,k3d,ns)
             else if (w .lt. 0.d0) then
                qzp(i,j,kc,ns) = q(i,j,k3d,ns) &
                     + flatn(i,j,k3d)*(Im(i,j,kc,3,2,ns) - q(i,j,k3d,ns))
             else
                qzp(i,j,kc,ns) = q(i,j,k3d,ns) &
                     + 0.5d0*flatn(i,j,k3d)*(Im(i,j,kc,3,2,ns) - q(i,j,k3d,ns))
             endif

             ! minus state on face k
             w = q(i,j,k3d-1,QW)
             if (w .gt. 0.d0) then
                qzm(i,j,kc,ns) = q(i,j,k3d-1,ns) &
                     + flatn(i,j,k3d-1)*(Ip(i,j,km,3,2,ns) - q(i,j,k3d-1,ns))
             else if (w .lt. 0.d0) then
                qzm(i,j,kc,ns) = q(i,j,k3d-1,ns)
             else
                qzm(i,j,kc,ns) = q(i,j,k3d-1,ns) &
                     + 0.5d0*flatn(i,j,k3d-1)*(Ip(i,j,km,3,2,ns) - q(i,j,k3d-1,ns))
             endif

          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

  end subroutine tracez_ppm

end module ppm_module
