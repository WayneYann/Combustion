      subroutine advance(nx,vel_old,vel_new,scal_old,scal_new,
     $                   Ydot_old,Ydot_new,press_old,press_new,
     $                   divu_old,divu_new,dsdt,beta_old,beta_new,
     $                   intra,dx,dt,time)

      implicit none
      include 'spec.h'
      integer nx
      real*8   vel_new(-1:nx  )
      real*8   vel_old(-1:nx  )
      real*8  scal_new(-1:nx  ,nscal)
      real*8  scal_old(-1:nx  ,nscal)
      real*8 press_new(0 :nx  )
      real*8 press_old(0 :nx  )
      real*8  Ydot_new(0 :nx-1,nspec)
      real*8  Ydot_old(0 :nx-1,nspec)
      real*8    macvel(0 :nx  )
      real*8   veledge(0 :nx  ,nscal)
      real*8      aofs(0 :nx-1,nscal)
      real*8  divu_old(0 :nx-1)
      real*8  divu_new(0 :nx-1)
      real*8  divu_tmp(0 :nx-1)
      real*8  beta_old(-1:nx,nscal)
      real*8  beta_new(-1:nx,nscal)
      real*8        mu(-1:nx)
      real*8      dsdt(0 :nx-1)
      real*8        gp(0 :nx-1)
      real*8   rhohalf(0 :nx-1)
      real*8    tforce(0 :nx-1,nscal)
      real*8      visc(0 :nx-1)
      real*8     intra(0 :nx-1,nscal)
      real*8        cp(0 :nx-1)
      real*8 x
      real*8 dx
      real*8 dt
      real*8 time
      real*8 be_cn_theta
      
      real*8    diff_old(0 :nx-1,nscal)
      real*8    diff_hat(0 :nx-1,nscal)
      real*8    diff_new(0 :nx-1,nscal)
      real*8   const_src(0 :nx-1,nscal)
      real*8 lin_src_old(0 :nx-1,nscal)
      real*8 lin_src_new(0 :nx-1,nscal)
      real*8    scal_tmp(0 :nx-1,nscal)
      real*8    sumh
      integer misdc
      
      integer i,n,ispec
      integer iunit
      
      real*8 divu_max
      real*8    alpha(0:nx-1)
      real*8      Rhs(0:nx-1)
      real*8    rhort(-1:nx  )
      real*8    Ydot_max
      
      print *,'advance: at start of time step'
c     
c*****************************************************************
c     Create MAC velocities.
c*****************************************************************
c     
      
      do i = 0,nx-1
         gp(i) = (press_old(i+1) - press_old(i)) / dx
      enddo
      
      call minmax_scal(nx,scal_old)
      print *,'... predict edge velocities'
      call pre_mac_predict(nx,vel_old,scal_old,gp,
     $                     macvel,dx,dt)
      
      call compute_rhort(nx,scal_old,rhort)

      do i = 0,nx-1
         divu_tmp(i) = divu_old(i) + 0.5d0 * dt * dsdt(i)
      enddo
      divu_max = ABS(divu_tmp(0))
      do i = 1,nx-1
         divu_max = MAX(divu_max,ABS(divu_tmp(i)))
      enddo
      print *,'DIVU norm old = ',divu_max 
      call add_dpdt(nx,rhort,divu_tmp,macvel,dx,dt)
      divu_max = ABS(divu_tmp(0))
      do i = 1,nx-1
         divu_max = MAX(divu_max,ABS(divu_tmp(i)))
      enddo
      print *,'DIVU norm new = ',divu_max 
      call macproj(nx,macvel,divu_tmp,dx)

      do n = 1,nscal
         do i = 0,nx-1
            tforce(i,n) = 0.d0
c           intra(i,n) = 0.d0
            diff_old(i,n) = 0.d0
            diff_hat(i,n) = 0.d0
            diff_new(i,n) = 0.d0
            aofs(i,n)     = 0.d0
         enddo
      enddo
c
c*****************************************************************
c     
      print *,'... creating the diffusive terms with old data'
      
c     First create the del dot lambda grad H terms 
      call calc_diffusivities(nx,scal_old,beta_old,mu)
      call rhoh_visc_terms(nx,scal_old,beta_old,visc,dx)
      do i = 0,nx-1
         diff_old(i,RhoH) = visc(i)
      enddo
      
      call compute_cp(nx,cp,scal_old)
      do i = 0,nx-1
         diff_old(i,Temp) = diff_old(i,RhoH) / 
     $        (cp(i)*scal_old(i,Density))
      enddo
      
c*****************************************************************
      
      print *,'... computing aofs with source = diff_old + intra'
      do n = 1,nscal
         do i = 0,nx-1
            tforce(i,n) = diff_old(i,n) + intra(i,n)
         enddo
      enddo
      
      call scal_aofs(nx,scal_old,macvel,aofs,tforce,dx,dt)
      
      do n = 1,nscal
         do i = 0,nx-1
            tforce(i,n) = intra(i,n)
         enddo
      enddo
      
c*****************************************************************
      
c     print *,'... update species and rho with advective terms only '
c     call update_spec(nx,scal_old,scal_new,aofs,tforce,dx,dt)
      
c*****************************************************************
      
c     print *,'... update to temp. to define new diff coeffs'
c     call update_temp(nx,scal_old,scal_new,aofs,
c     $                   alpha,beta_old,Rhs,dx,dt,be_cn_theta)
c     call cn_solve(nx,scal_new,alpha,beta_old,Rhs,
c     $                dx,dt,Temp,be_cn_theta)
      
c     call calc_diffusivities(nx,scal_new,beta_new)
      
      print *,'... set new beta = old beta'
      do n = 1,nscal
         do i = -1,nx
            beta_new(i,n) = beta_old(i,n)
         enddo
      enddo
      
c*****************************************************************

        print *,'... update to rhoH with new diff. coeffs'
        call update_rhoh(nx,scal_old,scal_new,beta_old,beta_new,
     $                   aofs,tforce,alpha,Rhs,dx,dt,be_cn_theta)
        call cn_solve(nx,scal_new,alpha,beta_new,Rhs,
     $                dx,dt,RhoH,be_cn_theta)

        print *,'... create new diff. terms for RhoH: diff_new'
        if (be_cn_theta .eq. 1.d0) then
         do i = 0,nx-1
           diff_new(i,RhoH) = 
     $               (scal_new(i,RhoH)-scal_old(i,RhoH))/dt 
     $               -aofs(i,RhoH)-tforce(i,RhoH)
         enddo
        else if (be_cn_theta .eq. 0.5d0) then
         do i = 0,nx-1
           diff_new(i,RhoH) = 2.d0 * (
     $               (scal_new(i,RhoH)-scal_old(i,RhoH))/dt 
     $               -aofs(i,RhoH)-tforce(i,RhoH)
     $               -0.5d0*diff_old(i,RhoH))
         enddo
        else 
          print *,'OOPS - BOGUS BE_CN_THETA ',be_cn_theta
          stop
        endif

        print *,'... create new temp from new RhoH'
        call rhoh_to_temp(nx,scal_new)

        print *,'... create new diff. terms for Temp: diff_new'
        call compute_cp(nx,cp,scal_new)
        do i = 0,nx-1
          diff_new(i,Temp) = diff_new(i,RhoH)/
     $                       (cp(i)*scal_new(i,Density))
        enddo

c*****************************************************************
c       Do chemistry update with const source = adv + diff_new
c*****************************************************************

        print *,'... advancing chem with const source'
        do n = 1,nscal
          do i = 0,nx-1
              const_src(i,n) =     aofs(i,n)
            lin_src_old(i,n) = diff_old(i,n)
            lin_src_new(i,n) = diff_new(i,n)
          enddo
        enddo

        call strang_chem(nx,scal_old,scal_new,
     $                   const_src,lin_src_old,lin_src_new,
     $                   intra,dt)

        print *,'MOVING ON TO MISDC LOOP '

c*****************************************************************
c       Now begin iterations.
c*****************************************************************
   
        do misdc = 1,4

c*****************************************************************

        do n=1,nscal
        do i=0,nx-1
          scal_tmp(i,n) = scal_new(i,n)
        enddo
        enddo

        call calc_diffusivities(nx,scal_new,beta_new)

        print *,'... create new diff. terms : diff_hat'
        call rhoh_visc_terms(nx,scal_new,beta_new,visc,dx)
        do i = 0,nx-1
          diff_hat(i,RhoH) = visc(i)
        enddo
        call compute_cp(nx,cp,scal_new)
        do i = 0,nx-1
          diff_hat(i,Temp) = diff_hat(i,RhoH) / 
     $                       (cp(i)*scal_new(i,Density))
        enddo

c*****************************************************************

        do n = 1,nscal
         do i = 0,nx-1
           tforce(i,n) = intra(i,n) + diff_old(i,n)
         enddo
        enddo

        print *,'... computing aofs with source = diff_old + intra'
        call scal_aofs(nx,scal_old,macvel,aofs,tforce,dx,dt)

c*****************************************************************

        do n = 1,nscal
          do i = 0,nx-1
            tforce(i,n) = intra(i,n)
     $                  + 0.5d0 * (diff_old(i,n) + diff_hat(i,n))
     $                  -          diff_hat(i,n)
          enddo
        enddo

        print *,'Updating species and rho with advective + intra terms'
        call update_spec(nx,scal_old,scal_new,aofs,tforce,dx,dt)

        print *,'... update to rhoH with new diff. coeffs'
        call update_rhoh(nx,scal_old,scal_new,beta_old,beta_new,
     $                   aofs,tforce,alpha,Rhs,dx,dt,be_cn_theta)
        call cn_solve(nx,scal_new,alpha,beta_new,Rhs,
     $                dx,dt,RhoH,be_cn_theta)

        print *,'... create new diff. terms for RhoH: diff_new'
        if (be_cn_theta .eq. 1.d0) then
         do i = 0,nx-1
           diff_new(i,RhoH) = 
     $               (scal_new(i,RhoH)-scal_old(i,RhoH))/dt 
     $               -aofs(i,RhoH)-tforce(i,RhoH)
         enddo
        else if (be_cn_theta .eq. 0.5d0) then
         do i = 0,nx-1
           diff_new(i,RhoH) = 2.d0 * (
     $               (scal_new(i,RhoH)-scal_old(i,RhoH))/dt 
     $               -aofs(i,RhoH)-tforce(i,RhoH)
     $               -0.5d0*diff_old(i,RhoH))
         enddo
        else 
          print *,'OOPS - BOGUS BE_CN_THETA ',be_cn_theta
          stop
        endif

        print *,'... create new temp from new rhoH'
        call rhoh_to_temp(nx,scal_new)

        print *,'... create new diff. terms for Temp: diff_new'
        call compute_cp(nx,cp,scal_new)
        do i = 0,nx-1
          diff_new(i,Temp) = diff_new(i,RhoH) / 
     $                       (cp(i)*scal_new(i,Density))
        enddo

c*****************************************************************
c       Do chemistry update with const source = adv + diff_new
c*****************************************************************

        print *,'... advancing chem with const and linear sources'
        do n = 1,nscal
          do i = 0,nx-1

              const_src(i,n) = aofs(i,n)
     $                       + diff_new(i,n) - diff_hat(i,n)
            lin_src_old(i,n) = diff_old(i,n)
            lin_src_new(i,n) = diff_hat(i,n)

          enddo
        enddo

        call strang_chem(nx,scal_old,scal_new,
     $                   const_src,lin_src_old,lin_src_new,
     $                   intra,dt)

        sumh = 0.d0
        do i=0,nx-1
           sumh = sumh+abs(scal_new(i,RhoH)-scal_tmp(i,RhoH))
        enddo
        print *,'SUMH ',misdc,sumh

        enddo

c*****************************************************************
c       Now end iterations and go on to update velocity.
c*****************************************************************
        
        do i = 0,nx-1
        do n = 1,nspec
          ispec = FirstSpec-1+n
          Ydot_new(i,n) = intra(i,ispec)
        enddo
        enddo

c       HACKING THIS OUT
c       call calc_diffusivities(nx,scal_new,beta_new)
        call calc_divu(nx,scal_new,beta_new,Ydot_new,divu_new,dx)

        do i = 0,nx-1
          rhohalf(i) = 0.5d0 * (scal_old(i,Density)+scal_new(i,Density))
        enddo
      
        Ydot_max = 0.d0
        do i = 0,nx-1
          dsdt(i) = (divu_new(i) - divu_old(i)) / dt
          Ydot_max = 
     $      max(abs(Ydot_new(i,1))/scal_new(i,Density),Ydot_max)
        enddo
        print *,'YDOT MAX ',Ydot_max

        print *,'... update velocities'
        call vel_edge_states(nx,vel_old,scal_old,gp,
     $                       macvel,veledge,dx,dt)
        call update_vel(nx,vel_old,vel_new,gp,rhohalf,
     $                  macvel,veledge,dx,dt)

        print *,'...project: levels = 0 0'
        call project(nx,vel_old,vel_new,rhohalf,divu_new,
     $               press_old,press_new,dx,dt)

        sumh = 0.d0
        do i=0,nx-1
           sumh = sumh+abs(scal_new(i,RhoH)-scal_old(i,RhoH))
        enddo
        print *,'SUMON ',sumh
        end
