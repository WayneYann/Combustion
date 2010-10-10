! ::: 
! ::: ------------------------------------------------------------------
! ::: 

     subroutine ca_estdt_burning(u,u_l1,u_l2,u_h1,u_h2, u_h3,u_h3, &
                                 enuc,e_l1,e_l2,e_h1,e_h2,e_h3,e_h3,&
                                 omegadot,o_l1,o_l2,o_h1,o_h2,o_h3,o_h3,lo,hi,dx,dt)

     use network, only : nspec, naux
     use eos_module
     use meth_params_module, only : NVAR, URHO, UMX, UMY, UEINT, UTEMP, UFS, UFX, &
                                    allow_negative_energy
     use eos_module 
     use burner_dt_module	

     implicit none

     integer          :: u_l1,u_l2,u_l3,u_h1,u_h2,u_h3
     integer          :: o_l1,o_l2,o_l3,o_h1,o_h2,o_h3
     integer          :: e_l1,e_l2,e_l3,e_h1,e_h2,e_h3
     integer          :: lo(3), hi(3)
     double precision ::u(u_l1:u_h1,u_l2:u_h2,u_l3:u_h3,NVAR)
     double precision ::enuc(e_l1:e_h1,e_l2:e_h2,e_l3:e_h3)
     double precision ::omegadot(o_l1:o_h1,o_l2:o_h2,o_l3:o_h3,nspec)
     double precision ::dx,dt,nucdt

     double precision :: p, gamc, c, T, dpdr, dpde, xin(nspec),xout(nspec), deltx(nspec)
     double precision ::  ein,eout,delte 
     double precision :: rhoInv,ux,uy,dt1,dt2
     integer          :: i,j,n

     logical, save    :: first = .true.
     !
     !    NOTE:  enuc contains change in enuc / dt
     !    NOTE:  omegadot_i contains only change in X_i, not scaled by dt.
     !
     do k = lo(3),hi(3)
        do j = lo(2),hi(2)
           do i = lo(1),hi(1)

              rho = u(i,j,k,URHO)
              rhoInv = 1.d0 / rho
              xin(1:nspec) = u(i,j,k,UFS:UFS+nspec-1) * rhoInv
              T = u(i,j,k,UTEMP)
              ux  = u(i,j,k,UMX)*rhoInv
              uy  = u(i,j,k,UMY)*rhoInv
              uz  = u(i,j,k,UMZ)*rhoInv
              ein  = u(i,j,k,UEINT)*rhoInv
              delte = enuc(i,j,k)*dt
	      eout = ein + delte

              deltx(1:nspec)  = omegadot(i,j,1:UFS+nspec-1) * rhoInv
              xout(1:nspec) =  xin(1:nspec)+deltx(1:nspec) 

              if (first) then
                 first = .false.
                 call ca_init_nuclear_esdt(rho,T,xin,dt)
              else
                 call ca_nuclear_esdt(xin,xout,ein,eout,dt, nucdt)
                 dt = nucdt
              endif

           enddo
        enddo
     enddo

     end subroutine ca_estdt_burning
