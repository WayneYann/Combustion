module slope_module
  
  implicit none

  private

  public uslope

contains

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine uslope(q,flatn,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                        dqx,dqy,dqz,qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
                        ilo1,ilo2,ihi1,ihi2,kc,k3d,nv)

      use meth_params_module

      implicit none

      integer,intent(in):: qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
      integer,intent(in):: qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3
      integer,intent(in):: ilo1,ilo2,ihi1,ihi2,kc,k3d,nv

      double precision,intent(in ):: q(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,nv)
      double precision,intent(in ):: flatn(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
      double precision,intent(out):: dqx(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,nv)
      double precision,intent(out):: dqy(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,nv)
      double precision,intent(out):: dqz(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,nv)

      integer i, j, k, n, ilo, ihi

      double precision dlft, drgt, slop, dq1
      double precision dm, dp, dc, ds, sl, dl, dfm, dfp

      double precision, allocatable::dsgn(:,:),dlim(:,:),df(:,:),dcen(:,:)

      double precision, parameter :: four3rd = 4.d0/3.d0, sixth = 1.d0/6.d0

      ilo = MIN(ilo1,ilo2)
      ihi = MAX(ihi1,ihi2)

      allocate (dsgn(ilo-2:ihi+2,ilo-2:ihi+2))
      allocate (dlim(ilo-2:ihi+2,ilo-2:ihi+2))
      allocate (  df(ilo-2:ihi+2,ilo-2:ihi+2))
      allocate (dcen(ilo-2:ihi+2,ilo-2:ihi+2))

      if(iorder.eq.1) then

         do n = 1, nv
            do j = ilo2-1, ihi2+1
               do i = ilo1-1, ihi1+1
                  dqx(i,j,kc,n) = 0.d0
                  dqy(i,j,kc,n) = 0.d0
                  dqz(i,j,kc,n) = 0.d0
               enddo
            enddo
         enddo

      else

         do n = 1, nv 

            ! Compute slopes in first coordinate direction
            !$OMP PARALLEL DO PRIVATE(i,j,dlft,drgt,slop,dq1)
            do j = ilo2-1, ihi2+1

               ! First compute Fromm slopes
               do i = ilo1-2, ihi1+2
                  dlft = 2.d0*(q(i ,j,k3d,n) - q(i-1,j,k3d,n))
                  drgt = 2.d0*(q(i+1,j,k3d,n) - q(i ,j,k3d,n))
                  dcen(i,j) = .25d0 * (dlft+drgt)
                  dsgn(i,j) = sign(1.d0, dcen(i,j))
                  slop = min( abs(dlft), abs(drgt) )
                  if (dlft*drgt .ge. 0.d0) then
                     dlim(i,j) = slop
                  else
                     dlim(i,j) = 0.d0
                  endif
                  df(i,j) = dsgn(i,j)*min( dlim(i,j), abs(dcen(i,j)) )
               enddo

               ! Now compute limited fourth order slopes
               do i = ilo1-1, ihi1+1
                  dq1       = four3rd*dcen(i,j) - sixth*(df(i+1,j) + df(i-1,j))
                  dqx(i,j,kc,n) = flatn(i,j,k3d)*dsgn(i,j)*min(dlim(i,j),abs(dq1))
               enddo

            enddo
            !$OMP END PARALLEL DO

            ! Compute slopes in second coordinate direction
            !$OMP PARALLEL DO PRIVATE(i,j,dlft,drgt,slop,dq1)
            do i = ilo1-1, ihi1+1
               ! First compute Fromm slopes for this column
               do j = ilo2-2, ihi2+2
                  dlft = 2.d0*(q(i,j ,k3d,n) - q(i,j-1,k3d,n))
                  drgt = 2.d0*(q(i,j+1,k3d,n) - q(i,j ,k3d,n))
                  dcen(i,j) = .25d0 * (dlft+drgt)
                  dsgn(i,j) = sign( 1.d0, dcen(i,j) )
                  slop = min( abs(dlft), abs(drgt) )
                  if (dlft*drgt .ge. 0.d0) then
                     dlim(i,j) = slop
                  else
                     dlim(i,j) = 0.d0
                  endif
                  df(i,j) = dsgn(i,j)*min( dlim(i,j),abs(dcen(i,j)) )
               enddo

               ! Now compute limited fourth order slopes
               do j = ilo2-1, ihi2+1
                  dq1 = four3rd*dcen(i,j) - sixth*( df(i,j+1) + df(i,j-1) )
                  dqy(i,j,kc,n) = flatn(i,j,k3d)*dsgn(i,j)*min(dlim(i,j),abs(dq1))
               enddo
            enddo
            !$OMP END PARALLEL DO

            ! Compute slopes in third coordinate direction
            !$OMP PARALLEL DO PRIVATE(i,j,k,dm,dp,dc,ds,sl,dl,dfm,dfp,dq1)
            do j = ilo2-1, ihi2+1
               do i = ilo1-1, ihi1+1

                  ! Compute Fromm slope on slab below
                  k = k3d-1
                  dm = 2.d0*(q(i,j,k ,n) - q(i,j,k-1,n))
                  dp = 2.d0*(q(i,j,k+1,n) - q(i,j,k, n))
                  dc = .25d0*(dm+dp)
                  ds = sign( 1.d0, dc )
                  sl = min( abs(dm), abs(dp) )
                  if (dm*dp .ge. 0.d0) then
                     dl = sl
                  else
                     dl = 0.d0
                  endif
                  dfm = ds*min(dl,abs(dc))

                  ! Compute Fromm slope on slab above
                  k = k3d+1
                  dm = 2.d0*(q(i,j,k ,n) - q(i,j,k-1,n))
                  dp = 2.d0*(q(i,j,k+1,n) - q(i,j,k, n))
                  dc = .25d0*(dm+dp)
                  ds = sign( 1.d0, dc )
                  sl = min( abs(dm), abs(dp) )
                  if (dm*dp .ge. 0.d0) then
                     dl = sl
                  else
                     dl = 0.d0
                  endif
                  dfp = ds*min(dl,abs(dc))

                  ! Compute Fromm slope on current slab
                  k = k3d
                  dm = 2.d0*(q(i,j,k ,n) - q(i,j,k-1,n))
                  dp = 2.d0*(q(i,j,k+1,n) - q(i,j,k, n))
                  dc = .25d0*(dm+dp)
                  ds = sign( 1.d0, dc )
                  sl = min( abs(dm), abs(dp) )
                  if (dm*dp .ge. 0.d0) then
                     dl = sl
                  else
                     dl = 0.d0
                  endif

                  ! Now compute limited fourth order slopes
                  dq1 = four3rd*dc - sixth*( dfp + dfm )
                  dqz(i,j,kc,n) = flatn(i,j,k3d)*ds*min(dl,abs(dq1))
               enddo
            enddo
            !$OMP END PARALLEL DO
         enddo

      endif

      deallocate(dsgn,dlim,df,dcen)

      end subroutine uslope

end module slope_module
