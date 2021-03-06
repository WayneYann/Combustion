
      subroutine ca_apply_metric(lo, hi, &
                                 rhs, rl1, rh1,  &
                                 ecx, ecxl1, ecxh1, dx, coord_type)

      implicit none
      integer lo(1),hi(1)
      integer rl1, rh1
      integer ecxl1, ecxh1
      integer coord_type
      double precision rhs(rl1:rh1)
      double precision ecx(ecxl1:ecxh1)
      double precision dx(1)

      double precision r,rlo,rhi
      integer i

      ! r-z
      if (coord_type .eq. 1) then

         ! At centers
         do i=lo(1),hi(1)
            r = (dble(i)+0.5d0) * dx(1)
            rhs(i) = rhs(i) * r
         enddo

         ! On edges
         do i=lo(1),hi(1)+1
            r = dble(i)*dx(1)
            ecx(i) = ecx(i) * r
         enddo

      ! spherical
      else if (coord_type .eq. 2) then

         ! At centers
         do i=lo(1),hi(1)
!           r = (dble(i)+0.5d0) * dx(1)
!           rhs(i) = rhs(i) * r**2
            rlo = dble(i) * dx(1)
            rhi = rlo + dx(1)
            rhs(i) = rhs(i) * (rhi**3 - rlo**3) / (3.d0 * dx(1))
         enddo

         ! Note that (rhi**3 - rlo**3) / (3 dr) = ( (r + dr/2)**3 - (r - dr/2)**3 ) / (3 dr)
         !                                 = r^2 + dr^2 / 12

         ! On edges
         do i=lo(1),hi(1)+1
            r = dble(i)*dx(1)
            ecx(i) = ecx(i) * r**2
         enddo

      else 
         print *,'Bogus coord_type in apply_metric ' ,coord_type
         call bl_error("Error:: MGutils_1d.f90 :: ca_apply_metric")
      end if

      end subroutine ca_apply_metric

!-----------------------------------------------------------------------

      subroutine ca_unweight_cc(lo, hi, &
                                cc, cl1, ch1,  &
                                dx, coord_type)

      implicit none
      integer lo(1),hi(1)
      integer cl1, ch1
      integer coord_type
      double precision cc(cl1:ch1)
      double precision dx(1)

      double precision r
      integer i

      ! r-z
      if (coord_type .eq. 1) then

         do i=lo(1),hi(1)
            r = (dble(i)+0.5d0) * dx(1)
            cc(i) = cc(i) / r
         enddo

      ! spherical
      else if (coord_type .eq. 2) then

         do i=lo(1),hi(1)
            r = (dble(i)+0.5d0) * dx(1)
            cc(i) = cc(i) / r**2
         enddo

      else 
         print *,'Bogus coord_type in unweight_cc ' ,coord_type
         call bl_error("Error:: MGutils_1d.f90 :: ca_unweight_cc")
      end if

      end subroutine ca_unweight_cc

!-----------------------------------------------------------------------

      subroutine ca_unweight_edges(lo, hi, &
                                   ecx, ecxl1, ecxh1, dx, coord_type)

      implicit none
      integer lo(1),hi(1)
      integer ecxl1, ecxh1
      integer coord_type
      double precision ecx(ecxl1:ecxh1)
      double precision dx(1)

      double precision r
      integer i

      ! r-z
      if (coord_type .eq. 1) then

         ! On edges
         do i=lo(1),hi(1)+1
            r = abs(dble(i))*dx(1)
            if (i.ne.0) ecx(i) = ecx(i) / r
         enddo

      ! spherical
      else if (coord_type .eq. 2) then

         ! On edges
         do i=lo(1),hi(1)+1
            r = dble(i)*dx(1)
            if (i.ne.0) ecx(i) = ecx(i) / r**2
         enddo

      else 
         print *,'Bogus coord_type in unweight_edges ' ,coord_type
         call bl_error("Error:: MGutils_1d.f90 :: ca_unweight_edges")
      end if

      end subroutine ca_unweight_edges
