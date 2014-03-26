! ::: -----------------------------------------------------------
! ::: This routine will tag high error cells based on the density
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: tag      <=  integer tag array
! ::: lo,hi     => index extent of tag array
! ::: set       => integer value to tag cell for refinement
! ::: clear     => integer value to untag cell
! ::: den       => density array
! ::: nd        => number of components in den array (should be 1)
! ::: domlo,hi  => index extent of problem domain
! ::: delta     => cell spacing
! ::: xlo       => physical location of lower left hand
! :::              corner of tag array
! ::: problo    => phys loc of lower left corner of prob domain
! ::: time      => problem evolution time
! ::: level     => refinement level of this array
! ::: -----------------------------------------------------------
      subroutine rns_denerror(tag,tagl1,tagl2,tagl3,tagh1,tagh2,tagh3, &
                             set,clear, &
                             den,denl1,denl2,denl3,denh1,denh2,denh3, &
                             lo,hi,nd,domlo,domhi, &
                             delta,xlo,problo,time,level)
      use probdata_module
      implicit none

      integer set, clear, nd, level
      integer tagl1,tagl2,tagl3,tagh1,tagh2,tagh3
      integer denl1,denl2,denl3,denh1,denh2,denh3
      integer lo(3), hi(3), domlo(3), domhi(3)
      integer tag(tagl1:tagh1,tagl2:tagh2,tagl3:tagh3)
      double precision den(denl1:denh1,denl2:denh2,denl3:denh3,nd)
      double precision delta(3), xlo(3), problo(3), time

      double precision ax,ay,az,xcen,ycen,zcen,rtag
      integer i, j, k

      if (level .eq. 0) then
         rtag = 0.25d0*Length(1) - 2.d0*delta(1)
      else
         rtag = 0.125d0*Length(1) - 2.d0*delta(1)
      end if

      !$omp parallel private(i,j,k,ax,ay,az,xcen,ycen,zcen)

!     Tag on regions of high density
      if (level .lt. max_denerr_lev) then
         !$omp do
         do k = lo(3), hi(3)
            zcen = abs(xlo(3) + delta(3)*(dble(k-lo(3)) + 0.5d0) - center(3))
            do j = lo(2), hi(2)
               ycen = abs(xlo(2) + delta(2)*(dble(j-lo(2)) + 0.5d0) - center(2))
               do i = lo(1), hi(1)
                  xcen = abs(xlo(1) + delta(1)*(dble(i-lo(1)) + 0.5d0) - center(1))
                  if (xcen<rtag .and. ycen<rtag .and. zcen<rtag) then
                     tag(i,j,k) = set
                  endif
               enddo
            enddo
         enddo
         !$omp end do nowait
      endif

!     Tag on regions of high density gradient
      if (level .lt. max_dengrad_lev) then
         !$omp do collapse(2)
         do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               ax = ABS(den(i+1,j,k,1) - den(i,j,k,1))
               ay = ABS(den(i,j+1,k,1) - den(i,j,k,1))
               az = ABS(den(i,j,k+1,1) - den(i,j,k,1))
               ax = MAX(ax,ABS(den(i,j,k,1) - den(i-1,j,k,1)))
               ay = MAX(ay,ABS(den(i,j,k,1) - den(i,j-1,k,1)))
               az = MAX(az,ABS(den(i,j,k,1) - den(i,j,k-1,1)))
               if ( sqrt(ax*ax+ay*ay+az*az) .ge. dengrad * den(i,j,k,1) * delta(1) ) then
                  tag(i,j,k) = set
               endif
            enddo
         enddo
         enddo
         !$omp end do
      endif

      !$omp end parallel
      
      end

! ::: -----------------------------------------------------------
! ::: This routine will tag high error cells based on the temperature
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: tag      <=  integer tag array
! ::: lo,hi     => index extent of tag array
! ::: set       => integer value to tag cell for refinement
! ::: clear     => integer value to untag cell
! ::: temp      => temperature array
! ::: np        => number of components in temp array (should be 1)
! ::: domlo,hi  => index extent of problem domain
! ::: delta     => cell spacing
! ::: xlo       => physical location of lower left hand
! :::              corner of tag array
! ::: problo    => phys loc of lower left corner of prob domain
! ::: time      => problem evolution time
! ::: level     => refinement level of this array
! ::: -----------------------------------------------------------

      subroutine rns_temperror(tag,tagl1,tagl2,tagl3,tagh1,tagh2,tagh3, &
                             set,clear, &
                             temp,templ1,templ2,templ3,temph1,temph2,temph3, &
                             lo,hi,nd,domlo,domhi, &
                             delta,xlo,problo,time,level)
      use probdata_module
      implicit none

      integer set, clear, nd, level
      integer tagl1,tagl2,tagl3,tagh1,tagh2,tagh3
      integer templ1,templ2,templ3,temph1,temph2,temph3
      integer lo(3), hi(3), domlo(3), domhi(3)
      integer tag(tagl1:tagh1,tagl2:tagh2,tagl3:tagh3)
      double precision temp(templ1:temph1,templ2:temph2,templ3:temph3,nd)
      double precision delta(3), xlo(3), problo(3), time

      double precision ax,ay,az
      integer i, j, k

      !$omp parallel private(i,j,k,ax,ay,az)

!     Tag on regions of high temperature
      if (level .lt. max_temperr_lev) then
         !$omp do collapse(2)
         do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               if (temp(i,j,k,1) .ge. temperr) then
                  tag(i,j,k) = set
               endif
            enddo
         enddo
         enddo
         !$omp end do nowait
      endif

!     Tag on regions of high temperature gradient
      if (level .lt. max_tempgrad_lev) then
         !$omp do collapse(2)
         do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               ax = ABS(temp(i+1,j,k,1) - temp(i,j,k,1))
               ay = ABS(temp(i,j+1,k,1) - temp(i,j,k,1))
               az = ABS(temp(i,j,k+1,1) - temp(i,j,k,1))
               ax = MAX(ax,ABS(temp(i,j,k,1) - temp(i-1,j,k,1)))
               ay = MAX(ay,ABS(temp(i,j,k,1) - temp(i,j-1,k,1)))
               az = MAX(az,ABS(temp(i,j,k,1) - temp(i,j,k-1,1)))
               if ( MAX(ax,ay,az) .ge. tempgrad) then
                  tag(i,j,k) = set
               endif
            enddo
         enddo
         enddo
         !$omp end do
      endif

      !$omp end parallel
      
      end

! ::: -----------------------------------------------------------
! ::: This routine will tag high error cells based on the pressure
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: tag      <=  integer tag array
! ::: lo,hi     => index extent of tag array
! ::: set       => integer value to tag cell for refinement
! ::: clear     => integer value to untag cell
! ::: press     => pressure array
! ::: np        => number of components in press array (should be 1)
! ::: domlo,hi  => index extent of problem domain
! ::: delta     => cell spacing
! ::: xlo       => physical location of lower left hand
! :::              corner of tag array
! ::: problo    => phys loc of lower left corner of prob domain
! ::: time      => problem evolution time
! ::: level     => refinement level of this array
! ::: -----------------------------------------------------------
      subroutine rns_presserror(tag,tagl1,tagl2,tagl3,tagh1,tagh2,tagh3, &
                               set,clear, &
                               press,pressl1,pressl2,pressl3, &
                                     pressh1,pressh2,pressh3, &
                               lo,hi,np,domlo,domhi, &
                               delta,xlo,problo,time,level)

      use probdata_module
      implicit none

      integer set, clear, np, level
      integer tagl1,tagl2,tagl3,tagh1,tagh2,tagh3
      integer pressl1,pressl2,pressl3,pressh1,pressh2,pressh3
      integer lo(3), hi(3), domlo(3), domhi(3)
      integer tag(tagl1:tagh1,tagl2:tagh2,tagl3:tagh3)
      double precision press(pressl1:pressh1,pressl2:pressh2,pressl3:pressh3,np)
      double precision delta(3), xlo(3), problo(3), time

      double precision ax,ay,az
      integer i, j, k

      !$omp parallel private(i,j,k,ax,ay,az)

!     Tag on regions of high pressure
      if (level .lt. max_presserr_lev) then
         !$omp do collapse(2)
         do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               if (press(i,j,k,1) .ge. presserr) then
                  tag(i,j,k) = set
               endif
            enddo
         enddo
         enddo
         !$omp end do nowait
      endif

!     Tag on regions of high pressure gradient
      if (level .lt. max_pressgrad_lev) then
         !$omp do collapse(2)
         do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               ax = ABS(press(i+1,j,k,1) - press(i,j,k,1))
               ay = ABS(press(i,j+1,k,1) - press(i,j,k,1))
               az = ABS(press(i,j,k+1,1) - press(i,j,k,1))
               ax = MAX(ax,ABS(press(i,j,k,1) - press(i-1,j,k,1)))
               ay = MAX(ay,ABS(press(i,j,k,1) - press(i,j-1,k,1)))
               az = MAX(az,ABS(press(i,j,k,1) - press(i,j,k-1,1)))
               if ( MAX(ax,ay,az) .ge. pressgrad) then
                  tag(i,j,k) = set
               endif
            enddo
         enddo
         enddo
         !$omp end do
      endif

      !$omp end parallel

      end

! ::: -----------------------------------------------------------
! ::: This routine will tag high error cells based on the velocity
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: tag      <=  integer tag array
! ::: lo,hi     => index extent of tag array
! ::: set       => integer value to tag cell for refinement
! ::: clear     => integer value to untag cell
! ::: vel       => velocity array
! ::: nv        => number of components in vel array (should be 1)
! ::: domlo,hi  => index extent of problem domain
! ::: delta     => cell spacing
! ::: xlo       => physical location of lower left hand
! :::              corner of tag array
! ::: problo    => phys loc of lower left corner of prob domain
! ::: time      => problem evolution time
! ::: level     => refinement level of this array
! ::: -----------------------------------------------------------
      subroutine rns_velerror(tag,tagl1,tagl2,tagl3,tagh1,tagh2,tagh3, &
                             set,clear, &
                             vel,vell1,vell2,vell3,velh1,velh2,velh3, &
                             lo,hi,nv,domlo,domhi, &
                             delta,xlo,problo,time,level)
      use probdata_module
      implicit none

      integer set, clear, nv, level
      integer tagl1,tagl2,tagl3,tagh1,tagh2,tagh3
      integer vell1,vell2,vell3,velh1,velh2,velh3
      integer lo(3), hi(3), domlo(3), domhi(3)
      integer tag(tagl1:tagh1,tagl2:tagh2,tagl3:tagh3)
      double precision vel(vell1:velh1,vell2:velh2,vell3:velh3,nv)
      double precision delta(3), xlo(3), problo(3), time

      double precision ax,ay,az
      integer i, j, k

      !$omp parallel private(i,j,k,ax,ay,az)

!     Tag on regions of high velocity gradient
      if (level .lt. max_velgrad_lev) then
         !$omp do collapse(2)
         do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               ax = ABS(vel(i+1,j,k,1) - vel(i,j,k,1))
               ay = ABS(vel(i,j+1,k,1) - vel(i,j,k,1))
               az = ABS(vel(i,j,k+1,1) - vel(i,j,k,1))
               ax = MAX(ax,ABS(vel(i,j,k,1) - vel(i-1,j,k,1)))
               ay = MAX(ay,ABS(vel(i,j,k,1) - vel(i,j-1,k,1)))
               az = MAX(az,ABS(vel(i,j,k,1) - vel(i,j,k-1,1)))
               if ( MAX(ax,ay,az) .ge. velgrad) then
                  tag(i,j,k) = set
               endif
            enddo
         enddo
         enddo
         !$omp end do
      endif

      !$omp end parallel

      end


! ::: -----------------------------------------------------------
! ::: This routine will tag high error cells based on the vorticity
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: tag      <=  integer tag array
! ::: lo,hi     => index extent of tag array
! ::: set       => integer value to tag cell for refinement
! ::: clear     => integer value to untag cell
! ::: vort      => vorticity array
! ::: nd        => number of components in vort array (should be 1)
! ::: domlo,hi  => index extent of problem domain
! ::: delta     => cell spacing
! ::: xlo       => physical location of lower left hand
! :::              corner of tag array
! ::: problo    => phys loc of lower left corner of prob domain
! ::: time      => problem evolution time
! ::: level     => refinement level of this array
! ::: -----------------------------------------------------------
      subroutine rns_vorterror(tag,tagl1,tagl2,tagl3,tagh1,tagh2,tagh3, &
                             set,clear, &
                             vort,vortl1,vortl2,vortl3,vorth1,vorth2,vorth3, &
                             lo,hi,nd,domlo,domhi, &
                             delta,xlo,problo,time,level)
      use probdata_module
      implicit none

      integer set, clear, nd, level
      integer tagl1,tagl2,tagl3,tagh1,tagh2,tagh3
      integer vortl1,vortl2,vortl3,vorth1,vorth2,vorth3
      integer lo(3), hi(3), domlo(3), domhi(3)
      integer tag(tagl1:tagh1,tagl2:tagh2,tagl3:tagh3)
      double precision vort(vortl1:vorth1,vortl2:vorth2,vortl3:vorth3,nd)
      double precision delta(3), xlo(3), problo(3), time

      double precision ax,ay,az
      integer i, j, k

      !$omp parallel private(i,j,k,ax,ay,az)

!     Tag on regions of high vorticity
      if (level .lt. max_vorterr_lev) then
         !$omp do collapse(2)
         do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               if (abs(vort(i,j,k,1)) .ge. vorterr) then
                  tag(i,j,k) = set
               endif
            enddo
         enddo
         enddo
         !$omp end do nowait
      endif

!     Tag on regions of high vorticity gradient
      if (level .lt. max_vortgrad_lev) then
         !$omp do collapse(2)
         do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               ax = ABS(vort(i+1,j,k,1) - vort(i,j,k,1))
               ay = ABS(vort(i,j+1,k,1) - vort(i,j,k,1))
               az = ABS(vort(i,j,k+1,1) - vort(i,j,k,1))
               ax = MAX(ax,ABS(vort(i,j,k,1) - vort(i-1,j,k,1)))
               ay = MAX(ay,ABS(vort(i,j,k,1) - vort(i,j-1,k,1)))
               az = MAX(az,ABS(vort(i,j,k,1) - vort(i,j,k-1,1)))
               if ( MAX(ax,ay,az) .ge. vortgrad) then
                  tag(i,j,k) = set
               endif
            enddo
         enddo
         enddo
         !$omp end do
      endif

      !$omp end parallel
      
      end


! ::: -----------------------------------------------------------
! ::: This routine will tag high error cells based on the flame tracer
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: tag      <=  integer tag array
! ::: lo,hi     => index extent of tag array
! ::: set       => integer value to tag cell for refinement
! ::: clear     => integer value to untag cell
! ::: trac      => flame trac array
! ::: nd        => number of components in trac array (should be 1)
! ::: domlo,hi  => index extent of problem domain
! ::: delta     => cell spacing
! ::: xlo       => physical location of lower left hand
! :::              corner of tag array
! ::: problo    => phys loc of lower left corner of prob domain
! ::: time      => problem evolution time
! ::: level     => refinement level of this array
! ::: -----------------------------------------------------------
      subroutine rns_tracerror(tag,tagl1,tagl2,tagl3,tagh1,tagh2,tagh3, &
                             set,clear, &
                             trac,tracl1,tracl2,tracl3,trach1,trach2,trach3, &
                             lo,hi,nd,domlo,domhi, &
                             delta,xlo,problo,time,level)
      use probdata_module
      implicit none

      integer set, clear, nd, level
      integer tagl1,tagl2,tagl3,tagh1,tagh2,tagh3
      integer tracl1,tracl2,tracl3,trach1,trach2,trach3
      integer lo(3), hi(3), domlo(3), domhi(3)
      integer tag(tagl1:tagh1,tagl2:tagh2,tagl3:tagh3)
      double precision trac(tracl1:trach1,tracl2:trach2,tracl3:trach3,nd)
      double precision delta(3), xlo(3), problo(3), time

      integer i, j, k

!     Tag on regions of high flame tracer
      if (level .lt. max_tracerr_lev) then
         !$omp parallel do private(i,j,k) collapse(2)
         do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               if (trac(i,j,k,1) .ge. tracerr) then
                  tag(i,j,k) = set
               endif
            enddo
         enddo
         enddo
         !$omp end parallel do
      endif
      
      end
