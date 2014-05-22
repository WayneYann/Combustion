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
      subroutine rns_denerror(tag,tagl1,tagl2,tagh1,tagh2, &
                             set,clear, &
                             den,denl1,denl2,denh1,denh2, &
                             lo,hi,nd,domlo,domhi, &
                             delta,xlo,problo,time,level)
      use probdata_module
      implicit none

      integer set, clear, nd, level
      integer tagl1,tagl2,tagh1,tagh2
      integer denl1,denl2,denh1,denh2
      integer lo(2), hi(2), domlo(2), domhi(2)
      integer tag(tagl1:tagh1,tagl2:tagh2)
      double precision den(denl1:denh1,denl2:denh2,nd)
      double precision delta(2), xlo(2), problo(2), time

      double precision ax,ay,xcen,ycen,rtag
      integer i, j

      if (level .eq. 0) then
         rtag = 0.25d0*Length(1) - 2.d0*delta(1)
      else
         rtag = 0.125d0*Length(1) - 2.d0*delta(1)
      end if

      !$omp parallel private(i,j,ax,ay,xcen,ycen)

!     Tag on regions of high density
      if (level .lt. max_denerr_lev) then
         !$omp do
         do j = lo(2), hi(2)
            ycen = abs(xlo(2) + delta(2)*(dble(j-lo(2)) + 0.5d0) - center(2))
            do i = lo(1), hi(1)
               xcen = abs(xlo(1) + delta(1)*(dble(i-lo(1)) + 0.5d0) - center(1))
               if (xcen<rtag .and. ycen<rtag) then
                  tag(i,j) = set
               endif
            enddo
         enddo
         !$omp end do nowait
      endif

!     Tag on regions of high density gradient
      if (level .lt. max_dengrad_lev) then
         !$omp do 
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               ax = ABS(den(i+1,j,1) - den(i,j,1))
               ay = ABS(den(i,j+1,1) - den(i,j,1))
               ax = MAX(ax,ABS(den(i,j,1) - den(i-1,j,1)))
               ay = MAX(ay,ABS(den(i,j,1) - den(i,j-1,1)))
               if ( MAX(ax,ay) .ge. dengrad) then
                  tag(i,j) = set
               endif
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

      subroutine rns_temperror(tag,tagl1,tagl2,tagh1,tagh2, &
                             set,clear, &
                             temp,templ1,templ2,temph1,temph2, &
                             lo,hi,nd,domlo,domhi, &
                             delta,xlo,problo,time,level)
      use probdata_module
      implicit none

      integer set, clear, nd, level
      integer tagl1,tagl2,tagh1,tagh2
      integer templ1,templ2,temph1,temph2
      integer lo(2), hi(2), domlo(2), domhi(2)
      integer tag(tagl1:tagh1,tagl2:tagh2)
      double precision temp(templ1:temph1,templ2:temph2,nd)
      double precision delta(2), xlo(2), problo(2), time

      double precision ax,ay
      integer i, j

      !$omp parallel private(i,j,ax,ay)

!     Tag on regions of high temperature
      if (level .lt. max_temperr_lev) then
         !$omp do 
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               if (temp(i,j,1) .ge. temperr) then
                  tag(i,j) = set
               endif
            enddo
         enddo
         !$omp end do nowait
      endif

!     Tag on regions of high temperature gradient
      if (level .lt. max_tempgrad_lev) then
         !$omp do 
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               ax = ABS(temp(i+1,j,1) - temp(i,j,1))
               ay = ABS(temp(i,j+1,1) - temp(i,j,1))
               ax = MAX(ax,ABS(temp(i,j,1) - temp(i-1,j,1)))
               ay = MAX(ay,ABS(temp(i,j,1) - temp(i,j-1,1)))
               if ( MAX(ax,ay) .ge. tempgrad) then
                  tag(i,j) = set
               endif
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
      subroutine rns_presserror(tag,tagl1,tagl2,tagh1,tagh2, &
                               set,clear, &
                               press,pressl1,pressl2, &
                                     pressh1,pressh2, &
                               lo,hi,np,domlo,domhi, &
                               delta,xlo,problo,time,level)

      use probdata_module
      implicit none

      integer set, clear, np, level
      integer tagl1,tagl2,tagh1,tagh2
      integer pressl1,pressl2,pressh1,pressh2
      integer lo(2), hi(2), domlo(2), domhi(2)
      integer tag(tagl1:tagh1,tagl2:tagh2)
      double precision press(pressl1:pressh1,pressl2:pressh2,np)
      double precision delta(2), xlo(2), problo(2), time

      double precision ax,ay
      integer i, j

      !$omp parallel private(i,j,ax,ay)

!     Tag on regions of high pressure
      if (level .lt. max_presserr_lev) then
         !$omp do 
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               if (press(i,j,1) .ge. presserr) then
                  tag(i,j) = set
               endif
            enddo
         enddo
         !$omp end do nowait
      endif

!     Tag on regions of high pressure gradient
      if (level .lt. max_pressgrad_lev) then
         !$omp do 
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               ax = ABS(press(i+1,j,1) - press(i,j,1))
               ay = ABS(press(i,j+1,1) - press(i,j,1))
               ax = MAX(ax,ABS(press(i,j,1) - press(i-1,j,1)))
               ay = MAX(ay,ABS(press(i,j,1) - press(i,j-1,1)))
               if ( MAX(ax,ay) .ge. pressgrad) then
                  tag(i,j) = set
               endif
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
      subroutine rns_velerror(tag,tagl1,tagl2,tagh1,tagh2, &
                             set,clear, &
                             vel,vell1,vell2,velh1,velh2, &
                             lo,hi,nv,domlo,domhi, &
                             delta,xlo,problo,time,level)
      use probdata_module
      implicit none

      integer set, clear, nv, level
      integer tagl1,tagl2,tagh1,tagh2
      integer vell1,vell2,velh1,velh2
      integer lo(2), hi(2), domlo(2), domhi(2)
      integer tag(tagl1:tagh1,tagl2:tagh2)
      double precision vel(vell1:velh1,vell2:velh2,nv)
      double precision delta(2), xlo(2), problo(2), time

      double precision ax,ay
      integer i, j

      !$omp parallel private(i,j,ax,ay)

!     Tag on regions of high velocity gradient
      if (level .lt. max_velgrad_lev) then
         !$omp do 
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               ax = ABS(vel(i+1,j,1) - vel(i,j,1))
               ay = ABS(vel(i,j+1,1) - vel(i,j,1))
               ax = MAX(ax,ABS(vel(i,j,1) - vel(i-1,j,1)))
               ay = MAX(ay,ABS(vel(i,j,1) - vel(i,j-1,1)))
               if ( MAX(ax,ay) .ge. velgrad) then
                  tag(i,j) = set
               endif
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
      subroutine rns_vorterror(tag,tagl1,tagl2,tagh1,tagh2, &
                             set,clear, &
                             vort,vortl1,vortl2,vorth1,vorth2, &
                             lo,hi,nd,domlo,domhi, &
                             delta,xlo,problo,time,level)
      use probdata_module
      implicit none

      integer set, clear, nd, level
      integer tagl1,tagl2,tagh1,tagh2
      integer vortl1,vortl2,vorth1,vorth2
      integer lo(2), hi(2), domlo(2), domhi(2)
      integer tag(tagl1:tagh1,tagl2:tagh2)
      double precision vort(vortl1:vorth1,vortl2:vorth2,nd)
      double precision delta(2), xlo(2), problo(2), time

      double precision ax,ay
      integer i, j

      !$omp parallel private(i,j,ax,ay)

!     Tag on regions of high vorticity
      if (level .lt. max_vorterr_lev) then
         !$omp do 
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               if (abs(vort(i,j,1)) .ge. vorterr) then
                  tag(i,j) = set
               endif
            enddo
         enddo
         !$omp end do nowait
      endif

!     Tag on regions of high vorticity gradient
      if (level .lt. max_vortgrad_lev) then
         !$omp do 
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               ax = ABS(vort(i+1,j,1) - vort(i,j,1))
               ay = ABS(vort(i,j+1,1) - vort(i,j,1))
               ax = MAX(ax,ABS(vort(i,j,1) - vort(i-1,j,1)))
               ay = MAX(ay,ABS(vort(i,j,1) - vort(i,j-1,1)))
               if ( MAX(ax,ay) .ge. vortgrad) then
                  tag(i,j) = set
               endif
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
      subroutine rns_tracerror(tag,tagl1,tagl2,tagh1,tagh2, &
                             set,clear, &
                             trac,tracl1,tracl2,trach1,trach2, &
                             lo,hi,nd,domlo,domhi, &
                             delta,xlo,problo,time,level)
      use probdata_module
      implicit none

      integer set, clear, nd, level
      integer tagl1,tagl2,tagh1,tagh2
      integer tracl1,tracl2,trach1,trach2
      integer lo(2), hi(2), domlo(2), domhi(2)
      integer tag(tagl1:tagh1,tagl2:tagh2)
      double precision trac(tracl1:trach1,tracl2:trach2,nd)
      double precision delta(2), xlo(2), problo(2), time

      integer i, j

!     Tag on regions of high flame tracer
      if (level .lt. max_tracerr_lev) then
         !$omp parallel do private(i,j) 
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               if (trac(i,j,1) .ge. tracerr) then
                  tag(i,j) = set
               endif
            enddo
         enddo
         !$omp end parallel do
      endif
      
      end
