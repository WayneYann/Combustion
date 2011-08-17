! ::: -----------------------------------------------------------
! ::: This routine will tag high error cells based on the Laplacian.
! :::
! ::: INPUTS/OUTPUTS:
! :::
! ::: tag      <=  integer tag array
! ::: lo,hi     => index extent of tag array
! ::: set       => integer value to tag cell for refinement
! ::: clear     => integer value to untag cell
! ::: var       => array of data
! ::: nd        => number of components in var array (should be 1)
! ::: domlo,hi  => index extent of problem domain
! ::: delta     => cell spacing
! ::: xlo       => physical location of lower left hand
! :::              corner of tag array
! ::: problo    => phys loc of lower left corner of prob domain
! ::: time      => problem evolution time
! ::: level     => refinement level of this array
! ::: -----------------------------------------------------------
      subroutine ca_laplac_error(tag,tagl1,tagl2,tagl3,tagh1,tagh2,tagh3, &
                                 set,clear, &
                                 var,varl1,varl2,varl3,varh1,varh2,varh3, &
                                 lo,hi,nd,domlo,domhi, &
                                 delta,xlo,problo,time,level)
      use probdata_module
      implicit none

      integer          :: set, clear, nd, level
      integer          :: tagl1,tagl2,tagl3,tagh1,tagh2,tagh3
      integer          :: varl1,varl2,varl3,varh1,varh2,varh3
      integer          :: lo(3), hi(3), domlo(3), domhi(3)
      integer          :: tag(tagl1:tagh1,tagl2:tagh2,tagl3:tagh3)
      double precision :: var(varl1:varh1,varl2:varh2,varl3:varh3,nd)
      double precision :: delta(3), xlo(3), problo(3), time
      integer          :: i,j,k

      ! This value is  taken from FLASH
      double precision, parameter :: ctore=0.8

      do k = lo(3),hi(3)
      do j = lo(2),hi(2)
      do i = lo(1),hi(1)
         if (var(i,j,k,1) .gt. ctore) tag(i,j,k)=set
      end do
      end do
      end do

      end subroutine ca_laplac_error

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
      subroutine ca_denerror(tag,tagl1,tagl2,tagl3,tagh1,tagh2,tagh3, &
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

      double precision ax,ay,az
      integer i, j, k

!     Tag on regions of high density
      if (level .lt. max_denerr_lev) then
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  if (den(i,j,k,1) .ge. denerr) then
                     tag(i,j,k) = set
                  endif
               enddo
            enddo
         enddo
      endif

!     Tag on regions of high density gradient
      if (level .lt. max_dengrad_lev) then
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  ax = ABS(den(i+1,j,k,1) - den(i,j,k,1))
                  ay = ABS(den(i,j+1,k,1) - den(i,j,k,1))
                  az = ABS(den(i,j,k+1,1) - den(i,j,k,1))
                  ax = MAX(ax,ABS(den(i,j,k,1) - den(i-1,j,k,1)))
                  ay = MAX(ay,ABS(den(i,j,k,1) - den(i,j-1,k,1)))
                  az = MAX(az,ABS(den(i,j,k,1) - den(i,j,k-1,1)))
                  if ( MAX(ax,ay,az) .ge. dengrad) then
                     tag(i,j,k) = set
                  endif
               enddo
            enddo
         enddo
      endif

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

      subroutine ca_temperror(tag,tagl1,tagl2,tagl3,tagh1,tagh2,tagh3, &
                              set,clear, &
                              temp,templ1,templ2,templ3,temph1, &
                              temph2,temph3, &
                              lo,hi,np,domlo,domhi, &
                              delta,xlo,problo,time,level)

      use probdata_module
      implicit none

      integer set, clear, np, level
      integer tagl1,tagl2,tagl3,tagh1,tagh2,tagh3
      integer templ1,templ2,templ3,temph1,temph2,temph3
      integer lo(3), hi(3), domlo(3), domhi(3)
      integer tag(tagl1:tagh1,tagl2:tagh2,tagl3:tagh3)
      double precision temp(templ1:temph1,templ2:temph2, &
                            templ3:temph3,np)
      double precision delta(3), xlo(3), problo(3), time

      double precision ax,ay,az
      integer i, j, k

!     Tag on regions of high temperature
      if (level .lt. max_temperr_lev) then
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  if (temp(i,j,k,1) .ge. temperr) then
                     tag(i,j,k) = set
                  endif
               enddo
            enddo
         enddo
      endif

!     Tag on regions of high temperature gradient
      if (level .lt. max_tempgrad_lev) then
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
      endif

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
      subroutine ca_presserror(tag,tagl1,tagl2,tagl3,tagh1,tagh2,tagh3, &
                             set,clear, &
                             press,pressl1,pressl2,pressl3,pressh1, &
                             pressh2,pressh3, &
                             lo,hi,np,domlo,domhi, &
                             delta,xlo,problo,time,level)

      use probdata_module
      implicit none

      integer set, clear, np, level
      integer tagl1,tagl2,tagl3,tagh1,tagh2,tagh3
      integer pressl1,pressl2,pressl3,pressh1,pressh2,pressh3
      integer lo(3), hi(3), domlo(3), domhi(3)
      integer tag(tagl1:tagh1,tagl2:tagh2,tagl3:tagh3)
      double precision press(pressl1:pressh1,pressl2:pressh2, &
                             pressl3:pressh3,np)
      double precision delta(3), xlo(3), problo(3), time

      double precision ax,ay,az
      integer i, j, k

!     Tag on regions of high pressure
      if (level .lt. max_presserr_lev) then
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  if (press(i,j,k,1) .ge. presserr) then
                     tag(i,j,k) = set
                  endif
               enddo
            enddo
         enddo
      endif

!     Tag on regions of high pressure gradient
      if (level .lt. max_pressgrad_lev) then
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
      endif

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
      subroutine ca_velerror(tag,tagl1,tagl2,tagl3,tagh1,tagh2,tagh3, &
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

!     Tag on regions of high velocity gradient
      if (level .lt. max_velgrad_lev) then
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
      endif

      end


! ::: -----------------------------------------------------------
! ::: This routine will tag high error cells based on the radiation
! ::: 
! ::: INPUTS/OUTPUTS:
! ::: 
! ::: tag      <=  integer tag array
! ::: lo,hi     => index extent of tag array
! ::: set       => integer value to tag cell for refinement
! ::: clear     => integer value to untag cell
! ::: rad       => radiation array
! ::: nr        => number of components in rad array (should be 1)
! ::: domlo,hi  => index extent of problem domain
! ::: delta     => cell spacing
! ::: xlo       => physical location of lower left hand
! :::              corner of tag array
! ::: problo    => phys loc of lower left corner of prob domain
! ::: time      => problem evolution time
! ::: level     => refinement level of this array
! ::: -----------------------------------------------------------
      subroutine ca_raderror(tag,tagl1,tagl2,tagl3,tagh1,tagh2,tagh3, &
                             set,clear, &
                             rad,radl1,radl2,radl3,radh1,radh2,radh3, &
                             lo,hi,nr,domlo,domhi, &
                             delta,xlo,problo,time,level)

      use probdata_module
      implicit none

      integer set, clear, nr, level
      integer tagl1,tagl2,tagl3,tagh1,tagh2,tagh3
      integer radl1,radl2,radl3,radh1,radh2,radh3
      integer lo(3), hi(3), domlo(3), domhi(3)
      integer tag(tagl1:tagh1,tagl2:tagh2,tagl3:tagh3)
      double precision rad(radl1:radh1,radl2:radh2,radl3:radh3,nr)
      double precision delta(3), xlo(3), problo(3), time

      double precision ax,ay,az
      integer i, j, k

!     Tag on regions of high radiation
      if (level .lt. max_raderr_lev) then
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  if (rad(i,j,k,1) .ge. raderr) then
                     tag(i,j,k) = set
                  endif
               enddo
            enddo
         enddo
      endif

!     Tag on regions of high radiation gradient
      if (level .lt. max_radgrad_lev) then
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  ax = ABS(rad(i+1,j,k,1) - rad(i,j,k,1))
                  ay = ABS(rad(i,j+1,k,1) - rad(i,j,k,1))
                  az = ABS(rad(i,j,k+1,1) - rad(i,j,k,1))
                  ax = MAX(ax,ABS(rad(i,j,k,1) - rad(i-1,j,k,1)))
                  ay = MAX(ay,ABS(rad(i,j,k,1) - rad(i,j-1,k,1)))
                  az = MAX(az,ABS(rad(i,j,k,1) - rad(i,j,k-1,1)))
                  if ( MAX(ax,ay,az) .ge. radgrad) then
                     tag(i,j,k) = set
                  endif
               enddo
            enddo
         enddo
      endif

      end
