        subroutine strang_chem(nx,scal_old,scal_new,
     $                         const_src,lin_src_old,lin_src_new,
     $                         intra,tstep)

        implicit none
 
        include 'nums.fi'
        include 'sndata.fi'
        include 'nkbrn.fi'

c       Quantities passed in.
        integer nx
        real*8     scal_old(-1:nx  ,nscal)
        real*8     scal_new(-1:nx  ,nscal)
        real*8    const_src( 0:nx-1,nscal)
        real*8  lin_src_old( 0:nx-1,nscal)
        real*8  lin_src_new( 0:nx-1,nscal)
        real*8        intra( 0:nx-1,nscal)
        real*8  tstep

c       Local variables
        real*8  scal_tmp(-1:nx  ,nscal)
        integer i,n
        integer ispec

        do i = -1,nx
          do n = 1,nscal
            scal_tmp(i,n) = scal_old(i,n)
          enddo
        enddo

        call sn_transient(nx,scal_tmp,const_src,
     $                    lin_src_old,lin_src_new,tstep)

        do i = -1,nx
          do n = 1,nscal
            scal_new(i,n) = scal_tmp(i,n)
          enddo
        enddo

c       Set outflow condition.
        do n = 1,nscal
          scal_new(nx,n) = scal_new(nx-1,n)
        enddo

c       Define change in state due to chemistry.
        do i = 0,nx-1
          do n = 1,nscal
             intra(i,n) =
     $          (scal_new(i,n)-scal_old(i,n)) / tstep
     $        - const_src(i,n)
     $        - 0.5d0*(lin_src_old(i,n)+lin_src_new(i,n))
          enddo
        enddo

        print *,' '
        print *,'strang_chem time '
        print *,' '

        end

