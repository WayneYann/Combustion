      subroutine dscal (n,da,dx,incx)
c
c     scales a vector by a constant.
c     uses unrolled loops for increment equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c
      double precision da,dx(*)
      integer i,incx,m,mp1,n,nincx
c
      if ( n.le.0 .or. incx.le.0 ) return
      if (incx.eq.1) goto 20
c
c     code for increment not equal to 1
c
      nincx = n*incx
      do i = 1,nincx,incx
        dx(i) = da*dx(i)
      enddo
      return
c
c     code for increment equal to 1
c
c
c     clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do i = 1,m
        dx(i) = da*dx(i)
      enddo
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do i = mp1,n,5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
      enddo
      end


      subroutine dcopy(n,dx,incx,dy,incy)
c
c     copies a vector, x, to a vector, y.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(*),dy(*)
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if (n.le.0) return
      if (incx.eq.1.and.incy.eq.1) goto 20
c
c     code for unequal increments or equal increments not equal to 1
c
      ix = 1
      iy = 1
      if (incx.lt.0) ix = (-n+1)*incx + 1
      if (incy.lt.0) iy = (-n+1)*incy + 1
      do i = 1,n
        dy(iy) = dx(ix)
        ix = ix + incx
        iy = iy + incy
      enddo
      return
c
c     code for both increments equal to 1
c
c
c     clean-up loop
c
   20 m = mod(n,7)
      if( m .eq. 0 ) goto 40
      do i = 1,m
        dy(i) = dx(i)
      enddo
      if( n .lt. 7 ) return
   40 mp1 = m + 1
      do i = mp1,n,7
        dy(i) = dx(i)
        dy(i + 1) = dx(i + 1)
        dy(i + 2) = dx(i + 2)
        dy(i + 3) = dx(i + 3)
        dy(i + 4) = dx(i + 4)
        dy(i + 5) = dx(i + 5)
        dy(i + 6) = dx(i + 6)
      enddo
      end


      double precision function vddot (n,dx,incx,dy,incy)
c
c     forms the dot product of two vectors.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      vddot = 0.0d0
      dtemp = 0.0d0
      if (n.le.0) return
      if (incx.eq.1.and.incy.eq.1) goto 20
c
c     code for unequal increments or equal increments not equal to 1
c
      ix = 1
      iy = 1
      if (incx.lt.0) ix = (-n+1)*incx + 1
      if (incy.lt.0) iy = (-n+1)*incy + 1
      do i = 1,n
        dtemp = dtemp + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
      enddo
      vddot = dtemp
      return
c
c     code for both increments equal to 1
c
c
c     clean-up loop
c
   20 m = mod(n,5)
      if ( m .eq. 0 ) goto 40
      do i = 1,m
        dtemp = dtemp + dx(i)*dy(i)
      enddo
      if( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do i = mp1,n,5
        dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) +
     *   dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
      enddo
   60 vddot = dtemp
      end

