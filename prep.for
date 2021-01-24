      program prep
      implicit none

      character*80 string
      integer i, ifrag, istop, np

      integer mxfrgg
      parameter (mxfrgg = 27)

      double precision x, y, z, chg(mxfrgg)

      double precision tobohr
      parameter (tobohr = 1.0d0/0.52917724924d0)


      open(unit=5,file='MMpart_x_y_z_charge',form='formatted')

      i     = 1
      ifrag = 0

      do
        istop = 0
        read(unit=5,fmt='(A80)',iostat=istop) string
        if (istop .eq. 0) then
          read(string,*,iostat=istop) x,y,z,chg(i)
        end if

        if (istop .eq. 0) then
          if (i .eq. 1) then
            ifrag = ifrag + 1
            write(*,fmt='('' $FRG'',I3.3)') ifrag
            write(*,*) ' fragment #', ifrag
            write(*,fmt='('' COORDINATES'')')
            write(*,fmt='('' X-ax '',f20.12,'' 0 0 1 0'')') tobohr
            write(*,fmt='('' Y-ax 0 '',f20.12,'' 0 1 0'')') tobohr
            write(*,fmt='('' Z-ax 0 0 '',f20.12,'' 1 0'')') tobohr
          end if

          write(*,fmt='('' CHG'',I3.3,3F20.12,'' 1 0'')') i,
     *          x*tobohr,y*tobohr,z*tobohr
        end if

        if (istop .eq. 0) i = i + 1

        if ((i .gt. mxfrgg) .or. (istop .ne. 0)) then
          np = i - 1
          write(*,fmt='('' STOP'')')
          write(*,fmt='('' MONOPOLES'')')
          do i=1,np
            write(*,fmt='('' CHG'',I3.3,F20.12)') i,chg(i)
          end do
          write(*,fmt='('' STOP'')')
          write(*,fmt='('' $END'')')
          i = 1
        end if

        if (istop .ne. 0) exit
      end do

      if (ifrag .ne. 0) then
        write(*,fmt='('' $EFRAG'')')
        write(*,fmt='('' COORD=CART'')')
        do i=1,ifrag
          write(*,fmt='('' FRAGNAME=FRG'',I3.3)') i
          write(*,fmt='('' X-ax 1 0 0'')')
          write(*,fmt='('' Y-ax 0 1 0'')')
          write(*,fmt='('' Z-ax 0 0 1'')')
        end do
        write(*,fmt='('' $END'')')
        write(*,fmt='('' $FRGRPL'')')
        write(*,fmt='('' $END'')')
      end if
      

      return
      end
