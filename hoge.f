
      double precision function ranran(ihaji)
ccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision(a-h,o-z)
      implicit integer(n)
      include 'inc/parameternano3.dat'
      include 'inc/variable4.dat'
ccccccccccccccccccccccccccccccccccccccccccccccc
      lconst = 843314861
      nconst = 453816693
      maxcont = 1073741824
      mucnt = 2147483647
      idumy = lconst*ihaji+nconst
      idumy = mod(idumy,mucnt)
      iransu = idumy
      rrout = dble(idumy/mucnt)
      ranran = rrout
      end function

      implicit double precision(a-h,o-z)
      implicit integer(n)
      nano=9999
      bb1 = ranran(nano)
      end