!!!!!!!!!!!!!!!!!!!!!!!!!温度分布と密度分布の導入!!!!!!!!!!!!!!!!!!!!!!!!!!
! example.f
! ナノ粒子の初期状態を分子動力学法で作成する．
! 必要ファイル:parate.dat，random0.dat,random1.dat

! 倍精度計算の宣言
      implicit double precision(a-h,o-z)
      implicit integer(n)

! パラメータファイルparate.datを読み込む
      include 'inc/parameternano3.dat'
! commonファイルを読み込む
      include 'inc/variable4.dat'
! 読み込み用乱数ファイル
      open(1,file='inc/random1.dat')
      open(2,file='inc/random0.dat')

! 各分子の位置データの出力
      open(3,file='result/nano3/posit_pt9.dat')
      open(4,file='result/nano3/posit_ar9.dat')
!　各分子の速度データの出力
      open(5,file='result/nano3/veloc_pt9.dat')
      open(6,file='result/nano3/veloc_ar9.dat')
! 各分子のエネルギーデータの出力
      open(7,file='result/nano3/energy_pt9.dat')
      open(8,file='result/nano3/energy_ar9.dat')
      open(9,file='result/nano3/total energy9.dat')
! 各分子の最終的な位置データの出力
      open(10,file='result/nano3/pos_vel9.dat')
! 全粒子の最大速度の出力
      open(12,file='result/nano3/max_velocity9.dat')
! 系の周期長さ
      open(13,file='result/nano3/syuuki9.dat')
! Ptの層ごとの温度データ
      open(42,file='result/nano3/temp_pt_bot9.dat')
      open(43,file='result/nano3/temp_pt_top9.dat')
      open(700,file='result/nano3/nhos.dat')
! Arの温度データ
      open(14,file='result/nano3/temp_ar9.dat')
! Maxwell分布データ
      open(20,file='result/nano3/maxwell9.dat')

! 分布の出力
      open(40,file='result/nano3/temp_bunpu9.dat')
      open(41,file='result/nano3/dens_bunpu9.dat')
      open(32,file='result/nano3/count9.dat')
      open(70,file='result/nano3/density_2d.dat')
! 圧力の測定  
      open(31,file='result/nano3/pressure9.dat')
      open(33,file='result/nano3/ave_pressure9.dat')
! 熱流束の計算
      open(50,file='result/nano3/heatflux_interface.dat')
      open(51,file='result/nano3/heatflux_langevin.dat')
      open(500,file='result/nano3/heatflux_test.dat')
      open(501,file='result/nano3/heatfluxtop_test.dat')
      open(600,file='result/nano3/nano-surface.dat')

! DDLの計算
      open(60,file='result/nano3/DDL.dat')
!	  open(600,file='result/nano3/DDL-test.dat')
      open(200,file='result/nano3/DDL-hisekibun.dat')
      open(210,file='result/nano3/posz.dat')

! pvch用のファイル
      ! open(16,file='animation/nano3/color.dat')
      open(17,file='animation/nano3/pos.dat')
      open(18,file='animation/nano3/mask.dat')
      open(19,file='animation/nano3/ball.dat')

      nowstp = 0
      iransu = 9999
      ncount = 0
      ncount_pt = 0
      ncountd = 0
      n_dens_tot = 0
      tot_temp_ar = 0
      tot_zdensi = 0
      tot_zdensi_pt = 0
      tot_q_bot = 0
      tot_q_top = 0
      ddl = 0
      tot_fff = 0
      tot_fave_top_z = 0
      tot_fave_bot_z = 0
      pt_bot2 = 0
      pt_bot3 = 0
      pt_bot4 = 0
      pt_bot5 = 0
      pt_bot6 = 0
      pt_bot7 = 0
      pt_bot8 = 0
      pt_bot9 = 0
      pt_bot10 = 0
      pt_bot11 = 0
      pt_bot12 = 0
      pt_top2 = 0
      pt_top3 = 0
      pt_top4 = 0
      ! tot_bot2 = 0
      !	tot_bot3 = 0
      !	tot_bot4 = 0
      !	tot_top2 = 0
      !	tot_top3 = 0
      !	tot_top4 = 0

! 各分子の初期位置，初期速度などの設定
      call seting
!	write(*,*)nkoss_pt,nkoss_pt_flat,nkoss_pt_nano
!	write(*,*)nkoss_pt_bot,nhos_nanoxy,nhos_nanoyz
!	write(*,*)xsyu10,ysyu10,zsyu10
! 系内の全分子の並進速度の補正
      call cortra1
      call cortra2
! 系内の全分子の温度の補正
      call scale2
      call langevin



! メインルーチン
      do 10 i=1,maxstep
        nowstp = i

      if (mod(nowstp,1000) .eq. 0) then
          write(*,*)nowstp
      endif
      !		 call scale1

      if(nowstp.le.n_time_t)then
! 系内の全分子の温度の補正
!           call scale1
          call scale2
      endif

!call langevin

      if(mod(nowstp,50).eq.1)then
        call book  !粒子登録法
      endif

! 各分子に働く力，速度，位置の分子動力学計算
        ! call calcu

! 境界条件の付与
      call boundary		 
      if(nowstp.ge.ns0)then
        call heatflux
        call heatflux_lange
      endif

      if(mod(nowstp, 1000) .eq. 1) then
! データの出力１
        !  call record
! データの出力２
        !  call record2  
        call pvch !アニメーションの作成  
        if(nowstp.gt.n_time_t)then
          call record_maxwell
        endif
      endif

      if(nowstp.gt.ns0)then
      !  call bunpu
      !  call density_2d
      endif

10     continue

!call calcu_ddl


!!!!!!!!!!!!!!!!!!緩和計算の終了!!!!!!!!!!!!!!!!
! データの出力３
      !call write_pos_vel
      !		write(*,*)slab_v,dz
      write(*,*) 'Calculation END'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      stop
      end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine seting
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision(a-h,o-z)
      implicit integer(n)
      include 'inc/parameternano3.dat'
      include 'inc/variable4.dat'
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      zmass_pt = 1.00D+26*bunsi_pt/avoga
      zmass_ar = 1.00D+26*bunsi_ar/avoga
      !	  write(*,*) bunsi_ar
      sig1   = sigpp
      eps1   = epspp
      sig2   = sigaa
      eps2   = epsaa
      !	  write(*,*) sigaa,sig2,epsaa,eps2
      sig3   = sigap
      eps3_bot = alpha_bot * epsap
      eps3_top = alpha_top * epsap
      cforce1 = 24.00D0*eps1/sig1
      cforce2 = 24.00D0*eps2/sig2
      !	  write(*,*) cforce2
      cforce3_bot = 24.00D0*eps3_bot/sig3
      !	  write(*,*)cforce3_bot
      cforce3_top = 24.00D0*eps3_top/sig3
      !      write(*,*)xsyul0
      !      write(*,*)ysyul0 
      !      write(*,*)zsyul0
cccccccccccccccccccccccccccccccccccccccccccccccccc
      num = 0
      x = 0.0000D0
      y = 0.0000D0
      z = 0.0000D0
      ofstx1 = 0.0000D0
      ofsty1 = 0.0000D0
      ofstz1 = 0.0000D0
      underx = ptkosi3 / 4.0D0
      undery = ptkosi3 / 4.0D0
      underz = ptkosi3 / 4.0D0	 
      upz = underz+dble(n1z_pt)*ptkosi3/2.0D0
      topz = zsyul0 - (dble(n2z_pt-1)*ptkosi3/2.00D0+underz)
      !      write(*,*)underx,undery,underz,topz	  
      ! Pt-Pt間のカットオフ
      xcutof1 = xsyul0 - cutoff33*sigpp
      ycutof1 = ysyul0 - cutoff33*sigpp
      zcutof1 = zsyul0 - cutoff33*sigpp
      !	  write(*,*)sigpp,sigap,sigaa
      ! Ar-Ar間のカットオフ
      xcutof2 = xsyul0 - cutoff33*sigaa
      ycutof2 = ysyul0 - cutoff33*sigaa
      zcutof2 = zsyul0 - cutoff33*sigaa
      !	  write(*,*)xcutof2,ycutof2,zcutof2
      ! Pt-Ar間のカットオフ
      xcutof3 = xsyul0 - cutoff33*sigap
      ycutof3 = ysyul0 - cutoff33*sigap
      zcutof3 = zsyul0 - cutoff33*sigap

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!Pt粒子並べる処理(下平滑面)
      stdist = ptkosi3
      do 100 k=1,n1z_pt
      z =underz + dble(k-1)*stdist/2.0D0
        do 101 i=1,n1x_pt
          x = underx + dble(i-1)*stdist/2.0D0
          do 102 j=1,n1y_pt
          if(mod(k,2) .eq. 0) then
            if(mod(i,2) .eq. 0) then
              y = undery + dble(j-1)*stdist
            else
              y = undery + dble(j-1)*stdist + stdist/2.0D0
          endif
          else
            if(mod(i,2) .eq. 0) then
              y = undery + dble(j-1)*stdist + stdist/2.0D0
            else
              y = undery + dble(j-1)*stdist
            endif
          endif
          num = num + 1
          posx_pt(num) = x
          posy_pt(num) = y
          posz_pt(num) = z
102         continue
101       continue
      write(210,'(2E15.7)')posz_pt(num)
100   continue

ccccccccccccccccccccccccccccccccccccccccccccccccccccc!ナノ構造
      do 500 k=1,n1x_pt
      if(k.ge.3)then
        if(k.le.10)then
          x = underx+dble(k-1)*stdist/2.0D0
          do 501 i=1,nnz_pt
            z = upz + dble(i-1)*stdist/2.0D0
            do 502 j=1,nny_pt
              if(mod(k,2).eq.0)then
                if(mod(i,2).eq.0)then
                  y = undery + dble(j-1)*stdist
                else
                  y = undery + dble(j-1)*stdist + stdist/2.0D0
                endif
              else
                if(mod(i,2) .eq. 0) then
                  y = undery + dble(j-1)*stdist + stdist/2.0D0
                else
                  y = undery + dble(j-1)*stdist
                endif
              endif    
              num = num+1
              posx_pt(num) = x
              posy_pt(num) = y
              posz_pt(num) = z
502           continue
      write(210,'(2e15.7)')posz_pt(num)
501         continue
        endif
      endif

!       if(k.ge.14)then
!         if(k.le.18)then
!           x = underx+dble(k-1)*stdist/2.0D0
!           do 505 i=1,nnz_pt
!             z = upz + dble(i-1)*stdist/2.0D0
!             do 506 j=1,nny_pt
!               if(mod(k,2).eq.0)then
!                 if(mod(i,2).eq.0)then
!                   y = undery + dble(j-1)*stdist
!                 else
!                   y = undery + dble(j-1)*stdist + stdist/2.0D0
!                 endif
!               else
!                 if(mod(i,2) .eq. 0) then
!                   y = undery + dble(j-1)*stdist + stdist/2.0D0
!                 else
!                   y = undery + dble(j-1)*stdist
!                 endif
!               endif    
!               num = num+1
!               posx_pt(num) = x
!               posy_pt(num) = y
!               posz_pt(num) = z
! 506           continue
! 505         continue
!         endif
!       endif
500   continue

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc!Pt粒子並べる処理(上の平板)
      do 199 k=1,n2z_pt
      z = topz + dble(k-1)*stdist/2.0D0
          do 1100 i=1, n2x_pt
            x = underx + dble(i-1)*stdist/2.0D0
            do 1105 j=1,n2y_pt
            if(mod(k,2) .eq. 0) then
              if(mod(i,2) .eq. 0) then
                  y = undery + dble(j-1)*stdist
              else
                  y = undery + dble(j-1)*stdist + stdist/2.0D0
              endif
              else
              if(mod(i,2) .eq. 0) then
                  y = undery + dble(j-1)*stdist + stdist/2.0D0
              else
                  y = undery + dble(j-1)*stdist
              endif
            endif
            num = num + 1
            posx_pt(num) = x
            posy_pt(num) = y
            posz_pt(num) = z
1105           continue
1100        continue
      write(210,'(2E15.7)')posz_pt(num)
199   continue
!      write(*,*)n2x_pt,n2y_pt,n2z_pt
!      write(*,*)topz
cccccccccccccccccccccccccccccccccccccccccccccccccc!Pt粒子に初期速度を与える処理
      num = 0
      cr = 1.00D-6
      do 200 i=1, nkoss_pt
      read(1,*)ran
      theta = pi*ran
      read(2,*)ran
      phi = 2.000D0*pi*ran
      vx = cr*dsin(theta)*dcos(phi)
      vy = cr*dsin(theta)*dsin(phi)
      vz = cr*dcos(theta)
      num = num + 1
      velx_pt(num) = vx
      vely_pt(num) = vy
      velz_pt(num) = vz
200   continue

!!!!!!!!!!!!!!!!!固定層について!!!!!!!!!!!!!!!
      do 201 i = 1, nhos
      velx_pt(i) = 0.00D0
      vely_pt(i) = 0.00D0
      velz_pt(i) = 0.00D0
201	  continue

      do 215 i = nkoss_pt-nhos+1, nkoss_pt
      velx_pt(i) = 0.00D0
      vely_pt(i) = 0.00D0
      velz_pt(i) = 0.00D0
215	  continue


ccccccccccccccccccccccccccccccccccccccccccccccc!Ar粒子の初期配置
      num = 0
      x = 0.0000D0
      y = 0.0000D0
      z = 0.0000D0
      arofstx1 = underx*2.00D0
      arofsty1 = arkosi3/4.00D0
      arofstz1 = underz*2.0D0 + dble(n1z_pt+nnz_pt)*ptkosi3/2.00D0

!	   write(*,*)arofstx1,arofsty1,arofstz1
      stdist = arkosi3
      stdist = arkosi3
!	   write(*,*)stdist,arkosi3
      do 299 k=1,nz_ar
        z = arofstz1 + dble(k-1)*stdist/2.0D0
          do 300 i=1,nx_ar
            x = arofstx1 + dble(i-1)*stdist/2.0D0
            do 305 j=1,ny_ar
            if(mod(k,2) .eq. 0) then
              if(mod(i,2) .eq. 0) then
                y = arofsty1 + dble(j-1)*stdist + stdist/2.0D0 
              else
                y = arofsty1 + dble(j-1)*stdist 
              endif
            else
              if(mod(i,2) .eq. 0) then
                y = arofsty1 + dble(j-1)*stdist 
              else
                y = arofsty1 + dble(j-1)*stdist + stdist/2.0D0
              endif
          endif
            num = num + 1
            posx_ar(num) = x
            posy_ar(num) = y
            posz_ar(num) = z
!			  write(*,*)posx_ar(num)
305           continue
300         continue
299      continue
cccccccccccccccccccccccccccccccccccccccccccccccccc
      num = 0
      cr = 1.00D-6
      do 250 i=1, nkoss_ar
        read(1,*)ran
        theta = pi*ran
        read(2,*)ran
        phi = 2.000D0*pi*ran
        vx = dsin(theta)*dcos(phi)*cr
        vy = dsin(theta)*dsin(phi)*cr
        vz = dcos(theta)*cr
        num = num + 1
        velx_ar(num) = vx
        vely_ar(num) = vy
        velz_ar(num) = vz
250     continue
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine cortra1	!Pt粒子について全体重視の並進速度を0にする処理
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision(a-h,o-z)
      implicit integer(n)
      include 'inc/parameternano3.dat'
      include 'inc/variable4.dat'
***********************************************************************
      trvx = 0.00000D0
      trvy = 0.00000D0
      trvz = 0.00000D0

      do 110 i=1,nkoss_pt
      trvx = trvx + velx_pt(i)
      trvy = trvy + vely_pt(i)
      trvz = trvz + velz_pt(i)
110   continue

      trvx = trvx/dble(nkoss_pt)
      trvy = trvy/dble(nkoss_pt)
      trvz = trvz/dble(nkoss_pt)

      do 120 j=1,nkoss_pt
      velx_pt(j) = velx_pt(j) - trvx
      vely_pt(j) = vely_pt(j) - trvy
      velz_pt(j) = velz_pt(j) - trvz
!		 write(*,*)velx_pt(j)
120   continue
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine cortra2 !Ar粒子について全体重視の並進速度を0にする処理
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision(a-h,o-z)
      implicit integer(n)
      include 'inc/parameternano3.dat'
      include 'inc/variable4.dat'
************************************************************************
      trvx = 0.00000D0
      trvy = 0.00000D0
      trvz = 0.00000D0
      do 20 i= 1, nkoss_ar
      trvx = trvx + velx_ar(i)
      trvy = trvy + vely_ar(i)
      trvz = trvz + velz_ar(i)
!			write(*,*)trvx,velx_ar(i)
!			write(*,*)velx_ar(i)
20     continue
      trvx = trvx/dble(nkoss_ar)
      trvy = trvy/dble(nkoss_ar)
      trvz = trvz/dble(nkoss_ar)
!			write(*,*)trvx
      do 30 j= 1, nkoss_ar
      velx_ar(j) = velx_ar(j) - trvx
      vely_ar(j) = vely_ar(j) - trvy
      velz_ar(j) = velz_ar(j) - trvz
!			write(*,*)velx_ar(j)
30    continue
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine scale1   !Pt粒子について行う速度スケーリング
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision(a-h,o-z)
      implicit integer(n)
      include 'inc/parameternano3.dat'
      include 'inc/variable4.dat'
************************************************************************
      temptp = 0.00D0
      do 601 i= 1, nkoss_pt/2
      vel_pt = velx_pt(i)**2 + vely_pt(i)**2 + velz_pt(i)**2
      temptp = temptp + vel_pt
601     continue
      temptp = temptp/dble(nkoss_pt/2)/1.000D+16
      aimtem = atemp_pt_bot
      aimnot = 3.00D0*boltz*aimtem/zmass_pt
      baiss = dsqrt(aimnot/temptp)
      do 602 i= 1, nkoss_pt/2
        velx_pt(i) = velx_pt(i)*baiss
        vely_pt(i) = vely_pt(i)*baiss
        velz_pt(i) = velz_pt(i)*baiss
!		  write(*,*)velx_pt(i200),vely_pt(i200),velz_pt(i200)
602     continue

      temptp = 0.00D0
      do 603 i= nkoss_pt/2+1, nkoss_pt
      vel_pt = velx_pt(i)**2 + vely_pt(i)**2 + velz_pt(i)**2
      temptp = temptp + vel_pt
603     continue
      temptp = temptp/dble(nkoss_pt/2)/1.000D+16
      aimtem = atemp_pt_top
      aimnot = 3.00D0*boltz*aimtem/zmass_pt
      baiss = dsqrt(aimnot/temptp)
      do 604 i= nkoss_pt/2+1, nkoss_pt
        velx_pt(i) = velx_pt(i)*baiss
        vely_pt(i) = vely_pt(i)*baiss
        velz_pt(i) = velz_pt(i)*baiss
!		  write(*,*)velx_pt(i200),vely_pt(i200),velz_pt(i200)
604     continue
      return        
      end		

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine scale2   !Ar粒子について行う速度スケーリング
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision(a-h,o-z)
      implicit integer(n)
      include 'inc/parameternano3.dat'
      include 'inc/variable4.dat'
************************************************************************
      temptp = 0.00D0
      do 700 i= 1, nkoss_ar
      vel_ar = velx_ar(i)**2 + vely_ar(i)**2 + velz_ar(i)**2
      temptp = temptp + vel_ar
700     continue
      temptp = temptp/dble(nkoss_ar)/1.000D+16
      aimtem = atemp_ar
      aimnot = 3.00D0*boltz*aimtem/zmass_ar
      baiss = dsqrt(aimnot/temptp)
      do 701 i200= 1, nkoss_ar
        velx_ar(i200) = velx_ar(i200)*baiss
        vely_ar(i200) = vely_ar(i200)*baiss
        velz_ar(i200) = velz_ar(i200)*baiss
!		  write(*,*)velx_ar(i200),vely_ar(i200),velz_ar(i200)
701     continue
      return        
      end

ccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine jyusin_ar
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision(a-h,o-z)
      implicit integer(n)
      include 'inc/parameternano3.dat'
      include 'inc/variable4.dat'

************************************************************************
      cmsx = xsyul0/2.00D0
      cmsy = ysyul0/2.00D0
      cmsz = zsyul0/2.00D0
      tcmsx = 0.0000D0
      tcmsy = 0.0000D0
      tcmsz = 0.0000D0
      do 100 i= 1, nkoss_ar
      tcmsx= tcmsx + posx_ar(i)
      tcmsy= tcmsy + posy_ar(i)
      tcmsz= tcmsz + posz_ar(i)
100     continue
      tcmsx = cmsx - tcmsx/dble(nkoss_ar)
      tcmsy = cmsy - tcmsy/dble(nkoss_ar)
      tcmsz = cmsz - tcmsz/dble(nkoss_ar)
      do 200 j= 1, nkoss_ar
      posx_ar(j) = posx_ar(j) + tcmsx 
      posy_ar(j) = posy_ar(j) + tcmsy
      posz_ar(j) = posz_ar(j) + tcmsz
200    continue
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine book
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision(a-h,o-z)
      implicit integer(n)
      include 'inc/parameternano3.dat'
      include 'inc/variable4.dat'

************************************************************************
!粒子の最大速度について計算
      v_max = 0.0D0
      do 900 i=1,nkoss_pt
        vel_w = velx_pt(i)**2+vely_pt(i)**2+velz_pt(i)**2
        vel_w = dsqrt(vel_w)
        if(vel_w.gt.v_max)then
          v_max = vel_w
        endif
900     continue

      do 901 i=1,nkoss_ar
        vel_w = velx_ar(i)**2+vely_ar(i)**2+velz_ar(i)**2
        vel_w = dsqrt(vel_w)
        if(vel_w.gt.v_max)then
          v_max = vel_w
        endif
901     continue

      if(v_max.gt.vmax_pt)then
        vmax_pt = v_max
      endif
      if(v_max.gt.vmax_ar)then
        vmax_ar = v_max
      endif

!!!!!!!!!!!!!!!!!!!!!!!!!Pt-Pt間の粒子登録法!!!!!!!!!!!!!!!!!!!!!!!
      bkr = safe*(cutoff33*sig1+vmax_pt*dt*50.0)
      bkr2 = bkr*bkr/sig1/sig1
      bkx = xsyul0 - bkr
      bky = ysyul0 - bkr
      bkz = zsyul0 - bkr
      np1 = 0

!粒子登録を行うかの判定
      do 902 i=1,nkoss_pt-1
        do 903 j=i+1,nkoss_pt
          divx = posx_pt(i) - posx_pt(j)
          !周期境界を考慮
          if(divx.lt.-xsyul0/2.00D0)then
            divx = divx + xsyul0
          else if(divx.gt.xsyul0/2.00D0)then
            divx = divx - xsyul0
          endif
          
          divy = posy_pt(i) - posy_pt(j)
          if(divy.lt.-ysyul0/2.00D0)then
            divy = divy + ysyul0
          else if(divy.gt.ysyul0/2.00D0)then
            divy = divy - ysyul0
          endif
          
          divz = posz_pt(i) - posz_pt(j)
      
  !粒子登録範囲外の粒子をはじく処理
        divx = divx/sig1
        if(divx.gt.(bkr/sig1))then
          go to 903
        endif
        if(divx.lt.-(bkr/sig1))then
          go to 903
        endif
        
        divy = divy/sig1
        if(divy.gt.(bkr/sig1))then
          go to 903
        endif
        if(divy.lt.-(bkr/sig1))then
          go to 903
        endif			

        divz = divz/sig1
        if(divz.gt.(bkr/sig1))then
          go to 903
        endif
        if(divz.lt.-(bkr/sig1))then
          go to 903
        endif

        dit2 = divx**2+divy**2+divz**2
  !粒子登録を行う粒子についての記憶
        if(dit2.lt.bkr2)then
          np1 =np1+1
          np11(np1) = i
          np12(np1) = j	
        endif
903       continue
902     continue		

!!!!!!!!!!!!!!!!!!!!!!!!!Ar-Ar間の粒子登録法!!!!!!!!!!!!!!!!!!!!!!!
      bkr = safe*(cutoff33*sig2+vmax_ar*dt*50.0)
      bkr2 = bkr*bkr/sig2/sig2
      bkx = xsyul0 - bkr
      bky = ysyul0 - bkr
      bkz = zsyul0 - bkr
      np2 = 0

!粒子登録を行うかの判定
      do 904 i=1,nkoss_ar-1
        do 905 j=i+1,nkoss_ar
          divx = posx_ar(i) - posx_ar(j)
          !周期境界を考慮
          if(divx.lt.-xsyul0/2.00D0)then
            divx = divx + xsyul0
          else if(divx.gt.xsyul0/2.00D0)then
            divx = divx - xsyul0
          endif
          
          divy = posy_ar(i) - posy_ar(j)
          if(divy.lt.-ysyul0/2.00D0)then
            divy = divy + ysyul0
          else if(divy.gt.ysyul0/2.00D0)then
            divy = divy - ysyul0
          endif
          
          divz = posz_ar(i) - posz_ar(j)
    
!粒子登録範囲外の粒子をはじく処理
      divx = divx/sig2
      if(divx.gt.(bkr/sig2))then
        go to 905
      endif
      if(divx.lt.-(bkr/sig2))then
        go to 905
      endif
      
      divy = divy/sig2
      if(divy.gt.(bkr/sig2))then
        go to 905
      endif
      if(divy.lt.-(bkr/sig2))then
        go to 905
      endif			

      divz = divz/sig2
      if(divz.gt.(bkr/sig2))then
        go to 905
      endif
      if(divz.lt.-(bkr/sig2))then
        go to 905
      endif

      dit2 = divx**2+divy**2+divz**2
!粒子登録を行う粒子についての記憶
      if(dit2.lt.bkr2)then
        np2 =np2+1
        np21(np2) = i
        np22(np2) = j	
      endif
905       continue
904     continue	

!!!!!!!!!!!!!!!!!!!!!!!!!Pt-Ar間の粒子登録法!!!!!!!!!!!!!!!!!!!!!!!
      bkr = safe*(cutoff33*sig3+v_max*dt*50.0)
      bkr2 = bkr*bkr/sig3/sig3
      bkx = xsyul0 - bkr
      bky = ysyul0 - bkr
      bkz = zsyul0 - bkr
      np3 = 0

!粒子登録を行うかの判定
      do 906 i=1,nkoss_pt
        do 907 j=1,nkoss_ar
          divx = posx_pt(i) - posx_ar(j)
          !周期境界を考慮
          if(divx.lt.-xsyul0/2.00D0)then
            divx = divx + xsyul0
          else if(divx.gt.xsyul0/2.00D0)then
            divx = divx - xsyul0
          endif
          
          divy = posy_pt(i) - posy_ar(j)
          if(divy.lt.-ysyul0/2.00D0)then
            divy = divy + ysyul0
          else if(divy.gt.ysyul0/2.00D0)then
            divy = divy - ysyul0
          endif
          
          divz = posz_pt(i) - posz_ar(j)
    
!粒子登録範囲外の粒子をはじく処理
      divx = divx/sig3
      if(divx.gt.(bkr/sig3))then
        go to 907
      endif
      if(divx.lt.-(bkr/sig3))then
        go to 907
      endif
      
      divy = divy/sig3
      if(divy.gt.(bkr/sig3))then
        go to 907
      endif
      if(divy.lt.-(bkr/sig3))then
        go to 907
      endif			

      divz = divz/sig3
      if(divz.gt.(bkr/sig3))then
        go to 907
      endif
      if(divz.lt.-(bkr/sig3))then
        go to 907
      endif

      dit2 = divx**2+divy**2+divz**2
!粒子登録を行う粒子についての記憶
      if(dit2.lt.bkr2)then
        np3 =np3+1
        np31(np3) = i
        np32(np3) = j
      endif
907       continue
906     continue
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine calcu
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision(a-h,o-z)
      implicit integer(n)
      include 'inc/parameternano3.dat'
      include 'inc/variable4.dat'
*******************************************************************!各粒子に働く力の初期化
      do 95 i = 1, nkoss_pt
        forx_pt(i) = 0.0000D0
        fory_pt(i) = 0.0000D0
        forz_pt(i) = 0.0000D0
        poten_pt(i) = 0.0000D0
        ukine_pt(i) = 0.0000D0
!		  write(*,*)poten_pt(i),ukine_pt(i)
        force_qx(i) = 0.0000D0
        force_qy(i) = 0.0000D0
        force_qz(i) = 0.0000D0
95      continue
    
  
      do 96 i=1, nkoss_ar
        forx_ar(i) = 0.0000D0
        fory_ar(i) = 0.0000D0
        forz_ar(i) = 0.0000D0
        poten_ar(i) = 0.0000D0
        ukine_ar(i) = 0.0000D0
96		 continue

*****************************************************************		
!!!!!!!!!!!!!!!!!Ar-Pt間の相互作用力の計算!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do 2200 i1= 1, np3
        n1 = np31(i1)
        n2 = np32(i1)		  
          divx = posx_pt(n1)-posx_ar(n2)
          divy = posy_pt(n1)-posy_ar(n2)
          divz = posz_pt(n1)-posz_ar(n2)
    
!			write(*,*)divx,divy,divz

ccccc
! 周期境界条件を考慮
      if (divx .lt. -xcutof3) then
        divx = divx + xsyul0
      else if(divx .gt. xcutof3) then
        divx = divx - xsyul0
      endif
ccccc
      if (divy .lt. -ycutof3) then
        divy = divy + ysyul0
      else if(divy .gt. ycutof3) then
        divy = divy - ysyul0
      endif
ccccc
!        if (divz .lt. -zcutof3) then
!          divz = divz + zsyul0
!        else if(divz .gt. zcutof3) then
!          divz = divz - zsyul0
!        endif
ccccc
! カットオフ距離より遠い位置に存在する粒子は省略
      divx = divx/sig3
!		write(*,*)divx
      if (divx .gt. cutoff33) then
        goto 2200
      endif
      if (divx .lt. -cutoff33) then
        goto 2200
      endif

c
      divy = divy/sig3
!		write(*,*)divy
      if (divy .gt. cutoff33) then
        goto 2200
      endif
      if (divy .lt. -cutoff33) then
        goto 2200
      endif
c

      divz = divz/sig3
!		write(*,*)divz,cutoff33
      if (divz .gt. cutoff33) then
        goto 2200
      endif
      if (divz .lt. -cutoff33) then
        goto 2200
      endif
ccccccccccccccccccccccccccccccccccccccccccccccccccc
      dit2 = divx**2 + divy**2 + divz**2
      dist = dsqrt(dit2)
!		write(*,*)dist
      if(dist .gt. cutoff33) then
        goto 2200
      endif
      dit4 = dit2*dit2
      dit6 = dit4*dit2
      dit8 = dit4*dit4
      dit12 = dit6*dit6
      dit14 = dit8*dit6
        if(i1.le.nkoss_pt_bot)then
          ppp = 4.00D0*eps3_bot*(1.00D0/dit12-1.00D0/dit6)
          force = cforce3_bot*(-2.00D0/dit14+1.00D0/dit8)
        elseif(i1.gt.nkoss_pt_bot)then
          ppp = 4.00D0*eps3_top*(1.00D0/dit12-1.00D0/dit6)
          force = cforce3_top*(-2.00D0/dit14+1.00D0/dit8)
      endif
!		write(*,*)force,cforce2
      forcex_pt = -force*divx
      forcex_ar = -force*divx/zmass_ar
      forcey_pt = -force*divy
      forcey_ar = -force*divy/zmass_ar
      forcez_pt = -force*divz
      forcez_ar = -force*divz/zmass_ar
!		write(*,*)zmass_ar
!		write(*,*)forcex,forcey,forcez
cccccccccccccccccccccccccccccccccccccccccccccccccc !熱流束の計算
      force_qx(n1) = force_qx(n1) + forcex_pt
      force_qy(n1) = force_qy(n1) + forcey_pt
      force_qz(n1) = force_qz(n1) + forcez_pt
!		write(*,*)force_qz(n1)
cccccccccccccccccccccccccccccccccccccccccccccccccc!平板にかかる力の足し合わせ
      if(nowstp.gt.ns0)then
        if(n1.gt.nkoss_pt_bot)then
          ftotx_top = ftotx_top + forcex_pt
          ftoty_top = ftoty_top + forcey_pt
          ftotz_top = ftotz_top + forcez_pt
        else
          ftotx_bot = ftotx_bot + forcex_pt
          ftoty_bot = ftoty_bot + forcey_pt
          ftotz_bot = ftotz_bot + forcez_pt
        endif
      endif		
cccccccccccccccccccccccccccccccccccccccccccccccccc  !粒子にかかる力
      forx_pt(n1) = forx_pt(n1) + forcex_pt
      forx_ar(n2) = forx_ar(n2) - forcex_ar
!		write(*,*)forx_ar(i1),forx_ar(i2)
      fory_pt(n1) = fory_pt(n1) + forcey_pt
      fory_ar(n2) = fory_ar(n2) - forcey_ar
      forz_pt(n1) = forz_pt(n1) + forcez_pt
      forz_ar(n2) = forz_ar(n2) - forcez_ar
!		write(*,*)forx_ar(i1),forx_ar(i2)
      poten_pt(n1) = poten_pt(n1) + ppp*0.500D0
      poten_ar(n2) = poten_ar(n2) + ppp*0.500D0
2200  continue

***************************************************!平板にかかる力から圧力を計算
      fave_top_x = fave_top_x + ftotx_top
      fave_top_y = fave_top_y + ftoty_top
      fave_top_z = fave_top_z + ftotz_top
      forw_top_x = forw_top_x + ftotx_top
      forw_top_y = forw_top_y + ftoty_top
      forw_top_z = forw_top_z + ftotz_top

      fave_bot_x = fave_bot_x + ftotx_bot
      fave_bot_y = fave_bot_y + ftoty_bot
      fave_bot_z = fave_bot_z + ftotz_bot
      forw_bot_x = forw_bot_x + ftotx_bot
      forw_bot_y = forw_bot_y + ftoty_bot
      forw_bot_z = forw_bot_z + ftotz_bot

      if((mod(nowstp,5000).eq.0))then
      if(nowstp.gt.ns0)then
      fave_top_x = fave_top_x/(ysyul0*wallbot)*1.00D14/dble(5000)
      fave_top_y = fave_top_y/(x_plate*wallbot)*1.0D14/dble(5000)
      fave_top_z = fave_top_z/(x_plate*ysyul0)*1.0D14/dble(5000)
          
      fave_bot_x = fave_bot_x/(ysyul0*wallbot)*1.00D14/dble(5000)
      fave_bot_y = fave_bot_y/(x_plate*wallbot)*1.0D14/dble(5000)
      fave_bot_z = fave_bot_z/(x_plate*ysyul0)*1.00D14/dble(5000)	

      fff = (fave_top_z + (-1.0D0)*fave_bot_z)*0.50D0
      tot_fff = tot_fff + fff
      tot_fave_top_z = tot_fave_top_z+fave_top_z
      tot_fave_bot_z = tot_fave_bot_z+(-1.0D0)*fave_bot_z
!		write(*,*) '1'
      write(31,'(2D15.7)')dble(nowstp-5000)*dt*1.0D-6, fff,tot_fff

      if(nowstp.eq.maxstep)then
      !		  fave_top_z = fave_top_z / dble(maxstep-ns0)
      !		  fave_bot_z = fave_bot_z / dble(maxstep-ns0)
        tot_fff = tot_fff / (dble(maxstep-ns0)/dble(5000))
        tot_fave_top_z=tot_fave_top_z/(dble(maxstep-ns0)/dble(5000))
        tot_fave_bot_z=tot_fave_bot_z/(dble(maxstep-ns0)/dble(5000))		  
        write(33,'(3D15.7)')tot_fff,tot_fave_top_z,tot_fave_bot_z
      endif

      fave_top_x = 0.00D0
      fave_top_y = 0.00D0
      fave_top_z = 0.00D0
      fave_bot_x = 0.00D0
      fave_bot_y = 0.00D0
      fave_bot_z = 0.00D0
      endif
      endif

!!!!!!!!!!!!!!!!!!!!!!!!Pt-Pt間の相互作用力の計算!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do 2211 i=1,np1  !粒子登録した粒子を呼び出す
      n1 = np11(i)
      n2 = np12(i)
          divx = posx_pt(n1)-posx_pt(n2)
          divy = posy_pt(n1)-posy_pt(n2)
          divz = posz_pt(n1)-posz_pt(n2)
    
!			write(*,*)divx,divy,divz
ccccc
! 周期境界条件を考慮
      if (divx .lt. -xcutof1) then
        divx = divx + xsyul0
      else if(divx .gt. xcutof1) then
        divx = divx - xsyul0
      endif
ccccc
      if (divy .lt. -ycutof1) then
        divy = divy + ysyul0
      else if(divy .gt. ycutof1) then
        divy = divy - ysyul0
      endif
ccccc
!        if (divz .lt. -zcutof1) then
!          divz = divz + zsyul0
!        else if(divz .gt. zcutof1) then
!          divz = divz - zsyul0
!        endif
ccccc
! カットオフ距離より遠い位置に存在する粒子は省略
      divx = divx/sig1
!		write(*,*)divx
      if (divx .gt. cutoff33) then
        goto 2211
      endif
      if (divx .lt. -cutoff33) then
        goto 2211
      endif
c
      divy = divy/sig1
!		write(*,*)divy
      if (divy .gt. cutoff33) then
        goto 2211
      endif
      if (divy .lt. -cutoff33) then
        goto 2211
      endif
c
      divz = divz/sig1
!		write(*,*)divz,cutoff33
      if (divz .gt. cutoff33) then
        goto 2211
      endif
      if (divz .lt. -cutoff33) then
        goto 2211
      endif
ccccccccccccccccccccccccccccccccccccccccccccccccccc
      dit2 = divx**2 + divy**2 + divz**2
      dist = dsqrt(dit2)
!		write(*,*)dist
      if(dist .gt. cutoff33) then
        goto 2211
      endif
      dit4 = dit2*dit2
      dit6 = dit4*dit2
      dit8 = dit4*dit4
      dit12 = dit6*dit6
      dit14 = dit8*dit6
      ppp = 4.00D0*eps1*(1.00D0/dit12-1.00D0/dit6)
      force = cforce1*(-2.00D0/dit14+1.00D0/dit8)
      forcex = -force*divx
      forcey = -force*divy
      forcez = -force*divz
      !		write(*,*)forcex,forcey,forcez
      forx_pt(n1) = forx_pt(n1) + forcex
      forx_pt(n2) = forx_pt(n2) - forcex
      fory_pt(n1) = fory_pt(n1) + forcey
      fory_pt(n2) = fory_pt(n2) - forcey
      forz_pt(n1) = forz_pt(n1) + forcez
      forz_pt(n2) = forz_pt(n2) - forcez
      poten_pt(n1) = poten_pt(n1) + ppp*0.500D0
      poten_pt(n2) = poten_pt(n2) + ppp*0.500D0
2211  continue


!!!!!!!!!!!!!!!!!Ar-Ar間の相互作用力の計算!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do 2212 i=1,np2  !粒子登録した粒子を呼び出す
      n1 = np21(i)
      n2 = np22(i)
          divx = posx_ar(n1)-posx_ar(n2)
          divy = posy_ar(n1)-posy_ar(n2)
          divz = posz_ar(n1)-posz_ar(n2)
    
!			write(*,*)divx,divy,divz
ccccc
! 周期境界条件を考慮
      if (divx .lt. -xcutof2) then
        divx = divx + xsyul0
      else if(divx .gt. xcutof2) then
        divx = divx - xsyul0
      endif
ccccc
      if (divy .lt. -ycutof2) then
        divy = divy + ysyul0
      else if(divy .gt. ycutof2) then
        divy = divy - ysyul0
      endif
ccccc
!        if (divz .lt. -zcutof2) then
!          divz = divz + zsyul0
!        else if(divz .gt. zcutof2) then
!          divz = divz - zsyul0
!        endif
ccccc
! カットオフ距離より遠い位置に存在する粒子は省略
      divx = divx/sig2
!		write(*,*)divx
      if (divx .gt. cutoff33) then
        goto 2212
      endif
      if (divx .lt. -cutoff33) then
        goto 2212
      endif
c
      divy = divy/sig2
!		write(*,*)divy
      if (divy .gt. cutoff33) then
        goto 2212
      endif
      if (divy .lt. -cutoff33) then
        goto 2212
      endif
c
      divz = divz/sig2
!		write(*,*)divz,cutoff33
      if (divz .gt. cutoff33) then
        goto 2212
      endif
      if (divz .lt. -cutoff33) then
        goto 2212
      endif
ccccccccccccccccccccccccccccccccccccccccccccccccccc
      dit2 = divx**2 + divy**2 + divz**2
      dist = dsqrt(dit2)
!		write(*,*)dist
      if(dist .gt. cutoff33) then
        goto 2212
      endif
      dit4 = dit2*dit2
      dit6 = dit4*dit2
      dit8 = dit4*dit4
      dit12 = dit6*dit6
      dit14 = dit8*dit6
      ppp = 4.00D0*eps2*(1.00D0/dit12-1.00D0/dit6)
      force = cforce2*(-2.00D0/dit14+1.00D0/dit8)
      !		write(*,*)force,cforce2
      forcex = -force*divx/zmass_ar
      forcey = -force*divy/zmass_ar
      forcez = -force*divz/zmass_ar
      !		write(*,*)zmass_ar
      !		write(*,*)forcex,forcey,forcez
      forx_ar(n1) = forx_ar(n1) + forcex
      forx_ar(n2) = forx_ar(n2) - forcex
      !		write(*,*)forx_ar(i1),forx_ar(i2)
      fory_ar(n1) = fory_ar(n1) + forcey
      fory_ar(n2) = fory_ar(n2) - forcey
      forz_ar(n1) = forz_ar(n1) + forcez
      forz_ar(n2) = forz_ar(n2) - forcez
      !		write(*,*)forx_ar(i1),forx_ar(i2)
      poten_ar(n1) = poten_ar(n1) + ppp*0.500D0
      poten_ar(n2) = poten_ar(n2) + ppp*0.500D0
2212  continue

*********************************************************
      do 200 i=1, nkoss_pt
        forx_pt(i) = forx_pt(i) + flanx(i)
        fory_pt(i) = fory_pt(i) + flany(i)
        forz_pt(i) = forz_pt(i) + flanz(i)
200   continue
************************************************************************
!計算された力からPt粒子の速度，位置について更新する(固定層は更新しない)
      do 5001 i=1, nkoss_pt
      forx_pt(i) = forx_pt(i)/zmass_pt
      fory_pt(i) = fory_pt(i)/zmass_pt
      forz_pt(i) = forz_pt(i)/zmass_pt
      velx_lan(i) = velx_pt(i)+ 0.500D0*forx_pt(i)*dt
      vely_lan(i) = vely_pt(i)+ 0.500D0*fory_pt(i)*dt
      velz_lan(i) = velz_pt(i)+ 0.500D0*forz_pt(i)*dt
      vxene = velx_pt(i) + 0.500D0*forx_pt(i)*dt
      vyene = vely_pt(i) + 0.500D0*fory_pt(i)*dt
      vzene = velz_pt(i) + 0.500D0*forz_pt(i)*dt
      !		write(*,*)velx_pt(i),vely_pt(i),velz_pt(i)
      vene = vxene*vxene + vyene*vyene + vzene*vzene
      ukine_pt(i) = 0.500D0*zmass_pt*vene
      !		write(*,*)ukine_pt(i)
5001  continue

      do 5111 i = 1,nhos
      velx_pt(i) = 0.00D0
      vely_pt(i) = 0.00D0
      velz_pt(i) = 0.00D0
      posx_pt(i) = posx_pt(i)
      posy_pt(i) = posy_pt(i)
      posz_pt(i) = posz_pt(i)
5111  continue

      do 5221 i = nhos+1, nkoss_pt_bot
!	    vx_pt(i) = velx_pt(i) + 0.500D0*forx_pt(i)*dt
!     vy_pt(i) = vely_pt(i) + 0.500D0*fory_pt(i)*dt
!     vz_pt(i) = velz_pt(i) + 0.500D0*forz_pt(i)*dt

      velx_pt(i) = velx_pt(i) + forx_pt(i)*dt
      vely_pt(i) = vely_pt(i) + fory_pt(i)*dt
      velz_pt(i) = velz_pt(i) + forz_pt(i)*dt
      !		write(*,*)velx_pt(i),vely_pt(i),velz_pt(i)
      posx_pt(i) = posx_pt(i) + velx_pt(i)*dt
      posy_pt(i) = posy_pt(i) + vely_pt(i)*dt
      posz_pt(i) = posz_pt(i) + velz_pt(i)*dt

5221  continue

      do 5222 i = nkoss_pt_bot+1, nkoss_pt-nhos
!	    vx_pt(i) = velx_pt(i) + 0.500D0*forx_pt(i)*dt
!     vy_pt(i) = vely_pt(i) + 0.500D0*fory_pt(i)*dt
!     vz_pt(i) = velz_pt(i) + 0.500D0*forz_pt(i)*dt

      velx_pt(i) = velx_pt(i) + forx_pt(i)*dt
      vely_pt(i) = vely_pt(i) + fory_pt(i)*dt
      velz_pt(i) = velz_pt(i) + forz_pt(i)*dt
      posx_pt(i) = posx_pt(i) + velx_pt(i)*dt
      posy_pt(i) = posy_pt(i) + vely_pt(i)*dt
      posz_pt(i) = posz_pt(i) + velz_pt(i)*dt

5222  continue

      do 5223 i=nkoss_pt-nhos+1,nkoss_pt
      velx_pt(i) = 0.00D0
      vely_pt(i) = 0.00D0
      velz_pt(i) = 0.00D0
      posx_pt(i) = posx_pt(i)
      posy_pt(i) = posy_pt(i)
      posz_pt(i) = posz_pt(i)
5223  continue	    

******************************************************************
!計算された力からAr粒子の速度，位置について更新する
      do 5002 i=1, nkoss_ar
        vxene = velx_ar(i) + forx_ar(i)*0.500D0*dt
        vyene = vely_ar(i) + fory_ar(i)*0.500D0*dt
        vzene = velz_ar(i) + forz_ar(i)*0.500D0*dt
      !		  write(*,*)velx_ar(i),forx_ar(i),vxene
        vene = vxene*vxene + vyene*vyene + vzene*vzene
        ukine_ar(i) = 0.500D0*zmass_ar*vene
      !		  write(*,*)ukine_ar(i)
5002   continue

      do 5102 i=1, nkoss_ar
        velx_ar(i) = velx_ar(i) + forx_ar(i)*dt
        vely_ar(i) = vely_ar(i) + fory_ar(i)*dt
        velz_ar(i) = velz_ar(i) + forz_ar(i)*dt
!          write(*,*)velx_ar(i)		  
      posx_ar(i) = posx_ar(i) + velx_ar(i)*dt
      posy_ar(i) = posy_ar(i) + vely_ar(i)*dt
      posz_ar(i) = posz_ar(i) + velz_ar(i)*dt
!		  write(*,*)posx_ar(i),posy_ar(i),posz_ar(i)
5102    continue
*********************************************************
      return
      end

************************************************************************

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine boundary
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision(a-h,o-z)
      implicit integer(n)
      include 'inc/parameternano3.dat'
      include 'inc/variable4.dat'
************************************************************************
!        write(*,*)xsyul0,ysyul0,zsyul0
      do 97 i=1, nkoss_pt
        if(posx_pt(i) .lt. 0.00D0) then
            posx_pt(i) = posx_pt(i) + xsyul0
        else if(posx_pt(i) .gt. xsyul0) then
            posx_pt(i) = posx_pt(i) - xsyul0
        endif
        if(posy_pt(i) .lt. 0.00D0) then
          posy_pt(i) = posy_pt(i) + ysyul0
        else if(posy_pt(i) .gt. ysyul0) then
            posy_pt(i) = posy_pt(i) - ysyul0
        endif
        if(posz_pt(i) .lt. 0.00D0) then
          posz_pt(i) = -1.0D0*posz_pt(i)
          velz_pt(i) = -1.0D0*velz_pt(i)
        else if(posz_pt(i) .gt. zsyul0) then
          posz_pt(i) = 2.0D0*zsyul0 - posz_pt(i)
          velz_pt(i) = -1.0D0*velz_pt(i)
        endif
97      continue

      do 950 i=1, nkoss_ar
        if(posx_ar(i) .lt. 0.00D0) then
            posx_ar(i) = posx_ar(i) + xsyul0
        else if(posx_ar(i) .gt. xsyul0) then
            posx_ar(i) = posx_ar(i) - xsyul0
        endif
        if(posy_ar(i) .lt. 0.00D0) then
            posy_ar(i) = posy_ar(i) + ysyul0
        else if(posy_ar(i) .gt. ysyul0) then
            posy_ar(i) = posy_ar(i) - ysyul0
        endif
        if(posz_ar(i) .lt. 0.00D0) then
            posz_ar(i) = -1.0D0*posz_pt(i)
            velz_pt(i) = -1.0D0*velz_pt(i)
        else if(posz_ar(i) .gt. zsyul0) then
            posz_ar(i) = 2.0D0*zsyul0 - posz_pt(i)
            velz_pt(i) = -1.0D0*velz_pt(i)
        endif
!		  write(*,*)xsyul0
950     continue
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine langevin
cccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision(a-h,o-z)
      implicit integer(n)
      include 'inc/parameternano3.dat'
      include 'inc/variable4.dat'
************************************************************************
************************************************************************
      do 82 i=1,nkoss_pt
        flanx(i) = 0.00D0
        flany(i) = 0.00D0
        flanz(i) = 0.00D0
82      continue
**************************************************************!底のLangevin層
      margin1 = nhos
      h1 = planck / (2.00D0*pi) !ディラック定数
      alpha_d0 = zmass_pt*1.00D-26*pi*debye*boltz/6.00D0/h1		!ダンパーの減衰係数
      alpha_d = alpha_d0/1.00D-26*1.00D-15
      a = 2.00D0*alpha_d0*boltz*atemp_pt_bot/(dt*1.00D-15)
      hyouj = dsqrt(a)
      !		write(*,*)a

      do 200 i = nhos+1, 2*nhos 
        !ランダム力のx方向成分の計算
        bb1 = ranran(iransu)
        if((-1.00D-10 .lt. bb1).and.(bb1 .lt. 1.00D-10))then
            bb1 = ranran(iransu + 1)
        endif
        bb1 = dabs(bb1)   !ランダム力を出すための乱数
        bb2 = ranran(iransu)*2.00D0*pi
        x1 = dsqrt(-2.0D0*dlog(bb1))
        ran2 = x1 * dsin(bb2)
        frandx = ran2*hyouj*1.00D6   !底面のPhantom粒子にかけるランダム力のx成分
      !		   write(*,*)frandx
        
        !ランダム力のy方向成分の計算
        bb1 = ranran(iransu)
        if((-1.00D-10 .lt. bb1).and.(bb1 .lt. 1.00D-10))then
            bb1 = ranran(iransu + 1)
        endif
        bb1 = dabs(bb1)   !ランダム力を出すための乱数
        bb2 = ranran(iransu)*2.00D0*pi
        y1 = dsqrt(-2.0D0*dlog(bb1))
        ran2 = y1 * dsin(bb2)
        frandy = ran2*hyouj*1.00D6   !底面のPhantom粒子にかけるランダム力のy成分

        !ランダム力のz方向成分の計算
        bb1 = ranran(iransu)
        if((-1.00D-10 .lt. bb1).and.(bb1 .lt. 1.00D-10))then
            bb1 = ranran(iransu + 1)
        endif
        bb1 = dabs(bb1)   !ランダム力を出すための乱数
        bb2 = ranran(iransu)*2.00D0*pi
        z1 = dsqrt(-2.0D0*dlog(bb1))
        ran2 = z1 * dsin(bb2)
        frandz = ran2*hyouj*1.00D6   !底面のPhantom粒子にかけるランダム力のz成分	
      !           write(*,*)frandz

        !Phantom粒子にかかるダンパーの力
        fdumx = -alpha_d*velx_pt(i)
        fdumy = -alpha_d*vely_pt(i)
        fdumz = -alpha_d*velz_pt(i)
      !		   write(*,*)fdumx,fdumy,fdumz
      !           write(*,*)velx_pt(i),vely_pt(i),velz_pt(i)
        
        !Phantom粒子にかける力の合計
        flanx(i) = frandx + fdumx
        flany(i) = frandy + fdumy
        flanz(i) = frandz + fdumz
        
      !		   write(*,*)fdumx,frandx,flanx(i),velx_pt(i)
      !           write(*,*)flanx(i),flany(i),flanz(i)
      !           write(*,*)flanx(i),frandx,fdumx
      !           write(*,*)alpha_d

200     continue

ccccccccccccccccccccccccccccccccccccccccccccccc!上面のLangevin層
      margin1 = nhos
      h1 = planck / (2.00D0*pi)
      alpha_d0 = zmass_pt*1.00D-26*pi*debye*boltz/6.00D0/h1 !ダンパーの減衰係数
      alpha_d = alpha_d0/1.00D-26*1.00D-15
      a = 2.00D0*alpha_d0*boltz*atemp_pt_top/(dt*1.00D-15)
      hyouj = dsqrt(a)
      !		write(*,*)a
      !		write(*,*)atemp_pt_bot,atemp_pt_top
      !		write(*,*)nkoss_pt-2*nhos+1,nkoss_pt-nhos
      do 220 i=nkoss_pt-2*nhos+1,nkoss_pt-nhos
          
  
      !ランダム力のx方向成分の計算
      bb1 = ranran(iransu)
      if((-1.00D-10 .lt. bb1).and.(bb1 .lt. 1.00D-10))then
          bb1 = ranran(iransu + 1)
      endif
      bb1 = dabs(bb1)   !ランダム力を出すための乱数
      bb2 = ranran(iransu)*2.00D0*pi
      x1 = dsqrt(-2.0D0*dlog(bb1))
      ran2 = x1 * dsin(bb2)
      frandx = ran2*hyouj*1.00D6   !底面のPhantom粒子にかけるランダム力のx成分
      
      
      !ランダム力のy方向成分の計算
      bb1 = ranran(iransu)
      if((-1.00D-10 .lt. bb1).and.(bb1 .lt. 1.00D-10))then
          bb1 = ranran(iransu + 1)
      endif
      bb1 = dabs(bb1)   !ランダム力を出すための乱数
      bb2 = ranran(iransu)*2.00D0*pi
      y1 = dsqrt(-2.0D0*dlog(bb1))
      ran2 = y1 * dsin(bb2)
      frandy = ran2*hyouj*1.00D6   !底面のPhantom粒子にかけるランダム力のy成分

      !ランダム力のz方向成分の計算
      bb1 = ranran(iransu)
      if((-1.00D-10 .lt. bb1).and.(bb1 .lt. 1.00D-10))then
          bb1 = ranran(iransu + 1)
      endif
      bb1 = dabs(bb1)   !ランダム力を出すための乱数
      bb2 = ranran(iransu)*2.00D0*pi
      z1 = dsqrt(-2.0D0*dlog(bb1))
      ran2 = z1 * dsin(bb2)
      frandz = ran2*hyouj*1.00D6   !底面のPhantom粒子にかけるランダム力のz成分	

      !Phantom粒子にかかるダンパーの力
      fdumx = -alpha_d*velx_pt(i) 
      fdumy = -alpha_d*vely_pt(i)
      fdumz = -alpha_d*velz_pt(i)
      
    !		   write(*,*)fdumx,flanx(i),velx_pt(i)
      !Phantom粒子にかける力の合計
      flanx(i) = frandx + fdumx
      flany(i) = frandy + fdumy 
      flanz(i) = frandz + fdumz
      
    !		   write(*,*)fdumx,flanx(i),velx_pt(i)
220     continue
      return 
      end		

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine record
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision(a-h,o-z)
      implicit integer(n)
      include 'inc/parameternano3.dat'
      include 'inc/variable4.dat'
************************************************************************
      num = 0

      do 500 i=1, nkoss_pt
        write(3,'(I6, 3E15.7)')i, posx_pt(i), posy_pt(i), posz_pt(i)
        write(5,'(I6, 3E15.7)')i, velx_pt(i), vely_pt(i), velz_pt(i)
500     continue
      do 505 i=1, nkoss_ar
        write(4,'(I6, 3E15.7)')i, posx_ar(i), posy_ar(i), posz_ar(i)
        write(6,'(I6, 3E15.7)')i, velx_ar(i), vely_ar(i), velz_ar(i)
505     continue
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine record2
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision(a-h,o-z)
      implicit integer(n)
      include 'inc/parameternano3.dat'
      include 'inc/variable4.dat' 
************************************************************************
      totene_pt = 0.00D0
      totpot_pt = 0.00D0
      totkin_pt = 0.00D0
      totpot_pt_top = 0.00D0		 
      totkin_pt_bot = 0.00D0
      totene_ar = 0.00D0
      totpot_ar = 0.00D0
      totkin_ar = 0.00D0

      do 5100 i=nhos+1,nkoss_pt_bot
        totpot_pt = totpot_pt + poten_pt(i)
        totkin_pt_bot = totkin_pt_bot + ukine_pt(i)
5100     continue

      do 801 i = nkoss_pt_bot+1, nkoss_pt-nhos
        totpot_pt = totpot_pt + poten_pt(i)
        totkin_pt_top = totkin_pt_top + ukine_pt(i)
801      continue

!         write(*,*)totkin_pt_bot,totkin_pt_top

      totpot_pt = totpot_pt/1.00D16
      totkin_pt_bot = totkin_pt_bot/1.00D16
      totkin_pt_top = totkin_pt_top/1.00D16
      totkin_pt = totkin_pt_bot + totkin_pt_top
      totene_pt = totpot_pt + totkin_pt
      temp_pt_bot = 2.000D0*totkin_pt_bot/
     &                     (3.00D0*dble(nkoss_pt_bot-nhos)*boltz)
      temp_pt_top = 2.000D0*totkin_pt_top/
     &                     (3.00D0*dble(3*nhos)*boltz)
      write(7,'(5E15.7)')totene_pt,totpot_pt,totkin_pt,
     &                      temp_pt_bot,temp_pt_top

      if(nowstp.ge.ns0)then
      ! 下壁面について
      ! 第2層目
        totkin_pt = 0.00D0
        do 6000 i=nhos+1,2*nhos
          totkin_pt = totkin_pt + ukine_pt(i)
      !			pt_bot2 = pt_bot2 + posz_pt(i)
6000      continue
      totkin_pt = totkin_pt/1.00D16
      temp_pt_bot2=2.00D0*totkin_pt/(3.00D0*dble(nhos)*boltz)
      tot_bot2 = tot_bot2 + temp_pt_bot2
    !		  pt_bot2 = pt_bot2/dble(nhos)
    !		  write(*,*)tot_bot2,temp_pt_bot2
    ! 第3層目 
      totkin_pt = 0.00D0
      do 6001 i=2*nhos+1,3*nhos
        totkin_pt = totkin_pt + ukine_pt(i)
      !		  pt_bot3 = pt_bot3 + posz_pt(i)
6001     continue
      totkin_pt = totkin_pt/1.00D16
      temp_pt_bot3=2.00D0*totkin_pt/(3.00D0*dble(nhos)*boltz)
      tot_bot3 = tot_bot3 + temp_pt_bot3
      !		 pt_bot3 = pt_bot3/dble(nhos)

      ! 第4層目
      totkin_pt = 0.00D0
      do 6002 i=3*nhos+1,4*nhos
        totkin_pt = totkin_pt + ukine_pt(i)
      !		  pt_bot4 = pt_bot4 + posz_pt(i)
6002     continue
      totkin_pt = totkin_pt/1.00D16
      temp_pt_bot4=2.00D0*totkin_pt/(3.00D0*dble(nhos)*boltz)	
      tot_bot4 = tot_bot4 + temp_pt_bot4
      !		 pt_bot4 = pt_bot4/dble(nhos)


      totkin_pt5 = 0.00D0
      totkin_pt6 = 0.00D0
      totkin_pt7 = 0.00D0
      totkin_pt8 = 0.00D0
      totkin_pt9 = 0.00D0
      totkin_pt10 = 0.00D0
      totkin_pt11 = 0.00D0
      totkin_pt12 = 0.00D0

      do 6100 i=nkoss_pt_flat+1,nkoss_pt_flat+nkoss_pt_nano
      ! 第5層目

        if(mod(i,nhos_nanoyz).ge.1)then
          if(mod(i,nhos_nanoyz).le.nny_pt)then
            totkin_pt5 = totkin_pt5+ukine_pt(i)
          endif
        endif
      !第6層目
        if(mod(i,nhos_nanoyz).ge.nny_pt+1)then
          if(mod(i,nhos_nanoyz).le.2*nny_pt)then
            totkin_pt6 = totkin_pt6+ukine_pt(i)
          endif
        endif
      !第7層目
        if(mod(i,nhos_nanoyz).ge.2*nny_pt+1)then
          if(mod(i,nhos_nanoyz).le.3*nny_pt)then
              totkin_pt7 = totkin_pt7+ukine_pt(i)
          endif
        endif
      !第8層目
        if(mod(i,nhos_nanoyz).ge.3*nny_pt+1)then
          if(mod(i,nhos_nanoyz).le.4*nny_pt)then
            totkin_pt8 = totkin_pt8+ukine_pt(i)
          endif
        endif

        if(mod(i,nhos_nanoyz).ge.4*nny_pt+1)then
          if(mod(i,nhos_nanoyz).le.5*nny_pt)then
            totkin_pt9 = totkin_pt9+ukine_pt(i)
          endif
        endif
        
        if(mod(i,nhos_nanoyz).ge.5*nny_pt+1)then
          if(mod(i,nhos_nanoyz).le.6*nny_pt)then
            totkin_pt10 = totkin_pt10+ukine_pt(i)
          endif
        endif
        
        if(mod(i,nhos_nanoyz).ge.6*nny_pt+1)then
          if(mod(i,nhos_nanoyz).le.7*nny_pt)then
            totkin_pt11 = totkin_pt11+ukine_pt(i)
          endif
        endif
        
        if(mod(i,nhos_nanoyz).ge.7*nny_pt+1)then
          if(mod(i,nhos_nanoyz).le.8*nny_pt)then
            totkin_pt12 = totkin_pt12+ukine_pt(i)
          endif
        endif
        
        if(mod(i,nhos_nanoyz).eq.0)then
            totkin_pt12 = totkin_pt12 + ukine_pt(i)
        endif		  
6100    continue

      totkin_pt5 = totkin_pt5/1.00D16
      totkin_pt6 = totkin_pt6/1.00D16
      totkin_pt7 = totkin_pt7/1.00D16
      totkin_pt8 = totkin_pt8/1.00D16
      totkin_pt9 = totkin_pt9/1.00D16
      totkin_pt10 = totkin_pt10/1.00D16
      totkin_pt11 = totkin_pt11/1.00D16
      totkin_pt12 = totkin_pt12/1.00D16
      temp_pt_bot5=2.00D0*totkin_pt5/
     &                   (3.00D0*dble(2*nhos_nanoxy)*boltz)
      temp_pt_bot6=2.00D0*totkin_pt6/
     &                   (3.00D0*dble(2*nhos_nanoxy)*boltz)
      temp_pt_bot7=2.00D0*totkin_pt7/
     &                   (3.00D0*dble(2*nhos_nanoxy)*boltz)
      temp_pt_bot8=2.00D0*totkin_pt8/
     &                   (3.00D0*dble(2*nhos_nanoxy)*boltz)
      temp_pt_bot9=2.00D0*totkin_pt9/
     &                   (3.00D0*dble(2*nhos_nanoxy)*boltz)
      temp_pt_bot10=2.00D0*totkin_pt10/
     &                   (3.00D0*dble(2*nhos_nanoxy)*boltz)
      temp_pt_bot11=2.00D0*totkin_pt11/
     &                   (3.00D0*dble(2*nhos_nanoxy)*boltz)
      temp_pt_bot12=2.00D0*totkin_pt12/
     &                   (3.00D0*dble(2*nhos_nanoxy)*boltz)
      tot_bot5 = tot_bot5+temp_pt_bot5
      tot_bot6 = tot_bot6+temp_pt_bot6
      tot_bot7 = tot_bot7+temp_pt_bot7
      tot_bot8 = tot_bot8+temp_pt_bot8
      tot_bot9 = tot_bot9+temp_pt_bot9
      tot_bot10 = tot_bot10+temp_pt_bot10
      tot_bot11 = tot_bot11+temp_pt_bot11
      tot_bot12 = tot_bot12+temp_pt_bot12

      !上壁面について
      ! 第2層目	
      totkin_pt = 0.00D0
      do 6003 i=nkoss_pt-2*nhos+1,nkoss_pt-nhos
        totkin_pt = totkin_pt + ukine_pt(i)
      !		  pt_top2 = pt_top2 + posz_pt(i)
6003     continue
      totkin_pt = totkin_pt/1.00D16
      temp_pt_top2 = 2.00D0*totkin_pt/(3.00D0*dble(nhos)*boltz)		
      tot_top2 = tot_top2 + temp_pt_top2
      !		 pt_top2 = pt_top2/dble(nhos)
      ! 第3層目
      totkin_pt = 0.00D0
      do 6004 i=nkoss_pt-3*nhos+1,nkoss_pt-2*nhos
        totkin_pt = totkin_pt + ukine_pt(i)
      !		  pt_top3 = pt_top3 + posz_pt(i)
6004     continue
      totkin_pt = totkin_pt/1.00D16
      temp_pt_top3 = 2.00D0*totkin_pt/(3.00D0*dble(nhos)*boltz)	
      tot_top3 = tot_top3 + temp_pt_top3
      !		 pt_top3 = pt_top3/dble(nhos)
      ! 第4層目
      totkin_pt = 0.00D0
      do 6005 i=nkoss_pt-4*nhos+1,nkoss_pt-3*nhos
        totkin_pt = totkin_pt + ukine_pt(i)
  
6005     continue
      totkin_pt = totkin_pt/1.00D16
      temp_pt_top4 = 2.00D0*totkin_pt/(3.00D0*dble(nhos)*boltz)	
      tot_top4 = tot_top4 + temp_pt_top4


      !		 write(42,'(5E15.7)')tot_bot2,tot_bot3,tot_bot4
      !	     write(43,'(5E15.7)')tot_top2,tot_top3,tot_top4

10       continue

! z座標の取得
      if(nowstp.eq.(maxstep-10000)+1)then
        do i=nhos+1,2*nhos
          pt_bot2=pt_bot2+posz_pt(i)
        enddo
        do i=2*nhos+1,3*nhos
          pt_bot3=pt_bot3+posz_pt(i) 
        enddo
        do i=3*nhos+1,nkoss_pt_flat
          pt_bot4=pt_bot4+posz_pt(i)
        enddo
        
        do i=nkoss_pt_flat+1,nkoss_pt_bot
          j= i-nkos_pt_flat
          if(mod(i,nhos_nanoyz).ge.1)then
            if(mod(i,nhos_nanoyz).le.nny_pt)then
              pt_bot5= pt_bot5+posz_pt(i)
            endif
          endif
          if(mod(i,nhos_nanoyz).ge.nny_pt+1)then
            if(mod(i,nhos_nanoyz).le.2*nny_pt)then
              pt_bot6= pt_bot6+posz_pt(i)
            endif
          endif
          if(mod(i,nhos_nanoyz).ge.2*nny_pt+1)then
            if(mod(i,nhos_nanoyz).le.3*nny_pt)then
              pt_bot7= pt_bot7+posz_pt(i)
            endif
          endif
          if(mod(i,nhos_nanoyz).ge.3*nny_pt+1)then
            if(mod(i,nhos_nanoyz).le.4*nny_pt)then
              pt_bot8= pt_bot8+posz_pt(i)
            endif
          endif
          if(mod(i,nhos_nanoyz).ge.4*nny_pt+1)then
            if(mod(i,nhos_nanoyz).le.5*nny_pt)then
              pt_bot9= pt_bot9+posz_pt(i)
            endif
          endif
          if(mod(i,nhos_nanoyz).ge.5*nny_pt+1)then
            if(mod(i,nhos_nanoyz).le.6*nny_pt)then
              pt_bot10= pt_bot10+posz_pt(i)
            endif
          endif
          if(mod(i,nhos_nanoyz).ge.6*nny_pt+1)then
            if(mod(i,nhos_nanoyz).le.7*nny_pt)then
              pt_bot11= pt_bot11+posz_pt(i)
            endif
          endif
          if(mod(i,nhos_nanoyz).ge.7*nny_pt+1)then
            if(mod(i,nhos_nanoyz).le.8*nny_pt)then
              pt_bot12= pt_bot12+posz_pt(i)
            endif
          endif
          if(mod(i,nhos_nanoyz).eq.0)then
              pt_bot12 = pt_bot12+posz_pt(i)
          endif
        enddo		  
        
        do i=nkoss_pt-2*nhos+1,nkoss_pt-nhos
          pt_top2=pt_top2+posz_pt(i)
        enddo
        do i=nkoss_pt-3*nhos+1,nkoss_pt-2*nhos
          pt_top3=pt_top3+posz_pt(i)
        enddo
        do i=nkoss_pt-4*nhos+1,nkoss_pt-3*nhos
          pt_top4=pt_top4+posz_pt(i)
        enddo
        
        pt_bot2=pt_bot2/dble(nhos)/dble(10)
        pt_bot3=pt_bot3/dble(nhos)/dble(10)
        pt_bot4=pt_bot4/dble(nhos)/dble(10)
        
        
        pt_bot5=pt_bot5/dble(2*nhos_nanoxy)/dble(10)
        pt_bot6=pt_bot6/dble(2*nhos_nanoxy)/dble(10)
        pt_bot7=pt_bot7/dble(2*nhos_nanoxy)/dble(10)
        pt_bot8=pt_bot8/dble(2*nhos_nanoxy)/dble(10)
        pt_bot9=pt_bot9/dble(2*nhos_nanoxy)/dble(10)
        pt_bot10=pt_bot10/dble(2*nhos_nanoxy)/dble(10)
        pt_bot11=pt_bot11/dble(2*nhos_nanoxy)/dble(10)
        pt_bot12=pt_bot12/dble(2*nhos_nanoxy)/dble(10)
        
        
        pt_top2=pt_top2/dble(nhos)/dble(10)
        pt_top3=pt_top3/dble(nhos)/dble(10)
        pt_top4=pt_top4/dble(nhos)/dble(10)
        
        
        tot_bot2=tot_bot2/(dble(maxstep-ns0)/dble(10000))
        tot_bot3=tot_bot3/(dble(maxstep-ns0)/dble(10000))
        tot_bot4=tot_bot4/(dble(maxstep-ns0)/dble(10000))
        tot_bot5=tot_bot5/(dble(maxstep-ns0)/dble(10000))
        tot_bot6=tot_bot6/(dble(maxstep-ns0)/dble(10000))
        tot_bot7=tot_bot7/(dble(maxstep-ns0)/dble(10000))
        tot_bot8=tot_bot8/(dble(maxstep-ns0)/dble(10000))
        tot_bot9=tot_bot9/(dble(maxstep-ns0)/dble(10000))
        tot_bot10=tot_bot10/(dble(maxstep-ns0)/dble(10000))
        tot_bot11=tot_bot11/(dble(maxstep-ns0)/dble(10000))
        tot_bot12=tot_bot12/(dble(maxstep-ns0)/dble(10000))
        
        tot_top2=tot_top2/(dble(maxstep-ns0)/dble(10000))
        tot_top3=tot_top3/(dble(maxstep-ns0)/dble(10000))
        tot_top4=tot_top4/(dble(maxstep-ns0)/dble(10000))
      !		  write(*,*)(dble(maxstep)/dble(1000))
        write(42,'(2E15.7)')pt_bot2,tot_bot2
        write(42,'(2E15.7)')pt_bot3,tot_bot3
        write(42,'(2E15.7)')pt_bot4,tot_bot4
        write(42,'(2E15.7)')pt_bot5,tot_bot5
        write(42,'(2E15.7)')pt_bot6,tot_bot6
        write(42,'(2E15.7)')pt_bot7,tot_bot7
        write(42,'(2E15.7)')pt_bot8,tot_bot8
        write(42,'(2E15.7)')pt_bot9,tot_bot9
        write(42,'(2E15.7)')pt_bot10,tot_bot10
        write(42,'(2E15.7)')pt_bot11,tot_bot11
        write(42,'(2E15.7)')pt_bot12,tot_bot12		   
        write(42,'(2E15.7)')pt_top4,tot_top4
        write(42,'(2E15.7)')pt_top3,tot_top3
        write(42,'(2E15.7)')pt_top2,tot_top2
        write(700,'(2E15.7)')nhos_nanoxy
        endif
      endif

                
      do 900 i900 = 1, nkoss_ar
        totpot_ar = totpot_ar + poten_ar(i900)
        totkin_ar = totkin_ar + ukine_ar(i900)
900      continue
        totpot_ar = totpot_ar/1.00D16
        totkin_ar = totkin_ar/1.00D16
        totene_ar = totpot_ar + totkin_ar
      temp_ar = 2.000D0*totkin_ar/(3.00D0*dble(nkoss_ar)*boltz)
      write(8,'(4E15.7)')totene_ar, totpot_ar, totkin_ar, temp_ar

      totene = totene_pt + totene_ar
      totpot = totpot_pt + totpot_ar
      totkin = totkin_pt + totkin_ar

      write(9,'(4E15.7)')totene,totpot,totkin
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine write_pos_vel
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision(a-h,o-z)
      implicit integer(n)
      include 'inc/parameternano3.dat'
      include 'inc/variable4.dat'
************************************************************************
      num = 0
      write(10,'(I14)')nowstp
      do 400 i=1, nkoss_pt
      write(10,'(I6, 6D15.7)')i,posx_pt(i), posy_pt(i), posz_pt(i),
     &                              velx_pt(i), vely_pt(i), velz_pt(i)
400     continue

      do 402 i=1, nkoss_ar
      write(10,'(I6, 6D15.7)')i,posx_ar(i), posy_ar(i), posz_ar(i),
     &                              velx_ar(i), vely_ar(i), velz_ar(i)
402     continue
      write(10,'(2E24.16)')vmax_pt,vmax_ar
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccc	  
      subroutine pvch
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision(a-h,o-z)
      implicit integer(n)
      include 'inc/parameternano3.dat'
      include 'inc/variable4.dat'
      real    pomx,pomy,pomz
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if (nowstp .eq. 1) then
      write(17,'(3I7)') moltype,nmol,ndat1
      write(17,'(3F15.5)') xsyul0,ysyul0,zsyul0
      write(17,'(2I7)') ntime0, ndt
      endif
cccccccccccccccccccccccccccccccccccccccccccccccc
      do 403 j=1, nkoss_pt
        pomx = posx_pt(j)
        pomy = posy_pt(j)
        pomz = posz_pt(j)
        write(17,'(3E15.7)') pomx, pomy, pomz
403     continue
!        write(*,*)nkoss_pt

      do 404 j=1, nkoss_ar
        pomx = posx_ar(j)
        pomy = posy_ar(j)
        pomz = posz_ar(j)
        write(17,'(3E15.7)') pomx, pomy, pomz
404     continue

      do 750 i75=1,nkoss_pt+nkoss_ar
        if (i75 .le. nkoss_pt)then
          write(18,'(I7)')nspt
      !			write(*,*)nspt
        else
          write(18,'(I7)')nsar
      !			write(*,*)nsar
        endif
750     continue

      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine bunpu     !温度分布と密度分布の算出
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision(a-h,o-z)
      implicit integer(n)
      include 'inc/parameternano3.dat'
      include 'inc/variable4.dat'
************************************************************************
      underz_ar = 0           !系内の下限
      topz_ar = zsyul0     !系内の上限
      z_ar = topz_ar - underz_ar    !系内のz方向長さ
      dzd = z_ar / dble(ndzd)  !密度分布のz方向の幅
      dzt = z_ar / dble(ndzt)  !温度分布のz方向の幅
      ncountd_ar = 0
      ncount_ar = 0
      ncount_pt = 0
      slab_v = xsyul0*ysyul0*dzd   !スラブの体積

      ! Arの一次元密度分布		
      do 20 i=1,ndzd
        z_up = i*dzd + underz_ar              !スラブのz方向の上限
        z_down = (i-1)*dzd + underz_ar        !スラブのz方向の下限
        do 21 j=1,nkoss_ar
          if(posz_ar(j).ge.z_down)then
            if(posz_ar(j).lt.z_up)then
              ncountd_ar(i) = ncountd_ar(i) + 1   !領域内のAr粒子の個数を求める
            endif
          endif
21        continue
!          totkin_ar_slab(i) = totkin_ar_slab(i) / 1.00D16
!		  temp_ar_slab(i) = 2.00D0*totkin_ar_slab(i)/
!     &                     (3.00D0*dble(ncount_ar(i))*boltz)
!	      tot_temp_ar(i) = tot_temp_ar(i) + temp_ar_slab(i)
      zdensi(i) = zmass_ar*dble(ncountd_ar(i))/slab_v
!		  write(*,*)zdensi(i)
      tot_zdensi(i) = tot_zdensi(i) + zdensi(i)
!		  write(*,*) tot_temp_ar(i),tot_zdensi(i)

! Ptの一次元密度分布
      do 22 k=1,nkoss_pt
        if(posz_pt(k).ge.z_down)then
          if(posz_pt(k).lt.z_up)then
            ncount_pt(i)=ncount_pt(i)+1
          endif
        endif
22        continue
      zdensi_pt(i)=zmass_pt*dble(ncount_pt(i))/slab_v
      tot_zdensi_pt(i)=tot_zdensi_pt(i)+zdensi_pt(i)

      if(nowstp.eq.maxstep)then
        dz(i) = dzd*dble(i)/10    !nm		
      !		  tot_temp_ar(i) = tot_temp_ar(i)/dble(maxstep-ns0)
        tot_zdensi(i) = tot_zdensi(i)*1.0D04/dble(maxstep-ns0)
        tot_zdensi_pt(i) = tot_zdensi_pt(i)*1.0D04/dble(maxstep-ns0)
        write(41,'(4E15.7)')dz(i),tot_zdensi(i),tot_zdensi_pt(i)
        
!		  write(32,*)dble(i)*dz,ncount_ar(i)
      endif		  
20      continue

      do 30 i=1,ndzt
        z_up = i*dzt + underz_ar              !スラブのz方向の上限
        z_down = (i-1)*dzt + underz_ar        !スラブのz方向の下限
        do 31 j=1,nkoss_ar
          if(posz_ar(j).ge.z_down)then
            if(posz_ar(j).lt.z_up)then
              ncount_ar(i) = ncount_ar(i) + 1   !領域内のAr粒子の個数を求める
              totkin_ar_slab(i) = totkin_ar_slab(i) + ukine_ar(j)
      !				write(*,*) '1'
            endif
          endif
31        continue
        totkin_ar_slab(i) = totkin_ar_slab(i) / 1.00D16
        temp_ar_slab(i) = 2.00D0*totkin_ar_slab(i)/
     &                     (3.00D0*dble(ncount_ar(i))*boltz)
        tot_temp_ar(i) = tot_temp_ar(i) + temp_ar_slab(i)
!	      zdensi(i) = zmass_ar*dble(ncount_ar(i))/slab_v
!          tot_zdensi(i) = tot_zdensi(i) + zdensi(i)
!		  write(*,*) tot_temp_ar(i),tot_zdensi(i)

      if(nowstp.eq.maxstep)then	
        dz(i) = dzt*dble(i)/10		  
        tot_temp_ar(i) = tot_temp_ar(i)/dble(maxstep-ns0)
    !		    tot_zdensi(i) = tot_zdensi(i)/dble(maxstep-ns0)
        write(40,'(3E15.7)')dz(i),tot_temp_ar(i)
    !		    write(32,*)dble(i)*dz,ncount_ar(i)
      endif		  
30      continue		
      return 
      end			
cccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine calcu_ddl
cccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision(a-h,o-z)
      implicit integer(n)
      include 'inc/parameternano3.dat'
      include 'inc/variable4.dat'
cccccccccccccccccccccccccccccccccccccccccccccccc      
      dens_ar = tot_zdensi(ndzd/2) !バルクの液体密度

      do 2 i2=1,600
      d=(1-(tot_zdensi_pt(i2)/dens_pt)
     &             -(tot_zdensi(i2)/dens_ar))
      z_d=(i2-1)*dzd
      ddl = ddl+(1-(tot_zdensi_pt(i2)/dens_pt)
     &             -(tot_zdensi(i2)/dens_ar))*dzd
      write(200,'(2E15.7)')z_d,d/(dble(10))
      !	  write(600,'(2D15.7)')tot_zdensi_pt(i2),tot_zdensi(i2)
2     continue



      ddl = ddl/dble(10)

      write(60,'(3E15.7)')ddl
      return
      end


ccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine record_maxwell
ccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision(a-h,o-z)
      implicit integer(n)
      include 'inc/parameternano3.dat'
      include 'inc/variable4.dat'
      double precision vel_maxwell(nkoss_ar)
ccccccccccccccccccccccccccccccccccccccccccccccc

      do 100 i=1,nkoss_ar
      vel_maxwell(i)=dsqrt(velx_ar(i)**2+vely_ar(i)**2+velz_ar(i)**2)
      vel_maxwell(i)=vel_maxwell(i) * (10**5)
100   continue

      do 200 i=1,nmaxwell
      do 201 j=1,nkoss_ar
        if(vel_maxwell(j).lt.dble(i)*dt_maxwell) then
          if(vel_maxwell(j).ge.dble(i-1)*dt_maxwell) then
            ncount(i) = ncount(i)+1
      !			  write(*,*)ncount(i)
          endif
        endif
201     continue
200   continue

      if(nowstp.ge.maxstep-1000)then
      write(20,*) (ncount(i),i=1,nmaxwell)
      endif

      return 
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine heatflux      !相互作用による熱流束の算出
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision(a-h,o-z)
      implicit integer(n)
      include 'inc/parameternano3.dat'
      include 'inc/variable4.dat'
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  初期化		  
      do 30 i=1,nkoss_pt
          qx_pt(i) = 0.00D0
          qy_pt(i) = 0.00D0
          qz_pt(i) = 0.00D0
          q_pt(i) = 0.00D0
30      continue

! 各Pt分子の熱量の総和
      do 40 i=1,nkoss_pt
        qx_pt(i) = qx_pt(i) + force_qx(i)*velx_lan(i)*dt		
        qy_pt(i) = qy_pt(i) + force_qy(i)*vely_lan(i)*dt
        qz_pt(i) = qz_pt(i) + force_qz(i)*velz_lan(i)*dt
        q_pt(i) = qx_pt(i) + qy_pt(i) + qz_pt(i)
      !		  write(*,*)q_pt(i)
40      continue

! 各壁面の伝熱量測定
      do 50 i=1,nkoss_pt_bot
      tot_q_bot = tot_q_bot + q_pt(i)
          if(mod(nowstp,1000).eq.0)then
          write(500,'(2D15.7)')tot_q_bot,q_pt(i)
        endif
50      continue

      do 51 j=nkoss_pt_bot+1,nkoss_pt
      tot_q_top = tot_q_top + q_pt(j)
        if(mod(nowstp,1000).eq.0)then
          write(501,'(2D15.7)')tot_q_top,q_pt(j)
        endif
51      continue

      q_bot = tot_q_bot*1.00D4/(xsyul0*ysyul0+2*surface)*(-1.0D0)
      q_top = tot_q_top*1.00D4/(xsyul0*ysyul0)*(-1.0D0)
      if(mod(nowstp,1000).eq.0)then
      tnowtime = dble(nowstp)*dt*1.0D-6 !現在時刻(ns)
      write(50,'(3E15.7)')tnowtime,q_bot,q_top
      endif		
      return
      end		

cccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine heatflux_lange
cccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision(a-h,o-z)
      implicit integer(n)
      include 'inc/parameternano3.dat'
      include 'inc/variable4.dat'
cccccccccccccccccccccccccccccccccccccccccccccccccc
!  初期化
      do 40 i=1,nkoss_pt
        qx_pt(i) = 0.00D0
        qy_pt(i) = 0.00D0
        qz_pt(i) = 0.00D0
        q_pt(i) = 0.00D0
40      continue

      do 10 i=1,nkoss_pt
      qx_pt(i) = qx_pt(i)+flanx(i)*velx_lan(i)*dt
      qy_pt(i) = qy_pt(i)+flany(i)*vely_lan(i)*dt
      qz_pt(i) = qz_pt(i)+flanz(i)*velz_lan(i)*dt
      q_pt(i) = qx_pt(i)+qy_pt(i)+qz_pt(i)
10      continue

      do 20 i=nhos+1,2*nhos
      qx_bot_L = qx_bot_L + qx_pt(i)
      qy_bot_L = qy_bot_L + qy_pt(i)
      qz_bot_L = qz_bot_L + qz_pt(i)
      q_bot_L = q_bot_L + q_pt(i)
20      continue

      do 21 j=nkoss_pt-2*nhos+1,nkoss_pt-nhos
      qx_top_L = qx_top_L + qx_pt(j)
      qy_top_L = qy_top_L + qy_pt(j)
      qz_top_L = qz_top_L + qz_pt(j)
      q_top_L = q_top_L + q_pt(j)
21      continue

      q_bot_phantom = q_bot_L*1.0D4/(xsyul0*ysyul0)
      q_top_phantom = q_top_L*1.0D4/(xsyul0*ysyul0)
      if(mod(nowstp,1000).eq.0)then
      tnowtime=dble(nowstp)*dt*1.0D-6 !現在時刻(ns)
      write(51,'(3E15.7)')tnowtime,q_bot_phantom,q_top_phantom
      endif
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine density_2d
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit double precision(a-h,o-z)
      implicit integer(n)
      include 'inc/parameternano3.dat'
      include 'inc/variable4.dat'            
ccccccccccccccccccccccccccccccccccccccccccccccccccccc		
      do 1 i=1,ndx
      do 2 j=1,ndz
      n_dens(i,j)=0
2       continue
1     continue

!パラメータ
      dzz = zsyul0/dble(ndz)
      dxx = xsyul0/dble(ndx)
      !	  n_dens=0
      !	  density_max=0
      !	  name_max_group = 15
      slab_2d_v = dxx*ysyul0*dzz

      ! 領域内のAr粒子の個数を求める	  
      do 10 i=1,nkoss_ar
      !	    name_b_group(i) = name_group(i)
      !		name_group(i) = 0
      do 11 j=1,ndx
      if(posx_ar(i).ge.dxx*dble(j-1))then
        if(posx_ar(i).lt.dxx*dble(j))then
          do 12 k=1,ndz
          if(posz_ar(i).ge.dzz*dble(k-1))then
            if(posz_ar(i).lt.dzz*dble(k))then
              n_dens(j,k) = n_dens(j,k)+1
              go to 10
            endif
          endif
12            continue
        endif
      endif
11      continue
10    continue

! 領域内のArの総数
      do 100 i=1,ndx
      do 101 j=1,ndz
      n_dens_tot(i,j)=n_dens_tot(i,j)+n_dens(i,j)
101     continue
100   continue

! メッシュ内のArの密度
      if(nowstp.eq.maxstep)then
      do 200 i=1,ndx
      do 201 j=1,ndz
        dens_2d(i,j)=dble(n_dens_tot(i,j))*zmass_ar*1.0D04/
     &                      slab_2d_v/dble(maxstep-ns0)
        x_2d(i) = (dble(i)-0.5)*dxx/dble(10)
        z_2d(j) = (dble(j)-0.5)*dzz/dble(10)
        time = dble(nowstp)*dt*1.0D-6
        write(70,'(4E15.7)')time,x_2d(i),z_2d(j),dens_2d(i,j)
201       continue
200     continue
      endif

      return 
      end
        