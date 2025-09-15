      SUBROUTINE FLUXAV_VMEC(nr_r,rhot_r,a0,f0,dvol_r,vp_r,shift_r,
     &                       elong_r,triang_r,gth_r,gph_r,b2_r,bm2_r,
     &                       grho2_r,gr2bm2_r,grth_r,ftrap_r,fm_r,f_r,
     &                       q_r,psi_r,iflag,message)
!***********************************************************************
!FLUXAV_VMEC calculates flux surface quantities from VMEC equilibria
!References:
!  Houlberg, Shaing, Hirshman, Zarnstorff, Phys Plasmas 4 (1997) 3230
!  Lin-Liu, Miller, Phys Plasmas 2 (1995) 1666
!  W.A.Houlberg 2/99
!Input:
!  nr_r-radial nodes in rhot_r grid (-)
!  rhot_r(i)-radial (cell interface) grid prop sqrt(toroidal flux) (m)
!  a0-minor radius (m)
!  f0-poloidal current external to plasma, f_r(nr_r)=f0 (A)
!Output:
!  dvol_r(i)-volume of cell i (m**3)
!  vp_r(i)-dV/d(rhot_r) at rhot_r(i) (m**2)
!  shift_r(i)-shift of surface rhot_r(i) (-)
!  elong_r(i)-elongation of surface rhot_r(i) (-)
!  triang_r(i)-triangularity of surface rhot_r(i) (-)
!  gth_r(i)-<gtt/sqrt(g)>-theta average (-)
!  gph_r(i)-<sqrt(g)/R**2>-theta average (-)
!  b2_r(i)-<B**2>-flux surface average (T**2)
!  bm2_r(i)-<1/B**2>-flux surface average (/T**2)
!  grho2_r(i)-<grad(rho)**2>-flux surface average (-)
!  gr2bm2_r(i)-<grad(rho)**2/B**2>-flux surface average (/T**2)
!  grth_r(i)-<n.grad(theta)>-flux surface average (/m)
!  ftrap_r(i)-trapped particle fraction (-)
!  fm_r(3,i)-poloidal moments of metric for viscosity (-)
!  f_r(i)-poloidal current external to rhot_r(i) (A)
!  q_r(i)-safety factor at surface rhot_r(i) (-)
!  psi_r(i)-poloidal flux at rhot_r(i) (Wb)
!  iflag-error flag for magnetic field calculation (-)
!       =0 or 1, normal exit
!       >1 no convergence in TRAFLX
!Comments:
!  Metrics in the scrape-off region are obtained by linear extrapolation
!  Trapped fraction uses approximation by Lin-Liu
!***********************************************************************
      IMPLICIT NONE
      INCLUDE 'pamxkm.inc'
      INCLUDE 'pamxmp.inc'
      INCLUDE 'pamxmt.inc'
      INCLUDE 'pamxnc.inc'
      INCLUDE 'pamxnx.inc'
      INCLUDE 'comm30.inc'
      INCLUDE 'comm31.inc'
!Declaration of input variables
      INTEGER        nr_r
      REAL           a0,                      f0
      REAL           rhot_r(*)
!Declaration of output variables
      CHARACTER*(*)  message
      INTEGER        iflag
      REAL           bm2_r(*),                b2_r(*),
     &               dvol_r(*),               elong_r(*),
     &               f_r(*),                  ftrap_r(*),
     &               gph_r(*),                gr2bm2_r(*),
     &               grho2_r(*),              grth_r(*),
     &               gth_r(*),                psi_r(*),
     &               q_r(*),                  shift_r(*),
     &               triang_r(*),             vp_r(*)
      REAL           fm_r(3,*)
!Declaration of local variables
      INTEGER        i,                       iph,
     &               ith,                     j,
     &               n
      REAL           ar2,                     as2,
     &               b2,                      br,
     &               btor,                    bz,
     &               con1,                    con2,
     &               dph,                     dth,
     &               ph,                      r,
     &               rg0,                     rmax,
     &               rmin,                    rp,
     &               rs1,                     rs2,
     &               rt,                      rx,
     &               scal,                    th,
     &               th1,                     th2,
     &               tolinv,                  xd,
     &               xx,                      z,
     &               zidum,                   zp,
     &               zs1,                     zs2,
     &               zt,                      zt1,
     &               zt2,                     zx
      REAL           bmax(mxnx+1),            el3d(mxnx+1),
     &               fpol(mxnx+1),            psij(mxnx+1),
     &               qj(mxnx+1),              rm3d(mxnx+1),
     &               sh3d(mxnx+1),            tr3d(mxnx+1),
     &               yd3d(mxnx+1),            zd3d(mxnx+1)
      REAL           r2bph(mxmp)
      REAL           b1(mxmt),                bth(mxmt),
     &               capth(mxmt)
      REAL           atg1(mxmp,18)
      REAL           g1(mxmt,18)
      REAL           atpg1(mxnx+1,18)
      REAL           zdum1(mxnc1),            zdum2(mxnc1)
      REAL           z_mu0,                   z_pi,
     &               z_precision
!Declaration of external functions
      REAL           RARRAY_XDOTY
!Initialization
!  Error flag
      iflag=0
      message=' '
!  Physical and conversion constants
      z_mu0=1.2566e-06
      z_pi=ACOS(-1.0)
      z_precision=2.0e-7
      tolinv=1.0e-5*a0
!Find the geometric center of plasma boundary, rg0
      CALL TRAGRZ(1.0,0.0,0.0,r,z,rx,rt,rp,zx,zt,zp)
      rg0=0.5*r
      CALL TRAGRZ(1.0,z_pi,0.0,r,z,rx,rt,rp,zx,zt,zp)
      rg0=rg0+0.5*r
!Set radial grid
      CALL RARRAY_COPY(nx3d,xm3d,1,rm3d,1)
      CALL RARRAY_SCALE(nx3d,a0,rm3d,1)
!  Check whether there is a scrape-off layer
      IF(rhot_r(nr_r).gt.a0*(1.0+10.0*z_precision)) THEN
        n=nx3d+1
        rm3d(n)=rhot_r(nr_r)
      ELSE
        n=nx3d
      ENDIF
!Compute metric coefficients on MHD grid
!  Set toroidal angle increment
      IF(mp3d.eq.1) THEN
        dph=0.0
      ELSE
        dph=2.0*z_pi/(nper*FLOAT(mp3d-1))
      ENDIF
!  Set poloidal angle increment
      dth=2.0*z_pi/FLOAT(mt3d-1)
!  Loop over radial nodes
      DO j=2,n
        xx=rm3d(j)/a0
        bmax(j)=0.0
!  Loop over toroidal angle
        DO iph=1,mp3d
          ph=(iph-1)*dph
!  Loop over poloidal angle
          capth(1)=0.0
          DO ith=1,mt3d
            th=(ith-1)*dth
!  Get the R and Z coordinates
            CALL TRAGRZ(xx,th,ph,r,z,rx,rt,rp,zx,zt,zp)
!  Geometric quantities
            g1(ith,1)=r*(rt*zx-rx*zt)/a0
            g1(ith,2)=r**2*(rt**2+zt**2)/g1(ith,1)
            g1(ith,3)=(rt**2+zt**2)/g1(ith,1)
            g1(ith,4)=g1(ith,1)/r**2
!    Get magnetic field components
            IF(xx.le.1.0+10.0*z_precision) THEN
!    Use pi/2 for poloidal angle so that d(lambda)/d(theta) is near zero
              CALL TRABFL(0,r,z,ph,tolinv,br,bz,btor,xx,th,rx,rt,
     &                    rp,zx,zt,zp,zidum,con1,con2,iflag,message)
              IF(iflag.eq.1) THEN
                GOTO 1000
              ELSEIF(iflag.eq.0) THEN
                b2=br**2+bz**2+btor**2
                b1(ith)=SQRT(b2)
                bth(ith)=r*(zx*br-rx*bz)/a0/g1(ith,1)
                g1(ith,5)=g1(ith,1)*b1(ith)
                g1(ith,6)=g1(ith,1)*b2
                g1(ith,7)=b1(ith)/bth(ith)
                g1(ith,8)=g1(ith,1)*bth(ith)
                g1(ith,17)=g1(ith,1)/b2
                g1(ith,18)=g1(ith,2)/b2
!    Get the maximum magnetic field on this surface
                IF(b1(ith).gt.bmax(j)) bmax(j)=b1(ith)
!    Get the transformation angle, Theta
                IF(ith.gt.1) capth(ith)=capth(ith-1)
     &                                  +(g1(ith-1,7)+g1(ith,7))*dth/2.0
              ENDIF
            ENDIF
          ENDDO
!    Scale transformation angle
          scal=2.0*z_pi/capth(mt3d)
          CALL RARRAY_SCALE(mt3d,scal,capth,1)  
          IF(xx.le.1.0+10.0*z_precision) THEN
!    Set quantities involving grad(B)
            DO ith=1,mt3d   
              th=(ith-1)*dth
              con1=b1(ith)/bmax(j)
!    Term for Lin-Liu expression for upper bound on ftrap
              IF((1.0-con1).gt.0.0) THEN
                g1(ith,10)=g1(ith,1)
     &                     *(1.0-SQRT(1.0-con1)*(1.0+con1/2.0))/con1**2
              ELSE
                g1(ith,10)=g1(ith,1)
              ENDIF
              IF(ith.eq.1.or.ith.eq.mt3d) THEN
!               Get d(B)/d(theta) at end points noting periodicity
                con2=(b1(2)-b1(mt3d-1))/(2.0*dth)
              ELSE
                con2=(b1(ith+1)-b1(ith-1))/(2.0*dth)
              ENDIF
              con1=bth(ith)*con2/b1(ith)
              g1(ith,9)=g1(ith,1)*con1**2
              g1(ith,11)=g1(ith,1)*SIN(capth(ith))*con1
              g1(ith,14)=g1(ith,11)*b1(ith)*scal
              g1(ith,12)=g1(ith,1)*SIN(2.0*capth(ith))*con1
              g1(ith,15)=g1(ith,12)*b1(ith)*scal
              g1(ith,13)=g1(ith,1)*SIN(3.0*capth(ith))*con1
              g1(ith,16)=g1(ith,13)*b1(ith)*scal
            ENDDO   
          ENDIF
!    Average over poloidal angle
          DO i=1,18
            atg1(iph,i)=RARRAY_XDOTY(mt3d,thfcn,1,g1(1,i),1)/thnorm
          ENDDO          
          IF(j.le.nx3d) THEN
            IF(j.eq.nx3d) THEN
!             Insure that evaluation is inside plasma boundary
              xd=1.0-2.0*tolinv/a0
            ELSE
              xd=xx
            ENDIF
!    Get R and Z coordinates of the point
            CALL TRAGRZ(xd,0.5*z_pi,ph,r,z,rx,rt,rp,zx,zt,zp)
!    Get toroidal field component of the point
!    Use pi/2 for poloidal angle so that d(lambda)/d(theta) is near zero
            iflag=0
            CALL TRABFL(0,r,z,ph,tolinv,br,bz,btor,xd,0.5*z_pi,rx,
     &                  rt,rp,zx,zt,zp,zidum,con1,con2,iflag,message)
            IF(iflag.eq.0) THEN
              r2bph(iph)=r*btor
            ELSE
!             Point not found inside plasma
              r2bph(iph)=0.0
              GOTO 1000
            ENDIF
          ELSE
            r2bph(iph)=0.0
          ENDIF
        ENDDO   
!  Average over toroidal angle
        DO i=1,18
          atpg1(j,i)=RARRAY_XDOTY(mp3d,phfcn,1,atg1(1,i),1)/phnorm
        ENDDO   
        IF(j.le.nx3d) THEN
          fpol(j)=2.0*z_pi*RARRAY_XDOTY(mp3d,phfcn,1,r2bph,1)
     &            /(phnorm*z_mu0)
        ELSE
          fpol(j)=fpol(nx3d)
        ENDIF
      ENDDO   
!Set values at origin
      DO i=1,18
        atpg1(1,i)=0.0
      ENDDO   
      fpol(1)=fpol(2)-xm3d(2)*(fpol(3)-fpol(2))/(xm3d(3)-xm3d(2))
      g1(1,2)=g1(2,2)
      g1(1,5)=g1(2,5)
      g1(1,6)=g1(2,6)
      g1(1,10)=g1(2,10)
      g1(1,17)=g1(2,17)
      g1(1,18)=g1(2,18)
!Set values at edge if SOL
      IF(n.gt.nx3d) THEN
        bmax(n)=bmax(n-1)
      ENDIF
!Normalize external poloidal current to tf coil current
      CALL RARRAY_SCALE(n,f0/fpol(nx3d),fpol,1)
!Other dependent variables and metrics
      psij(1)=0.0
      qj(1)=1.0/fiota(1,1)
      DO j=2,n
        xx=rm3d(j)/a0
!  Local minor radius, geometric center, and shift
        CALL TRAGRZ(xx,0.0,0.0,r,z,rx,rt,rp,zx,zt,zp)
        rmax=r
        rmin=r
        CALL TRAGRZ(xx,z_pi,0.0,r,z,rx,rt,rp,zx,zt,zp)
        IF(rmax.lt.r) rmax=r
        IF(rmin.gt.r) rmin=r
        as2=0.5*(rmax-rmin)
        ar2=0.5*(rmax+rmin)   
        sh3d(j)=(ar2-rg0)/a0
!  Find where dz/d(theta)=0 to get elongation and triangularity
!  Calculation allows for vertical asymmetry
!       Bottom 
        th1=0.6*z_pi
        CALL TRAGRZ(xx,th1,0.0,rs1,zs1,rx,rt,rp,zx,zt1,zp)
        th2=th1-0.1*z_pi
        CALL TRAGRZ(xx,th2,0.0,rs2,zs2,rx,rt,rp,zx,zt2,zp)
!       Normally 2 iterations is sufficient
        DO WHILE(ABS(zt2*2.0*z_pi/a0).gt.0.3*tolinv)
          con1=th2-zt2*(th2-th1)/(zt2-zt1)
          th1=th2
          rs1=rs2
          zs1=zs2
          zt1=zt2
          th2=con1
          CALL TRAGRZ(xx,th2,0.0,rs2,zs2,rx,rt,rp,zx,zt2,zp)
        ENDDO
        el3d(j)=-zs2/as2
        tr3d(j)=(ar2-rs2)/as2
!       Top
        th1=-0.6*z_pi
        CALL TRAGRZ(xx,th1,0.0,rs1,zs1,rx,rt,rp,zx,zt1,zp)
        th2=th1+0.1*z_pi
        CALL TRAGRZ(xx,th2,0.0,rs2,zs2,rx,rt,rp,zx,zt2,zp)
!       Normally 2 iterations is sufficient
        DO WHILE(ABS(zt2*2.0*z_pi/a0).gt.0.3*tolinv)
          con1=th2-zt2*(th2-th1)/(zt2-zt1)
          th1=th2
          rs1=rs2
          zs1=zs2
          zt1=zt2
          th2=con1
          CALL TRAGRZ(xx,th2,0.0,rs2,zs2,rx,rt,rp,zx,zt2,zp)
        ENDDO   
        el3d(j)=0.5*ABS(el3d(j)+zs2/as2)
        tr3d(j)=0.5*(tr3d(j)+(ar2-rs2)/as2)
!  Safety factor and poloidal flux
        IF(j.le.nx3d) THEN
!         Inside plasma
          qj(j)=1.0/fiota(1,j)
          psij(j)=psij(j-1)
     &            +(fiota(1,j-1)*xm3d(j-1)+fiota(1,j)*xx)*(xx-xm3d(j-1))
     &            *2.0*phitot
        ELSE
!         Outside plasma
          qj(j)=1.0/fiota(1,nx3d)
          psij(j)=psij(j-1)+(fiota(1,j-1)*xm3d(j-1)+fiota(1,nx3d)*xx)
     &            *(xx-xm3d(j-1))*2.0*phitot
        ENDIF
      ENDDO   
      CALL TRAGRZ(0.0,0.0,0.0,r,z,rx,rt,rp,zx,zt,zp)
      sh3d(1)=(r-rg0)/a0
      el3d(1)=el3d(2)
      tr3d(1)=0.0
!Interpolate onto input grid
!  shift_r        
      CALL W_LIN_INTERP(n,rm3d,sh3d,nr_r,rhot_r,shift_r,iflag,message)
!  elong_r
      CALL W_LIN_INTERP(n,rm3d,el3d,nr_r,rhot_r,elong_r,iflag,message)
!  triang_r
      CALL W_LIN_INTERP(n,rm3d,tr3d,nr_r,rhot_r,triang_r,iflag,message)
!  f_r
      CALL W_LIN_INTERP(n,rm3d,fpol,nr_r,rhot_r,f_r,iflag,message)
!  q_r
      CALL W_LIN_INTERP(n,rm3d,qj,nr_r,rhot_r,q_r,iflag,message)
!  psi_r
      CALL W_LIN_INTERP(n,rm3d,psij,nr_r,rhot_r,psi_r,iflag,message)
!  vp_r
      CALL W_LIN_INTERP(n,rm3d,atpg1(1,1),nr_r,rhot_r,vp_r,iflag,
     &                  message)
      CALL RARRAY_SCALE(nr_r,(2.0*z_pi)**2,vp_r,1)
!  grho2_r
      DO i=2,n
        zd3d(i)=atpg1(i,2)/atpg1(i,1)
      ENDDO
      zd3d(1)=zd3d(2)
      CALL W_LIN_INTERP(n,rm3d,zd3d,nr_r,rhot_r,grho2_r,iflag,message)
!  gth_r
      CALL W_LIN_INTERP(n,rm3d,atpg1(1,3),nr_r,rhot_r,gth_r,iflag,
     &                  message)
!  gph_r
      CALL W_LIN_INTERP(n,rm3d,atpg1(1,4),nr_r,rhot_r,gph_r,iflag,
     &                  message)
!  b2_r
      DO i=2,n
        zd3d(i)=atpg1(i,6)/atpg1(i,1)
      ENDDO
      zd3d(1)=zd3d(2)
      CALL W_LIN_INTERP(n,rm3d,zd3d,nr_r,rhot_r,b2_r,iflag,message)
!  grth_r
      DO i=2,n
        zd3d(i)=1.0/atpg1(i,7)
      ENDDO
      zd3d(1)=0.0
      CALL W_LIN_INTERP(n,rm3d,zd3d,nr_r,rhot_r,grth_r,iflag,message)
!  ftrap_r
!  yd3d-Lin-Liu and Miller lower limit
!  zd3d-Lin-Liu and Miller upper limit
!  Normalize to (r/R)**0.5 for interpolation - important near axis
      DO i=2,n
        con1=atpg1(i,5)/bmax(i)/atpg1(i,1)
        con2=atpg1(i,6)/bmax(i)**2/atpg1(i,1)
        yd3d(i)=1.0-con2*atpg1(i,10)/atpg1(i,1)
        yd3d(i)=yd3d(i)/SQRT(rm3d(i)/a0)
        zd3d(i)=1.0-(con2/con1**2)*(1.0-SQRT(1.0-con1)*(1.0+con1/2.0))
        zd3d(i)=zd3d(i)/SQRT(rm3d(i)/a0)
      ENDDO
      yd3d(1)=3.0*SQRT(2.0)/z_pi
      zd3d(1)=1.5
      CALL W_LIN_INTERP(n,rm3d,yd3d,nr_r,rhot_r,zdum1,iflag,message)
      CALL W_LIN_INTERP(n,rm3d,zd3d,nr_r,rhot_r,zdum2,iflag,message)
      DO i=1,nr_r
        ftrap_r(i)=(0.25*zdum1(i)+0.75*zdum2(i))*SQRT(rhot_r(i)/a0)
      ENDDO
!  fm_r(1)
      DO i=2,n
        zd3d(i)=2.0*atpg1(i,11)*atpg1(i,14)/atpg1(i,6)/atpg1(i,8)
      ENDDO
      zd3d(1)=zd3d(2)
      CALL W_LIN_INTERP(n,rm3d,zd3d,nr_r,rhot_r,zdum1,iflag,message)
      DO i=1,nr_r
        fm_r(1,i)=zdum1(i)
      ENDDO
!  fm_r(2)
      DO i=2,n
        zd3d(i)=2.0*atpg1(i,12)*atpg1(i,15)/atpg1(i,6)/atpg1(i,8)
      ENDDO
      zd3d(1)=zd3d(2)
      CALL W_LIN_INTERP(n,rm3d,zd3d,nr_r,rhot_r,zdum1,iflag,message)
      DO i=1,nr_r
        fm_r(2,i)=zdum1(i)
      ENDDO
!  fm_r(3)
      DO i=2,n
        zd3d(i)=2.0*atpg1(i,13)*atpg1(i,16)/atpg1(i,6)/atpg1(i,8)
      ENDDO
      zd3d(1)=zd3d(2)
      CALL W_LIN_INTERP(n,rm3d,zd3d,nr_r,rhot_r,zdum1,iflag,message)
      DO i=1,nr_r
        fm_r(3,i)=zdum1(i)
      ENDDO
!  bm2_r
      DO i=2,n
        zd3d(i)=atpg1(i,17)/atpg1(i,1)
      ENDDO
      zd3d(1)=zd3d(2)
      CALL W_LIN_INTERP(n,rm3d,zd3d,nr_r,rhot_r,bm2_r,iflag,message)
!  gr2bm2_r
      DO i=2,n
        zd3d(i)=atpg1(i,18)/atpg1(i,1)
      ENDDO
      zd3d(1)=zd3d(2)
      CALL W_LIN_INTERP(n,rm3d,zd3d,nr_r,rhot_r,gr2bm2_r,iflag,message)
!  dvol_r
      dvol_r(1)=vp_r(2)*rhot_r(2)*0.5
      DO i=2,nr_r-1
        dvol_r(i)=(vp_r(i)/rhot_r(i)+vp_r(i+1)/rhot_r(i+1))*0.5
     &         *(rhot_r(i+1)**2-rhot_r(i)**2)*0.5
      ENDDO   
 1000 RETURN
      END
