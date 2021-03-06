       SUBROUTINE USER_RHS_CHAN_PHYSICAL
       INCLUDE 'header'
! Here, you can add terms to the right hand side
! of the momentum and scalar equations.  
! The right hand side forcing arrays, CF2, CF2, CF3, CFTH
! are in Fourier space.  The velocity and scalars are available 
! in physical space.
! S1 is available as a working variable

       INTEGER I,J,K,N
#ifdef SOLAR
       REAL*8 CSolar,srflx,COEF,Dangle,Hangle
       REAL*8 e_sat,LatRad,cff1,cff2,cff,zenith,lonr,Tair,Hair
       REAL*8 yday,RSolar,vap_p
       REAL*8 acoef,bcoef,zeta1coef,zeta2coef
       REAL*8 analyticalflag,FRCTIME,FRCTIME2
       REAL*8 FRCSOLAR,FRCSOLAR2
! Type II Jerlov water, Paulson and Simpson Table 2
       acoef=0.77d0
!       zeta1coef=1.5d0
! SWITCH TO 5.0 here because of numerical issues with all the heat accumulating in the top grid cell
       zeta1coef=5.0d0
       bcoef=1.0d0-acoef
       zeta2coef=14.0d0
! IRENE: use analytical flag
       analyticalflag=0.0d0
! g*alpha/rho*c
! buoyancy units
!       COEF=9.81d0*2.6d-4/(1025d0*3850)
! temperature units
       COEF=1.0d0/(1022.8d0*4000.0d0)
       IF (analyticalflag.ge.1.0d0) THEN
! declination angle (radians) following roms
! yearday mid june, june 14 
       yday=165d0
       Dangle=23.44d0*COS((172.0d0-yday)*2.0d0*3.14159d0/365.25d0)
       Dangle=Dangle*2.0d0*3.14159d0/360.0d0
!
!  Compute hour angle (radians).
!
      Hangle=(12.0d0-MOD(TIME/3600.0d0,24d0))*3.14159d0/12.0d0
! Latitude
      LatRad=18.0d0*2.0d0*3.14159d0/360.0d0
      cff1=SIN(LatRad)*SIN(Dangle)
      cff2=COS(LatRad)*COS(Dangle)
      CSolar=1355.0d0
      RSolar=CSolar*COEF
      srflx=0.0d0
! -17.0d0 is the longitude
      lonr=-17d0
      Tair=24.0d0
! humidity fraction, relative to 1
      Hair=0.8d0 
! the above 2 are estimated from loading heatfluxes_dbw.mat
      srflx=0.0d0
      zenith=cff1+cff2*COS(Hangle-(lonr*0.017453d0)/15.0d0)
          IF (zenith.gt.0.0d0) THEN
            cff=(0.7859d0+0.03477d0*Tair)/                              &
     &          (1.0d0+0.00412d0*Tair)
            e_sat=10.0d0**cff    ! saturation vapor pressure (hPa=mbar)
            vap_p=e_sat*Hair ! water vapor pressure (hPa=mbar)
            srflx=0.85d0*RSolar*zenith*zenith/                          &
     &                 ((zenith+2.7d0)*vap_p*0.001d0+                   &
     &                  1.085d0*zenith+0.1d0)
          END IF
      ELSE
      open(201,file='ocean_frc_time.txt',
     &  form='formatted',status='unknown')
      open(202,file='ocean_frc_swrad.txt',
     &  form='formatted',status='unknown')
      READ(201,*),FRCTIME
      READ(202,*),FRCSOLAR
      FRCTIME2=FRCTIME
      FRCSOLAR2=FRCSOLAR
      DO WHILE (FRCTIME2.LT.TIME)
      READ(201,*),FRCTIME2
      READ(202,*),FRCSOLAR2
      IF (FRCTIME2.GE.TIME) THEN
      srflx=COEF*(1.0d0-(TIME-FRCTIME)/(FRCTIME2-FRCTIME))*FRCSOLAR+
     & COEF*(1.0d0-(FRCTIME2-TIME)/(FRCTIME2-FRCTIME))*FRCSOLAR2
      ELSE
      FRCTIME=FRCTIME2
      FRCSOLAR=FRCSOLAR2
      END IF
      END DO
      CLOSE(201)
      CLOSE(202)
      END IF
      IF (RANK.eq.0) THEN
      WRITE(*,*),'SRFLX',FRCTIME,FRCTIME2,FRCSOLAR,FRCSOLAR2,srflx
      END IF
!
#endif
!       REAL*8 ALPHA 
!
!
! SW solar radiation (3d buoyancy source and light distribution)
!
! biogeochemical reactions
#ifdef BIO
      IF (FLAVOR.EQ.'WLT_BIO') THEN
         CALL WLT_BIO_RHS
      END IF
#endif
!! dynamical forcing 
#ifdef SOLAR
! add diurnal cycle penetrating shortwave radiation to RHS of TH1
       DO J=JSTART_TH(1),JEND_TH(1)
         DO K=0,NZP-1
           DO I=0,NXM
           S1(I,K,J)=srflx*acoef/zeta1coef*EXP((GYF(J)-LY)/zeta1coef)+  &
     &               srflx*bcoef/zeta2coef*EXP((GYF(J)-LY)/zeta2coef)
           END DO
         END DO
       END DO
! Convert to Fourier space
       CALL FFT_XZ_TO_FOURIER(S1,CS1,0,NY+1)
       DO J=JSTART_TH(1),JEND_TH(1)
         DO K=0,TNKZ
           DO I=0,NXP-1
             CFTH(I,K,J,1)=CFTH(I,K,J,1)+CS1(I,K,J)
           END DO
         END DO
       END DO
#endif

! sponges

! For example, to add a linear damping term (e.g. -alpha*U) to the RHS:
!       alpha=-0.1d0
!       DO J=JSTART,JEND
!         DO K=0,NZP-1
!           DO I=0,NXM
!             S1(I,K,J)=-alpha*U1(I,K,J)
!           END DO
!         END DO
!       END DO
! Convert to Fourier space
!       CALL FFT_XZ_TO_FOURIER(S1,CS1,0,NY+1)
!       DO J=JSTART,JEND
!         DO K=0,TNKZ 
!           DO I=0,NXP-1
!             CF1(I,K,J)=CF1(I,K,J)+CS1(I,K,J)
!           END DO
!         END DO
!       END DO

! For U2 do this...
! Note that the only thing that changes are the bounds of the J index
!       DO J=2,NY
!         DO K=0,NZP-1
!           DO I=0,NXM
!             S1(I,K,J)=-alpha*U2(I,K,J)
!           END DO
!         END DO
!       END DO
! Convert to Fourier space
!       CALL FFT_XZ_TO_FOURIER(S1,CS1,0,NY+1)
!       DO J=2,NY
!         DO K=0,TNKZ
!           DO I=0,NXP-1
!             CF2(I,K,J)=CF2(I,K,J)+CS1(I,K,J)
!           END DO
!         END DO
!       END DO

! For scalars, do this...
! Loop over all scalars
!       DO N=1,N_TH
!       DO J=JSTART,JEND
!         DO K=0,NZP-1
!           DO I=0,NXM
!             S1(I,K,J)=-alpha*TH(I,K,J,N)
!           END DO
!         END DO
!       END DO
! Convert to Fourier space
!       CALL FFT_XZ_TO_FOURIER(S1,CS1,0,NY+1)
!       DO J=JSTART,JEND
!         DO K=0,TNKZ
!           DO I=0,NXP-1
!             CFTH(I,K,J,N)=CFTH(I,K,J,N)+CS1(I,K,J)
!           END DO
!         END DO
!       END DO
!      END DO 

      RETURN 
      END


       SUBROUTINE USER_RHS_CHAN_FOURIER
       include 'header'
! Here, you can add terms to the right hand side
! of the momentum and scalar equations.  
! The right hand side forcing arrays, CF2, CF2, CF3, CFTH
! are in Fourier space.  The velocity and scalars are available 
! in Fourier space.
! S1 is available as a working variable

       integer i,j,k,n

       real*8 alpha

! For example, to add a linear damping term (e.g. -alpha*U) to the RHS:
!       alpha=-0.1d0
!       DO J=JSTART,JEND
!         DO K=0,TNKZ 
!           DO I=0,NXP-1
!             CF1(I,K,J)=CF1(I,K,J)-alpha*CU1(I,K,J)
!           END DO
!         END DO
!       END DO

! For U2 do this...
! Note that the only thing that changes are the bounds of the J index
!       DO J=2,NY
!         DO K=0,TNKZ
!           DO I=0,NXP-1
!             CF2(I,K,J)=CF2(I,K,J)-alpha*CU2(I,K,J)
!           END DO
!         END DO
!       END DO

! Advection owing to thermal wind
      IF (I_RO.ne.0.d0) THEN
      DO N=1,N_TH
      
      do j=JSTART,JEND
      do k=0,TNKZ
      do i=0,NXP-1
        CF1(I,K,J)=CF1(I,K,J)
     &           -(DRHODX(N)*RI(N)*GYF(J)/I_RO-
     &           0.5d0*DRHODX(N)*RI(N)*LY/I_RO)
     &                            *CIKZ(K)*CU1(I,K,J)
     &           -(-1.d0*DRHODZ(N)*RI(N)*GYF(J)/I_RO
     &             +0.5d0*DRHODZ(N)*RI(N)*LY/I_RO)
     &                            *CIKX(I)*CU1(I,K,J)
     &                   -(-1.d0*DRHODZ(N)*RI(N)/I_RO)
     &                      *0.5d0*(CU2(I,K,J)+CU2(I,K,J+1))
        CF3(I,K,J)=CF3(I,K,J)
     &           -(DRHODX(N)*RI(N)*GYF(J)/I_RO-
     &             0.5d0*DRHODX(N)*RI(N)*LY/I_RO)
     &                           *CIKZ(K)*CU3(I,K,J)
     &           -(-1.d0*DRHODZ(N)*RI(N)*GYF(J)/I_RO+
     &              0.5d0*DRHODZ(N)*RI(N)*LY/I_RO)
     &                           *CIKX(I)*CU3(I,K,J)
     &                   -(DRHODX(N)*RI(N)/I_RO)
     &                      *0.5d0*(CU2(I,K,J)+CU2(I,K,J+1))

      end do
      end do
      end do

      do j=2,NY
      do k=0,TNKZ
      do i=0,NXP-1
        CF2(I,K,J)=CF2(I,K,J)
     &            -(DRHODX(N)*RI(N)*GY(J)/I_RO-
     &              0.5d0*DRHODX(N)*RI(N)*LY/I_RO)
     &                     *CIKZ(K)*CU2(I,K,J)
     &            -(-1.d0*DRHODZ(N)*RI(N)*GY(J)/I_RO+
     &               0.5d0*DRHODZ(N)*RI(N)*LY/I_RO)
     &                     *CIKX(I)*CU2(I,K,J)
      end do
      end do
      end do

! Do for each scalar

      DO J=JSTART_TH(N),JEND_TH(N)
        DO K=0,TNKZ
          DO I=0,NXP-1
! Advection by the thermal wind flow that depends on
! one scalar e.g. buoyancy only (TH(1))
! NEED TO CHANGE THIS IF MULTIPLE TRACERS MODIFY DENSITY
            CFTH(I,K,J,N)=CFTH(I,K,J,N)
     &     -(RI(1)/I_RO*DRHODX(1)*GYF(J)
     &       -0.5d0*RI(1)/I_RO*DRHODX(1)*LY)
     &                 *CIKZ(K)*CTH(I,K,J,N)
     &     -(RI(1)/I_RO*-1.d0*DRHODZ(1)*GYF(J)
     &       +RI(1)/I_RO*0.5d0*DRHODZ(1)*LY)
     &                 *CIKX(I)*CTH(I,K,J,N)
          END DO
        END DO
      END DO

      END DO
      end if


! For scalars, do this...
!       DO J=JSTART,JEND
!         DO K=0,TNKZ
!           DO I=0,NXP-1
!             CFTH(I,K,J,N)=CFTH(I,K,J,N)-alpha*CTH(I,K,J,N)
!           END DO
!         END DO
!       END DO
!      END DO 
! Add sponge layer forcing
#ifdef SOLAR
! IRENE commeted out:
!      CALL SPONGE_TH(1)
!      CALL SPONGE_VEL
#else
      CALL SPONGE_TH(1)
      CALL SPONGE_VEL
#endif
      RETURN 
      END


      SUBROUTINE USER_RHS_PER_PHYSICAL
      include 'header'
C Optionally, add forcing terms to the right hand side
C of the momentum and scalar evolution equations
C Here, CF1, CF2, CF3, and CFTH are the forcing arrays in Fourier space
C U1, U2, U3, and TH are the velocity and scalar variables in Physical space

      integer i,j,k,n

C For forced homogeneous turbulence:
C Add some forcing to the system to keep the Batchelor scale fixed
      EK=0.d0
      DO J=0,NYM
        DO K=0,NZM
          DO I=0,NXM
              EK=EK+U1(I,K,J)**2.d0+U2(I,K,J)**2.d0+U3(I,K,J)**2.d0
          END DO
        END DO
      END DO
C Note, that each cell has the same volume, so we can just average over all points
      EK=EK/dble(NX*NY*NZ)
! Scale EK by an amount to compensate for dissipation from 2/3 de-aliasing:
      EK=0.8d0*EK
      DO J=0,NYM
        DO K=0,NZM
          DO I=0,NXM
            S1(I,K,J)=(EPSILON_TARGET/EK)*U1(I,K,J)
          END DO
        END DO
      END DO
      CALL FFT_XZY_TO_FOURIER(S1,CS1)
      DO J=0,TNKY
        DO K=0,TNKZ
          DO I=0,NKX
            CF1(I,K,J)=CF1(I,K,J)+CS1(I,K,J)
          END DO
        END DO
      END DO
      DO J=0,NYM
        DO K=0,NZM
          DO I=0,NXM
            S1(I,K,J)=(EPSILON_TARGET/EK)*U2(I,K,J)
          END DO
        END DO
      END DO
      CALL FFT_XZY_TO_FOURIER(S1,CS1)
      DO J=0,TNKY
        DO K=0,TNKZ
          DO I=0,NKX
            CF2(I,K,J)=CF2(I,K,J)+CS1(I,K,J)
          END DO
        END DO
      END DO
      DO J=0,NYM
        DO K=0,NZM
          DO I=0,NXM
            S1(I,K,J)=(EPSILON_TARGET/EK)*U3(I,K,J)
          END DO
        END DO
      END DO
      CALL FFT_XZY_TO_FOURIER(S1,CS1)
      DO J=0,TNKY
        DO K=0,TNKZ
          DO I=0,NKX
            CF3(I,K,J)=CF3(I,K,J)+CS1(I,K,J)
          END DO
        END DO
      END DO

      RETURN
      END

      SUBROUTINE USER_RHS_DUCT_PHYSICAL
      include 'header'
      RETURN
      END

      SUBROUTINE USER_RHS_CAVITY_PHYSICAL
      include 'header'
      RETURN
      END 


      SUBROUTINE SPONGE_TH(N)
! This subroutine applies a sponge relaxation (Rayleigh damping) towards a
! specified background state for the temperature field
! The intention is to allow an open boundary
      include 'header'
      integer i,j,k,n
      real*8 L_sponge,L_bottom
      real*8 SPONGE_AMP

! The following variables will store the background state
      real*8 TH_0(-1:NY+1)

      real*8 RI_B(0:NY+1)

! This variable will hold the forcing rate
      real*8 SPONGE_SIGMA(0:NY+1)

! Set the amplitude of the sponge
      SPONGE_AMP=5.0d0*I_RO
! Set the top of the sponge layer in physical units
      L_sponge=18.d0
! Set the bottom of the computational domain in physical units
      L_bottom=0.d0
        DO J=0,NY+1
! Quadratic damping at lower wall
         if (GYF(J).lt.L_sponge) then
           SPONGE_SIGMA(j)=SPONGE_AMP*((L_sponge-GYF(J))
     &       /(L_sponge-L_bottom))**2.d0
         else
           SPONGE_SIGMA(j)=0.d0
         end if
        END DO

!     if (n.eq.1) then
      TH_0(0)=0.d0
      do j=1,NY+1
!        RI_B(J)=20.d0
             IF (GYF(J).LT.30.d0) THEN
! integral form
             TH_0(J)=TH_BC_YMIN_C1(1)*GYF(J)+
     & (TH_BC_YMIN_C2(1)-TH_BC_YMIN_C1(1))/2.0d0*GYF(J)+
     & 3.0d0*(TH_BC_YMIN_C2(1)-TH_BC_YMIN_C1(1))/2.0d0*
     & (LOG(COSH((GYF(J)-20.0d0)/3.0d0))-LOG(COSH(-20.0d0/3.0d0)))
!             TH_0(J)=TH_0(J-1)+(GY(J)-GY(J-1))*
!     & ((TH_BC_YMIN_C2(1)-TH_BC_YMIN_C1(1))
!     &   *(0.5d0*TANH((GYF(J)-20.0d0)/3.0d0)+0.5d0)
!     &   +TH_BC_YMIN_C1(1))
             ELSE
         TH_0(J)=0.0d0
! differential form fails for multiple vertical splits
!             TH_0(J)=TH_0(J-1)+(GY(J)-GY(J-1))*
!     & ((TH_BC_YMIN_C2(1)-TH_BC_YMIN_C3(1))
!     &   *(0.5d0*TANH((-GYF(J)+40.0d0)/3.0d0)+0.5d0)
!     &    +TH_BC_YMIN_C3(1))
             END IF
!     &                    RI_B(J)*(RI(N)*DRHODX(N))**2.d0
!     &                    /I_RO**2.d0/RI(N)
      end do
!      else
!        do j=0,NY+1
!          TH_0(j)=0.d0
!        end do
!      end if

! Add damping to R-K terms
! Damp the perturbations towards 0
      do k=0,TNKZ
         do i=0,NXP-1
           if ((RANKZ.ne.0).or.(i.ne.0).or.(k.ne.0)) then
           do j=JSTART_TH(N),JEND_TH(N)
              CFTH(i,k,j,n)=CFTH(i,k,j,n)
     &                 -SPONGE_SIGMA(j)*(CTH(i,k,j,n)-0.)
           end do
           end if
         end do
      end do
! Damp the mean gradient towards TH_0
      if (RANKZ.eq.0) then
        do j=JSTART_TH(N),JEND_TH(N)
          CFTH(0,0,j,n)=CFTH(0,0,j,n)-SPONGE_SIGMA(j)
     &          *(CTH(0,0,j,n)-TH_0(J))
        end do
      end if

      return
      end
  
      SUBROUTINE SPONGE_VEL
! This subroutine applies a sponge relaxation (Rayleigh damping) towards a
! specified background state
! The intention is to allow an open boundary
      include 'header'
      integer i,j,k

      real*8 L_sponge,L_bottom
      real*8 SPONGE_AMP

! The following variables will store the background state
      real*8 U1_0(-1:NY+1), U2_0(0:NY+1), U3_0(-1:NY+1)

! This variable will hold the forcing rate
      real*8 SPONGE_SIGMA(0:NY+1)

! Set the amplitude of the sponge
      SPONGE_AMP=5.0d0*I_RO
! Set the top of the sponge layer in physical units
      L_sponge=18.d0
! Set the bottom of the computational domain in physical units
      L_bottom=0.d0
        DO J=0,NY+1
! Quadratic damping at lower wall
         if (GYF(J).lt.L_sponge) then
           SPONGE_SIGMA(j)=SPONGE_AMP*((L_sponge-GYF(J))
     &       /(L_sponge-L_bottom))**2.d0
         else
           SPONGE_SIGMA(j)=0.d0
         end if
        END DO

! Set the background state
! Here, set the background to be geostrophic, with a linear temperature profile
      do j=0,NY+1
        U1_0(j)=0.d0
        U3_0(j)=0.d0
      end do
      do j=0,NY+1
        U2_0(j)=0.d0
      end do

! Add damping function to explicit R-K 
       do k=0,TNKZ
         do i=0,NXP-1
           if ((RANKZ.ne.0).or.(i.ne.0).or.(k.ne.0)) then
           do j=jstart,jend
             CF1(I,K,J)=CF1(I,K,J)-SPONGE_SIGMA(j)*(CU1(i,k,j)-0.d0)
             CF3(I,K,J)=CF3(I,K,J)-SPONGE_SIGMA(j)*(CU3(i,k,j)-0.d0)
           end do
           do j=2,NY
             CF2(I,K,J)=CF2(I,K,J)-
     &        0.5*(SPONGE_SIGMA(j)+SPONGE_SIGMA(j+1))*(CU2(i,k,j)-0.d0)
           end do
           end if   
         end do
      end do
! Damp mean flow
      if (RANKZ.eq.0) then
      do j=jstart,jend
        CF1(0,0,j)=CF1(0,0,j)-SPONGE_SIGMA(j)*(CU1(0,0,j)-U1_0(j))
        CF3(0,0,j)=CF3(0,0,j)-SPONGE_SIGMA(j)*(CU3(0,0,j)-U3_0(j))
      end do
      do j=2,NY
        CF2(0,0,j)=CF2(0,0,j)-SPONGE_SIGMA(j)*(CU2(0,0,j)-U2_0(j))
      end do
      end if

      return
      end






