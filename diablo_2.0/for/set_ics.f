C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE CREATE_FLOW_CHAN
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'


      INTEGER I,J,K,N
      REAL*8 RNUM1,RNUM2,RNUM3,XN,ZN,A1,WNUM,LNX,LNZ
      INTEGER, DIMENSION(:), ALLOCATABLE :: seed

C Initialize the random number generator
      CALL RANDOM_SEED(SIZE = K)
      Allocate (seed(1:K))
      do k=1,K
        seed(k)=RANK+k+999
      end do
      CALL RANDOM_SEED(PUT = seed)

C UBULK0 and KICK should be set in input.dat

C IC_TYPE is set in input_chan.dat and can be used to easily
C control which initial condition is used.  A few examples
C are given here. These can be modified, or new types can be
C added 

       IF (IC_TYPE.eq.0) then
C Parabolic profile for laminar closed channel flow
       DO J=0,NY
         DO K=0,NZP-1
           DO I=0,NXM
             U1(I,K,J)=(3./2.)*UBULK0*(1.d0-GYF(J)**2.)
             U2(I,K,J)=0.
             U3(I,K,J)=0.
           END DO
         END DO
      END DO
      else if (IC_TYPE.eq.1) then
C Laminar profile for open channel flow :
       DO K=0,NZP-1
         DO I=0,NXM
           DO J=1,NY
             U1(I,K,J)=-(3./2.)*UBULK0*GYF(J)**2.+3.*UBULK0*GYF(J)
             U2(I,K,J)=0.
             U3(I,K,J)=0.
           END DO
           U1(I,K,0)=0.
           U3(I,K,0)=0.
           U1(I,K,NY+1)=0.
           U3(I,K,NY+1)=0.
         END DO
      END DO
      else if (IC_TYPE.eq.2) then
C Linear profile for laminar Couette flow:
       DO J=0,NY
         DO K=0,NZP-1
           DO I=0,NXM
             U1(I,K,J)=gyf(j)
             U2(I,K,J)=0.
             U3(I,K,J)=0.
           END DO
         END DO
      END DO
      else if (IC_TYPE.eq.3) then
C Tanh shear layer
       DO J=0,NY
         DO K=0,NZP-1
           DO I=0,NXM
             U1(I,K,J)=TANH(GYF(J))
             U2(I,K,J)=0.d0
             U3(I,K,J)=0.d0
           END DO
         END DO
       END DO
      else if ((IC_TYPE.eq.4).or.(IC_TYPE.eq.5)) then
       DO J=0,NY
         DO K=0,NZP-1
           DO I=0,NXM
             U1(I,K,J)=0.d0
             U2(I,K,J)=0.d0
             U3(I,K,J)=0.d0
           END DO
         END DO
       END DO
      else if (IC_TYPE.eq.6) then
C EADY/STORM/FIXED RED NOISE 
       DO J=0,NY
         DO K=0,NZP-1
           DO I=0,NXM
!             IF (J.GE.(NY-37)) THEN
       XN=REAL(I,8)
       LNX=REAL(NX,8)
       ZN=REAL(RANK,8)*REAL(NZP,8)+REAL(K,8)
       LNZ=REAL(NZ,8)
!       A1=0.0010d0
       A1=0.0d0
       WNUM=1.0d0
       U1(I,K,J)=U1(I,K,J)+A1/WNUM*SIN(6.2832d0*(WNUM*XN/LNX+0.796d0))
       U1(I,K,J)=U1(I,K,J)+A1/WNUM*SIN(6.2832d0*(WNUM*ZN/LNZ+0.099d0))
       WNUM=2.0d0
       U1(I,K,J)=U1(I,K,J)+A1/WNUM*SIN(6.2832d0*(WNUM*XN/LNX+0.262d0))
       U1(I,K,J)=U1(I,K,J)+A1/WNUM*SIN(6.2832d0*(WNUM*ZN/LNZ+0.335d0))
       WNUM=3.0d0
       U1(I,K,J)=U1(I,K,J)+A1/WNUM*SIN(6.2832d0*(WNUM*XN/LNX+0.680d0))
       U1(I,K,J)=U1(I,K,J)+A1/WNUM*SIN(6.2832d0*(WNUM*ZN/LNZ+0.137d0))
       WNUM=4.0d0
       U1(I,K,J)=U1(I,K,J)+A1/WNUM*SIN(6.2832d0*(WNUM*XN/LNX+0.721d0))
       U1(I,K,J)=U1(I,K,J)+A1/WNUM*SIN(6.2832d0*(WNUM*ZN/LNZ+0.107d0))
       WNUM=5.0d0
       U1(I,K,J)=U1(I,K,J)+A1/WNUM*SIN(6.2832d0*(WNUM*XN/LNX+0.654d0))
       U1(I,K,J)=U1(I,K,J)+A1/WNUM*SIN(6.2832d0*(WNUM*ZN/LNZ+0.494d0))
       WNUM=6.0d0
       U1(I,K,J)=U1(I,K,J)+A1/WNUM*SIN(6.2832d0*(WNUM*XN/LNX+0.779d0))
       U1(I,K,J)=U1(I,K,J)+A1/WNUM*SIN(6.2832d0*(WNUM*ZN/LNZ+0.715d0))
       WNUM=7.0d0
       U1(I,K,J)=U1(I,K,J)+A1/WNUM*SIN(6.2832d0*(WNUM*XN/LNX+0.904d0))
       U1(I,K,J)=U1(I,K,J)+A1/WNUM*SIN(6.2832d0*(WNUM*ZN/LNZ+0.891d0))
       WNUM=8.0d0
       U1(I,K,J)=U1(I,K,J)+A1/WNUM*SIN(6.2832d0*(WNUM*XN/LNX+0.334d0))
       U1(I,K,J)=U1(I,K,J)+A1/WNUM*SIN(6.2832d0*(WNUM*ZN/LNZ+0.699d0))
       XN=REAL(I,8)
       LNX=REAL(NX,8)
       ZN=REAL(RANK,8)*REAL(NZP,8)+REAL(K,8)
       LNZ=REAL(NZ,8)
!
       WNUM=1.0d0
       U3(I,K,J)=U3(I,K,J)+A1/WNUM*SIN(6.2832d0*(WNUM*XN/LNX+0.796d0))
       U3(I,K,J)=U3(I,K,J)+A1/WNUM*SIN(6.2832d0*(WNUM*ZN/LNZ+0.099d0))
       WNUM=2.0d0
       U3(I,K,J)=U3(I,K,J)+A1/WNUM*SIN(6.2832d0*(WNUM*XN/LNX+0.262d0))
       U3(I,K,J)=U3(I,K,J)+A1/WNUM*SIN(6.2832d0*(WNUM*ZN/LNZ+0.335d0))
       WNUM=3.0d0
       U3(I,K,J)=U3(I,K,J)+A1/WNUM*SIN(6.2832d0*(WNUM*XN/LNX+0.680d0))
       U3(I,K,J)=U3(I,K,J)+A1/WNUM*SIN(6.2832d0*(WNUM*ZN/LNZ+0.137d0))
       WNUM=4.0d0
       U3(I,K,J)=U3(I,K,J)+A1/WNUM*SIN(6.2832d0*(WNUM*XN/LNX+0.721d0))
       U3(I,K,J)=U3(I,K,J)+A1/WNUM*SIN(6.2832d0*(WNUM*ZN/LNZ+0.107d0))
       WNUM=5.0d0
       U3(I,K,J)=U3(I,K,J)+A1/WNUM*SIN(6.2832d0*(WNUM*XN/LNX+0.654d0))
       U3(I,K,J)=U3(I,K,J)+A1/WNUM*SIN(6.2832d0*(WNUM*ZN/LNZ+0.494d0))
       WNUM=6.0d0
       U3(I,K,J)=U3(I,K,J)+A1/WNUM*SIN(6.2832d0*(WNUM*XN/LNX+0.779d0))
       U3(I,K,J)=U3(I,K,J)+A1/WNUM*SIN(6.2832d0*(WNUM*ZN/LNZ+0.715d0))
       WNUM=7.0d0
       U3(I,K,J)=U3(I,K,J)+A1/WNUM*SIN(6.2832d0*(WNUM*XN/LNX+0.904d0))
       U3(I,K,J)=U3(I,K,J)+A1/WNUM*SIN(6.2832d0*(WNUM*ZN/LNZ+0.891d0))
       WNUM=8.0d0
       U3(I,K,J)=U3(I,K,J)+A1/WNUM*SIN(6.2832d0*(WNUM*XN/LNX+0.334d0))
       U3(I,K,J)=U3(I,K,J)+A1/WNUM*SIN(6.2832d0*(WNUM*ZN/LNZ+0.699d0))
!
       U2(I,K,J)=0.d0
!       ELSE
!       U1(I,K,J)=0.d0
!       U2(I,K,J)=0.d0
!       U3(I,K,J)=0.d0
!       END IF
! multiply by tanh to remove strong shear at ML base
       U3(I,K,J)=U3(I,K,J)*
     &     REAL(0.5d0*TANH(REAL((GYF(J)-40.0d0)/10.0d0))+0.5d0)
       U1(I,K,J)=U1(I,K,J)*
     &     REAL(0.5d0*TANH(REAL((GYF(J)-40.0d0)/10.0d0))+0.5d0)
             END DO
         END DO
       END DO
       else if (IC_TYPE.eq.7) then
       DO J=0,NY
         DO K=0,NZP-1
           DO I=0,NXM
       U1(I,K,J)=0.d0
       U2(I,K,J)=0.d0
       U3(I,K,J)=0.d0
           END DO
         END DO
       END DO
      else
        WRITE(*,*) 'WARNING, unsupported IC_TYPE in CREATE_FLOW'
      end if
C Add random noise in physical space
      CALL RANDOM_NUMBER(RNUM1)
      CALL RANDOM_NUMBER(RNUM1)
      CALL RANDOM_NUMBER(RNUM1)

      DO J=0,NY+1
       DO K=0,NZP-1
         DO I=0,NXM
           CALL RANDOM_NUMBER(RNUM1)
           U1(I,K,J)=U1(I,K,J)+KICK*(RNUM1-0.5d0)
           CALL RANDOM_NUMBER(RNUM1)
           U2(I,K,J)=U2(I,K,J)+KICK*(RNUM1-0.5d0)
           CALL RANDOM_NUMBER(RNUM1)
           U3(I,K,J)=U3(I,K,J)+KICK*(RNUM1-0.5d0)
         END DO
        END DO
      END DO

C Zero the ghost cells
       IF (.NOT.USE_MPI) THEN
       DO K=0,NZM
         DO I=0,NXM
           U1(I,K,0)=0.
           U2(I,K,0)=0.
           U3(I,K,0)=0.
           U1(I,K,NY+1)=0.
           U2(I,K,NY+1)=0.
           U3(I,K,NY+1)=0.
         END DO
      END DO
      END IF

C Convert to Fourier space      
      CALL FFT_XZ_TO_FOURIER(U1,CU1,0,NY+1)
      CALL FFT_XZ_TO_FOURIER(U2,CU2,0,NY+1)
      CALL FFT_XZ_TO_FOURIER(U3,CU3,0,NY+1)
      CALL FFT_XZ_TO_FOURIER(P,CP,0,NY+1)

! Optionally, add random noise in Fourier space instead
!      DO I=0,NXP-1
!        DO J=1,NY
!          DO K=0,TNKZ
C Now, give the velocity field a random perturbation
!            CALL RANDOM_NUMBER(RNUM1)
!            CALL RANDOM_NUMBER(RNUM2)
!            CU1(I,K,J)=CU1(I,K,J)
!     &           +CMPLX((RNUM1-0.5d0),(RNUM2-0.5d0))*KICK
!            CALL RANDOM_NUMBER(RNUM1)
!            CALL RANDOM_NUMBER(RNUM2)
!            CU2(I,K,J)=CU2(I,K,J)
!     &           +CMPLX((RNUM1-0.5d0),(RNUM2-0.5d0))*KICK
!            CALL RANDOM_NUMBER(RNUM1)
!            CALL RANDOM_NUMBER(RNUM2)
!            CU3(I,K,J)=CU3(I,K,J)
!     &           +CMPLX((RNUM1-0.5d0),(RNUM2-0.5d0))*KICK
!          END DO
!          IF (TNKZ.EQ.0) THEN
! Here, In the 2d case we want to add a kick to the mean in z
!            K=0         
!            CALL RANDOM_NUMBER(RNUM1)
!            CALL RANDOM_NUMBER(RNUM2)
!            CALL RANDOM_NUMBER(RNUM3)

!            IF (IC_TYPE.eq.3) THEN
!              CU1(I,K,J)=CU1(I,K,J)
!     &             +(RNUM1-0.5)*KICK*EXP(-(GYF(J)*20.d0)**2.d0)
!              CU2(I,K,J)=CU2(I,K,J)
!     &             +(RNUM1-0.5)*KICK*EXP(-(GYF(J)*20.d0)**2.d0)
!              CU3(I,K,J)=CU3(I,K,J)
!     &             +(RNUM1-0.5)*KICK*EXP(-(GYF(J)*20.d0)**2.d0)
!            ELSE
!              CU1(I,K,J)=CU1(I,K,J)+(RNUM1-0.5)*KICK
!              CU2(I,K,J)=CU2(I,K,J)+(RNUM2-0.5)*KICK
!              CU3(I,K,J)=CU3(I,K,J)+(RNUM3-0.5)*KICK
!            END IF
!          END IF 
!        END DO
!      END DO

      IF (USE_MPI) THEN
        CALL GHOST_CHAN_MPI
      END IF
C Apply Boundary conditions to velocity field
      IF (USE_MPI) THEN
        CALL APPLY_BC_VEL_MPI
      ELSE
        CALL APPLY_BC_VEL_LOWER
        CALL APPLY_BC_VEL_UPPER
      END IF

C Remove the divergence of the velocity field
      CALL REM_DIV_CHAN

      IF (USE_MPI) THEN
        CALL GHOST_CHAN_MPI
      END IF

C Get the pressure from the poisson equation
!      CALL POISSON_P_CHAN
! Fix for the pressure
!      IF (USE_MPI) THEN
!        CALL GHOST_CHAN_MPI
!      END IF

C Save various statistics to keep track of the initial condition
      CALL SAVE_STATS_CHAN(.FALSE.)

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE CREATE_TH_CHAN
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
C Initialize the scalar fields
C In this subroutine, you should initialize each scalar field for the
C particular problem of interest

      INCLUDE 'header'
      INTEGER I,J,K,N

      DO N=1,N_TH
        IF (CREATE_NEW_TH(N)) THEN

      IF (IC_TYPE.eq.0) THEN
! As an example, initialize TH1 with a sine in x
       DO K=0,NZP-1
         DO I=0,NXM
           DO J=1,NY
             TH(I,K,J,N)=sin(2.d0*PI*GX(I)/LX)/(4.d0*PI**2.d0)
           END DO
         END DO
       END DO
       ELSE IF ((IC_TYPE.eq.1).or.(IC_TYPE.eq.2)) THEN
! Initialize with a linear profile using the bcs
       DO K=0,NZP-1
         DO I=0,NXM
           IF ((TH_BC_YMIN(N).EQ.0).AND.(TH_BC_YMAX(N).EQ.0)) THEN
               DO J=1,NY
               IF (GYF(J).LE.2.0) THEN
                 TH(I,K,J,N)=(TH_BC_YMAX_C1(N)-TH_BC_YMIN_C1(N))
     &                *(GYF(J)+1.)/2.0+TH_BC_YMIN_C1(N)
               ELSE
                 TH(I,K,J,N)=TH_BC_YMAX_C1(N)
               END IF
             END DO
           ELSE IF ((TH_BC_YMIN(N).EQ.1)
     &            .AND.(TH_BC_YMAX(N).EQ.1)) THEN
             DO J=1,NY
! Linear profile with slope corresponding to upper value
                TH(I,K,J,N)=TH_BC_YMIN_C1(N)*GYF(J)
              END DO
           ELSE
             IF (RANK.EQ.0) then
                WRITE(*,*) 'WARNING, THETA INITIALIZED TO ZERO ...'
                WRITE(*,*) 'CREATE AN INITIAL VALUE IN CREATE_FLOW_CHAN'
             end if
           END IF
         END DO
        END DO
        ELSE IF (IC_TYPE.eq.3) then
! Tanh profile 
       DO K=0,NZP-1
         DO I=0,NXM
           DO J=1,NY
             TH(I,K,J,N)=TANH(GYF(J))
           END DO
         END DO
       END DO
       ELSE IF (IC_TYPE.eq.4) THEN
       DO K=0,NZP-1
         DO I=0,NXM
           DO J=1,NY
             IF ((J.GE.(NY-9)).AND.(RANKY.eq.NPROCY-1)) THEN
             TH(I,K,J,N)=TH_BC_YMIN_C1(N)*GYF(NY-9)
             ELSE
             TH(I,K,J,N)=TH_BC_YMIN_C1(N)*GYF(J)
             END IF
           END DO
         END DO
       END DO
       ELSE IF (IC_TYPE.eq.5) THEN
       DO K=0,NZP-1
         DO I=0,NXM
           DO J=1,NY
             IF (N.EQ.1) THEN
             IF ((J.GE.(NY-9)).AND.(RANKY.eq.NPROCY-1)) THEN
             TH(I,K,J,N)=TH_BC_YMIN_C1(N)*GYF(NY-9)
             ELSE
             TH(I,K,J,N)=TH_BC_YMIN_C1(N)*GYF(J)
             END IF
             ELSE IF (N.EQ.2) THEN
             TH(I,K,J,N)=2.0d0
             ELSE IF ((N.EQ.3).OR.(N.EQ.4).OR.(N.EQ.5)) THEN
             TH(I,K,J,N)=2.0d0
             ELSE IF (N.EQ.6) THEN
             IF (GYF(J).GE.(LY/1.5d0)) THEN
             TH(I,K,J,N)=0.0d0
             ELSE
             TH(I,K,J,N)=1.0d0
             END IF
             END IF
           END DO
         END DO
       END DO
       ELSE IF (IC_TYPE.eq.6) THEN
       DO K=0,NZP-1
         DO I=0,NXM
           DO J=1,NY
             IF (N.EQ.1) THEN
! TANH N2
             IF (GYF(J).LT.30.d0) THEN
              TH(I,K,J,N)=TH_BC_YMIN_C1(N)*GYF(J)/2.0d0-
     & 1.5d0*TH_BC_YMIN_C1(N)*LOG(COSH((GYF(J)-20.0d0)/3.0d0))+
     & TH_BC_YMIN_C2(N)*GYF(J)/2.0d0 +
     & 1.5d0*TH_BC_YMIN_C2(N)*LOG(COSH((GYF(J)-20.0d0)/3.0d0))+
     & 1.5d0*TH_BC_YMIN_C1(N)*LOG(COSH((-20.0d0)/3.0d0))-
     & 1.5d0*TH_BC_YMIN_C2(N)*LOG(COSH((-20.0d0)/3.0d0))
!             TH(I,K,J,N)=TH(I,K,J-1,N)+(GY(J)-GY(J-1))* 
!     & ((TH_BC_YMIN_C2(N)-TH_BC_YMIN_C1(N))              
!     &   *(0.5d0*TANH((GYF(J)-20.0d0)/3.0d0)+0.5d0)
!     &   +TH_BC_YMIN_C1(N))
             ELSE 
              TH(I,K,J,N)=1.5d0*(TH_BC_YMIN_C3(N)-TH_BC_YMIN_C2(N))*
     & LOG(COSH((-GYF(J)+40.0d0)/3.0d0))+0.5d0*(TH_BC_YMIN_C2(N)+
     & TH_BC_YMIN_C3(N))*GYF(J)-
     & 1.5d0*(TH_BC_YMIN_C3(N)-TH_BC_YMIN_C2(N))*
     & LOG(COSH((10.0d0)/3.0d0))-0.5d0*(TH_BC_YMIN_C2(N)+
     & TH_BC_YMIN_C3(N))*30.0d0+
     & 15.0d0*(TH_BC_YMIN_C1(N)+TH_BC_YMIN_C2(N))+
     & 1.5d0*(TH_BC_YMIN_C1(N)-TH_BC_YMIN_C2(N))*
     & LOG(COSH((20.0d0)/3.0d0)/COSH(10.0d0/3.0d0))
!             TH(I,K,J,N)=TH(I,K,J-1,N)+(GY(J)-GY(J-1))* 
!     & ((TH_BC_YMIN_C2(N)-TH_BC_YMIN_C3(N))  
!     &   *(0.5d0*TANH((-GYF(J)+40.0d0)/3.0d0)+0.5d0)
!     &    +TH_BC_YMIN_C3(N))
             END IF
!STEP N2
!             IF (J.LT.(NY-55)) THEN
!             TH(I,K,J,N)=TH_BC_YMIN_C1(N)*GYF(J)
!             ELSE IF ((J.LT.(NY-37)).AND.(J.GE.(NY-55))) THEN
!             TH(I,K,J,N)=TH_BC_YMIN_C1(N)*GYF(NY-56)+
!     &                   TH_BC_YMIN_C2(N)*(GYF(J)-GYF(NY-56))
!             ELSE
!             TH(I,K,J,N)=TH_BC_YMIN_C1(N)*GYF(NY-56)+
!     &                   TH_BC_YMIN_C2(N)*(GYF(NY-38)-GYF(NY-56))
!     &                   +TH_BC_YMIN_C3(N)*(GYF(J)-GYF(NY-38))
!             END IF
! TANH TRACER (bad)
             ELSE IF (N.EQ.2) THEN
             IF ((J.LT.2).and.(RANKY.eq.0)) THEN
             TH(I,K,J,N)=2.0d0
             ELSE IF
     &    (((GYF(J)).LT.(30.0d0)).AND.
     &      (((J.GT.1).AND.(RANKY.eq.0)).or.(RANKY.gt.0))) THEN
             TH(I,K,J,N)=TH(I,K,J-1,N)+(GY(J)-GY(J-1))*
     & ((TH_BC_YMIN_C2(N)-TH_BC_YMIN_C1(N))
     &   *(0.5d0*TANH((GYF(J)-20.0d0)/3.0d0)+0.5d0)
     &   +TH_BC_YMIN_C1(N))
             ELSE
             TH(I,K,J,N)=TH(I,K,J-1,N)+(GY(J)-GY(J-1))*
     & ((TH_BC_YMIN_C2(N)-TH_BC_YMIN_C3(N))
     &   *(0.5d0*TANH((-GYF(J)+40.0d0)/3.0d0)+0.5d0)
     &    +TH_BC_YMIN_C3(N))
             END IF
!! END TANH TRACER
!
!             IF (J.LT.(NY-55)) THEN
!             TH(I,K,J,N)=2.0d0
!             ELSE IF ((J.LT.(NY-37)).AND.(J.GE.(NY-55))) THEN
!             TH(I,K,J,N)=2.0d0-
!     &        1.9d0*(GYF(J)-GYF(NY-56))/(GYF(NY-38)-GYF(NY-56))
!             ELSE 
!             TH(I,K,J,N)=0.1d0
!             END IF
             ELSE
             TH(I,K,J,N)=0.5d0
             END IF
           END DO
         END DO
       END DO
! MAURITANIA
       ELSE IF (IC_TYPE.eq.7) THEN
       DO K=0,NZP-1
         DO I=0,NXM
           DO J=1,NY
!             IF (N.EQ.1) THEN
!! stretching type 4 coeff: 1.0, NY=101
!             IF ((J.LT.(NY-12)).AND.(J.GE.(NY-59)).AND.
!     &           (RANKY.eq.NPROCY-1)) THEN
!             TH(I,K,J,N)=TH_BC_YMIN_C1(N)*GYF(NY-59)+TH_BC_YMIN_C2(N)*
!     &                   (GYF(J)-GYF(NY-59))
!             ELSEIF ((J.GE.(NY-12)).AND.(RANKY.eq.NPROCY-1)) THEN
!             TH(I,K,J,N)=TH(I,K,NY-13,N)
!             ELSE
!             TH(I,K,J,N)=TH_BC_YMIN_C1(N)*GYF(J)
!             END IF
!             ELSE IF (N.EQ.2) THEN
!             TH(I,K,J,N)=2.0d0
!             ELSE IF ((N.EQ.3).OR.(N.EQ.4).OR.(N.EQ.5)) THEN
!             TH(I,K,J,N)=2.0d0
!             ELSE IF (N.EQ.6) THEN
!             IF (GYF(J).GE.(LY/1.5d0)) THEN
!             TH(I,K,J,N)=0.0d0
!             ELSE
             TH(I,K,J,N)=1.0d0
!             END IF
!             END IF
           END DO
         END DO
       END DO
       ELSE
        WRITE(*,*) 'WARNING, unsupported IC_TYPE in CREATE_FLOW'
       END IF


      S1(:,:,:)=TH(:,:,:,N)
      CALL FFT_XZ_TO_FOURIER(S1,CS1,0,NY+1)
      CTH(:,:,:,N)=CS1(:,:,:)

      END IF
      END DO

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE CREATE_FLOW_PER
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      INTEGER I, J, K
      REAL*8 RNUM1,RNUM2,RNUM3
      REAL*8 K0

C For an initial vortex, define the location of the centerline
      REAL*8 XC(0:NY+1),ZC(0:NY+1)

      WRITE(6,*) 'Creating new flow from scratch.'

C Initialize random number generator
      CALL RANDOM_SEED

      IF (IC_TYPE.eq.0) THEN
C Initizlize the flow using a Taylor-Green vortex
C Nondimensionalize with U0 and 1/kappa
      DO J=0,NYM
        DO K=0,NZM
          DO I=0,NXM
            U1(I,K,J)=cos(2*pi*(GY(J))/LY)
     &               *cos(2*pi*(GX(I))/LX)
     &               *SIN(2*pi*(GZ(K))/LZ)
            U2(I,K,J)=0.d0
            U3(I,K,J)=-cos(2*pi*(GY(J))/LY)
     &               *sin(2*pi*(GX(I))/LX)
     &               *COS(2*pi*(GZ(K))/LZ)
          END DO
        END DO
      END DO
      ELSE IF (IC_TYPE.eq.1) THEN
C Start with an ideal vortex centered in the domain
      DO J=0,NYM
        XC(J)=LX/2.
        ZC(J)=LZ/2.
        DO K=0,NZM
          DO I=0,NXM
            IF ((GX(I)-XC(j))**2.+(GZ(K)-ZC(j))**2..gt.0.1) then
! If we aren't too close to the vortex center
              U1(I,K,J)=-1.d0*(GZ(K)-ZC(j))
     &                /((GX(I)-XC(j))**2.+(GZ(K)-ZC(j))**2.)
              U3(I,K,J)=1.d0*(GX(I)-XC(j))
     &                /((GX(I)-XC(j))**2.+(GZ(K)-ZC(j))**2.)
              U2(I,K,J)=0.d0
            ELSE
! Otherwise:
              U1(I,K,J)=-1.d0*(GZ(K)-ZC(j))
     &                /0.1
              U3(I,K,J)=1.d0*(GX(I)-XC(j))
     &                /0.1
              U2(I,K,J)=0.d0
            END IF
          END DO
        END DO
      END DO
      END IF
! Add random noise in Fourier space

      CALL FFT_XZY_TO_FOURIER(U1,CU1)
      CALL FFT_XZY_TO_FOURIER(U2,CU2)
      CALL FFT_XZY_TO_FOURIER(U3,CU3)
      DO J=1,TNKY
        DO K=1,TNKZ
          DO I=1,NKX
            CALL RANDOM_NUMBER(RNUM1)
            CALL RANDOM_NUMBER(RNUM2)
            CALL RANDOM_NUMBER(RNUM3)
            K0=sqrt(KX(I)**2.d0+KY(J)**2.d0+KZ(K)**2.d0)
     &        /sqrt(KX(1)**2.d0+KY(1)**2.d0+KZ(1)**2.d0)
            CU1(I,K,J)=CU1(I,K,J)+(RNUM1-0.5)*KICK/K0
            CU2(I,K,J)=CU2(I,K,J)+(RNUM1-0.5)*KICK/K0
            CU3(I,K,J)=CU3(I,K,J)+(RNUM1-0.5)*KICK/K0
          end do
        end do
      end do
! get the initial energy in low wavenumbers
      CALL FFT_XZY_TO_PHYSICAL(CU1,U1)
      CALL FFT_XZY_TO_PHYSICAL(CU2,U2)
      CALL FFT_XZY_TO_PHYSICAL(CU3,U3)
      EK0=0.d0
      DO J=0,NYM
        DO K=0,NZM
          DO I=0,NXM
              EK0=EK0+U1(I,K,J)**2.d0+U2(I,K,J)**2.d0+U3(I,K,J)**2.d0
          END DO
        END DO
      END DO
      write(*,*) 'EK0: ',EK0
      IF (N_TH.gt.0) THEN
!      EPSILON_TARGET=((1.d0/DX(1))**4.d0)*(NU**3.d0)*(PR(1))**(-2.d0)
      EPSILON_TARGET=((1.d0/DX(1))**4.d0)*(NU**3.d0)*(100.d0)**(-2.d0)
      write(*,*) 'EPSILON_TARGET: ',EPSILON_TARGET
      write(*,*) 'TARGET KOLMOGOROV SCALE: ',
     &         (NU**3.d0/epsilon_target)**(0.25d0)
      END IF
      CALL FFT_XZY_TO_FOURIER(U1,CU1)
      CALL FFT_XZY_TO_FOURIER(U2,CU2)
      CALL FFT_XZY_TO_FOURIER(U3,CU3)


      CALL REM_DIV_PER
      CALL POISSON_P_PER

      CALL SAVE_STATS_PER(.FALSE.)

      RETURN
      END
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE CREATE_TH_PER
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      INCLUDE 'header'
      INTEGER I,J,K,N

C Note, Since stratification is not permitted in the periodic flow field
C Any background stratification must be added to the governing equations

      DO N=1,N_TH
      IF (CREATE_NEW_TH(N)) THEN
        DO J=0,NYM
          DO K=0,NZM
            DO I=0,NXM
C Example: Gaussian patch centered in the domain
         TH(I,K,J,N)=EXP(-((GX(I)-LX/2)*10.d0)**2.d0
     &                   -((GY(J)-LY/2)*10.d0)**2.d0
     &                   -((GZ(K)-LZ/2)*10.d0)**2.d0)
            END DO
          END DO
        END DO
       CALL FFT_XZY_TO_FOURIER(TH(0,0,0,N),CTH(0,0,0,N))

       END IF
       END DO

       RETURN
       END

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE CREATE_FLOW_DUCT
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      RETURN
      END


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE CREATE_FLOW_CAV
C----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      RETURN
      END





