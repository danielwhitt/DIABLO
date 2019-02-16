!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE SET_BCS 
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! Here we set time-dependent boundary conditions analytically
      INCLUDE 'header'
          IF (NUM_PER_DIR.EQ.3) THEN
          ELSEIF (NUM_PER_DIR.EQ.2) THEN
          IF (ANA_FRC) THEN
            CALL ANA_BCS_CHAN
          ELSE
            CALL READ_BCS_CHAN
          END IF
          ELSEIF (NUM_PER_DIR.EQ.1) THEN
          ELSEIF (NUM_PER_DIR.EQ.0) THEN
          END IF
      RETURN
      END


!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE ANA_BCS_CHAN
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! Here we set time-dependent boundary conditions analytically
      INCLUDE 'header'
      INTEGER N
      REAL*8 TGT_UBC,TGT_WBC,AST,BST,CST,DST
      REAL*8 THBCPEAK
! if bc_type does not equal one then BC's are not updated and remain
! set to initial constants
! if bc_type is equal to one then we impose an oscillatory surface stress
! with frequency given by W_BC_YMAX_C2 and U_BC_YMAX_C2 and amplitude
! IRENE:
!      U_BC_YMIN_VAL=-0.006d0*SQRT(CU1(0,0,2)**2.0d0 +CU3(0,0,2)**2.0d0)
!     &              *CU1(0,0,2)/NU
!      W_BC_YMIN_VAL=-0.006d0*SQRT(CU1(0,0,2)**2.0d0 +CU3(0,0,2)**2.0d0)
!     &              *CU3(0,0,2)/NU
! DAN's default
      U_BC_YMIN_VAL=U_BC_YMIN_C1
      W_BC_YMIN_VAL=W_BC_YMIN_C1
! W_BC_YMAX_C1 and U_BC_YMAX_C1
! BUT WE RAMP UP THE STRESS OVER THE FIRST 10h  of simulation time
! SPECIFY ANALYTIC FUNCTION FOR THE SURFACE MOMENTUM FLUX
      IF (BC_TYPE.EQ.1) THEN ! 3x (Fix THBC should not be reset in 2x)
      AST = 259200.0d0 ! storm start time
      BST = 432000.0d0 ! storm mid time
      CST = 518400.0d0 ! storm end time
      DST = 414.0d0 ! storm stress target (0.4 -> 390 tau/rho) 
      THBCPEAK = -0.06d0
      IF (TIME.LT.AST) THEN
      U_BC_YMAX_VAL=0.0d0
      ELSE IF ((TIME.GE.AST).AND.(TIME.LT.BST)) THEN
         TGT_UBC = (TIME-AST)/(BST-AST)
     &              *DST
         U_BC_YMAX_VAL=TGT_UBC
!     & 0.995d0*(U_BC_YMAX_VAL-TGT_UBC)
!     & +(1.0d0-0.995d0**2.0d0)*r4_normal_01(0)
      ELSE IF ((TIME.GE.BST).AND.(TIME.LT.CST)) THEN
         TGT_UBC = DST*(1.0d0-(TIME-BST)/             
     &     (CST-BST))
         U_BC_YMAX_VAL=TGT_UBC
!     & 0.995d0*(U_BC_YMAX_VAL-TGT_UBC)
!     & +(1.0d0-0.995d0**2.0d0)*r4_normal_01(0)
      ELSE
       U_BC_YMAX_VAL=0.0d0
      END IF
      IF (TIME.LT.AST) THEN
      W_BC_YMAX_VAL=0.0d0
      ELSE IF ((TIME.GE.AST).AND.(TIME.LT.BST)) THEN
         IF (DRHODX(1).GT.0.0d0) THEN
         TGT_WBC = (TIME-AST)/(BST-AST)
     &              *DST
         ELSE
         TGT_WBC = -(TIME-AST)/(BST-AST)
     &              *DST
         END IF
         W_BC_YMAX_VAL=TGT_WBC
!     & 0.995d0*(W_BC_YMAX_VAL-TGT_WBC)
!     & +(1.0d0-0.995d0**2.0d0)*r4_normal_01(0)
      ELSE IF ((TIME.GE.BST).AND.(TIME.LT.CST)) THEN
! replace minus sign 
         IF (DRHODX(1).GT.0.0d0) THEN
         TGT_WBC = DST*(1.0d0-(TIME-BST)/                  
     &     (CST-BST))
         ELSE 
         TGT_WBC = -DST*(1.0d0-(TIME-BST)/                  
     &     (CST-BST))
         END IF
         W_BC_YMAX_VAL=TGT_WBC
!     & 0.995d0*(W_BC_YMAX_VAL-TGT_WBC)
!     & +(1.0d0-0.995d0**2.0d0)*r4_normal_01(0)
      ELSE
       W_BC_YMAX_VAL=0.0d0
      END IF
      DO N=1,N_TH
         IF (N.EQ.1) THEN
          IF ((TIME.GE.AST).AND.(TIME.LT.BST)) THEN
          TH_BC_YMAX_VAL(N)= (TIME-AST)/(BST-AST)
     &              *THBCPEAK
          ELSEIF ((TIME.GE.BST).AND.(TIME.LT.CST)) THEN
          TH_BC_YMAX_VAL(N)=THBCPEAK*(1.0d0-(TIME-BST)/
     &     (CST-BST))
          ELSE
          TH_BC_YMAX_VAL(N)=TH_BC_YMAX_C1(N)
          END IF
         ELSE
         TH_BC_YMAX_VAL(N)=TH_BC_YMAX_C1(N)
         END IF
      END DO
      ELSEIF (BC_TYPE.EQ.8) THEN
! buoyancy flux equivalent to EBF of W1024_*2x
      AST = 259200.0d0 ! storm start time
      BST = 432000.0d0 ! storm mid time
      CST = 518400.0d0 ! storm end time
      THBCPEAK = -0.207d0
      DO N=1,N_TH
         IF (N.EQ.1) THEN
          IF ((TIME.GE.AST).AND.(TIME.LT.BST)) THEN
          TH_BC_YMAX_VAL(N)= (TIME-AST)/(BST-AST)
     &              *THBCPEAK+TH_BC_YMAX_C1(N)
          ELSEIF ((TIME.GE.BST).AND.(TIME.LT.CST)) THEN
          TH_BC_YMAX_VAL(N)=THBCPEAK*(1.0d0-(TIME-BST)/
     &     (CST-BST))+TH_BC_YMAX_C1(N)
          ELSE
          TH_BC_YMAX_VAL(N)=TH_BC_YMAX_C1(N)
          END IF
         ELSE
         TH_BC_YMAX_VAL(N)=TH_BC_YMAX_C1(N)
         END IF
      END DO
      W_BC_YMAX_VAL=W_BC_YMAX_C1
      U_BC_YMAX_VAL=U_BC_YMAX_C1
      ELSEIF (BC_TYPE.EQ.9) THEN
! buoyancy flux equivalent to HALF of EBF of W1024_*2x
      AST = 259200.0d0 ! storm start time
      BST = 432000.0d0 ! storm mid time
      CST = 518400.0d0 ! storm end time
      THBCPEAK = -0.1035d0
      DO N=1,N_TH
         IF (N.EQ.1) THEN
          IF ((TIME.GE.AST).AND.(TIME.LT.BST)) THEN
          TH_BC_YMAX_VAL(N)= (TIME-AST)/(BST-AST)
     &              *THBCPEAK+TH_BC_YMAX_C1(N)
          ELSEIF ((TIME.GE.BST).AND.(TIME.LT.CST)) THEN
          TH_BC_YMAX_VAL(N)=THBCPEAK*(1.0d0-(TIME-BST)/
     &     (CST-BST))+TH_BC_YMAX_C1(N)
          ELSE
          TH_BC_YMAX_VAL(N)=TH_BC_YMAX_C1(N)
          END IF
         ELSE
         TH_BC_YMAX_VAL(N)=TH_BC_YMAX_C1(N)
         END IF
      END DO
      W_BC_YMAX_VAL=W_BC_YMAX_C1
      U_BC_YMAX_VAL=U_BC_YMAX_C1
      ELSEIF (BC_TYPE.EQ.10) THEN
! buoyancy flux equivalent to EBF of W1024_*2x, NW256 reset time
      AST = 2500.0d0 ! storm start time
      BST = 175300.0d0 ! storm mid time
      CST = 261700.0d0 ! storm end time
      THBCPEAK = -0.207d0
      DO N=1,N_TH
         IF (N.EQ.1) THEN
          IF ((TIME.GE.AST).AND.(TIME.LT.BST)) THEN
          TH_BC_YMAX_VAL(N)= (TIME-AST)/(BST-AST)
     &              *THBCPEAK+TH_BC_YMAX_C1(N)
          ELSEIF ((TIME.GE.BST).AND.(TIME.LT.CST)) THEN
          TH_BC_YMAX_VAL(N)=THBCPEAK*(1.0d0-(TIME-BST)/
     &     (CST-BST))+TH_BC_YMAX_C1(N)
          ELSE
          TH_BC_YMAX_VAL(N)=TH_BC_YMAX_C1(N)
          END IF
         ELSE
         TH_BC_YMAX_VAL(N)=TH_BC_YMAX_C1(N)
         END IF
      END DO
      W_BC_YMAX_VAL=W_BC_YMAX_C1
      U_BC_YMAX_VAL=U_BC_YMAX_C1
      ELSEIF (BC_TYPE.EQ.7) THEN
! No buoyancy flux start at day 0 
      AST = 2500.0d0 ! storm start time
      BST = 175300.0d0 ! storm mid time
      CST = 261700.0d0 ! storm end time
      DST = 414.0d0 ! storm stress target (0.4 -> 390 tau/rho) 
      IF (TIME.LT.AST) THEN
      U_BC_YMAX_VAL=0.0d0
      ELSE IF ((TIME.GE.AST).AND.(TIME.LT.BST)) THEN
         TGT_UBC = (TIME-AST)/(BST-AST)
     &              *DST
         U_BC_YMAX_VAL=TGT_UBC
      ELSE IF ((TIME.GE.BST).AND.(TIME.LT.CST)) THEN
         TGT_UBC = DST*(1.0d0-(TIME-BST)/             
     &     (CST-BST))
         U_BC_YMAX_VAL=TGT_UBC
      ELSE
       U_BC_YMAX_VAL=0.0d0
      END IF
      IF (TIME.LT.AST) THEN
      W_BC_YMAX_VAL=0.0d0
      ELSE IF ((TIME.GE.AST).AND.(TIME.LT.BST)) THEN
         IF (DRHODX(1).GT.0.0d0) THEN
         TGT_WBC = (TIME-AST)/(BST-AST)
     &              *DST
         ELSE
         TGT_WBC = -(TIME-AST)/(BST-AST)
     &              *DST
         END IF
         W_BC_YMAX_VAL=TGT_WBC
      ELSE IF ((TIME.GE.BST).AND.(TIME.LT.CST)) THEN
! replace minus sign 
         IF (DRHODX(1).GT.0.0d0) THEN
         TGT_WBC = DST*(1.0d0-(TIME-BST)/                  
     &     (CST-BST))
         ELSE 
         TGT_WBC = -DST*(1.0d0-(TIME-BST)/                  
     &     (CST-BST))
         END IF
         W_BC_YMAX_VAL=TGT_WBC
      ELSE
       W_BC_YMAX_VAL=0.0d0
      END IF
      DO N=1,N_TH
         TH_BC_YMAX_VAL(N)=TH_BC_YMAX_C1(N)
      END DO
      ELSEIF (BC_TYPE.EQ.11) THEN
! No buoyancy flux start at day 0 - small stress 1/8 of default WT17 
      AST = 2500.0d0 ! storm start time
      BST = 175300.0d0 ! storm mid time
      CST = 261700.0d0 ! storm end time
      DST = 51.75d0 ! storm stress target (0.4 -> 390 tau/rho) 
      IF (TIME.LT.AST) THEN
      U_BC_YMAX_VAL=0.0d0
      ELSE IF ((TIME.GE.AST).AND.(TIME.LT.BST)) THEN
         TGT_UBC = (TIME-AST)/(BST-AST)
     &              *DST
         U_BC_YMAX_VAL=TGT_UBC
      ELSE IF ((TIME.GE.BST).AND.(TIME.LT.CST)) THEN
         TGT_UBC = DST*(1.0d0-(TIME-BST)/
     &     (CST-BST))
         U_BC_YMAX_VAL=TGT_UBC
      ELSE
       U_BC_YMAX_VAL=0.0d0
      END IF
      IF (TIME.LT.AST) THEN
      W_BC_YMAX_VAL=0.0d0
      ELSE IF ((TIME.GE.AST).AND.(TIME.LT.BST)) THEN
         IF (DRHODX(1).GT.0.0d0) THEN
         TGT_WBC = (TIME-AST)/(BST-AST)
     &              *DST
         ELSE
         TGT_WBC = -(TIME-AST)/(BST-AST)
     &              *DST
         END IF
         W_BC_YMAX_VAL=TGT_WBC
      ELSE IF ((TIME.GE.BST).AND.(TIME.LT.CST)) THEN
! replace minus sign 
         IF (DRHODX(1).GT.0.0d0) THEN
         TGT_WBC = DST*(1.0d0-(TIME-BST)/
     &     (CST-BST))
         ELSE
         TGT_WBC = -DST*(1.0d0-(TIME-BST)/
     &     (CST-BST))
         END IF
         W_BC_YMAX_VAL=TGT_WBC
      ELSE
       W_BC_YMAX_VAL=0.0d0
      END IF
      DO N=1,N_TH
         TH_BC_YMAX_VAL(N)=TH_BC_YMAX_C1(N)
      END DO
      ELSEIF (BC_TYPE.EQ.12) THEN
! No buoyancy flux start at day 0 -  same as WT17 as published 
      AST = 2500.0d0 ! storm start time
      BST = 175300.0d0 ! storm mid time
      CST = 261700.0d0 ! storm end time
      DST = 207.0d0 ! storm stress target (0.4 -> 390 tau/rho) 
      IF (TIME.LT.AST) THEN
      U_BC_YMAX_VAL=0.0d0
      ELSE IF ((TIME.GE.AST).AND.(TIME.LT.BST)) THEN
         TGT_UBC = (TIME-AST)/(BST-AST)
     &              *DST
         U_BC_YMAX_VAL=TGT_UBC
      ELSE IF ((TIME.GE.BST).AND.(TIME.LT.CST)) THEN
         TGT_UBC = DST*(1.0d0-(TIME-BST)/
     &     (CST-BST))
         U_BC_YMAX_VAL=TGT_UBC
      ELSE
       U_BC_YMAX_VAL=0.0d0
      END IF
      IF (TIME.LT.AST) THEN
      W_BC_YMAX_VAL=0.0d0
      ELSE IF ((TIME.GE.AST).AND.(TIME.LT.BST)) THEN
         IF (DRHODX(1).GT.0.0d0) THEN
         TGT_WBC = (TIME-AST)/(BST-AST)
     &              *DST
         ELSE
         TGT_WBC = -(TIME-AST)/(BST-AST)
     &              *DST
         END IF
         W_BC_YMAX_VAL=TGT_WBC
      ELSE IF ((TIME.GE.BST).AND.(TIME.LT.CST)) THEN
! replace minus sign 
         IF (DRHODX(1).GT.0.0d0) THEN
         TGT_WBC = DST*(1.0d0-(TIME-BST)/
     &     (CST-BST))
         ELSE
         TGT_WBC = -DST*(1.0d0-(TIME-BST)/
     &     (CST-BST))
         END IF
         W_BC_YMAX_VAL=TGT_WBC
      ELSE
       W_BC_YMAX_VAL=0.0d0
      END IF
      DO N=1,N_TH
         TH_BC_YMAX_VAL(N)=TH_BC_YMAX_C1(N)
      END DO
      ELSEIF (BC_TYPE.EQ.13) THEN
! buoyancy flux + wind WT17 
      AST = 2500.0d0 ! storm start time
      BST = 175300.0d0 ! storm mid time
      CST = 261700.0d0 ! storm end time
      DST = 207.0d0 ! storm stress target (0.4 -> 390 tau/rho) 
      IF (TIME.LT.AST) THEN
      U_BC_YMAX_VAL=0.0d0
      ELSE IF ((TIME.GE.AST).AND.(TIME.LT.BST)) THEN
         TGT_UBC = (TIME-AST)/(BST-AST)
     &              *DST
         U_BC_YMAX_VAL=TGT_UBC
      ELSE IF ((TIME.GE.BST).AND.(TIME.LT.CST)) THEN
         TGT_UBC = DST*(1.0d0-(TIME-BST)/
     &     (CST-BST))
         U_BC_YMAX_VAL=TGT_UBC
      ELSE
       U_BC_YMAX_VAL=0.0d0
      END IF
      IF (TIME.LT.AST) THEN
      W_BC_YMAX_VAL=0.0d0
      ELSE IF ((TIME.GE.AST).AND.(TIME.LT.BST)) THEN
         IF (DRHODX(1).GT.0.0d0) THEN
         TGT_WBC = (TIME-AST)/(BST-AST)
     &              *DST
         ELSE
         TGT_WBC = -(TIME-AST)/(BST-AST)
     &              *DST
         END IF
         W_BC_YMAX_VAL=TGT_WBC
      ELSE IF ((TIME.GE.BST).AND.(TIME.LT.CST)) THEN
! replace minus sign 
         IF (DRHODX(1).GT.0.0d0) THEN
         TGT_WBC = DST*(1.0d0-(TIME-BST)/
     &     (CST-BST))
         ELSE
         TGT_WBC = -DST*(1.0d0-(TIME-BST)/
     &     (CST-BST))
         END IF
         W_BC_YMAX_VAL=TGT_WBC
      ELSE
       W_BC_YMAX_VAL=0.0d0
      END IF
      THBCPEAK = -0.1035d0
      DO N=1,N_TH
         IF (N.EQ.1) THEN
          IF ((TIME.GE.AST).AND.(TIME.LT.BST)) THEN
          TH_BC_YMAX_VAL(N)= (TIME-AST)/(BST-AST)
     &              *THBCPEAK+TH_BC_YMAX_C1(N)
          ELSEIF ((TIME.GE.BST).AND.(TIME.LT.CST)) THEN
          TH_BC_YMAX_VAL(N)=THBCPEAK*(1.0d0-(TIME-BST)/
     &     (CST-BST))+TH_BC_YMAX_C1(N)
          ELSE
          TH_BC_YMAX_VAL(N)=TH_BC_YMAX_C1(N)
          END IF
         ELSE
         TH_BC_YMAX_VAL(N)=TH_BC_YMAX_C1(N)
         END IF
      END DO
      ELSEIF (BC_TYPE.EQ.6) THEN
! No buoyancy flux test
      AST = 259200.0d0 ! storm start time
      BST = 432000.0d0 ! storm mid time
      CST = 518400.0d0 ! storm end time
      DST = 414.0d0 ! storm stress target (0.4 -> 390 tau/rho) 
      IF (TIME.LT.AST) THEN
      U_BC_YMAX_VAL=0.0d0
      ELSE IF ((TIME.GE.AST).AND.(TIME.LT.BST)) THEN
         TGT_UBC = (TIME-AST)/(BST-AST)
     &              *DST
         U_BC_YMAX_VAL=TGT_UBC
      ELSE IF ((TIME.GE.BST).AND.(TIME.LT.CST)) THEN
         TGT_UBC = DST*(1.0d0-(TIME-BST)/             
     &     (CST-BST))
         U_BC_YMAX_VAL=TGT_UBC
      ELSE
       U_BC_YMAX_VAL=0.0d0
      END IF
      IF (TIME.LT.AST) THEN
      W_BC_YMAX_VAL=0.0d0
      ELSE IF ((TIME.GE.AST).AND.(TIME.LT.BST)) THEN
         IF (DRHODX(1).GT.0.0d0) THEN
         TGT_WBC = (TIME-AST)/(BST-AST)
     &              *DST
         ELSE
         TGT_WBC = -(TIME-AST)/(BST-AST)
     &              *DST
         END IF
         W_BC_YMAX_VAL=TGT_WBC
      ELSE IF ((TIME.GE.BST).AND.(TIME.LT.CST)) THEN
! replace minus sign 
         IF (DRHODX(1).GT.0.0d0) THEN
         TGT_WBC = DST*(1.0d0-(TIME-BST)/                  
     &     (CST-BST))
         ELSE 
         TGT_WBC = -DST*(1.0d0-(TIME-BST)/                  
     &     (CST-BST))
         END IF
         W_BC_YMAX_VAL=TGT_WBC
      ELSE
       W_BC_YMAX_VAL=0.0d0
      END IF
      DO N=1,N_TH
         TH_BC_YMAX_VAL(N)=TH_BC_YMAX_C1(N)
      END DO
      ELSEIF (BC_TYPE.EQ.5) THEN
      AST = 259200.0d0 ! storm start time
      BST = 432000.0d0 ! storm mid time
      CST = 518400.0d0 ! storm end time
      DST = 506.9d0 ! storm stress target (0.4 -> 390 tau/rho) 
      THBCPEAK = -0.06d0
      IF (TIME.LT.AST) THEN
      U_BC_YMAX_VAL=0.0d0
      ELSE IF ((TIME.GE.AST).AND.(TIME.LT.BST)) THEN
         TGT_UBC = (TIME-AST)/(BST-AST)
     &              *DST
         U_BC_YMAX_VAL=TGT_UBC/1.732d0
      ELSE IF ((TIME.GE.BST).AND.(TIME.LT.CST)) THEN
         TGT_UBC = DST*(1.0d0-(TIME-BST)/             
     &     (CST-BST))
         U_BC_YMAX_VAL=TGT_UBC/1.732d0
      ELSE
       U_BC_YMAX_VAL=0.0d0
      END IF
      IF (TIME.LT.AST) THEN
      W_BC_YMAX_VAL=0.0d0
      ELSE IF ((TIME.GE.AST).AND.(TIME.LT.BST)) THEN
         IF (DRHODX(1).GT.0.0d0) THEN
         TGT_WBC = (TIME-AST)/(BST-AST)
     &              *DST
         ELSE
         TGT_WBC = -(TIME-AST)/(BST-AST)
     &              *DST
         END IF
         W_BC_YMAX_VAL=TGT_WBC
      ELSE IF ((TIME.GE.BST).AND.(TIME.LT.CST)) THEN
! replace minus sign 
         IF (DRHODX(1).GT.0.0d0) THEN
         TGT_WBC = DST*(1.0d0-(TIME-BST)/                  
     &     (CST-BST))
         ELSE 
         TGT_WBC = -DST*(1.0d0-(TIME-BST)/                  
     &     (CST-BST))
         END IF
         W_BC_YMAX_VAL=TGT_WBC
      ELSE
       W_BC_YMAX_VAL=0.0d0
      END IF
      DO N=1,N_TH
         IF (N.EQ.1) THEN
          IF ((TIME.GE.AST).AND.(TIME.LT.BST)) THEN
          TH_BC_YMAX_VAL(N)= (TIME-AST)/(BST-AST)
     &              *THBCPEAK
          ELSEIF ((TIME.GE.BST).AND.(TIME.LT.CST)) THEN
          TH_BC_YMAX_VAL(N)=THBCPEAK*(1.0d0-(TIME-BST)/
     &     (CST-BST))
          ELSE
          TH_BC_YMAX_VAL(N)=TH_BC_YMAX_C1(N)
          END IF
         ELSE
         TH_BC_YMAX_VAL(N)=TH_BC_YMAX_C1(N)
         END IF
      END DO
      ELSEIF (BC_TYPE.EQ.3) THEN 
! WIND ANGLE TEST 45 degree
      AST = 259200.0d0 ! storm start time
      BST = 432000.0d0 ! storm mid time
      CST = 518400.0d0 ! storm end time
      DST = 414.0d0 ! storm stress target (0.4 -> 390 tau/rho) 
      THBCPEAK = -0.06d0
      IF (TIME.LT.AST) THEN
      U_BC_YMAX_VAL=0.0d0
      ELSE IF ((TIME.GE.AST).AND.(TIME.LT.BST)) THEN
         TGT_UBC = (TIME-AST)/(BST-AST)
     &              *DST
         U_BC_YMAX_VAL=TGT_UBC
!     & 0.995d0*(U_BC_YMAX_VAL-TGT_UBC)
!     & +(1.0d0-0.995d0**2.0d0)*r4_normal_01(0)
      ELSE IF ((TIME.GE.BST).AND.(TIME.LT.CST)) THEN
         TGT_UBC = DST*(1.0d0-(TIME-BST)/
     &     (CST-BST))
         U_BC_YMAX_VAL=TGT_UBC
!     & 0.995d0*(U_BC_YMAX_VAL-TGT_UBC)
!     & +(1.0d0-0.995d0**2.0d0)*r4_normal_01(0)
      ELSE
       U_BC_YMAX_VAL=0.0d0
      END IF
      IF (TIME.LT.AST) THEN
      W_BC_YMAX_VAL=0.0d0
      ELSE IF ((TIME.GE.AST).AND.(TIME.LT.BST)) THEN
         IF (DRHODX(1).GT.0.0d0) THEN
         TGT_WBC = -(TIME-AST)/(BST-AST)
     &              *DST
         ELSE
         TGT_WBC = (TIME-AST)/(BST-AST)
     &              *DST
         END IF
         W_BC_YMAX_VAL=TGT_WBC
!     & 0.995d0*(W_BC_YMAX_VAL-TGT_WBC)
!     & +(1.0d0-0.995d0**2.0d0)*r4_normal_01(0)
      ELSE IF ((TIME.GE.BST).AND.(TIME.LT.CST)) THEN
         IF (DRHODX(1).GT.0.0d0) THEN
         TGT_WBC = -DST*(1.0d0-(TIME-BST)/
     &     (CST-BST))
         ELSE
         TGT_WBC = DST*(1.0d0-(TIME-BST)/
     &     (CST-BST))
         END IF
         W_BC_YMAX_VAL=TGT_WBC
!     & 0.995d0*(W_BC_YMAX_VAL-TGT_WBC)
!     & +(1.0d0-0.995d0**2.0d0)*r4_normal_01(0)
      ELSE
       W_BC_YMAX_VAL=0.0d0
      END IF
      DO N=1,N_TH
         IF (N.EQ.1) THEN
          IF ((TIME.GE.AST).AND.(TIME.LT.BST)) THEN
          TH_BC_YMAX_VAL(N)= (TIME-AST)/(BST-AST)
     &              *THBCPEAK
          ELSEIF ((TIME.GE.BST).AND.(TIME.LT.CST)) THEN
          TH_BC_YMAX_VAL(N)=THBCPEAK*(1.0d0-(TIME-BST)/
     &     (CST-BST))
          ELSE
          TH_BC_YMAX_VAL(N)=TH_BC_YMAX_C1(N)
          END IF
         ELSE
         TH_BC_YMAX_VAL(N)=TH_BC_YMAX_C1(N)
         END IF
      END DO
      ELSEIF (BC_TYPE.EQ.4) THEN
! WIND ANGLE TEST 0 degree
      AST = 259200.0d0 ! storm start time
      BST = 432000.0d0 ! storm mid time
      CST = 518400.0d0 ! storm end time
      DST = 585.40d0 ! storm stress target (0.4 -> 390 tau/rho) 
      THBCPEAK = -0.06d0
      IF (TIME.LT.AST) THEN
      U_BC_YMAX_VAL=0.0d0
      ELSE IF ((TIME.GE.AST).AND.(TIME.LT.BST)) THEN
         TGT_UBC = (TIME-AST)/(BST-AST)
     &              *DST
         U_BC_YMAX_VAL=TGT_UBC
      ELSE IF ((TIME.GE.BST).AND.(TIME.LT.CST)) THEN
         TGT_UBC = DST*(1.0d0-(TIME-BST)/
     &     (CST-BST))
         U_BC_YMAX_VAL=TGT_UBC
      ELSE
       U_BC_YMAX_VAL=0.0d0
      END IF
      W_BC_YMAX_VAL=0.0d0
      DO N=1,N_TH
         IF (N.EQ.1) THEN
          IF ((TIME.GE.AST).AND.(TIME.LT.BST)) THEN
          TH_BC_YMAX_VAL(N)= (TIME-AST)/(BST-AST)
     &              *THBCPEAK
          ELSEIF ((TIME.GE.BST).AND.(TIME.LT.CST)) THEN
          TH_BC_YMAX_VAL(N)=THBCPEAK*(1.0d0-(TIME-BST)/
     &     (CST-BST))
          ELSE
          TH_BC_YMAX_VAL(N)=TH_BC_YMAX_C1(N)
          END IF
         ELSE
         TH_BC_YMAX_VAL(N)=TH_BC_YMAX_C1(N)
         END IF
      END DO
      ELSE IF (BC_TYPE.EQ.2) THEN
! MAURITANIA
      IF (TIME.LT.U_BC_YMAX_C3) THEN
         U_BC_YMAX_VAL=TIME/U_BC_YMAX_C3*U_BC_YMAX_C1
      ELSE
!         U_BC_YMAX_VAL=COS((TIME-U_BC_YMAX_C3)
!         &
!     &                 *U_BC_YMAX_C2)*U_BC_YMAX_C1
          U_BC_YMAX_VAL=U_BC_YMAX_C1
      END IF
      IF (TIME.LT.W_BC_YMAX_C3) THEN
         W_BC_YMAX_VAL=TIME/W_BC_YMAX_C3*W_BC_YMAX_C1
      ELSE
!         W_BC_YMAX_VAL=COS((TIME-W_BC_YMAX_C3)
!     &   *W_BC_YMAX_C2)*W_BC_YMAX_C1
          W_BC_YMAX_VAL=W_BC_YMAX_C1
      END IF
      DO N=1,N_TH
         IF (N.EQ.1) THEN
       TH_BC_YMAX_VAL(N)=COS((TIME-TH_BC_YMAX_C3(N))*TH_BC_YMAX_C2(N))
     &                     *TH_BC_YMAX_C1(N)
         ELSE
         TH_BC_YMAX_VAL(N)=TH_BC_YMAX_C1(N)
         END IF
      END DO
      ELSEIF (BC_TYPE.EQ.21) THEN
! REVISED MAURITANIA January 24, 2018 
      AST = 2500.0d0 ! wind start time
      BST = 175300.0d0 ! wind mid time
      DST = -130.0d0 ! storm mid time
      IF (TIME.LT.AST) THEN
      W_BC_YMAX_VAL=0.0d0
      ELSE IF ((TIME.GE.AST).AND.(TIME.LT.BST)) THEN
         TGT_WBC = (TIME-AST)/(BST-AST)
     &              *DST
         W_BC_YMAX_VAL=TGT_WBC
      ELSE IF ((TIME.GE.BST)) THEN
         TGT_WBC = DST
         W_BC_YMAX_VAL=TGT_WBC
      ELSE
       W_BC_YMAX_VAL=0.0d0
      END IF
      TH_BC_YMAX_VAL(1)=TH_BC_YMAX_C1(1)
      U_BC_YMAX_VAL=U_BC_YMAX_C1
      ELSEIF (BC_TYPE.EQ.14) THEN ! WT17 CS reverse wind 
      AST = 259200.0d0 ! storm start time
      BST = 432000.0d0 ! storm mid time
      CST = 518400.0d0 ! storm end time
      DST = 414.0d0 ! storm stress target (0.4 -> 390 tau/rho) 
      IF (TIME.LT.AST) THEN
      U_BC_YMAX_VAL=0.0d0
      ELSE IF ((TIME.GE.AST).AND.(TIME.LT.BST)) THEN
         TGT_UBC = (TIME-AST)/(BST-AST)
     &              *DST
         U_BC_YMAX_VAL=-TGT_UBC ! minus sign for reverse wind
      ELSE IF ((TIME.GE.BST).AND.(TIME.LT.CST)) THEN
         TGT_UBC = DST*(1.0d0-(TIME-BST)/             
     &     (CST-BST))
         U_BC_YMAX_VAL=-TGT_UBC ! minus sign for reverse wind
      ELSE
       U_BC_YMAX_VAL=0.0d0
      END IF
      IF (TIME.LT.AST) THEN
      W_BC_YMAX_VAL=0.0d0
      ELSE IF ((TIME.GE.AST).AND.(TIME.LT.BST)) THEN
         TGT_WBC = -(TIME-AST)/(BST-AST)
     &              *DST
         W_BC_YMAX_VAL=-TGT_WBC ! minus sign, reversed wind
      ELSE IF ((TIME.GE.BST).AND.(TIME.LT.CST)) THEN
         TGT_WBC = -DST*(1.0d0-(TIME-BST)/                  
     &     (CST-BST))
         W_BC_YMAX_VAL=-TGT_WBC ! minus sign, reverse wind
      ELSE
       W_BC_YMAX_VAL=0.0d0
      END IF
      DO N=1,N_TH
         TH_BC_YMAX_VAL(N)=TH_BC_YMAX_C1(N)
      END DO
      ELSEIF (BC_TYPE.EQ.15) THEN
! No buoyancy flux start at day 0 
      AST = 2500.0d0 ! storm start time
      BST = 175300.0d0 ! storm mid time
      CST = 261700.0d0 ! storm end time
      DST = 414.0d0 ! storm stress target (0.4 -> 390 tau/rho) 
      THBCPEAK = -0.06d0
      IF (TIME.LT.AST) THEN
      U_BC_YMAX_VAL=0.0d0
      ELSE IF ((TIME.GE.AST).AND.(TIME.LT.BST)) THEN
         TGT_UBC = (TIME-AST)/(BST-AST)
     &              *DST
         U_BC_YMAX_VAL=TGT_UBC
      ELSE IF ((TIME.GE.BST).AND.(TIME.LT.CST)) THEN
         TGT_UBC = DST*(1.0d0-(TIME-BST)/
     &     (CST-BST))
         U_BC_YMAX_VAL=TGT_UBC
      ELSE
       U_BC_YMAX_VAL=0.0d0
      END IF
      IF (TIME.LT.AST) THEN
      W_BC_YMAX_VAL=0.0d0
      ELSE IF ((TIME.GE.AST).AND.(TIME.LT.BST)) THEN
         IF (DRHODX(1).GT.0.0d0) THEN
         TGT_WBC = (TIME-AST)/(BST-AST)
     &              *DST
         ELSE
         TGT_WBC = -(TIME-AST)/(BST-AST)
     &              *DST
         END IF
         W_BC_YMAX_VAL=TGT_WBC
      ELSE IF ((TIME.GE.BST).AND.(TIME.LT.CST)) THEN
! replace minus sign 
         IF (DRHODX(1).GT.0.0d0) THEN
         TGT_WBC = DST*(1.0d0-(TIME-BST)/
     &     (CST-BST))
         ELSE
         TGT_WBC = -DST*(1.0d0-(TIME-BST)/
     &     (CST-BST))
         END IF
         W_BC_YMAX_VAL=TGT_WBC
      ELSE
       W_BC_YMAX_VAL=0.0d0
      END IF
      DO N=1,N_TH
         IF (N.EQ.1) THEN
          IF ((TIME.GE.AST).AND.(TIME.LT.BST)) THEN
          TH_BC_YMAX_VAL(N)= (TIME-AST)/(BST-AST)
     &              *THBCPEAK
          ELSEIF ((TIME.GE.BST).AND.(TIME.LT.CST)) THEN
          TH_BC_YMAX_VAL(N)=THBCPEAK*(1.0d0-(TIME-BST)/
     &     (CST-BST))
          ELSE
          TH_BC_YMAX_VAL(N)=TH_BC_YMAX_C1(N)
          END IF
         ELSE
         TH_BC_YMAX_VAL(N)=TH_BC_YMAX_C1(N)
         END IF
      END DO
      ELSEIF (BC_TYPE.EQ.16) THEN ! WT17 CS +45 degrees 
      AST = 259200.0d0 ! storm start time
      BST = 432000.0d0 ! storm mid time
      CST = 518400.0d0 ! storm end time
      DST = 414.0d0 ! storm stress target (0.4 -> 390 tau/rho) 
      IF (TIME.LT.AST) THEN
      U_BC_YMAX_VAL=0.0d0
      ELSE IF ((TIME.GE.AST).AND.(TIME.LT.BST)) THEN
         TGT_UBC = (TIME-AST)/(BST-AST)
     &              *DST
         U_BC_YMAX_VAL=TGT_UBC
      ELSE IF ((TIME.GE.BST).AND.(TIME.LT.CST)) THEN
         TGT_UBC = DST*(1.0d0-(TIME-BST)/
     &     (CST-BST))
         U_BC_YMAX_VAL=TGT_UBC
      ELSE
       U_BC_YMAX_VAL=0.0d0
      END IF
      IF (TIME.LT.AST) THEN
      W_BC_YMAX_VAL=0.0d0
      ELSE IF ((TIME.GE.AST).AND.(TIME.LT.BST)) THEN
         TGT_WBC = -(TIME-AST)/(BST-AST)
     &              *DST
         W_BC_YMAX_VAL=-TGT_WBC
      ELSE IF ((TIME.GE.BST).AND.(TIME.LT.CST)) THEN
         TGT_WBC = -DST*(1.0d0-(TIME-BST)/
     &     (CST-BST))
         W_BC_YMAX_VAL=-TGT_WBC
      ELSE
       W_BC_YMAX_VAL=0.0d0
      END IF
      DO N=1,N_TH
         TH_BC_YMAX_VAL(N)=TH_BC_YMAX_C1(N)
      END DO
      ELSEIF (BC_TYPE.EQ.17) THEN ! WT17 CS +225 degrees 
      AST = 259200.0d0 ! storm start time
      BST = 432000.0d0 ! storm mid time
      CST = 518400.0d0 ! storm end time
      DST = 414.0d0 ! storm stress target (0.4 -> 390 tau/rho) 
      IF (TIME.LT.AST) THEN
      U_BC_YMAX_VAL=0.0d0
      ELSE IF ((TIME.GE.AST).AND.(TIME.LT.BST)) THEN
         TGT_UBC = (TIME-AST)/(BST-AST)
     &              *DST
         U_BC_YMAX_VAL=-TGT_UBC
      ELSE IF ((TIME.GE.BST).AND.(TIME.LT.CST)) THEN
         TGT_UBC = DST*(1.0d0-(TIME-BST)/
     &     (CST-BST))
         U_BC_YMAX_VAL=-TGT_UBC
      ELSE
       U_BC_YMAX_VAL=0.0d0
      END IF
      IF (TIME.LT.AST) THEN
      W_BC_YMAX_VAL=0.0d0
      ELSE IF ((TIME.GE.AST).AND.(TIME.LT.BST)) THEN
         TGT_WBC = -(TIME-AST)/(BST-AST)
     &              *DST
         W_BC_YMAX_VAL=TGT_WBC
      ELSE IF ((TIME.GE.BST).AND.(TIME.LT.CST)) THEN
         TGT_WBC = -DST*(1.0d0-(TIME-BST)/
     &     (CST-BST))
         W_BC_YMAX_VAL=TGT_WBC
      ELSE
       W_BC_YMAX_VAL=0.0d0
      END IF
      DO N=1,N_TH
         TH_BC_YMAX_VAL(N)=TH_BC_YMAX_C1(N)
      END DO
      ELSE
!ANOTHER BC 
         U_BC_YMAX_VAL=U_BC_YMAX_C1
         W_BC_YMAX_VAL=W_BC_YMAX_C1
      DO N=1,N_TH
         TH_BC_YMAX_VAL(N)=TH_BC_YMAX_C1(N)
      END DO
      END IF

      RETURN
      END
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE READ_BCS_CHAN
!----*|--.---------.---------.---------.---------.---------.---------.-|-------|
! Here we set time-dependent boundary conditions by reading from a 
! forcing file
      INCLUDE 'header'
      INTEGER I,N
      REAL*8 FRCTIME,FRCU,FRCW,FRCTH1,FRCTH2
      REAL*8 FRCTIME2,FRCU2,FRCW2,FRCTH12,FRCTH22
! IRENE:
      IF ((RANKZ.EQ.0).AND.(RANKY.EQ.0)) THEN
! coefficient = (karman/log(DY/roughness length))^2
      U_BC_YMIN_VAL=((0.41d0/LOG(DY(2)/0.002d0))**2.0d0)
     &                             *SQRT(dble(CU1(0,0,2))**2.0d0 
     &                            +dble(CU3(0,0,2))**2.0d0)
     &              *dble(CU1(0,0,2))/NU
      W_BC_YMIN_VAL=((0.41d0/LOG(DY(2)/0.002d0))**2.0d0)
     &                       *SQRT(dble(CU1(0,0,2))**2.0d0 
     &                            +dble(CU3(0,0,2))**2.0d0)
     &              *dble(CU3(0,0,2))/NU
         CALL MPI_BCAST(U_BC_YMIN_VAL,1,MPI_DOUBLE_PRECISION,0,
     &        MPI_COMM_WORLD,IERROR)
         CALL MPI_BCAST(W_BC_YMIN_VAL,1,MPI_DOUBLE_PRECISION,0,
     &        MPI_COMM_WORLD,IERROR)
      ELSE
         CALL MPI_BCAST(U_BC_YMIN_VAL,1,MPI_DOUBLE_PRECISION,0,
     &        MPI_COMM_WORLD,IERROR)
         CALL MPI_BCAST(W_BC_YMIN_VAL,1,MPI_DOUBLE_PRECISION,0,
     &        MPI_COMM_WORLD,IERROR)
      END IF
      IF ((RANKY.EQ.0).AND.(RANKZ.EQ.0)) THEN
      WRITE(*,*) "U_BC_YMIN_VAL",RANKY,RANKZ,CU1(0,0,2),U_BC_YMIN_VAL
      WRITE(*,*) "W_BC_YMIN_VAL",RANKY,RANKZ,CU3(0,0,2),W_BC_YMIN_VAL
      END IF
! DAN's default
!      U_BC_YMIN_VAL=U_BC_YMIN_C1
!      W_BC_YMIN_VAL=W_BC_YMIN_C1
      open(101,file='ocean_frc_time.txt',
     &  form='formatted',status='unknown')
      open(102,file='ocean_frc_sustr.txt', 
     &  form='formatted',status='unknown')
      open(103,file='ocean_frc_svstr.txt',
     &  form='formatted',status='unknown')
      open(104,file='ocean_frc_shfluxminusswrad.txt',
     &  form='formatted',status='unknown')
      open(105,file='ocean_frc_EminusP.txt',
     &  form='formatted',status='unknown')
      READ(101,*),FRCTIME
      READ(102,*),FRCU
      READ(103,*),FRCW
      READ(104,*),FRCTH1
      READ(105,*),FRCTH2
      FRCTIME=FRCTIME2
      FRCU=FRCU2
      FRCW=FRCW2
      FRCTH1=FRCTH12
      FRCTH2=FRCTH22
      DO WHILE (FRCTIME2.LT.TIME)
      READ(101,*),FRCTIME2
      READ(102,*),FRCU2
      READ(103,*),FRCW2
      READ(104,*),FRCTH12
      READ(105,*),FRCTH22
      IF (FRCTIME2.GE.TIME) THEN
      U_BC_YMAX_VAL=(1.0d0-(TIME-FRCTIME)/(FRCTIME2-FRCTIME))
     &              *FRCU/NU/1022.8d0+
     & (1.0d0-(FRCTIME2-TIME)/(FRCTIME2-FRCTIME))
     &              *FRCU2/NU/1022.8d0
      W_BC_YMAX_VAL=(1.0d0-(TIME-FRCTIME)/(FRCTIME2-FRCTIME))
     &              *FRCW/NU/1022.8d0+ 
     & (1.0d0-(FRCTIME2-TIME)/(FRCTIME2-FRCTIME))
     &              *FRCW2/NU/1022.8d0
      TH_BC_YMAX_VAL(1)=(1.0d0-(TIME-FRCTIME)/(FRCTIME2-FRCTIME))
     &                  *FRCTH1/NU*PR(1)/(1022.8d0*4000.0d0)+ 
     & (1.0d0-(FRCTIME2-TIME)/(FRCTIME2-FRCTIME))
     &                  *FRCTH12/NU*PR(1)/(1022.8d0*4000.0d0)
! EMINUSP is forcing in units of kg fresh water per meter squared per
! second 
      TH_BC_YMAX_VAL(2)=(1.0d0-(TIME-FRCTIME)/(FRCTIME2-FRCTIME))
     &                  *FRCTH2/NU*PR(1)/(1000.0d0/31.5d0)+ 
     & (1.0d0-(FRCTIME2-TIME)/(FRCTIME2-FRCTIME))
     &                  *FRCTH22/NU*PR(1)/(1000.0d0/31.5d0)
      ELSE
      FRCTIME=FRCTIME2
      FRCU=FRCU2
      FRCW=FRCW2
      FRCTH1=FRCTH12
      FRCTH2=FRCTH22
      END IF
      END DO
      CLOSE(101)
      CLOSE(102)
      CLOSE(103)
      CLOSE(104)
      CLOSE(105)
      RETURN
      END
