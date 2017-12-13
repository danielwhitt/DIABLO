!  subroutines necessary for the WLT_BIO
#ifdef BIO
      SUBROUTINE INPUT_WLT_BIO
! read in reaction rate constants for the Whitt Levy Taylor 2016 NPZD model
      INCLUDE 'header'

      OPEN (12,file='input_wlt.dat',form='formatted',status='old')

      READ(12,*)
      READ(12,*)
      READ(12,*)
      READ(12,*) KWLIGHT, ALPHALIGHT, I0LIGHT
      READ(12,*)
      READ(12,*) VMNUT, KNNUT, SIGMADPHY
      READ(12,*)
      READ(12,*) IVLEVZOO, RZOO, GAMMANZOO, ZETAZOO
      READ(12,*)
      READ(12,*) ZETA2ZOO, DELTADET, WDET, DEEPN, NUDGETS 
      READ(12,*)

      CLOSE(12)

      END
!
      SUBROUTINE DEFINE_PAR
! define the 3-d par field for the biogeochemical model
! can be updated at each time step like the BCs if need be
! or initialized and held fixed
      INCLUDE 'header'
      INTEGER I,J,K
      
!       IF (RANKZ.eq.0) THEN
       WRITE(6,*) '**** DEFINING PAR USING ANALYTICAL FUNCTION ****'
!       END IF
       DO K=0,NZP-1
         DO I=0,NXM
           DO J=1,NY
             PAR(I,K,J)=I0LIGHT*EXP((GYF(J)-LY)*KWLIGHT)
!             IF ((RANKZ.eq.0)                                           &
!     &           .AND.(K.EQ.1).AND.(I.EQ.1)) THEN
!             WRITE(6,*) 'DEPTH:',GY(J),'PAR: ', PAR(I,K,J)
!             END IF
           END DO
         END DO
       END DO
!       IF (RANKZ.eq.0) THEN
!       WRITE(6,*) '**** END PAR PROFILE ****'
!       END IF
      END
!
      SUBROUTINE DEFINE_U2DET
! define the 3-d detrital sinking velocity for the biogeochemical model
! can be updated at each time step like the BCs if need be
! or initialized and held fixed
      INCLUDE 'header'
      INTEGER I,J,K

       WRITE(6,*) '** DEFINING U2DET USING ANALYTICAL FUNCTION **'
!       IF (RANKZ.eq.0) THEN
!       WRITE(6,*) '**** U2DET PROFILE ****'
!       END IF
       DO K=0,NZP-1
         DO I=0,NXM
           DO J=0,NY-1
             U2DET(I,K,J)=WDET
!             IF ((RANKZ.eq.0)                                           &
!     &           .AND.(K.EQ.1).AND.(I.EQ.1)) THEN
!             WRITE(6,*) 'DEPTH:',GYF(J),'U2DET: ', U2DET(I,K,J)
!             END IF
           END DO
         END DO
       END DO
       IF (RANKY.EQ.NPROCY-1) THEN
       DO K=0,NZP-1
         DO I=0,NXM
           DO J=NY,NY+1 
             U2DET(I,K,J)=0.0d0
!             IF ((RANKY.eq.(NPROCY-1)).AND.(RANKZ.eq.0)                 &
!     &           .AND.(K.EQ.1).AND.(I.EQ.1)) THEN
!             WRITE(6,*) 'DEPTH:',GYF(J),'U2DET: ', U2DET(I,K,J)
!             END IF
           END DO
         END DO
       END DO
       ELSE IF (RANKY.LT.(NPROCY-1)) THEN
       DO K=0,NZP-1
         DO I=0,NXM
           DO J=NY,NY+1
             U2DET(I,K,J)=WDET
           END DO
         END DO
       END DO
       END IF        
      END
!
      SUBROUTINE WLT_BIO_RHS
! Reaction terms only (sinking/swimming terms appear in advection operator)
       INCLUDE 'header'

       INTEGER I,J,K,N,CTFILL
       INTEGER NIX,PIX,ZIX,DIX
! N_TH indices of the nutrient (NIX), phytoplankton (PIX), zooplankton (ZIX), and detritus (DIX)
       NIX=2
       PIX=3
       ZIX=4
       DIX=5
       CTFILL=0
! remove/"fill" negative values
       DO J=JSTART_TH(NIX),JEND_TH(NIX)
         DO K=0,NZP-1
           DO I=0,NXM
            DO N=2,5
             IF (TH(I,K,J,N).LT.0.0d0) THEN
                 TH(I,K,J,N)=0.0d0
                 CTFILL=CTFILL+1
             END IF
            END DO
           END DO
         END DO
       END DO
!       IF (RANK.eq.0) THEN
!       WRITE(*,*) 'RANKY,RANKZ,FILLED TH:',RANKY,RANKZ,CTFILL
!       END IF
! RHS of nutrient equation
! DELTADET*D + GAMMANZOO*G*Z - U*P
!
       DO J=JSTART_TH(NIX),JEND_TH(NIX)
         DO K=0,NZP-1
           DO I=0,NXM
             S1(I,K,J)=DELTADET*TH(I,K,J,DIX)+GAMMANZOO*                &
     &                 RZOO*(1.0d0 - EXP(-IVLEVZOO*TH(I,K,J,PIX)))*     &
     &                 TH(I,K,J,ZIX) - TH(I,K,J,PIX)*                   &
     &                (VMNUT*TH(I,K,J,NIX)/(KNNUT + TH(I,K,J,NIX)))*    &
     &                  ALPHALIGHT*PAR(I,K,J)/                          &
     &        SQRT(VMNUT**2.0d0 + (ALPHALIGHT*PAR(I,K,J))**2.0d0)
              IF ((J.LT.3).AND.(RANKY.EQ.0)) THEN
              S1(I,K,J)=S1(I,K,J)-NUDGETS*(TH(I,K,J,NIX)-DEEPN)
              END IF
           END DO
         END DO
       END DO
! Convert to Fourier space
       CALL FFT_XZ_TO_FOURIER(S1,CS1,0,NY+1)
       DO J=JSTART_TH(NIX),JEND_TH(NIX)
         DO K=0,TNKZ
           DO I=0,NXP-1
             CFTH(I,K,J,NIX)=CFTH(I,K,J,NIX)+CS1(I,K,J)
           END DO
         END DO
       END DO
! Second RHS of P eqn
! U*P - G*Z - SIGMAD*P
!
       DO J=JSTART_TH(PIX),JEND_TH(PIX)
         DO K=0,NZP-1
           DO I=0,NXM
             S1(I,K,J)=-SIGMADPHY*TH(I,K,J,PIX)                         &
     &                -RZOO*(1.0d0 - EXP(-IVLEVZOO*TH(I,K,J,PIX)))*     &
     &                 TH(I,K,J,ZIX) + TH(I,K,J,PIX)*                   &
     &                (VMNUT*TH(I,K,J,NIX)/(KNNUT + TH(I,K,J,NIX)))*    &
     &                  ALPHALIGHT*PAR(I,K,J)/                          &
     &        SQRT(VMNUT**2.0d0 + (ALPHALIGHT*PAR(I,K,J))**2.0d0)
           END DO
         END DO
       END DO
! Convert to Fourier space
       CALL FFT_XZ_TO_FOURIER(S1,CS1,0,NY+1)
       DO J=JSTART_TH(PIX),JEND_TH(PIX)
         DO K=0,TNKZ
           DO I=0,NXP-1
             CFTH(I,K,J,PIX)=CFTH(I,K,J,PIX)+CS1(I,K,J)
           END DO
         END DO
       END DO
! Third RHS of Z eqn
! (1-GAMMANZOO)*G*Z - ZETAZOO*Z-ZETA2ZOO*Z^2
!
       DO J=JSTART_TH(ZIX),JEND_TH(ZIX)
         DO K=0,NZP-1
           DO I=0,NXM
              S1(I,K,J)=                                                &
     &   -ZETAZOO*TH(I,K,J,ZIX)-ZETA2ZOO*(TH(I,K,J,ZIX)**2.0d0)         &
     &   +(1.0d0-GAMMANZOO)*RZOO*(1.0d0-EXP(-IVLEVZOO*TH(I,K,J,PIX)))   &
     &                 *TH(I,K,J,ZIX)
           END DO
         END DO
       END DO
! Convert to Fourier space
       CALL FFT_XZ_TO_FOURIER(S1,CS1,0,NY+1)
       DO J=JSTART_TH(ZIX),JEND_TH(ZIX)
         DO K=0,TNKZ
           DO I=0,NXP-1
             CFTH(I,K,J,ZIX)=CFTH(I,K,J,ZIX)+CS1(I,K,J)
           END DO
         END DO
       END DO
! Fourth RHS of D eqn
! + ZETAZOO*Z+ZETA2ZOO*Z^2+SIGMADPHY*P-DELTADET*D
!     
       DO J=JSTART_TH(DIX),JEND_TH(DIX)
         DO K=0,NZP-1
           DO I=0,NXM
              S1(I,K,J)= SIGMADPHY*TH(I,K,J,PIX)                        &
     &   +ZETAZOO*TH(I,K,J,ZIX)+ZETA2ZOO*(TH(I,K,J,ZIX)**2.0d0)         &
     &   -DELTADET*TH(I,K,J,DIX)
           END DO
         END DO
       END DO
! Convert to Fourier space
       CALL FFT_XZ_TO_FOURIER(S1,CS1,0,NY+1)
       DO J=JSTART_TH(DIX),JEND_TH(DIX)
         DO K=0,TNKZ
           DO I=0,NXP-1
             CFTH(I,K,J,DIX)=CFTH(I,K,J,DIX)+CS1(I,K,J)
           END DO
         END DO
       END DO
! Add sinking term to the RHS of the detritus equation
! this is now added in RK_CHAN
!      DO J=JSTART_TH(DIX),JEND_TH(DIX)
!        DO K=0,NZP-1
!          DO I=0,NXM
!            S1(I,K,J)=
!     & (TH(I,K,J+1,DIX)*U2DET(I,K,J+1)                                  &
!     & +TH(I,K,J,DIX)*U2DET(I,K,J+1)                                    &
!     & -TH(I,K,J,DIX)*U2DET(I,K,J)                                      &
!     & -TH(I,K,J-1,DIX)*U2DET(I,K,J))/(2.d0*DYF(J))
!          END DO
!        END DO
!      END DO
!      CALL FFT_XZ_TO_FOURIER(S1,CS1,0,NY+1)
!      DO J=JSTART_TH(DIX),JEND_TH(DIX)
!        DO K=0,TNKZ
!          DO I=0,NXP-1
!            CFTH(I,K,J,DIX)=CFTH(I,K,J,DIX) - CS1(I,K,J)
!          END DO
!        END DO
!      END DO
      END
#endif
