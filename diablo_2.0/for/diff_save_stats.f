9c9
<       CHARACTER*10 GNAME
---
>       CHARACTER*20 GNAME
16d15
<  
19a19,21
> ! These variable are used for HDF5 writing
>       real*8 Diag(1:NY)
> 
38a41,69
> #ifdef HDF5
>         FNAME='stats.h5'
>         if (USE_MPI) then
>           call mpi_barrier(MPI_COMM_WORLD,ierror)
>         end if
>         IF (RANKZ.EQ.0) THEN
>           Diag=GYF(1:NY)
>           gname='GYF'
>           call WriteStatH5(FNAME,gname,Diag)
>           Diag=UBAR(1:NY)
>           gname='UBAR'
>           call WriteStatH5(FNAME,gname,Diag)
>           Diag=VBAR(1:NY)
>           gname='VBAR'
>           call WriteStatH5(FNAME,gname,Diag)
>           Diag=WBAR(1:NY)
>           gname='WBAR'
>           call WriteStatH5(FNAME,gname,Diag)
>           do n=1,N_TH
>             Diag=THBAR(1:NY,n)
>             gname='THBAR'
>      &           //CHAR(MOD(N,100)/10+48)
>      &           //CHAR(MOD(N,10)+48)
>             call WriteStatH5(FNAME,gname,Diag)
>           end do
>         END IF
> 
> #else
> ! Here we aren't using HDF5, so save to text files
45d75
< 
58a89
> #endif
225,258d255
< ! Calculate the mean viscous fluxes F_x = <(NU+NU_T)*du/dy>
< ! F_y = <(NU+NU_T)*dw/dy> 
<       do j=1,NY+1
<         nutdudy(j)=0.d0
<         nutdugdy(j)=0.d0
<         nutdwdy(j)=0.d0
<         do k=0,NZP-1
<           do i=0,NXM
<             nutdudy(j)=nutdudy(j)
<      &       +(NU+NU_T(i,k,j))*(U1(i,k,j)-U1(i,k,j-1))/(DY(j))
<             nutdugdy(j)=nutdugdy(j)
<      &       +(NU+NU_T(i,k,j))*(-DRHODZ(1)*RI(1)/I_RO)
<             nutdwdy(j)=nutdwdy(j)
<      &       +(NU+NU_T(i,k,j))*(U3(i,k,j)-U3(i,k,j-1))/(DY(j))
<           end do
<         end do
<       end do
<       call mpi_allreduce(mpi_in_place,nutdudy,NY+2,MPI_DOUBLE_PRECISION,
<      &     MPI_SUM,MPI_COMM_Z,ierror)
<       do j=0,NY+1
<         nutdudy(j)=nutdudy(j)/dble(NX*NZ)
<       end do
<       call mpi_allreduce(mpi_in_place,nutdugdy,NY+2,
<      &                          MPI_DOUBLE_PRECISION,
<      &     MPI_SUM,MPI_COMM_Z,ierror)
<       do j=0,NY+1
<         nutdugdy(j)=nutdugdy(j)/dble(NX*NZ)
<       end do
<       call mpi_allreduce(mpi_in_place,nutdwdy,NY+2,MPI_DOUBLE_PRECISION,
<      &     MPI_SUM,MPI_COMM_Z,ierror)
<       do j=0,NY+1
<         nutdwdy(j)=nutdwdy(j)/dble(NX*NZ)
<       end do
< 
343a341,420
> #ifdef HDF5
>       FNAME='mean.h5'
> 
>       gname='time'
>       call WriteHDF5_real(FNAME,gname,TIME)
> 
>       IF (RANKZ.eq.0) then
> 
>       gname='gyf'
>       Diag=gyf(1:NY)      
>       call WriteStatH5(FNAME,gname,Diag)
> 
>       gname='ume'
>       Diag=ume(1:NY)
>       call WriteStatH5(FNAME,gname,Diag)
>     
>       gname='vme'
>       Diag=vme(1:NY)
>       call WriteStatH5(FNAME,gname,Diag)
> 
>       gname='wme'
>       Diag=wme(1:NY)
>       call WriteStatH5(FNAME,gname,Diag)
> 
>       gname='urms'
>       Diag=urms(1:NY)
>       call WriteStatH5(FNAME,gname,Diag)
> 
>       gname='vrms'
>       Diag=vrms(1:NY)
>       call WriteStatH5(FNAME,gname,Diag)
> 
>       gname='wrms'
>       Diag=wrms(1:NY)
>       call WriteStatH5(FNAME,gname,Diag)
> 
>       gname='uv'
>       Diag=uv(1:NY)
>       call WriteStatH5(FNAME,gname,Diag)
> 
>       gname='uw'
>       Diag=uw(1:NY)
>       call WriteStatH5(FNAME,gname,Diag)
> 
>       gname='wv'
>       Diag=wv(1:NY)
>       call WriteStatH5(FNAME,gname,Diag)
> 
>       gname='dudy'
>       Diag=dudy(1:NY)
>       call WriteStatH5(FNAME,gname,Diag)
> 
>       gname='dwdy'
>       Diag=dwdy(1:NY)
>       call WriteStatH5(FNAME,gname,Diag)
> 
>       gname='cp'
>       Diag=dble(cp(0,0,1:NY))
>       call WriteStatH5(FNAME,gname,Diag)
> 
>       gname='shear'
>       Diag=shear(1:NY)
>       call WriteStatH5(FNAME,gname,Diag)
> 
>       gname='omega_x'
>       Diag=omega_x(1:NY)
>       call WriteStatH5(FNAME,gname,Diag)
> 
>       gname='omega_y'
>       Diag=omega_y(1:NY)
>       call WriteStatH5(FNAME,gname,Diag)
> 
>       gname='omega_z'
>       Diag=omega_y(1:NY)
>       call WriteStatH5(FNAME,gname,Diag)
> 
>       END IF
> 
> #else
> ! Here are aren't using HDF5, so write mean statistcs to text files
360d436
<      &      ,nutdudy(j),nutdugdy(j),nutdwdy(j)
362a439,440
> 401   format(I3,' ',17(F30.20,' '))
> #endif
364,373c442,444
< 401   format(I3,' ',20(F30.20,' '))
< 
< 
< ! DAN Write out tkebudget BEFORE writing over CR1,2,3 in the TH loop
<       call mpi_barrier(MPI_COMM_WORLD,ierror)
<       CALL tkebudget_chan
<       IF (LES) THEN 
<       CALL tkebudget_chan_les
<       END IF
<       call mpi_barrier(MPI_COMM_WORLD,ierror)
---
> ! Calculate the dissipation rate
>       call tkebudget_chan
>       IF (LES) call tkebudget_chan_les
440,455d510
< ! Compute mean diffusive flux
<       do j=1,NY+1
<         thsum(j)=0.
<       do k=0,NZP-1
<       do i=0,NXM
<        thsum(j)=thsum(j)+(KAPPA_T(i,k,j,n)+
<      &                    +NU/PR(n))
<      &    *(TH(i,k,j,n)-TH(i,k,j-1,n))/DY(j)
<       end do
<       end do
<       end do
<       call mpi_allreduce(mpi_in_place,thsum,(NY+2),
<      &     MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_Z,ierror)
<       do j=0,NY+1
<       kappadthdy(j,n)=thsum(j)/(float(NZ)*float(NX))
<       end do
551,655d605
<       else if (n.eq.3) then
<          if (USE_MPI) then
<          call mpi_barrier(MPI_COMM_WORLD,ierror)
<          end if
<          IF (RANKZ.EQ.RANKZMOVIE) THEN
<             do I=0,NXM
<             do J=1,NY
<                varxy(i,j)=TH(i,NzMovie,j,n)
<             end do
<             end do
<             GNAME='th3_xy'
<             call writeHDF5_xyplane(FNAME,GNAME,varxy)
<          END IF
<          if (USE_MPI) then
<          call mpi_barrier(MPI_COMM_WORLD,ierror)
<          end if
<          IF (RANKY.EQ.RANKYMOVIE) THEN
<             do I=0,NXM
<             do J=0,NZP-1
<                varxz(i,j)=TH(i,j,NyMovie,n)
<             end do
<             end do
<             GNAME='th3_xz'
<             call writeHDF5_xzplane(FNAME,GNAME,varxz)
<          END IF
<          if (USE_MPI) then
<          call mpi_barrier(MPI_COMM_WORLD,ierror)
<          end if
<          do I=0,NZP-1
<          do J=1,NY
<             varzy(i,j)=TH(NxMovie,i,j,n)
<          end do
<          end do
<          GNAME='th3_zy'
<          call writeHDF5_zyplane(FNAME,GNAME,varzy)
<       else if (n.eq.4) then
<          if (USE_MPI) then
<          call mpi_barrier(MPI_COMM_WORLD,ierror)
<          end if
<          IF (RANKZ.EQ.RANKZMOVIE) THEN
<             do I=0,NXM
<             do J=1,NY
<                varxy(i,j)=TH(i,NzMovie,j,n)
<             end do
<             end do
<             GNAME='th4_xy'
<             call writeHDF5_xyplane(FNAME,GNAME,varxy)
<          END IF
<          if (USE_MPI) then
<          call mpi_barrier(MPI_COMM_WORLD,ierror)
<          end if
<          IF (RANKY.EQ.RANKYMOVIE) THEN
<             do I=0,NXM
<             do J=0,NZP-1
<                varxz(i,j)=TH(i,j,NyMovie,n)
<             end do
<             end do
<             GNAME='th4_xz'
<             call writeHDF5_xzplane(FNAME,GNAME,varxz)
<          END IF
<          if (USE_MPI) then
<          call mpi_barrier(MPI_COMM_WORLD,ierror)
<          end if
<          do I=0,NZP-1
<          do J=1,NY
<             varzy(i,j)=TH(NxMovie,i,j,n)
<          end do
<          end do
<          GNAME='th4_zy'
<          call writeHDF5_zyplane(FNAME,GNAME,varzy)
<       else if (n.eq.5) then
<          if (USE_MPI) then
<          call mpi_barrier(MPI_COMM_WORLD,ierror)
<          end if
<          IF (RANKZ.EQ.RANKZMOVIE) THEN
<             do I=0,NXM
<             do J=1,NY
<                varxy(i,j)=TH(i,NzMovie,j,n)
<             end do
<             end do
<             GNAME='th5_xy'
<             call writeHDF5_xyplane(FNAME,GNAME,varxy)
<          END IF
<          if (USE_MPI) then
<          call mpi_barrier(MPI_COMM_WORLD,ierror)
<          end if
<          IF (RANKY.EQ.RANKYMOVIE) THEN
<             do I=0,NXM
<             do J=0,NZP-1
<                varxz(i,j)=TH(i,j,NyMovie,n)
<             end do
<             end do
<             GNAME='th5_xz'
<             call writeHDF5_xzplane(FNAME,GNAME,varxz)
<          END IF
<          if (USE_MPI) then
<          call mpi_barrier(MPI_COMM_WORLD,ierror)
<          end if
<          do I=0,NZP-1
<          do J=1,NY
<             varzy(i,j)=TH(NxMovie,i,j,n)
<          end do
<          end do
<          GNAME='th5_zy'
<          call writeHDF5_zyplane(FNAME,GNAME,varzy)
656a607
> 
666a618,663
> 
> 
> #ifdef HDF5
>  
>       FNAME='mean.h5' 
>   
>       IF (RANKZ.eq.0) THEN
>  
>       do n=1,N_TH
>  
>         Diag=thme(1:NY,n)
>         gname='thme'
>      &           //CHAR(MOD(N,100)/10+48)
>      &           //CHAR(MOD(N,10)+48)
>         call WriteStatH5(FNAME,gname,Diag)
> 
>         Diag=dthdy(1:NY,n)
>         gname='dthdy'
>      &           //CHAR(MOD(N,100)/10+48)
>      &           //CHAR(MOD(N,10)+48)
>         call WriteStatH5(FNAME,gname,Diag)
> 
>         Diag=thrms(1:NY,n)
>         gname='thrms'
>      &           //CHAR(MOD(N,100)/10+48)
>      &           //CHAR(MOD(N,10)+48)
>         call WriteStatH5(FNAME,gname,Diag)
> 
>         Diag=thv(1:NY,n)
>         gname='thv'
>      &           //CHAR(MOD(N,100)/10+48)
>      &           //CHAR(MOD(N,10)+48)
>         call WriteStatH5(FNAME,gname,Diag)
> 
>         Diag=pe_diss(1:NY,n)
>         gname='pe_diss'
>      &           //CHAR(MOD(N,100)/10+48)
>      &           //CHAR(MOD(N,10)+48)
>         call WriteStatH5(FNAME,gname,Diag)
>  
>       end do
> 
>       END IF
> 
> #else
> ! Here we aren't using HDF5, so write to a text file
680c677
<      +      ,dthdy(j,n),thrms(j,n),thv(j,n),pe_diss(j,n),kappadthdy(j,n)
---
>      +      ,dthdy(j,n),thrms(j,n),thv(j,n),pe_diss(j,n)
684,685c681,682
< 
< 402   format(I3,' ',7(F30.20,' '))
---
> 402   format(I3,' ',6(F30.20,' '))
> #endif
758,769d754
<          if (USE_MPI) then
<          call mpi_barrier(MPI_COMM_WORLD,ierror)
<          end if
<          IF (RANKZ.EQ.RANKZMOVIE) THEN
<             do I=0,NXM
<             do J=1,NY
<                varxy(i,j)=KAPPA_T(i,NzMovie,j,1)
<             end do
<             end do
<             GNAME='kappa_t_xy'
<             call writeHDF5_xyplane(FNAME,GNAME,varxy)
<          END IF
824,835d808
<          if (USE_MPI) then
<          call mpi_barrier(MPI_COMM_WORLD,ierror)
<          end if
<          IF (RANKY.EQ.RANKYMOVIE) THEN
<             do I=0,NXM
<             do J=0,NZP-1
<               varxz(i,j)=KAPPA_T(i,j,NyMovie,1)
<             end do
<             end do
<             GNAME='kappa_t_xz'
<             call writeHDF5_xzplane(FNAME,GNAME,varxz)
<          end if
871d843
<          IF (LES) then
874a847
>          IF (LES) then
882,892d854
< 
<          if (USE_MPI) then
<          call mpi_barrier(MPI_COMM_WORLD,ierror)
<          end if
<          do I=0,NZP-1
<          do J=1,NY
<             varzy(i,j)=KAPPA_T(NxMovie,i,j,1)
<          end do
<          end do
<          GNAME='kappa_t_zy'
<          call writeHDF5_zyplane(FNAME,GNAME,varzy)
906a869
> ! END IF FINAL
918a882,959
>       subroutine tkebudget_chan_les
> ! Calculate the componet of th SGS dissipation rate 
> ! only includes the terms timestepped implicitly
>       include 'header'
>       include 'header_les'
> 
>       character*35 FNAME
>       CHARACTER*20 GNAME
>       real*8 epsilon_sgs(NY)
>       real*8 Diag(NY)
>       integer i,j,k
> 
> ! Compute the turbulent dissipation rate, epsilon=nu*<du_i/dx_j du_i/dx_j>
> 
>       DO J=1,NY
>         DO K=0,NZP-1
>           DO I=0,NXM
>             TEMP(I,K,J)=U1(I,K,J)*
>      &        (  (NU_T(I,K,J+1) * (U1(I,K,J+1) - U1(I,K,J)) / DY(J+1)
>      &         -  NU_T(I,K,J) * (U1(I,K,J)   - U1(I,K,J-1)) / DY(J))
>      &               /DYF(J)  )
>      &           +U3(I,K,J)*
>      &        (  (NU_T(I,K,J+1) * (U3(I,K,J+1) - U3(I,K,J)) / DY(J+1)
>      &        - NU_T(I,K,J) * (U3(I,K,J)   - U3(I,K,J-1)) / DY(J))
>      &              /DYF(J)  )
>      &           +U2(I,K,J)*
>      &     ((0.5d0*(NU_T(I,K,J)+NU_T(I,K,J+1))*(U2(I,K,J+1)-U2(I,K,J))
>      &                                              / DYF(J)
>      &    -0.5d0*(NU_T(I,K,J)+NU_T(I,K,J-1))*(U2(I,K,J)-U2(I,K,J-1))
>      &                                          / DYF(J-1))   /DY(J)  )
>           END DO
>         END DO
>       END DO
>  
> ! Now calculate the horizontal average
>         do j=1,NY
>           epsilon_sgs(j)=0.d0
>           do i=0,NXM
>           do k=0,NZP-1
>             epsilon_sgs(j)=epsilon_sgs(j)+TEMP(I,K,J)
>           end do
>           end do
>         end do
> 
>       call mpi_allreduce(mpi_in_place,epsilon_sgs,NY+2
>      &    ,MPI_DOUBLE_PRECISION,
>      &     MPI_SUM,MPI_COMM_Z,ierror)        
> 
> 
> #ifdef HDF5
>       FNAME='tke.h5'
> 
>       IF (RANKZ.eq.0) THEN
>         gname='epsilon_sgs'
>         Diag=epsilon_sgs(1:NY)
>         call WriteStatH5(FNAME,gname,Diag)
>       END IF
> #else
> ! Here are aren't using HDF5, so write to text files
>       IF (RANKZ.EQ.0) THEN
>       IF (USE_MPI) THEN
>         FNAME='tke_les'//trim(MPI_IO_NUM)//'.txt'
>       ELSE
>         FNAME='tke_les.txt'
>       END IF
>       open(46,file=FNAME,form='formatted',status='unknown')
> 
>       write(46,*) TIME_STEP,TIME,DELTA_T
>         do j=1,NY
>           write(46,460) j,GYF(J),epsilon_sgs(J)
>         end do
>       END IF
> 460     format(I3,' ',2(F30.20,' '))
> #endif
> 
>       END
> 
> 
924a966,967
>       character*20 GNAME
>       real*8 Diag(1:NY)
927a971
> 
1089c1133
<         epsilon(j)=NU*epsilon(j)/(float(NX)*float(NZ))
---
>         epsilon(j)=epsilon(j)/float(NX*NZ)
1093a1138,1150
> #ifdef HDF5
>       FNAME='tke.h5'
> 
>       IF (RANKZ.eq.0) THEN
>         epsilon(1)=0.d0
>         epsilon(NY)=0.d0      
>         gname='epsilon' 
>         Diag=epsilon(1:NY)
>         call WriteStatH5(FNAME,gname,Diag)
>       END IF
> 
> #else
> ! Here we aren't using HDF5, so write to a text file
1106c1163
< 401   format(I3,' ',2(F30.20,' '))
---
> 401   format(I3,' ',2(F20.9,' '))
1107a1165
> #endif
