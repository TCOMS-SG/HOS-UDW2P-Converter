!-------------------------------------------------------------------------
!   Copyright (C) 2021 -Technology Centre for Offshore and Marine, Singapore, 12 Prince George's Park, Singapore 118411
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program.  If not, see <https://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------    
PROGRAM  HOS_UDW2P_Converter
!-------------------------------------------------------------------------
! PROGRAM :     HOS_UDW2P_Converter
!-------------------------------------------------------------------------
!> @brief This program convert the HOS data to UDW2P file format
!> @details This program convert the HOS data to UDW2P file format
!> @remarks: 
!>  1.unit conversion: 
!>           The HOS_Ocean_SWENSE file output data is dimensional data
!>           The UDW2P data file prefer non-dimensionalized data
!>              -. Length : normalzied by water depth h
!>              -. Time   : normalzied by sqrt(h/g), where g is gravity.
!>           During the computation, all the data is dimensional data, will only no-dimensioalzied when output
!>           
!>  2.grid conversion:
!>           The UDW2P grid in x direction can be user defined in input_converter.dat. 
!>           The UDW2P grid in z direction is fixed to 20.
!>  3.time conversion:
!>           The UDW2P and HOS time step is the same.
!>
!> @author Xu Haihua, Technology Centre for Offshore and Marine, Singapore (TCOMS)
!> @date 17 Feb 2022
!>
!>
!> The HOS UDW2P Converter program was contributed by 
!> Dr Haihua Xu, Mr Yun Zhi Law and Dr Harrif Santo under the guidance of Prof Allan Magee. 
!> This research is supported by A∗STAR Science and Engineering Research Council with grant number 172 19 00089 under the Marine & Offshore Strategic Research Programme (M&O SRP). 
!> The source code is available at: https://github.com/ittcoms/HOS-UDW2P-Converter
!---------------------------------------------------------------------------------------------------------------------    
USE UDW2p_TPNWT  
USE Sigma_Grid_MOD
USE HOS_Ocean_SWENSE_MOD ,HOS_xlen=>xlen, HOS_ylen=>ylen,HOS_T_Stop=>T_Stop;
USE HOS_NWT_SWENSE_MOD,HOS_NWT_nx=>n1;   
USE OMP_LIB

IMPLICIT NONE
!--sigmaGrid will use dimensional data.
TYPE(TSigmaGrid2D):: sigmaGrid
!iHOS SWENSE type: 
!     =1 data is HOS Ocean
!     =2 data is HOS NWT(TBD)
INTEGER,PARAMETER:: HOS_OCEAN =1;
INTEGER,PARAMETER:: HOS_NWT =2;

INTEGER(4) :: iHOS =HOS_OCEAN;
!iOMPCores: number of cores for OpenMP computation.
INTEGER(4) :: iOMPCores = 6;


REAL(RK)    :: HOS_g,HOS_depth;

!nx,nz: number of grid in x and z direction.
INTEGER(4)  :: nx=1024,nz=20;
REAL(RK)    :: zp(100)
!xmin,xmax: minimum and maximum coordiante in UDW2P file format
!ymin,ymax: not used
!depth    : water depth (unit m)  
!tStart,tStop: start and end time to extract data from HOS to UDW2P
!dtOut    : UDW2P output time step which is them same as HOS
REAL(RK)    :: xmin,xmax,ymin,ymax,depth;
REAL(RK)    :: tStart,tStop; !start and stop convert time
REAL(RK)    :: dtOut         
REAL(RK)    :: dx, xCalib !xCalib is unclear assume it is zero
REAL(RK)    :: rTime; 
INTEGER(4)  :: iTime,iRecord;  
CHARACTER(400) :: filename="modes_HOS_SWENSE.dat"
INTEGER(4) :: iDebug = 1

REAL(RK)    :: xScale,tScale,vScale;

INTEGER(4)  :: i;
REAL(RK)    :: tmp;
REAL(RK),DIMENSION(:),ALLOCATABLE   ::elev,phi,leftPhi,rightPhi;
REAL(RK),DIMENSION(:),ALLOCATABLE   ::elevDt;
REAL(RK),DIMENSION(:),ALLOCATABLE   ::elevDt2;
REAL(RK),DIMENSION(:,:),ALLOCATABLE ::vel,leftVel,rightVel
REAL(RK)   :: ptmp,leftElev,leftXp,rightElev,rightXp,ztmp



INTEGER(4) :: unit,iUnit  
!-------CPU time varialbes
CHARACTER(LEN=6) ::cTime;
REAL(RK)   :: tCostS,tCostE;


!---read input data from file "input_converter.dat"
CALL read_input("input_converter.dat")


!Check wheter nz is equal to 20? current version nz==20 for UDW2P
IF( nz /=20) THEN 
    WRITE(*,"(A,I4)") "In current version nz = 20 is needed, nz =", nz;
    CALL sleep(10);
    STOP;
ENDIF

!--init the HOS Ocean Mode--
IF(iHOS == HOS_OCEAN) THEN
    CALL HOS_Ocean_init(TRIM(filename),dtOut,HOS_depth,HOS_g)  
ELSEIF(iHOS == HOS_NWT) THEN
    CALL HOS_NWT_init(TRIM(filename),dtOut,HOS_depth,HOS_g)
!--current HOS version has the following limition    
    xmin = 0.0d0;      !minimum set to 0
    xmax = HOS_xlen;   !maximum set to domain max
!    nx   = HOS_NWT_nx;         !nx set to total number of grid
ENDIF




!----compute the length scale, time scale in UDW2P
xScale = HOS_depth;             !length scale
tScale = SQRT(HOS_depth/HOS_g); !time scale
vScale = SQRT(HOS_depth*HOS_g)  !velocity scale 


!--check the input data 
IF(xmax > HOS_xlen)THEN
    WRITE(*,"(A,2F10.5)")"UDW2P xmax size larger than HOS domain size, xmax, HOS_xlen",xmax, HOS_xlen;
    CALL sleep(10);
    STOP;    
ENDIF

!--check the input data 
IF(tStop + tiny  > HOS_T_Stop)THEN
    WRITE(*,"(A,2F15.6)")"UDW2P maximum time larger than HOS maxium output time tStop, HOS_T_Stop",tStop, HOS_T_Stop;
    WRITE(*,"(A,2F15.6)") "Set maximum output time to HOS_T_Stop = ",HOS_T_Stop  ;
    tStop=  HOS_T_Stop
    !CALL sleep(10);
    !STOP;    
ENDIF

!--writhe number of grid for testing
WRITE(*,*) "Nx = ",nx;
WRITE(*,*) "Nz = ",nz;


!--set up the number of cores for OpenMP computation.
CALL omp_set_num_threads(iOMPCores);


!---get the vertical points distribution
CALL UDW2p_getGQPoints(zp)

!---start to record the time cost
tCostS = omp_get_wtime();


!--water depth in UDW2P, which is the same as HOS
depth = HOS_depth

!(nx,nz,xmin,xmax,h,zp)
!--create sigma coordinate in cartesian coordinate
CALL sigmaGrid%init(nx,nz,xmin,xmax,depth,zp)


ALLOCATE(elev(nx),elevDt(nx),vel(3,nx),leftVel(3,nz),rightVel(3,nz), &
         phi(nx),leftPhi(nz),rightPhi(nz))
         
!--elev_dt2 is used for test purpose
!ALLOCATE(elev_dt2(nx))
!WRITE(*,*)dt_out_star,T_stop_star,xlen_star,ylen_star,depth_star,g_star,L,T


!Now we start to process HOS-SWENSE results step by step.
!---------------set up the UDW2P file
dx     = sigmaGrid%dx   ;  !unit value
xCalib = 0.0d0

!CALL UDW2p_setUp(Nx_, Np_, Nt_, dx_, dt_, xCalib_, xScale_, tScale_)
CALL UDW2p_setUp(nx, nz, Nt_=0, dx_=dx/xScale, dt_=dtOut/tScale, xCalib_=xCalib, xScale_=xScale, tScale_=tScale)

!--first time to extract the data.
rTime = tStart;  !assign the time to convert.

itime = NINT(rTime/dtOut)+1 !

!-The main do loop.
rTime = (itime-1)*dtOut ;   !re-compute the time to HOS Output time

iRecord = 1;
DO WHILE (rTime < tStop )


!--read HOS data (modes) at T=rtime from the HOS file
    IF(iHOS == HOS_OCEAN) THEN
        CALL read_HOS_Ocean_mod(TRIM(filename),rTime);
    ELSEIF(iHOS == HOS_NWT) THEN
        CALL read_HOS_NWT_mod(TRIM(filename),rTime);    
    ENDIF
!----compute elevation and velocity at free surface. 
!$OMP PARALLEL DO     
    DO i=1,sigmaGrid%nx
!----compute the elevation at sigmaGrid
        IF(iHOS == HOS_OCEAN) THEN
            elev(i)   = get_HOS_Ocean_pt_eta(2,sigmaGrid%xp(i),0.0d0);
        
!---compute the time derivative of elevation.
            elevDt(i) = get_HOS_Ocean_pt_deta_dt(2,sigmaGrid%xp(i),0.0d0); 
        ELSEIF(iHOS == HOS_NWT) THEN
            elev(i)   = get_HOS_NWT_pt_eta(2,sigmaGrid%xp(i),0.0d0);   
            elevDt(i) = get_HOS_NWT_pt_deta_dt(2,sigmaGrid%xp(i),0.0d0);             
        ENDIF
!----compute the velocity and pressure at free surface for test purpose.
        IF(iDebug ==1) THEN  ! if debug mode, compute the velocity.
            IF(iHOS == HOS_OCEAN) THEN        
                CALL get_HOS_Ocean_pt_uvwp(2,sigmaGrid%xp(i),0.0d0,elev(i), vel(1,i),vel(2,i),vel(3,i),ptmp);
            ELSEIF(iHOS == HOS_NWT) THEN
                CALL get_HOS_NWT_pt_uvwp  (2,sigmaGrid%xp(i),0.0d0,elev(i), vel(1,i),vel(2,i),vel(3,i),ptmp);  

            ENDIF
        ENDIF
    ENDDO
!$OMP END PARALLEL DO
    
!----compute elevation and velocity at left boundary 
    leftXp   = sigmaGrid%xp(1)
    leftElev = elev(1)
    
!$OMP PARALLEL DO PRIVATE(ztmp)    
    DO i=1,sigmaGrid%nz
!---compute z coordiante     
        ztmp = sigmaGrid%S2Z(leftElev,sigmaGrid%zp(i));
        
!---compute phi (need to develop)        
!        leftPhi(i) = get_HOS_Ocean_pt_phi(2, leftXp,0.0d0,ztmp);     
        
!---compute velocity and pressure at (leftxp,0,ztmp) 
        IF(iHOS == HOS_OCEAN) THEN             
            CALL get_HOS_Ocean_pt_uvwp(2,leftXp,0.0d0,ztmp, leftVel(1,i),leftVel(2,i),leftVel(3,i),ptmp); 
        ELSEIF(iHOS == HOS_NWT) THEN
            CALL get_HOS_NWT_pt_uvwp  (2,leftXp,0.0d0,ztmp, leftVel(1,i),leftVel(2,i),leftVel(3,i),ptmp);          
        ENDIF            
    ENDDO
!$OMP END PARALLEL DO 
    
!----compute elevation and velocity at right boundary 
    rightXp   = sigmaGrid%xp(sigmaGrid%nx)
    rightElev = elev(sigmaGrid%nx)    
!$OMP PARALLEL DO PRIVATE(ztmp)      
    DO i=1,sigmaGrid%nz 
!---compute z coordiante       
        ztmp = sigmaGrid%S2Z(rightElev,sigmaGrid%zp(i));
        
!---compute phi (need to develop)           
!        rightPhi(i) = get_HOS_Ocean_pt_phi(2, rightXp,0.0d0,ztmp);   
        
!---compute velocity and pressure at (rightXp,0,ztmp) 
        IF(iHOS == HOS_OCEAN) THEN               
            CALL get_HOS_Ocean_pt_uvwp(2,rightXp,0.0d0,ztmp, rightVel(1,i),rightVel(2,i),rightVel(3,i),ptmp); 
        ELSEIF(iHOS == HOS_NWT) THEN
            CALL get_HOS_NWT_pt_uvwp  (2,rightXp,0.0d0,ztmp, rightVel(1,i),rightVel(2,i),rightVel(3,i),ptmp);       
        ENDIF            
    ENDDO   
!$OMP END PARALLEL DO 
    
!-----debug output the data  
!---export the unit value
    IF(iDebug ==1 .AND. mod(iRecord ,100) ==0) THEN
!     IF(iDebug ==1) THEN
        WRITE(cTime,"(I6.6)") iRecord
        
        OPEN(NEWUNIT = iUnit,FILE="debug_data/Debug_SF_"//TRIM(cTime)//".dat", ACTION="WRITE",STATUS="UNKNOWN")
        WRITE(iUnit,"(8A20)")  "x(m)","eta(m)","deta_dt(m/s)","u(m/s)","v(m/s)","w(m/s)"
        DO i=1,sigmaGrid%nx
            WRITE(iUnit,"(8E20.10)") sigmaGrid%xp(i),elev(i),elevDt(i),vel(:,i)
        ENDDO
        CLOSE(iUnit)
        
        OPEN(NEWUNIT = iUnit,FILE="debug_data/Debug_LF_"//TRIM(cTime)//".dat", ACTION="WRITE",STATUS="UNKNOWN")
        WRITE(iUnit,"(8A20)")  "x(m)","z(m)","u(m/s)","v(m/s)","w(m/s)"
        DO i=1,nz
            ztmp = sigmaGrid%S2Z(leftElev,sigmaGrid%zp(i));        
            WRITE(iUnit,"(6E20.10)") leftXp,ztmp,leftVel(:,i)
        ENDDO
        CLOSE(iUnit)
        
        OPEN(NEWUNIT = iUnit,FILE="debug_data/Debug_RF_"//TRIM(cTime)//".dat", ACTION="WRITE",STATUS="UNKNOWN")
        WRITE(iUnit,"(8A20)")  "x(m)","z(m)","u(m/s)","v(m/s)","w(m/s)"        
        DO i=1,nz
            ztmp = sigmaGrid%S2Z(rightElev,sigmaGrid%zp(i));        
            WRITE(iUnit,"(6E20.10)") rightXp,ztmp,rightVel(:,i)
        ENDDO
        CLOSE(iUnit)        
    ENDIF
    
    
!!$OMP END PARALLEL DO    
!----Write the data to UDW2P file .
!---normalize the data for output, add 28th-May-2020
    elev    = elev/xScale
    elevDt  = elevDt/vScale;
    leftVel = leftVel/vScale;
    rightVel= rightVel/vScale;
    
    CALL UDW2p_write_V2(elev,elevDt,leftVel(1,:),rightVel(1,:));

!update counter for the next record                
    itime = itime + 1;
    rTime = (itime-1)*dtOut ;  
    iRecord = iRecord +1;  
    
    IF(mod(iRecord ,100) ==0) THEN
        tCostE =  omp_get_wtime();
        WRITE(*,"(A,I6,A,F10.5)")"Pocressed itime = ",itime, " CPU time cost = ",   tCostE - tCostS
    ENDIF
    

ENDDO
!--close the file
CALL UDW2p_Close();

tCostE =  omp_get_wtime();
WRITE(*,"(A,F10.2,A)") "HOS data convert finished, total time cost is ",  tCostE - tCostS,"s"
CALL sleep(2);

STOP;

CONTAINS

SUBROUTINE read_input(inFile)

USE HOS_Ocean_SWENSE_MOD
IMPLICIT NONE
! Input variables
CHARACTER(LEN=*), INTENT(IN)  :: inFile
! Local variables
INTEGER  :: unit
!
unit = 100
OPEN(unit, FILE=inFile)
line_counter = 0

!CALL read_blank_line(unit)             ! --- HOS Type 
!CALL read_datum(unit, iHOS)            ! Minimum x in UDW2P grid

!
CALL read_blank_line(unit)             ! --- UDW2P Grid info
CALL read_datum(unit, xmin)            ! Minimum x in UDW2P grid
CALL read_datum(unit, xmax)            ! Maximum x in UDW2P grid
CALL read_datum(unit,nx)               ! Number of grid in x 
CALL read_datum(unit,nz)               ! Number of grid in z 
CALL read_datum(unit,tStart)           ! Start time convert 
CALL read_datum(unit,tStop)            ! End time convert 

WRITE(*,*)
CALL read_blank_line(unit)             ! ---  Input files 
CALL read_datum(unit, iHOS)            ! Type of HOS Data  
CALL read_datum(unit, filename)        ! Name of modal description

WRITE(*,*)
CALL read_blank_line(unit)             ! ---  Input files 
CALL read_datum(unit, iOMPCores)       ! Numer of cores used  
CALL read_datum(unit, iDebug)          ! Whether to output the data for debug 


CLOSE(unit)
!
END SUBROUTINE read_input

END PROGRAM


