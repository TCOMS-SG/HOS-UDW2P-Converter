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
MODULE  HOS_Ocean_SWENSE_MOD
!This module contains the subroutine to read HOS SWENSE file and compute the surface elevation and velocity 
!from the file.
!    
!This subroutine is referenced and revised from HOS Ocean software  https://github.com/LHEEA/HOS-ocean
!The Copyright/License of HOS-Ocean is present here 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    Copyright (C) 2014 - LHEEA Lab., Ecole Centrale de Nantes, UMR CNRS 6598
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!USE utilities_mod    
!USE variables_MOD
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------    
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


!HOSFile Name         : The file name of HOS SWESEN, default is HOSFileName ="modes_HOS_SWENSE.dat"
!read_HOS_Ocean_mod   : read FFT modes from HOS Ocean SWENSE file.
!HOS_Ocean_init       : init the HOS Ocen SWENSE file 
!get_HOS_Ocean_pt_eta : compute eta (free surface elevation) for a given point.
!get_HOS_Ocean_pt_uvwp: compute u,v,w and p (velocity and pressure) for a given point.
!get_HOS_Ocean_pt_deta_dt: compute d(eta)/dt
PUBLIC :: HOSFileName,get_HOS_Ocean_pt_eta, get_HOS_Ocean_pt_uvwp,read_HOS_Ocean_mod,HOS_Ocean_init  
PUBLIC :: get_HOS_Ocean_pt_deta_dt
PUBLIC :: int2str
PUBLIC :: Tiny
!depth :water depth of the domain
!T_Stop:stop time. 
!xlen,ylen :domain length
PUBLIC :: xlen,ylen,T_stop ; 

!PUBLIC :: get_HOS_Ocean_pt_phi
PUBLIC :: read_datum, read_blank_line
PUBLIC :: line_counter


INTEGER(4) :: line_counter=0;


!iHOSModeMethod: the method to use HOS modes
!=1: use all the modes kx,ky in HOS model
!=2: use the maximum nmodes in HOS model
INTEGER,PUBLIC       :: iHOSModeMethod = 1;
!if iHOSModeMethod=2, will only use maximum nHOSModeUsed in the ele, u,v,w and p computation. 
INTEGER,PUBLIC       :: nHOSModeUsed =128 !

!--The 
!=.TRUE. The HOS mode is selected, hence
LOGICAL,PRIVATE      ::lHOSModeSelected =.FALSE.

PRIVATE

INTERFACE get_HOS_Ocean_pt_eta
    MODULE PROCEDURE get_pt_eta
END INTERFACE

!INTERFACE get_HOS_Ocean_pt_phi
!    MODULE PROCEDURE get_pt_phi
!END INTERFACE
INTERFACE get_HOS_Ocean_pt_deta_dt
    MODULE PROCEDURE get_pt_deta_dt
END INTERFACE


INTERFACE get_HOS_Ocean_pt_uvwp
    MODULE PROCEDURE get_pt_uvwp
END INTERFACE 

INTERFACE read_HOS_Ocean_mod
    MODULE PROCEDURE read_mod
END INTERFACE

INTERFACE HOS_Ocean_init
    MODULE PROCEDURE recons_HOS_init
END INTERFACE

!PUBLIC :: recons_HOS_init, get_pt_eta,get_pt_uvwp,read_mod 
!iModeIndexNx: modes index selected for X direction
INTEGER,DIMENSION(:,:),ALLOCATABLE::iModeUsedIndex

CHARACTER(LEN=40)        :: HOSFileName ="modes_HOS_SWENSE.dat"
INTEGER,PARAMETER,PRIVATE:: ntn=2;

REAL(8),PRIVATE :: rPFResultScale =1.0d0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    Copyright (C) 2014 - LHEEA Lab., Ecole Centrale de Nantes, UMR CNRS 6598
!
!    This program is part of HOS-ocean
!
!    HOS-ocean is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Definition of symbols for real types (RP) and complex ones (CP)
!
! Real numbers are simple or double precision
INTEGER, PARAMETER :: SP = SELECTED_REAL_KIND(6,   37)     ! REAL32
INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(15, 307)     ! REAL64
!
! Complex numbers are simple or double precision
INTEGER, PARAMETER :: SPC = SP
INTEGER, PARAMETER :: DPC = DP
!
! Current types
INTEGER, PARAMETER :: RP = DP
INTEGER, PARAMETER :: CP = DPC
!
! Define usual mathematical constants i, pi, 2pi, pi/2, square root of 2.
COMPLEX(CP), PARAMETER :: i     = ((0.0_rp, 1.0_rp))
REAL(RP), PARAMETER    :: PI    = 3.141592653589793238462643383279502888419_rp
!REAL(RP), PARAMETER    :: g     = 9.81_rp
REAL(RP), PARAMETER    :: PIO2  = 1.570796326794896619231321691639751442098_rp
REAL(RP), PARAMETER    :: TWOPI = 6.283185307179586476925286766559005768394_rp
REAL(RP), PARAMETER    :: SQ2   = 1.414213562373095048801688724209698078569_rp
!
! For comparison of real numbers
!REAL(RK),PARAMETER:: EPSILON    =1.0E-6_RK
REAL(RP), PARAMETER    :: tiny = 1.0E-8_RP ! epsilon(1.0_rp)
!

INTEGER :: n1,n2
REAL(RP), ALLOCATABLE, DIMENSION(:)   :: x,y
REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: eta, phis
!
REAL(RP), ALLOCATABLE, DIMENSION(:)      :: kx,ky_n2
REAL(RP), ALLOCATABLE, DIMENSION(:,:)    :: kth
COMPLEX(CP), ALLOCATABLE, DIMENSION(:,:) :: ikx,iky
!
! Input file
INTEGER               :: i_ana, i_card, tecplot, i_zvect

CHARACTER(LEN=100)    :: file_3d, file_mod
!
INTEGER, PARAMETER    :: n_hdr = 34 ! Headerlines in '3d.dat' including line variables
REAL(RP), PARAMETER   :: HfoHs = 2.0_rp ! Freak wave height on Hs threshold for detection

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
INTEGER :: i_xvect,i_yvect,nzvect

INTEGER :: n1o2p1,n2o2p1
! test for fourier
INTEGER :: m1,m2,Nd1,Nd2,Nd1o2p1,m1o2p1,md1o2p1,md1,md2


INTEGER :: i_Test=0;   ! =1 compute vertical profile on point for test
INTEGER :: imode =1    ! =1: use all modes to compute eta,u,v,w and p; =2: use maximum ncomp modes only
INTEGER :: ncomp=128;
INTEGER :: ncores=12;
!

COMPLEX(CP), DIMENSION(:,:,:),ALLOCATABLE,TARGET :: modesspecx_a,modesspecy_a,modesspecz_a,modesspect_a,modesFS_a,modesFSt_a
!REAL(RP),    ALLOCATABLE, DIMENSION(m1o2p1,m2) :: anglex, angley, anglez, anglet, angleut, anglevt, anglewt
REAL(RP),    DIMENSION(:,:,:),ALLOCATABLE,TARGET  :: anglex_a, angley_a, anglez_a, anglet_a, angleut_a, anglevt_a, anglewt_a

COMPLEX(CP), DIMENSION(:,:),POINTER  :: modesspecx,modesspecy,modesspecz,modesspect,modesFS,modesFSt
!REAL(RP),    ALLOCATABLE, DIMENSION(m1o2p1,m2) :: anglex, angley, anglez, anglet, angleut, anglevt, anglewt
REAL(RP),    DIMENSION(:,:),POINTER  :: anglex, angley, anglez, anglet, angleut, anglevt, anglewt


!---additional variables.
!*_star: non-dimensional variables.
!
REAL(RP)        :: dt_out_star,T_stop_star,xlen_star,ylen_star,depth_star,g_star,L,T;
REAL(RP)        :: T_start, T_stop, x_min, x_max, y_min, y_max, z_min, z_max
REAL(RP)        :: density=1000.0_RP;
REAL(RP)        :: depth,g,xlen,ylen;

!----bubble sort with index as return value
INTERFACE sortBubble
    MODULE PROCEDURE R8SORT_BUBBLE
END INTERFACE sortBubble

INTEGER, PARAMETER            :: N_descr = 40
INTEGER, PARAMETER            :: N_tot   = 52
INTEGER, PARAMETER            :: len_form_read = 7
INTEGER, PARAMETER            :: len_form_write = 25

INTERFACE read_datum
     MODULE PROCEDURE read_datum_i, read_datum_r, read_datum_c
END INTERFACE

CHARACTER(LEN=len_form_read)  :: format_read(0:4)
CHARACTER(LEN=N_tot)          :: description
CONTAINS
!int2str: convert int to string
!-------------------------------------------------------------
!SUBROUTINE         int2str
!-------------------------------------------------------------
!> @brief convert integer data to string type data
!> @details convert integer data to string type data
!> @param[in] int integer to convert
!> @results The converted string results
!
!> @author Xu Haihua, Technology Centre for Offshore and Marine, Singapore (TCOMS)
!> @date 11 Nov 2019 
!> @remarks None
!-------------------------------------------------------------
CHARACTER(LEN=12) FUNCTION  int2str(int)
!CHARACTER(LEN=*) :: int2str
INTEGER          :: int
WRITE(int2str,'(I12)') int
int2str = TRIM(ADJUSTL(int2str))
END FUNCTION int2str


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    Copyright (C) 2014 - LHEEA Lab., Ecole Centrale de Nantes, UMR CNRS 6598
!
!    This program is part of HOS-ocean
!
!    HOS-ocean is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!USE HOS_type
!


!CONTAINS
!
LOGICAL FUNCTION iseven(n)
!
IMPLICIT NONE
!
INTEGER :: n
!
iseven = (MOD(n,2) == 0)
!
END FUNCTION iseven
!
!-------------------------------------------------------------
!SUBROUTINE         init_read_mod
!-------------------------------------------------------------
!> @brief initial the module and read the parameters from the HOS SWESEN file
!> @details  initial the module and read the parameters from the HOS SWESEN file
!> @param[in] fileName: The HOS SWESEN file name
!> @param[OUT] n1: Number of modes in x direction
!> @param[OUT] n2: Number of modes in y direction
!> @param[OUT] dt_out: Time step for output (non-dimensional)
!> @param[OUT] T_stop: Maximum time in the SWESEN file (non-dimensional)
!> @param[OUT] xlen: Length in x direction (non-dimensional)
!> @param[OUT] ylen: Length in y direction (non-dimensional)
!> @param[OUT] depth: Water depth (non-dimensional)
!> @param[OUT] g: Gravity (non-dimensional)
!> @param[OUT] L: Non-dimensional length value
!> @param[OUT] T: Non-dimensional time value 
!
!> @author Xu Haihua, Technology Centre for Offshore and Marine, Singapore (TCOMS)
!> @date 11 Nov 2019 
!> @remarks None
!-------------------------------------------------------------
SUBROUTINE init_read_mod(filename,n1,n2,dt_out,T_stop,xlen,ylen,depth,g,L,T)
!
! Initialize data from volumic mode description generated by HOS-ocean
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(IN) :: filename

!
INTEGER, INTENT(OUT)         :: n1,n2
REAL(RP), INTENT(OUT)        :: dt_out,T_stop,xlen,ylen,depth,g,L,T
!
INTEGER        :: i_unit

! Local variables
REAL(RP) :: x1, x2
!
! We will look at first eight variables written on 18 characters
OPEN(NEWUNIT=i_unit,file=filename,status='OLD', FORM='FORMATTED', ACCESS='DIRECT',RECL=18*10)
READ(i_unit,'(10(ES17.10,1X))',REC=1) x1, x2, dt_out, T_stop, xlen, ylen, depth, g, L, T
!
n1 = NINT(x1)
n2 = NINT(x2)
!
CLOSE(i_unit)


!
END SUBROUTINE init_read_mod
!
!
!
!---READ HOS MODES from the file
!-------------------------------------------------------------
!SUBROUTINE         read_mod
!-------------------------------------------------------------
!> @brief Read HOS modes from the HOS SWESEN file 
!> @details Read HOS modes from the HOS SWESEN file 
!> @param[in] fileName: The HOS SWESEN file name
!> @param[in] time: the time to read (dimensional time)
!
!> @author Xu Haihua, Technology Centre for Offshore and Marine, Singapore (TCOMS)
!> @date 11 Nov 2019 
!> @remarks None
!-------------------------------------------------------------
SUBROUTINE read_mod(filename,time)
!
! Initialize data from volumic mode description generated by HOS-ocean
!
!
IMPLICIT NONE
!---------Arguments----------------------
CHARACTER(LEN=*), INTENT(IN) :: filename
REAL(RP), INTENT(IN)         :: time

!---------Local Variables----------------------
!COMPLEX(CP), INTENT(OUT), DIMENSION(n1o2p1,n2) :: modesspecx,modesspecy,modesspecz,modesspect,modesFS,modesFSt
!
INTEGER                      :: i_unit
! Local variables
INTEGER                      :: i1, i2, it,ix

REAL(RP)                     :: time_star,dt_out
REAL(RP)                     :: mode_abs(n1o2p1,n2)
REAL(RP)                     :: rATotal,rASelect
CHARACTER(100)               :: cBuff
INTEGER,SAVE                 :: iread=1;
!
! We read the specific records corresponding to time
!
!it = NINT(time/dt_out)+1


modesspecx_a(:,:,1)  =  modesspecx_a(:,:,2)
modesspecy_a(:,:,1)  =  modesspecy_a(:,:,2)
modesspecz_a(:,:,1)  =  modesspecz_a(:,:,2)
modesspect_a(:,:,1)  =  modesspect_a(:,:,2)

modesFS_a(:,:,1)   =  modesFS_a(:,:,2)
modesFSt_a(:,:,1)  =  modesFSt_a(:,:,2);

anglex_a(:,:,1)  =  anglex_a(:,:,2)
angley_a(:,:,1)  =  angley_a(:,:,2)
anglez_a(:,:,1)  =  anglez_a(:,:,2)
anglet_a(:,:,1)  =  anglet_a(:,:,2)

angleut_a(:,:,1)  =  angleut_a(:,:,2)
anglevt_a(:,:,1)  =  anglevt_a(:,:,2)
anglewt_a(:,:,1)  =  anglewt_a(:,:,2)

!--assgin working array
CALL assign_pointer_arrays(2)

time_star = time / T;
it = NINT(time_star/dt_out_star)+1
!
OPEN(NEWUNIT=i_unit,file=filename,status='OLD', FORM='FORMATTED', ACCESS='DIRECT',RECL=18*(2*n1o2p1))
!
DO i2=1,n2
    READ(i_unit,'(5000(ES17.10,1X))',REC=((it)*n2*6)+1+6*(i2-1)) (modesspecx(i1,i2), i1=1,n1o2p1)
    READ(i_unit,'(5000(ES17.10,1X))',REC=((it)*n2*6)+2+6*(i2-1)) (modesspecy(i1,i2), i1=1,n1o2p1)
    READ(i_unit,'(5000(ES17.10,1X))',REC=((it)*n2*6)+3+6*(i2-1)) (modesspecz(i1,i2), i1=1,n1o2p1)
    READ(i_unit,'(5000(ES17.10,1X))',REC=((it)*n2*6)+4+6*(i2-1)) (modesspect(i1,i2), i1=1,n1o2p1)
    READ(i_unit,'(5000(ES17.10,1X))',REC=((it)*n2*6)+5+6*(i2-1)) (modesFS(i1,i2)   , i1=1,n1o2p1)
    READ(i_unit,'(5000(ES17.10,1X))',REC=((it)*n2*6)+6+6*(i2-1)) (modesFSt(i1,i2)  , i1=1,n1o2p1)
ENDDO
!

!ENDDO
DO i1=1,n1o2p1
    DO i2=1,n2
        anglex(i1,i2)  = ATAN2(AIMAG(modesspecx(i1,i2)),REAL(modesspecx(i1,i2),RP))
        angley(i1,i2)  = ATAN2(AIMAG(modesspecy(i1,i2)),REAL(modesspecy(i1,i2),RP))
        anglez(i1,i2)  = ATAN2(AIMAG(modesspecz(i1,i2)),REAL(modesspecz(i1,i2),RP))
        anglet(i1,i2)  = ATAN2(AIMAG(modesspect(i1,i2)),REAL(modesspect(i1,i2),RP))
        angleut(i1,i2) = ATAN2(AIMAG(ikx(i1,i2)*modesspect(i1,i2)),REAL(ikx(i1,i2)*modesspect(i1,i2),RP))
        anglevt(i1,i2) = ATAN2(AIMAG(iky(i1,i2)*modesspect(i1,i2)),REAL(iky(i1,i2)*modesspect(i1,i2),RP))
        anglewt(i1,i2) = ATAN2(AIMAG(kth(i1,i2)*modesspect(i1,i2)),REAL(kth(i1,i2)*modesspect(i1,i2),RP))
    ENDDO
ENDDO

CLOSE(i_unit)
!

!--HOS reconstruction only use selected mode,then
!--order the modes by surface spectrum and get the index.
!IF(iHOSModeMethod ==2 .AND. lHOSModeSelected ==.FALSE.) THEN
IF(iHOSModeMethod ==2 .AND. MOD(iread,10)== 1 ) THEN
! IF (lHOSModeSelected ==.FALSE.) THEN
!    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i2)

    DO i1=1,n1o2p1
        DO i2=1,n2
            mode_abs(i1,i2) = ABS(modesFS(i1,i2)) !only consider the y=1
        ENDDO
    ENDDO
!    !$OMP END PARALLEL DO  

    
!sort the index for each row and the constan mode (first mode) is always in.    
    DO i2=1,n2
        iModeUsedIndex(1,i2)=0;
!--sort the mode from second model to the last one. the firt mode is constant mode and must contain        
        CALL sortBubble(n1o2p1-1, mode_abs(2:n1o2p1,i2), iModeUsedIndex(2:n1o2p1,i2), -1)
        
        iModeUsedIndex(:,i2) = iModeUsedIndex(:,i2) + 1
        iModeUsedIndex(1,i2) = 1;
        
    ENDDO

!--compute the total amplitude in HOS    
    rATotal = 0.0d0
    DO i1=1,n1o2p1
        DO i2=1,n2
            rATotal = rATotal + mode_abs(i1,i2)
        ENDDO
    ENDDO
    
!--compute the selecte mode amplitude in HOSs    
    rASelect =0.0
    DO i2=1,n2
        DO ix =1,nHOSModeUsed !loop the selected modes.
            i1 = iModeUsedIndex(ix,i2);  !get the modes index
!            get_pt_eta= get_pt_eta  + 1.0_rp*REAL(modesFS(i1,i2) * EXP(i*kx(i1)*xx)*EXP(i*ky_n2(i2)*yy))
            rASelect = rASelect + mode_abs(i1,i2)
        ENDDO
    ENDDO 
    
    WRITE(*,"(A,F10.6,A)") "The total energy selected modes contains : " , rASelect/rATotal * 100, "%";
    
    lHOSModeSelected =.TRUE.
!ENDIF
ENDIF
END SUBROUTINE read_mod



!----assgin pointer to the working array.
!-------------------------------------------------------------
!SUBROUTINE         assign_pointer_arrays
!-------------------------------------------------------------
!> @brief Assign working array (modesspecx ... modesFSt) 
!> @details Assign working array (modesspecx ... modesFSt) to the storage array.  
!> @param[in] it: the array index
!
!> @author Xu Haihua, Technology Centre for Offshore and Marine, Singapore (TCOMS)
!> @date 11 Nov 2019 
!> @remarks if it=1, modesspecx ... modesFSt point to modesspecx_a(:,:,1)...modesFSt(:,:,1)
!>          if it=2, modesspecx ... modesFSt point to modesspecx_a(:,:,2)...modesFSt(:,:,2)
!-------------------------------------------------------------

SUBROUTINE assign_pointer_arrays(it)
!---------Arguments----------------------
INTEGER(4),  INTENT(IN)  :: it;
!---------Local Variable----------------------

modesspecx => modesspecx_a(:,:,it)
modesspecy => modesspecy_a(:,:,it)
modesspecz => modesspecz_a(:,:,it)
modesspect => modesspect_a(:,:,it)

modesFS  => modesFS_a(:,:,it)
modesFSt => modesFSt_a(:,:,it);

anglex => anglex_a(:,:,it)
angley => angley_a(:,:,it)
anglez => anglez_a(:,:,it)
anglet => anglet_a(:,:,it)

angleut => angleut_a(:,:,it)
anglevt => anglevt_a(:,:,it)
anglewt => anglewt_a(:,:,it)



END SUBROUTINE

!END MODULE variables_3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    Copyright (C) 2014 - LHEEA Lab., Ecole Centrale de Nantes, UMR CNRS 6598
!
!    This program is part of HOS-ocean
!
!    HOS-ocean is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!





!PRIVATE:: get_pt_eta;
!PRIVATE:: get_pt_uvwp

!IMPLICIT NONE
!CONTAINS
!init the HOS model and return the output time step in the SWENSE file
!-------------------------------------------------------------
!SUBROUTINE recons_HOS_init
!-------------------------------------------------------------
!> @brief init the HOS model and return the output time step in the SWENSE file
!> @details init the HOS model and return the output time step in the SWENSE file
!> @param[in] filename : HOS Ocean SWENSE file name;
!> @param[in] dt_out   : output time step unit is second
!> @param[in] HOS_depth: water depth, unit is m
!> @param[in] HOS_g    : acceleration unit is m/s^2
!
!> @remarks the returned value is non-dimensional
!>
!> @author Xu Haihua 
!> @date 10 March 2020 
!REAL(RP) FUNCTION get_pt_eta(modesFS,xx,yy)
!-------------------------------------------------------------
SUBROUTINE recons_HOS_init(filename,dt_out,HOS_depth,HOS_g)
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(IN) :: filename
REAL(RP),INTENT(OUT)         :: dt_out
REAL(RP),INTENT(OUT)         :: HOS_depth,HOS_g

INTEGER                      :: i_unit,info
!
! Local variables
!INTEGER  :: n1,n2
!REAL(RP), INTENT(OUT) :: dt_out_star,T_stop_star,xlen_star,ylen_star,depth_star,g_star,L,T
!
!COMPLEX(CP), ALLOCATABLE, DIMENSION(:,:), INTENT(OUT) :: modesspecx,modesspecy,modesspecz,modesspect,modesFS,modesFSt
!

HOSFileName = TRIM(filename)   ;

! Initialize variables reading filename
! In the file, everything is non-dimensional with L and T length and time scales
CALL init_read_mod(filename,n1,n2,dt_out_star,T_stop_star,xlen_star,ylen_star,depth_star,g_star,L,T)

lHOSModeSelected = .FALSE.

dt_out  = dt_out_star * T;
T_stop  = T_stop_star * T;
depth   = depth_star  * L;
g       = g_star * L / (T*T)
xlen    = xlen_star * L;
ylen    = ylen_star * L;

!
!n1o2p1 = n1/2+1

n1o2p1 = n1/2+1
n2o2p1 = n2/2+1
    
    
m1      = n1
m2      = n2
Nd1     = 1
Nd2     = 1
Nd1o2p1 = 1
m1o2p1  = n1o2p1
md1o2p1 = 1
md1     = 1
md2     = 1
    
IF(ALLOCATED(modesspecx_a)) DEALLOCATE(modesspecx_a)
IF(ALLOCATED(modesspecy_a)) DEALLOCATE(modesspecy_a)
IF(ALLOCATED(modesspecz_a)) DEALLOCATE(modesspecz_a)
IF(ALLOCATED(modesspect_a)) DEALLOCATE(modesspect_a)
IF(ALLOCATED(modesFS_a))    DEALLOCATE(modesFS_a)
IF(ALLOCATED(modesFSt_a))   DEALLOCATE(modesFSt_a)

ALLOCATE(modesspecx_a(n1o2p1,n2,ntn),modesspecy_a(n1o2p1,n2,ntn),modesspecz_a(n1o2p1,n2,ntn),modesspect_a(n1o2p1,n2,ntn),modesFS_a(n1o2p1,n2,ntn), &
         modesFSt_a(n1o2p1,n2,ntn),STAT=info);
IF(info /=0) THEN
    WRITE(*,*)"Failed to allocate memory modesspecx_a, modesspecy_a, ... , modesFSt_a"
ENDIF

IF(ALLOCATED(anglex_a)) DEALLOCATE(anglex_a)
IF(ALLOCATED(angley_a)) DEALLOCATE(angley_a)
IF(ALLOCATED(anglez_a)) DEALLOCATE(anglez_a)
IF(ALLOCATED(anglet_a)) DEALLOCATE(anglet_a)
IF(ALLOCATED(angleut_a)) DEALLOCATE(angleut_a)
IF(ALLOCATED(anglevt_a)) DEALLOCATE(anglevt_a)
IF(ALLOCATED(anglewt_a)) DEALLOCATE(anglewt_a)
        
ALLOCATE(anglex_a(m1o2p1,m2,ntn),angley_a(m1o2p1,m2,ntn),anglez_a(m1o2p1,m2,ntn),anglet_a(m1o2p1,m2,ntn),angleut_a(m1o2p1,m2,ntn), &
         anglevt_a(m1o2p1,m2,ntn) , anglewt_a(m1o2p1,m2,ntn),STAT=info);
IF(info /=0) THEN
    WRITE(*,*)"Failed to allocate memory anglex_a, angley_a, ... , anglewt_a"
ENDIF

IF(iHOSModeMethod ==2) THEN
    IF(ALLOCATED(iModeUsedIndex)) DEALLOCATE(iModeUsedIndex)
  
    ALLOCATE(iModeUsedIndex(n1o2p1,n2),STAT=info);
    IF(info /=0) THEN
        WRITE(*,*)"Failed to allocate memory iModeUsedIndexNx, iModeUsedIndexNy"
    ENDIF
!---limit the maximum modes select to n1 or n2
    nHOSModeUsed = MIN(nHOSModeUsed,n1o2p1);
 
ENDIF

!IF(ALLOCATED(x)) DEALLOCATE(x)
!IF(ALLOCATED(y)) DEALLOCATE(y)
!IF(ALLOCATED(kx)) DEALLOCATE(kx)
!IF(ALLOCATED(ky_n2)) DEALLOCATE(ky_n2)
!IF(ALLOCATED(ikx)) DEALLOCATE(ikx)
!IF(ALLOCATED(iky)) DEALLOCATE(iky)
!IF(ALLOCATED(kth)) DEALLOCATE(kth)
!
!ALLOCATE(x(n1),y(n2),kx(n1o2p1),ky_n2(n2),ikx(n1o2p1,n2),iky(n1o2p1,n2),kth(n1o2p1,n2))
    
CALL build_mesh_global()
!
! Read time=0
CALL read_mod(filename,0.0_rp)
!
HOS_depth = depth
HOS_g     = g;

END SUBROUTINE recons_HOS_init
!


!
!SUBROUTINE build_mesh_global(xlen_star,ylen_star,depth_star,n1,n2,x,y,kx,ky_n2,ikx,iky,kth)
SUBROUTINE build_mesh_global()
!
IMPLICIT NONE
!
!REAL(RP), INTENT(IN) :: xlen_star, ylen_star,depth_star
!INTEGER, INTENT(IN)  :: n1, n2
!
!REAL(RP), DIMENSION(n1), INTENT(OUT)           :: x
!REAL(RP), DIMENSION(n1o2p1), INTENT(OUT)       :: kx
!REAL(RP), DIMENSION(n2), INTENT(OUT)           :: y, ky_n2
!REAL(RP), DIMENSION(n1o2p1,n2), INTENT(OUT)    :: kth
!COMPLEX(CP), DIMENSION(n1o2p1,n2), INTENT(OUT) :: ikx, iky
! Local variables
REAL(RP) :: pioxlen, pioylen, delx, dely, k2
INTEGER  :: N_der(2)
INTEGER  :: i1,i2

IF(ALLOCATED(x)) DEALLOCATE(x)
IF(ALLOCATED(y)) DEALLOCATE(y)
IF(ALLOCATED(kx)) DEALLOCATE(kx)
IF(ALLOCATED(ky_n2)) DEALLOCATE(ky_n2)
IF(ALLOCATED(ikx)) DEALLOCATE(ikx)
IF(ALLOCATED(iky)) DEALLOCATE(iky)
IF(ALLOCATED(kth)) DEALLOCATE(kth)
ALLOCATE(x(n1),y(n2),kx(n1o2p1),ky_n2(n2),ikx(n1o2p1,n2),iky(n1o2p1,n2),kth(n1o2p1,n2))
!
! Specify temporary number of points
N_der(1) = n1o2p1
N_der(2) = n2o2p1
!
! Specify length of domain
!
pioxlen = TWOPI / xlen_star
!
IF (n2 == 1) THEN
    pioylen = 0.0_rp
ELSE
    pioylen = TWOPI / ylen_star
ENDIF
!
!   mesh generation
!
delx = xlen_star / n1
DO i1 = 1,n1
    x(i1) = (i1 - 1) * delx
ENDDO
!
IF (n2 == 1) THEN
    dely = 0.0_rp
ELSE
    dely = ylen_star / n2
ENDIF
DO i2 = 1,n2
    y(i2) = (i2 - 1) * dely
ENDDO
!
!   wave numbers
DO i1 = 1, n1o2p1
    kx(i1)  = REAL(i1 - 1,RP) * pioxlen
ENDDO
!  y-wave numbers (on n2 modes)
DO i2 = 1, n2o2p1
    ky_n2(i2) = REAL(i2 - 1,RP) * pioylen
ENDDO
DO i2 = 2,n2o2p1
    ky_n2(n2-i2+2) = - REAL(i2 - 1,RP) * pioylen
ENDDO
!
IF (iseven(n2)) ky_n2(n2o2p1) = REAL(n2o2p1 - 1,RP) * pioylen
! Storage for derivatives
ikx = 0.0_cp
!  x-derivative on n1 points (i.e. n1o2p1 modes)
DO i2 = 1, n2o2p1
    ikx(1:MIN(N_der(1),n1o2p1),i2) = i * kx(1:MIN(N_der(1),n1o2p1))
ENDDO
! Last mode contains cos information only and must not be part of the differentiation.
IF (iseven(n1)) ikx(n1o2p1,:) = 0.0_cp
! negative ky
DO i2 = 2, n2o2p1
    ikx(1:MIN(N_der(1),n1o2p1),n2-i2+2) = i * kx(1:MIN(N_der(1),n1o2p1))
ENDDO
!
iky = 0.0_cp
! y-derivative on n1 points (i.e. n1o2p1 modes)
DO i1 = 1, n1o2p1
    iky(i1,1:MIN(N_der(2),n2o2p1)) = i * ky_n2(1:MIN(N_der(2),n2o2p1))
    ! negative ky
    DO i2 = 2, MIN(N_der(2), n2o2p1)
        iky(i1,n2-i2+2) = - i * ky_n2(i2)
    ENDDO
    IF (iseven(n2) .AND. N_der(2)>=n2o2p1) iky(i1,n2o2p1) = 0.0_cp
ENDDO
!
! HOS modal coefficients of the vertical derivatives
DO i2 = 1, n2
    DO i1 = 1, n1o2p1
      k2         = kx(i1) * kx(i1) + ky_n2(i2) * ky_n2(i2)
      kth(i1,i2) = SQRT(k2)*TANH(SQRT(k2) * depth_star)
    ENDDO
ENDDO
!
END SUBROUTINE build_mesh_global


!-------------------------------------------------------------
!SUBROUTINE get_pt_eta
!-------------------------------------------------------------
!> @brief Get dimensionalal surface elevation at one point (x,y) in horizontal 2D
!> @details Get dimensionalal surface elevation at one point (x,y)  in horizontal 2D
!> @param[in] it       : time index =1: get eta at t0;
!                                   =2: get eta at t1;
!> @param[in] xx       : dimensional coordinate in x direction.
!> @param[in] yy       : dimensional coordinate in y direction.
!> @return             : dimensional surface elevation
!
!> @remarks the returned value is non-dimensional
!>
!> @author Xu Haihua 
!> @date 10 March 2020 
!REAL(RP) FUNCTION get_pt_eta(modesFS,xx,yy)
REAL(RP) FUNCTION get_pt_eta(it,xxl,yyl)
IMPLICIT NONE
!COMPLEX(CP),  DIMENSION(:,:),INTENT(IN)            :: modesFS
INTEGER,INTENT(IN)    :: it;            
REAL(RP),INTENT(IN)   :: xxl,yyl; 

!---Local variables    
INTEGER             :: i1,i2,ix; 
REAL(RP)            :: xx,yy,rtmp; 

CALL assign_pointer_arrays(it);

xx = xxl/L
yy = yyl/L


!use all the modes in the results.
IF(iHOSModeMethod ==1) THEN
    get_pt_eta = 0.0_RP
    DO i2=1,n2 
        DO i1=1,n1o2p1
!            rtmp =  1.0_rp*REAL(modesFS(i1,i2) * EXP(i*kx(i1)*xx)*EXP(i*ky_n2(i2)*yy));
            get_pt_eta= get_pt_eta  + 1.0_rp*REAL(modesFS(i1,i2) * EXP(i*kx(i1)*xx)*EXP(i*ky_n2(i2)*yy))
        ENDDO
    ENDDO
   
ELSEIF(iHOSModeMethod ==2) THEN
!IF(iHOSModeMethod ==2) THEN
    get_pt_eta = 0.0_RP
    DO i2=1,n2
        DO ix =1,nHOSModeUsed !loop the selected modes.
            i1 = iModeUsedIndex(ix,i2);  !get the modes index
            get_pt_eta= get_pt_eta  + 1.0_rp*REAL(modesFS(i1,i2) * EXP(i*kx(i1)*xx)*EXP(i*ky_n2(i2)*yy))
        ENDDO
    ENDDO 

ENDIF

 
!get_pt_eta = get_pt_eta * L
get_pt_eta = rPFResultScale * get_pt_eta * L              

END FUNCTION


!-------------------------------------------------------------
!SUBROUTINE get_pt_deta_dt
!-------------------------------------------------------------
!> @brief Get dimensionalal surface elevation time derivative at one point (x,y) in horizontal 2D
!> @details Get dimensionalal surface elevation time derivative at one point (x,y)  in horizontal 2D
!> @param[in] it       : time index =1: get eta at t0;
!                                   =2: get eta at t1;
!> @param[in] xx       : dimensional coordinate in x direction.
!> @param[in] yy       : dimensional coordinate in y direction.
!> @return             : dimensional surface elevation
!
!> @remarks the returned value is non-dimensional
!>
!> @author Xu Haihua 
!> @date 10 March 2020 
!REAL(RP) FUNCTION get_pt_eta(modesFS,xx,yy)
REAL(RP) FUNCTION get_pt_deta_dt(it,xxl,yyl)
IMPLICIT NONE
!COMPLEX(CP),  DIMENSION(:,:),INTENT(IN)            :: modesFS
INTEGER,INTENT(IN)    :: it;            
REAL(RP),INTENT(IN)   :: xxl,yyl; 

!---Local variables    
INTEGER             :: i1,i2,ix; 
REAL(RP)            :: xx,yy,rtmp; 

CALL assign_pointer_arrays(it);

xx = xxl/L
yy = yyl/L


!use all the modes in the results.
IF(iHOSModeMethod ==1) THEN
    get_pt_deta_dt = 0.0_RP
    DO i2=1,n2 
        DO i1=1,n1o2p1
!            rtmp =  1.0_rp*REAL(modesFS(i1,i2) * EXP(i*kx(i1)*xx)*EXP(i*ky_n2(i2)*yy));
            get_pt_deta_dt= get_pt_deta_dt  + 1.0_rp*REAL(modesFSt(i1,i2) * EXP(i*kx(i1)*xx)*EXP(i*ky_n2(i2)*yy))
        ENDDO
    ENDDO
   
ELSEIF(iHOSModeMethod ==2) THEN
!IF(iHOSModeMethod ==2) THEN
    get_pt_deta_dt = 0.0_RP
    DO i2=1,n2
        DO ix =1,nHOSModeUsed !loop the selected modes.
            i1 = iModeUsedIndex(ix,i2);  !get the modes index
            get_pt_deta_dt= get_pt_deta_dt  + 1.0_rp*REAL(modesFSt(i1,i2) * EXP(i*kx(i1)*xx)*EXP(i*ky_n2(i2)*yy))
        ENDDO
    ENDDO 

ENDIF

 
!get_pt_eta = get_pt_eta * L
get_pt_deta_dt = rPFResultScale * get_pt_deta_dt * L / T             

END FUNCTION


!
!!-------------------------------------------------------------
!!SUBROUTINE get_pt_eta
!!-------------------------------------------------------------
!!> @brief Get dimensionalal surface elevation at one point (x,y) in horizontal 2D
!!> @details Get dimensionalal surface elevation at one point (x,y)  in horizontal 2D
!!> @param[in] it       : time index =1: get eta at t0;
!!                                   =2: get eta at t1;
!!> @param[in] xx       : dimensional coordinate in x direction.
!!> @param[in] yy       : dimensional coordinate in y direction.
!!> @return             : dimensional surface elevation
!!
!!> @remarks the returned value is non-dimensional
!!>
!!> @author Xu Haihua 
!!> @date 10 March 2020 
!!REAL(RP) FUNCTION get_pt_eta(modesFS,xx,yy)
!REAL(RP) FUNCTION get_pt_phi(it,xxl,yyl,zzl)
!IMPLICIT NONE
!!COMPLEX(CP),  DIMENSION(:,:),INTENT(IN)            :: modesFS
!INTEGER,INTENT(IN)    :: it;            
!REAL(RP),INTENT(IN)   :: xxl,yyl,zzl; 
!
!!---Local variables    
!INTEGER             :: i1,i2,ix; 
!REAL(RP)            :: xx,yy,zz,rtmp; 
!COMPLEX(CP), DIMENSION(n1o2p1,n2):: modesPhi
!REAL(RP)                         :: xvect,yvect,zvect;    !non dimensional
!!-------Local variables
!
!REAL(RP)                        :: phi; !phit,dudt,dvdt,dwdt;
!COMPLEX(CP)                     :: phi_1; !vitx_l, vity_l, vitz_l, phit_l, dudt_l, dvdt_l, dwdt_l
!COMPLEX(CP)                     :: coeff
!!INTEGER                         :: nTotal,iTotal
!REAL(RP) :: k_n2
!
!
!CALL assign_pointer_arrays(it);
!
!xvect = xxl/L
!yvect = yyl/L
!zvect = zzl/L
!
!i1 =1
!i2 =1
!modesPhi(i1,i2) =  modesFS(i1,i2)
!DO i2=1,n2 
!    DO i1=2,n1o2p1
!        modesPhi(i1,i2) =  modesFS(i1,i2) / (ikx(i1,i2))   
!    ENDDO
!ENDDO
!
!i1 = 1
!i2 = 1
!
!phi = 0.0d0
!
!DO i2=1,n2 
!    DO i1=2,n1o2p1
!    k_n2 = SQRT(kx(i1)**2+ky_n2(i2)**2)
!    IF ((k_n2*(zvect+depth_star).LT.50.).AND.(k_n2*depth_star.LT.50.)) THEN
!        coeff = COSH(k_n2*(zvect+depth_star))/COSH(k_n2*depth_star) !* EXP(i*ky_n2(i2)*y(ii2))
!!        coeff2= SINH(k_n2*(zvect+depth_star))/SINH(k_n2*depth_star) !* EXP(i*ky_n2(i2)*y(ii2))
!    ELSE
!        coeff = EXP(k_n2*zvect)
!!        coeff2= coeff
!    ENDIF
!    phi  = phi  + 2.0_rp*REAL(modesPhi(i1,i2) * coeff * EXP(i*kx(i1)*xvect)*EXP(i*ky_n2(i2)*yvect))
!    !vitx = vitx + 2.0_rp*ABS(modesspecx(i1,i2) * coeff) *COS(ky_n2(i2)*yvect+anglex(i1,i2))
!    !vity = vity + 2.0_rp*ABS(modesspecy(i1,i2) * coeff) *COS(ky_n2(i2)*yvect+angley(i1,i2))
!    !vitz = vitz + 2.0_rp*ABS(modesspecz(i1,i2) * coeff2)*COS(ky_n2(i2)*yvect+anglez(i1,i2))
!    !phit = phit + 2.0_rp*ABS(modesspect(i1,i2) * coeff) *COS(ky_n2(i2)*yvect+anglet(i1,i2))
!    !dudt= dudt+ 2.0_rp*ABS(ikx(i1,i2)*modesspect(i1,i2) * coeff) &
!    !    *COS(ky_n2(i2)*yvect+angleut(i1,i2))
!    !dvdt = dvdt + 2.0_rp*ABS(iky(i1,i2)*modesspect(i1,i2) * coeff) &
!    !    *COS(ky_n2(i2)*yvect+anglevt(i1,i2))
!    !dwdt = dwdt + 2.0_rp*ABS(kth(i1,i2)*modesspect(i1,i2) * coeff2) &
!    !    *COS(ky_n2(i2)*yvect+anglewt(i1,i2))
!    ENDDO
!ENDDO
!
! 
!!get_pt_eta = get_pt_eta * L
!get_pt_phi = rPFResultScale * phi * L* L / T              
!
!END FUNCTION

!-------------------------------------------------------------
!SUBROUTINE get_pt_uvwp
!-------------------------------------------------------------
!> @brief Get dimensional velocity and pressure and one 3D point
!> @details Get dimensional surface elevation at one point (x,y)
!> @param[in] it       : time index =1: get eta at t0;
!                                   =2: get eta at t1;
!> @param[in] xvectl       : dimensional coordinate in x direction.
!> @param[in] yvectl       : dimensional coordinate in y direction.
!> @param[in] zvectl       : dimensional coordinate in y direction.
!> @param[out] vitx       : dimensional velocity in x direction.
!> @param[out] vity       : dimensional velocity in y direction.
!> @param[out] vity       : dimensional velocity in y direction.
!> @param[out] pit        : dimensional total pressure .
!
!> @remarks the returned value is dimensional
!>
!> @author Xu Haihua 
!> @date 10 March 2020 
!SUBROUTINE get_pt_uvwp(modesspecx,modesspecy,modesspecz,modesspect,depth_star,g_star,xvect,yvect,zvect,vitx,vity,vitz,pit)
SUBROUTINE get_pt_uvwp(it,xvectl,yvectl,zvectl,vitx,vity,vitz,pit)
!
IMPLICIT NONE
!% INPUT VARIABLES
!COMPLEX(CP), DIMENSION(m1o2p1,m2), INTENT(IN) :: modesspecx,modesspecy,modesspecz,modesspect
!REAL(RP)                         , INTENT(IN) :: depth_star,g_star
!% INPUT VARIABLES
INTEGER(4),  INTENT(IN) :: it;
REAL(RP)   , INTENT(IN) :: xvectl,yvectl,zvectl;    !non dimensional
REAL(RP)   , INTENT(OUT):: vitx,vity,vitz,pit
!----------------------------------------------------------

INTEGER                        :: i1,i2,ii1,ii2,ix
!
REAL(RP)                       :: xvect,yvect,zvect;    !non dimensional
!-------Local variables

REAL(RP)                        :: phit,dudt,dvdt,dwdt;
COMPLEX(CP)                     :: vitx_l, vity_l, vitz_l, phit_l, dudt_l, dvdt_l, dwdt_l
COMPLEX(CP)                     :: coeff, coeff2
!INTEGER                         :: nTotal,iTotal
REAL(RP) :: k_n2

CALL assign_pointer_arrays(it);

xvect = xvectl / L 
yvect = yvectl / L
zvect = zvectl / L
!
! constant mode
i1 = 1
i2 = 1
!
vitx = REAL(modesspecx(i1,i2),RP)
vity = REAL(modesspecy(i1,i2),RP)
vitz = REAL(modesspecz(i1,i2),RP)
phit = REAL(modesspect(i1,i2),RP)
dudt = REAL(ikx(i1,i2)*modesspect(i1,i2),RP)
dvdt = REAL(iky(i1,i2)*modesspect(i1,i2),RP)
dwdt = REAL(kth(i1,i2)*modesspect(i1,i2),RP)
! i1=1 and all i2
DO i2=2,n2o2p1
    k_n2 = SQRT(kx(i1)**2+ky_n2(i2)**2)
    IF ((k_n2*(zvect+depth_star).LT.50.).AND.(k_n2*depth_star.LT.50.)) THEN
        coeff = COSH(k_n2*(zvect+depth_star))/COSH(k_n2*depth_star) !* EXP(i*ky_n2(i2)*y(ii2))
        coeff2= SINH(k_n2*(zvect+depth_star))/SINH(k_n2*depth_star) !* EXP(i*ky_n2(i2)*y(ii2))
    ELSE
        coeff = EXP(k_n2*zvect)
        coeff2= coeff
    ENDIF
    vitx = vitx + 2.0_rp*ABS(modesspecx(i1,i2) * coeff) *COS(ky_n2(i2)*yvect+anglex(i1,i2))
    vity = vity + 2.0_rp*ABS(modesspecy(i1,i2) * coeff) *COS(ky_n2(i2)*yvect+angley(i1,i2))
    vitz = vitz + 2.0_rp*ABS(modesspecz(i1,i2) * coeff2)*COS(ky_n2(i2)*yvect+anglez(i1,i2))
    phit = phit + 2.0_rp*ABS(modesspect(i1,i2) * coeff) *COS(ky_n2(i2)*yvect+anglet(i1,i2))
    dudt= dudt+ 2.0_rp*ABS(ikx(i1,i2)*modesspect(i1,i2) * coeff) &
        *COS(ky_n2(i2)*yvect+angleut(i1,i2))
    dvdt = dvdt + 2.0_rp*ABS(iky(i1,i2)*modesspect(i1,i2) * coeff) &
        *COS(ky_n2(i2)*yvect+anglevt(i1,i2))
    dwdt = dwdt + 2.0_rp*ABS(kth(i1,i2)*modesspect(i1,i2) * coeff2) &
        *COS(ky_n2(i2)*yvect+anglewt(i1,i2))
ENDDO
! FIXME: add the case n2 even
IF (iseven(n2)) THEN
    i1=1
    i2=n2o2p1
    vitx = vitx - 1.0_rp*ABS(modesspecx(i1,i2) * coeff) *COS(ky_n2(i2)*yvect+anglex(i1,i2))
    vity = vity - 1.0_rp*ABS(modesspecy(i1,i2) * coeff) *COS(ky_n2(i2)*yvect+angley(i1,i2))
    vitz = vitz - 1.0_rp*ABS(modesspecz(i1,i2) * coeff2)*COS(ky_n2(i2)*yvect+anglez(i1,i2))
    phit = phit - 1.0_rp*ABS(modesspect(i1,i2) * coeff) *COS(ky_n2(i2)*yvect+anglet(i1,i2))
    dudt= dudt- 1.0_rp*ABS(ikx(i1,i2)*modesspect(i1,i2) * coeff) &
        *COS(ky_n2(i2)*yvect+angleut(i1,i2))
    dvdt = dvdt - 1.0_rp*ABS(iky(i1,i2)*modesspect(i1,i2) * coeff) &
        *COS(ky_n2(i2)*yvect+anglevt(i1,i2))
    dwdt = dwdt - 1.0_rp*ABS(kth(i1,i2)*modesspect(i1,i2) * coeff2) &
        *COS(ky_n2(i2)*yvect+anglewt(i1,i2))
ENDIF
! i2 and i1 =/ 1

IF(iHOSModeMethod ==1) THEN

    DO i1=2,n1o2p1
        DO i2=1,n2
            k_n2 = SQRT(kx(i1)**2+ky_n2(i2)**2)
            IF ((k_n2*(zvect+depth_star).LT.50.).AND.(k_n2*depth_star.LT.50.)) THEN
                coeff = COSH(k_n2*(zvect+depth_star))/COSH(k_n2*depth_star) * EXP(i*ky_n2(i2)*yvect)
                coeff2= SINH(k_n2*(zvect+depth_star))/SINH(k_n2*depth_star) * EXP(i*ky_n2(i2)*yvect)
            ELSE
                coeff = EXP(k_n2*zvect) * EXP(i*ky_n2(i2)*yvect)
                coeff2= coeff
            ENDIF
            !
            vitx_l = modesspecx(i1,i2) * coeff
            vity_l = modesspecy(i1,i2) * coeff
            vitz_l = modesspecz(i1,i2) * coeff2
            phit_l = modesspect(i1,i2) * coeff
            dudt_l = ikx(i1,i2)*phit_l
            dvdt_l = iky(i1,i2)*phit_l
            dwdt_l = kth(i1,i2)*modesspect(i1,i2) * coeff2
            !
            vitx = vitx + 1.0_rp*ABS(vitx_l)*COS(kx(i1)*xvect+ATAN2(AIMAG(vitx_l),REAL(vitx_l,RP)))
            vity = vity + 1.0_rp*ABS(vity_l)*COS(kx(i1)*xvect+ATAN2(AIMAG(vity_l),REAL(vity_l,RP)))
            vitz = vitz + 1.0_rp*ABS(vitz_l)*COS(kx(i1)*xvect+ATAN2(AIMAG(vitz_l),REAL(vitz_l,RP)))
            phit = phit + 1.0_rp*ABS(phit_l)*COS(kx(i1)*xvect+ATAN2(AIMAG(phit_l),REAL(phit_l,RP)))
            dudt = dudt  + 1.0_rp*ABS(dudt_l)*COS(kx(i1)*xvect+ATAN2(AIMAG(dudt_l),REAL(dudt_l,RP)))
            dvdt = dvdt + 1.0_rp*ABS(dvdt_l)*COS(kx(i1)*xvect+ATAN2(AIMAG(dvdt_l),REAL(dvdt_l,RP)))
            dwdt = dwdt + 1.0_rp*ABS(dwdt_l)*COS(kx(i1)*xvect+ATAN2(AIMAG(dwdt_l),REAL(dwdt_l,RP)))
        ENDDO
    ENDDO

ELSEIF(iHOSModeMethod ==2) THEN

       DO i2=1,n2
        DO ix =2,nHOSModeUsed  !only start with the 2nd mode
            i1 = iModeUsedIndex(ix,i2);  !get the modes index
         
            k_n2 = SQRT(kx(i1)**2+ky_n2(i2)**2)
            IF ((k_n2*(zvect+depth_star).LT.50.).AND.(k_n2*depth_star.LT.50.)) THEN
                coeff = COSH(k_n2*(zvect+depth_star))/COSH(k_n2*depth_star) * EXP(i*ky_n2(i2)*yvect)
                coeff2= SINH(k_n2*(zvect+depth_star))/SINH(k_n2*depth_star) * EXP(i*ky_n2(i2)*yvect)
            ELSE
                coeff = EXP(k_n2*zvect) * EXP(i*ky_n2(i2)*yvect)
                coeff2= coeff
            ENDIF
            !
            vitx_l = modesspecx(i1,i2) * coeff
            vity_l = modesspecy(i1,i2) * coeff
            vitz_l = modesspecz(i1,i2) * coeff2
            phit_l = modesspect(i1,i2) * coeff
            dudt_l = ikx(i1,i2)*phit_l
            dvdt_l = iky(i1,i2)*phit_l
            dwdt_l = kth(i1,i2)*modesspect(i1,i2) * coeff2
            !
            vitx = vitx + 1.0_rp*ABS(vitx_l)*COS(kx(i1)*xvect+ATAN2(AIMAG(vitx_l),REAL(vitx_l,RP)))
            vity = vity + 1.0_rp*ABS(vity_l)*COS(kx(i1)*xvect+ATAN2(AIMAG(vity_l),REAL(vity_l,RP)))
            vitz = vitz + 1.0_rp*ABS(vitz_l)*COS(kx(i1)*xvect+ATAN2(AIMAG(vitz_l),REAL(vitz_l,RP)))
            phit = phit + 1.0_rp*ABS(phit_l)*COS(kx(i1)*xvect+ATAN2(AIMAG(phit_l),REAL(phit_l,RP)))
            dudt = dudt  + 1.0_rp*ABS(dudt_l)*COS(kx(i1)*xvect+ATAN2(AIMAG(dudt_l),REAL(dudt_l,RP)))
            dvdt = dvdt + 1.0_rp*ABS(dvdt_l)*COS(kx(i1)*xvect+ATAN2(AIMAG(dvdt_l),REAL(dvdt_l,RP)))
            dwdt = dwdt + 1.0_rp*ABS(dwdt_l)*COS(kx(i1)*xvect+ATAN2(AIMAG(dwdt_l),REAL(dwdt_l,RP)))
        ENDDO
    ENDDO

ENDIF


pit = - g_star*zvect - 0.5_rp*(vitx**2+vity**2+vitz**2)-phit;


!vitx,vity,vitz,pit

!vitx =vitx * L / T;    !velocity U convert to SI unit
!vity =vity * L / T;    !velocity V convert to SI unit
!vitz =vitz * L / T;    !velocity W convert to SI uint 

!pit = density * pit *L**2/T**2

!rFPResultScale
vitx =rPFResultScale * vitx * L / T;    !velocity U convert to SI unit
vity =rPFResultScale * vity * L / T;    !velocity V convert to SI unit
vitz =rPFResultScale * vitz * L / T;    !velocity W convert to SI uint 

pit =rPFResultScale * density * pit *L**2/T**2

!
!compute z_star

!zvect(ii1,ii2)=zmin+(-zmin+eta_tmp(ii1+imin-1,ii2+jmin-1))*SIN(pio2*REAL(ii3-1,RP)/REAL(i_zvect-1,RP))


!-----conver the results to SI UNIT
!ee = eta_star * L  ;    !z coordinate, convert to SI unit

!uu =vitx * L / T;    !velocity U convert to SI unit
!vv =vity * L / T;    !velocity V convert to SI unit
!ww =vitz * L / T;    !velocity W convert to SI uint 
!pit = - g_star*zvect - 0.5_rp*(vitx**2+vity**2+vitz**2)-phit;
!PP = - g_star*zvect - 0.5_rp*(vitx**2+vity**2+vitz**2)-phit;
!pp = density*PP *L**2/T**2
END SUBROUTINE

!-------------------------------------------------------------
!SUBROUTINE get_pt_uvwp_omp
!-------------------------------------------------------------
!> @brief Get dimensional velocity and pressure and one 3D point OpenMP version
!> @details Get dimensional surface elevation at one point (x,y) OpenMP version
!> @param[in] it       : time index =1: get eta at t0;
!                                   =2: get eta at t1;
!> @param[in] xvectl       : dimensional coordinate in x direction.
!> @param[in] yvectl       : dimensional coordinate in y direction.
!> @param[in] zvectl       : dimensional coordinate in y direction.
!> @param[out] vitx       : dimensional velocity in x direction.
!> @param[out] vity       : dimensional velocity in y direction.
!> @param[out] vity       : dimensional velocity in y direction.
!> @param[out] pit        : dimensional total pressure .
!
!> @remarks the returned value is dimensional
!>
!> @author Xu Haihua 
!> @date 10 March 2020 
!SUBROUTINE get_pt_uvwp(modesspecx,modesspecy,modesspecz,modesspect,depth_star,g_star,xvect,yvect,zvect,vitx,vity,vitz,pit)
SUBROUTINE get_pt_uvwp_omp(it,xvectl,yvectl,zvectl,vitx,vity,vitz,pit)
!
IMPLICIT NONE
!% INPUT VARIABLES
!COMPLEX(CP), DIMENSION(m1o2p1,m2), INTENT(IN) :: modesspecx,modesspecy,modesspecz,modesspect
!REAL(RP)                         , INTENT(IN) :: depth_star,g_star
!% INPUT VARIABLES
INTEGER(4),  INTENT(IN) :: it;
REAL(RP)   , INTENT(IN) :: xvectl,yvectl,zvectl;    !non dimensional
REAL(RP)   , INTENT(OUT):: vitx,vity,vitz,pit
!----------------------------------------------------------

INTEGER                        :: i1,i2,ii1,ii2,ix
!
REAL(RP)                       :: xvect,yvect,zvect;    !non dimensional
!-------Local variables

REAL(RP)                        :: phit,dudt,dvdt,dwdt;
COMPLEX(CP)                     :: vitx_l, vity_l, vitz_l, phit_l, dudt_l, dvdt_l, dwdt_l
COMPLEX(CP)                     :: coeff, coeff2
INTEGER                         :: nTotal,iTotal
REAL(RP) :: k_n2

CALL assign_pointer_arrays(it);

xvect = xvectl / L 
yvect = yvectl / L
zvect = zvectl / L
!
! constant mode
i1 = 1
i2 = 1
!
vitx = REAL(modesspecx(i1,i2),RP)
vity = REAL(modesspecy(i1,i2),RP)
vitz = REAL(modesspecz(i1,i2),RP)
phit = REAL(modesspect(i1,i2),RP)
dudt = REAL(ikx(i1,i2)*modesspect(i1,i2),RP)
dvdt = REAL(iky(i1,i2)*modesspect(i1,i2),RP)
dwdt = REAL(kth(i1,i2)*modesspect(i1,i2),RP)
! i1=1 and all i2

!NOTE: i2 is prive in defalut.
!$OMP PARALLEL DO PRIVATE(k_n2,coeff,coeff2) REDUCTION(+:vitx) REDUCTION(+:vity) REDUCTION(+:vitz) &
!$OMP & REDUCTION(+:phit) REDUCTION(+:dudt) REDUCTION(+:dvdt) REDUCTION(+:dwdt) 
DO i2=2,n2o2p1
    k_n2 = SQRT(kx(i1)**2+ky_n2(i2)**2)
    IF ((k_n2*(zvect+depth_star).LT.50.).AND.(k_n2*depth_star.LT.50.)) THEN
        coeff = COSH(k_n2*(zvect+depth_star))/COSH(k_n2*depth_star) !* EXP(i*ky_n2(i2)*y(ii2))
        coeff2= SINH(k_n2*(zvect+depth_star))/SINH(k_n2*depth_star) !* EXP(i*ky_n2(i2)*y(ii2))
    ELSE
        coeff = EXP(k_n2*zvect)
        coeff2= coeff
    ENDIF
    vitx = vitx + 2.0_rp*ABS(modesspecx(i1,i2) * coeff) *COS(ky_n2(i2)*yvect+anglex(i1,i2))
    vity = vity + 2.0_rp*ABS(modesspecy(i1,i2) * coeff) *COS(ky_n2(i2)*yvect+angley(i1,i2))
    vitz = vitz + 2.0_rp*ABS(modesspecz(i1,i2) * coeff2)*COS(ky_n2(i2)*yvect+anglez(i1,i2))
    phit = phit + 2.0_rp*ABS(modesspect(i1,i2) * coeff) *COS(ky_n2(i2)*yvect+anglet(i1,i2))
    dudt= dudt+ 2.0_rp*ABS(ikx(i1,i2)*modesspect(i1,i2) * coeff) &
        *COS(ky_n2(i2)*yvect+angleut(i1,i2))
    dvdt = dvdt + 2.0_rp*ABS(iky(i1,i2)*modesspect(i1,i2) * coeff) &
        *COS(ky_n2(i2)*yvect+anglevt(i1,i2))
    dwdt = dwdt + 2.0_rp*ABS(kth(i1,i2)*modesspect(i1,i2) * coeff2) &
        *COS(ky_n2(i2)*yvect+anglewt(i1,i2))
ENDDO
!$OMP END PARALLEL DO

! FIXME: add the case n2 even
IF (iseven(n2)) THEN
    i1=1
    i2=n2o2p1
    vitx = vitx - 1.0_rp*ABS(modesspecx(i1,i2) * coeff) *COS(ky_n2(i2)*yvect+anglex(i1,i2))
    vity = vity - 1.0_rp*ABS(modesspecy(i1,i2) * coeff) *COS(ky_n2(i2)*yvect+angley(i1,i2))
    vitz = vitz - 1.0_rp*ABS(modesspecz(i1,i2) * coeff2)*COS(ky_n2(i2)*yvect+anglez(i1,i2))
    phit = phit - 1.0_rp*ABS(modesspect(i1,i2) * coeff) *COS(ky_n2(i2)*yvect+anglet(i1,i2))
    dudt= dudt- 1.0_rp*ABS(ikx(i1,i2)*modesspect(i1,i2) * coeff) &
        *COS(ky_n2(i2)*yvect+angleut(i1,i2))
    dvdt = dvdt - 1.0_rp*ABS(iky(i1,i2)*modesspect(i1,i2) * coeff) &
        *COS(ky_n2(i2)*yvect+anglevt(i1,i2))
    dwdt = dwdt - 1.0_rp*ABS(kth(i1,i2)*modesspect(i1,i2) * coeff2) &
        *COS(ky_n2(i2)*yvect+anglewt(i1,i2))
ENDIF
! i2 and i1 =/ 1


!$OMP PARALLEL DO PRIVATE(n1,n2,k_n2,coeff,coeff2,vitx_l,vity_l,vitz_l,phit_l,dudt_l,dvdt_l,dwdt_l) &
!$OMP & REDUCTION(+:vitx) REDUCTION(+:vity) REDUCTION(+:vitz) &
!$OMP & REDUCTION(+:phit) REDUCTION(+:dudt) REDUCTION(+:dvdt) REDUCTION(+:dwdt) 
DO iTotal = 1, n1o2p1*n2
    n1  = MOD(iTotal-1,n1o2p1)   +1;
    n2  = INT((iTotal-1)/n1o2p1) +1;
!DO i1=2,n1o2p1
!    DO i2=1,n2
        k_n2 = SQRT(kx(i1)**2+ky_n2(i2)**2)
        IF ((k_n2*(zvect+depth_star).LT.50.).AND.(k_n2*depth_star.LT.50.)) THEN
            coeff = COSH(k_n2*(zvect+depth_star))/COSH(k_n2*depth_star) * EXP(i*ky_n2(i2)*yvect)
            coeff2= SINH(k_n2*(zvect+depth_star))/SINH(k_n2*depth_star) * EXP(i*ky_n2(i2)*yvect)
        ELSE
            coeff = EXP(k_n2*zvect) * EXP(i*ky_n2(i2)*yvect)
            coeff2= coeff
        ENDIF
        !
        vitx_l = modesspecx(i1,i2) * coeff
        vity_l = modesspecy(i1,i2) * coeff
        vitz_l = modesspecz(i1,i2) * coeff2
        phit_l = modesspect(i1,i2) * coeff
        dudt_l = ikx(i1,i2)*phit_l
        dvdt_l = iky(i1,i2)*phit_l
        dwdt_l = kth(i1,i2)*modesspect(i1,i2) * coeff2
        !
        vitx = vitx + 1.0_rp*ABS(vitx_l)*COS(kx(i1)*xvect+ATAN2(AIMAG(vitx_l),REAL(vitx_l,RP)))
        vity = vity + 1.0_rp*ABS(vity_l)*COS(kx(i1)*xvect+ATAN2(AIMAG(vity_l),REAL(vity_l,RP)))
        vitz = vitz + 1.0_rp*ABS(vitz_l)*COS(kx(i1)*xvect+ATAN2(AIMAG(vitz_l),REAL(vitz_l,RP)))
        phit = phit + 1.0_rp*ABS(phit_l)*COS(kx(i1)*xvect+ATAN2(AIMAG(phit_l),REAL(phit_l,RP)))
        dudt = dudt  + 1.0_rp*ABS(dudt_l)*COS(kx(i1)*xvect+ATAN2(AIMAG(dudt_l),REAL(dudt_l,RP)))
        dvdt = dvdt + 1.0_rp*ABS(dvdt_l)*COS(kx(i1)*xvect+ATAN2(AIMAG(dvdt_l),REAL(dvdt_l,RP)))
        dwdt = dwdt + 1.0_rp*ABS(dwdt_l)*COS(kx(i1)*xvect+ATAN2(AIMAG(dwdt_l),REAL(dwdt_l,RP)))
!    ENDDO
ENDDO
!$OMP END PARALLEL DO




pit = - g_star*zvect - 0.5_rp*(vitx**2+vity**2+vitz**2)-phit;


!vitx,vity,vitz,pit

!vitx =vitx * L / T;    !velocity U convert to SI unit
!vity =vity * L / T;    !velocity V convert to SI unit
!vitz =vitz * L / T;    !velocity W convert to SI uint 

!pit = density * pit *L**2/T**2

!rFPResultScale
vitx =rPFResultScale * vitx * L / T;    !velocity U convert to SI unit
vity =rPFResultScale * vity * L / T;    !velocity V convert to SI unit
vitz =rPFResultScale * vitz * L / T;    !velocity W convert to SI uint 

pit =rPFResultScale * density * pit *L**2/T**2

!
!compute z_star

!zvect(ii1,ii2)=zmin+(-zmin+eta_tmp(ii1+imin-1,ii2+jmin-1))*SIN(pio2*REAL(ii3-1,RP)/REAL(i_zvect-1,RP))


!-----conver the results to SI UNIT
!ee = eta_star * L  ;    !z coordinate, convert to SI unit

!uu =vitx * L / T;    !velocity U convert to SI unit
!vv =vity * L / T;    !velocity V convert to SI unit
!ww =vitz * L / T;    !velocity W convert to SI uint 
!pit = - g_star*zvect - 0.5_rp*(vitx**2+vity**2+vitz**2)-phit;
!PP = - g_star*zvect - 0.5_rp*(vitx**2+vity**2+vitz**2)-phit;
!pp = density*PP *L**2/T**2
END SUBROUTINE
!---------------------------------------------------------------------------------------------------------------------
!  Subroutine   :                   R8SORT_BUBBLE
!---------------------------------------------------------------------------------------------------------------------
!
!  Purpose      : sort the integer array in increasing order or decreasing order using
!                      selection method. 
!
!  Input        :  N : number of iterms in the array,
!                  X : array hold the number to be sorted.
!                  ORDER : 1 increasing order
!                         -1 decreasing order
!                  P : return index of the array;
!                 
!  Input/output : 
!
!  Output       : 
!
!  Routines     : 
!                 
!
!  Remarks      :
!
!  References   :
!
!  Revisions    :
!---------------------------------------------------------------------------------------------------------------------
SUBROUTINE R8SORT_BUBBLE(N,X,P,ORDER)

IMPLICIT NONE
INTEGER(KIND=4)                 ,INTENT(IN)    ::N;
REAL(KIND=8)    ,DIMENSION(*)   ,INTENT(IN)    ::X;
INTEGER(KIND=4) ,DIMENSION(*)   ,INTENT(  OUT) ::P;
INTEGER(KIND=4)                 ,OPTIONAL      ::ORDER;

!-------------------------------------------------
!Local varialbe.
!-------------------------------------------------
INTEGER(KIND=4)           :: iOrder,i,j,TEMP2;
REAL(KIND=8)              :: TEMP;
IF( PRESENT(ORDER) == .FALSE.) THEN
    iOrder = 1;
ELSE
    iOrder = order;
ENDIF

DO i =1,N
    P(i) = i;
ENDDO

!---------increasing order-------------
IF( iOrder ==1) THEN
    DO I=1,N
    !Now i need to check the numbers
    !and swap their location according the condition
        DO J=1,N-1 
            IF (X(P(J))>X(P(J+1))) THEN
!                TEMP = X(J)
!                X(J)=X(J+1)
!                X(J+1)=TEMP
                
                TEMP2  = P(J)
                P(J)  = P(J+1);
                P(J+1)= TEMP2;
            END IF
        END DO !  DO J=1,N-1 
    END DO !DO I=1,N
    
ELSEIF( iOrder ==-1) THEN

    DO I=1,N
    !Now i need to check the numbers
    !and swap their location according the condition
        DO J=1,N-1 
            IF (X(P(J))<X(P(J+1))) THEN
!                TEMP = X(J)
!                X(J)=X(J+1)
!                X(J+1)=TEMP
                
                TEMP2  = P(J)
                P(J)   = P(J+1);
                P(J+1) = TEMP2;
                
            END IF
        END DO !  DO J=1,N-1 
    END DO !DO I=1,N

ENDIF

END SUBROUTINE


!
SUBROUTINE read_datum_i(unit, input)
!
IMPLICIT NONE
! Input variables
INTEGER, INTENT(IN)  :: unit
! Output variables
INTEGER, INTENT(OUT) :: input
!
line_counter = line_counter + 1
CALL build_read_format('I')
READ(unit,format_read,ERR=101,END=101) description, input
WRITE(*,'(A," : ",I5)') description(1:N_descr), input
RETURN
101 CALL error_message(description)
!
END SUBROUTINE read_datum_i
!
!
!
SUBROUTINE read_datum_r(unit, input)
!
IMPLICIT NONE
! Input variables
INTEGER, INTENT(IN)   :: unit
! Output variables
REAL(RP), INTENT(OUT) :: input
!
line_counter = line_counter + 1
CALL build_read_format('R')
READ(unit,format_read,ERR=101,END=101) description, input
WRITE(*,'(A," : ",ES25.16)') description(1:N_descr), input
IF (ABS(input) > tiny .AND. ABS(input) * 1.0E+16_rp < 1.0E+5_rp) THEN
    WRITE(*,'(A,A)') 'Numeric point is probably missing in current input ',description
    STOP 1
ENDIF
RETURN
101 CALL error_message(description)
!
END SUBROUTINE read_datum_r
!
!
!
SUBROUTINE read_datum_c(unit, input)
!
IMPLICIT NONE
! Input variables
INTEGER, INTENT(IN)   :: unit
! Output variables
CHARACTER(LEN=*), INTENT(OUT) :: input
!
line_counter = line_counter + 1
CALL build_read_format('A')
READ(unit,format_read,ERR=101,END=101) description, input
WRITE(*,'(A," : ",A)') description(1:N_descr), input
RETURN
101 CALL error_message(description)
!
END SUBROUTINE read_datum_c
!
!
!
SUBROUTINE read_blank_line(unit)
!
IMPLICIT NONE
! Input variables
INTEGER, INTENT(IN)   :: unit
!
line_counter = line_counter + 1
READ(unit,*)
!
END SUBROUTINE read_blank_line
!
!
SUBROUTINE build_read_format(code)
!
IMPLICIT NONE
!
CHARACTER(LEN=*) :: code
!
format_read(0) = '(A'
format_read(1) = int2str(N_tot)
format_read(2) = ','
SELECT CASE (code)
    CASE('I')          ! Integer
        format_read(3) = 'I5'
    CASE('F','R')      ! Real number
        format_read(3) = 'ES25.16'
    CASE('S','C','A')  ! Character string
        format_read(3) = 'A'
END SELECT
format_read(4) = ')'
! WRITE(*,*) format_read
!
END SUBROUTINE build_read_format

SUBROUTINE error_message(description)
!
IMPLICIT NONE
! Input variables
CHARACTER(LEN=N_tot) :: description
!
WRITE(*,'(A,I2)') 'Error while reading the input file on line: ', line_counter
WRITE(*,'(A)') description
STOP 1
!
END SUBROUTINE error_message
END MODULE
    