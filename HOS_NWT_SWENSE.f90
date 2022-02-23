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
MODULE  HOS_NWT_SWENSE_MOD
!This module contains the subroutine to read HOS NWT SWENSE file and compute the surface elevation and velocity 
!from the file.
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

!This subroutine is referenced and revised from HOS NWT software  https://github.com/LHEEA/HOS-NWT/wiki
!The Copyright/License of HOS-NWT is present here 
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
USE HOS_Ocean_SWENSE_MOD
IMPLICIT NONE

PUBLIC :: get_HOS_NWT_pt_eta, get_HOS_NWT_pt_deta_dt,get_HOS_NWT_pt_uvwp
PUBLIC :: read_HOS_NWT_mod,HOS_NWT_init
PUBLIC :: n1




!--The 
!=.TRUE. The HOS mode is selected, hence
LOGICAL,PRIVATE      ::lHOSModeSelected =.FALSE.

PRIVATE
INTERFACE get_HOS_NWT_pt_eta
    MODULE PROCEDURE get_pt_eta
END INTERFACE


INTERFACE get_HOS_NWT_pt_deta_dt
    MODULE PROCEDURE get_pt_deta_dt
END INTERFACE


INTERFACE get_HOS_NWT_pt_uvwp
    MODULE PROCEDURE get_pt_uvwp
END INTERFACE 

INTERFACE read_HOS_NWT_mod
    MODULE PROCEDURE read_NWT_mod
END INTERFACE

INTERFACE HOS_NWT_init
    MODULE PROCEDURE recons_HOS_NWT_init
END INTERFACE

!PUBLIC :: recons_HOS_init, get_pt_eta,get_pt_uvwp,read_mod 
!iModeIndexNx: modes index selected for X direction
INTEGER,DIMENSION(:,:),ALLOCATABLE::iModeUsedIndex

CHARACTER(LEN=40)        :: HOSFileName ="modes_HOS_SWENSE.dat"
INTEGER,PARAMETER,PRIVATE:: ntn=2;

REAL(8),PRIVATE          :: rPFResultScale =1.0d0
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
REAL(RP), PARAMETER    :: density=1000.0_RP;
!
! For comparison of real numbers
!REAL(RK),PARAMETER:: EPSILON    =1.0E-6_RK
REAL(RP), PARAMETER    :: tiny = 1.0E-8_RP ! epsilon(1.0_rp)
!

INTEGER :: n1,n2,n3,n3_add
REAL(RP), ALLOCATABLE, DIMENSION(:)      :: x,y
REAL(RP), ALLOCATABLE, DIMENSION(:,:)    :: eta, phis

!---HOS-NWT addtional mode
REAL(RP), ALLOCATABLE, DIMENSION(:)      :: kx,ky
REAL(RP), ALLOCATABLE, DIMENSION(:,:)    :: k,kth
COMPLEX(CP), ALLOCATABLE, DIMENSION(:,:) :: ikx,iky


!---HOS-NWT addtional mode. 
REAL(RP), ALLOCATABLE, DIMENSION(:)     :: kx_add,kx2_add
REAL(RP), ALLOCATABLE, DIMENSION(:,:)   :: k_add,k_add_2,k_add_thk_add
REAL(RP), ALLOCATABLE, DIMENSION(:,:,:) :: csh_add_x,k_add_sh_add_x,kycsh_add_x,kx_add_csh_add_x
!
! Input file
INTEGER               :: i_ana, i_card, tecplot, i_zvect
REAL(RP)              :: T_start, x_min, x_max, y_min, y_max, z_min, z_max
CHARACTER(LEN=100)    :: file_mod

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
INTEGER :: i_unit
INTEGER :: l_add
! test for fourier
!INTEGER :: m1,m2,m3,m3_add,Nd1,Nd2,md1,md2
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!INTEGER :: i_xvect,i_yvect,nzvect

INTEGER :: n1o2p1,n2o2p1
! test for fourier
INTEGER :: m1,m2,m3,Nd1,Nd2,Nd1o2p1,m1o2p1,md1o2p1,md1,md2
!---HOS-NWT added mode
INTEGER :: m3_add

!INTEGER :: nz
REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: modesspecx,modesspecy,modesspecz,modesspect,modesFS,modesFSt,modesadd,modesaddt
REAL(RP), ALLOCATABLE, DIMENSION(:)   :: xvect,yvect,zvect
REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: vitx,vity,vitz,phit,dudt,dvdt,dwdt,zlocal
REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: phixadd,phizadd,phiyadd,phitadd,phixtadd,phiytadd,phiztadd
!
REAL(RP) :: dt_out_star,T_stop_star,xlen_star,ylen_star,depth
!INTEGER  :: imin,imax,jmin,jmax,i_xvect,i_yvect,nzvect
INTEGER  :: i1,i2,i3,i_test
REAL(RP) :: g,g_star
!Normalzie 
REAL(RP) :: T_adim, L_adim,L, T;
REAL(RP) :: tiny_sp = epsilon(1.0);
!REAL(RP), PARAMETER    :: tiny = 1.0E-8_RP ! epsilon(1.0_rp)

REAL(RP) :: FNZ_VALUE =2.0d0;
!
CONTAINS
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

!CALL assign_pointer_arrays(it);

xx = xxl/L
yy = yyl/L


!use all the modes in the results.
IF(iHOSModeMethod ==1) THEN
    get_pt_deta_dt = 0.0_RP
    DO i2=1,n2 
        DO i1=1,n1
!            rtmp =  1.0_rp*REAL(modesFS(i1,i2) * EXP(i*kx(i1)*xx)*EXP(i*ky_n2(i2)*yy));
            get_pt_deta_dt= get_pt_deta_dt  + 1.0_rp*REAL(modesFSt(i1,i2) * EXP(i*kx(i1)*xx)*EXP(i*ky(i2)*yy))
        ENDDO
    ENDDO
   
ELSEIF(iHOSModeMethod ==2) THEN
!IF(iHOSModeMethod ==2) THEN
    get_pt_deta_dt = 0.0_RP
    DO i2=1,n2
        DO ix =1,nHOSModeUsed !loop the selected modes.
            i1 = iModeUsedIndex(ix,i2);  !get the modes index
            get_pt_deta_dt= get_pt_deta_dt  + 1.0_rp*REAL(modesFSt(i1,i2) * EXP(i*kx(i1)*xx)*EXP(i*ky(i2)*yy))
        ENDDO
    ENDDO 

ENDIF

 
!get_pt_eta = get_pt_eta * L
get_pt_deta_dt = rPFResultScale * get_pt_deta_dt * L / T             

END FUNCTION
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

!CALL assign_pointer_arrays(it);

xx = xxl/L
yy = yyl/L
!xx = xxl
!yy = yyl

!use all the modes in the results.
IF(iHOSModeMethod ==1) THEN
    get_pt_eta = 0.0_RP
    DO i2=1,n2 
        DO i1=1,n1
!            rtmp =  1.0_rp*REAL(modesFS(i1,i2) * EXP(i*kx(i1)*xx)*EXP(i*ky_n2(i2)*yy));
            get_pt_eta= get_pt_eta  + 1.0_rp*REAL(modesFS(i1,i2) * EXP(i*kx(i1)*xx)*EXP(i*ky(i2)*yy))
        ENDDO
    ENDDO
   
ELSEIF(iHOSModeMethod ==2) THEN
!IF(iHOSModeMethod ==2) THEN
    get_pt_eta = 0.0_RP
    DO i2=1,n2
        DO ix =1,nHOSModeUsed !loop the selected modes.
            i1 = iModeUsedIndex(ix,i2);  !get the modes index
            get_pt_eta= get_pt_eta  + 1.0_rp*REAL(modesFS(i1,i2) * EXP(i*kx(i1)*xx)*EXP(i*ky(i2)*yy))
        ENDDO
    ENDDO 

ENDIF

 
!get_pt_eta = get_pt_eta * L
get_pt_eta = rPFResultScale * get_pt_eta * L              

END FUNCTION


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
SUBROUTINE get_pt_uvwp(it,xvectl,yvectl,zvectl,vitxl,vityl,vitzl,pitl)
!
IMPLICIT NONE
!% INPUT VARIABLES
!COMPLEX(CP), DIMENSION(m1o2p1,m2), INTENT(IN) :: modesspecx,modesspecy,modesspecz,modesspect
!REAL(RP)                         , INTENT(IN) :: depth_star,g_star
!% INPUT VARIABLES
INTEGER(4),  INTENT(IN) :: it;
REAL(RP)   , INTENT(IN) :: xvectl,yvectl,zvectl;    !non dimensional
REAL(RP)   , INTENT(OUT):: vitxl,vityl,vitzl,pitl
!----------------------------------------------------------
INTEGER                 :: i1,i2,ii1,ii2,ii3,ix
INTEGER                 :: imin,jmin;
!
REAL(RP)                :: xvect(1),yvect(1),zvect(1,1);    !non dimensional
!-------Local variables
REAL(RP)                :: phit(1,1),dudt(1,1),dvdt(1,1),dwdt(1,1);
REAL(RP)                :: coeff, coeff2, coscos, cossin, sincos
!INTEGER                :: nTotal,iTotal
REAL(RP) :: k_n2
!--added component
REAL(RP), DIMENSION(1,1):: phixadd,phiyadd,phizadd,phitadd,phixtadd,phiytadd,phiztadd
REAL(RP)                :: coskx_add,sinkx_add,cosky,coeff1,coeff3,coeff4
REAL(RP)                :: vitx(1,1),vity(1,1),vitz(1,1),pit(1,1)

REAL(RP) :: k_add_X_max, k_add_x, k_add_x1, k_add_x2, expon1, expon2, expon12, expon3, expon3p1, csh_add, sh_add
REAL(RP), DIMENSION(n3_add) :: cs_add, s_add, k_add_s_add ! check if it is useful those values???
REAL(RP), DIMENSION(1,n2,n3_add) :: csh_add_x,k_add_sh_add_x,kycsh_add_x,kx_add_csh_add_x
!CALL assign_pointer_arrays(it);

xvect = xvectl / L 
yvect = yvectl / L
zvect = zvectl / L
!xvect = xvectl  
!yvect = yvectl 
!zvect = zvectl 
!
! constant mode
i1 = 1
i2 = 1

imin = 1;
jmin = 1;
!
! Constant mode
ii1=1
ii2=1
!!
coscos = COS(kx(ii1)*xvect(i1))*COS(ky(ii2)*yvect(i2))
cossin = COS(kx(ii1)*xvect(i1))*SIN(ky(ii2)*yvect(i2))
sincos = SIN(kx(ii1)*xvect(i1))*COS(ky(ii2)*yvect(i2))


!
vitx(i1,i2) = modesspecx(ii1,ii2)*sincos
vity(i1,i2) = modesspecy(ii1,ii2)*cossin
vitz(i1,i2) = modesspecz(ii1,ii2)*coscos
phit(i1,i2) = modesspect(ii1,ii2)*coscos
dudt(i1,i2) = -kx(ii1)*modesspect(ii1,ii2)*sincos
dvdt(i1,i2) = -ky(ii2)*modesspect(ii1,ii2)*cossin
dwdt(i1,i2) = kth(ii1,ii2)*modesspect(ii1,ii2)*coscos
!
ii1=1
DO ii2=2,n2
	IF ((k(ii1,ii2)*(zvect(i1,i2)+1.0_rp).LT.50.).AND.(k(ii1,ii2).LT.50.)) THEN
		coeff = COSH(k(ii1,ii2)*(zvect(i1,i2)+1.0_rp))/COSH(k(ii1,ii2))
		coeff2= SINH(k(ii1,ii2)*(zvect(i1,i2)+1.0_rp))/SINH(k(ii1,ii2))
	ELSE
		coeff = EXP(k(ii1,ii2)*zvect(i1,i2))
		coeff2= coeff
	ENDIF  
    
    if (coeff.ge.FNZ_VALUE) then
        coeff = FNZ_VALUE
    endif

    if (coeff2.ge.FNZ_VALUE) then
        coeff2 = FNZ_VALUE
    endif    
    
	!
	coscos = COS(kx(ii1)*xvect(i1))*COS(ky(ii2)*yvect(i2))
	cossin = COS(kx(ii1)*xvect(i1))*SIN(ky(ii2)*yvect(i2))
	sincos = SIN(kx(ii1)*xvect(i1))*COS(ky(ii2)*yvect(i2))
	!
	vitx(i1,i2) = vitx(i1,i2)+modesspecx(ii1,ii2)*coeff*sincos
	vity(i1,i2) = vity(i1,i2)+modesspecy(ii1,ii2)*coeff*cossin
	vitz(i1,i2) = vitz(i1,i2)+modesspecz(ii1,ii2)*coeff2*coscos
	phit(i1,i2) = phit(i1,i2)+modesspect(ii1,i2)*coeff*coscos
	dudt(i1,i2) = dudt(i1,i2)-kx(ii1)*modesspect(ii1,i2)*coeff*sincos
	dvdt(i1,i2) = dvdt(i1,i2)-ky(ii2)*modesspect(ii1,i2)*coeff*cossin
	dwdt(i1,i2) = dwdt(i1,i2)+kth(ii1,ii2)*modesspect(ii1,ii2)*coeff2*coscos
ENDDO
!
DO ii1=2,n1-n1/16
	DO ii2=1,n2
		IF ((k(ii1,ii2)*(zvect(i1,i2)+1.0_rp).LT.50.).AND.(k(ii1,ii2).LT.50.)) THEN
			coeff = COSH(k(ii1,ii2)*(zvect(i1,i2)+1.0_rp))/COSH(k(ii1,ii2))
			coeff2= SINH(k(ii1,ii2)*(zvect(i1,i2)+1.0_rp))/SINH(k(ii1,ii2))
		ELSE
			coeff = EXP(k(ii1,ii2)*zvect(i1,i2))
			coeff2= coeff
		ENDIF 
        
        if (coeff.ge.FNZ_VALUE) then
            coeff = FNZ_VALUE
        endif

        if (coeff2.ge.FNZ_VALUE) then
            coeff2 = FNZ_VALUE
        endif           
		!
		coscos = COS(kx(ii1)*xvect(i1))*COS(ky(ii2)*yvect(i2))
		cossin = COS(kx(ii1)*xvect(i1))*SIN(ky(ii2)*yvect(i2))
		sincos = SIN(kx(ii1)*xvect(i1))*COS(ky(ii2)*yvect(i2))
		!
		vitx(i1,i2) = vitx(i1,i2)+modesspecx(ii1,ii2)*coeff*sincos
		vity(i1,i2) = vity(i1,i2)+modesspecy(ii1,ii2)*coeff*cossin
		vitz(i1,i2) = vitz(i1,i2)+modesspecz(ii1,ii2)*coeff2*coscos
		phit(i1,i2) = phit(i1,i2)+modesspect(ii1,i2)*coeff*coscos
		dudt(i1,i2) = dudt(i1,i2)-kx(ii1)*modesspect(ii1,i2)*coeff*sincos
		dvdt(i1,i2) = dvdt(i1,i2)-ky(ii2)*modesspect(ii1,i2)*coeff*cossin
		dwdt(i1,i2) = dwdt(i1,i2)+kth(ii1,ii2)*modesspect(ii1,ii2)*coeff2*coscos
	ENDDO
ENDDO	

!------------------compute the wave mode part for current postion xvect,yvect,zvect
!NOTE: this is importat, which is not in HOS-NWT post-processing
!Add Xu Haihua 2020-06-04
!	    
!CALL compute_add_csh_k_kycsh_kx(xlen_star);
! wavemaker modes part
k_add_X_max = 700.0_rp
!
DO i3 = 1, n3_add
   cs_add(i3)   = COS(kx_add(i3))
   s_add(i3)    = SIN(kx_add(i3))
   k_add_s_add(i3) = kx_add(i3) * SIN(kx_add(i3))
   !
   DO i2 = 1, n2
      k_add_x = k_add_2(i3,i2) * xlen_star
      !
      IF (k_add_x <= k_add_X_max) THEN
         expon3 = EXP(-k_add_2(i3,i2) * xlen_star)
      ELSE
         expon3 = 0.0_rp  
      END IF
      !
      expon3p1  = expon3 + 1.d0
!      DO i1 = 1, n1	                             ! on the free surface
!         k_add_x1 = k_add(i3,i2) * x(i1)   ! Change to following Line Xu Haihua 2020-06-03
         k_add_x1 = k_add(i3,i2) * xvect(i1)
         !
         IF (k_add_x1 <= k_add_X_max) THEN
            expon1    = EXP(-k_add_x1) / expon3p1
         ELSE
            expon1 = 0.0_rp  
         END IF
         !
!         k_add_x2 = k_add_2(i3,i2) * (xlen - x(i1))   ! Change to following Line Xu Haihua 2020-06-03
         k_add_x2 = k_add_2(i3,i2) * (xlen_star -xvect(i1))   
         IF (k_add_x2 <= k_add_X_max) THEN
            expon2    = EXP(- k_add_x2)
         ELSE
            expon2 = 0.0_rp  
         END IF
         !
         IF (k_add_x1 + k_add_x2 <= k_add_x_max) THEN
            expon12   = expon1 * expon2
         ELSE
            expon12   = 0.0_rp
         END IF
         !
         csh_add     = expon1 + expon12
         sh_add      = expon1 - expon12
         csh_add_x(i1,i2,i3)     = csh_add
         k_add_sh_add_x(i1,i2,i3)   = sh_add  * k_add(i3,i2)
         kycsh_add_x(i1,i2,i3)      = csh_add * ky(i2)
         kx_add_csh_add_x(i1,i2,i3) = csh_add * kx_add(i3)
!      ENDDO
   ENDDO
ENDDO
i1 =1 
i2 =1

phixadd(i1,i2)  = 0.0_rp
phiyadd(i1,i2)  = 0.0_rp
phizadd(i1,i2)  = 0.0_rp
phitadd(i1,i2)  = 0.0_rp
phixtadd(i1,i2) = 0.0_rp
phiytadd(i1,i2) = 0.0_rp
phiztadd(i1,i2) = 0.0_rp
!
DO ii2=1,n2
	DO ii3=1,n3_add
		coskx_add = cos(kx_add(ii3) * (zvect(i1,i2)+1.0_rp))
		sinkx_add = sin(kx_add(ii3) * (zvect(i1,i2)+1.0_rp))
		cosky     = cos(ky(ii2)*y(i2+jmin-1))
		!
		coeff1 = coskx_add * csh_add_x(i1+imin-1,ii2,ii3)        * cosky
		coeff2 = sinkx_add * kx_add_csh_add_x(i1+imin-1,ii2,ii3) * cosky
		coeff3 = coskx_add * k_add_sh_add_x(i1+imin-1,ii2,ii3)   * cosky
		coeff4 = coskx_add * kycsh_add_x(i1+imin-1,ii2,ii3)      * sin(ky(ii2)*y(i2+jmin-1))
		!
		phitadd(i1,i2)    =  phitadd(i1,i2)  +  modesaddt(ii3,ii2) * coeff1 
		phixadd(i1,i2)    =  phixadd(i1,i2)  -  modesadd(ii3,ii2)  * coeff3 
		phiyadd(i1,i2)    =  phiyadd(i1,i2)  -  modesadd(ii3,ii2)  * coeff4 
		phizadd(i1,i2)    =  phizadd(i1,i2)  -  modesadd(ii3,ii2)  * coeff2 
		!
		phixtadd(i1,i2)    =  phixtadd(i1,i2)  -  modesaddt(ii3,ii2) * coeff3 
		phiytadd(i1,i2)    =  phiytadd(i1,i2)  -  modesaddt(ii3,ii2) * coeff4
		phiztadd(i1,i2)    =  phiztadd(i1,i2)  -  modesaddt(ii3,ii2) * coeff2 
	ENDDO
ENDDO

vitx = vitx + phixadd
vity = vity + phiyadd
vitz = vitz + phizadd
phit = phit + phitadd
! dudt, dvdt, dwdt : total quantities
dudt = dudt + phixtadd
dvdt = dvdt + phiytadd
dwdt = dwdt + phiztadd

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

vitxl = vitx(1,1)
vityl = vity(1,1)
vitzl = vitz(1,1)
pitl  = pit(1,1)

CONTAINS

SUBROUTINE compute_add_csh_k_kycsh_kx(xlen)  

REAL(RP),INTENT(IN)  :: xlen

! wavemaker modes part
k_add_X_max = 700.0_rp
!
DO i3 = 1, n3_add
   cs_add(i3)   = COS(kx_add(i3))
   s_add(i3)    = SIN(kx_add(i3))
   k_add_s_add(i3) = kx_add(i3) * SIN(kx_add(i3))
   !
   DO i2 = 1, n2
      k_add_x = k_add_2(i3,i2) * xlen
      !
      IF (k_add_x <= k_add_X_max) THEN
         expon3 = EXP(-k_add_2(i3,i2) * xlen)
      ELSE
         expon3 = 0.0_rp  
      END IF
      !
      expon3p1  = expon3 + 1.d0
!      DO i1 = 1, n1	                             ! on the free surface
!         k_add_x1 = k_add(i3,i2) * x(i1)   ! Change to following Line Xu Haihua 2020-06-03
         k_add_x1 = k_add(i3,i2) * xvect(i1)
         !
         IF (k_add_x1 <= k_add_X_max) THEN
            expon1    = EXP(-k_add_x1) / expon3p1
         ELSE
            expon1 = 0.0_rp  
         END IF
         !
!         k_add_x2 = k_add_2(i3,i2) * (xlen - x(i1))   ! Change to following Line Xu Haihua 2020-06-03
         k_add_x2 = k_add_2(i3,i2) * (xlen -xvect(i1))   
         IF (k_add_x2 <= k_add_X_max) THEN
            expon2    = EXP(- k_add_x2)
         ELSE
            expon2 = 0.0_rp  
         END IF
         !
         IF (k_add_x1 + k_add_x2 <= k_add_x_max) THEN
            expon12   = expon1 * expon2
         ELSE
            expon12   = 0.0_rp
         END IF
         !
         csh_add     = expon1 + expon12
         sh_add      = expon1 - expon12
         csh_add_x(i1,i2,i3)     = csh_add
         k_add_sh_add_x(i1,i2,i3)   = sh_add * k_add(i3,i2)
         kycsh_add_x(i1,i2,i3)      = csh_add * ky(i2)
         kx_add_csh_add_x(i1,i2,i3) = csh_add * kx_add(i3)
!      ENDDO
   ENDDO
ENDDO


END SUBROUTINE 
END SUBROUTINE



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




!-------------------------------------------------------------
!SUBROUTINE         init_read_NWT_mod
!-------------------------------------------------------------
!> @brief initial the module and read the parameters from the HOS NWT SWESEN file
!> @details  initial the module and read the parameters from the HOS  NWT SWESEN file
!> @param[in] fileName: The HOS SWESEN file name
!> @param[in] i_unit: the file unit
!> @param[OUT] n1: Number of modes in x direction
!> @param[OUT] n2: Number of modes in y direction
!> @param[OUT] n3_add: Additional terms 
!> @param[OUT] dt_out: Time step for output (non-dimensional)
!> @param[OUT] T_stop: Maximum time in the SWESEN file (non-dimensional)
!> @param[OUT] xlen: Length in x direction (non-dimensional)
!> @param[OUT] ylen: Length in y direction (non-dimensional)
!> @param[OUT] depth: Water depth (non-dimensional)
!
!> @date 02 JUN 2020 
!> @remarks None
!-------------------------------------------------------------
SUBROUTINE init_read_NWT_mod(filename,n1,n2,n3_add,dt_out,T_stop,xlen,ylen,depth)
!
! Initialize data from volumic mode description generated by HOS-NWT
! 
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(IN) :: filename
!
INTEGER, INTENT(OUT)         :: n1,n2,n3_add
REAL(RP), INTENT(OUT)        :: dt_out,T_stop,xlen,ylen,depth
!
! Local variables
REAL(RP) :: x1, x2, x3
INTEGER :: i_unit;
!
! We will look at first eight variables written on 18 characters
OPEN(newUnit=i_unit,file=filename,status='OLD', FORM='FORMATTED', ACCESS='DIRECT',RECL=18*8)
READ(i_unit,'(8(ES17.10,1X))',REC=1) x1, x3, x2, dt_out, T_stop, xlen, ylen, depth
!
n1     = NINT(x1)
n2     = NINT(x2)
n3_add = NINT(x3)
!
CLOSE(i_unit)
!
END SUBROUTINE init_read_NWT_mod



!-------------------------------------------------------------
!SUBROUTINE         read_NWT_mod
!-------------------------------------------------------------
!> @brief Read HOS NWT modes from the HOS NWT SWESEN file 
!> @details Read HOS NWT modes from the HOS NWT SWESEN file 
!> @param[in] fileName: The HOS SWESEN file name
!> @param[in] time: the time to read (uit s)
!> @author Xu Haihua, Technology Centre for Offshore and Marine, Singapore (TCOMS)
!> @date 11 Nov 2019 
!> @remarks None
!-------------------------------------------------------------
SUBROUTINE read_NWT_mod(filename,time)
!
! Initialize data from volumic mode description generated by HOS-NWT
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(IN) :: filename

REAL(RP), INTENT(IN)         :: time !, dt_out
!
!REAL(RP), INTENT(OUT), DIMENSION(n1,n2)     :: modesspecx,modesspecy,modesspecz,modesspect,modesFS,modesFSt
!REAL(RP), INTENT(OUT), DIMENSION(n3_add,n2) :: modesadd,modesaddt
!
! Local variables
REAL(RP):: time_star;
INTEGER :: i1, i2, it
INTEGER :: i_unit
INTEGER :: n1, n2, n3_add
!
! We read the specific records corresponding to time
!
n1 = m1;
n2 = m2;
n3_add = m3_add;

time_star = time/T !compute 

it = NINT(time_star/dt_out_star)+1
!
OPEN(newUnit=i_unit,file=filename,status='OLD', FORM='FORMATTED', ACCESS='DIRECT',RECL=18*n1)
!
DO i2=1,n2
	READ(i_unit,'(5000(ES17.10,1X))',REC=((it)*n2*8)+1+8*(i2-1)) (modesspecx(i1,i2), i1=1,n1)
	READ(i_unit,'(5000(ES17.10,1X))',REC=((it)*n2*8)+2+8*(i2-1)) (modesspecy(i1,i2), i1=1,n1)
	READ(i_unit,'(5000(ES17.10,1X))',REC=((it)*n2*8)+3+8*(i2-1)) (modesspecz(i1,i2), i1=1,n1)
	READ(i_unit,'(5000(ES17.10,1X))',REC=((it)*n2*8)+4+8*(i2-1)) (modesspect(i1,i2), i1=1,n1)
	READ(i_unit,'(5000(ES17.10,1X))',REC=((it)*n2*8)+5+8*(i2-1)) (modesFS(i1,i2)   , i1=1,n1)
	READ(i_unit,'(5000(ES17.10,1X))',REC=((it)*n2*8)+6+8*(i2-1)) (modesFSt(i1,i2)  , i1=1,n1)
	READ(i_unit,'(5000(ES17.10,1X))',REC=((it)*n2*8)+7+8*(i2-1)) (modesadd(i1,i2)  , i1=1,n3_add)
	READ(i_unit,'(5000(ES17.10,1X))',REC=((it)*n2*8)+8+8*(i2-1)) (modesaddt(i1,i2) , i1=1,n3_add)
ENDDO
!
CLOSE(i_unit)
!
END SUBROUTINE read_NWT_mod






SUBROUTINE recons_HOS_NWT_init(filename,dt_out,HOS_depth,HOS_g)
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(IN) :: filename
REAL(RP), INTENT(OUT)        :: dt_out,HOS_depth,HOS_g

INTEGER                     :: i_unit
!INTEGER                     :: n1,n2,n3=0,n3_add

dt_out =0.0d0
!-int HOS NWT g is pre-defined.
g     = 9.81_rp;

!
! Initialize variables reading filename
! In the file, everything is non-dimensional with L and T length and time scales
CALL init_read_NWT_mod(filename,n1,n2,n3_add,dt_out_star,T_stop_star,xlen_star,ylen_star,depth)


!
ALLOCATE(modesspecx(n1,n2),modesspecy(n1,n2),modesspecz(n1,n2),modesspect(n1,n2),modesFS(n1,n2),modesFSt(n1,n2))
ALLOCATE(modesadd(n3_add,n2),modesaddt(n3_add,n2))



!
! Initialize Fourier (global variables have to be defined)
m1      = n1
m2      = n2
m3_add  = n3_add
Nd1     = 1
Nd2     = 1
md1     = 1
md2     = 1

n1o2p1 = n1/2+1
n2o2p1 = n2/2+1

!----In HOS the L and T are normalzied as following, which is different from HOS-Ocean.
L_adim = depth
T_adim = SQRT(depth/g)

L = L_adim
T = T_adim
!---
xlen = xlen_star * L
ylen = ylen_star * L;
!--convert the output time to unit value.
dt_out = dt_out_star * T_adim 
T_stop = T_stop_star * T_adim
xlen = xlen_star *L_adim
ylen = ylen_star *L_adim 
!depth_star = 1.0_RP;
!
!depth_star = 1.0d0;
g_star = g /(L_adim/(T_adim*T_adim));
!

!
ALLOCATE(x(n1),y(n2),ky(n2),kx(n1),k(n1,n2),kth(n1,n2),eta(n1,n2))
!
ALLOCATE(kx_add(n3_add),kx2_add(n3_add),k_add(n3_add,n2),k_add_2(n3_add,n2),k_add_thk_add(n3_add,n2))
ALLOCATE(csh_add_x(n1,n2,n3_add),k_add_sh_add_x(n1,n2,n3_add),kycsh_add_x(n1,n2,n3_add),kx_add_csh_add_x(n1,n2,n3_add))


! Read time=0
!CALL read_NWT_mod(filename,0.0_rp,dt_out_star,n1,n2,n3_add)
 CALL read_NWT_mod(filename,0.0_rp);

CALL build_mesh_global_NWT(xlen_star,ylen_star,n1,n2,n3,n3_add,x,y,kx,ky,k,kth,&
	kx_add, kx2_add,k_add, k_add_2, k_add_thk_add,k_add_sh_add_x, kycsh_add_x, kx_add_csh_add_x)


HOS_depth = depth;
HOS_g     = g;
!
!--
END SUBROUTINE recons_HOS_NWT_init


!
SUBROUTINE build_mesh_global_NWT(xlen,ylen,n1,n2,n3,n3_add,x,y,kx,ky,k,kth,&
	kx_add, kx2_add,k_add, k_add_2, k_add_thk_add,k_add_sh_add_x, kycsh_add_x, kx_add_csh_add_x)
!
IMPLICIT NONE
!
REAL(RP), INTENT(IN) :: xlen, ylen
INTEGER, INTENT(IN)  :: n1,n2,n3,n3_add
!!
REAL(RP), DIMENSION(n1,n2), INTENT(OUT)        :: k,kth
REAL(RP), DIMENSION(n1), INTENT(OUT)           :: x,kx
REAL(RP), DIMENSION(n2), INTENT(OUT)           :: y,ky
!
REAL(RP), DIMENSION(n3_add), INTENT(OUT)       :: kx_add, kx2_add
REAL(RP), DIMENSION(n3_add,n2), INTENT(OUT)    :: k_add, k_add_2, k_add_thk_add
!
REAL(RP), DIMENSION(n1,n2,n3_add), INTENT(OUT) :: k_add_sh_add_x, kycsh_add_x, kx_add_csh_add_x
! Local variables
REAL(RP) :: pixlen, piylen, delx, dely, delz
INTEGER  :: i1,i2,i3
!
REAL(RP), DIMENSION(n1) :: kx2
REAL(RP), DIMENSION(n2) :: ky2
REAL(RP), DIMENSION(n3) :: cx_add, cz
!
REAL(RP) :: pixlen_add, xlen_add
REAL(RP) :: k_add_X_max, k_add_x, k_add_x1, k_add_x2, expon1, expon2, expon12, expon3, expon3p1, csh_add, sh_add
REAL(RP), DIMENSION(n3_add) :: cs_add, s_add, k_add_s_add ! check if it is useful those values???
!
INTEGER :: l_add
!
! FIXME: l_add should be put in mod_file probably... so that you do not assume a given size
l_add    = 2
xlen_add = 2.d0 * REAL(l_add)
!
!   mesh generation
delx = xlen / (n1 - 1)
delz = 1.d0 / (n3 - 1)
!
if (n2.ne.1) then
   dely = ylen / (n2 - 1)
else
   dely = 0.d0
endif
!
do i2 = 1, n2
   y(i2) = (i2 - 1) * dely
end do
!
do i1 = 1,n1
   x(i1) = (i1 - 1) * delx
end do
!
do i3=1,n3
   cx_add(i3) = delz*(i3-1)
   cz(i3)   = cx_add(i3)-1.d0
end do
!     __________________________________________________________________
!
!% SPECTRAL RESOLUTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     __________________________________________________________________
!
!	wave numbers
!
pixlen=pi/xlen
!
IF (n2 /= 1) THEN
   piylen=pi/ylen
ELSE
   piylen=0.d0
ENDIF
!
DO i2=1,n2
   ky(i2)  = (i2-1) * piylen
   ky2(i2) = ky(i2) * ky(i2)
ENDDO
!
DO i1 = 1, n1
   kx(i1)  = (i1-1) * pixlen
   kx2(i1) = kx(i1) * kx(i1)
   DO i2 = 1, n2
      k(i1,i2)     = SQRT(kx2(i1) + ky2(i2))
      kth(i1,i2)   = k(i1,i2)*TANH(k(i1,i2))
   ENDDO
ENDDO
!
! on the wavemaker surface
pixlen_add = PI / xlen_add
!
! wavemaker modes part
DO i3 = 1, n3_add
   kx_add(i3)  = (2*i3-1) * pixlen_add
   kx2_add(i3) = kx_add(i3) * kx_add(i3)
   DO i2 = 1, n2
      k_add(i3,i2)  = SQRT(kx2_add(i3) + ky2(i2))
      k_add_2(i3,i2) = 2.0_rp * k_add(i3,i2) ! GD : to check
      k_add_thk_add(i3,i2) = k_add(i3,i2) * dtanh(k_add(i3,i2) * xlen) ! on the wavemaker
   ENDDO
ENDDO
!	    
! wavemaker modes part
k_add_X_max = 700.0_rp
!
DO i3 = 1, n3_add
   cs_add(i3)   = COS(kx_add(i3))
   s_add(i3)    = SIN(kx_add(i3))
   k_add_s_add(i3) = kx_add(i3) * SIN(kx_add(i3))
   !
   DO i2 = 1, n2
      k_add_x = k_add_2(i3,i2) * xlen
      !
      IF (k_add_x <= k_add_X_max) THEN
         expon3 = EXP(-k_add_2(i3,i2) * xlen)
      ELSE
         expon3 = 0.0_rp  
      END IF
      !
      expon3p1  = expon3 + 1.d0
      DO i1 = 1, n1	                             ! on the free surface
         k_add_x1 = k_add(i3,i2) * x(i1)
         !
         IF (k_add_x1 <= k_add_X_max) THEN
            expon1    = EXP(-k_add_x1) / expon3p1
         ELSE
            expon1 = 0.0_rp  
         END IF
         !
         k_add_x2 = k_add_2(i3,i2) * (xlen - x(i1))
         IF (k_add_x2 <= k_add_X_max) THEN
            expon2    = EXP(- k_add_x2)
         ELSE
            expon2 = 0.0_rp  
         END IF
         !
         IF (k_add_x1 + k_add_x2 <= k_add_x_max) THEN
            expon12   = expon1 * expon2
         ELSE
            expon12   = 0.0_rp
         END IF
         !
         csh_add     = expon1 + expon12
         sh_add      = expon1 - expon12
         csh_add_x(i1,i2,i3)     = csh_add
         k_add_sh_add_x(i1,i2,i3)   = sh_add * k_add(i3,i2)
         kycsh_add_x(i1,i2,i3)      = csh_add * ky(i2)
         kx_add_csh_add_x(i1,i2,i3) = csh_add * kx_add(i3)
      ENDDO
   ENDDO
ENDDO
!
END SUBROUTINE build_mesh_global_NWT



END MODULE
    