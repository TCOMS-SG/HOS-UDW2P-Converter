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
MODULE Sigma_Grid_MOD
!-------------------------------------------------------------------------
! MODULE :     Sigma_Grid_MOD
!-------------------------------------------------------------------------
!> @brief This module contains the 2D sigma grid type.
!> @details This module contains the 2D sigma grid type. The sigma sigmae coordinate system is defined as
!> @details sigma = z + h /(eta + h), where
!           z: is vertical coordiante with origin at free surface, upward is +, downward is -
!         eta: free surface elevation in z coordinate system.
!           h: water depth (positive all the time).
!  Public :: TSigmaGrid2D , create_Chebyshev_nodes
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

!---------------------------------------------------------------------------------------------------------------------    
!USE HOS_Ocean_SWENSE_MOD
!USE UDW2p_TPNWT
IMPLICIT NONE
!-----------------------------------------------------
!grid variablese dimensional
!-----------------------------------------------------
!xs: start point 
!xe: end point;
!
INTEGER,PARAMETER ::IK =4;
INTEGER,PARAMETER ::RK =8;
REAL(RK),PARAMETER::PI =3.141592653589793238462643383279502888419_RK

PUBLIC ::TSigmaGrid2D,create_Chebyshev_nodes 

!---2D sigma grid type
!In x direction the grid is evenly distributed.
!In z direction the grid can be unevely distributed and zp always from 0-1
!The sigma coordinate value between [0,1] 0 is at still water level, upward is positive, downward is negative.
!
TYPE TSigmaGrid2D
    REAL(RK)           :: xmin;
    REAL(RK)           :: xmax;
    REAL(RK)           :: dx
    REAL(RK)           :: h   ;    !h is water depth, positive.
    REAL(RK),DIMENSION(:),ALLOCATABLE :: xp,zp  !zp alway from 0 to 1, bottom is 0, free surface is 1.
    INTEGER(IK)::nx,nz
CONTAINS
    PROCEDURE,PASS(this) :: init  =>TSigmaGrid2D_INIT;  
    PROCEDURE,PASS(this) :: S2Z   =>TSigmaGrid2D_transS2Z
    PROCEDURE,PASS(this) :: Z2S   =>TSigmaGrid2D_transZ2S 
END TYPE

PRIVATE:: TSigmaGrid2D_INIT,TSigmaGrid2D_transS2Z,TSigmaGrid2D_transZ2S;

CONTAINS









!-------------------------------------------------------------------------
! SUBROUTINE :    create_Chebyshev_nodes
!> @brief create the chebyshev nodes for later boundary
!> @detail create the chebyshev nodes for later boundary
!> @param[in]  np  number of points.
!> @param[out] xpts coordiante of the Chebyshev_nodes
!> @author Xu Haihua, Technology Centre for Offshore and Marine, Singapore (TCOMS)
!> @date 20 May 2020
!> @reference:  https://en.wikipedia.org/wiki/Chebyshev_nodes
!>
!-------------------------------------------------------------------------
SUBROUTINE create_Chebyshev_nodes(np, xpts)
!-------------------------------------------------------------------------------------------------------------
!  Modules 
!---------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE;
!---------------------------------------------------------------------------------------------------------------------
!  Arguments     
!---------------------------------------------------------------------------------------------------------------------
INTEGER(IK),          INTENT(IN):: np;
REAL(RK),DIMENSION(:),INTENT(OUT)::xpts;
!---------------------------------------------------------------------------------------------------------------------
!  Local Variables.     
!---------------------------------------------------------------------------------------------------------------------
INTEGER(IK) ::k;
REAL(RK)    ::a, b;

a =0.0d0;
b =1.0d0
DO k=np+1,2*np
    xpts(k-np) = 0.5d0*(a+b) + 0.5d0*(b-a)*cos((2.0d0*k-1.0d0)/(2.0d0*np)*PI)
ENDDO

END SUBROUTINE

!-------------------------------------------------------------------------
! SUBROUTINE :     TSigmaGrid2D_INIT
!-------------------------------------------------------------------------
!> @brief Init the type TSigmaGrid2D_INIT
!> @detail Init the type TSigmaGrid2D_INIT
!> @param[in] this    The TSigmaGrid2D_INIT
!> @param[in] nx      Number of grid in x direction.
!> @param[in] nz      Number of grid in z direciton.
!> @param[in] xmin    Minmum coordinate of the system.
!> @param[in] xmax    Maximum coordinate of the system.
!> @param[in] h       Water depth (positive value)
!> @param[in] zp      sigma coordiante distribution array between [0, 1],0 is bottom, 1 is free surface.
!
!>                  
!> @author Xu Haihua, Technology Centre for Offshore and Marine, Singapore (TCOMS)
!> @date 20 May 2020
!>
!  Remarks      : 
!
!  References   :
!
!  Revisions    :
!---------------------------------------------------------------------------------------------------------------------    
SUBROUTINE  TSigmaGrid2D_INIT(this,nx,nz,xmin,xmax,h,zp)
!-------------------------------------------------------------------------------------------------------------
!  Modules 
!---------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE;
!---------------------------------------------------------------------------------------------------------------------
!  Arguments     
!---------------------------------------------------------------------------------------------------------------------
CLASS(TSigmaGrid2D),          INTENT(INOUT) :: this
INTEGER(IK),                   INTENT(IN)   :: nx,nz;
REAL(RK),                     INTENT(IN)    :: xmin,xmax,h
REAL(RK),DIMENSION(:),        INTENT(IN)    :: zp

!---------------------------------------------------------------------------------------------------------------------
!  Local variables 
!---------------------------------------------------------------------------------------------------------------------
INTEGER(IK) :: i,info;


this%nx = nx;
this%nz = nz;
this%xmin = xmin;
this%xmax = xmax
this%h = h;

IF(ALLOCATED(this%xp)) DEALLOCATE(this%xp)
IF(ALLOCATED(this%zp)) DEALLOCATE(this%zp)

ALLOCATE(this%xp(nx), this%zp(nz),STAT=info)
IF(info /=0) THEN
    WRITE(*,*) "Failed to allocate array xp, zp "
ENDIF

this%dx = (xmax - xmin)/(nx-1);

DO i=1,nx
    this%xp(i) = xmin+(i-1)*this%dx;
ENDDO

DO i=1,nz
    this%zp(i) = zp(i);
ENDDO

END SUBROUTINE


!-------------------------------------------------------------
!SUBROUTINE transZ2S
!-------------------------------------------------------------
!> @brief transform cartesian coordinate Z to sigma coordinate sigma
!> @details transform cartesian coordinate Z to sigma coordinate sigma
!> @param[in] eta: surface elevation
!> @param[in] z  : z coordinate (in Cartesian coordinate system)
!> @return       : sigma coordiant $\sigma$ in sigma coordinate system
!
!> @author Xu Haihua, Technology Centre for Offshore and Marine, Singapore (TCOMS)
!> @date 20 May 2020
!> NOTE:after transform sigma = 0: bottom;
!                       sigma = 1: free sufrace;
!                       <0 or > 1: not in the Potential flow domain
REAL(RK) FUNCTION TSigmaGrid2D_transZ2S(this,eta,z) RESULT(sigma)
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------
!  Arguments     
!---------------------------------------------------------------------------------------------------------------------
CLASS(TSigmaGrid2D)     :: this 
REAL(RK), INTENT(IN)    :: eta,z
!---------------------------------------------------------------------------------------------------------------------
!  Local Variables     
!---------------------------------------------------------------------------------------------------------------------
sigma = (this%h+z)/(eta+this%h)

END FUNCTION

!-------------------------------------------------------------
!SUBROUTINE transS2Z
!-------------------------------------------------------------
!> @brief transform sigma coordinate sigma to Z cartesian coordinate Z 
!> @details transform sigma coordinate sigma to Z cartesian coordinate Z
!> @param[in]  h   : water depth 
!> @param[in] eta  : surface elevation
!> @param[in] sigam: sigam (sigam coordinate)
!> @return         : z coordinate
!
!> @author Xu Haihua, Technology Centre for Offshore and Marine, Singapore (TCOMS)
!> @date 20 May 2020
REAL(RK) FUNCTION TSigmaGrid2D_transS2Z(this,eta,sigma) RESULT(z)
IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------
!  Arguments     
!---------------------------------------------------------------------------------------------------------------------
CLASS(TSigmaGrid2D)     :: this 
REAL(RK), INTENT(IN)    :: eta,sigma
!---------------------------------------------------------------------------------------------------------------------
!  Local Variables     
!---------------------------------------------------------------------------------------------------------------------
z = sigma * (eta+this%h) - this%h

END FUNCTION


END MODULE