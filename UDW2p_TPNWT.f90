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
MODULE UDW2p_TPNWT
IMPLICIT NONE
PRIVATE
!---------------------------------------------------------------------------------------------------------------------    
!> @author Xu Haihua, Technology Centre for Offshore and Marine, Singapore (TCOMS)
!> @date 17 Feb 2022
!>
!>
!> The HOS UDW2P Converter program was contributed by 
!> Dr Haihua Xu, Mr Yun Zhi Law and Dr Harrif Santo under the guidance of Prof Allan Magee. 
!> This research is supported by A∗STAR Science and Engineering Research Council with grant number 172 19 00089 under the Marine & Offshore Strategic Research Programme (M&O SRP). 
!> The source code is available at: https://github.com/ittcoms/HOS-UDW2P-Converter
!---------------------------------------------------------------------------------------------------------------------    

!Reference:
! 1. JIP on Reproducible CFD Modeling Practice for Offshore Applications ¨C TAC Monthly Meeting, March 2020
!
!-----------------------------------------------
!IW: The unit number for UNW2P file.
!fileName: UNW2P file name
!Nx: Number of points in x direction.
!Np: Number fo point of left side and right side boundary points.
!Nt: Don't know the meaning, according to the Ref 1, page 11, Nt=0
!recSize: record length  recSize = ( Nx * 2 + Np * 2 )
!         Nx*2 = elev, fphi        (Real Value)
!         Ny*2 = leftPhi,rightPhi  (Real Value)
!headerSize: =1 or
!         headerSize = int( headerSize / recSize) + 1
!dx : grid size in x direction (UDW2P file grid)
!dt : time step size 
!xCalibration: X axis calibration location 
!xScale: non dimensional scale in length
!tScale: non dimensional scale in time 
!potScale: potential scale  = xScale**2 / tScale
!velScale: velocity scale  = xScale / tScale
!-----------------------------------------------
INTEGER, PARAMETER :: IW = 1002
!fileType: the type of UDW2P file type
!         =1: file contain elevation and velocity potentional
!         =2: file contain elevation and velocity

INTEGER            :: fileType = 2;
CHARACTER(*), PARAMETER :: fileName = 'PNWT.rec'            !*needs to be changed as variable
INTEGER,PRIVATE :: Nx, Np, Nt, recSize
INTEGER,PRIVATE:: irec, headerSize 

REAL(8) :: dx, dt, xCalibration
REAL(8) :: xScale, tScale, potScale, velScale

Logical, PUBLIC :: toSave = .false.

!--UDW2p_setUp: crate the UDW2P file and write the header
!--             the file 'PNWT.rec' is open after call this subroutine
!--UDW2p_write: write the data to the file 'PNWT.rec'. The data is :
!        1. free surface elvation (elev) and velocity potential (fphi), data size is Nx
!        2. left and right side boundary velocity potential (leftPhi,rightPhi), data size is Np
PUBLIC :: UDW2p_setUp, UDW2p_write ,UDW2p_write_V2,UDW2p_Close,UDW2p_getGQPoints;

INTERFACE UDW2p_getGQPoints
    MODULE PROCEDURE getGQPoints
END INTERFACE 

CONTAINS



SUBROUTINE UDW2p_setUp(Nx_, Np_, Nt_, dx_, dt_, xCalib_, xScale_, tScale_)      !*add fileName argument
INTEGER, INTENT(IN) :: Nx_, Np_, Nt_
REAL(8), INTENT(IN) :: dx_, dt_, xCalib_, xScale_, tScale_
REAL, ALLOCATABLE :: datBlock(:)
!INTEGER :: headerSize

INTEGER :: isize

Nx = Nx_
Np = Np_
Nt = Nt_
dx = dx_
dt = dt_
xCalibration = xCalib_

xScale = xScale_
tScale = tScale_
potScale = xScale**2 / tScale
velScale = xScale / tScale

CALL UDW2p_writeHeader( 1)
CLOSE(IW)
! check the head size to continue.
open(IW, file = trim(fileName), form='unformatted', access='stream', status='OLD')
INQUIRE(IW, size = headerSize)
CLOSE(IW)

!   Figure out record length (recSize)
IF(fileType ==1) THEN
    recSize = ( Nx * 2 + Np * 2 )
ELSEIF(fileType ==2) THEN
!elev(nx) + U(nx) + W(nx)  =Nx*3
!leftU(np) + leftW(np) + rightU(np) + rightW(np) = Np*4
!    recSize = ( Nx * 3 + Np * 4 )
    recSize = ( Nx * 2 + Np * 2 )
ENDIF

ALLOCATE( datBlock(recSize))
INQUIRE(iolength=recSize) datBlock
DEALLOCATE( datBlock)

headerSize = int( headerSize / recSize) + 1

write(*,*) 'headerSize, recSize:', headerSize, recSize

CALL UDW2p_writeHeader( headerSize)
CLOSE(IW)

open(IW, file = trim(fileName), form='unformatted', access='stream', status='OLD')
INQUIRE(IW, size = headerSize)
CLOSE(IW)

!   Figure out record length (recSize)
!recSize = ( Nx * 2 + Np * 2 )
!   Figure out record length (recSize)
IF(fileType ==1) THEN
    recSize = ( Nx * 2 + Np * 2 )
ELSEIF(fileType ==2) THEN
!elev(nx)  + DelevDt =Nx*2
!leftU(np) + rightU(np) = Np*2
    recSize = ( Nx * 2 + Np * 2 )
ENDIF

ALLOCATE( datBlock(recSize))
INQUIRE(iolength=recSize) datBlock
DEALLOCATE( datBlock)

headerSize = int( headerSize / recSize) + 1

write(*,*) 'headerSize, recSize:', headerSize, recSize

CALL UDW2p_writeHeader( headerSize)
CLOSE(IW)



open(IW, file = trim(fileName), form='unformatted', &
        access='direct', status='OLD', recl=recSize)

irec = headerSize + 1

END SUBROUTINE




SUBROUTINE UDW2p_writeHeader( headerSize)
INTEGER, INTENT(IN) :: headerSize
CHARACTER(72) :: txt
INTEGER :: i

!open(IW, file = trim(fileName) )
!
!rewind(15)
!txt = " "
!do while(txt(1:1).ne."#")
!    read(15,*) txt
!end do
!backspace(15)
!
!do i = 1, 100
!     	read(15,'(a72)',end=999) txt
!        CALL UDW2p_writeLine(txt)
!end do
!
!       
!999 continue
!---Write UDW2P comments
open(IW, file = trim(fileName) )

CALL Write_HOS_Ocean_Comment()

close(IW)

open(IW, file = trim(fileName), access='append' )
    
CALL UDW2p_writeParameter('PRECISION', 'single')
CALL UDW2p_writeParameter('Software', 'HOS')
CALL UDW2p_writeIntParameter('Nheader', headerSize)
CALL UDW2p_writeIntParameter('Nx', Nx)
CALL UDW2p_writeIntParameter('Np', Np)
CALL UDW2p_writeIntParameter('Nt', Nt)
CALL UDW2p_writeDoubleParameter('dx', dx)
CALL UDW2p_writeDoubleParameter('dt', dt)
CALL UDW2p_writeDoubleParameter('xCalib', xCalibration)
CALL UDW2p_writeDoubleParameter('xScale', xScale)
CALL UDW2p_writeDoubleParameter('tScale', tScale)
CALL UDW2p_writeDoubleParameter('WD', xScale)

CONTAINS

SUBROUTINE Write_HOS_Ocean_Comment()

    txt="######################################"
    CALL UDW2p_writeLine(txt)
    
    txt="#                                     "
    CALL UDW2p_writeLine(txt)     
    
    txt="#  Comments on Wave data"
    CALL UDW2p_writeLine(txt)
    
    txt="#  Wind Seed XX, Hs=XX, Tp= XX  (HOS) "
    CALL UDW2p_writeLine(txt)
    
    txt="#                                     "
    CALL UDW2p_writeLine(txt)
    
    txt="######################################"
    CALL UDW2p_writeLine(txt)

END SUBROUTINE

END SUBROUTINE

SUBROUTINE UDW2p_writeLine( key)
CHARACTER(*), INTENT(IN) :: key
 

write(IW,'(a)') key

END SUBROUTINE


SUBROUTINE UDW2p_writeParameter( key, value)
CHARACTER(*), INTENT(IN) :: key, value
 

write(IW,*) trim(key),' = ', adjustl(trim(value))

END SUBROUTINE


SUBROUTINE UDW2p_writeIntParameter( key, value)
CHARACTER(*), INTENT(IN) :: key
INTEGER, INTENT(IN) :: value
CHARACTER(20) :: txt

write(txt,'(i20)') value
write(IW,*) trim(key),' = ', adjustl(trim(txt))

END SUBROUTINE

SUBROUTINE UDW2p_writeDoubleParameter( key, value)
CHARACTER(*), INTENT(IN) :: key
REAL(8), INTENT(IN) :: value
CHARACTER(40) :: txt

write(txt,*) value
write(IW,*) trim(key),' = ', adjustl( trim(txt))

END SUBROUTINE



SUBROUTINE UDW2p_write( elev, fphi, leftPhi, rightPhi )
REAL(8), INTENT(IN) :: elev(Nx), fphi(Nx), leftPhi(Np), rightPhi(Np)

write(IW, rec=irec) real(elev), real(fphi), real(leftPhi), real(rightPhi)
irec = irec + 1
if( mod(irec,100).eq.0) write(*,*) ' irec:', irec
END SUBROUTINE

!-------------------------------------------------------
!write the data to the UDW2P file
!For free surface write eta and d(eta)/dt
!For side boundary write left BC U and right BC U
!-------------------------------------------------------
SUBROUTINE UDW2p_write_V2(elev, elevDt,leftU,rightU)
REAL(8), INTENT(IN) :: elev(Nx),elevDt(Nx), leftU(Np),rightU(Np)

write(IW, rec=irec) real(elev), real(elevDt),real(leftU), real(rightU)
irec = irec + 1
if( mod(irec,100).eq.0) write(*,*) ' irec:', irec
END SUBROUTINE


SUBROUTINE UDW2p_write_vel( elev, U,W, leftU,leftW,rightU,RightW)
!REAL(8), INTENT(IN) :: elev(Nx), fphi(Nx), leftPhi(Np), rightPhi(Np)
REAL(8), INTENT(IN) :: elev(Nx), U(Nx),W(Nx), leftU(Np),leftW(Np),rightU(Np),RightW(Np)

write(IW, rec=irec) real(elev), real(U),real(W), real(leftU), real(leftW), real(rightU), real(rightW)
irec = irec + 1
if( mod(irec,100).eq.0) write(*,*) ' irec:', irec
END SUBROUTINE


SUBROUTINE UDW2p_Close( )
!REAL(8), INTENT(IN) :: elev(Nx), fphi(Nx), leftPhi(Np), rightPhi(Np)
CLOSE(IW);
END SUBROUTINE

! 
! Return 20 sampling points (GQ20) for the vertical profile for UDW2p data 
!   in sigma coordinate normalized as 
!   sigma = 0 is sea bottom 
!   sigma = 1 is free surface 
! 
! 

SUBROUTINE getGQPoints( GQ20 ) 
IMPLICIT NONE 
! SPG: 20-point Gauss quadrature: sampling points and weight in the normalized interval (-1, 1) 
! Because of symmetry, only sampling points in the half internal (positive) are tabulated. 
REAL(8), PARAMETER, dimension(10) :: & 
spg = (/ 0.076526521133497333755d0,& 
         0.227785851141645078080d0,& 
         0.373706088715419560673d0,& 
         0.510867001950827098004d0,& 
         0.636053680726515025453d0,& 
         0.746331906460150792614d0,& 
         0.839116971822218823395d0,& 
         0.912234428251325905868d0,& 
         0.963971927277913791268d0,& 
         0.993128599185094924786d0/) 
  
REAL(8), INTENT(OUT) :: GQ20(20) 
    
INTEGER :: isp 
INTEGER, PARAMETER :: Ng20 = 10 

do isp = 1,Ng20 
    GQ20(isp) = .5d0 - .5d0 *spg(11 - isp) 
    GQ20(10 + isp) = .5d0 + .5d0 *spg(isp) 
end do 
  
END SUBROUTINE 
  





END MODULE
