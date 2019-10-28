!!****if* source/Simulation/SimulationMain/magnetoHD/OrszagTang/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer(in) :: blockID) 
!!                       
!!
!!
!!
!! DESCRIPTION
!!
!!  Initializes fluid data (density, pressure, velocity, etc.) for
!!  a specified block.
!!
!!  Reference: Orszag and Tang, J. Fluid Mech., 90:129--143, 1979
!!
!! 
!! ARGUMENTS
!!
!!  blockID -          the number of the block to update
!!
!! 
!!
!!***

!!REORDER(4): solnData, face[xy]Data

subroutine Simulation_initBlock(blockID)

  use Simulation_data

  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
                             Grid_getCellCoords,     &
                             Grid_getBlkPtr,         &
                             Grid_releaseBlkPtr

  implicit none

#include "constants.h"
#include "Flash.h"

  !!$ Arguments -----------------------
  integer, intent(in) :: blockID
  !!$ ---------------------------------

  integer :: i, j, k, n, istat, sizeX, sizeY, sizeZ
  integer :: axis(MDIM)
  real,allocatable,dimension(:) :: xCoord,yCoord,zCoord
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  real :: enerZone, ekinZone, eintZone
  real, pointer, dimension(:,:,:,:) :: solnData, facexData,faceyData,facezData
  real :: tradActual
  real :: rho, tele, trad, tion, zbar, abar
  integer :: species
#ifndef CHAM_SPEC
  integer :: CHAM_SPEC = 1, TARG_SPEC = 2
#endif

  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

  sizeX = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
  sizeY = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
  sizeZ = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1

  allocate(xCoord(sizeX),stat=istat)
  allocate(yCoord(sizeY),stat=istat)
  allocate(zCoord(sizeZ),stat=istat)

  xCoord = 0.0
  yCoord = 0.0
  zCoord = 0.0

  if (NDIM == 3) call Grid_getCellCoords(KAXIS,blockID,CENTER,.true.,zCoord,sizeZ)
  if (NDIM >= 2) call Grid_getCellCoords(JAXIS,blockID,CENTER,.true.,yCoord,sizeY)
  call Grid_getCellCoords(IAXIS,blockID,CENTER,.true.,xCoord,sizeX)
  !------------------------------------------------------------------------------

  call Grid_getBlkPtr(blockID,solnData,CENTER)

#if NFACE_VARS > 0
  if (sim_killdivb) then
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
     facexData(:,:,:,:)=0.0
     faceyData(:,:,:,:)=0.0

        call Grid_getBlkPtr(blockID,facezData,FACEZ)
        facezData(:,:,:,:) = 0.

  endif
#endif

  ! Loop over cells in the block.
  do k = blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
     do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
        do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
           axis(IAXIS) = i
           axis(JAXIS) = j
           axis(KAXIS) = k

           species = CHAM_SPEC
           if ( yCoord(j) <= sim_targetHeight .and. &
                      yCoord(j) >= sim_vacuumHeight ) then
                    species = TARG_SPEC
           end if
           if(species == TARG_SPEC) then
              rho = sim_rhoTarg
              tele = sim_teleTarg
              tion = sim_tionTarg
              trad = sim_tradTarg
           else
              rho = sim_rhoCham
              tele = sim_teleCham
              tion = sim_tionCham
              trad = sim_tradCham
           end if

           call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, axis, rho)
           call Grid_putPointData(blockId, CENTER, TEMP_VAR, EXTERIOR, axis, tele)

#ifdef FLASH_3T
           call Grid_putPointData(blockId, CENTER, TION_VAR, EXTERIOR, axis, tion)
           call Grid_putPointData(blockId, CENTER, TELE_VAR, EXTERIOR, axis, tele)

           ! Set up radiation energy density:
           call RadTrans_mgdEFromT(blockId, axis, trad, tradActual)
           call Grid_putPointData(blockId, CENTER, TRAD_VAR, EXTERIOR, axis, tradActual)
#endif

if (NSPECIES > 0) then
              ! Fill mass fractions in solution array if we have any SPECIES defined.
              ! We put nearly all the mass into either the Xe material if XE_SPEC is defined,
              ! or else into the first species.
              do n = SPECIES_BEGIN,SPECIES_END
                 if (n==species) then
                    call Grid_putPointData(blockID, CENTER, n, EXTERIOR, axis, 1.0e0-(NSPECIES-1)*sim_smallX)
                 else
                    call Grid_putPointData(blockID, CENTER, n, EXTERIOR, axis, sim_smallX)
                 end if
              enddo
           end if

#ifdef BDRY_VAR
           call Grid_putPointData(blockId, CENTER, BDRY_VAR, EXTERIOR, axis, -1.0)
#endif

           solnData(MAGX_VAR,i,j,k)=  200000.
           solnData(MAGY_VAR,i,j,k)=  0.
           solnData(MAGZ_VAR,i,j,k)=  0.

           solnData(MAGP_VAR,i,j,k)= .5*dot_product(solnData(MAGX_VAR:MAGZ_VAR,i,j,k),&
                                                    solnData(MAGX_VAR:MAGZ_VAR,i,j,k))

           ! Cell face-centered variables for StaggeredMesh scheme
#if NFACE_VARS > 0
           if (sim_killdivb) then
              facexData(MAG_FACE_VAR,i,j,k)= 200000.
              faceyData(MAG_FACE_VAR,i,j,k)= 0.
              facezData(MAG_FACE_VAR,i,j,k)= 0.
           endif
#endif
        enddo
     enddo
  enddo



#if NFACE_VARS > 0
  if (sim_killdivb) then

     do k=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
        i=blkLimitsGC(HIGH,IAXIS)+1
        do j=blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
           facexData(MAG_FACE_VAR,i,j,k)= 200000.
        enddo

        j=blkLimitsGC(HIGH,JAXIS)+1
        do i=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
           faceyData(MAG_FACE_VAR,i,j,k)= 0.
        enddo
     enddo

  endif

     facezData(MAG_FACE_VAR,:,:,:) = 0.
#endif


  ! Release pointer

#if NFACE_VARS > 0
  if (sim_killdivb) then
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)  
endif
#endif

  call Grid_releaseBlkPtr(blockID,solnData,CENTER)

  deallocate(xCoord)
  deallocate(yCoord)
  deallocate(zCoord)

end subroutine Simulation_initBlock






