!===============================================================
! mods_common.f90
!   Replacement for common1.inc using Fortran modules
!   - Keep original COMMON grouping
!   - Make NNP/NND dynamic and allocate arrays at runtime
!===============================================================

!-------------------------
! MPI state (formerly common/mpi/)
!-------------------------
module mod_mpi_state
  implicit none
  public
  integer :: nprocs=1, myrank=0, ista=1, iend=0, ista2=1, iend2=0, ised=0
end module mod_mpi_state


!-------------------------
! Domain (formerly common/domain/)
!   tx, ty, tstart, tz, depmax
!-------------------------
module mod_domain
  use mod_kinds, only: dp
  implicit none
  private
  public :: tx, ty, tstart, tz, depmax, alloc_domain, free_domain, &
       mx, my, mz, ibndtyp, npoly
  
  real(dp), allocatable :: tx(:), ty(:)
  real(dp) :: tstart = 0.0_dp, tz = 0.0_dp, depmax = 0.0_dp
  integer:: mx=0, my=0, mz=0, ibndtyp=-1, npoly=-1
contains
  subroutine alloc_domain()
    if (npoly.le.0) then
      write(*,*) "ERROR: alloc_domain: sizes not set"
      stop
    end if
    if (.not. allocated(tx)) allocate(tx(npoly+1))
    if (.not. allocated(ty)) allocate(ty(npoly+1))
    tx = 0.0_dp
    ty = 0.0_dp
  end subroutine alloc_domain

  subroutine free_domain()
    if (allocated(tx)) deallocate(tx)
    if (allocated(ty)) deallocate(ty)
  end subroutine free_domain
end module mod_domain


!-------------------------
! Data (formerly common/data/)
!   xx, yy, zz, zmg, zdp
!-------------------------
module mod_eq_data
  use mod_kinds, only: dp
  implicit none
  private
  public :: xx, yy, zz, zmg, zdp, alloc_eq_data, free_eq_data,ind, &
    xa,ya, indat, nn, nnc

  real(dp), allocatable :: xx(:), yy(:), zz(:), zmg(:), zdp(:)
  real(dp), allocatable :: xa(:), ya(:)
  integer,  allocatable :: ind(:)
  integer,  allocatable :: indat(:)
  integer:: nn=-1, nnc=-1
  
contains
  subroutine alloc_eq_data()
    if (nn .le.0) then
      write(*,*) "ERROR: alloc_eq_data: sizes not set"
      stop
    end if
    if (.not. allocated(xx)) allocate(xx(nn))
    if (.not. allocated(yy)) allocate(yy(nn))
    if (.not. allocated(zz)) allocate(zz(nn))
    if (.not. allocated(zmg)) allocate(zmg(nn))
    if (.not. allocated(zdp)) allocate(zdp(nn))
    if (.not. allocated(xa))    allocate(xa(nn))
    if (.not. allocated(ya))    allocate(ya(nn))
    if (.not. allocated(indat)) allocate(indat(nn))
    if (.not. allocated(ind)) allocate(ind(nn))
    xx = 0.0_dp; yy = 0.0_dp; zz = 0.0_dp; zmg = 0.0_dp; zdp = 0.0_dp
    xa = 0.0_dp; ya = 0.0_dp; indat = 0
  end subroutine alloc_eq_data

  subroutine free_eq_data()
    if (allocated(xx))  deallocate(xx)
    if (allocated(yy))  deallocate(yy)
    if (allocated(zz))  deallocate(zz)
    if (allocated(zmg)) deallocate(zmg)
    if (allocated(zdp)) deallocate(zdp)
    if (allocated(xa))    deallocate(xa)
    if (allocated(ya))    deallocate(ya)
    if (allocated(indat)) deallocate(indat)
    if (allocated(ind)) deallocate(ind)
 end subroutine free_eq_data
  
end module mod_eq_data


!-------------------------
! Back (formerly common/back/)
!   zprob, zbandw, zbkgd, xint0
!-------------------------
module mod_med_result
  use mod_kinds, only: dp
  use mod_eq_data, only: nn
  implicit none
  private
  public :: zprob, zbandw, zbkgd, xint0, lambdas, &
        alloc_med_result, free_med_result

  real(dp), allocatable :: zprob(:), zbandw(:), zbkgd(:),lambdas(:)
  real(dp) :: xint0 = 0.0_dp

contains
  subroutine alloc_med_result()
    if (nn.lt. 1) then
      write(*,*) "ERROR: alloc_back: sizes not set"
      stop
    end if
    if (.not. allocated(zprob))  allocate(zprob(nn))
    if (.not. allocated(zbandw)) allocate(zbandw(nn))
    if (.not. allocated(zbkgd))  allocate(zbkgd(nn))
    if (.not. allocated(lambdas))  allocate(lambdas(nn))
    zprob = 0.0_dp; zbandw = 0.0_dp; zbkgd = 0.0_dp
    lambdas= 0.0_dp
  end subroutine alloc_med_result

  subroutine free_med_result()
    if (allocated(zprob))  deallocate(zprob)
    if (allocated(zbandw)) deallocate(zbandw)
    if (allocated(zbkgd))  deallocate(zbkgd)
    if (allocated(lambdas))  deallocate(lambdas)
  end subroutine free_med_result
  
end module mod_med_result


!-------------------------
! Smoothing (formerly common/smoothing/)
!   bwm, bwh, npp
!-------------------------
module mod_smoothing
  use mod_kinds, only: dp
  implicit none
  public
  real(dp) :: bwm = 0.0_dp, bwh = 0.0_dp
  integer  :: npp = 0
end module mod_smoothing

!-------------------------
! Model (formerly common/model/)
!   b(9), x(9), n
!-------------------------
module mod_model
  use mod_kinds, only: dp
  implicit none
  public

  real(dp) :: b(9) = 0.0_dp
  real(dp) :: x(9) = 0.0_dp
  integer  :: n    = 0

end module mod_model
