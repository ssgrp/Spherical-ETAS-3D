

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

