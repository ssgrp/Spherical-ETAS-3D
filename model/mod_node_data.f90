
module mod_node_data
  use mod_kinds,only: dp,i8
  use mod_eq_data, only: nn
  use mod_mpi_state
  use mpi_f08
  use iso_c_binding
  implicit none
  private
  public :: hav_deltad, xis, idx_ut, inv_idx_ut,idx_lt, inv_idx_lt, &
            allocate_node_data, free_node_data, &
            win_pair1, win_pair2
  
  type(mpi_win):: win_pair1, win_pair2
  real(dp),pointer:: hav_deltad(:), xis(:)
contains
    
    integer(kind=8) function  idx_ut(i,j) result(k)
     integer,intent(in):: i, j
      k= (i-1_8)*(2_8*nn-i)/2_8+(j-i)
    end function idx_ut

        
    integer(kind=8) function  idx_lt(i,j) result(k)
     integer,intent(in):: i, j
      k= (i-1_8)*(i-2_8)/2_8+j
    end function idx_lt

    subroutine inv_idx_lt(k,i,j)
      integer(i8),intent(in):: k
      integer::i,j
      real(dp):: ri

      ri = (1.0_dp+ sqrt(1.0_dp+8.0_dp * k))/2.0_dp
      i = int(ceiling(ri))
      j= k-(i-1)*(i-2)/2
      
    end subroutine inv_idx_lt
    
  pure subroutine inv_idx_ut(k, i, j)
    integer(kind=8), intent(in) :: k
    integer, intent(out) :: i, j
    integer(kind=8) :: lo, hi, mid, s, len,n

     n =nn
     if (k < 1_8 .or. k > n*(n-1_8)/2_8) error stop "inv_idx_ut: k out of range"
     lo = 1_8
     hi = n-1_8

     do while (lo <= hi)
         mid = (lo + hi)/2_8
         s   = (mid-1_8) * (2_8*n - mid) / 2_8
         len = n - mid
         if (k < s + 1_8) then
           hi = mid - 1_8
         else if (k > s + len) then
           lo = mid + 1_8
         else
           i = int(mid)
           j = int(mid + (k - (s + 1_8)) + 1_8)
           return
         end if
     end do
     error stop "inv_idx_ut: logic error"
  end subroutine inv_idx_ut

    subroutine allocate_node_data()
      type(c_ptr)::baseptr1, baseptr2, qptr1, qptr2
      integer(kind=MPI_ADDRESS_KIND) :: qsize, winsize
      integer :: disp_unit
      
      integer ::ierr
      integer(i8)::npairs

      npairs= nn*(nn-1)/2
      
      disp_unit = 1   ! byte addressing（最简单）

      if (shmrank == 0) then
         winsize = int(npairs, MPI_ADDRESS_KIND) * c_sizeof(0.0_dp) 
      else 
         winsize = 0_MPI_ADDRESS_KIND
      end if

      call MPI_Win_allocate_shared(winsize, disp_unit, MPI_INFO_NULL, &
           comm_shm, baseptr1, win_pair1, ierr)

      call mpi_win_shared_query(win_pair1, 0, qsize, disp_unit, qptr1, ierr)

      call c_f_pointer(qptr1, hav_deltad, [npairs])
      
      call MPI_Win_allocate_shared(winsize, disp_unit, MPI_INFO_NULL, &
           comm_shm, baseptr2, win_pair2, ierr)
      
      call mpi_win_shared_query(win_pair2, 0, qsize, disp_unit, qptr2, ierr)
      
      call c_f_pointer(qptr2, xis, [npairs])

      if(shmrank ==0) then
         
         xis=0.0_dp
         hav_deltad = 0.0
      endif
      call MPI_win_sync(win_pair1, ierr)
      call mpi_win_sync(win_pair2, ierr)
      call mpi_barrier(comm_shm, ierr)
      
    end subroutine allocate_node_data
    
    subroutine free_node_data()
      integer:: ierr      
      call mpi_win_free(win_pair1, ierr)
      call mpi_win_free(win_pair2, ierr)
    end subroutine free_node_data
    
   

end module mod_node_data

