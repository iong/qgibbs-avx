module xyz
contains
subroutine dump_xyz(r, fname, molname)
    implicit none
    real*8, intent(in) :: r(:,:)
    character(len=*), intent(in) :: fname, molname
    integer :: j, N

    N = size(r, 2)
    open(33,file=trim(fname),status='REPLACE')
    write(33,'(I10)') N
    write(33,"(A)") trim(molname)
    do j=1,N
        write(33,"('Ne ', 3(F14.8, ' '))") r(:, j)
    end do
    close(33)
end subroutine

  
subroutine load_xyz(fname, r)
    implicit none
    real*8, intent(out),allocatable :: r(:,:)
    character(LEN=*), intent(in) :: fname
    character(256):: waste
    integer :: j, N

    open(33,file=fname,STATUS='OLD')
    read(33,*) N
    if (.not.allocated(r)) then
        allocate (r(3,N))
    end if

!    if (N /= size(r, 2)) then
!        write (*,*) "N /= size(r, 2)", N, size(r, 2)
!        stop
!    end if
 
   read(33, *)
    read(33, *) (waste, r(:, j), j=1,N)
    close(33)
end subroutine
end module xyz
