module vgw
    use utils
    use iso_c_binding
    implicit none
    private
    public :: vgwinit, vgw0, vgwcleanup
    
    integer :: Natom, Nmax, maxthreads
    real*8 :: BL, rfullmatsq
    real*8, dimension(10) :: LJA, LJC
    integer :: NGAUSS
    integer, allocatable :: NBIDX(:,:), NNB(:)
    
    real*8 :: U, TRUXXG
    real*8, allocatable :: UPV(:,:), UPM(:,:)
    
    real*8 :: invmass, RC, mass, dt0, dtmax, dtmin, vgw_atol(3)
    logical :: finished
    integer :: tid=0, nthr=1, thread_start, thread_stop, nnbmax
!$OMP THREADPRIVATE(tid, thread_start, thread_stop, nnbmax)

    
    real(RP), allocatable, dimension(:) :: DETA, invDETAG, qZq, expav, v0
    real(RP), allocatable, dimension(:, :) :: Zq
    real(RP), allocatable, dimension(:, :) :: GC, A, AG, Z

contains

subroutine vgwinit(Nmax_, species, M, rcutoff, massx)
    implicit none
    integer, intent(in) :: Nmax_
    character(*), intent(in) :: species
    real*8, intent(in), optional :: M, rcutoff, massx

    Nmax = Nmax_
    allocate(NNB(Nmax), NBIDX(Nmax,Nmax), upv(3, Nmax), &
        upm(6, Nmax))
    
    
include 'species.f90'
end subroutine


subroutine vgwcleanup()
    deallocate(NNB, NBIDX, UPV, UPM)
end subroutine

subroutine unpack_g(y, g)
    double precision, intent(in) :: y(:)
    double precision, intent(out) :: g(:,:,:)
    integer :: i, k

    i=0
    do k=3*Natom+1,9*Natom,6
        i=i+1
        g(:,1,i) = y(k : k+2)
        g(1, 2, i) = y(k+1)
        g(2:3,2,i) = y(k+3 : k+4)
        g(1, 3, i) = y(k+2)
        g(2, 3, i) = y(k+4)
        g(3, 3, i) = y(k+5)
    end do
end subroutine

include 'interaction_lists.f90'
include 'vgw0.f90'
!include 'vgw1.f90'
include 'rhss0.f90'
!include 'rhss1.f90'

end module vgw

