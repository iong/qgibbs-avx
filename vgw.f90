module vgw
    use utils
    use iso_c_binding
    implicit none
    private
    public :: vgwinit, vgw0, vgwcleanup
    
    integer :: Natom, Nmax, maxthreads
    real(c_double) :: BL, rfullmatsq
    real(c_float) :: LJA(16), LJC(16)
    integer(c_int) :: NGAUSS

    integer, allocatable :: NBIDX(:,:), NNB(:)
    integer(c_int) :: nnbmax
    
    real*8 :: U, TRUXXG, UXXlrc, Ulrc
    real*8, allocatable :: UPV(:,:), UPM(:,:)
    
    real(c_float) :: invmass, RC, mass, dt0, dtmax, dtmin, vgw_atol(3)
    logical :: dlsode_done=.FALSE., rhss_done, pbc=.FALSE.
!$omp threadprivate(rhss_done)

    
    real(RP), allocatable, dimension(:) :: qZq, expav, v0, DETA, invDETAG
    real(RP), allocatable, dimension(:,:) :: Zq, GC, A, AG, Z

    interface
        subroutine gaussian_average_avx_init(Natom, nnb, nbidx, nnbmax, LJA, LJC, NGAUSS, bl) bind(c)
            use iso_c_binding
            integer(c_int), value :: Natom, nnbmax, NGAUSS
            integer(c_int) :: nnb(*), nbidx(*)
            real(c_float) :: LJA(*), LJC(*)
            real(c_double), value :: BL
        end subroutine
        subroutine gaussian_average_avx(y, U, UPV, UPM) bind(c)
            use iso_c_binding
            real(c_double), intent(in) :: y(*)
            real(c_double), intent(out) :: U, UPV(*), UPM(*)
        end subroutine
        subroutine gaussian_average_avx_cleanup() bind(c)
        end subroutine
    end interface

contains


subroutine vgwinit(Nmax_, species, M, rcutoff)
!$  use omp_lib
    implicit none
    integer, intent(in) :: Nmax_
    character(*), intent(in) :: species
    real*8, intent(in), optional :: M, rcutoff
    integer :: nmax_threads = 1

    Nmax = Nmax_
!$  nmax_threads = omp_get_max_threads()
    allocate(NNB(Nmax), NBIDX(Nmax,Nmax), upv(3, Nmax*nmax_threads), &
        upm(6, Nmax*nmax_threads))
include 'species.f90'
end subroutine


subroutine vgwcleanup()
    deallocate(NNB, NBIDX, UPV, UPM)
end subroutine

include 'interaction_lists.f90'
include 'vgw0.f90'
include 'rhss0.f90'

end module vgw

