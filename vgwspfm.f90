module vgwspfm
    use utils
    use iso_c_binding
    implicit none
    private
    public :: vgwspfminit, vgw0spfm,vgwspfmcleanup
    
    integer :: Natom, Nmax, nnzmax, nnz, nrows, n3rows16
    real*8 :: BL, rfullmatsq = 5.5d0**2
    real*8, dimension(10) :: LJA, LJC
    integer :: NGAUSS
    integer, allocatable :: NBIDX(:,:), NNB(:),  &
            Gbia(:), Gbja(:), Giia(:), Gbiia(:)
    
    integer(C_INT), allocatable, target :: Gia(:), Gja(:)
    real(C_DOUBLE), allocatable, target :: G(:)

    real*8 :: gama, gamap, U
    real*8, allocatable :: Q(:,:), Gb(:,:,:), Gbdiag(:,:,:), &
            UXYdiag(:,:,:), UX(:,:)

    real*8, allocatable, target :: GU(:,:), UXYf(:,:)

    real*8 :: invmass, RC, mass, dt0, dtmax, dtmin, vgw_atol(3)
    logical :: finished
    integer :: nnbmax

    interface
        real(C_DOUBLE) function cholmod_logdet(G, ia, ja, N) BIND(C)
            use, intrinsic :: iso_c_binding
            type(C_PTR), value :: G, ia, ja
            integer(C_INT), value :: N
        end function
    end interface


contains

subroutine vgwspfminit(Np, species, M, rcutoff, massx)
    use omp_lib
    implicit none
    integer, intent(in) :: Np
    character(*), intent(in) :: species
    real*8, intent(in), optional :: M, rcutoff, massx

    Natom = Np
    n3rows16 = ((3*Natom+1)/2)*2 ! align first dimension on 16byte boundary

    allocate(NNB(Natom), NBIDX(Natom,Natom), Gbia(Natom+1), Gbja(1), &
        Gia(3*Natom+1), Gja(1), Giia(3*Natom), Gbiia(Natom), &
        Gb(3,3,1), Gbdiag(3,3,Natom), &
        UXYdiag(3,3,Natom), &
        UX(3,Natom), UXYf(n3rows16,3*Natom), &
        GU(n3rows16,3*Natom))
    
    
    include 'species.f90'
end subroutine


subroutine vgwspfmcleanup()
    deallocate(NNB, NBIDX, Gbia, Gbja, &
        Gia, Gja, Giia, Gbiia, &
        Gb, Gbdiag, &
        UXYdiag, &
        UX, UXYf, GU)
end subroutine

include 'vgw0spfm.f90'


include 'rhssspfm.f90'

end module vgwspfm
