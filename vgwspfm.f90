module vgwspfm
    use utils
    use iso_c_binding
    implicit none
    private
    public :: vgwspfminit, vgw0spfm,vgwspfmcleanup
    
    integer :: Natom, Nmax, nnzbmax, nnzb
    real*8 :: BL, rfullmatsq = 5.5d0**2
    real*8, dimension(10) :: LJA, LJC
    integer :: NGAUSS
    integer, allocatable :: NBIDX(:,:), NNB(:),  &
            Gria(:), Grja(:), fmdiag(:)
    
    integer(C_INT), allocatable, target :: Gbia(:), Gbja(:)
    real(C_DOUBLE), allocatable, target :: Gb(:,:,:)

    real*8 :: T, gama, gamap, U
    real*8, allocatable :: Q(:,:), Gcsr(:), Gbdiag(:,:,:), &
            UXY(:,:,:), UXYdiag(:,:,:), UXYr(:), &
            UX(:,:), GPb(:,:,:),  QP(:,:)

    real*8, allocatable, target :: GUT(:,:)

    real*8 :: invmass, RC, TAUMIN, mass, dt0, dtmax, dtmin, vgw_atol(3)
    logical :: finished
    integer :: nnbmax

    interface
        real(C_DOUBLE) function cholmod_logdet(G, ia, ja, N) BIND(C)
            use, intrinsic :: iso_c_binding
            type(C_PTR), value :: G, ia, ja
            integer(C_INT), value :: N
        end function

        subroutine init_cholmod(N, nnz, x, i, p) BIND(C)
            use, intrinsic :: iso_c_binding
            integer, value :: N, nnz
            type(c_ptr) :: x, i, p
        end subroutine

        subroutine destroy_cholmod() BIND(C)
        end subroutine
    end interface


contains

subroutine vgwspfminit(Np, species, M, rcutoff, massx)
    use omp_lib
    implicit none
    integer, intent(in) :: Np
    character(*), intent(in) :: species
    real*8, intent(in), optional :: M, rcutoff, massx

    Natom = Np
    allocate(NNB(Natom), NBIDX(Natom,Natom), Gbia(Natom+1), Gbja(Natom), &
        Gria(3*Natom+1), Grja(Natom), fmdiag(Natom), &
        Q(3,Natom), Gb(3,3,4), Gbdiag(3,3,Natom), Gcsr(4), &
        UXY(3,3,4), UXYdiag(3,3,Natom), UXYr(4), &
        UX(3,Natom), QP(3,Natom), GPb(3,3,4), GUT(3*Natom,3*Natom))
    
    
    include 'species.f90'
end subroutine


subroutine vgwspfmcleanup()
    deallocate(NNB, NBIDX, Gbia, Gbja, &
        Gria, Grja, fmdiag, &
        Q, Gb, Gbdiag, Gcsr, &
        UXY, UXYdiag, UXYr, &
        UX, QP, GPb, GUT)
end subroutine

subroutine unpack_y(y, Q, Gb, gama)
    implicit none
    double precision, intent(out) :: Q(:,:), Gb(:,:,:), gama
    double precision, intent(in) :: y(:)

    Q = reshape(y(1:3*Natom), (/3, Natom/) )
    Gb(:,:,1:nnzb) = reshape(y(3*Natom+1 : 3*Natom + 9*nnzb), (/3,3,nnzb/) )
    gama = y(3*Natom + 9*nnzb + 1)
end subroutine

subroutine pack_y(Q, Gb, gama, y)
    implicit none
    double precision, intent(in) :: Q(:,:), Gb(:,:,:), gama
    double precision, intent(out) :: y(:)

    y(1:3*Natom) = reshape(Q, (/3*Natom/) )
    y(3*Natom+1 : 3*Natom + 9*nnzb) = reshape(Gb(:,:,1:nnzb), (/9*nnzb/) )
    y(3*Natom + 9*nnzb + 1) = gama
end subroutine

include 'vgw0spfm.f90'


include 'rhssspfm.f90'

end module vgwspfm
