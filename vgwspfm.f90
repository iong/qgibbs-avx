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

    real*8 :: invmass, RC, TAUMIN, mass
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

subroutine vgwspfminit(Np, species, M, rcutoff)
    use omp_lib
    implicit none
    integer, intent(in) :: Np
    character(*), intent(in) :: species
    real*8, intent(in), optional :: M, rcutoff

    Natom = Np
    allocate(NNB(Natom), NBIDX(Natom,Natom), Gbia(Natom+1), Gbja(Natom), &
        Gria(3*Natom+1), Grja(Natom), fmdiag(Natom), &
        Q(3,Natom), Gb(3,3,4), Gbdiag(3,3,Natom), Gcsr(4), &
        UXY(3,3,4), UXYdiag(3,3,Natom), UXYr(4), &
        UX(3,Natom), QP(3,Natom), GPb(3,3,4), GUT(3*Natom,3*Natom))
    
    
    if (species=='pH2-3g') then
        NGAUSS=3
        LJA(1:3) = (/ 0.669311, 0.199426, 0.092713/)
        LJC(1:3) = (/ 29380.898517, -303.054026, -40.574585 /)
        mass = 2.0
        rc = 8.0
        TAUMIN=1d-4
    else if (species=='pH2-4g') then
        NGAUSS=4
        LJA(1:4) = (/ 1.038252215127D0, 0.5974039109464D0, 0.196476572277834D0, &
                    0.06668611771781D0 /)
        LJC(1:4) = (/ 96609.488289873d0, 14584.62075507514d0, -365.460614956589d0, &
                    -19.5534697800036d0 /)
        mass = 2.0
        rc = 8.0
        TAUMIN=1d-4
    end if
    
    if (present(M)) then
        mass = M
    end if
    if (present(rcutoff)) then
        rc = rcutoff
    end if
        
    mass = mass*0.020614788876D0
    invmass = 1.0/mass
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
