module vgwfm
    use utils
    implicit none
    private
    public :: vgwfminit, vgw0fm,vgwfmcleanup
    
    integer :: Natom, Nmax
    real*8 :: BL
    real*8, dimension(10) :: LJA, LJC
    integer :: NGAUSS
    
    real*8 :: T, gama, gamap, U
    real*8, allocatable :: Q(:), G(:,:), QP(:), GP(:,:), GU(:,:)
    real*8, allocatable :: UX(:), UXY(:,:)
    
    real*8 :: invmass, RC, TAUMIN, mass
    logical :: finished
    
contains

subroutine vgwfminit(Nmax_, species, M, rcutoff)
    use omp_lib
    implicit none
    integer, intent(in) :: Nmax_
    character(*), intent(in) :: species
    real*8, intent(in), optional :: M, rcutoff

    Natom = Nmax_
    allocate(Q(3*Natom), G(3*Natom,3*Natom), QP(3*Natom), GP(3*Natom, 3*Natom), &
        GU(3*Natom, 3*Natom), UX(3*Natom), UXY(3*Natom, 3*Natom))
    
    
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

subroutine vgwfmcleanup()
    deallocate(Q, G, QP, GP, GU, UX, UXY)
end subroutine

subroutine Upot_tau0(Q)
    IMPLICIT NONE
    REAL*8, intent(in) :: Q(:,:)
    INTEGER  I,J,N
    real*8 :: rsq,QIJ(3)

    N = size(Q, 2)

    U=0d0
    DO I=1,N-1
        DO J=I+1,N
                qij = Q(:,I) - Q(:,J)
                rsq = sum(min_image(qij, BL)**2)
                U = U + sum(LJC(1:NGAUSS)*EXP(-LJA(1:NGAUSS)*rsq))
        ENDDO
    ENDDO
end subroutine Upot_tau0

subroutine init_gaussians(q0, tau)
    REAL*8, intent(in) :: Q0(:,:), tau
    integer :: i
    
    call Upot_tau0(Q0)

    gama = -tau*U
    Q = reshape(Q0, (/ 3*Natom /) )

    G = 0
    do i=1,3*Natom
        G(i,i) = tau*invmass
    end do
end subroutine

subroutine get_gaussian(Qtau, Gtau, gamatau)
    double precision, intent(out), optional :: Qtau(:,:), Gtau(:,:), gamatau
    
    if (present(Qtau)) then
        Qtau = reshape(Q, (/ 3, Natom /) )
    end if
    
    if (present(Gtau)) then
        Gtau = G
    end if
    
    if (present(gamatau)) then
        gamatau = gama
    end if
end subroutine

include 'vgw0fm.f90'
include 'rhssfm.f90'

end module vgwfm