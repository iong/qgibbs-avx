module vgwfm
    use utils
    implicit none
    private
    public :: vgwfminit, vgw0fm,vgwfmcleanup, classical_Utot
    
    integer :: Natom, Nmax
    real*8 :: BL, rfullmatsq = 36d0
    real*8, dimension(10) :: LJA, LJC
    integer :: NGAUSS
    
    real*8 :: T, gama, gamap, U
    real*8, allocatable :: Q(:), G(:,:), QP(:), GP(:,:), GU(:,:)
    real*8, allocatable :: UX(:), UXY(:,:)
    
    real*8 :: invmass, RC, TAUMIN, mass
    logical :: finished
    
contains

subroutine vgwfminit(Nmax_, species, M, rcutoff, massx)
    use omp_lib
    implicit none
    integer, intent(in) :: Nmax_
    character(*), intent(in) :: species
    real*8, intent(in), optional :: M, rcutoff, massx

    Natom = Nmax_
    allocate(Q(3*Natom), G(3*Natom,3*Natom), QP(3*Natom), GP(3*Natom, 3*Natom), &
        GU(3*Natom, 3*Natom), UX(3*Natom), UXY(3*Natom, 3*Natom))
    
    
    if (species=='pH2-3g') then
        NGAUSS=3
        LJA(1:3) = (/ 0.669311, 0.199426, 0.092713/)
        LJC(1:3) = (/ 29380.898517, -303.054026, -40.574585 /)
        mass = 2.0*0.020614788876D0
        rc = 8.0
        TAUMIN=1d-4
    else if (species=='pH2-4g') then
        NGAUSS=4
        LJA(1:4) = (/ 1.038252215127D0, 0.5974039109464D0, 0.196476572277834D0, &
                    0.06668611771781D0 /)
        LJC(1:4) = (/ 96609.488289873d0, 14584.62075507514d0, -365.460614956589d0, &
                    -19.5534697800036d0 /)
        mass = 2.0*0.020614788876D0
        rc = 8.0
        TAUMIN=1d-4
    else if (species == 'LJ') then
        NGAUSS = 3
        LJA(1:3) = (/ 6.65, 0.79, 2.6 /)
        LJC(1:3) = (/ 1840d0, -1.48d0, -23.2d0 /)
        mass = 1.0
        rc = 2.5
        taumin=1d-4
    end if
    
    if (present(M)) then
        mass = M
    end if
    if (present(rcutoff)) then
        rc = rcutoff
    end if

    if (present(massx)) then
        mass = massx
    end if

    invmass = 1.0/mass
end subroutine

subroutine vgwfmcleanup()
    deallocate(Q, G, QP, GP, GU, UX, UXY)
end subroutine

subroutine unpack_y(y, Q, G, gama)
    implicit none
    double precision, intent(out) :: Q(:), G(:,:), gama
    double precision, intent(in) :: y(:)

    integer :: ypos, i

    Q = y(1:3*Natom)
    
    ypos = 3*Natom+1
    do i = 1,3*Natom
    	G(i:3*Natom, i) = y(ypos : ypos + 3*Natom - i)
    	G(i, i:3*Natom) = y(ypos : ypos + 3*Natom - i)
    	ypos = ypos + 3*Natom - i + 1
    end do
    
    gama = y(ypos)
end subroutine

subroutine pack_y(Q, G, gama, y)
    implicit none
    double precision, intent(in) :: Q(:), G(:,:), gama
    double precision, intent(out) :: y(:)
    
    integer :: ypos, i

    y(1:3*Natom) = Q
    
    ypos = 3*Natom+1
    do i = 1,3*Natom
    	y(ypos : ypos + 3*Natom - i) = G(i:3*Natom, i)
    	ypos = ypos + 3*Natom - i + 1
    end do
    
    y(ypos) = gama
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
