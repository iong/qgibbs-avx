!INCLUDE 'mkl_pardiso.f90'
module vgwspfm
    use utils
    use iso_c_binding
!    use mkl_pardiso
    implicit none
    private
    public :: vgwspfminit, vgw0spfm,vgwspfmcleanup
    
    integer :: Natom, Nmax
    real*8 :: BL, rfullmatsq = 36d0
    real*8, dimension(10) :: LJA, LJC
    integer :: NGAUSS
    integer, allocatable :: NBIDX(:,:), NNB(:),  &
            Gbia(:), Gbja(:), fmdiag(:)
    
    integer(C_INT), allocatable, target :: Gria(:), Grja(:)
    real(C_DOUBLE), allocatable, target :: Gcsr(:)
    
    real*8 :: T, gama, gamap, U
    real*8, allocatable :: Q(:,:), Gb(:,:,:), Gbdiag(:,:,:), &
            UXY(:,:,:), UXYdiag(:,:,:), UXYr(:), &
            UX(:,:), GPb(:,:,:),  QP(:,:), GUT(:,:)
    
    real*8 :: invmass, RC, TAUMIN, mass
    logical :: finished
    integer :: nnbmax

    !integer*8 :: pt(64) = 0
    !double precision :: dparm(64)
    !integer :: solver, mtype, iparm(64)


!    TYPE(MKL_PARDISO_HANDLE), ALLOCATABLE  :: pt(:)
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

!    solver = 0
!    mtype = -2
!    call pardisoinit(pt, mtype, solver, iparm, dparm, ierr)
!    iparm(3) = 1
!    if (ierr /= 0) then
!        IF (ierr == -10 ) WRITE(*,*) 'No license file found'
!        IF (ierr == -11 ) WRITE(*,*) 'License is expired'
!        IF (ierr == -12 ) WRITE(*,*) 'Wrong username or hostname'
!        STOP
!    end if


end subroutine


subroutine vgwspfmcleanup()
    deallocate(NNB, NBIDX, Gbia, Gbja, &
        Gria, Grja, fmdiag, &
        Q, Gb, Gbdiag, Gcsr, &
        UXY, UXYdiag, UXYr, &
        UX, QP, GPb)
end subroutine


subroutine init_gaussians(q0, tau)
    REAL*8, intent(in) :: Q0(:,:), tau
    integer :: i, k
    
    call Upot_tau0(Q0)

    gama = -tau*U
    Q = Q0

    Gb = 0
    do i=1,Natom
        do k=1,3
            Gb(k,k,fmdiag(i)) = tau*invmass
        end do
    end do
end subroutine

subroutine Upot_tau0(Q)
    IMPLICIT NONE
    REAL*8, intent(in) :: Q(:,:)
    INTEGER  I,J,N
    real*8 :: rsq,BL2,QIJ(3)

    N = size(Q, 2)
    BL2=BL/2.0

    U=0d0

    DO I=1,N-1
        DO J=1,NNB(I)
                qij = Q(:,I) - Q(:,NBIDX(J, I))
                rsq = sum(min_image(qij, BL)**2)
                U = U + sum(LJC(1:NGAUSS)*EXP(-LJA(1:NGAUSS)*rsq))
        ENDDO
    ENDDO
end subroutine Upot_tau0

subroutine interaction_lists(Q)
    implicit none
    real*8, intent(in) :: Q(:,:)
    integer :: N,I,J, NN
    real*8 rsq,rc2,qij(3)

    N = size(Q, 2)
    rc2=rc**2


!$OMP DO SCHEDULE(DYNAMIC) PRIVATE(I, J, QIJ, RSQ)
    do I=1,N-1
        NN = 0
        do J=I+1,N
            qij=Q(:,I)-Q(:,J)
            rsq = sum(min_image(qij, BL)**2)
            if(rsq <= rc2) then
                NN = NN + 1
                NBIDX(NN, I) = J
            endif
        enddo
        NNB(i) = NN
    enddo
!$OMP END DO
    nnbmax = maxval(nnb)
end subroutine interaction_lists


subroutine init_sparse_pattern(Q0)
    real*8, intent(in) :: Q0(:,:)
    integer :: I, J, Gbptr
    double precision :: qij(3), rsq

    Gbptr = 1
!$OMP DO SCHEDULE(DYNAMIC) PRIVATE(I, J, QIJ, RSQ)
    do I=1,Natom
        Gbia(I) = Gbptr
        do J=1,Natom
            qij=Q0(:,I)-Q0(:,J)
            rsq = sum(min_image(qij, BL)**2)
            if (rsq <= rfullmatsq) then
                Gbja(Gbptr) = J
                if (I==J) then
                    fmdiag(I) = Gbptr
                end if
                Gbptr = Gbptr + 1
            end if
        enddo
    enddo
!$OMP END DO
    Gbia(Natom + 1) = Gbptr
end subroutine


include 'vgw0spfm.f90'
include 'rhssspfm.f90'

end module vgwspfm
