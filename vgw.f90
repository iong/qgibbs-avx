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
        subroutine gaussian_average_acc_init(Natom, nnb, nbidx, nnbmax, LJA, LJC, NGAUSS, bl) bind(c)
            use iso_c_binding
            integer(c_int), value :: Natom, nnbmax, NGAUSS
            integer(c_int) :: nnb(*), nbidx(*)
            real(c_float) :: LJA(*), LJC(*)
            real(c_double), value :: BL
        end subroutine
        subroutine gaussian_average_acc(y, U, UPV, UPM) bind(c)
            use iso_c_binding
            real(c_double), intent(in) :: y(*)
            real(c_double), intent(out) :: U, UPV(*), UPM(*)
        end subroutine
        subroutine gaussian_average_acc_cleanup() bind(c)
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

subroutine interaction_lists(r)
    implicit none
    real*8, intent(in) :: r(3,Natom)
    integer :: I,J, NN
    real*8 :: rsq,rcsq, dr(3), bl2
    

    rcsq=rc**2
    bl2=0.5*bl
!$OMP DO SCHEDULE(DYNAMIC) PRIVATE(J, RSQ, dr, NN)
    do I=1,Natom-1
        NN = 0
        do J=I+1,Natom
            dr = r(:,j) - r(:,i)

            where (abs(dr) > bl2)
                dr  = dr - sign(bl, dr)
            end where

            rsq = sum(dr**2)
 
            if(rsq <= rcsq) then
                NN = NN + 1
                NBIDX(NN, I) = J
            endif
        enddo
        NNB(i) = NN
    enddo
!$OMP ENDDO

!$omp master
    NNB(Natom) = 0
    nnbmax = maxval(nnb)
    nnbmax = (nnbmax/8 + 1) * 8
!$omp end master
end subroutine interaction_lists

SUBROUTINE vgw0(Q0, BL_, beta,Ueff)
!$  use omp_lib
    IMPLICIT NONE
    double precision, intent(in) :: Q0(:,:), beta, BL_
    double precision, intent(out) :: Ueff
    real*8 :: LOGDET, logrho, T, TSTOP
    integer :: i, j, ncalls

    double precision, allocatable :: Y(:), RWORK(:), ATOL(:)
    integer, allocatable :: IWORK(:)

    integer :: tid = 0
    integer :: NEQ, ITOL, ITASK, IOPT, MF, ISTATE, LRW, LIW
    double precision :: RTOL

    Natom = size(Q0, 2)
    BL = BL_

    NEQ = 9*Natom + 1

    LRW = 20 + 16*NEQ
    LIW = 30
    allocate(Y(NEQ), ATOL(NEQ), RWORK(LRW), IWORK(LIW))

    ITOL=2
    RTOL=0
    ATOL(1:3*Natom) = vgw_atol(1)
    ATOL(3*Natom+1:3*Natom+6*Natom)=vgw_atol(2)
    ! tolerance for Qnkp and gamakp
    ATOL(NEQ) = vgw_atol(3)
    ITASK=1
    ISTATE=1
    IOPT = 1
    MF=10


    IWORK=0
    IWORK(6) = 50000 !MXSTEP

    RWORK(5)=dt0
    RWORK(6)=dtmax
    RWORK(7)=dtmin
    
    T = 0
    TSTOP = 0.5d0*beta

    y(1 : 3*Natom) = reshape(q0, (/3*Natom/))
    y(3*Natom+1:) = 0d0

    call presort_ppc(y(1:3*Natom), 8)
    if (pbc) then
        Ulrc = Ulrc * Natom**2 / BL**3
        UXXlrc = UXXlrc * Natom / BL**3
    else
        Ulrc = 0.0
        UXXlrc = 0.0
    end if

    dlsode_done=.FALSE.
    !print *, 'OAK',   omp_get_thread_num(), omp_get_wtime()
    tid = 0
!$omp parallel private(tid)

    call interaction_lists(y(1:3*Natom))

!$    tid = omp_get_thread_num()

    rhss_done = .FALSE.

    if (tid == 0) then
        !print *, 'master entry', omp_get_thread_num(), omp_get_wtime()
        call gaussian_average_acc_init(Natom, nnb, nbidx, nnbmax, LJA, LJC, NGAUSS, bl)
        CALL DLSODE(RHSS0,NEQ,Y,T,TSTOP,ITOL,RTOL,ATOL,ITASK,ISTATE,IOPT,&
            RWORK,LRW,IWORK,LIW,JAC,MF)
        dlsode_done = .TRUE.
        call rhss0(NEQ,0d0,Y,Y)
        !print *, 'master exit', omp_get_thread_num(), omp_get_wtime()
    else
        !print *, 'Entering tid =', omp_get_thread_num(), dlsode_done, omp_get_wtime()
        do
            call rhss0(NEQ,0d0,Y,Y)
            if (rhss_done) then
                !print *, 'Exiting tid=', tid, rhss_done, omp_get_wtime()
                exit
            end if
        end do
    end if

!$omp end parallel
    !print *, 'if exit tid=',  omp_get_thread_num(), omp_get_wtime()

    LOGDET=0d0
    DO j=3*Natom+1,9*Natom,6
        LOGDET = LOGDET + LOG( DETM_S(y(j : j+5)) )
    ENDDO
    !print *, 'p', LOGDET, y(NEQ)

    logrho = 2.0*Natom*y(NEQ) - 0.5*LOGDET - 1.5*Natom*log(4.0*M_PI)
    ncalls = IWORK(12)
    !write (*,*) IWORK(11), 'steps,', IWORK(12), ' RHSS calls, logdet =', logdet

    deallocate(y, RWORK, IWORK, ATOL)
    call gaussian_average_acc_cleanup()

    Ueff = -logrho/beta
END SUBROUTINE


subroutine presort_ppc(r, nppc)
    use utils, only: index_sort
    implicit none
    double precision, intent(inout) :: r(3,Natom)
    integer, intent(in) :: nppc

    integer(c_size_t) :: i, N
    integer(c_int) :: idx(Natom), nunits
    real(c_double) :: z_(Natom), r_(3,Natom)
    double precision :: bu, cbl, ll(3), ur(3), sysbox

    ll = minval(r, 2)
    ur = maxval(r, 2)
    sysbox = maxval(ur - ll)

    cbl = bl

    ! cluster
    if (bl > 2*sysbox) then
        nunits = nint(2.0 * (0.75*Natom/(M_PI*nppc))**(1.0/3.0))
        cbl = sysbox
        pbc = .FALSE.
    else
        nunits = nint((real(Natom)/real(nppc))**(1.0/3.0))
        pbc = .TRUE.
    end if

    bu = cbl / nunits

    
    forall (i=1:Natom) idx(i) = i
    
    z_ = ( floor((r(1,:) - ll(1)) / bu) * nunits + floor( (r(2,:) - ll(2)) / bu ) ) * cbl + r(3,:) - ll(3)    
    
    N = Natom
    call index_sort(N, idx, z_)

    r_ = r(:,idx)
    r = r_
end subroutine presort_ppc

subroutine adjust_nbidx_for_c()
    integer :: i
    do i=1,Natom
        nbidx(1:nnb(i),i) =  nbidx(1:nnb(i),i) - 1
    end do
end subroutine


subroutine JAC()
end subroutine


SUBROUTINE RHSS0(NEQ, T, YMASTER, YPMASTER)
    use utils, only: RP
    IMPLICIT NONE
    integer, intent(in) :: NEQ
    double precision, intent(in) :: T
    double precision, intent(in) , target :: YMASTER(NEQ)
    double precision, intent(out), target :: YPMASTER(NEQ)

    double precision, pointer, save :: Y(:), YP(:)

    INTEGER :: I1, J

    real(8) :: GU(6), GUG(6), QP_(3), G(6), UPM1(6)
    !real(8), pointer :: G(:)

!$omp master
    Y => YMASTER
    YP => YPMASTER
    TRUXXG = 0d0
!$omp end master
!$omp barrier
    if (dlsode_done) then
        rhss_done=.TRUE.
        return
    end if

    if (y(3*Natom+1) == 0.0) then
        call rhss_zero_time(NEQ, y, yp)
        return
    end if

    call gaussian_average_acc(y, U, UPV, UPM)

!$omp do schedule(static) reduction(+:TRUXXG)
    do i1=1,Natom
        G = y(3*Natom + 6*I1-5 : 3*Natom + 6*I1)
        UPM1 = UPM(:,I1)

        UPM1( (/ 1, 4, 6/) ) = UPM1( (/ 1, 4, 6/) ) + UXXlrc

        TRUXXG = TRUXXG +     G(1)*UPM1(1) + 2d0*G(2)*UPM1(2) &
                        + 2d0*G(3)*UPM1(3) +     G(4)*UPM1(4) &
                        + 2d0*G(5)*UPM1(5) +     G(6)*UPM1(6)

        QP_(1) = -G(1)*UPV(1,I1) - G(2)*UPV(2,I1) - G(3)*UPV(3,I1)
        QP_(2) = -G(2)*UPV(1,I1) - G(4)*UPV(2,I1) - G(5)*UPV(3,I1)
        QP_(3) = -G(3)*UPV(1,I1) - G(5)*UPV(2,I1) - G(6)*UPV(3,I1)

        GU = mm_ss(G, UPM1)
        GUG = -mm_ss(GU, G)

        GUG(1) = GUG(1) + invmass
        GUG(4) = GUG(4) + invmass
        GUG(6) = GUG(6) + invmass

        UPM(:,I1) = UPM1
        yp(3*I1-2 : 3*I1) = QP_
        yp(3*Natom + 6*I1-5 : 3*Natom + 6*I1) = GUG
    end do
!$omp end do

!$omp master
    U = U + Ulrc
    yp(NEQ) = -(0.25d0*TRUXXG + U)/real(Natom)
!$omp end master
END SUBROUTINE RHSS0

subroutine gaussian_average(y, U, UPV, UPM)
    double precision, intent(in) :: y(:)
    double precision, intent(out) :: U, UPV(:,:), UPM(:,:)
    INTEGER :: J,I1, NN1
    real(RP), dimension(nnbmax) :: x12, y12, z12
    real(RP), dimension(nnbmax, 3) :: UX0
    real(RP), dimension(nnbmax, 6) :: UXX0


    UPM = 0.0d0
    UPV = 0.0d0

    U = 0d0
    do I1=1,Natom-1
        NN1 = NNB(I1)
        if (NN1 == 0) cycle

        x12(1:NN1) = min_image(y(          I1) - y(          nbidx(1:NN1,I1)), bl)
        y12(1:NN1) = min_image(y(  Natom + I1) - y(  Natom + nbidx(1:NN1,I1)), bl)
        z12(1:NN1) = min_image(y(2*Natom + I1) - y(2*Natom + nbidx(1:NN1,I1)), bl)

        do J=1,6
            GC(1:NN1,J) = y((J+2)*Natom + nbidx(1:NN1,I1)) + y((J+2)*Natom + I1)
        end do

        call stream_kernel(NN1, x12, y12, z12, GC, U, UX0, UXX0)

        do J=1,NN1
            UPV(:,I1) = UPV(:, I1) + UX0 (J,:)
            UPV(:,nbidx(J,I1)) = UPV(:,nbidx(J,I1)) - UX0(J,:)

            UPM(:,I1) = UPM(:,I1) + UXX0(J,:)
            UPM(:,nbidx(J,I1)) = UPM(:,nbidx(J,I1)) + UXX0(J,:)
        end do
    ENDDO ! I
end subroutine



subroutine stream_kernel(NN1, x12, y12, z12, GC, U, UX0, UXX0)
    implicit none
    integer, intent(in) :: NN1
    real(RP), intent(in) :: x12(:), y12(:), z12(:), GC(:,:)
    double precision, intent(inout) :: U
    real(RP), intent(out) :: UX0(:,:), UXX0(:,:)


    integer :: IG, J

    call pdetminvm_sg(NN1, GC, DETA, A)

    DO IG=1,NGAUSS      ! BEGIN SUMMATION OVER GAUSSIANS
        AG(1:NN1,1) = LJA(IG) + A(1:NN1,1)
        AG(1:NN1,2) =           A(1:NN1,2)
        AG(1:NN1,3) =           A(1:NN1,3)
        AG(1:NN1,4) = LJA(IG) + A(1:NN1,4)
        AG(1:NN1,5) =           A(1:NN1,5)
        AG(1:NN1,6) = LJA(IG) + A(1:NN1,6)

        call pdetminvm_sg(NN1, AG, invDETAG, Z)
        do J=1,6
            Z(1:NN1,J) = - (LJA(IG)**2) * Z(1:NN1,J)
        end do

        Z(1:NN1,1)=LJA(IG) + Z(1:NN1,1)
        Z(1:NN1,4)=LJA(IG) + Z(1:NN1,4)
        Z(1:NN1,6)=LJA(IG) + Z(1:NN1,6)

        !Zq = matmul(Z, Q12) ! R = -2.0*Zq
        ZQ(1:NN1,1) = Z(1:NN1,1)*x12(1:NN1) + Z(1:NN1,2) * y12(1:NN1) + Z(1:NN1,3)*z12(1:NN1)
        ZQ(1:NN1,2) = Z(1:NN1,2)*x12(1:NN1) + Z(1:NN1,4) * y12(1:NN1) + Z(1:NN1,5)*z12(1:NN1)
        ZQ(1:NN1,3) = Z(1:NN1,3)*x12(1:NN1) + Z(1:NN1,5) * y12(1:NN1) + Z(1:NN1,6)*z12(1:NN1)

        !qZq = dot_product(Q12, Zq) 
        qZq(1:NN1) = Zq(1:NN1,1)*x12(1:NN1) + Zq(1:NN1,2)*y12(1:NN1) + Zq(1:NN1,3)*z12(1:NN1)

        EXPAV(1:NN1)=EXP(-qZq(1:NN1)) * SQRT(DETA(1:NN1)*invDETAG(1:NN1))

        v0(1:NN1) = 2.0_RP*expav(1:NN1)*LJC(IG)
        U = U + 0.5_RP *sum(v0(1:NN1))

        if (IG == 1) then
            DO J=1,3
                UX0(1:NN1,J) = - v0(1:NN1)*Zq(1:NN1,J)
            END DO

            UXX0(1:NN1,1) = v0(1:NN1)*(2.0_RP*Zq(1:NN1,1)*Zq(1:NN1,1) - Z(1:NN1,1))
            UXX0(1:NN1,2) = v0(1:NN1)*(2.0_RP*Zq(1:NN1,2)*Zq(1:NN1,1) - Z(1:NN1,2))
            UXX0(1:NN1,3) = v0(1:NN1)*(2.0_RP*Zq(1:NN1,3)*Zq(1:NN1,1) - Z(1:NN1,3))
            UXX0(1:NN1,4) = v0(1:NN1)*(2.0_RP*Zq(1:NN1,2)*Zq(1:NN1,2) - Z(1:NN1,4))
            UXX0(1:NN1,5) = v0(1:NN1)*(2.0_RP*Zq(1:NN1,3)*Zq(1:NN1,2) - Z(1:NN1,5))
            UXX0(1:NN1,6) = v0(1:NN1)*(2.0_RP*Zq(1:NN1,3)*Zq(1:NN1,3) - Z(1:NN1,6))
            cycle
        end if
        DO J=1,3
            UX0(1:NN1,J) = UX0(1:NN1,J) - v0(1:NN1)*Zq(1:NN1,J)
        END DO

        UXX0(1:NN1,1) = UXX0(1:NN1,1) + v0(1:NN1)*(2.0_RP*Zq(1:NN1,1)*Zq(1:NN1,1) - Z(1:NN1,1))
        UXX0(1:NN1,2) = UXX0(1:NN1,2) + v0(1:NN1)*(2.0_RP*Zq(1:NN1,2)*Zq(1:NN1,1) - Z(1:NN1,2))
        UXX0(1:NN1,3) = UXX0(1:NN1,3) + v0(1:NN1)*(2.0_RP*Zq(1:NN1,3)*Zq(1:NN1,1) - Z(1:NN1,3))
        UXX0(1:NN1,4) = UXX0(1:NN1,4) + v0(1:NN1)*(2.0_RP*Zq(1:NN1,2)*Zq(1:NN1,2) - Z(1:NN1,4))
        UXX0(1:NN1,5) = UXX0(1:NN1,5) + v0(1:NN1)*(2.0_RP*Zq(1:NN1,3)*Zq(1:NN1,2) - Z(1:NN1,5))
        UXX0(1:NN1,6) = UXX0(1:NN1,6) + v0(1:NN1)*(2.0_RP*Zq(1:NN1,3)*Zq(1:NN1,3) - Z(1:NN1,6))
    ENDDO
end subroutine stream_kernel

!    11 22 33
!    21 42 53
!    31 52 63
!    22 44 55
!    32 54 65
!    33 55 66
pure function mm_ss(A, B) result(C)
    double precision, intent(in) :: A(6), B(6)
    double precision :: C(6)

    double precision :: x22, x33, x53, x55

    x22 = A(2) * B(2)
    x33 = A(3) * B(3)
    x53 = A(5) * B(3)
    x55 = A(5) * B(5)
    C(1) = A(1)*B(1) + x22 + x33
    C(2) = A(2)*B(1) + A(4)*B(2) + x53
    C(3) = A(3)*B(1) + A(5)*B(2) + A(6)*B(3)
    C(4) = x22 + A(4)*B(4) + A(5)*B(5)
    C(5) = A(3)*B(2) + A(5)*B(4) + A(6)*B(5)
    C(6) = x33 + x55 + A(6)*B(6)
end function mm_ss


subroutine rhss_zero_time(NEQ, y, yp)
    use utils, only: RP
    implicit none
    integer, intent(in) :: NEQ
    double precision, intent(in) :: y(:)
    double precision, intent(out) :: yp(:)

    double precision :: dr(3), rsq, bl2
    integer :: i, j

!$omp workshare
    yp(1:3*Natom) = 0.0d0
!$omp end workshare

!$omp master
    U=0.0d0
!$omp end master

!$omp do schedule(static)
    do i=3*Natom+1,9*Natom,6
        yp(i : i+5) = (/invmass, 0.0, 0.0, invmass, 0.0, invmass/)
    end do
!$omp end do

    bl2 = 0.5*bl

!$omp do schedule(dynamic) reduction(+:U) private(rsq, dr, j)
    DO I=1,Natom-1
        if (nnb(i) == 0) cycle
        do j=1,nnb(i)
            dr = y(3*nbidx(j,i) - 2 :  3*nbidx(j,i)) - y(3*i-2:3*i)

            where (abs(dr) > bl2)
                dr  = dr - sign(bl, dr)
            end where

            rsq = sum(dr**2)

            U = U + sum(LJC(1:NGAUSS)*EXP(-LJA(1:NGAUSS)*rsq))
        end do
    ENDDO
!$omp end do

!$omp master
    !print *, 'U =', U
    U = U + Ulrc

    yp(NEQ) = -U/real(Natom )
!$omp end master
end subroutine

end module vgw

