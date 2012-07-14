module vgwref
    use utils
    implicit none
    private
    public :: vgwrefinit, get_fx, vgwref0, vgwrefcleanup
    
    integer :: Natom, Nmax, maxthreads
    real*8 :: BL, rfullmatsq
    real*8, dimension(10) :: LJA, LJC
    integer :: NGAUSS
    integer, allocatable :: NBIDX(:,:), NNB(:)
    
    real*8 :: T, gama, gamap, U, TRUXXG
    real*8, allocatable :: Q(:,:), G(:,:,:), Qnk(:,:,:), gamak(:,:), &
                            QP(:,:), GP(:,:,:)
                            
    real*8, allocatable :: UPV(:,:), UPM(:,:,:)
    
    real*8 :: invmass, RC, mass, dt0, dtmax, dtmin, vgw_atol(3)
    logical :: finished
    integer :: tid=0, nthr=1, thread_start, thread_stop, nnbmax
!$OMP THREADPRIVATE(tid, thread_start, thread_stop, nnbmax)

contains

subroutine vgwrefinit(Nmax_, species, M, rcutoff, massx)
    implicit none
    integer, intent(in) :: Nmax_
    character(*), intent(in) :: species
    real*8, intent(in), optional :: M, rcutoff, massx

    Nmax = Nmax_
    allocate(NNB(Nmax), NBIDX(Nmax,Nmax), upv(3,Nmax), &
        upm(3,3,Nmax), g(3,3,Nmax))
    
    
include 'species.f90'
end subroutine


subroutine vgwrefcleanup()
    deallocate(NNB, NBIDX, UPV, UPM, g)
end subroutine

function get_fx() result(fx)
    real*8 :: fx(3, Natom)
    fx = - gamak(3,Natom)/T
end function


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
    NNB(N) = 0
    nnbmax = maxval(nnb)
end subroutine interaction_lists

SUBROUTINE vgwref0(Q0, BL_, beta,Ueff, rt)
    IMPLICIT NONE
    double precision, intent(in) :: Q0(:,:), beta, BL_
    double precision, intent(out) :: Ueff
    double precision, intent(out), optional :: rt
    real*8 :: LOGDET, logrho, TSTOP, start_time, stop_time
    integer :: i, j, ncalls

    double precision, allocatable :: Y(:), RWORK(:), YP(:), ATOL(:)
    integer, allocatable :: IWORK(:)

    integer :: NEQ, ITOL, ITASK, IOPT, MF, ISTATE, LRW, LIW
    double precision :: RTOL

    Natom = size(Q0, 2)
    BL = BL_

    call interaction_lists(Q0)

    NEQ = 9*Natom + 1
    !if (present(WX)) NEQ = NEQ + 12*Natom

    LRW = 20 + 16*NEQ
    LIW = 30

    allocate(Y(NEQ), YP(NEQ), ATOL(NEQ), RWORK(LRW), IWORK(LIW))

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
    y = 0d0
    y(1:3*Natom) = reshape(Q0, (/ 3*Natom /) )

    call cpu_time(start_time)
    TSTOP = 0.5d0*beta
    CALL DLSODE(RHSS0,NEQ,Y,T,TSTOP,ITOL,RTOL,ATOL,ITASK,ISTATE,IOPT,&
        RWORK,LRW,IWORK,LIW,JAC,MF)

    LOGDET=0d0
    DO j=1,Natom
        LOGDET = LOGDET + LOG( DETM_S(y(3*Natom + 6*j - 5 : 3*Natom + 6*j)) )
    ENDDO
    call cpu_time(stop_time)

    logrho = 2.0*Natom*y(NEQ) - 0.5*LOGDET - 1.5*Natom*log(4.0*M_PI)
    ncalls = IWORK(12)
    !write (*,*) IWORK(11), 'steps,', IWORK(12), ' RHSS calls, logdet =', logdet

    deallocate(y, yp, RWORK, IWORK, ATOL)

    Ueff = -logrho/beta
    if (present(rt)) then
        rt = (stop_time - start_time) / real(ncalls)
     end if
END SUBROUTINE


SUBROUTINE vgw0gs(Q0, BL_, beta,Havg, rt)
    IMPLICIT NONE
    double precision, intent(in) :: Q0(:,:), beta, BL_
    double precision, intent(out) :: Havg
    double precision, intent(out), optional :: rt
    real*8 :: LOGDET, logrho(2), dbeta, TSTOP, start_time, stop_time
    integer :: i, j, ncalls

    double precision, allocatable :: Y(:), RWORK(:), YP(:), ATOL(:)
    integer, allocatable :: IWORK(:)

    integer :: NEQ, ITOL, ITASK, IOPT, MF, ISTATE, LRW, LIW
    double precision :: RTOL

    Natom = size(Q0, 2)
    BL = BL_

    call interaction_lists(Q0)

    NEQ = 9*Natom + 1
    !if (present(WX)) NEQ = NEQ + 12*Natom

    LRW = 20 + 16*NEQ
    LIW = 30

    allocate(Y(NEQ), YP(NEQ), ATOL(NEQ), RWORK(LRW), IWORK(LIW))

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
    y = 0d0
    y(1:3*Natom) = reshape(Q0, (/ 3*Natom /) )
    dbeta = 0.1*beta

    call cpu_time(start_time)
    do i=1,2
        TSTOP = 0.5d0*(beta - (2-i)*dbeta)
        CALL DLSODE(RHSS0,NEQ,Y,T,TSTOP,ITOL,RTOL,ATOL,ITASK,ISTATE,IOPT,&
        RWORK,LRW,IWORK,LIW,JAC,MF)

        LOGDET=0d0
        DO j=1,Natom
            LOGDET = LOGDET + LOG( DETM_S(y(3*Natom + 6*j - 5 : 3*Natom + 6*j)) )
        ENDDO

        logrho(i) = 2.0*Natom*y(NEQ) - 0.5*LOGDET - 1.5*Natom*log(4.0*M_PI)
    end do
    call cpu_time(stop_time)
    ncalls = IWORK(12)
    !write (*,*) IWORK(11), 'steps,', IWORK(12), ' RHSS calls, logdet =', logdet

    deallocate(y, yp, RWORK, IWORK, ATOL)

    Havg = -(logrho(2) - logrho(1)) / dbeta
    if (present(rt)) then
        rt = (stop_time - start_time) / real(ncalls)
     end if
END SUBROUTINE

subroutine JAC()
end subroutine

SUBROUTINE RHSS0(NEQ, T, Y, YP)
    IMPLICIT NONE
    integer, intent(in) :: NEQ!, IPAR(:)
    double precision, intent(in) :: T, Y(NEQ)!, RPAR(:)
    double precision, intent(out) :: YP(NEQ)

    INTEGER :: J,I1,I2,IG, NN1, k1, k2, k3, k4
    REAL*8 AG(3,3),GU(3,3), &
            DETA,DETAG,GUG(6),QZQ,EXPAV, &
            G12(6),A(3,3), &
            Zq(3), Z(3,3),Q12(3), v0, QP_(3), Q1(3), G1(6)
    real*8 :: UXX0(3,3), UX0(3), UPV_I1(3), UPM_I1(3,3)
    real*8 :: Ulocal
    real*8, allocatable :: UX(:,:),UXX(:,:, :), GC(:,:), QC(:, :)

    !print *, 'rhss', nnbmax
    allocate(UX(3,nnbmax),UXX(3,3, nnbmax), GC(6,nnbmax), QC(3, nnbmax))
    UPM = 0d0
    UPV = 0d0
    Ulocal = 0d0

    if (y(3*Natom+1) == 0d0) then
        call rhss_zero_time(NEQ, y, yp)
        return
    end if
        

    do I1=1,Natom-1
        NN1 = NNB(I1)
        if (NN1 == 0) cycle

        q1 = y(3*I1 - 2 : 3*I1)
        G1 = y(3*Natom + 6*I1 - 5 : 3*Natom + 6*I1)
        UPV_I1 = UPV(:,I1) 
        UPM_I1 = UPM(:,:,I1) 
        
        do I2=1,NN1
            qc(:,I2) = y(3*nbidx(I2,I1) - 2 : 3*nbidx(I2,I1))
            GC(:,I2) = y(3*Natom + 6*nbidx(I2,I1) - 5 : 3*Natom + 6*nbidx(I2,I1))
        end do

        DO I2=1,NN1
            Q12 = Q1 - QC(:,I2)
            Q12 = min_image(Q12, bl)
            G12=G1+GC(:,I2)

            call detminvm_sg(G12, DETA, A)
            DETA = 1.0d0/DETA

            UX0 = 0d0; UXX0 = 0d0
            DO IG=1,NGAUSS      ! BEGIN SUMMATION OVER GAUSSIANS
                AG = A
                do J=1,3
                    AG(J,J)=LJA(IG)+AG(J,J)
                end do

                call detminvm(AG, DETAG, Z)
                Z = - LJA(IG)**2 * Z

                do J=1,3
                    Z(J,J) = Z(J,J) + LJA(IG)
                end do

                Zq = matmul(Z, Q12) ! R = -2.0*Zq
                qZq = dot_product(Q12, Zq) 

                EXPAV=SQRT(DETA/DETAG)*EXP(-qZq)
                Ulocal=Ulocal+EXPAV*LJC(IG)

                v0 = 2d0*expav*LJC(IG)

                UX0 = UX0 - v0*Zq
                do J=1,3
                    UXX0(:,J) = UXX0(:,J) + v0*(2d0*Zq*Zq(J) - Z(:,J))
                end do
            ENDDO ! IG
! Avoid load and store as much as possbile. Store now, process as a stream later. Much faster.
            UPV_I1 = UPV_I1 + UX0
            if (I2>size(UX,2)) print *, size(UX, 2), I2
            UX(:,I2) = - UX0
            UPM_I1 = UPM_I1 + UXX0
            UXX(:,:,I2) = UXX0
        ENDDO ! I2
        UPV(:,I1) = UPV_I1
        UPM(:,:,I1) = UPM_I1
        UPV(:,nbidx(1:NN1,I1)) = UPV(:,nbidx(1:NN1,I1)) + UX(:,1:NN1)
        UPM(:,:,nbidx(1:NN1,I1)) = UPM(:,:,nbidx(1:NN1,I1)) + UXX(:,:,1:NN1)
    ENDDO ! I

    call unpack_g(y, g)
    TRUXXG = sum(UPM*G)
    U = Ulocal

    do i1=1,Natom
        k1=             3*(i1 - 1) + 1
        k2 =  3*Natom + 6*(i1 - 1) + 1
        k3 =  9*Natom + 9*(i1 - 1) + 1
        k4 = 18*Natom + 3*(i1 - 1) + 1

        QP_ = - matmul(G(:,:,I1), UPV(:,I1))

        GU = matmul(G(:,:,I1), UPM(:,:,I1))
        GUG = -matmul_sgs(GU, G(:,:,I1))
        GUG(1) = GUG(1) + invmass
        GUG(4) = GUG(4) + invmass
        GUG(6) = GUG(6) + invmass

        yp(k1 : k1+2) = QP_
        yp(k2 : k2+5) = GUG

        if (NEQ > 21*Natom ) then
            yp(k3 : k3+8) = - reshape(matmul(GU, Qnk(:,:,i1)), (/ 9 /) )

            yp(k4 : k4 + 2 ) = - matmul(transpose(Qnk(:,:,i1)), UPV(:,I1))
        end if
    end do

    yp(NEQ) = -(0.25d0*TRUXXG + U)/real(Natom)
    print *, maxval(abs(yp(1:9*Natom))), yp(NEQ)
    deallocate(UX,UXX, GC, QC)
END SUBROUTINE RHSS0

function matmul_sgs(A, B) result (C)
    double precision :: A(3,3), B(3,3)
    double precision :: C(6)

    C(1:3) = matmul(A, B(:,1))
    C(4) = sum(A(2,:)*B(:,2))
    C(5) = sum(A(3,:)*B(:,2))
    C(6) = sum(A(3,:)*B(:,3))
end function

subroutine rhss_zero_time(NEQ, y, yp)
    integer, intent(in) :: NEQ
    double precision, intent(in) :: y(:)
    double precision, intent(out) :: yp(:)

    double precision :: qij(3), qi(3), qj(3), rsq
    integer :: i, j

    yp = 0d0

    do i=1,Natom
        yp(3*Natom + 6*i - 5 : 3*Natom + 6*i) = invmass*(/1d0, 0d0, 0d0, 1d0, 0d0, 1d0/)
    end do

    if (NEQ > 18*Natom) then
        do i=1,Natom
            yp(9*Natom + 9*(i - 1) + 1) = 1d0
            yp(9*Natom + 9*(i - 1) + 5) = 1d0
            yp(9*Natom + 9*i          ) = 1d0
        end do
    end if

    U=0d0

    DO I=1,Natom-1
        qi = y(3*I-2:3*I)
        DO J=1,NNB(I)
                qj = y(3*NBIDX(J,I)-2 : 3*NBIDX(J,I))
                qij = qi - qj
                rsq = sum(min_image(qij, BL)**2)
                U = U + sum(LJC(1:NGAUSS)*EXP(-LJA(1:NGAUSS)*rsq))
        ENDDO
    ENDDO

    yp(NEQ) = -U/real(Natom)
end subroutine

end module vgwref

