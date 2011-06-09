SUBROUTINE vgw0spfm(Q0, BL_, TAUMAX, W)
    IMPLICIT NONE
    double precision, intent(in) :: Q0(:,:), TAUMAX, BL_
    double precision, intent(out) :: W
    double precision :: logdetg
    double precision :: DT, next_stop
    integer :: j, I1, I2, J2, info, s
    logical :: mm = .FALSE.

    double precision, allocatable :: Y(:), RWORK(:), YP(:), ATOL(:)
    integer, allocatable :: IWORK(:)

    integer :: NEQ, IPAR(10), ITOL, ITASK, IOPT, MF, ISTATE, LRW, LIW
    double precision :: RPAR(10), RTOL

    BL = BL_

    call interaction_lists(Q0)
    call hessian(Q0)

    nnzbmax = 2*sum(nnb) + Natom

    if (size(Gb, 3) < nnzbmax) then
        deallocate(Gb, Gcsr, UXY, UXYr, Gbja, Grja, GPb)
        allocate(Gb(3,3,nnzbmax), Gcsr(9*nnzbmax), UXY(3,3,nnzbmax), &
            UXYr(9*nnzbmax), Gbja(nnzbmax), Grja(9*nnzbmax), &
            GPb(3,3,nnzbmax))
    end if

    call init_sparse_pattern(Q0)


    NEQ = 3*Natom + 9*nnzb + 1

    !TAUMIN=5d-4

    !T=TAUMIN
    !call init_gaussians(Q0, T)

    !next_stop = 0.5d0*TAUMAX
    !do
    !    DT = 1d-2*sqrt(T)
    !    if (T+DT > next_stop) then
    !        DT = next_stop - T
    !        T = next_stop
    !    else
    !        T = T + DT
    !    end if
    !    call RHSSspfm(DT, mm)
    !    if (T == next_stop) exit
    !end do

    LRW = 20 + 16*NEQ
    LIW = 30
    allocate(Y(NEQ), YP(NEQ), ATOL(NEQ), RWORK(LRW), IWORK(LIW))

    T=0
    y=0d0
    y(1:3*Natom) = reshape(Q0, (/ 3*Natom /) )

    !call pack_y(Q, Gb(:,:,1:nnzb), gama, y)
    !call RHSSspFM(NEQ, 0d0, Y, YP)


    ITOL=2
    RTOL=0
    ATOL(1:3*Natom) = 1d-2
    ATOL(3*Natom+1:3*Natom+9*nnzb)=1d-4
    ATOL(3*Natom+9*nnzb+1) = 1000
    ITASK=1
    ISTATE=1
    IOPT = 1
    MF=10

    RWORK(5)=5d-4
    RWORK(6)=2d-3
    RWORK(7)=1d-5

    !CALL DVODE(RHSSspFM,NEQ,Y,T,0.5*TAUMAX,ITOL,RTOL,ATOL,ITASK,ISTATE,IOPT,&
    !    RWORK,LRW,IWORK,LIW,JAC,MF,RPAR,IPAR)
    CALL DLSODE(RHSSspFM,NEQ,Y,T,0.5*TAUMAX,ITOL,RTOL,ATOL(1),ITASK,ISTATE,IOPT,&
        RWORK,LRW,IWORK,LIW,JAC,MF)

    write (*,*) IWORK(11), 'steps,', IWORK(12), ' RHSS calls', 'last_dt =', RWORK(12)

    call unpack_y(y, Q, Gb(:,:,1:nnzb), gama)

    !rewind(50)
    !write (50, '(F16.8)') YP
    !rewind(51)
    !call RHSSspFM(NEQ, 0d0, Y, YP)
    !write (51, '(F16.8)') YP

    deallocate(y, yp, RWORK, IWORK, ATOL)

    gama = gama * Natom

    logdetg = cholmod_logdet(C_LOC(Gb), C_LOC(Gbia), C_LOC(Gbja), Natom)

    W=-(1/TAUMAX)*(2.0*gama - 0.5*logdetg)! - 3.0*Natom*log(2.0*sqrt(M_PI)))
    !write (*,*) gama
END SUBROUTINE


subroutine init_gaussians(q0, tau)
    REAL*8, intent(in) :: Q0(:,:), tau
    integer :: i, k
    
    call Upot_tau0(Q0)

    gama = -tau*U/real(Natom)
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


    NNB = 0
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
    nnzb = GBptr - 1
end subroutine

subroutine JAC()
end subroutine

function logdetg_gemm() result(logdetg)
    double precision :: logdetg
    double precision :: GF(3*Natom, 3*Natom)
    integer :: s, j, info

    call bsrdense(Gb, GBia, Gbja, GF)
    call dpotrf('U', 3*Natom, GF, 3*Natom, info)

    logdetg = 0.0
    s = 1
    DO j=1,3*Natom
        if (GF(j,j) < 0d0 ) s = -s
        logdetg = logdetg + 2.0*LOG(ABS(GF(j,j)))
    ENDDO

    if (s<0) then
        write(*,*) '||G|| is negative!'
        stop
    end if
end function

subroutine hessian(Q)
    double precision, intent(in) :: Q(:,:)
    double precision :: UXYf(3*Natom, 3*Natom), UX0(3), UXY0(3,3), U12, Zq(3), &
        qZq, Q12(3)
    integer :: I1, I2, J2, IG, J

    integer, save :: kkk=0

    U = 0; UX = 0; UXYf = 0
    do I1=1,Natom-1

        DO J2=1,NNB(I1)
            I2 = NBIDX(J2, I1)

            Q12 = Q(:, I1) - Q(:, I2)
            Q12 = min_image(Q12, bl)

            UX0 = 0d0; UXY0 = 0d0
            DO IG=1,NGAUSS ! BEGIN SUMMATION OVER GAUSSIANS
                Zq = LJA(IG) * Q12 ! R = -2.0*Zq
                qZq = LJA(IG)*sum(Q12**2)

                U12 = EXP(-qZq)*LJC(IG)
                U = U + U12

                UX0 = UX0 - 2d0*U12*Zq
                do J=1,3
                    UXY0(:,J) = UXY0(:,J) + 2d0*U12*2d0*Zq*Zq(J)
                    UXY0(J,J) = UXY0(J,J) - 2d0*U12 * LJA(IG)
                end do
            end do ! IG

            UX(:, I1) = UX(:, I1) + UX0
            UX(:, I2) = UX(:, I2) - UX0

            UXYf(3*I1-2 : 3*I1, 3*I1-2 : 3*I1) = UXYf(3*I1-2 : 3*I1, 3*I1-2 : 3*I1) + UXY0
            UXYf(3*I2-2 : 3*I2, 3*I2-2 : 3*I2) = UXYf(3*I2-2 : 3*I2, 3*I2-2 : 3*I2) + UXY0
            UXYf(3*I1-2 : 3*I1, 3*I2-2 : 3*I2) = -UXY0
            UXYf(3*I2-2 : 3*I2, 3*I1-2 : 3*I1) = -UXY0
        end do ! I2
    end do ! I1

    write(39+kkk, '(F16.8)') UXYf
    kkk=kkk+1
end subroutine
