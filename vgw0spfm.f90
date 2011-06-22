SUBROUTINE vgw0spfm(Q0, BL_, TAUMAX, Havg, rt, logfd)
    IMPLICIT NONE
    double precision, intent(in) :: Q0(:,:), TAUMAX, BL_
    double precision, intent(out) :: Havg
    double precision, intent(out), optional :: rt
    integer, intent(in), optional :: logfd
    real*8 :: LOGDET, logrho(2), dbeta, TSTOP, start_time, stop_time, T
    integer :: i, ncalls, nsteps

    double precision, allocatable :: Y(:), RWORK(:), YP(:), ATOL(:)
    integer, allocatable :: IWORK(:)

    integer :: NEQ, ITOL, ITASK, IOPT, MF, ISTATE, LRW, LIW
    double precision ::  RTOL

    Natom = size(Q0, 2)
    BL = BL_

    call interaction_lists(Q0)

    nnzmax = 9*(2*sum(nnb) + Natom)

    if (size(Gja) < nnzmax) then

        deallocate(Gb, Gbja, Gja)

        allocate(Gb(3,3,nnzmax/9), Gbja(nnzmax/9), Gja(nnzmax))
    end if

    call init_sparse_pattern(Q0)
    call test_ia_ja(Gia, Gja)

    NEQ = 3*Natom + nnz + 1


    LRW = 20 + 16*NEQ
    LIW = 30
    allocate(Y(NEQ), YP(NEQ), ATOL(NEQ), RWORK(LRW), IWORK(LIW))

    ITOL=2
    RTOL=0
    ATOL(1:3*Natom) = vgw_atol(1)
    ATOL(3*Natom+1:3*Natom+nnz)=vgw_atol(2)
    ATOL(NEQ) = vgw_atol(3)
    ITASK=1
    ISTATE=1
    IOPT = 1
    MF=10
    IWORK=0

    IWORK(6) = 50000 ! MXSTEP

    RWORK(5)=dt0
    RWORK(6)=dtmax
    RWORK(7)=dtmin

    y=0d0
    y(1:3*Natom) = reshape(Q0, (/ 3*Natom /) )
    dbeta = min(0.1*0.5*TAUMAX, 1d0)
    nsteps = ceiling(0.5*TAUMAX/dbeta) 
    dbeta = 0.5*TAUMAX/real(nsteps)

    call cpu_time(start_time)

    T=0
    logrho=0
    do i=1,nsteps
        logrho(1) = logrho(2)

        TSTOP = T  + dbeta
        CALL DLSODE(RHSSspFM,NEQ,Y,T,TSTOP,ITOL,RTOL,ATOL,ITASK,ISTATE,IOPT,&
            RWORK,LRW,IWORK,LIW,JAC,MF)

        gama = Y(NEQ) * real(Natom)

        logdet = cholmod_logdet(C_LOC(y(3*Natom+1)), C_LOC(Gia), &
            C_LOC(Gja), 3*Natom)

        logrho(2) = 2.0*gama - 0.5*logdet - 1.5*Natom*log(4.0*M_PI)

        Havg = -(logrho(2) - logrho(1)) *0.5d0/ dbeta
        if (present(logfd)) then
            write (logfd,*) 2.0*T, Havg
        end if
    end do

    call cpu_time(stop_time)
    ncalls = IWORK(12)
    write (*,*) IWORK(11), 'steps,', IWORK(12), ' RHSS calls, logdet =', logdet

    deallocate(y, yp, RWORK, IWORK, ATOL)

    if (present(rt)) then
        rt = (stop_time - start_time) / real(ncalls)
     end if
END SUBROUTINE


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
    integer :: I, J, p, p0, ib, jb, nz
    double precision :: qij(3), rsq

    p = 1
    do ib=1,Natom
        i = 3*(ib - 1) + 1

        p0 = p

        do jb=1,Natom
            j = 3*(jb - 1) + 1

            qij=Q0(:,ib)-Q0(:,jb)
            rsq = sum(min_image(qij, BL)**2)

            if (rsq <= rfullmatsq) then
                Gja(p : p+2) = (/ j, j+1, j+2/)

                if (ib == jb) then
                    Giia(i) = p
                    Gbiia(ib) = (p0-1)/9 + (p-p0) / 3 + 1
                end if

                p = p + 3
            end if
        enddo
        nz = p - p0

        Gia(i : i+2) = (/ p0, p, p + nz/)

        Gja(p : p + nz - 1) = Gja(p0 : p - 1) 
        Gja(p + nz : p + 2*nz - 1) = Gja(p0 : p - 1) 

        p = p + 2*nz

        Giia(i+1) = Giia(i)   + nz + 1
        Giia(i+2) = Giia(i+1) + nz + 1
    enddo
    Gia(3*Natom + 1) = p
    nnz = p - 1
end subroutine

subroutine JAC()
end subroutine

subroutine hessian(Q)
    double precision, intent(in) :: Q(:,:)
    double precision :: UXYf(3*Natom, 3*Natom), UX0(3), UXY0(3,3), U12, Zq(3), &
        qZq, Q12(3)
    integer :: I1, I2, J2, IG, J

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
end subroutine

subroutine test_ia_ja(ia, ja)
    integer, intent(in) :: ia(:), ja(:)
    integer :: P(size(ia)-1, size(ia)-1), i, ptr, N

    N = size(ia) - 1

    P = 0
    do i=1,N
        do ptr=ia(i), ia(i+1)-1
            P(i, ja(ptr)) = 1
        end do
    end do

    if (sum(P-transpose(P)) /= 0) then
        write (*,*) 'The pattern is not symmetric!'
        stop
    end if
end  subroutine
