SUBROUTINE vgw0spfm(Q0, BL_, beta, Ueff, rt)
    IMPLICIT NONE
    double precision, intent(in) :: Q0(:,:), beta, BL_
    double precision, intent(out) :: Ueff
    double precision, intent(out), optional :: rt
    real*8 :: LOGDET, logrho, TSTOP, start_time, stop_time, T
    integer :: i, ncalls

    double precision, allocatable :: Y(:), RWORK(:), YP(:), ATOL(:)
    integer, allocatable :: IWORK(:)

    integer :: NEQ, ITOL, ITASK, IOPT, MF, ISTATE, LRW, LIW
    double precision ::  RTOL

    Natom = size(Q0, 2)
    BL = BL_

    call interaction_lists(Q0)

    if (size(Gja) < nnz) then

        deallocate(Gb, Gbja, Gja)

        allocate(Gb(3,3,nnz/9), Gbja(nnz/9), Gja(nnz))
    end if

    call init_sparse_pattern(Q0)
    call test_ia_ja(Gia, Gja)

    NEQ = 3*Natom + nnz + 1


    LRW = 20 + 16*NEQ
    LIW = 30
    allocate(Y(NEQ), YP(NEQ), ATOL(NEQ), RWORK(LRW), IWORK(LIW))

    ITOL=2
    RTOL=1d-5
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

    call cpu_time(start_time)

    T=0
    TSTOP = 0.5d0*beta
    CALL DLSODE(RHSSspFM,NEQ,Y,T,TSTOP,ITOL,RTOL,ATOL,ITASK,ISTATE,IOPT,&
        RWORK,LRW,IWORK,LIW,JAC,MF)

    gama = Y(NEQ) * real(Natom)

    logdet = cholmod_logdet(C_LOC(y(3*Natom+1)), C_LOC(Gia), C_LOC(Gja), &
        3*Natom)
    if (logdet == 1234.1234d0) then
        print *, 'The G matrix is not positive definite!'
        open(10,file='nonposdef.xyz', status='REPLACE')
        write(10, '(I5/2F16.8/("Ne",3F14.8))') Natom, bl, beta, q0
        close(10)
        stop
    end if
        

    logrho = 2.0*gama - 0.5*logdet - 1.5*Natom*log(4.0*M_PI)

    call cpu_time(stop_time)

    ncalls = IWORK(12)
    write (*,*) IWORK(11), 'steps,', IWORK(12), ' RHSS calls, logdet =', logdet

    deallocate(y, yp, RWORK, IWORK, ATOL)

    Ueff = -logrho/beta
    if (present(rt)) then
        rt = (stop_time - start_time) / real(ncalls)
     end if
END SUBROUTINE


SUBROUTINE vgw0spfmgs(Q0, BL_, TAUMAX, Havg, rt, logfd, cfgname)
    IMPLICIT NONE
    double precision, intent(in) :: Q0(:,:), TAUMAX, BL_
    double precision, intent(out) :: Havg
    double precision, intent(out), optional :: rt
    integer, intent(in), optional :: logfd
    character(*), intent(in), optional :: cfgname
    real*8 :: LOGDET, logrho(2), dbeta, TSTOP, start_time, stop_time, T
    integer :: i, ncalls, nsteps

    double precision, allocatable :: Y(:), RWORK(:), YP(:), ATOL(:)
    integer, allocatable :: IWORK(:)

    integer :: NEQ, ITOL, ITASK, IOPT, MF, ISTATE, LRW, LIW
    double precision ::  RTOL

    Natom = size(Q0, 2)
    BL = BL_

    call interaction_lists(Q0)

    if (size(Gja) < nnz) then

        deallocate(Gb, Gbja, Gja)

        allocate(Gb(3,3,nnz/9), Gbja(nnz/9), Gja(nnz))
    end if

    call init_sparse_pattern(Q0)
    call test_ia_ja(Gia, Gja)

    NEQ = 3*Natom + nnz + 1


    LRW = 20 + 16*NEQ
    LIW = 30
    allocate(Y(NEQ), YP(NEQ), ATOL(NEQ), RWORK(LRW), IWORK(LIW))

    ITOL=2
    RTOL=1d-5
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
        if (present(cfgname)) then
            open(57,file=trim(cfgname),status='REPLACE')
            write (57,*) Natom
            write(57, *) 'vgw intermediate'
            write(57, '(("Ne ", 3F13.8))') y(1:3*Natom)
            close(57)
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
    nnz = 0
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
            if (rsq <= rfullmatsq) then
                nnz = nnz + 1
            end if
        enddo
        NNB(i) = NN
    enddo
!$OMP END DO
    nnbmax = maxval(nnb)
    nnz = 9*(2*nnz + Natom)
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
    if (nnz /= p - 1) then
        write (*,*) 'nnz error!', nnz, p-1
        stop
    end if
end subroutine

subroutine JAC()
end subroutine

subroutine test_ia_ja(ia, ja)
    integer, intent(in) :: ia(:), ja(:)
    integer,allocatable :: P(:,:)
    integer :: i, ptr, N
    
    allocate(P(size(ia)-1, size(ia)-1))

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

    deallocate(P)
end  subroutine
