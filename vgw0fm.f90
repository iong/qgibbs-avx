SUBROUTINE vgw0fm(Q0, BL_, TAUMAX,Havg, rt)
    IMPLICIT NONE
    double precision, intent(in) :: Q0(:,:), TAUMAX, BL_
    double precision, intent(out) :: Havg
    double precision, intent(out), optional :: rt
    real*8 :: logdet, logrho(2), dbeta, TSTOP, start_time, stop_time
    integer :: i, j, ncalls, info

    double precision, allocatable :: Y(:), RWORK(:), YP(:), ATOL(:)
    integer, allocatable :: IWORK(:)

    integer :: NEQ, ITOL, ITASK, IOPT, MF, ISTATE, LRW, LIW
    double precision :: RTOL

    Natom = size(Q0, 2)
    BL = BL_

    nlg = (9*Natom**2 + 3*Natom)/2
    NEQ = 3*Natom + nlg + 1

    LRW = 20 + 16*NEQ
    LIW = 30
    allocate(Y(NEQ), YP(NEQ), ATOL(NEQ), RWORK(LRW), IWORK(LIW))

    ITOL=2
    RTOL=0
    ATOL(1:3*Natom) = vgw_atol(1)
    ATOL(3*Natom+1:3*Natom+nlg)=vgw_atol(2)
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
    dbeta = 0.1*TAUMAX

    call cpu_time(start_time)
    do i=1,2
        TSTOP = 0.5d0*(TAUMAX - (2-i)*dbeta)
        CALL DLSODE(RHSSFM,NEQ,Y,T,TSTOP,ITOL,RTOL,ATOL,ITASK,ISTATE,IOPT,&
            RWORK,LRW,IWORK,LIW,JAC,MF)

        call unpack_y(y, Q, G, gama)
        gama = gama * real(Natom)

        GU = G
        call dpotrf('U', 3*Natom, GU, 3*Natom, info)
        logdet=0.0
        DO j=1,3*Natom
            logdet = logdet + LOG(ABS( GU(j,j) ))
        ENDDO
        logdet = 2d0* logdet

        logrho(i) = 2.0*gama - 0.5*logdet - 1.5*Natom*log(4.0*M_PI)
    end do
    call cpu_time(stop_time)
    ncalls = IWORK(12)
	write (*,*) IWORK(11), 'steps,', IWORK(12), ' RHSS calls, logdet =', logdet

	deallocate(y, yp, RWORK, IWORK, ATOL)

    Havg = -(logrho(2) - logrho(1)) / dbeta
    if (present(rt)) then
        rt = (stop_time - start_time) / real(ncalls)
     end if
END SUBROUTINE

subroutine JAC()
end subroutine

function classical_Utot(Q0, BLX) result(U)
    double precision, intent(in) :: Q0(:,:), BLX
    double precision :: U
    INTEGER  I,J,N
    real*8 :: rsq, QIJ(3)

    N = size(Q0, 2)

    U=0d0
    DO I=1,N-1
        DO J=I+1,N
                qij = Q0(:,I) - Q0(:,J)
                rsq = sum(min_image(qij, BLX)**2)
                U = U + sum(LJC(1:NGAUSS)*EXP(-LJA(1:NGAUSS)*rsq))
        ENDDO
    ENDDO
end function
