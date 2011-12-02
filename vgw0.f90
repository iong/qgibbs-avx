SUBROUTINE vgw0(Q0, BL_, beta,Ueff, rt, bl)
    IMPLICIT NONE
    double precision, intent(in) :: Q0(:,:), beta, BL_
    double precision, intent(out) :: Ueff
    double precision, intent(out), optional :: rt
    logical, intent(in), optional :: scale2bl
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
    if (present(scale2bl)) then
        y(1:3*Natom) = bl*reshape(Q0, (/ 3*Natom /) )
    else
        y(1:3*Natom) = reshape(Q0, (/ 3*Natom /) )
    end if

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
    write (*,*) IWORK(11), 'steps,', IWORK(12), ' RHSS calls, logdet =', logdet

    deallocate(y, yp, RWORK, IWORK, ATOL)

    Havg = -(logrho(2) - logrho(1)) / dbeta
    if (present(rt)) then
        rt = (stop_time - start_time) / real(ncalls)
     end if
END SUBROUTINE

subroutine JAC()
end subroutine
