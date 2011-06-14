SUBROUTINE vgw0(Q0, BL_, TAUMAX, W, WX)
    use omp_lib
    IMPLICIT NONE
    REAL*8, intent(in) :: Q0(:,:), TAUMAX, BL_
    REAL*8, intent(out) :: W
    real*8, intent(out), optional :: WX(:,:)
    real*8 :: LOGDET
    real*8 :: DT, next_stop
    integer :: i, j, chunk_size
    logical :: mm = .FALSE.

    double precision, allocatable :: Y(:), RWORK(:), YP(:), ATOL(:)
    integer, allocatable :: IWORK(:)

    integer :: NEQ, IPAR(10), ITOL, ITASK, IOPT, MF, ISTATE, LRW, LIW
    double precision :: RPAR(10), RTOL

    Natom = size(Q0, 2)
    BL = BL_

    NEQ = 9*Natom + 1
    if (present(WX)) NEQ = NEQ + 12*Natom

    LRW = 20 + 16*NEQ
    LIW = 30

    allocate(Y(NEQ), YP(NEQ), ATOL(NEQ), RWORK(LRW), IWORK(LIW))

    ITOL=2
    RTOL=0
    ATOL(1:3*Natom) = vgw_atol(1)
    ATOL(3*Natom+1:3*Natom+6*Natom)=vgw_atol(2)
    ATOL(NEQ-nthr + 1 : NEQ) = vgw_atol(3)
    ITASK=1
    ISTATE=1
    IOPT = 1
    MF=10
    IWORK=0

    IWORK(6) = 5000 !MXSTEP

    RWORK(5)=dt0
    RWORK(6)=dtmax
    RWORK(7)=dtmin

    call interaction_lists(Q0) !Determine which particles
    
    T = 0
    y = 0d0
    y(1:3*Natom) = reshape(Q0, (/ 3*NAtom /) )

    CALL DLSODE(RHSS0,NEQ,Y,T,0.5*TAUMAX,ITOL,RTOL,ATOL,ITASK,ISTATE,IOPT,&
        RWORK,LRW,IWORK,LIW,JAC,MF)

    LOGDET=0d0
    DO j=1,Natom
        LOGDET = LOGDET + LOG( DETM_S(y(3*Natom + 6*j - 5 : 3*Natom + 6*j)) )
    ENDDO


    W=-(1/TAUMAX)*(2.0*y(NEQ) - 0.5*LOGDET)
    deallocate(y, yp, RWORK, IWORK, ATOL)
END SUBROUTINE

subroutine JAC()
end subroutine
! vim:et
