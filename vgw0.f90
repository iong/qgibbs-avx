SUBROUTINE vgw0(Q0, BL_, TAUMAX,Havg)
    IMPLICIT NONE
    double precision, intent(in) :: Q0(:,:), TAUMAX, BL_
    double precision, intent(out) :: Havg
    real*8 :: LOGDET, logrho(2), dbeta, TSTOP
    integer :: i, j

    double precision, allocatable :: Y(:), RWORK(:), YP(:), ATOL(:)
    integer, allocatable :: IWORK(:)

    integer :: NEQ, IPAR(10), ITOL, ITASK, IOPT, MF, ISTATE, LRW, LIW
    double precision :: RPAR(10), RTOL

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
    dbeta = 0.1*TAUMAX

    do i=1,2
        TSTOP = 0.5d0*(TAUMAX - (2-i)*dbeta)
        CALL DLSODE(RHSS0,NEQ,Y,T,TSTOP,ITOL,RTOL,ATOL,ITASK,ISTATE,IOPT,&
        RWORK,LRW,IWORK,LIW,JAC,MF)

    	LOGDET=0d0
    	DO j=1,Natom
        	LOGDET = LOGDET + LOG( DETM_S(y(3*Natom + 6*j - 5 : 3*Natom + 6*j)) )
    	ENDDO

		logrho(i) = 2.0*Natom*y(NEQ) - 0.5*LOGDET - 1.5*Natom*log(4.0*M_PI)
    end do

	write (*,*) IWORK(11), 'steps,', IWORK(12), ' RHSS calls, logdet =', logdet

	deallocate(y, yp, RWORK, IWORK, ATOL)

    !W=-(1/TAUMAX)*logrho(2)
    Havg = -(logrho(2) - logrho(1)) / dbeta
END SUBROUTINE

subroutine JAC()
end subroutine
