SUBROUTINE vgw0fm(Q0, BL_, TAUMAX,W)
    use omp_lib
    IMPLICIT NONE
    REAL*8, intent(in) :: Q0(:,:), TAUMAX, BL_
    REAL*8, intent(out) :: W
    real*8 :: LOGDET
    real*8 :: DT, next_stop
    integer :: j, info, nlg
    logical :: mm = .FALSE.

    double precision, allocatable :: Y(:), RWORK(:), YP(:), ATOL(:)
    integer, allocatable :: IWORK(:)

    integer :: NEQ, IPAR(10), ITOL, ITASK, IOPT, MF, ISTATE, LRW, LIW
    double precision :: RPAR(10), RTOL

    Natom = size(Q0, 2)
    BL = BL_

    nlg = (9*Natom**2 + 3*Natom)/2
    NEQ = 3*Natom + nlg + 1

    call init_gaussians(Q0, TAUMIN)

    T = TAUMIN

    LRW = 20 + 16*NEQ
    LIW = 30
    allocate(Y(NEQ), YP(NEQ), ATOL(NEQ), RWORK(LRW), IWORK(LIW))

    call pack_y(Q, G, gama, y)
    !call RHSSFM(NEQ, 0d0, Y, YP)

    ITOL=2
    RTOL=0
    ATOL(1:3*Natom) = vgw_atol(1)
    ATOL(3*Natom+1:3*Natom+nlg)=vgw_atol(2)
    ATOL(3*Natom+nlg+1) = vgw_atol(3)
    ITASK=1
    ISTATE=1
    IOPT = 1
    MF=10
    IWORK=0

    IWORK(6) = 50000 !MXSTEP

    RWORK(5)=dt0
    RWORK(6)=dtmax
    RWORK(7)=dtmin
    !CALL DVODE(RHSSspFM,NEQ,Y,T,0.5*TAUMAX,ITOL,RTOL,ATOL,ITASK,ISTATE,IOPT,&
    !    RWORK,LRW,IWORK,LIW,JAC,MF,RPAR,IPAR)
    CALL DLSODE(RHSSFM,NEQ,Y,T,0.5*TAUMAX,ITOL,RTOL,ATOL,ITASK,ISTATE,IOPT,&
        RWORK,LRW,IWORK,LIW,JAC,MF)

    call unpack_y(y, Q, G, gama)
    gama = gama * real(Natom)

    !rewind(50)
    !write (50, '(F16.8)'), YP
    !rewind(51)
    !call RHSSFM(NEQ, 0d0, Y, YP)
    !write (51, '(F16.8)'), YP

    deallocate(y, yp, ATOL, RWORK, IWORK)

    GU = G
    call dpotrf('U', 3*Natom, GU, 3*Natom, info)
    LOGDET=0.0
    DO j=1,3*Natom
        LOGDET = LOGDET + LOG(ABS( GU(j,j) ))
    ENDDO
    LOGDET = 2d0* LOGDET

    W=-(1/TAUMAX)*(2.0*gama - 0.5*LOGDET - 1.5*Natom*log(4.0*M_PI))
    !write (*,*) gama
END SUBROUTINE


subroutine init_gaussians(q0, tau)
    REAL*8, intent(in) :: Q0(:,:), tau
    integer :: i
    
    call Upot_tau0(Q0)

    gama = -tau*U/real(Natom)
    Q = reshape(Q0, (/ 3*Natom /) )

    G = 0
    do i=1,3*Natom
        G(i,i) = tau*invmass
    end do
end subroutine


subroutine Upot_tau0(Q)
    IMPLICIT NONE
    REAL*8, intent(in) :: Q(:,:)
    INTEGER  I,J,N
    real*8 :: rsq, QIJ(3)

    N = size(Q, 2)

    U=0d0
    DO I=1,N-1
        DO J=I+1,N
                qij = Q(:,I) - Q(:,J)
                rsq = sum(min_image(qij, BL)**2)
                U = U + sum(LJC(1:NGAUSS)*EXP(-LJA(1:NGAUSS)*rsq))
        ENDDO
    ENDDO
end subroutine Upot_tau0


function classical_Utot(Q0, BLX)
    double precision, intent(in) :: Q0(:,:), BLX
    double precision :: classical_Utot

    BL = BLX
    call Upot_tau0(Q0)
    classical_Utot = U
end function

subroutine JAC()
end subroutine

