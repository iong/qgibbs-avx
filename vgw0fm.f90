SUBROUTINE vgw0fm(Q0, BL_, TAUMAX,W)
    use omp_lib
    IMPLICIT NONE
    REAL*8, intent(in) :: Q0(:,:), TAUMAX, BL_
    REAL*8, intent(out) :: W
    real*8 :: LOGDET
    real*8 :: DT, next_stop
    integer :: j, info
    logical :: mm = .FALSE.

    double precision, allocatable :: Y(:), RWORK(:), YP(:)
    integer, allocatable :: IWORK(:)

    integer :: NEQ, IPAR(10), ITOL, ITASK, IOPT, MF, ISTATE, LRW, LIW
    double precision :: RPAR(10), RTOL, ATOL

    Natom = size(Q0, 2)
    BL = BL_

    NEQ = 3*Natom + (9*Natom**2 + 3*Natom/2) + 1

    call init_gaussians(Q0, TAUMIN)

    T = TAUMIN
    !next_stop = 0.5d0*TAUMAX
    !do
    !    DT = 1d-2*sqrt(T)
    !    if (T+DT > next_stop) then
    !        DT = next_stop - T
    !        T = next_stop
    !    else
    !        T = T + DT
    !    end if
    !    call RHSSfm(DT, mm)
    !    if (T == next_stop) exit
    !end do

    LRW = 20 + 16*NEQ
    LIW = 30
    allocate(Y(NEQ), YP(NEQ), RWORK(LRW), IWORK(LIW))

    call pack_y(Q, G, gama, y)
    !call RHSSFM(NEQ, 0d0, Y, YP)

    ITOL=1
    RTOL=1d-3
    ATOL=1d-4
    ITASK=1
    ISTATE=1
    IOPT = 0
    MF = 10
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

    deallocate(y, yp, RWORK, IWORK)

    GU = G
    call dpotrf('U', 3*Natom, GU, 3*Natom, info)
    LOGDET=0.0
    DO j=1,3*Natom
        LOGDET = LOGDET + LOG(ABS( GU(j,j) ))
    ENDDO
    LOGDET = 2d0* LOGDET

    W=-(1/TAUMAX)*(2.0*gama - 0.5*LOGDET)! - 3.0*Natom*log(2.0*sqrt(M_PI)))
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
    real*8 :: rsq,QIJ(3)

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

subroutine JAC()
end subroutine

