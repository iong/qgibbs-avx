module vgwfm
    use utils
    implicit none
    private
    public :: vgwfminit, vgw0fm,vgwfmcleanup, classical_Utot
    
    integer :: Natom, Nmax, nlg
    real*8 :: BL, rfullmatsq = 36d0
    real*8, dimension(10) :: LJA, LJC
    integer :: NGAUSS
    
    real*8 :: T, gama, gamap, U
    real*8, allocatable :: Q(:), G(:,:), QP(:), GP(:,:), GU(:,:)
    real*8, allocatable :: UX(:), UXY(:,:)
    
    real*8 :: invmass, RC, mass, dt0, dtmax, dtmin, vgw_atol(3)
    logical :: finished
    
contains

subroutine vgwfminit(Nmax_, species, M, rcutoff, massx)
    use omp_lib
    implicit none
    integer, intent(in) :: Nmax_
    character(*), intent(in) :: species
    real*8, intent(in), optional :: M, rcutoff, massx

    Natom = Nmax_
    allocate(Q(3*Natom), G(3*Natom,3*Natom), QP(3*Natom), GP(3*Natom, 3*Natom), &
        GU(3*Natom, 3*Natom), UX(3*Natom), UXY(3*Natom, 3*Natom))
    
include 'species.f90'
end subroutine

subroutine vgwfmcleanup()
    deallocate(Q, G, QP, GP, GU, UX, UXY)
end subroutine

subroutine unpack_y(y, Q, G, gama)
    implicit none
    double precision, intent(out) :: Q(:), G(:,:), gama
    double precision, intent(in) :: y(:)

    integer :: ypos, i

    Q = y(1:3*Natom)
    
    ypos = 3*Natom+1
    do i = 1,3*Natom
    	G(i:3*Natom, i) = y(ypos : ypos + 3*Natom - i)
    	G(i, i:3*Natom) = y(ypos : ypos + 3*Natom - i)
    	ypos = ypos + 3*Natom - i + 1
    end do
    
    gama = y(ypos)
end subroutine

subroutine pack_y(Q, G, gama, y)
    implicit none
    double precision, intent(in) :: Q(:), G(:,:), gama
    double precision, intent(out) :: y(:)
    
    integer :: ypos, i

    y(1:3*Natom) = Q
    
    ypos = 3*Natom+1
    do i = 1,3*Natom
    	y(ypos : ypos + 3*Natom - i) = G(i:3*Natom, i)
    	ypos = ypos + 3*Natom - i + 1
    end do
    
    y(ypos) = gama
end subroutine


subroutine get_gaussian(Qtau, Gtau, gamatau)
    double precision, intent(out), optional :: Qtau(:,:), Gtau(:,:), gamatau
    
    if (present(Qtau)) then
        Qtau = reshape(Q, (/ 3, Natom /) )
    end if
    
    if (present(Gtau)) then
        Gtau = G
    end if
    
    if (present(gamatau)) then
        gamatau = gama
    end if
end subroutine


SUBROUTINE vgw0fm(Q0, BL_, beta, Ueff, rt)
    IMPLICIT NONE
    double precision, intent(in) :: Q0(:,:), beta, BL_
    double precision, intent(out) :: Ueff
    double precision, intent(out), optional :: rt
    real*8 :: logdet, logrho, TSTOP, start_time, stop_time
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

    call cpu_time(start_time)
    TSTOP = 0.5d0*beta
    CALL DLSODE(RHSSFM,NEQ,Y,T,TSTOP,ITOL,RTOL,ATOL,ITASK,ISTATE,IOPT,&
        RWORK,LRW,IWORK,LIW,JAC,MF)

    call unpack_y(y, Q, G, gama)
    gama = gama * real(Natom)

     write(79,'(F16.10)') GU

    GU = G
    call dpotrf('U', 3*Natom, GU, 3*Natom, info)
    logdet=0.0
    DO j=1,3*Natom
        logdet = logdet + LOG(ABS( GU(j,j) ))
    ENDDO
    logdet = 2d0* logdet

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


SUBROUTINE vgw0fmgs(Q0, BL_, beta,Havg, rt)
    IMPLICIT NONE
    double precision, intent(in) :: Q0(:,:), beta, BL_
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
    dbeta = 0.1*beta

    call cpu_time(start_time)
    do i=1,2
        TSTOP = 0.5d0*(beta - (2-i)*dbeta)
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


SUBROUTINE RHSSFM(NEQ, T, Y, YP)
!    use omp_lib
    IMPLICIT NONE
    integer, intent(in) :: NEQ!, IPAR(:)
    double precision, intent(in) :: T, Y(NEQ)!, RPAR(:)
    double precision, intent(out) :: YP(NEQ)
    INTEGER :: J,I1,I2,IG, I1_3, I2_3
    REAL*8 AG(3,3), &
            DETA,DETAG,QZQ,U12, &
            G12(3,3),A(3,3), &
            Zq(3), Z(3,3),Q12(3)
    real*8 :: UXY0(3,3), UX0(3)

    if (y(3*Natom+1) == 0d0) then
        call rhss_zero_time(NEQ, y, yp)
        return
    end if

    call unpack_y(y, Q, G, gama)

!    do I1=1,Natom-1
	call unpack_y(y, Q, G, gama) 

    U = 0; UX = 0; UXY = 0;
    do I1=1,Natom-1
        I1_3 = 3*(I1-1) + 1
        DO I2=I1+1,Natom
            I2_3 = 3*(I2-1) + 1
            Q12 = Q(I1_3:I1_3+2) - Q(I2_3:I2_3+2)
            Q12 = min_image(Q12, bl)
            G12=G(I1_3 : I1_3 + 2, I1_3 : I1_3 + 2) &
                + G(I2_3 : I2_3 + 2, I2_3 : I2_3 + 2) &
                - G(I2_3 : I2_3 + 2, I1_3 : I1_3 + 2) &
                - G(I1_3 : I1_3 + 2, I2_3 : I2_3 + 2)

            call detminvm(G12, DETA, A)
            DETA = 1.0d0/DETA

            UX0 = 0d0; UXY0 = 0d0
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

                U12 = SQRT(DETA/DETAG)*EXP(-qZq)*LJC(IG)
                U = U + U12

                UX0 = UX0 - 2d0*U12*Zq
                do J=1,3
                    UXY0(:,J) = UXY0(:,J) + 2d0*U12*(2d0*Zq*Zq(J) - Z(:,J))
                end do
            end do ! IG
            
! Avoid load and store as much as possbile. Store now, process as a stream later. Much faster.
            UX(I1_3 : I1_3 + 2) = UX(I1_3 : I1_3 + 2) + UX0
            UX(I2_3 : I2_3 + 2) = UX(I2_3 : I2_3 + 2) - UX0
            UXY(I1_3 : I1_3 + 2, I1_3 : I1_3 + 2) = UXY(I1_3 : I1_3 + 2, I1_3 : I1_3 + 2) + UXY0
            UXY(I2_3 : I2_3 + 2, I2_3 : I2_3 + 2) = UXY(I2_3 : I2_3 + 2, I2_3 : I2_3 + 2) + UXY0
            
            !if (sum(q12**2) <= rfullmatsq) then
                UXY(I1_3 : I1_3 + 2, I2_3 : I2_3 + 2) = -UXY0
                UXY(I2_3 : I2_3 + 2, I1_3 : I1_3 + 2) = -UXY0
            !end if
        end do ! I2
    end do ! I1

    
    !QP = -matmul(G, UX)
    call dsymm('L', 'L', 3*Natom, 1, -1d0, G, 3*Natom, UX, 3*Natom, &
        0d0, QP, 3*Natom)

    !GU = matmul(G, UXY)
    call dsymm('L', 'L', 3*Natom, 3*Natom, 1d0, G, 3*Natom, UXY, 3*Natom, &
        0d0, GU, 3*Natom)

    !GP = -matmul(GU, G)
    call dsymm('R', 'L', 3*Natom, 3*Natom, -1d0, G, 3*Natom, GU, 3*Natom, &
        0d0, GP, 3*Natom)

    do J=1,3*Natom
        GP(J,J) = GP(J,J) + invmass
    end do

!    do I1=1,Natom-1
!        I1_3 = 3*(I1-1) + 1
!        DO I2=I1+1,Natom
!            I2_3 = 3*(I2-1) + 1
!            Q12 = Q(I1_3:I1_3+2) - Q(I2_3:I2_3+2)
!            Q12 = min_image(Q12, bl)
!            if (sum(q12**2) > rfullmatsq) then
!                GP(I1_3 : I1_3 + 2, I2_3 : I2_3 + 2) = 0
!                GP(I2_3 : I2_3 + 2, I1_3 : I1_3 + 2) = 0
!           end if
!        end do
!    end do
    
    gamap = -(0.25d0 * sum(UXY*G) + U)/real(Natom)

	call pack_y(QP, GP, gamap, yp)
END SUBROUTINE RHSSFM

subroutine rhss_zero_time(NEQ, y, yp)
    integer, intent(in) :: NEQ
    double precision, intent(in) :: y(:)
    double precision, intent(out) :: yp(:)

    double precision :: qij(3), qi(3), qj(3), rsq
    integer :: i, j

    yp = 0d0

    j = 3*Natom + 1
    do i=1,3*Natom
        yp(j) = invmass
        j = j + 3*Natom - i + 1
    end do

    if (NEQ > 3*Natom + nlg + 1) then
        j = 3*Natom + nlg
        do i=1,Natom
            yp(j + 9*(i - 1) + 1) = 1d0
            yp(j + 9*(i - 1) + 5) = 1d0
            yp(j + 9*i          ) = 1d0
        end do
    end if

    U=0d0

    DO I=1,Natom-1
        qi = y(3*I-2:3*I)
        DO J=I+1,Natom
                qj = y(3*J-2 : 3*J)
                qij = qi - qj
                rsq = sum(min_image(qij, BL)**2)
                U = U + sum(LJC(1:NGAUSS)*EXP(-LJA(1:NGAUSS)*rsq))
        ENDDO
    ENDDO

    yp(NEQ) = -U/real(Natom)
end subroutine
end module vgwfm
