SUBROUTINE RHSS0(NEQ, T, Y, YP)
    use utils, only: RP
    IMPLICIT NONE
    integer, intent(in) :: NEQ
    double precision, intent(in) :: T
    double precision, intent(in) , target :: Y(NEQ)
    double precision, intent(out), target :: YP(NEQ)

    INTEGER :: I1, J, qidx(3), gidx(6)

    real(8) :: GU(6), GUG(6), QP_(3), G(6)


    if (y(3*Natom+1) == 0.0) then
        call rhss_zero_time(NEQ, y, yp)
        return
    end if

    call gaussian_average_avx(y, U, UPV, UPM)

    TRUXXG = 0d0
!$omp parallel do schedule(static) reduction(+:TRUXXG) private(QP_,qidx,gidx,G,GU,GUG)
    do i1=1,Natom
        qidx=   (/ (I1 + J*Natom, J=0,2) /)
        gidx=   (/ (I1 + J*Natom, J=3,8) /)

        G = y(gidx)

        TRUXXG = TRUXXG +     G(1)*UPM(1,I1) + 2d0*G(2)*UPM(2,I1) &
	                    + 2d0*G(3)*UPM(3,I1) +     G(4)*UPM(4,I1) &
                        + 2d0*G(5)*UPM(5,I1) +     G(6)*UPM(6,I1)

        QP_(1) = -G(1)*UPV(1,I1) - G(2)*UPV(2,I1) - G(3)*UPV(3,I1)
        QP_(2) = -G(2)*UPV(1,I1) - G(4)*UPV(2,I1) - G(5)*UPV(3,I1)
        QP_(3) = -G(3)*UPV(1,I1) - G(5)*UPV(2,I1) - G(6)*UPV(3,I1)

        GU = mm_ss(G, UPM(:,I1))
        GUG = -mm_ss(GU, G)

        GUG(1) = GUG(1) + invmass
        GUG(4) = GUG(4) + invmass
        GUG(6) = GUG(6) + invmass

        yp(qidx) = QP_
        yp(gidx) = GUG
    end do
!$omp end parallel do

    yp(NEQ) = -(0.25d0*TRUXXG + U)/real(Natom)
END SUBROUTINE RHSS0

subroutine gaussian_average(y, U, UPV, UPM)
    double precision, intent(in) :: y(:)
    double precision, intent(out) :: U, UPV(:,:), UPM(:,:)
    INTEGER :: J,I1, NN1
    real(RP), dimension(nnbmax) :: x12, y12, z12
    real(RP), dimension(nnbmax, 3) :: UX0
    real(RP), dimension(nnbmax, 6) :: UXX0


    UPM = 0.0d0
    UPV = 0.0d0

    U = 0d0
!!    !$omp parallel private(x12, y12, z12, DETA, invDETAG, qZq, expav, v0, Zq, GC, A, AG, Z, UX0, UXX0, UPV1, UPM1, J, IG, NN1)
!!    !$omp do schedule(dynamic) reduction(+:U,UPV, UPM)
    do I1=1,Natom-1
        NN1 = NNB(I1)
        if (NN1 == 0) cycle

        x12(1:NN1) = min_image(y(          I1) - y(          nbidx(1:NN1,I1)), bl)
        y12(1:NN1) = min_image(y(  Natom + I1) - y(  Natom + nbidx(1:NN1,I1)), bl)
        z12(1:NN1) = min_image(y(2*Natom + I1) - y(2*Natom + nbidx(1:NN1,I1)), bl)

        do J=1,6
            GC(1:NN1,J) = y((J+2)*Natom + nbidx(1:NN1,I1)) + y((J+2)*Natom + I1)
        end do

        call stream_kernel(NN1, x12, y12, z12, GC, U, UX0, UXX0)

        do J=1,NN1
            UPV(:,I1) = UPV(:, I1) + UX0 (J,:)
            UPV(:,nbidx(J,I1)) = UPV(:,nbidx(J,I1)) - UX0(J,:)

            UPM(:,I1) = UPM(:,I1) + UXX0(J,:)
            UPM(:,nbidx(J,I1)) = UPM(:,nbidx(J,I1)) + UXX0(J,:)
        end do
    ENDDO ! I
!!    !$omp end do
end subroutine



subroutine stream_kernel(NN1, x12, y12, z12, GC, U, UX0, UXX0)
    implicit none
    integer, intent(in) :: NN1
    real(RP), intent(in) :: x12(:), y12(:), z12(:), GC(:,:)
    double precision, intent(inout) :: U
    real(RP), intent(out) :: UX0(:,:), UXX0(:,:)


    integer :: IG, J

    call pdetminvm_sg(NN1, GC, DETA, A)

    DO IG=1,NGAUSS      ! BEGIN SUMMATION OVER GAUSSIANS
        AG(1:NN1,1) = LJA(IG) + A(1:NN1,1)
        AG(1:NN1,2) =           A(1:NN1,2)
        AG(1:NN1,3) =           A(1:NN1,3)
        AG(1:NN1,4) = LJA(IG) + A(1:NN1,4)
        AG(1:NN1,5) =           A(1:NN1,5)
        AG(1:NN1,6) = LJA(IG) + A(1:NN1,6)

        call pdetminvm_sg(NN1, AG, invDETAG, Z)
        do J=1,6
            Z(1:NN1,J) = - (LJA(IG)**2) * Z(1:NN1,J)
        end do

        Z(1:NN1,1)=LJA(IG) + Z(1:NN1,1)
        Z(1:NN1,4)=LJA(IG) + Z(1:NN1,4)
        Z(1:NN1,6)=LJA(IG) + Z(1:NN1,6)

        !Zq = matmul(Z, Q12) ! R = -2.0*Zq
        ZQ(1:NN1,1) = Z(1:NN1,1)*x12(1:NN1) + Z(1:NN1,2) * y12(1:NN1) + Z(1:NN1,3)*z12(1:NN1)
        ZQ(1:NN1,2) = Z(1:NN1,2)*x12(1:NN1) + Z(1:NN1,4) * y12(1:NN1) + Z(1:NN1,5)*z12(1:NN1)
        ZQ(1:NN1,3) = Z(1:NN1,3)*x12(1:NN1) + Z(1:NN1,5) * y12(1:NN1) + Z(1:NN1,6)*z12(1:NN1)

        !qZq = dot_product(Q12, Zq) 
        qZq(1:NN1) = Zq(1:NN1,1)*x12(1:NN1) + Zq(1:NN1,2)*y12(1:NN1) + Zq(1:NN1,3)*z12(1:NN1)

        EXPAV(1:NN1)=EXP(-qZq(1:NN1)) * SQRT(DETA(1:NN1)*invDETAG(1:NN1))

        v0(1:NN1) = 2.0_RP*expav(1:NN1)*LJC(IG)
        U = U + 0.5_RP *sum(v0(1:NN1))

        if (IG == 1) then
            DO J=1,3
                UX0(1:NN1,J) = - v0(1:NN1)*Zq(1:NN1,J)
            END DO

            UXX0(1:NN1,1) = v0(1:NN1)*(2.0_RP*Zq(1:NN1,1)*Zq(1:NN1,1) - Z(1:NN1,1))
            UXX0(1:NN1,2) = v0(1:NN1)*(2.0_RP*Zq(1:NN1,2)*Zq(1:NN1,1) - Z(1:NN1,2))
            UXX0(1:NN1,3) = v0(1:NN1)*(2.0_RP*Zq(1:NN1,3)*Zq(1:NN1,1) - Z(1:NN1,3))
            UXX0(1:NN1,4) = v0(1:NN1)*(2.0_RP*Zq(1:NN1,2)*Zq(1:NN1,2) - Z(1:NN1,4))
            UXX0(1:NN1,5) = v0(1:NN1)*(2.0_RP*Zq(1:NN1,3)*Zq(1:NN1,2) - Z(1:NN1,5))
            UXX0(1:NN1,6) = v0(1:NN1)*(2.0_RP*Zq(1:NN1,3)*Zq(1:NN1,3) - Z(1:NN1,6))
            cycle
        end if
        DO J=1,3
            UX0(1:NN1,J) = UX0(1:NN1,J) - v0(1:NN1)*Zq(1:NN1,J)
        END DO

        UXX0(1:NN1,1) = UXX0(1:NN1,1) + v0(1:NN1)*(2.0_RP*Zq(1:NN1,1)*Zq(1:NN1,1) - Z(1:NN1,1))
        UXX0(1:NN1,2) = UXX0(1:NN1,2) + v0(1:NN1)*(2.0_RP*Zq(1:NN1,2)*Zq(1:NN1,1) - Z(1:NN1,2))
        UXX0(1:NN1,3) = UXX0(1:NN1,3) + v0(1:NN1)*(2.0_RP*Zq(1:NN1,3)*Zq(1:NN1,1) - Z(1:NN1,3))
        UXX0(1:NN1,4) = UXX0(1:NN1,4) + v0(1:NN1)*(2.0_RP*Zq(1:NN1,2)*Zq(1:NN1,2) - Z(1:NN1,4))
        UXX0(1:NN1,5) = UXX0(1:NN1,5) + v0(1:NN1)*(2.0_RP*Zq(1:NN1,3)*Zq(1:NN1,2) - Z(1:NN1,5))
        UXX0(1:NN1,6) = UXX0(1:NN1,6) + v0(1:NN1)*(2.0_RP*Zq(1:NN1,3)*Zq(1:NN1,3) - Z(1:NN1,6))
    ENDDO
end subroutine stream_kernel

!    11 22 33
!    21 42 53
!    31 52 63
!    22 44 55
!    32 54 65
!    33 55 66
pure function mm_ss(A, B) result(C)
    double precision, intent(in) :: A(6), B(6)
    double precision :: C(6)

    double precision :: x22, x33, x53, x55

    x22 = A(2) * B(2)
    x33 = A(3) * B(3)
    x53 = A(5) * B(3)
    x55 = A(5) * B(5)
    C(1) = A(1)*B(1) + x22 + x33
    C(2) = A(2)*B(1) + A(4)*B(2) + x53
    C(3) = A(3)*B(1) + A(5)*B(2) + A(6)*B(3)
    C(4) = x22 + A(4)*B(4) + A(5)*B(5)
    C(5) = A(3)*B(2) + A(5)*B(4) + A(6)*B(5)
    C(6) = x33 + x55 + A(6)*B(6)
end function mm_ss


subroutine rhss_zero_time(NEQ, y, yp)
    use utils, only: RP
    implicit none
    integer, intent(in) :: NEQ
    double precision, intent(in) :: y(:)
    double precision, intent(out) :: yp(:)

    double precision :: qij(3), qi(3), qj(3)
    double precision :: rsq(nnbmax)
    integer :: i, j

    yp = 0.0d0

    yp(3*Natom + 1 : 4*Natom) = invmass
    yp(6*Natom + 1 : 7*Natom) = invmass
    yp(8*Natom + 1 : 9*Natom) = invmass

    U=0.0d0

    DO I=1,Natom-1
        if (nnb(i) == 0) cycle
        rsq(1:NNB(i)) = min_image(y(nbidx(1:NNB(i),i)) - y(          i), bl)**2 &
            + min_image(y(  Natom + nbidx(1:NNB(i),i)) - y(  Natom + i), bl)**2 &
            + min_image(y(2*Natom + nbidx(1:NNB(i),i)) - y(2*Natom + i), bl)**2
        DO J=1,NGAUSS
            U = U + LJC(J)*sum(EXP(-LJA(J)*rsq(1:NNB(i))))
        END DO
    ENDDO

    yp(NEQ) = -U/real(Natom )
end subroutine

