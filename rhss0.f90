SUBROUTINE RHSS0(NEQ, T, Y, YP)
    use utils, only: RP
    IMPLICIT NONE
    integer, intent(in) :: NEQ
    double precision, intent(in) :: T
    double precision, intent(in) , target :: Y(NEQ)
    double precision, intent(out), target :: YP(NEQ)

    INTEGER :: J,I1,I2,IG, NN1
    real(RP), dimension(nnbmax) :: x12, y12, z12, DETA, invDETAG, qZq, expav, v0
    real(RP), dimension(nnbmax, 3) :: Zq, UX0
    real(RP), dimension(nnbmax, 6) :: GC, A, AG, Z, UXX0
    real(8) :: Ulocal, GU(Natom, 6), UPV1(3), UPM1(6)
    type parray
        real(8), pointer :: p(:)
    end type
    type(parray) :: G(6), GUG(6)


    UPM = 0.0d0
    UPV = 0.0d0

    Ulocal = 0.0

    if (y(3*Natom+1) == 0.0) then
        call rhss_zero_time(NEQ, y, yp)
        return
    end if

    do J=1,6
        G(J)%p => y((3+J-1)*Natom + 1 : (3 + J) * Natom)
    end do
        
    do I1=1,Natom-1
        NN1 = NNB(I1)
        if (NN1 == 0) cycle

        x12(1:NN1) = min_image(y(          I1) - y(          nbidx(1:NN1,I1)), bl)
        y12(1:NN1) = min_image(y(  Natom + I1) - y(  Natom + nbidx(1:NN1,I1)), bl)
        z12(1:NN1) = min_image(y(2*Natom + I1) - y(2*Natom + nbidx(1:NN1,I1)), bl)

        do J=1,6
            GC(1:NN1,J) = G(J)%p(nbidx(1:NN1,I1)) + G(J)%p(I1)
        end do

        call pdetminvm_sg(GC(1:NN1,:), DETA, A)

        DO IG=1,NGAUSS      ! BEGIN SUMMATION OVER GAUSSIANS
            AG(1:NN1,1) = LJA(IG) + A(1:NN1,1)
            AG(1:NN1,2) =           A(1:NN1,2)
            AG(1:NN1,3) =           A(1:NN1,3)
            AG(1:NN1,4) = LJA(IG) + A(1:NN1,4)
            AG(1:NN1,5) =           A(1:NN1,5)
            AG(1:NN1,6) = LJA(IG) + A(1:NN1,6)

            call pdetminvm_sg(AG(1:NN1,:), invDETAG, Z)
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
            Ulocal=Ulocal + 0.5_RP *sum(v0)

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

        ! I2
        UPV1 = UPV(I1,:)
        UPM1 = UPM(I1,:)
        do J=1,NN1
            UPV1 = UPV1 + UX0 (J,:)
            UPV(nbidx(J,I1),:) = UPV(nbidx(J,I1),:) - UX0(J,:)

            UPM1 = UPM1 + UXX0(J,:)
            UPM(nbidx(J,I1),:) = UPM(nbidx(J,I1),:) + UXX0(J,:)
        end do
        UPV(I1,:) = UPV1
        UPM(I1,:) = UPM1
    ENDDO ! I

    TRUXXG =  sum(UPM(:,1)*G(1)%p) + 2*sum(UPM(:,2)*G(2)%p) &
        + 2.0*sum(UPM(:,3)*G(3)%p) +   sum(UPM(:,4)*G(4)%p) &
        + 2.0*sum(UPM(:,5)*G(5)%p) +   sum(UPM(:,6)*G(6)%p)
    U = Ulocal

    yp(        1:  Natom) = - G(1)%p*UPV(:,1) - G(2)%p*UPV(:,2) - G(3)%p*UPV(:,3)
    yp(  Natom+1:2*Natom) = - G(2)%p*UPV(:,1) - G(4)%p*UPV(:,2) - G(5)%p*UPV(:,3)
    yp(2*Natom+1:3*Natom) = - G(3)%p*UPV(:,1) - G(5)%p*UPV(:,2) - G(6)%p*UPV(:,3)

    do J=1,6
        GUG(J)%p => yp(3*Natom + (J-1)*Natom + 1 : 3*Natom + J*Natom)
    end do
    call pmm_sg(reshape(y(3*Natom+1:9*Natom), (/Natom, 6/)), UPM, GU)
    call pmm_sg(GU, y(3*Natom+1:9*Natom), yp(3*Natom+1:9*Natom))
    yp(3*Natom+1:9*Natom) = -yp(3*Natom+1:9*Natom)
    GUG(1)%p = GUG(1)%p + invmass
    GUG(4)%p = GUG(4)%p + invmass
    GUG(6)%p = GUG(6)%p + invmass
!    do i1=1,Natom
!        k1=             3*(i1 - 1) + 1
!        k2 =  3*Natom + 6*(i1 - 1) + 1
!        k3 =  9*Natom + 9*(i1 - 1) + 1
!        k4 = 18*Natom + 3*(i1 - 1) + 1
!
!        if (NEQ > 21*Natom ) then
!            yp(k3 : k3+8) = - reshape(matmul(GU, Qnk(:,:,i1)), (/ 9 /) )
!
!            yp(k4 : k4 + 2 ) = - matmul(transpose(Qnk(:,:,i1)), UPV(:,I1))
!        end if
!    end do

    yp(NEQ) = -(0.25_RP*TRUXXG + U)/real(Natom)
END SUBROUTINE RHSS0

!    11 22 33
!    21 42 53
!    31 52 63
!    22 44 55
!    32 54 65
!    33 55 66
subroutine pmm_sg(A, B, C)
    double precision, intent(in) :: A(:,:), B(size(A,1),6)
    double precision, intent(out) :: C(size(A,1),6)

    double precision, dimension(size(A, 1)) :: x22, x33, x53, x55
    
    !call atom_range(i1, i2)

    x22 = A(:,2) * B(:,2)
    x33 = A(:,3) * B(:,3)
    x53 = A(:,5) * B(:,3)
    x55 = A(:,5) * B(:,5)
    C(:,1) = A(:,1)*B(:,1) + x22 + x33
    C(:,2) = A(:,2)*B(:,1) + A(:,4)*B(:,2) + x53
    C(:,3) = A(:,3)*B(:,1) + A(:,5)*B(:,2) + A(:,6)*B(:,3)
    C(:,4) = x22 + A(:,4)*B(:,4) + A(:,5)*B(:,5)
    C(:,5) = A(:,3)*B(:,2) + A(:,5)*B(:,4) + A(:,6)*B(:,5)
    C(:,6) = x33 + x55 + A(:,6)*B(:,6)
end subroutine

function matmul_sgs(A, B) result (C)
    double precision :: A(3,3), B(3,3)
    double precision :: C(6)

    C(1:3) = matmul(A, B(:,1))
    C(4) = sum(A(2,:)*B(:,2))
    C(5) = sum(A(3,:)*B(:,2))
    C(6) = sum(A(3,:)*B(:,3))
end function

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
        rsq(1:NNB(i)) = min_image(y(nbidx(1:NNB(i),i)) - y(i), bl)**2 &
            + min_image(y(Natom + nbidx(1:NNB(i),i)) - y(Natom + i), bl)**2 &
            + min_image(y(2*Natom + nbidx(1:NNB(i),i)) - y(2*Natom + i), bl)**2
        DO J=1,NGAUSS
                U = U + LJC(J)*sum(EXP(-LJA(J)*rsq(1:NNB(i))))
        END DO
    ENDDO

    yp(NEQ) = -U/real(Natom )
end subroutine

