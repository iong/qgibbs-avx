SUBROUTINE RHSSspFM(NEQ, T, Y, YP)!, RPAR, IPAR)
    IMPLICIT NONE
    integer, intent(in) :: NEQ!, IPAR(:)
    double precision, intent(in) :: T!, RPAR(:)
    double precision, intent(in), target :: Y(NEQ)
    double precision, intent(out), target :: YP(NEQ)
    INTEGER :: J,I1,I2,IG, J2, Gb_ptr, ptrb(3*Natom+1), ptre(3*Natom+1), p
    double precision :: AG(3,3), DETA,DETAG,QZQ,U12, G12(3,3),A(3,3), &
            Zq(3), Z(3,3),Q12(3), UXY0(3,3), UX0(3), TrUXYG

    integer :: job(10), info
    character(6) :: matdescra='GxxFxx'

    double precision, pointer :: Gptr(:), GPptr(:)

    write (*,*) T

    ! first call, G=0
    if (y(3*Natom+1)==0d0) then
        call rhss_zero_time(NEQ, y, yp)
        return
    end if

    Gptr => y(3*Natom+1 : 3*Natom + nnz)
    GPptr => yp(3*Natom+1 : 3*Natom + nnz)

    job(1:6) = (/ 0, 1, 1, 0, 0, 1 /)
    call mkl_dcsrbsr(job, 3*Natom, 3, 9, Gptr, Gja, Gia, Gb, Gbja, Gbia, info)
    Gbdiag = Gb(:,:,Gbiia)

    U = 0; UX = 0; UXYdiag=0; UXYf=0


    do I1=1,Natom-1
        Gb_ptr =  Gbiia(I1) + 1
        ! what if there is no other interaction on the right hand side?
        if (Gbiia(I1) + 1 == Gbia(I1+1)) then
            Gb_ptr = Gbiia(i1)
        end if

        DO J2=1,NNB(I1)
            I2 = NBIDX(J2, I1)
            
            Q12 = y(3*I1-2 : 3*I1) - y(3*I2-2 : 3*I2)
            Q12 = min_image(Q12, bl)

            G12=Gbdiag(:,:, I1) + Gbdiag(:,:, I2)
            if ( Gbja(Gb_ptr) == I2) then
                G12 = G12 - Gb(:,:,Gb_ptr)- transpose(Gb(:,:,Gb_ptr))
            end if

            call detminvm(G12, DETA, A)
            DETA = 1.0d0/DETA

            UX0 = 0d0; UXY0 = 0d0
            DO IG=1,NGAUSS ! BEGIN SUMMATION OVER GAUSSIANS
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

                if (DETA*DETAG <= 0.0d0 ) then
                    write (69, '(F16.8)'), Y
                    write (*,*) DETA, DETAG
                    stop
                end if
                U12 = SQRT(DETA/DETAG)*EXP(-qZq)*LJC(IG)
                U = U + U12

                UX0 = UX0 - 2d0*U12*Zq
                do J=1,3
                    UXY0(:,J) = UXY0(:,J) + 2d0*U12*(2d0*Zq*Zq(J) - Z(:,J))
                end do
            end do ! IG

! Avoid load and store as much as possbile. Store now, process as a stream
! later. Much faster.
            UX(:, I1) = UX(:, I1) + UX0
            UX(:, I2) = UX(:, I2) - UX0

            UXYdiag(:,:, I1) = UXYdiag(:,:, I1) + UXY0
            UXYdiag(:,:, I2) = UXYdiag(:,:, I2) + UXY0

            UXYf(3*I1-2 : 3*I1, 3*I2-2 : 3*I2) = -UXY0
            UXYf(3*I2-2 : 3*I2, 3*I1-2 : 3*I1) = -transpose(UXY0)

            if (Gbja(Gb_ptr) == I2) then
                Gb_ptr = Gb_ptr + 1
            end if

        end do ! I2
    end do ! I1

    do I1=1,Natom
        UXYf(3*I1-2 : 3*I1, 3*I1-2 : 3*I1) = UXYdiag(:,:,I1)
    end do

    call mkl_dcsrgemv('N', 3*Natom, Gptr, Gia, Gja, UX, yp)
    yp(1:3*Natom) = -yp(1:3*Natom)

    GU = 0
    ptrb(1:3*Natom) = Gia(1:3*Natom)
    ptre(1:3*Natom) = Gia(2:3*Natom+1)
    call mkl_dcsrmm('N', 3*Natom, 3*Natom, 3*Natom, -1d0, matdescra, &
        Gptr, Gja, ptrb, ptre, UXYf, n3rows16, 0d0, GU, n3rows16)

    TrUXYG = 0d0
    do I1=1,3*Natom
        TrUXYG = TrUXYG - GU(I1,I1)
    end do

    call mmcsrsym(3*Natom, 3*Natom, 3*Natom, GU, Gptr, Gia, Gja, GPptr)
    call csr_symmetrize(Gia, Gja, GPptr)

    do I1=1,3*Natom
        GPptr(Giia(I1)) = GPptr(Giia(I1)) + invmass
    end do


    yp(NEQ) = -(0.25d0 * TrUXYG + U)/real(Natom)
end subroutine RHSSspFM

subroutine rhss_zero_time(NEQ, y, yp)
    integer, intent(in) :: NEQ
    double precision, intent(in) :: y(:)
    double precision, intent(out) :: yp(:)

    double precision :: qij(3), rsq
    integer :: i, j

    yp = 0

    do i=1,3*Natom
        yp(3*Natom + Giia(i)) = invmass
    end do

    U=0d0

    DO I=1,Natom-1
        DO J=1,NNB(I)
                qij = y(3*I - 2 : 3*I) - y(3*NBIDX(J, I) - 2 : 3*NBIDX(J, I))
                rsq = sum(min_image(qij, BL)**2)
                U = U + sum(LJC(1:NGAUSS)*EXP(-LJA(1:NGAUSS)*rsq))
        ENDDO
    ENDDO

    yp(NEQ) = -U/real(Natom)
end subroutine

subroutine mmcsrsym(nrowsA, ncolsC, nrowsB, A, B, ia, ja, C)
    implicit none
    integer, intent(in) :: nrowsA, ncolsC, nrowsB, ia(:), ja(:)
    double precision, intent(in) :: A(:,:), B(:)
    double precision, intent(out) :: C(:)
    integer :: j, k, p, k_j_p

    double precision :: x(size(A, 1))
    C = 0d0
    do k=1,nrowsB
        do k_j_p=ia(k), ia(k+1)-1
            j = ja(k_j_p)
            !do i=1,nrowsA
                !c(i,j) = c(i, j) + a(i, k) * b(k, j)
            !end do
            !x = a(:,k)*b(k_j_p)


            do p=ia(j),ia(j+1)-1
                c(p) = c(p) + a(ja(p), k) * b(k_j_p)
            end do
        end do
    end do
end subroutine

subroutine csr_dense(ia, ja, x, A)
    integer, intent(in) :: ia(:), ja(:)
    double precision, intent(in) :: x(:)
    double precision, intent(out) :: A(:,:)

    integer :: i, p, N

    N = size(ia) - 1

    A = 0d0
    do i=1,N
        do p=ia(i),ia(i+1)-1
            A(i, ja(p)) = x(p)
        end do
    end do
end subroutine

subroutine csr_symmetrize(ia, ja, x)
    integer, intent(in) :: ia(:), ja(:)
    double precision, intent(inout) :: x(:)

    double precision :: xavg

    integer :: w(size(ia) - 1), i, j, p, N

    N = size(ia) - 1

    w = 0

    do i=1,N
        do p=ia(i), ia(i+1)-1
            j = ja(p)
            if (j .le. i) cycle

            xavg = 0.5d0* (x(p) + x(ia(j) + w(j)) )
            x(p) = xavg
            x(ia(j) + w(j)) = xavg

            w(j) = w(j) + 1
        end do
    end do
end subroutine
