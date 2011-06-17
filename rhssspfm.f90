SUBROUTINE RHSSspFM(NEQ, T, Y, YP)!, RPAR, IPAR)
!    use omp_lib
    IMPLICIT NONE
    integer, intent(in) :: NEQ!, IPAR(:)
    double precision, intent(in) :: T, Y(NEQ)!, RPAR(:)
    double precision, intent(out) :: YP(NEQ)
    INTEGER :: J,I1,I2,IG, J2, Gb_ptr, info, nGbcols, ptrb(Natom+1), ptre(Natom+1), k, I2_3
    double precision :: AG(3,3), DETA,DETAG,QZQ,U12, G12(3,3),A(3,3), &
            Zq(3), Z(3,3),Q12(3), UXY0(3,3), UX0(3), TrUXYG
    double precision, pointer :: UG_slab(:,:)

    integer :: job(6) = (/ 1, 1, 1, 0, 0, 1 /)
    character(6) :: matdescra='GxxFxx'

    ! first call, G=0
    if (y(3*Natom+1)==0d0) then
        call rhss_zero_time(NEQ, y, yp)
        return
    end if

    call unpack_y(y, Q, Gb(:,:,1:nnzb), gama)
 

    U = 0; UX = 0; UXY = 0;  UXYdiag=0; UXYf=0

    Gbdiag = Gb(:,:, FMdiag)

    do I1=1,Natom-1
        Gb_ptr =  FMdiag(I1) + 1
        ! what if there is no other interaction on the right hand side?
        if (fmdiag(I1) + 1 == Gbia(I1+1)) then
            Gb_ptr = fmdiag(i1)
        end if

        DO J2=1,NNB(I1)
            I2 = NBIDX(J2, I1)
            
            Q12 = Q(:, I1) - Q(:, I2)
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
!                UXY(:,:,Gb_ptr) = -UXY0
                Gb_ptr = Gb_ptr + 1
            end if
            
        end do ! I2
    end do ! I1

    do I1=1,Natom
        UXYf(3*I1-2 : 3*I1, 3*I1-2 : 3*I1) = UXYdiag(:,:,I1)
    end do


!    UXY(:,:,FMDIAG) = UXYdiag
!    call bsr_copy_up_lo(Gbia, Gbja, UXY)




    !call mkl_dcsrbsr(job, Natom, 3, 9, Gcsr, Grja, Gria, Gb, Gbja, Gbia, info)
    !call mkl_dcsrbsr(job, Natom, 3, 9, UXYr, Grja, Gria, UXY, Gbja, Gbia, info)

    ! dG/dt = - G U G = -GU G = (G^T GU^T)^T = (G^T UG)^T
    !call mkl_dcsrmultd('N', 3*Natom, 3*Natom, 3*Natom, UXYr, Grja, Gria, Gcsr, &
    !    Grja, Gria, UG, 3*Natom)



    ptrb(1:Natom) = Gbia(1:Natom)
    ptre(1:Natom) = Gbia(2:Natom+1)
    call mkl_dbsrmm('N', Natom, 3*Natom, Natom, 3, 1d0, matdescra, Gb, Gbja, &
        ptrb, ptre, UXYf, 3*Natom, 0d0, GU, 3*Natom)

    UG = transpose(GU)

    do I1=1,Natom
        ptrb(1) = 1
        ptre(1) = Gbia(I1+1) - Gbia(I1) + 1
        do J2=fmdiag(I1),Gbia(I1+1)-1
            I2_3=3*Gbja(J2) - 2
            UG_slab => UG(:, I2_3 : I2_3 + 2)
            call mkl_dbsrmm('N', 1, 3, Natom, 3, -1d0, matdescra, &
                Gb(:,:,Gbia(I1):Gbia(I1+1)-1), Gbja(Gbia(I1):Gbia(I1+1)-1), &
                ptrb, ptre, UG_slab, 3*Natom, 0d0, &
                GPb(:,:,J2), 3)
        end do
        
        do k=1,3
            GPb(k,k,FMDIAG(I1)) = GPb(k,k,FMDIAG(I1)) + invmass
        end do
        
    end do
    call bsr_copy_up_lo(Gbia, Gbja, GPb)

    TrUXYG = 0d0
    do I1=1,3*Natom
        TrUXYG = TrUXYG + UG(I1,I1)
    end do

!    call bsrdense(Gb, Gbia, Gbja, GU)
!    call dgemm('N', 'N', 3*Natom, 3*Natom, 3*Natom, 1d0, UXYf, 3*Natom, &
!        GU, 3*Natom, 0d0, UG, 3*Natom)
!    call dgemm('N', 'N', 3*Natom, 3*Natom, 3*Natom, -1d0, GU, 3*Natom, &
!        UG, 3*Natom, 0d0, UXYf, 3*Natom)
!    do I1=1,3*Natom
!        UXYf(I1,I1) = UXYf(I1,I1) + invmass
!    end do
!
!    do I1=1,Natom-1
!        DO I2=I1+1,Natom
!            Q12 = Q(:,I1) - Q(:,I2)
!            Q12 = min_image(Q12, bl)
!            if (sum(q12**2) > rfullmatsq) then
!                UXYf(3*I1-2 : 3*I1, 3*I2-2 : 3*I2) = 0
!                UXYf(3*I2-2 : 3*I2, 3*I1-2 : 3*I1) = 0
!           end if
!        end do
!    end do
!
!    call bsrdense(GPb, Gbia, Gbja, GU)
!
!    if (sum(abs(GU-UXYf)) > 1d-5) then
!        write (*,*) 'bs'
!    else if (abs( sum(UXY(:,:,1:nnzb) * Gb(:,:,1:nnzb)) -TrUXYG ) > 1d-5) then
!        write (*,*) 'bs squared'
!    end if

    call mkl_dbsrgemv('N', Natom, 3, Gb, Gbia, Gbja, UX, QP)
    QP = -1d0*QP


    !gamap = -(0.25d0 * sum(UXY(:,:,1:nnzb) * Gb(:,:,1:nnzb)) + U)/real(Natom)
    gamap = -(0.25d0 * TrUXYG + U)/real(Natom)

    call pack_y(QP, GPb(:,:,1:nnzb), gamap, yp)
end subroutine RHSSspFM

subroutine rhss_zero_time(NEQ, y, yp)
    integer, intent(in) :: NEQ
    double precision, intent(in) :: y(:)
    double precision, intent(out) :: yp(:)

    double precision :: qij(3), rsq
    integer :: i, j, k

    Q = reshape(y(1:3*Natom), (/ 3, Natom /) )

    QP = 0d0

    GPb = 0
    do i=1,Natom
        do k=1,3
            GPb(k,k,fmdiag(i)) = invmass
        end do
    end do

    U=0d0

    DO I=1,Natom-1
        DO J=1,NNB(I)
                qij = Q(:,I) - Q(:,NBIDX(J, I))
                rsq = sum(min_image(qij, BL)**2)
                U = U + sum(LJC(1:NGAUSS)*EXP(-LJA(1:NGAUSS)*rsq))
        ENDDO
    ENDDO

    gamap = -U/real(Natom)

    call pack_y(QP, GPb(:,:,1:nnzb), gamap, yp)
end subroutine

subroutine bsr_copy_up_lo(p, i, x)
    integer, intent(in) :: p(:), i(:)
    double precision, intent(inout) :: x(:,:,:)

    integer :: j, k, r, nstored(size(p) - 1)

    nstored = 0
    do j=1,size(p) - 1
        ! skip the diagonal
        nstored(j) = nstored(j) + 1

        do k=p(j)+nstored(j), p(j+1)-1
            r = i(k)
            nstored(r) = nstored(r)+1
            x(:,:,p(r) + nstored(r) - 1) = transpose(x(:,:,k))
        end do
   end do
end subroutine


subroutine bsr_copy_lo_up(p, i, x, d)
    integer, intent(in) :: p(:), i(:), d(:)
    double precision, intent(inout) :: x(:,:,:)

    integer :: j, k, c, nstored(size(p) - 1), N

    N = size(p) - 1

    nstored = d - p(1:N) + 1
    do j=1,size(p) - 1
        do k=p(j), d(j) - 1
            c = i(k)
            nstored(c) = nstored(c)+1
            x(:,:,p(c) + nstored(c) - 1) = transpose(x(:,:,k))
        end do
   end do
end subroutine

subroutine bsr_chop(ia, slices, ptrb, ptre)
    integer, intent(in) :: ia(:), slices(:)
    integer, intent(out) :: ptrb(:), ptre(:)

    integer :: k

    do k=1,size(slices)
        ptrb(k) = ia(slices(k))
        ptre(k) = ia(slices(k)+1)
    end do
end subroutine

subroutine bsrdense(x, ia, ja, A)
    double precision, intent(in) :: x(:,:,:)
    integer, intent(in) :: ia(:), ja(:)
    double precision, intent(out) :: A(:,:)
    
    integer :: i, j, jit, N, mblk

    A = 0
    N = size(ia) - 1
    mblk = size(x, 1)
    do i=1,N
        do jit=ia(i),ia(i+1)-1
            j = ja(jit)
            
            A(mblk*(i-1) + 1 : mblk*i, mblk*(j-1) + 1 : mblk*j) = x(:,:,jit)
        end do
    end do
end subroutine
