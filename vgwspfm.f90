module vgwspfm
    use utils
    use iso_c_binding
    implicit none
    private
    public :: vgwspfminit, vgw0spfm,vgwspfmcleanup
    
    integer :: Natom, Nmax, nnz, n3rows16
    real*8 :: BL, rfullmatsq = 5.5d0**2
    real*8, dimension(10) :: LJA, LJC
    integer :: NGAUSS
    integer, allocatable :: NBIDX(:,:), NNB(:),  &
            Gbia(:), Gbja(:), Giia(:), Gbiia(:)
    
    integer(C_INT), allocatable, target :: Gia(:), Gja(:)
    real(C_DOUBLE), allocatable, target :: G(:)

    real*8 :: gama, gamap, U
    real*8, allocatable :: Gb(:,:,:), Gbdiag(:,:,:), &
            UXYdiag(:,:,:), UX(:,:)

    real*8, allocatable, target :: GU(:,:), UXYf(:,:)

    real*8 :: invmass, RC, mass, dt0, dtmax, dtmin, vgw_atol(3)
    logical :: finished
    integer :: nnbmax

    interface
        real(C_DOUBLE) function cholmod_logdet(G, ia, ja, N) BIND(C)
            use, intrinsic :: iso_c_binding
            type(C_PTR), value :: G, ia, ja
            integer(C_INT), value :: N
        end function
    end interface


contains

subroutine vgwspfminit(Np, species, M, rcutoff, massx)
    use omp_lib
    implicit none
    integer, intent(in) :: Np
    character(*), intent(in) :: species
    real*8, intent(in), optional :: M, rcutoff, massx

    Natom = Np
    n3rows16 = ((3*Natom+1)/2)*2 ! align first dimension on 16byte boundary

    allocate(NNB(Natom), NBIDX(Natom,Natom), Gbia(Natom+1), Gbja(1), &
        Gia(3*Natom+1), Gja(1), Giia(3*Natom), Gbiia(Natom), &
        Gb(3,3,1), Gbdiag(3,3,Natom), &
        UXYdiag(3,3,Natom), &
        UX(3,Natom), UXYf(n3rows16,3*Natom), &
        GU(n3rows16,3*Natom))
    
    
    include 'species.f90'
end subroutine


subroutine vgwspfmcleanup()
    deallocate(NNB, NBIDX, Gbia, Gbja, &
        Gia, Gja, Giia, Gbiia, &
        Gb, Gbdiag, &
        UXYdiag, &
        UX, UXYf, GU)
end subroutine


SUBROUTINE vgw0spfm(Q0, BL_, beta, Ueff, rt)
    IMPLICIT NONE
    double precision, intent(in) :: Q0(:,:), beta, BL_
    double precision, intent(out) :: Ueff
    double precision, intent(out), optional :: rt
    real*8 :: LOGDET, logrho, TSTOP, start_time, stop_time, T
    integer :: i, ncalls

    double precision, allocatable, target :: Y(:)
    double precision, allocatable :: RWORK(:), YP(:), ATOL(:)
    integer, allocatable :: IWORK(:)

    integer :: NEQ, ITOL, ITASK, IOPT, MF, ISTATE, LRW, LIW
    double precision ::  RTOL

    Natom = size(Q0, 2)
    BL = BL_

    call interaction_lists(Q0)

    if (size(Gja) < nnz) then

        deallocate(Gb, Gbja, Gja)

        allocate(Gb(3,3,nnz/9), Gbja(nnz/9), Gja(nnz))
    end if

    call init_sparse_pattern(Q0)
    call test_ia_ja(Gia, Gja)

    NEQ = 3*Natom + nnz + 1


    LRW = 20 + 16*NEQ
    LIW = 30
    allocate(Y(NEQ), YP(NEQ), ATOL(NEQ), RWORK(LRW), IWORK(LIW))

    ITOL=2
    RTOL=1d-5
    ATOL(1:3*Natom) = vgw_atol(1)
    ATOL(3*Natom+1:3*Natom+nnz)=vgw_atol(2)
    ATOL(NEQ) = vgw_atol(3)
    ITASK=1
    ISTATE=1
    IOPT = 1
    MF=10
    IWORK=0

    IWORK(6) = 50000 ! MXSTEP

    RWORK(5)=dt0
    RWORK(6)=dtmax
    RWORK(7)=dtmin

    y=0d0
    y(1:3*Natom) = reshape(Q0, (/ 3*Natom /) )

    call cpu_time(start_time)

    T=0
    TSTOP = 0.5d0*beta
    CALL DLSODE(RHSSspFM,NEQ,Y,T,TSTOP,ITOL,RTOL,ATOL,ITASK,ISTATE,IOPT,&
        RWORK,LRW,IWORK,LIW,JAC,MF)

    gama = Y(NEQ) * real(Natom)

    logdet = cholmod_logdet(C_LOC(y(3*Natom+1)), C_LOC(Gia), C_LOC(Gja), &
        3*Natom)
    if (logdet == 1234.1234d0) then
        print *, 'The G matrix is not positive definite!'
        open(10,file='nonposdef.xyz', status='REPLACE')
        write(10, '(I5/2F16.8/("Ne",3F14.8))') Natom, bl, beta, q0
        close(10)
        stop
    end if
        

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


SUBROUTINE vgw0spfmgs(Q0, BL_, TAUMAX, Havg, rt, logfd, cfgname)
    IMPLICIT NONE
    double precision, intent(in) :: Q0(:,:), TAUMAX, BL_
    double precision, intent(out) :: Havg
    double precision, intent(out), optional :: rt
    integer, intent(in), optional :: logfd
    character(*), intent(in), optional :: cfgname
    real*8 :: LOGDET, logrho(2), dbeta, TSTOP, start_time, stop_time, T
    integer :: i, ncalls, nsteps

    double precision, allocatable, target :: Y(:)
    double precision, allocatable :: RWORK(:), YP(:), ATOL(:)
    integer, allocatable :: IWORK(:)

    integer :: NEQ, ITOL, ITASK, IOPT, MF, ISTATE, LRW, LIW
    double precision ::  RTOL

    Natom = size(Q0, 2)
    BL = BL_

    call interaction_lists(Q0)

    if (size(Gja) < nnz) then

        deallocate(Gb, Gbja, Gja)

        allocate(Gb(3,3,nnz/9), Gbja(nnz/9), Gja(nnz))
    end if

    call init_sparse_pattern(Q0)
    call test_ia_ja(Gia, Gja)

    NEQ = 3*Natom + nnz + 1


    LRW = 20 + 16*NEQ
    LIW = 30
    allocate(Y(NEQ), YP(NEQ), ATOL(NEQ), RWORK(LRW), IWORK(LIW))

    ITOL=2
    RTOL=1d-5
    ATOL(1:3*Natom) = vgw_atol(1)
    ATOL(3*Natom+1:3*Natom+nnz)=vgw_atol(2)
    ATOL(NEQ) = vgw_atol(3)
    ITASK=1
    ISTATE=1
    IOPT = 1
    MF=10
    IWORK=0

    IWORK(6) = 50000 ! MXSTEP

    RWORK(5)=dt0
    RWORK(6)=dtmax
    RWORK(7)=dtmin

    y=0d0
    y(1:3*Natom) = reshape(Q0, (/ 3*Natom /) )
    dbeta = min(0.1*0.5*TAUMAX, 1d0)
    nsteps = ceiling(0.5*TAUMAX/dbeta) 
    dbeta = 0.5*TAUMAX/real(nsteps)

    call cpu_time(start_time)

    T=0
    logrho=0
    do i=1,nsteps
        logrho(1) = logrho(2)

        TSTOP = T  + dbeta
        CALL DLSODE(RHSSspFM,NEQ,Y,T,TSTOP,ITOL,RTOL,ATOL,ITASK,ISTATE,IOPT,&
            RWORK,LRW,IWORK,LIW,JAC,MF)

        gama = Y(NEQ) * real(Natom)

        logdet = cholmod_logdet(C_LOC(y(3*Natom+1)), C_LOC(Gia), &
            C_LOC(Gja), 3*Natom)

        logrho(2) = 2.0*gama - 0.5*logdet - 1.5*Natom*log(4.0*M_PI)

        Havg = -(logrho(2) - logrho(1)) *0.5d0/ dbeta
        if (present(logfd)) then
            write (logfd,*) 2.0*T, Havg
        end if
        if (present(cfgname)) then
            open(57,file=trim(cfgname),status='REPLACE')
            write (57,*) Natom
            write(57, *) 'vgw intermediate'
            write(57, '(("Ne ", 3F13.8))') y(1:3*Natom)
            close(57)
        end if
    end do

    call cpu_time(stop_time)
    ncalls = IWORK(12)
    write (*,*) IWORK(11), 'steps,', IWORK(12), ' RHSS calls, logdet =', logdet

    deallocate(y, yp, RWORK, IWORK, ATOL)

    if (present(rt)) then
        rt = (stop_time - start_time) / real(ncalls)
     end if
END SUBROUTINE


subroutine interaction_lists(Q)
    implicit none
    real*8, intent(in) :: Q(:,:)
    integer :: N,I,J, NN
    real*8 rsq,rc2,qij(3)

    N = size(Q, 2)
    rc2=rc**2


    NNB = 0
    nnz = 0
!$OMP DO SCHEDULE(DYNAMIC) PRIVATE(I, J, QIJ, RSQ)
    do I=1,N-1
        NN = 0
        do J=I+1,N
            qij=Q(:,I)-Q(:,J)
            rsq = sum(min_image(qij, BL)**2)
            if(rsq <= rc2) then
                NN = NN + 1
                NBIDX(NN, I) = J
            endif
            if (rsq <= rfullmatsq) then
                nnz = nnz + 1
            end if
        enddo
        NNB(i) = NN
    enddo
!$OMP END DO
    nnbmax = maxval(nnb)
    nnz = 9*(2*nnz + Natom)
end subroutine interaction_lists


subroutine init_sparse_pattern(Q0)
    real*8, intent(in) :: Q0(:,:)
    integer :: I, J, p, p0, ib, jb, nz
    double precision :: qij(3), rsq

    p = 1
    do ib=1,Natom
        i = 3*(ib - 1) + 1

        p0 = p

        do jb=1,Natom
            j = 3*(jb - 1) + 1

            qij=Q0(:,ib)-Q0(:,jb)
            rsq = sum(min_image(qij, BL)**2)

            if (rsq <= rfullmatsq) then
                Gja(p : p+2) = (/ j, j+1, j+2/)

                if (ib == jb) then
                    Giia(i) = p
                    Gbiia(ib) = (p0-1)/9 + (p-p0) / 3 + 1
                end if

                p = p + 3
            end if
        enddo
        nz = p - p0

        Gia(i : i+2) = (/ p0, p, p + nz/)

        Gja(p : p + nz - 1) = Gja(p0 : p - 1) 
        Gja(p + nz : p + 2*nz - 1) = Gja(p0 : p - 1) 

        p = p + 2*nz

        Giia(i+1) = Giia(i)   + nz + 1
        Giia(i+2) = Giia(i+1) + nz + 1
    enddo
    Gia(3*Natom + 1) = p
    if (nnz /= p - 1) then
        write (*,*) 'nnz error!', nnz, p-1
        stop
    end if
end subroutine

subroutine JAC()
end subroutine

subroutine test_ia_ja(ia, ja)
    integer, intent(in) :: ia(:), ja(:)
    integer,allocatable :: P(:,:)
    integer :: i, ptr, N
    
    allocate(P(size(ia)-1, size(ia)-1))

    N = size(ia) - 1

    P = 0
    do i=1,N
        do ptr=ia(i), ia(i+1)-1
            P(i, ja(ptr)) = 1
        end do
    end do

    if (sum(P-transpose(P)) /= 0) then
        write (*,*) 'The pattern is not symmetric!'
        stop
    end if

    deallocate(P)
end  subroutine


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

!    write (*,*) T

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
    double precision, intent(in) :: y(NEQ)
    double precision, intent(out) :: yp(NEQ)

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
end module vgwspfm
