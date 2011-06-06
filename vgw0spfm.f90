SUBROUTINE vgw0spfm(Q0, BL_, TAUMAX, W)
    IMPLICIT NONE
    double precision, intent(in) :: Q0(:,:), TAUMAX, BL_
    double precision, intent(out) :: W
    double precision :: logdetg, GF(3*Natom, 3*Natom)
    double precision :: DT, next_stop
    integer :: j, nnzb, I1, I2, J2, info, s
    logical :: mm = .FALSE.

    BL = BL_

    call interaction_lists(Q0)

    nnzb = 2*sum(nnb) + Natom

    if (size(Gb, 3) < nnzb) then
        deallocate(Gb, Gcsr, UXY, UXYr, Gbja, Grja, GPb)
        allocate(Gb(3,3,nnzb), Gcsr(9*nnzb), UXY(3,3,nnzb), UXYr(9*nnzb), &
            Gbja(nnzb), Grja(9*nnzb), GPb(3,3,nnzb))
    end if

    call init_sparse_pattern(Q0)

    call init_gaussians(Q0, TAUMIN)

    T=TAUMIN
    next_stop = 0.5d0*TAUMAX
    do
        DT = 1d-2*sqrt(T)
        if (T+DT > next_stop) then
            DT = next_stop - T
            T = next_stop
        else
            T = T + DT
        end if
        call RHSSspfm(DT, mm)
        if (T == next_stop) exit
    end do

    GF = 0
    do I1=1,Natom
        do J2=Gbia(I1),Gbia(I1+1)-1
            I2 = Gbja(J2)
            GF(3*(I1-1) + 1 : 3*I1, 3*(I2-1) + 1 : 3*I2) = Gb(:,:,J2)
        end do
    end do

    !write (45, '(F16.12)') reshape(GF, (/ 9*Natom**2 /) )
    call dpotrf('U', 3*Natom, GF, 3*Natom, info)

    logdetg = 0.0
    s = 1
    DO j=1,3*Natom
        if (GF(j,j) < 0d0 ) s = -s
        logdetg = logdetg + 2.0*LOG(ABS(GF(j,j)))
    ENDDO

    if (s<0) then
        write(*,*) '||G|| is negative!'
        stop
    end if
!    logdetg = ll()
!    if (logdetg == 0.d0) then
!        write (45, '(F16.12)') reshape(GF, (/ 9*Natom**2 /) )
!        stop
!    end if
!
    write (*,*) gama
    W=-(1/TAUMAX)*(2.0*gama - 0.5*logdetg)! - 3.0*Natom*log(2.0*sqrt(M_PI)))
    !write (*,*) gama
END SUBROUTINE

!function ll()
!    double precision :: ll
!    interface
!        real(C_DOUBLE) function det_sparse_g(G, NNBFM, FMIDX, FMDIAG, Natom, &
!                                        NNBFMMAX) BIND(C)
!            use, intrinsic :: iso_c_binding
!            !real(C_DOUBLE) :: G(3,3,NNBFMMAX,Natom)
!            !integer(C_INT) :: NNBFM(Natom), FMIDX(NNBFMMAX, Natom), FMDIAG(Natom)
!            type(C_PTR), value :: G, NNBFM, FMIDX, FMDIAG
!            integer(C_INT), value :: Natom, NNBFMMAX
!        end function
!    end interface
!
!    ll = det_sparse_g(C_LOC(G), C_LOC(NNBFM), C_LOC(FMIDX), C_LOC(FMDIAG), &
!        Natom, size(G, 3))
!end function

!function logdet_g()
!    double precision :: logdet_g
!    double precision :: ddum(3*Natom)
!    integer :: idum(3*Natom)
!    integer :: maxfct, mnum, phase, nrhs, msglvl, ierr
!    integer :: I1, J2, k, j, nnz, curpos
!
!    interface
!        subroutine get_cholmod_ptrs(N, nnz, c_ia, c_ja, c_A), BIND(C)
!            integer(C_INT), value :: N, nnz
!            TYPE(C_PTR) :: c_ia, c_ja, c_A
!        end subroutine
!    end interface
!            
!    nnz = 9*(sum(nnbfm) - Natom)/2 + 6*Natom
!    
!    call get_cholmod_ptrs(3*Natom, nnz, c_ia, c_ja, c_A)
!    call c_f_pointer(c_ia, ia, shape=[3*Natom+1])
!    call c_f_pointer(c_ja, ja, shape=[nnz])
!    call c_f_pointer(c_A, A, shape=[nnz])
!    
!    curpos = 1
!    do I1=1,Natom
!        do k=1,3
!            ia(3*(I1-1) + k) = curpos
!
!            A(curpos:curpos+3-k) = G(k,k:3,FMdiag(I1), I1)
!            ja(curpos:curpos+3-k) = 3*(I1-1) + (/ (j, j=k,3) /)
!            curpos = curpos + 4 - k
!
!            do J2=FMdiag(I1)+1, NNBFM(I1)
!                A(curpos:curpos+2) = G(k,:,J2, I1)
!                ja(curpos:curpos+2) = 3*(FMIDX(J2,I1)-1) + (/ (j, j=1,3) /)
!                curpos = curpos + 3
!            end do
!        end do
!    end do
!    ia(3*Natom + 1) = curpos
!
!    ia = ia - 1
!    ja = ja - 1
!    
!    logdet_g = det_sparse_g(c_ia, c_ja, c_A)
!
!
!!    if (curpos /= nnz + 1) then
!!        write (*,*) 'Copy error: curpos != nnz + 1'
!!        stop
!!    end if
!!
!!    maxfct = 1
!!    mnum = 1
!!    phase = 12
!!    nrhs = 1
!!    msglvl  = 1
!!
!!    CALL pardiso_chkmatrix  (mtype, 3*Natom, a, ia, ja, ierr);
!!      IF (ierr .NE. 0) THEN
!!         WRITE(*,*) 'The following ERROR was detected: ', ierr
!!        STOP
!!    ENDIF
!!
!!    iparm(33) = 1
!!    CALL pardiso (pt, maxfct, mnum, mtype, phase, 3*Natom, a(1:nnz), ia, ja(1:nnz),  &
!!                   idum, nrhs, iparm, msglvl, ddum, ddum, ierr, dparm)
!!    logdet_sparse_g = dparm(33)
!!    if (dparm(33) < 0d0) then
!!        logdet_sparse_g = - abs(dparm(33))
!!    end if
!!    write (*,*) 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx', maxval(A(1:nnz))
!!    phase = -1
!!    CALL pardiso (pt, maxfct, mnum, mtype, phase, 3*Natom, a(1:nnz), ia, ja(1:nnz),  &
!!                   idum, nrhs, iparm, msglvl, ddum, ddum, ierr, dparm)
!end function

