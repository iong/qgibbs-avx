subroutine Upot_tau0(Q)
    use utils
    IMPLICIT NONE
    REAL*8, intent(in) :: Q(:,:)
    INTEGER  I,J,N
    real*8 :: RSQ,BL2,QIJ(3)

    N = size(Q, 2)
    BL2=BL/2.0

!$OMP SINGLE
        U=0d0
!$OMP END SINGLE
!$OMP DO SCHEDULE(DYNAMIC, 32) PRIVATE(I, J, QIJ, RSQ)
    DO I=1,N-1
        DO J=1,NNB(I)
                qij = min_image(Q(:,I) - Q(:,NBIDX(J, I)), BL)
                rsq = sum(qij**2)
                U(tid) = U(tid) + sum(LJC(1:NGAUSS)*EXP(-LJA(1:NGAUSS)*Rsq))
        ENDDO
    ENDDO
!$OMP END DO
end subroutine Upot_tau0


!subroutine U_Ux_tau0(Q, U, UX)
!    IMPLICIT NONE
!    REAL*8, intent(in) :: Q(:,:)
!    real*8, intent(in) :: U, UX(:)
!    real*8 :: UXI(3), UXJ(3,size(Q,2))
!    INTEGER  I,J,K,P, N, cnt
!    REAL*8  RSQ,BL2,VEC(3),GRAD
!
!    N = size(Q, 2)
!
!    BL2=BL/2
!
!    cnt = 1
!    UX=0.0
!    DO I=1,N-1
!        UXI = 0d0
!        UXJ = 0d0
!        DO J=I+1,N
!            if (QRC(J, I)) then
!                qij = min_image(Q(:,I) - Q(:,J), BL)
!                rsq = sum(qij**2)
!                LJC_EXPAV=LJC(1:NGAUSS)*EXP(-LJA(1:NGAUSS)*Rsq)
!                U = U + sum(LJC_EXPAV)
!                UX0 = qij*(2.0*sum(LJA(1:NGAUSS)*LJC_EXPAV))
!                UXI = UXI + UX0
!                UXJ(:,J) = -UX0
!            endif
!        ENDDO
!        UX(cnt:cnt+2) = UX(cnt:cnt+2) + UXI
!        UX(cnt+3:N3) = UX(cnt+3:N3) + reshape(UXJ(:,I+1:N), (/ 3*(N-I) /))
!        cnt = cnt + 3
!    ENDDO
!end subroutine U_Ux_tau0

! The writes the result is into U and UPV, which are shared across all threads
subroutine Upot_Ux_tau0(Q)
    use utils
    IMPLICIT NONE
    REAL*8, intent(in) :: Q(:,:)
    real*8 :: UX0(3), UXI(3), UXJ(3,size(Q,2)), QC(3,size(Q,2))
    INTEGER :: I,J, N, cnt
    REAL*8 :: RSQ,QIJ(3), LJC_EXPAV(NGAUSS)

    N = size(Q, 2)

    cnt = 1
!$OMP SINGLE
        U=0d0
        UPV=0d0
!$OMP END SINGLE
!$OMP DO SCHEDULE(DYNAMIC, 32) PRIVATE(I, J, UXI, UX0, UXJ, QC, QIJ, RSQ, LJC_EXPAV, CNT)
    DO I=1,N-1
        UXI =  UPV(:,I,tid)
        qc(:,1:NNB(I)) = q(:,nbidx(1:NNB(I), I))
        DO J=1,NNB(I)
           qij = min_image(Q(:,I) - QC(:,J), BL)
           rsq = sum(qij**2)
           LJC_EXPAV=LJC(1:NGAUSS)*EXP(-LJA(1:NGAUSS)*Rsq)
           U(tid) = U(tid) + sum(LJC_EXPAV)
           UX0 = qij*(2.0*sum(LJA(1:NGAUSS)*LJC_EXPAV))
           UXI = UXI + UX0
           UXJ(:,J) = -UX0
        ENDDO
        UPV(:,I,tid) =UXI
        UPV(:,nbidx(1:NNB(I),I),tid) = UPV(:,nbidx(1:NNB(I),I),tid) + UXJ(:,1:NNB(I))
    ENDDO
!$OMP END DO

    do i=1,nthreads-1
        UPV(:,thread_start:thread_stop,0) = UPV(:,thread_start:thread_stop,0) + UPV(:,thread_start:thread_stop,i)
    end do
!$OMP BARRIER
end subroutine Upot_Ux_tau0
