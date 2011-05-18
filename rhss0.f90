SUBROUTINE RHSS0(DT, mm)
    use omp_lib
    IMPLICIT NONE
    REAL*8, intent(in) :: DT
    logical, intent(in) :: mm
    INTEGER :: J,I1,I2,IG,CNT,CNT2, NN1
    REAL*8 AG(3,3),GU(3,3), &
            DETA,DETAG,GUG(3,3),UX(3,nnbmax),UXX(3,3, nnbmax),QZQ,EXPAV, &
            G12(3,3),A(3,3), &
            Zq(3), Z(3,3),Q12(3), v0, QP_(3)
    real*8 :: GC(3,3,nnbmax), &
            QC(3, nnbmax), UXX0(3,3), UX0(3), UPV_I1(3), UPM_I1(3,3)
    real*8 :: UPV_local(3,thread_start:thread_stop), &
            UPM_local(3,3,thread_start:thread_stop), Ulocal

    UPM(:,:,:,tid) = 0d0
    UPV(:,:,tid) = 0d0
    Ulocal = 0d0
!$OMP DO SCHEDULE(DYNAMIC)
    do I1=1,Natom-1
        NN1 = NNB(I1)
        if (NN1 == 0) cycle

        UPV_I1 = UPV(:,I1,tid) 
        UPM_I1 = UPM(:,:,I1,tid) 
        
        qc(:,1:NN1) = q(:,nbidx(1:NN1,I1))
        GC(:,:,1:NN1) = G(:,:,nbidx(1:NN1,I1))
        DO I2=1,NN1
            Q12 = Q(:,I1) - QC(:,I2)
            Q12 = min_image(Q12, bl)
            G12=G(:,:,I1)+GC(:,:,I2)

            call detminvm(G12, DETA, A)
            DETA = 1.0d0/DETA

            UX0 = 0d0; UXX0 = 0d0
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

                EXPAV=SQRT(DETA/DETAG)*EXP(-qZq)
                Ulocal=Ulocal+EXPAV*LJC(IG)

                v0 = 2d0*expav*LJC(IG)

                UX0 = UX0 - v0*Zq
                do J=1,3
                    UXX0(:,J) = UXX0(:,J) + v0*(2d0*Zq*Zq(J) - Z(:,J))
                end do
            ENDDO ! IG
! Avoid load and store as much as possbile. Store now, process as a stream later. Much faster.
            UPV_I1 = UPV_I1 + UX0
            UX(:,I2) = - UX0
            UPM_I1 = UPM_I1 + UXX0
            UXX(:,:,I2) = UXX0
        ENDDO ! I2
        UPV(:,I1,tid) = UPV_I1
        UPM(:,:,I1,tid) = UPM_I1
        UPV(:,nbidx(1:NN1,I1),tid) = UPV(:,nbidx(1:NN1,I1),tid) + UX(:,1:NN1)
        UPM(:,:,nbidx(1:NN1,I1),tid) = UPM(:,:,nbidx(1:NN1,I1),tid) + UXX(:,:,1:NN1)
    ENDDO ! I
!$OMP END DO

    UPV_local(:,thread_start:thread_stop) = UPV(:,thread_start:thread_stop,0)
    UPM_local(:,:,thread_start:thread_stop) = UPM(:,:,thread_start:thread_stop,0)
    do i1=1,nthr-1
        UPV_local(:,thread_start:thread_stop) = UPV_local(:,thread_start:thread_stop) + UPV(:,thread_start:thread_stop,i1)
        UPM_local(:,:,thread_start:thread_stop) = UPM_local(:,:,thread_start:thread_stop) + UPM(:,:,thread_start:thread_stop,i1)
    end do
    TRUXXG(tid) = sum(UPM_local*G(:,:,thread_start:thread_stop))
    U(tid) = Ulocal   
    do i1=thread_start,thread_stop
        QP_ = - matmul(G(:,:,I1), UPV_local(:,I1))

        GU = matmul(G(:,:,I1), UPM_local(:,:,I1))
        GUG = -matmul(GU, G(:,:,I1))
        do J=1,3
            GUG(J,J) = GUG(J,J) + invmass
        end do
        
        
        Q(:,I1) = Q(:,I1) + DT * QP_
        QP(:,I1) = QP_
        
        GP(:,:,I1) = GUG
        G(:,:,I1) = G(:,:,I1) + DT*GUG
        
        if (mm) then
            gamak(:, I1) = gamak(:, I1) - DT * matmul(transpose(Qnk(:,:,i1)), UPV_local(:,I1))

    	    Qnk(:,:,I1) = Qnk(:,:,I1) - DT * matmul(GU, Qnk(:,:,i1))
	    end if
    end do

!$OMP BARRIER
!$OMP MASTER
    gamap = -0.25d0 * sum(TRUXXG(0:nthr-1)) - sum(U(0:nthr-1))
    gama = gama + DT * gamap
!$OMP END MASTER
END SUBROUTINE RHSS0
