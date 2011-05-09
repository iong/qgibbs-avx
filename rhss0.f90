SUBROUTINE RHSS0(NEQ, DT, y, yp)
    use omp_lib
    IMPLICIT NONE
    integer, intent(in) :: NEQ
    REAL*8, intent(in) :: DT
    REAL*8, intent(inout) :: y(NEQ)
    REAL*8, intent(out) :: yp(NEQ)
    INTEGER I, J,I1,I2,IG,CNT,CNT2, NN1
    REAL*8 AG(3,3),GU(3,3), &
            DETA,DETAG,GUG(3,3),UX(3,nnbmax),UXX(3,3, nnbmax),QZQ,EXPAV, &
            BL2, G12(3,3),A(3,3), &
            Zq(3), Z(3,3),Q12(3), v0
    real*8 :: Q(3,Natom), G(3,3,Natom), GC(3,3,nnbmax), QC(3, nnbmax), UXX0(3,3), UX0(3), UPV_I1(3), UPM_I1(3,3)
    real*8 :: UPV_local(3,thread_start:thread_stop), UPM_local(3,3,thread_start:thread_stop), Ulocal

    if (NEQ /= (1+9*Natom) ) then
        write (*,*) 'NEQ != 1+9*Natom', NEQ, 1+9*Natom
        stop
    endif

    Q = reshape(y(2:1+3*Natom), (/3, Natom/) )
    call unpack_g(y, G)
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
            Q12 = min_image(Q(:,I1) - QC(:,I2), bl)
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
    do i1=1,nthreads-1
        UPV_local(:,thread_start:thread_stop) = UPV_local(:,thread_start:thread_stop) + UPV(:,thread_start:thread_stop,i1)
        UPM_local(:,:,thread_start:thread_stop) = UPM_local(:,:,thread_start:thread_stop) + UPM(:,:,thread_start:thread_stop,i1)
    end do
    TRUXXG(tid) = sum(UPM_local*G(:,:,thread_start:thread_stop))
    U(tid) = Ulocal
    do i1=thread_start,thread_stop
        cnt = 3*(i1-1)+2
        cnt2 = 2 + 6*(i1-1) + 3*Natom

        yp(CNT:CNT+2) = - matmul(G(:,:,I1), UPV_local(:,I1))

        GU = matmul(G(:,:,I1), UPM_local(:,:,I1))
        GUG = -matmul(GU, G(:,:,I1))
        do J=1,3
            GUG(J,J) = GUG(J,J) + invmass
        end do
        yp(cnt2  :cnt2+2) = GUG(:,1)
        yp(cnt2+3:cnt2+4) = GUG(2:3,2)
        yp(cnt2+5) = GUG(3,3)
        y(CNT:CNT+2) = y(CNT:CNT+2) + DT*yp(CNT:CNT+2)
        y(CNT2:CNT2+5) = y(CNT2:CNT2+5) + DT*yp(CNT2:CNT2+5)
    ENDDO

!$OMP BARRIER
!$OMP MASTER
            yp(1) = -0.25d0 * sum(TRUXXG) - sum(U)
            y(1) = y(1)  + yp(1) *DT
!$OMP END MASTER
END SUBROUTINE RHSS0
