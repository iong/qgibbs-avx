SUBROUTINE vgw0v(Q0, BL_, TAUMAX,TAUI, W, WX)
    use omp_lib
    IMPLICIT NONE
    REAL*8, intent(in) :: Q0(:,:), TAUMAX(:), TAUI, BL_
    REAL*8, intent(out) :: W(:)
    real*8, intent(out), optional :: WX(:,:)
    real*8 :: LOGDET
    real*8 :: DT, next_stop
    integer :: i, j, chunk_size
    logical :: mm = .FALSE.


    Natom = size(Q0, 2)
    BL = BL_
    if (present(WX)) mm = .TRUE.

!$OMP PARALLEL PRIVATE(I,J, T, DT, next_stop)
    if (TAUI <= 0.0d0) then
        T = TAUMIN
    else
        T = 0.5*TAUI
    endif
 
    tid = OMP_get_thread_num()
    nthr = OMP_get_num_threads()

    chunk_size = Natom/nthr
    thread_start=tid*chunk_size + 1
    thread_stop=(tid+1)*chunk_size
    if (tid==nthr-1) then
       thread_stop = Natom
    end if


    if (TAUI <= 0.0d0) then
        call interaction_lists(Q0) !Determine which particles
        call init_gaussians(Q0, TAUMIN, mm)
    end if
    
!$OMP BARRIER
    do i=1,size(TAUMAX)
        next_stop = 0.5d0*TAUMAX(i)
        do
            DT = 1d-4!1d-2*sqrt(T)
            if (T+DT > next_stop) then
                DT = next_stop - T
                T = next_stop
            else
                T = T + DT
            end if
            call RHSS0(DT, mm)
            if (T == next_stop) exit
        end do

!$OMP SINGLE
    LOGDET=0d0
!$OMP END SINGLE
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:LOGDET)
        DO j=1,Natom
            LOGDET = LOGDET + LOG(DETM(G(:,:,j)))
        ENDDO
!$OMP END DO

!$OMP MASTER
        W(i)=-(1/TAUMAX(i))*(2.0*gama - 0.5*LOGDET)! - 3.0*Natom*log(2.0*sqrt(M_PI)))
!$OMP END MASTER
    end do
!$OMP END PARALLEL

    if (present(WX)) then
        WX(:,1:Natom) =  gamak(:,1:Natom) / T
    end if
END SUBROUTINE


SUBROUTINE vgw0(Q0, BL_, TAUMAX,TAUI, W, WX)
IMPLICIT NONE
REAL*8, intent(in) :: Q0(:,:), TAUMAX, TAUI, BL_
REAL*8, intent(out) :: W
real*8, intent(out), optional :: WX(:,:)

real*8, dimension(1) :: W_, TAUMAX_

    TAUMAX_(1) = TAUMAX
    if (present(WX)) then
        call vgw0v(Q0, BL_, TAUMAX_, TAUI, W_, WX)
    else
        call vgw0v(Q0, BL_, TAUMAX_, TAUI, W_)
    end if
    W = W_(1)
END SUBROUTINE

! vim:et
