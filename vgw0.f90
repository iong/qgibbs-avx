SUBROUTINE vgw0v(Q0, W, TAUMAX,TAUI, Y)
    use propagation
    use utils
    use omp_lib
    IMPLICIT NONE
    REAL*8, intent(in) :: Q0(3,Natom), TAUMAX(:), TAUI
    REAL*8, intent(inout) :: Y(1+9*Natom)
    REAL*8, intent(out) :: W(:)
    real*8 :: G(3,3), T, LOGDET, taustops(1), waste(1)
    real*8 :: DT, YP(1+9*Natom), next_stop
    integer :: i, j, NEQ, chunk_size


    NEQ = 1 + 9*Natom

!$OMP PARALLEL PRIVATE(I,J, G, T, DT, next_stop)
    if (TAUI <= 0.0d0) then
        T = TAUMIN
    else
        T = 0.5*TAUI
    endif
 
    tid = OMP_get_thread_num()
    nthreads = OMP_get_num_threads()

    chunk_size = Natom/nthreads
    thread_start=tid*chunk_size + 1
    thread_stop=(tid+1)*chunk_size
    if (tid==nthreads-1) then
       thread_stop = Natom
    end if

!$OMP SINGLE
        if (.not.allocated(upv)) then
            allocate(upv(3,Natom,0:nthreads-1), upm(3,3,Natom,0:nthreads-1), &
            TRUXXG(0:nthreads-1), U(0:nthreads-1))
        end if
!$OMP END SINGLE

    if (TAUI <= 0.0d0) then
        call interaction_lists(Q0) !Determine which particles
        call init_gaussians(Q0, TAUMIN, y)
    endif
!$OMP BARRIER
    do i=1,size(TAUMAX)
        next_stop = 0.5d0*TAUMAX(i)
        do
            DT = 1d-2*sqrt(T)
            if (T+DT > next_stop) then
                DT = next_stop - T
                T = next_stop
            else
                T = T + DT
            end if
            call RHSS0(NEQ, DT, Y, YP)
            if (T == next_stop) exit
        end do

!$OMP SINGLE
    LOGDET=0d0
!$OMP END SINGLE
!$OMP DO SCHEDULE(STATIC) REDUCTION(+:LOGDET)
        DO j=1,Natom
            call unpack_g_i(y, j, G)
            LOGDET = LOGDET + LOG(DETM(G))
        ENDDO
!$OMP END DO


        !The constant should be irrelevant for MC and FX. don't forget to
        !add it to the Cv
        !W(i)=-(1/TAUMAX(i))*(2.0*y(1) - 0.5*LOGDET - 3.0*Natom*log(2.0*sqrt(M_PI)))
!$OMP MASTER
        W(i)=-(1/TAUMAX(i))*(2.0*y(1) - 0.5*LOGDET)
!$OMP END MASTER
    end do
!$OMP END PARALLEL
END SUBROUTINE


SUBROUTINE vgw0(Q0, W, TAUMAX,TAUI, Y)
IMPLICIT NONE
REAL*8, intent(in) :: Q0(3,Natom), TAUMAX, TAUI
REAL*8, intent(inout) :: Y(1+9*Natom)
REAL*8, intent(out) :: W
real*8, dimension(1) :: W_, TAUMAX_

    TAUMAX_(1) = TAUMAX
    call vgw0v(Q0, W_, TAUMAX_, TAUI, Y)
    W = W_(1)
END SUBROUTINE

! vim:et
