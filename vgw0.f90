SUBROUTINE vgw0(Q0, BL_, beta,Ueff)
!$  use omp_lib
    IMPLICIT NONE
    double precision, intent(in) :: Q0(:,:), beta, BL_
    double precision, intent(out) :: Ueff
    real*8 :: LOGDET, logrho, T, TSTOP
    integer :: i, j, ncalls

    double precision, allocatable :: Y(:), RWORK(:), ATOL(:)
    integer, allocatable :: IWORK(:)

    integer :: tid = 0
    integer :: NEQ, ITOL, ITASK, IOPT, MF, ISTATE, LRW, LIW
    double precision :: RTOL

    Natom = size(Q0, 2)
    BL = BL_

    NEQ = 9*Natom + 1

    LRW = 20 + 16*NEQ
    LIW = 30
    allocate(Y(NEQ), ATOL(NEQ), RWORK(LRW), IWORK(LIW))

    ITOL=2
    RTOL=0
    ATOL(1:3*Natom) = vgw_atol(1)
    ATOL(3*Natom+1:3*Natom+6*Natom)=vgw_atol(2)
    ! tolerance for Qnkp and gamakp
    ATOL(NEQ) = vgw_atol(3)
    ITASK=1
    ISTATE=1
    IOPT = 1
    MF=10


    IWORK=0
    IWORK(6) = 50000 !MXSTEP

    RWORK(5)=dt0
    RWORK(6)=dtmax
    RWORK(7)=dtmin
    
    T = 0
    TSTOP = 0.5d0*beta

    y(1 : 3*Natom) = reshape(q0, (/3*Natom/))
    y(3*Natom+1:) = 0d0

    call presort_ppc(y(1:3*Natom), 8)
    if (pbc) then
        Ulrc = Ulrc * Natom**2 / BL**3
        UXXlrc = Ulrc * Natom**2 / BL**3
    else
        Ulrc = 0.0
        UXXlrc = 0.0
    end if

    dlsode_done=.FALSE.
    !print *, 'OAK',   omp_get_thread_num(), omp_get_wtime()
    tid = 0
!$omp parallel private(tid)

    call interaction_lists(y(1:3*Natom))

!$    tid = omp_get_thread_num()

    rhss_done = .FALSE.

    if (tid == 0) then
        !print *, 'master entry', omp_get_thread_num(), omp_get_wtime()
        call gaussian_average_avx_init(Natom, nnb, nbidx, nnbmax, LJA, LJC, NGAUSS, bl)
        CALL DLSODE(RHSS0,NEQ,Y,T,TSTOP,ITOL,RTOL,ATOL,ITASK,ISTATE,IOPT,&
            RWORK,LRW,IWORK,LIW,JAC,MF)
        dlsode_done = .TRUE.
        call rhss0(NEQ,0d0,Y,Y)
        !print *, 'master exit', omp_get_thread_num(), omp_get_wtime()
    else
        !print *, 'Entering tid =', omp_get_thread_num(), dlsode_done, omp_get_wtime()
        do
            call rhss0(NEQ,0d0,Y,Y)
            if (rhss_done) then
                !print *, 'Exiting tid=', tid, rhss_done, omp_get_wtime()
                exit
            end if
        end do
    end if

!$omp end parallel
    !print *, 'if exit tid=',  omp_get_thread_num(), omp_get_wtime()

    LOGDET=0d0
    DO j=3*Natom+1,9*Natom,6
        LOGDET = LOGDET + LOG( DETM_S(y(j : j+5)) )
    ENDDO
    print *, 'p', LOGDET, y(NEQ)

    logrho = 2.0*Natom*y(NEQ) - 0.5*LOGDET - 1.5*Natom*log(4.0*M_PI)
    ncalls = IWORK(12)
    !write (*,*) IWORK(11), 'steps,', IWORK(12), ' RHSS calls, logdet =', logdet

    deallocate(y, RWORK, IWORK, ATOL)
    call gaussian_average_avx_cleanup()

    Ueff = -logrho/beta
END SUBROUTINE


subroutine presort_ppc(r, nppc)
    use utils, only: index_sort
    implicit none
    double precision, intent(inout) :: r(3,Natom)
    integer, intent(in) :: nppc

    integer(c_size_t) :: i, N
    integer(c_int) :: idx(Natom), nunits
    real(c_double) :: z_(Natom), r_(3,Natom)
    double precision :: bu, cbl, ll(3), ur(3), sysbox

    ll = minval(r, 2)
    ur = maxval(r, 2)
    sysbox = maxval(ur - ll)

    cbl = bl

    ! cluster
    if (bl > 2*sysbox) then
        nunits = nint(2.0 * (0.75*Natom/(M_PI*nppc))**(1.0/3.0))
        cbl = sysbox
        pbc = .FALSE.
    else
        nunits = nint((real(Natom)/real(nppc))**(1.0/3.0))
        pbc = .TRUE.
    end if

    bu = cbl / nunits

    
    forall (i=1:Natom) idx(i) = i
    
    z_ = ( floor((r(1,:) - ll(1)) / bu) * nunits + floor( (r(2,:) - ll(2)) / bu ) ) * cbl + r(3,:) - ll(3)    
    
    N = Natom
    call index_sort(N, idx, z_)

    r_ = r(:,idx)
    r = r_
end subroutine presort_ppc

subroutine adjust_nbidx_for_c()
    integer :: i
    do i=1,Natom
        nbidx(1:nnb(i),i) =  nbidx(1:nnb(i),i) - 1
    end do
end subroutine


subroutine JAC()
end subroutine
