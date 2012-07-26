SUBROUTINE vgw0(Q0, BL_, beta,Ueff, rt)
    IMPLICIT NONE
    double precision, intent(in) :: Q0(:,:), beta, BL_
    double precision, intent(out) :: Ueff
    double precision, intent(out), optional :: rt
    real*8 :: LOGDET, logrho, TSTOP, start_time, stop_time
    integer :: i, j, ncalls

    double precision, allocatable :: Y(:), RWORK(:), YP(:), ATOL(:)
    integer, allocatable :: IWORK(:)

    integer :: NEQ, ITOL, ITASK, IOPT, MF, ISTATE, LRW, LIW
    double precision :: RTOL
    double precision, pointer, dimension(:, :) :: G

    Natom = size(Q0, 2)
    BL = BL_

    NEQ = 9*Natom + 1
    !if (present(WX)) NEQ = NEQ + 12*Natom

    LRW = 20 + 16*NEQ
    LIW = 30

    allocate(Y(NEQ), YP(NEQ), ATOL(NEQ), RWORK(LRW), IWORK(LIW))

    y(        1 :   Natom) = Q0(1,:)
    y(  Natom+1 : 2*Natom) = Q0(2,:)
    y(2*Natom+1 : 3*Natom) = Q0(3,:)

    y(3*Natom+1:) = 0d0

    call       presort_ppc(y(1:Natom), y(Natom+1:2*Natom), y(2*Natom+1:3*Natom), 4)
    call interaction_lists(y(1:Natom), y(Natom+1:2*Natom), y(2*Natom+1:3*Natom))
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

    call cpu_time(start_time)
    TSTOP = 0.5d0*beta
    CALL DLSODE(RHSS0,NEQ,Y,T,TSTOP,ITOL,RTOL,ATOL,ITASK,ISTATE,IOPT,&
        RWORK,LRW,IWORK,LIW,JAC,MF)


    LOGDET=sum(log(pdetm_s(reshape(y(3*Natom+1:9*Natom), (/Natom, 6/)))))
    call cpu_time(stop_time)

    logrho = 2.0*Natom*y(NEQ) - 0.5*LOGDET - 1.5*Natom*log(4.0*M_PI)
    ncalls = IWORK(12)
    !write (*,*) IWORK(11), 'steps,', IWORK(12), ' RHSS calls, logdet =', logdet

    deallocate(y, yp, RWORK, IWORK, ATOL)

    Ueff = -logrho/beta
    if (present(rt)) then
        rt = (stop_time - start_time) / real(ncalls)
     end if
END SUBROUTINE


subroutine presort_ppc(x, y, z, nppc)
    use utils, only: index_sort
    implicit none
    double precision, intent(inout) :: x(:), y(:), z(:)
    integer, intent(in) :: nppc

    integer(c_size_t) :: i, N
    integer(c_int) :: idx(size(x)), nunits
    real(c_double) :: z_(size(x))
    double precision :: bu, cbl, minx, miny, minz, maxx, maxy, maxz, sysbox

    minx = minval(x); miny = minval(y); minz = minval(z) 
    maxx = maxval(x); maxy = maxval(y); maxz = maxval(z) 
    sysbox = max(maxx-minx, maxy-miny, maxz-minz)

    N = size(x)
    cbl = bl

    ! cluster
    if (bl > 2*sysbox) then
        nunits = nint(2.0 * (0.75*N/(M_PI*nppc))**(1.0/3.0))
        cbl = sysbox
    else
        nunits = nint((real(N)/real(nppc))**(1.0/3.0))
    end if

    bu = cbl / nunits

    
    forall (i=1:N) idx(i) = i
    
    z_ = ( floor((x - minx) / bu) * nunits + floor( (y - miny) / bu ) ) * cbl + z - minz    
    
    call index_sort(N, idx, z_)

    z_ = x(idx); x = z_
    z_ = y(idx); y = z_
    z_ = z(idx); z = z_
end subroutine presort_ppc


subroutine JAC()
end subroutine
