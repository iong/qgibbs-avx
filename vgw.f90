module vgw
    use utils
    implicit none
    integer, private :: Natom
    real*8, private :: BL
    real*8, dimension(10) :: LJA, LJC
    integer :: NGAUSS
    integer, allocatable :: NBIDX(:,:), NNB(:)
    real*8, allocatable :: UPV(:,:,:), UPM(:,:,:,:), TRUXXG(:), U(:)
    real*8, pointer :: yp_ptr(:)
    real*8 :: invmass, RC, ATOL, RTOL, TAUMIN, mass
    logical :: finished
    integer :: tid=0, nthreads=1, thread_start, thread_stop, nnbmax
!$OMP THREADPRIVATE(tid, nthreads, thread_start, thread_stop, nnbmax)

contains

subroutine vgwinit(N, boxlen)
    integer, intent(in) :: N
    real*8, intent(in) :: boxlen
    Natom = N
    BL = boxlen
    invmass=1.0D0/mass

    allocate(NNB(Natom), NBIDX(natom,natom))
end subroutine



subroutine unpack_q(y, q)
    real*8, intent(in) :: y(:)
    real*8, intent(out) :: q(:,:)
    integer :: Natom

    Natom = size(q,2)
    q = reshape(y(2:1+3*Natom), (/3, Natom/) )
end subroutine

subroutine unpack_g(y, G)
    real*8, intent(in) :: y(:)
    real*8, intent(out) :: G(:,:,:)
    integer :: i, j, k, cnt, N

    N = size(G, 3)
    cnt = 1+3*N+1
    do i=1,N
        do j=1,3
            do k=j,3
                G(k, j, i) = y(cnt)
                G(j, k, i) = y(cnt)
                cnt = cnt + 1
            enddo
        enddo
    enddo
end subroutine

subroutine unpack_g_i(y, i, G)
    real*8, intent(in) :: y(:)
    integer, intent(in) :: i
    real*8, intent(out) :: G(3,3)
    integer :: cnt

    cnt = 2 + 3*Natom + 6*(i-1)
    G(:,1) = y(cnt:cnt+2)
    G(2:3,2) = y(cnt+3:cnt+4)
    G(3,3) = y(cnt+5)

    G(1,2) = G(2,1)
    G(1:2,3) = G(3,1:2)
end subroutine

subroutine unpack_Qnk(y, Qnk)
    real*8, intent(in) :: y(:)
    real*8, intent(out) :: Qnk(:,:,:)
    integer :: i, j, k, CNT, N

    N = size(Qnk, 3)

    CNT = 2+9*N
    DO I=1,N
        DO J=1,3
            DO K=1,3
                QNK(J,K,I)=Y(CNT)
                CNT=CNT+1
            ENDDO
        ENDDO
    ENDDO
end subroutine

function unpack_f(y) result(f)
    real*8, intent(in) :: y(:)
    real*8 :: f(3,Natom)

    f = reshape(y(2+18*Natom:1+21*Natom), (/3, Natom/) )
end function

subroutine init_gaussians(q0, tau, y)
    REAL*8, intent(in) :: Q0(3,Natom), tau
    REAL*8, intent(inout) :: Y(:)
    real*8 :: G0(6)
    integer :: NY, J

    NY = size(y)

    if (NY == 1 + 21*Natom) then
        call Upot_UX_tau0(Q0)
    else
        call Upot_tau0(Q0)
    end if

!$OMP SINGLE
    Y(1) = -tau*sum(U)
    Y(2:1+3*Natom)=reshape(Q0, (/3*Natom/))

    G0=tau*invmass*(/1.0, 0.0, 0.0, 1.0, 0.0, 1.0/)
    Y(2+3*Natom:1+9*Natom) = reshape(spread(G0, 2, Natom), (/ 6*Natom /) )
    if (NY == 1 + 21*Natom) then
        do j=2+9*Natom,1+18*Natom,9
            Y(j:j+8) = (/1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0/)
        end do
        Y(2+18*Natom:1+21*Natom) = tau * reshape(UPV(:,:,0), (/3*Natom/) )
    end if
!$OMP END SINGLE
end subroutine

include 'interaction_lists.f90'
include 'potential_energy.f90'
include 'vgw0.f90'
!include 'vgw1.f90'
include 'rhss0.f90'
!include 'rhss1.f90'

end module vgw

