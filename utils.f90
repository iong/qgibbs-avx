module utils
    double precision :: M_PI = 3.14159265358979323846264338327950288d0
    integer, parameter :: RP = 4
    interface min_image
        module procedure min_image1
        module procedure min_image2
    end interface

    interface
        subroutine index_sort(N, idx, data) bind(c)
            use iso_c_binding
            integer(c_size_t), value :: N
            integer(c_int) :: idx(*)
            real(c_double) :: data(*)
        end subroutine
        subroutine pinvdetinv_avx(N, stride, a, inva, detinva) bind(c)
            use iso_c_binding
            integer(c_int), value :: N, stride
            real(c_float) :: a(*)
            real(c_float) :: inva(*), detinva(*)
        end subroutine
        subroutine posix_memalign(ptr, align, nbytes) bind(c)
            use iso_c_binding
            type(c_ptr) :: ptr
            integer(c_size_t), value :: align, nbytes
        end subroutine
        subroutine posix_free(ptr) bind(c, name='free')
            use iso_c_binding
            type(c_ptr), value :: ptr
        end subroutine
    end interface
contains
subroutine seed_rng()
        integer :: seed(128), seed_size, c, crate, cmax

        call random_seed(SIZE=seed_size)
!        write(*,*) 'Seed size =', seed_size
        call random_seed(GET=seed(1:seed_size))
        call system_clock(c, crate, cmax)
        seed(1) = mod(c,crate)
        call random_seed(PUT=seed(1:seed_size))
end subroutine

double precision function gaussran(sigma, x0) result(y)
    implicit none
    double precision, intent(in) :: sigma, x0
    double precision :: x(2)
    do
    call random_number(x)
        if (x(1) /= 0.0d0) exit
    enddo
    !   write (*,*) x
    y = sqrt( -2.0 * log(x(1))) * cos(2.0*M_PI*x(2))

    y = y*sigma + x0
end function gaussran

subroutine pol2cart(pol, cart)
    double precision, intent(in) :: pol(3)
    double precision, intent(out) :: cart(3)
    double precision :: rxy

    rxy = pol(1) * sin(pol(2))
    cart(1) = rxy * cos(pol(3))
    cart(2) = rxy * sin(pol(3))
    cart(3) = pol(1) * cos(pol(2))
end subroutine pol2cart


subroutine int2strz(n, w, str0)
    implicit none
    integer, intent(in) :: n, w
    character(len=*), intent(out) :: str0
    integer :: z, n2, i

    if (len(str0) < w) then
        write (*,*) 'Number too large for string'
        stop
    end if

    n2 = n
    z=iachar('0')
    do i=w,1,-1
        str0(i:i) = achar(z + mod(n2,10))
        n2 = n2/10
    end do
end subroutine int2strz

subroutine linspace(xmin, xmax, N, xout)
    double precision, intent(in) :: xmin, xmax
    integer, intent(in) :: N
    double precision, intent(out) :: xout(N)
    integer :: i
    double precision:: dx

    dx = (xmax - xmin) / (N-1)
    xout = xmin + dx * (/(i,i=0,N-1)/)
end subroutine linspace

subroutine replace_char(str, a, b)
    character(len=*), intent(inout) :: str
    character, intent(in) :: a, b
    integer :: i

    do i=1,len_trim(str)
        if (str(i:i)==a) then
            str(i:i) = b
        end if
    end do
end subroutine replace_char

function fliplr(v) result(y)
    double precision, intent(in) :: v(:)
    double precision :: y(size(v))
    integer :: i, N

    N = size(v)
    do i=1,N
        y(N-i+1) = v(i)
    end do
end function

subroutine detminvm(A, DETA, INVA)
    implicit none
    double precision, intent(in) :: A(3,3)
    double precision, intent(out) :: DETA, INVA(3,3)
    double precision :: INVDET

    INVA(1,1) = A(2,2)*A(3,3)-A(2,3)*A(3,2)
    INVA(2,1) = -A(1,2)*A(3,3)+A(1,3)*A(3,2)
    INVA(3,1) = A(1,2)*A(2,3)-A(1,3)*A(2,2)
    INVA(1,2) = INVA(2,1)
    INVA(2,2) = A(1,1)*A(3,3)-A(1,3)*A(3,1)
    INVA(3,2) = -A(1,1)*A(2,3)+A(1,3)*A(2,1)
    INVA(1,3) = INVA(3,1)
    INVA(2,3) = INVA(3,2)
    INVA(3,3) = A(1,1)*A(2,2)-A(1,2)*A(2,1)

    DETA = INVA(1,1)*A(1,1)+INVA(2,1)*A(2,1)+INVA(3,1)*A(3,1)
    INVDET = 1.0d0/DETA
    INVA = INVA * INVDET
end subroutine detminvm

subroutine detminvm_sg(A, DETA, INVA)
    implicit none
    double precision, intent(in) :: A(6)
    double precision, intent(out) :: DETA, INVA(3,3)
    double precision :: INVDET

    INVA(1,1) = A(4)*A(6)-A(5)**2
    INVA(2,1) = -A(2)*A(6)+A(3)*A(5)
    INVA(3,1) = A(2)*A(5)-A(3)*A(4)
    INVA(1,2) = INVA(2,1)
    INVA(2,2) = A(1)*A(6)-A(3)**2
    INVA(3,2) = -A(1)*A(5)+A(3)*A(2)
    INVA(1,3) = INVA(3,1)
    INVA(2,3) = INVA(3,2)
    INVA(3,3) = A(1)*A(4)-A(2)**2

    DETA = INVA(1,1)*A(1)+INVA(2,1)*A(2)+INVA(3,1)*A(3)
    INVDET = 1.0d0/DETA
    INVA = INVA * INVDET
end subroutine detminvm_sg

subroutine pdetminvm_sg(N, A, DETINVA, INVA)
    implicit none
    integer, intent(in) :: N
    real(RP), intent(in) :: A(:,:)
    real(RP), intent(out) :: DETINVA(:), INVA(:,:)
    integer :: I

    INVA(1:N,1) =  A(1:N,4)*A(1:N,6) - A(1:N,5)**2
    INVA(1:N,2) = -A(1:N,2)*A(1:N,6) + A(1:N,3)*A(1:N,5)
    INVA(1:N,3) =  A(1:N,2)*A(1:N,5) - A(1:N,3)*A(1:N,4)
    INVA(1:N,4) =  A(1:N,1)*A(1:N,6) - A(1:N,3)**2
    INVA(1:N,5) = -A(1:N,1)*A(1:N,5) + A(1:N,3)*A(1:N,2)
    INVA(1:N,6) =  A(1:N,1)*A(1:N,4) - A(1:N,2)**2

    DETINVA(1:N) = 1.0_RP/(INVA(1:N,1)*A(1:N,1) +INVA(1:N,2)*A(1:N,2) +INVA(1:N,3)*A(1:N,3))
    do I=1,6
        INVA(1:N,I) = INVA(1:N,I) * DETINVA(1:N)
    end do
end subroutine pdetminvm_sg

subroutine pdetminvm_sg2(N, A, DETINVA, INVA)
    implicit none
    integer, intent(in) :: N
    real(RP), intent(in) :: A(:,:)
    real(RP), intent(out) :: DETINVA(:), INVA(:,:)
    integer :: I, J, N8
    real(RP) :: DETIA(8), IA(8,6)

    N8 = (N/8)*8
    do i=1,N8,8
        IA(:,1) =  A(i:i+7,4)*A(i:i+7,6) - A(i:i+7,5)**2
        IA(:,2) = -A(i:i+7,2)*A(i:i+7,6) + A(i:i+7,3)*A(i:i+7,5)
        IA(:,3) =  A(i:i+7,2)*A(i:i+7,5) - A(i:i+7,3)*A(i:i+7,4)
        IA(:,4) =  A(i:i+7,1)*A(i:i+7,6) - A(i:i+7,3)**2
        IA(:,5) = -A(i:i+7,1)*A(i:i+7,5) + A(i:i+7,3)*A(i:i+7,2)
        IA(:,6) =  A(i:i+7,1)*A(i:i+7,4) - A(i:i+7,2)**2

        DETIA = 1.0_RP / ( IA(:,1)*A(i:i+7,1) +IA(:,2)*A(i:i+7,2) &
            + IA(:,3)*A(i:i+7,3) )

        do J=1,6
            INVA(i:i+7,J) = IA(:,J) * DETIA
        end do
        DETINVA(i:i+7) = DETIA
    end do

    if (N8 == N) return

    INVA(N8+1:N,1) =  A(N8+1:N,4)*A(N8+1:N,6) - A(N8+1:N,5)**2
    INVA(N8+1:N,2) = -A(N8+1:N,2)*A(N8+1:N,6) + A(N8+1:N,3)*A(N8+1:N,5)
    INVA(N8+1:N,3) =  A(N8+1:N,2)*A(N8+1:N,5) - A(N8+1:N,3)*A(N8+1:N,4)
    INVA(N8+1:N,4) =  A(N8+1:N,1)*A(N8+1:N,6) - A(N8+1:N,3)**2
    INVA(N8+1:N,5) = -A(N8+1:N,1)*A(N8+1:N,5) + A(N8+1:N,3)*A(N8+1:N,2)
    INVA(N8+1:N,6) =  A(N8+1:N,1)*A(N8+1:N,4) - A(N8+1:N,2)**2

    DETINVA(N8+1:N) = 1.0_RP/(INVA(N8+1:N,1)*A(N8+1:N,1) +INVA(N8+1:N,2)*A(N8+1:N,2) +INVA(N8+1:N,3)*A(N8+1:N,3))
    do I=1,6
        INVA(N8+1:N,I) = INVA(N8+1:N,I) * DETINVA(N8+1:N)
    end do
end subroutine



subroutine detminvm_ss(A, DETA, INVA)
    implicit none
    double precision, intent(in) :: A(6)
    double precision, intent(out) :: DETA, INVA(6)
    double precision :: INVDET

    INVA(1) = A(4)*A(6)-A(5)**2
    INVA(2) = -A(2)*A(6)+A(3)*A(5)
    INVA(3) = A(2)*A(5)-A(3)*A(4)
    INVA(4) = A(1)*A(6)-A(3)**2
    INVA(5) = -A(1)*A(5)+A(3)*A(2)
    INVA(6) = A(1)*A(4)-A(2)**2

    DETA = sum(INVA(1:3) * A(1:3) )
    INVDET = 1.0d0/DETA
    INVA = INVA * INVDET
end subroutine detminvm_ss

pure function detm(A)
    implicit none
    double precision, intent(in) :: A(3,3)
    double precision :: DETM

    DETM = (A(2,2)*A(3,3) - A(2,3)*A(3,2)) * A(1,1)&
          + ( -A(1,2)*A(3,3) + A(1,3)*A(3,2) ) * A(2,1) &
          + ( A(1,2)*A(2,3) - A(1,3)*A(2,2) ) * A(3,1)
end function detm

pure function detm_s(A)
    implicit none
    double precision, intent(in) :: A(6)
    double precision :: DETM_S

    DETM_S = (A(4)*A(6) - A(5)**2) * A(1)&
          + ( -A(2)*A(6) + A(3)*A(5) ) * A(2) &
          + ( A(2)*A(5) - A(3)*A(4) ) * A(3)
end function detm_s

pure function pdetm_s(A)
    implicit none
    double precision, intent(in) :: A(:,:)
    double precision :: PDETM_S(size(A, 1))

    PDETM_S = (A(:,4)*A(:,6) - A(:,5)**2) * A(:,1)&
          + ( -A(:,2)*A(:,6) + A(:,3)*A(:,5) ) * A(:,2) &
          + ( A(:,2)*A(:,5) - A(:,3)*A(:,4) ) * A(:,3)
end function pdetm_s


function outer_product(l, r) result(m)
    double precision, intent(in) :: l(:), r(:)
    double precision :: m(size(l), size(r))
    integer :: i

    forall (i=1:size(r))
        m(:,i) = l*r(i)
    end forall
end function outer_product

pure function outer_product3(l, r) result(m)
    double precision, intent(in) :: l(3), r(3)
    double precision :: m(3,3)
    integer :: i

    forall (i=1:3)
        m(:,i) = l*r(i)
    end forall
end function outer_product3

pure function min_image1(r)
    implicit none
    real(8), intent(in) :: r(:)
    real(8) :: min_image1(size(r))
    integer :: i

    where (abs(r) > 0.5)
            min_image1 = r - sign(1.0d0, r)
    elsewhere
            min_image1 = r
    endwhere

!    do i=1,3
!        if (r(i) > 0.5) then
!            min_image1(i) = r(i) - 1.0
!        elseif (r(i) < -0.5) then
!            min_image1(i) = r(i) + 1.0
!        else
!            min_image1(i) = r(i)
!        end if
!    end do
end function

pure function min_image2(r, bl)
    implicit none
    real(8), intent(in) :: r(:), bl
    real(8) :: min_image2(size(r))
    real(8) :: bl2

    bl2=0.5d0*bl
    where (abs(r) > bl2)
        min_image2  = r - sign(bl, r)
    elsewhere
        min_image2 = r
    end where
end function

pure subroutine find(l, i, ni)
        logical, intent(in) :: l(:)
        integer, intent(out) :: i(:), ni
        integer :: j

        ni = 0
        do j=1,size(l)
            if (l(j)) then
                ni = ni + 1
                i(ni) = j
            endif
        end do
end subroutine find

function allocate_aligned(nbytes) result(cptr)
    use iso_c_binding
    integer, intent(in) :: nbytes
    type(c_ptr) :: cptr
    call posix_memalign(cptr, 32_8, int(nbytes, c_size_t))
end function

end module utils
