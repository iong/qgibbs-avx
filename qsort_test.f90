program qsort_test
    use utils
    use iso_c_binding
    integer :: i
    integer(c_size_t) :: nel, width
    integer(c_int), allocatable :: idx(:)
    real(c_double), allocatable :: z_(:)

    nel = 1024
    width = c_int

    allocate(idx(nel), z_(nel))

    forall (i=1:nel)
        idx(i) = i
        z_(i) = mod(i-1, nel/2)
    end forall
    !call index_sort(nel, idx, z_)

    print '(I6, F12.6)', (idx(i), z_(idx(i)), i=1,nel)
contains

end program
