module BeamElement
    !use mkl
    implicit none
    private
    integer, parameter :: ndofs = 12
    real(8) :: constParams(1:ndofs, 1:ndofs)
    public :: Segment
    type :: Segment
        private
        integer :: m_eID
        integer :: m_localToGlobal(1:ndofs)
        integer :: m_globalToLocal(1:ndofs)
        real(8) :: m_stifMat(1:ndofs, 1:ndofs)
        real(8) :: m_coefMat(1:ndofs, 1:ndofs)
    contains
        procedure :: init => Segment_init
        procedure :: Multiply => Segment_Multiply
        procedure :: UpdateMatrix => Segment_UpdateMatrix
    end type Segment

  contains
    subroutine Segment_init(this, eId, p0Id, p1Id, x)
        class(Segment), intent(inout) :: this
        real(8) :: x(1:ndofs)
        integer :: eID, p0Id, p1Id
        !symm(m_coefMat, x, b)
        this%m_localToGlobal(1) = ndofs * p0Id
    end subroutine Segment_init

    subroutine Segment_Multiply(this, x, b)
        class(Segment), intent(in) :: this
        real(8) :: x(1:ndofs), b(1:ndofs)
        !symm(m_coefMat, x, b)
        b = matmul(this%m_coefMat, x)
    end subroutine Segment_Multiply

    subroutine Segment_UpdateMatrix(this, x)
        class(Segment) :: this
        real(8) :: x(1:ndofs)
        integer :: i, j
        !update m_coefMat
        do i=1,ndofs
            do j=1,ndofs
                this%m_coefMat(i, j) = 0.
            enddo
            this%m_coefMat(i, i) = 1.
        enddo
    end subroutine Segment_UpdateMatrix
end module BeamElement

program test
    use BeamElement
    type(Segment) :: beam(2)
    real(8) :: x(1:12), b(1:12)
    do i = 1, 12
        x(i) = dble(i)
    enddo
    call beam(1)%UpdateMatrix(x)
    call beam(1)%Multiply(x, b)
    write(*, *) b
end