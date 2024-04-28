module BeamStrucutre
    !use mkl
    implicit none
    private
    integer, parameter :: nElmtDofs = 12
    real(8) :: coeffs(0:7)
    integer :: m_npts, m_nelmts, m_nmaterials, g_ndofs
    character(LEN=20) :: m_meshfile
    public :: Segment
    type :: Segment
        private
        integer :: m_localToGlobal(1:nElmtDofs)
        real(8) :: x0(1:nElmtDofs), len0
        real(8) :: m_property(1:8)
        real(8) :: m_coefMat(1:nElmtDofs, 1:nElmtDofs)
        real(8) :: m_tanMat(1:nElmtDofs, 1:nElmtDofs)
        real(8) :: m_stfMat(1:nElmtDofs, 1:nElmtDofs)
        real(8) :: m_masMat(1:nElmtDofs, 1:nElmtDofs)
        real(8) :: m_geoMat(1:nElmtDofs, 1:nElmtDofs)
    contains
        procedure :: init => Segment_init
        procedure :: Multiply => Segment_Multiply
        procedure :: UpdateMatrix => Segment_UpdateMatrix
        procedure :: LocToGlobal => Segment_LocToGlobal
        procedure :: GlobalToLoc => Segment_GlobalToLoc
        procedure :: FormMassMatrix => Segment_FormMassMatrix
        procedure :: FormStiffMatrix => Segment_FormStiffMatrix
        procedure :: FormGeomMatrix => Segment_FormGeomMatrix
        procedure :: FormDampMatrix => Segment_FormDampMatrix
        procedure :: RotateMatrix => Segment_RotateMatrix
    end type Segment
    type(Segment), allocatable :: m_elements(:)

  contains
    subroutine Beam_initialise(filename, gamma, beta, dt)
        implicit none
        character (LEN=20):: filename
        real(8) :: gamma, beta, dt
        integer :: fileiD = 996
        character (LEN=1000):: buffer
        real(8), allocatable :: xyz(:, :), material(:, :)
        ! coefficients in the Newmark-beta method
        coeffs(2) = 1. / (beta * dt)
        coeffs(1) = gamma * coeffs(2)
        coeffs(0) = coeffs(2) / dt
        coeffs(3) = 0.5 / beta - 1.
        coeffs(4) = gamma / beta - 1.
        coeffs(5) = dt * (0.5 * gamma / beta - 1.)
        coeffs(6) = dt * (1. - gamma)
        coeffs(7) = dt * gamma
        ! load mesh information
        m_meshfile = filename
        open(unit=fileiD, file = filename, status = 'old')
            read(fileiD,*) buffer
            !if(buffer(1:3).eq.'00')
            read(buffer, *) m_npts, m_nelmts, m_nmaterials
        close(fileiD)
        g_ndofs = m_npts * 6
        ! load points data
        allocate(xyz(1:3, 1:m_npts))
        call Beam_ReadPoints(xyz)
        ! load matieral data
        allocate(material(1:8, 1:m_nmaterials))
        call Beam_ReadMaterials(material)
        ! load and build elements 
        allocate(m_elements(1:m_nelmts))
        call Beam_ReadBuildElements(xyz, material)
    end subroutine Beam_initialise

    subroutine Beam_UpdateMatrix(disp)
        implicit none
        real(8), intent(inout) :: disp(1:6, 1:m_npts)
        integer :: i
        do i = 1, m_nelmts
            call m_elements(i)%UpdateMatrix(disp)
        enddo
    end subroutine

    subroutine Beam_Solve(disp, vel, acc)
        implicit none
        real(8), intent(inout) :: disp(1:6, 1:m_npts), vel(1:6, 1:m_npts), acc(1:6, 1:m_npts)
        real(8) :: dispn(1:6, 1:m_npts), rhs(1:6, 1:g_ndofs)
        integer :: maxNewtonRaphson, maxCG, i, j, n
        logical :: testconvergence
        ! solve the next dispalce, velocity and acceleration using CG method
        do i = 1, maxNewtonRaphson
            call Beam_UpdateMatrix(disp)
            call Beam_UpdateRHS(rhs)
            do j = 1, maxCG
                call CG_Solve(dispn, rhs)
                if(testconvergence) exit
            enddo
            if(testconvergence) exit
        enddo
        call Beam_UpdateVelAcc(dispn, disp, vel, acc)
    end subroutine

    subroutine CG_Solve(dispn, rhs)
        implicit none
        real(8) :: dispn(1:6, 1:m_npts), rhs(1:6, 1:g_ndofs), b(1:6, 1:g_ndofs) ! todo move local arrary to module variable
        !call Beam_MatrixMultipy(dispn, b)
    end subroutine

    subroutine Beam_UpdateVelAcc(dispn, disp, vel, acc)
        implicit none
        real(8), intent(inout) :: dispn(1:6, 1:m_npts), disp(1:6, 1:m_npts), vel(1:6, 1:m_npts), acc(1:6, 1:m_npts)
    end subroutine

    subroutine Beam_ReadMaterials(material)
        implicit none
        real(8) :: material(1:8, 1:m_nmaterials)
        ! to do, fill material
    end subroutine Beam_ReadMaterials

    subroutine Beam_readPoints(xyz)
        implicit none
        integer :: npts
        real(8) :: xyz(1:3, 1:m_npts)
        integer :: fileiD = 996, tmpid
        real(8) :: tmpx(1:3)

        open(unit=fileiD, file = m_meshfile, status = 'old')
            read(fileiD,*)
            read(fileiD,*)tmpid,tmpx(1),tmpx(2),tmpx(3)
            read(fileiD,*)
        close(fileiD)
    end subroutine Beam_ReadPoints
    subroutine Beam_ReadBuildElements(xyz, material)
        implicit none
        integer :: npts
        real(8) :: xyz(1:3, 1:m_npts), material(1:8, m_nmaterials)
        integer :: fileiD = 996, tmpid, n, i, j, k, im, itype
        integer(8) :: points(1:3)
        character(LEN=1000) :: buffer

        open(unit=fileiD, file = m_meshfile, status = 'old')
            do i=1,10000
                read(fileiD,*) buffer
                buffer = trim(buffer)
                if(buffer(1:3).eq. 'END') exit
            enddo
            do n=1,m_nelmts
                read(fileiD,*) tmpid, i, j, k, itype, im
                if(1.le.tmpid .and. tmpid.le.m_nelmts) then
                    call m_elements(tmpid)%init(i, j, xyz, material(1:8, im))
                endif
            enddo
        close(fileiD)
    end subroutine Beam_ReadBuildElements

    subroutine  Beam_UpdateRHS(b)
        implicit none
        real(8) :: b(1:g_ndofs)

    end subroutine Beam_UpdateRHS

    subroutine Beam_MatrixMultipy(x, b)
        implicit none
        real(8) :: x(1:g_ndofs), b(1:g_ndofs)
        integer :: i
        do i=1,m_nelmts
            call m_elements(i)%Multiply(x, b)
        enddo
    end subroutine Beam_MatrixMultipy

    subroutine Segment_init(this, p0Id, p1Id, xyz, material)
        class(Segment), intent(inout) :: this
        real(8), intent(in) :: xyz(1:3, 1:m_npts), material(1:8)
        integer :: p0Id, p1Id, i, offset0, offset1
        ! material property
        this%m_property(1:8) = material(1:8)
        !local dof to global dof mapping
        offset0 = (p0Id - 1) * 6
        offset1 = (p1Id - 1) * 6
        do i=1,6
            this%m_localToGlobal(i) = offset0 + i
            this%m_localToGlobal(i+6) = offset1 + i
        enddo
        this%x0(1:3) = xyz(1:3, p0Id)
        this%x0(7:9) = xyz(1:3, p1Id)
        do i=1,3
            !initialise angle dofs
        enddo
    end subroutine Segment_init

    subroutine Segment_GlobalToLoc(this, x, lx)
        class(Segment), intent(in) :: this
        real(8) :: x(1:g_ndofs)
        real(8) :: lx(1:nElmtDofs)
        integer :: i
        do i=1,nElmtDofs
            lx(i) = x(this%m_localToGlobal(i))
        enddo
    end subroutine Segment_GlobalToLoc

    subroutine Segment_LocToGlobal(this, lx, x)
        class(Segment), intent(in) :: this
        real(8) :: x(1:g_ndofs)
        real(8) :: lx(1:nElmtDofs)
        integer :: i
        do i=1,nElmtDofs
            x(this%m_localToGlobal(i)) = x(this%m_localToGlobal(i)) + lx(i)
        enddo
    end subroutine Segment_LocToGlobal

    subroutine Segment_Multiply(this, x, b)
        class(Segment), intent(in) :: this
        real(8) :: x(1:g_ndofs), b(1:g_ndofs)
        real(8) :: lx(1:nElmtDofs), lb(1:nElmtDofs)
        integer :: isdamp
        !symm(m_coefMat, x, b)
        call this%GlobalToLoc(x, lx)
        lb = matmul(this%m_coefMat, lx)
        call this%LocToGlobal(lb, b)
    end subroutine Segment_Multiply

    subroutine Segment_UpdateMatrix(this, x)
        class(Segment) :: this
        real(8) :: x(1:nElmtDofs), rotmat(1:3,1:3)
        integer :: i, j
        !update m_coefMat
        call this%FormMassMatrix
        call this%FormStiffMatrix
        call this%FormGeomMatrix
        call this%FormDampMatrix
        this%m_coefMat = this%m_stfMat + coeffs(0) * this%m_masMat ! to do consider the damping matrix
        call this%RotateMatrix()
    end subroutine Segment_UpdateMatrix

    subroutine Segment_FormMassMatrix(this)
        implicit none
        class(Segment), intent(in) :: this
    end subroutine Segment_FormMassMatrix

    subroutine Segment_FormStiffMatrix(this)
        implicit none
        class(Segment), intent(in) :: this
    end subroutine Segment_FormStiffMatrix

    subroutine Segment_FormGeomMatrix(this)
        implicit none
        class(Segment), intent(in) :: this
    end subroutine Segment_FormGeomMatrix

    subroutine Segment_FormDampMatrix(this)
        implicit none
        class(Segment), intent(in) :: this
    end subroutine Segment_FormDampMatrix

    subroutine Segment_RotateMatrix(this)
        implicit none
        class(Segment), intent(in) :: this
    end subroutine Segment_RotateMatrix
end module BeamStrucutre

program test
    use BeamStrucutre
    implicit none
    character*1 ::     trans
    double precision :: gamma = 0.75, beta = 0.25
    double precision :: A(12, 12), x(12), y(12)
    integer(8) :: i, j, m=12, n=12, lda=12, icx=1, icy=1
    external         DGEMV
    do i = 1, 12
        do j = 1, 12
            A(i, j) = 0.
        enddo
        A(i, i) = 1.
        x(i) = dble(i)
        y(i) = 0.
    enddo
    trans(1:1) = "T"
    write(*, *)x, trans
    !call gemv(M, x, y, 1., 0., 'N')
    call DGEMV(trans, m, m, gamma, A, lda, x, icx, beta, y, icy)
    !call dsymv('L', 12, 1., M, 12, x, 1, 0., y, 1)
    write(*, *)y

    !call Beam_initialise(filename, gamma, beta, dt)
    !scall Beam_Solve(disp, vel, acc)
end