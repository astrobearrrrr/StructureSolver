module BeamStrucutre
    !use mkl
    implicit none
    private
    integer, parameter :: nElmtDofs = 12
    real(8) :: coeffs(0:7),dampM,dampK,gamma,alphaf
    integer :: m_npts, m_nelmts, m_nmaterials, g_ndofs
    character(LEN=20) :: m_meshfile
    public :: Segment, Beam_initialise, Beam_Solve
    type :: Segment
        private
        integer :: m_localToGlobal(1:nElmtDofs)
        real(8) :: x0(1:nElmtDofs),x1(1:nElmtDofs)
        real(8) :: dx0, dy0, dz0, dx1, dy1, dz1, xll0, xmm0, xnn0, xll1, xmm1, xnn1, len0, len1
        real(8) :: bc(1:nElmtDofs), geoFRM,stf(1:nElmtDofs, 1:nElmtDofs)
        real(8) :: triad_ee(3,3),triad_n1(3,3),triad_n2(3,3)
        real(8) :: m_property(1:8)
        real(8) :: m_coefMat(1:nElmtDofs, 1:nElmtDofs)
        real(8) :: m_tanMat(1:nElmtDofs, 1:nElmtDofs)
        real(8) :: m_stfMat(1:nElmtDofs, 1:nElmtDofs)
        real(8) :: m_masMat(1:nElmtDofs, 1:nElmtDofs)
        real(8) :: m_geoMat(1:nElmtDofs, 1:nElmtDofs)
        real(8) :: m_rotMat(1:3, 1:3)
    contains
        procedure :: init => Segment_init
        procedure :: Multiply => Segment_Multiply
        procedure :: UpdateMatrix => Segment_UpdateMatrix
        procedure :: UpdateLoad => Segment_UpdateLoad
        procedure :: LocToGlobal => Segment_LocToGlobal
        procedure :: GlobalToLoc => Segment_GlobalToLoc
        procedure :: FormMassMatrix => Segment_FormMassMatrix
        procedure :: FormStiffMatrix => Segment_FormStiffMatrix
        procedure :: FormGeomMatrix => Segment_FormGeomMatrix
        procedure :: RotateMatrix => Segment_RotateMatrix
        procedure :: RKR => Segment_RKR
        procedure :: BoundaryCond => Segment_BoundaryCond
        procedure :: BodyStress => Segment_BodyStress_D
        procedure :: InitTriad_D => Segment_InitTriad_D
        procedure :: UpdateTriad_D => Segment_UpdateTriad_D
        procedure :: MakeTriad_ee => Segment_MakeTriad_ee
    end type Segment
    type(Segment), allocatable :: m_elements(:)
    real(8), allocatable ::lodInte(:),lodExte(:),lodEffe(:),vBC(:)

  contains
    subroutine Beam_initialise(filename, Newmarkgamma, Newmarkbeta, dt, dampM1, dampK1, gamma1, alphaf1)
        implicit none
        character (LEN=20):: filename
        real(8) :: Newmarkgamma, Newmarkbeta, dt, dampM1, dampK1, gamma1, alphaf1
        integer :: fileiD = 996
        character (LEN=1000):: buffer
        real(8), allocatable :: xyz(:, :), material(:, :), boundary(:, :)
        ! coefficients in the Newmark-beta method
        coeffs(2) = 1. / (Newmarkbeta * dt)
        coeffs(1) = Newmarkgamma * coeffs(2)
        coeffs(0) = coeffs(2) / dt
        coeffs(3) = 0.5 / Newmarkbeta - 1.
        coeffs(4) = Newmarkgamma / Newmarkbeta - 1.
        coeffs(5) = dt * (0.5 * Newmarkgamma / Newmarkbeta - 1.)
        coeffs(6) = dt * (1. - Newmarkgamma)
        coeffs(7) = dt * Newmarkgamma
        dampM = dampM1
        dampK = dampK1
        gamma = gamma1
        alphaf = alphaf1
        ! load mesh information
        m_meshfile = filename
        open(unit=fileiD, file = filename )
            read(fileiD,*) buffer
            read(fileiD,*) m_npts, m_nelmts, m_nmaterials
        close(fileiD)
        g_ndofs = m_npts * 6
        allocate(vBC(1:g_ndofs))
        vBC(1:g_ndofs) = 0.0d0
        ! load points data
        allocate(xyz(1:3, 1:m_npts))
        call Beam_ReadPoints(xyz)
        ! load matieral data
        allocate(material(1:8, 1:m_nmaterials))
        call Beam_ReadMaterials(material)
        ! load boundary condition
        allocate(boundary(1:6, 1:m_npts))
        call Beam_ReadBoundary(boundary)
        ! load and build elements 
        allocate(m_elements(1:m_nelmts))
        call Beam_ReadBuildElements(xyz, material, boundary)
    end subroutine Beam_initialise
    
    subroutine Beam_ReadPoints(xyz)
        implicit none
        real(8) :: xyz(1:3, 1:m_npts)
        integer :: fileiD = 996, tmpid, i
        character(LEN=1000) :: buffer

        open(unit=fileiD, file = m_meshfile )
            do i=1,10000
                read(fileiD,*) buffer
                buffer = trim(buffer)
                if(buffer(1:5) .eq. 'POINT') exit
            enddo
            read(fileiD,*) m_npts
            do i = 1,m_npts
                read(fileiD,*)tmpid,xyz(1,i),xyz(2,i),xyz(3,i)
            enddo
        close(fileiD)
    end subroutine Beam_ReadPoints

    subroutine Beam_ReadMaterials(material)
        implicit none
        real(8) :: material(1:8, 1:m_nmaterials)
        integer :: fileiD = 996, tmpid, i
        character(LEN=1000) :: buffer

        open(unit=fileiD, file = m_meshfile )
            do i=1,10000
                read(fileiD,*) buffer
                buffer = trim(buffer)
                if(buffer(1:8) .eq. 'MATERIAL') exit
            enddo
            read(fileiD,*) m_nmaterials
            do i = 1,m_nmaterials
                read(fileiD,*) tmpid, material(1:8,i)
            enddo
        close(fileiD)
    end subroutine Beam_ReadMaterials
    
    subroutine Beam_ReadBoundary(boundary)
        implicit none
        real(8) :: boundary(1:6, 1:m_npts)
        integer :: fileiD = 996, tmpid, i
        character(LEN=1000) :: buffer

        open(unit=fileiD, file = m_meshfile )
            do i=1,10000
                read(fileiD,*) buffer
                buffer = trim(buffer)
                if(buffer(1:9) .eq. 'CONSTRAIN') exit
            enddo
            read(fileiD,*) m_npts
            do i = 1,m_npts
                read(fileiD,*) tmpid, boundary(1:6,i)
            enddo
        close(fileiD)
    end subroutine Beam_ReadBoundary

    subroutine Beam_ReadBuildElements(xyz, material, boundary)
        implicit none
        real(8) :: xyz(1:3, 1:m_npts), material(1:8, m_nmaterials), boundary(1:6, 1:m_npts)
        integer :: fileiD = 996, tmpid, n, i, j, k, imat, itype
        character(LEN=1000) :: buffer

        open(unit=fileiD, file = m_meshfile )
            do n=1,10000
                read(fileiD,*) buffer
                buffer = trim(buffer)
                if(buffer(1:7).eq. 'ELEMENT') exit
            enddo
            read(fileiD,*) m_nelmts
            do n=1,m_nelmts
                read(fileiD,*) tmpid, i, j, k, itype, imat
                if(1.le.tmpid .and. tmpid.le.m_nelmts) then
                    call m_elements(tmpid)%init(i, j, xyz, material(1:8, imat), boundary)
                endif
            enddo
        close(fileiD)
    end subroutine Beam_ReadBuildElements

    subroutine Beam_InitLoad
        implicit none
        allocate(lodInte(1:g_ndofs),lodExte(1:g_ndofs),lodEffe(1:g_ndofs))
        lodInte(1:g_ndofs) = 0.0d0
        call Beam_ReadExternalLoad
        lodEffe(1:g_ndofs) = 0.0d0
    end subroutine Beam_InitLoad
    subroutine Beam_ReadExternalLoad
        implicit none
        real(8) :: load(1:6, 1:m_npts)
        integer :: fileiD = 996, tmpid, i
        character(LEN=1000) :: buffer

        open(unit=fileiD, file = m_meshfile )
            do i=1,10000
                read(fileiD,*) buffer
                buffer = trim(buffer)
                if(buffer(1:12) .eq. 'EXTERNALLOAD') exit
            enddo
            read(fileiD,*) m_npts
            do i = 1,m_npts
                read(fileiD,*) tmpid, load(1:6,i)
            enddo
            do i = 1,m_npts
                lodExte((i-1)*6+1:(i-1)*6+6) = load(1:6,i)
            enddo
        close(fileiD)
    end subroutine Beam_ReadExternalLoad


    subroutine Beam_Solve(maxDynamic, maxNewtonRaphson, dtol)
        implicit none
        real(8) :: dspO(1:6, 1:m_npts), velO(1:6, 1:m_npts), accO(1:6, 1:m_npts)
        real(8) :: dsp(1:6, 1:m_npts), vel(1:6, 1:m_npts), acc(1:6, 1:m_npts)
        real(8) :: dspn(1:g_ndofs)
        integer :: maxDynamic, maxNewtonRaphson, i, j
        real(8) :: dnorm, dtol
        ! solve the next dispalce, velocity and acceleration using CG method
        call Beam_InitTriadANDFormMass
        call Beam_InitLoad
        call Beam_InitDspVelAcc(dsp, vel, acc)
        do i = 1, maxDynamic
            call Beam_InitDspVelAccATTimeT(dspO, velO, accO, dsp, vel, acc)
            dnorm=1.0d0
            do j = 1, maxNewtonRaphson
                call Beam_UpdateMatrixANDLoad(j,dspO,dsp,vel,acc)
                call CG_Solve(dspn, lodEffe)
                call Beam_UpdateDspANDTride(j, dspn, dsp, dnorm)
                if(dnorm .le. dtol) exit
            enddo
            call Beam_UpdateVelAcc(dspO, velO, accO, dsp, vel, acc)
        enddo
        return
    end subroutine Beam_Solve

    subroutine Beam_InitDspVelAcc(dsp, vel, acc)
        implicit none
        real(8), intent(inout) :: dsp(1:6, 1:m_npts), vel(1:6, 1:m_npts), acc(1:6, 1:m_npts)
        dsp(1:6, 1:m_npts) = 0.0d0
        vel(1:6, 1:m_npts) = 0.0d0
        acc(1:6, 1:m_npts) = 0.0d0
    end subroutine Beam_InitDspVelAcc
    subroutine Beam_InitDspVelAccATTimeT(dspO, velO, accO, dsp, vel, acc)
        implicit none
        real(8), intent(inout) :: dspO(1:6, 1:m_npts), velO(1:6, 1:m_npts), accO(1:6, 1:m_npts)
        real(8), intent(inout) :: dsp(1:6, 1:m_npts), vel(1:6, 1:m_npts), acc(1:6, 1:m_npts)
        dspO(1:6, 1:m_npts) = dsp(1:6, 1:m_npts)
        velO(1:6, 1:m_npts) = vel(1:6, 1:m_npts)
        accO(1:6, 1:m_npts) = acc(1:6, 1:m_npts)
    end subroutine Beam_InitDspVelAccATTimeT

    subroutine CG_Solve(x, b)
        implicit none
        real(8) :: x(1:g_ndofs), b(1:g_ndofs), r(1:g_ndofs), p(1:g_ndofs), Ap(1:g_ndofs)
        integer :: iter, max_iter
        double precision :: alpha, beta, rsold, rsnew, err
        ! maximum iteration count.
        max_iter = 10000
        err = 1d-10
        ! initialize the solution vector x to 0.
        x(1:g_ndofs) = 0.0d0
        ! initialize the residual vector. r = b - matmul(A, x)
        call Beam_MatrixMultipy(x, Ap)
        r = b-Ap
        p = r

        rsold = dot_product(r, r)
        do iter = 1, max_iter
            ! matrix free method, Calculate the matrix-vector product. Ap = matmul(A, p)
            call Beam_MatrixMultipy(p, Ap)
            alpha = rsold / dot_product(p, Ap)
            x = x + alpha * p
            r = r - alpha * Ap
            rsnew = dot_product(r, r)
            if(sqrt(rsnew) < err) exit
            beta = rsnew/rsold
            p = r + beta * p
            rsold = rsnew
        enddo
        return
    end subroutine CG_Solve

    subroutine Beam_UpdateDspANDTride(iter, dspn, dsp, dnorm)
        implicit none
        real(8) :: dspn(1:g_ndofs), dspnn(1:6, 1:m_npts)
        real(8), intent(inout) :: dsp(1:6, 1:m_npts)
        real(8) :: beta0,beta,zi,z0,maxramp,dnorm
        integer :: iter,i,node0,node1
        beta0 = 1.0d0
        maxramp = 0
        if    (iter <= maxramp) then
            zi=2.0d0**(iter)
            z0=2.0d0**(maxramp)
            beta=zi/z0*beta0
        else
            beta=1.0d0*beta0
        endif
        do i = 1, m_npts
            dspnn(1:6,i)= beta*dspn((i-1)*6+1:(i-1)*6+6)
        enddo
        dsp(1:6, 1:m_npts) = dsp(1:6, 1:m_npts) + dspnn(1:6, 1:m_npts)
        do i = 1, m_nelmts
            node0 = (m_elements(i)%m_localToGlobal(1)+5)/6
            node1 = (m_elements(i)%m_localToGlobal(7)+5)/6
            ! UpdateElementPosition
            m_elements(i)%x1(1:6) = m_elements(i)%x0(1:6) + dsp(1:6, node0)
            m_elements(i)%x1(7:12) = m_elements(i)%x0(7:12) + dsp(1:6, node1)
            m_elements(i)%dx1  = m_elements(i)%x1(7) - m_elements(i)%x1(1)
            m_elements(i)%dy1  = m_elements(i)%x1(8) - m_elements(i)%x1(2)
            m_elements(i)%dz1  = m_elements(i)%x1(9) - m_elements(i)%x1(3)
            m_elements(i)%len1 = dsqrt(m_elements(i)%dx1*m_elements(i)%dx1+m_elements(i)%dy1*m_elements(i)%dy1+m_elements(i)%dz1*m_elements(i)%dz1)
            m_elements(i)%xll1 = m_elements(i)%dx1/m_elements(i)%len1
            m_elements(i)%xmm1 = m_elements(i)%dy1/m_elements(i)%len1
            m_elements(i)%xnn1 = m_elements(i)%dz1/m_elements(i)%len1
            ! UpdateNodeTriad
            call m_elements(i)%UpdateTriad_D(dspnn)
            ! MakeElementTriad
            call m_elements(i)%MakeTriad_ee
        enddo
        dnorm=dabs(maxval(dspnn(1:6, 1:m_npts)**2))
        return
    end subroutine Beam_UpdateDspANDTride

    subroutine Beam_UpdateVelAcc(dspO, velO, accO, dsp, vel, acc)
        implicit none
        real(8), intent(inout) :: dspO(1:6, 1:m_npts), velO(1:6, 1:m_npts), accO(1:6, 1:m_npts)
        real(8), intent(inout) :: dsp(1:6, 1:m_npts), vel(1:6, 1:m_npts), acc(1:6, 1:m_npts)
            acc(1:6, 1:m_npts)  = coeffs(0)*(dsp(1:6, 1:m_npts) - dspO(1:6, 1:m_npts)) -coeffs(2)*velO(1:6, 1:m_npts) - coeffs(3)*accO(1:6, 1:m_npts)
            vel(1:6, 1:m_npts)  = velO(1:6, 1:m_npts) + coeffs(6)*accO(1:6, 1:m_npts) + coeffs(7)*acc(1:6, 1:m_npts)
        return
    end subroutine

    subroutine Beam_UpdateMatrixANDLoad(iter,dspO,dsp,vel,acc)
        implicit none
        real(8), intent(inout) :: dspO(1:6, 1:m_npts)
        real(8), intent(inout) :: dsp(1:6, 1:m_npts), vel(1:6, 1:m_npts), acc(1:6, 1:m_npts)
        integer :: i,iter
        
        lodInte = 0.0d0
        do i = 1, m_nelmts
            call m_elements(i)%UpdateMatrix
            call m_elements(i)%BodyStress
        enddo
        lodEffe = lodExte - lodInte
        do i = 1, m_nelmts
            call m_elements(i)%UpdateLoad(lodEffe,dspO,dsp,vel,acc)
            call m_elements(i)%BoundaryCond(iter)
        enddo
    end subroutine

    subroutine Segment_UpdateLoad(this,b,dspO,dsp,vel,acc)
        implicit none
        class(Segment), intent(inout) :: this
        real(8) :: b(1:g_ndofs)
        real(8) :: dspO(1:6, 1:m_npts)
        real(8) :: dsp(1:6, 1:m_npts), vel(1:6, 1:m_npts), acc(1:6, 1:m_npts)
        real(8) :: mss(1:nElmtDofs),wk1(1:nElmtDofs),wk2(1:nElmtDofs)
        integer :: i,node0,node1
        node0 = (this%m_localToGlobal(1)+5)/6
        node1 = (this%m_localToGlobal(7)+5)/6
        do i = 1,6
            mss(i) = (coeffs(0)*(dspO(i,node0)-dsp(i,node0))+coeffs(2)*vel(i,node0)+coeffs(3)*acc(i,node0))*this%m_masMat(i,i)   &
                                    +(coeffs(1)*(dspO(i,node0)-dsp(i,node0))+coeffs(4)*vel(i,node0)+coeffs(5)*acc(i,node0))*dampM*this%m_masMat(i,i)
            mss(i+6) = (coeffs(0)*(dspO(i,node1)-dsp(i,node1))+coeffs(2)*vel(i,node1)+coeffs(3)*acc(i,node1))*this%m_masMat(i+6,i+6)   &
                                    +(coeffs(1)*(dspO(i,node1)-dsp(i,node1))+coeffs(4)*vel(i,node1)+coeffs(5)*acc(i,node1))*dampM*this%m_masMat(i+6,i+6)
        enddo
        call this%LocToGlobal(mss, b)
        if (dampK .gt. 0.0d0) then
            do i = 1,6
                wk1(i)= coeffs(1)*(dspO(i,node0)-dsp(i,node0)) + coeffs(4)*vel(i,node0) +coeffs(5)*acc(i,node0)
                wk1(i+6)= coeffs(1)*(dspO(i,node1)-dsp(i,node1)) + coeffs(4)*vel(i,node1) +coeffs(5)*acc(i,node1)
            enddo
            wk2 = dampK*matmul(this%m_stfMat, wk1)
            call this%LocToGlobal(wk2, b)
        endif
        return
    end subroutine Segment_UpdateLoad


    subroutine Segment_init(this, p0Id, p1Id, xyz, material, boundary)
        class(Segment), intent(inout) :: this
        real(8), intent(in) :: xyz(1:3, 1:m_npts), material(1:8), boundary(1:6, 1:m_npts)
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
            ! initialise angle dofs
            this%x0(4:6) = 0.0d0
            this%x0(10:12) = 0.0d0
        enddo
        this%x1(1:12) = this%x0(1:12)
        this%bc(1:6) = boundary(1:6, p0Id)
        this%bc(7:12) = boundary(1:6, p1Id)
        this%dx0  = this%x0(7) - this%x0(1)
        this%dy0  = this%x0(8) - this%x0(2)
        this%dz0  = this%x0(9) - this%x0(3)
        this%len0 = dsqrt(this%dx0*this%dx0+this%dy0*this%dy0+this%dz0*this%dz0)
        this%xll0 = this%dx0/this%len0
        this%xmm0 = this%dy0/this%len0
        this%xnn0 = this%dz0/this%len0
        this%dx1  = this%x1(7) - this%x1(1)
        this%dy1  = this%x1(8) - this%x1(2)
        this%dz1  = this%x1(9) - this%x1(3)
        this%len1 = dsqrt(this%dx1*this%dx1+this%dy1*this%dy1+this%dz1*this%dz1)
        this%xll1 = this%dx1/this%len1
        this%xmm1 = this%dy1/this%len1
        this%xnn1 = this%dz1/this%len1
    end subroutine Segment_init

    subroutine Beam_MatrixMultipy(x, b)
        implicit none
        real(8) :: x(1:g_ndofs), b(1:g_ndofs)
        integer :: i
        b(1:g_ndofs) = 0.0d0
        do i=1,m_nelmts
            call m_elements(i)%Multiply(x, b)
        enddo
    end subroutine Beam_MatrixMultipy

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
        call this%GlobalToLoc(x, lx)
        lb = matmul(this%m_coefMat, lx)
        call this%LocToGlobal(lb, b)
    end subroutine Segment_Multiply

    subroutine Segment_UpdateMatrix(this)
        class(Segment), intent(inout) :: this
        !update m_coefMat
        call this%FormStiffMatrix
        call this%FormGeomMatrix
        this%stf = this%m_stfMat
        call this%RotateMatrix(this%xll1,this%xmm1,this%xnn1)

        call this%RKR(this%m_stfMat)

        call this%RKR(this%m_geoMat)

        this%m_tanMat = this%m_stfMat + gamma * this%m_geoMat
        this%m_coefMat = this%m_tanMat + coeffs(0) * this%m_masMat + coeffs(1) * dampM * this%m_masMat
        if(dampK .gt. 0.0d0) then
            this%m_coefMat = this%m_coefMat + coeffs(1) * dampK * this%m_stfMat
        endif
        return
    end subroutine Segment_UpdateMatrix

    subroutine Segment_FormMassMatrix(this)
        implicit none
        class(Segment), intent(inout) :: this
        real(8):: area,rho,zix,length
        real(8):: em(12,12)
        real(8):: roal

        area = this%m_property(3)
        rho  = this%m_property(4)
        zix  = this%m_property(6)
        length  = this%len0

        em(1:12,1:12) = 0.0d0
        roal = rho*area*length/2.0d0
        em(1,1)     = roal
        em(2,2)     = roal
        em(3,3)     = roal
        em(4,4)     = roal*zix/area
        em(5,5)     = roal*length*length*alphaf
        em(6,6)     = roal*length*length*alphaf
        em(7,7)     = em(1,1)
        em(8,8)     = em(2,2)
        em(9,9)     = em(3,3)
        em(10,10)   = em(4,4)
        em(11,11)   = em(5,5)
        em(12,12)   = em(6,6)

        this%m_masMat = 0.0d0
        this%m_masMat = em
        return
    end subroutine Segment_FormMassMatrix

    subroutine Segment_FormStiffMatrix(this)
        implicit none
        class(Segment), intent(inout) :: this
        real(8):: emod,gmod,area,zix,ziy,ziz,length
        real(8):: Invlength
        real(8):: ek(12,12)
        integer:: i,j
        real(8):: emlen,emlen2,emlen3
        
        emod = this%m_property(1)
        gmod = this%m_property(2)
        area = this%m_property(3)
        zix  = this%m_property(6)
        ziy  = this%m_property(7)
        ziz  = this%m_property(8)
        length  = this%len0

        ! initialize all ek elements to zero
        ek(1:12,1:12)=0.0
        
        ! STIFFNESS matrix in local coordinates
        ! emod is Modulus of elasticity
        Invlength = 1/length
        emlen  = emod*Invlength
        emlen2 = emlen*Invlength
        emlen3 = emlen2*Invlength

        ek(1,1)   =   area*emlen
        ek(2,2)   =   12.0*emlen3*ziz
        ek(3,3)   =   12.0*emlen3*ziy
        ek(4,4)   =   gmod*zix*Invlength
        ek(5,5)   =   4.0*emlen*ziy
        ek(6,6)   =   4.0*emlen*ziz
        
        ek(2,6)   =   6.0*emlen2*ziz
        ek(3,5)   =  -6.0*emlen2*ziy
        
        ek(7,7)   =   ek(1,1)
        ek(8,8)   =   ek(2,2)
        ek(9,9)   =   ek(3,3)
        ek(10,10) =   ek(4,4)
        ek(11,11) =   ek(5,5)
        ek(12,12) =   ek(6,6)
        
        ek(1,7)   =   -ek(1,1)
        ek(2,8)   =   -ek(2,2)
        ek(2,12)  =    ek(2,6)
        ek(3,9)   =   -ek(3,3)
        ek(3,11)  =    ek(3,5)
        ek(4,10)  =   -ek(4,4)
        ek(5,9)   =   -ek(3,5)
        ek(5,11)  =    ek(5,5)/2.0
        ek(6,8)   =   -ek(2,6)
        ek(6,12)  =    ek(6,6)/2.0
        
        ek(8,12)  =   -ek(2,6)
        ek(9,11)  =   -ek(3,5)
        
        ! impose the symmetry
        do  i= 1, 12
            do  j= i, 12
                ek(j,i) = ek(i,j)
            enddo
        enddo
        
        this%m_stfMat = 0.0d0
        this%m_stfMat = ek
        return
    end subroutine Segment_FormStiffMatrix

    subroutine Segment_FormGeomMatrix(this)
        implicit none
        class(Segment), intent(inout) :: this
        real(8):: length
        real(8):: eg(12,12),s!s is axial forces
        real(8):: emlenz,alpha
        integer:: i,j

        s = this%geoFRM

        length  = this%len0

        ! initialize all eg elements to zero
        eg(1:12,1:12)=0.0d0
        alpha   =   (s/length)*1.0d-6
        emlenz  =   s/(30.0*length)
        if (abs(alpha) .lt. 1.0e-10)then
            alpha = 1.0e-10
        endif
        ! Frame
        eg(1,1)   =   alpha
        eg(2,2)   =   36*emlenz 
        eg(3,3)   =   36*emlenz   
        eg(4,4)   =   alpha
        eg(5,5)   =   4.0*emlenz*length*length
        eg(6,6)   =   4.0*emlenz*length*length

        eg(2,6)   =   3.0*emlenz*length
        eg(3,5)   =  -3.0*emlenz*length

        ! Truss
        ! eg(1,1)   =   alpha
        ! eg(2,2)   =   emlenz
        ! eg(3,3)   =   emlenz
        ! eg(4,4)   =   alpha
        ! eg(5,5)   =   0.0d0
        ! eg(6,6)   =   0.0d0

        eg(7,7)   =   eg(1,1)
        eg(8,8)   =   eg(2,2)
        eg(9,9)   =   eg(3,3)
        eg(10,10) =   eg(4,4)
        eg(11,11) =   eg(5,5)
        eg(12,12) =   eg(6,6)

        eg(1,7)   =   -eg(1,1)
        eg(2,8)   =   -eg(2,2)
        eg(2,12)  =    eg(2,6)
        eg(3,9)   =   -eg(3,3)
        eg(3,11)  =    eg(3,5)
        eg(4,10)  =   -eg(4,4)
        eg(5,9)   =   -eg(3,5)
        eg(5,11)  =   -eg(5,5)/4.0
        eg(6,8)   =   -eg(2,6)
        eg(6,12)  =   -eg(6,6)/4.0

        eg(8,12)  =   -eg(2,6)
        eg(9,11)  =   -eg(3,5)

        ! impose the symmetry
        do  i= 1, 12
            do  j= i, 12
                eg(j,i) = eg(i,j)
            enddo
        enddo

        ! check diagonal terms
        do    i=1,12
            if (dabs(eg(i,i)) .lt. 1.0d-20) eg(i,i)=1.0d-20
        enddo

        this%m_geoMat = 0.0d0
        this%m_geoMat = eg
        return
    end subroutine Segment_FormGeomMatrix

    subroutine Segment_RotateMatrix(this,l,m,n)
        implicit none
        class(Segment), intent(inout) :: this
        real(8):: r(3,3)
        real(8):: l,m,n,beta,pi,sb,cb,d,Invd

        beta = this%m_property(5)

        pi=4.0*datan(1.0d0)
        !
        sb=dsin(beta*pi/180)
        cb=dcos(beta*pi/180)
        d=dsqrt(1.0-n**2)
        Invd=1/d
        ! if (abs(l).ge. 0.995 .and. abs(beta).le. 0.01) return
        if (abs(n).gt.0.995) then
            r(1,1)  =  0.0
            r(1,2)  =  0.0
            r(1,3)  =  n
            r(2,1)  = -n*sb
            r(2,2)  =  cb
            r(2,3)  =  0.0
            r(3,1)  = -n*cb
            r(3,2)  = -sb
            r(3,3)  =  0.0
        else
            r(1,1)  =  l
            r(1,2)  =  m
            r(1,3)  =  n
            if (abs(beta) .le. 0.01) then
                r(2,1)  =  -m*Invd
                r(2,2)  =  l*Invd
                r(2,3)  =  0.0
                r(3,1)  =  -l*n*Invd
                r(3,2)  =  -m*n*Invd
                r(3,3)  =  d
            else
                r(2,1)  =  -(m*cb+l*n*sb)*Invd
                r(2,2)  =  (l*cb-m*n*sb)*Invd
                r(2,3)  =  d*sb
                r(3,1)  =  (m*sb-l*n*cb)*Invd
                r(3,2)  =  -(l*sb+m*n*cb)*Invd
                r(3,3)  =  d*cb
            endif
        endif

        this%m_rotMat = r
        return
    end subroutine Segment_RotateMatrix

    subroutine Segment_RKR(this,ek)
        implicit none
        class(Segment), intent(inout) :: this
        real(8):: r(3,3),rt(3,3),ktemp(12,12),ek(12,12)
        integer:: i,j,k,j1,j2,ii,jj,in,jn

        r = this%m_rotMat
        do  in=1,3
        do  jn=1,3
            rt(jn,in)=this%m_rotMat(in,jn)
        enddo
        enddo
        ! take [Rtrans][K][R] using the nature of [R] to speed computation.
        ! k is sectioned off into 3x3s then multiplied [rtrans][k][r]
        !
        do  i=0,3
        do  j=0,3
            do    k=1,3
            do    ii=1,3
                j1=i*3
                j2=j*3
                ktemp(j1+k,j2+ii)=0.0
                do     jj=1,3
                ktemp(j1+k,j2+ii)=ktemp(j1+k,j2+ii)+ek(j1+k,j2+jj)*r(jj,ii)
                enddo
            enddo
            enddo
            do  k=1,3
            do  ii=1,3
                ek(j1+k,j2+ii)=0.0
                do  jj=1,3
                    ek(j1+k,j2+ii)=ek(j1+k,j2+ii)+rt(k,jj)*ktemp(j1+jj,j2+ii)
                enddo
            enddo
            enddo
        enddo
        enddo
        return
    end subroutine Segment_RKR

    subroutine Segment_BoundaryCond(this,iter)
        implicit none
        class(Segment), intent(inout) :: this
        integer :: i,iter
        do i=1,nElmtDofs
            if (this%bc(i).gt.0.0d0)then
                this%m_coefMat(i,i) = this%m_coefMat(i,i) * 1.0d20
                if(iter==0)then
                    lodEffe(this%m_localToGlobal(i))=this%m_coefMat(i,i)*vBC(this%m_localToGlobal(i)) &
                                                     + lodEffe(this%m_localToGlobal(i))
                else
                    lodEffe(this%m_localToGlobal(i))=0.0d0
                endif
            endif
        enddo
    end subroutine Segment_BoundaryCond

    subroutine Segment_InitTriad_D(this)
        implicit none
        class(Segment), intent(inout) :: this
        real(8):: xll,xmm,xnn
        real(8):: dd

        xll=this%xll0
        xmm=this%xmm0
        xnn=this%xnn0
        dd=dsqrt(xll*xll+xmm*xmm)
        this%triad_n1(1,1)=xll
        this%triad_n1(2,1)=xmm
        this%triad_n1(3,1)=xnn
        
        if    (dd .lt. 0.001d0) then
            this%triad_n1(1,2)=0.0d0
            this%triad_n1(2,2)=1.0d0
            this%triad_n1(3,2)=0.0d0
            this%triad_n1(1,3)=-xnn
            this%triad_n1(2,3)=0.00d0
            this%triad_n1(3,3)=0.00d0
        else
            this%triad_n1(1,2)=-xmm/dd
            this%triad_n1(2,2)=+xll/dd
            this%triad_n1(3,2)=0.0d0

            this%triad_n1(1,3)=-xll*xnn/dd
            this%triad_n1(2,3)=-xmm*xnn/dd                                    
            this%triad_n1(3,3)= dd 
        endif
        ! all element triads have same initial orientation
        this%triad_n2(1:3,1:3)=this%triad_n1(1:3,1:3)
        this%triad_ee(1:3,1:3)=this%triad_n1(1:3,1:3)
        return
    end subroutine Segment_InitTriad_D

    subroutine Segment_BodyStress_D(this)
        implicit none
        class(Segment), intent(inout) :: this
        real(8) :: triad_00(3,3),triad_11(3,3),triad_22(3,3)
        real(8) :: ub(12),dl
        real(8) :: rr(3,3)
        real(8) :: force(12),forceb(12)
        real(8) :: du,dv,dw
        real(8) :: tx1,tx2,ty1,ty2,tz1,tz2,tx,ty,tz
        real(8) :: fxx,emod,area
        integer :: i,j

        du = this%dx1 - this%dx0
        dv = this%dy1 - this%dy0
        dw = this%dz1 - this%dz0
        ! Define the the xyz directional distance of the two nodes of the beam element after deformation
        ! as well as the length, xl, of the element after deformation
        dl = ( (this%dx0+this%dx1)*du +(this%dy0+this%dy1)*dv +(this%dz0+this%dz1)*dw )/ (this%len0+this%len1)
        ! get twisting angles
        do    i=1,3
        do    j=1,3
            triad_00(i,j)=this%triad_ee(i,j)
            triad_11(i,j)=this%triad_ee(i,j)
            triad_22(i,j)=this%triad_n1(i,j)
        enddo
        enddo
        call Segment_get_angle_triad(triad_11,triad_22,tx,ty,tz)
        call Segment_global_to_local(triad_00,tx,ty,tz,tx1,ty1,tz1)
        !
        do    i=1,3
        do    j=1,3
            triad_11(i,j)=this%triad_ee(i,j)
            triad_22(i,j)=this%triad_n2(i,j)
        enddo
        enddo
        call Segment_get_angle_triad(triad_11,triad_22,tx,ty,tz)
        call Segment_global_to_local(triad_00,tx,ty,tz,tx2,ty2,tz2)

        ! non-zero ty1 tz1 u2 tx2 ty2 tz2
        ub(1)=0.0
        ub(2)=0.0
        ub(3)=0.0
        ub(4)=tx1
        ub(5)=ty1
        ub(6)=tz1
        !
        ub(7)=dl
        ub(8)=0.0
        ub(9)=0.0
        ub(10)=tx2
        ub(11)=ty2
        ub(12)=tz2
        !
        ! compute axial force
        emod = this%m_property(1)
        area = this%m_property(3)
        fxx=dl*emod*area/this%len0
        ! save local force for geo stiff
        ! write(igeoFRM,'(D25.15)') fxx
        this%geoFRM=fxx
        ! nodal forces in local coords
        ! {F}=[k]{u}
        forceb(1:12) = matmul(this%stf,ub)
        ! transform to global
        do  i=1,3
        do  j=1,3
            rr(i,j)=this%triad_ee(i,j)
        enddo
        enddo
        do  i=1,12
            force(i)=0.0d0
        enddo
        do    i=1,3
        do    j=1,3
            force(0+i) = force(0+i) + rr(i,j)*forceb(0+j)
            force(3+i) = force(3+i) + rr(i,j)*forceb(3+j)
            force(6+i) = force(6+i) + rr(i,j)*forceb(6+j)
            force(9+i) = force(9+i) + rr(i,j)*forceb(9+j)
        enddo
        enddo
        call this%LocToGlobal(force, lodInte)
        return
    end subroutine Segment_BodyStress_D

    subroutine Segment_get_angle_triad(triad_11,triad_22,tx,ty,tz)
        implicit none
        real(8):: triad_11(3,3),triad_22(3,3)
        real(8):: rr(3,3)
        real(8):: tx,ty,tz, dtx,dty,dtz,c1,tt,sint
        integer:: i,j,k
        !
        ! get angle between two triads
        do    i=1,3
            do    j=1,3
                rr(i,j)=0.0d0
                do    k=1,3
                    rr(i,j)=rr(i,j) + triad_22(i,k)*triad_11(j,k)
                enddo
            enddo
        enddo
    
        dtx = (rr(3,2)-rr(2,3))/2.0d0
        dty = (rr(1,3)-rr(3,1))/2.0d0
        dtz = (rr(2,1)-rr(1,2))/2.0d0
    
        c1=1.0d0
        sint = dsqrt(dtx*dtx+dty*dty+dtz*dtz)
    
        if (sint .gt. 1.0d0) sint=1.0d0
        tt = dasin(sint)
        if ( sint .lt. 1.0d-6) then
             c1=1.0d0
        else
             c1 = tt/sint
        endif
    
        tx=c1*dtx
        ty=c1*dty
        tz=c1*dtz
    
        return
    end subroutine Segment_get_angle_triad
    
    subroutine Segment_global_to_local(triad,tx,ty,tz,tx2,ty2,tz2)
        implicit none
        real(8):: triad(3,3)
        real(8):: tx,ty,tz,tx2,ty2,tz2
        tx2 = triad(1,1)*tx+triad(2,1)*ty+triad(3,1)*tz
        ty2 = triad(1,2)*tx+triad(2,2)*ty+triad(3,2)*tz
        tz2 = triad(1,3)*tx+triad(2,3)*ty+triad(3,3)*tz
        return
    end subroutine Segment_global_to_local

    subroutine Segment_UpdateTriad_D(this,dspnn)
        implicit none
        class(Segment), intent(inout) :: this
        real(8):: dspnn(1:6,1:m_npts)
        real(8):: rr(3,3)
        real(8):: dtx1,dty1,dtz1
        integer:: node0,node1

        node0 = (this%m_localToGlobal(1)+5)/6
        node1 = (this%m_localToGlobal(7)+5)/6
        ! n1 node
        dtx1=dspnn(4,node0)
        dty1=dspnn(5,node0)
        dtz1=dspnn(6,node0)
        call Segment_FiniteRot(dtx1,dty1,dtz1,rr)
        this%triad_n1=matmul(rr,this%triad_n1)
        ! n2 node
        dtx1=dspnn(4,node1)
        dty1=dspnn(5,node1)
        dtz1=dspnn(6,node1)
        call Segment_FiniteRot(dtx1,dty1,dtz1,rr)
        this%triad_n2=matmul(rr,this%triad_n2)

        return
    end subroutine Segment_UpdateTriad_D

    subroutine Segment_MakeTriad_ee(this)
        implicit none
        class(Segment), intent(inout) :: this
        real(8):: triad_aa(3,3)
        real(8):: rr(3,3),tx,ty,tz
        real(8):: triad_11(3,3),triad_22(3,3)
        real(8):: xll,xmm,xnn,dd,r2e1,r3e1
        integer:: i,j,k

        xll=this%xll1
        xmm=this%xmm1
        xnn=this%xnn1
        dd =dsqrt(xll*xll+xmm*xmm)
        do    i=1,3
            do    j=1,3
                this%triad_ee(i,j)=0.0d0
            enddo
        enddo
        this%triad_ee(1,1)=xll
        this%triad_ee(2,1)=xmm
        this%triad_ee(3,1)=xnn
        !
        ! get angle between two triads
        do    i=1,3
            do    j=1,3
                triad_11(i,j)=this%triad_n1(i,j)
                triad_22(i,j)=this%triad_n2(i,j)
            enddo
        enddo
        call Segment_get_angle_triad(triad_11,triad_22,tx,ty,tz)
        !
        ! rotate n1 to intermediate
        tx=tx/2.0d0
        ty=ty/2.0d0
        tz=tz/2.0d0
        call Segment_FiniteRot(tx,ty,tz,rr)
        triad_aa(1:3,1:3)=matmul(rr(1:3,1:3),this%triad_n1(1:3,1:3))
        !
        ! vectors e2 e3
        r2e1 = 0.0d0
        r3e1 = 0.0d0
        do    k=1,3
            r2e1 = r2e1 + triad_aa(k,2)*this%triad_ee(k,1)
            r3e1 = r3e1 + triad_aa(k,3)*this%triad_ee(k,1)
        enddo
        do    j=1,3
            this%triad_ee(j,2)=triad_aa(j,2) - r2e1*(triad_aa(j,1)+this%triad_ee(j,1))/2.0d0
            this%triad_ee(j,3)=triad_aa(j,3) - r3e1*(triad_aa(j,1)+this%triad_ee(j,1))/2.0d0
        enddo
        return
    end subroutine Segment_MakeTriad_ee
    
    subroutine Segment_FiniteRot(t1,t2,t3,rr)
        implicit none
        real(8):: rr(3,3),rr1(3,3),rr2(3,3),rr3(3,3)
        real(8):: t1,t2,t3,tt,ss,cc,c1,c2
        !
        integer:: i,j
        !
        tt=dsqrt( t1**2 + t2**2 + t3**2 )
        ss=dsin(tt)
        cc=dcos(tt)
        !
        rr1(1,1)=1.0d0
        rr1(1,2)=0.0d0
        rr1(1,3)=0.0d0
        rr1(2,1)=0.0d0
        rr1(2,2)=1.0d0
        rr1(2,3)=0.0d0
        rr1(3,1)=0.0d0
        rr1(3,2)=0.0d0
        rr1(3,3)=1.0d0
        !
        rr2(1,1)=0.0d0
        rr2(1,2)=-t3
        rr2(1,3)= t2
        rr2(2,1)= t3
        rr2(2,2)=0.0d0
        rr2(2,3)=-t1
        rr2(3,1)=-t2
        rr2(3,2)= t1
        rr2(3,3)=0d0
        !
        rr3(1,1)=-t3*t3-t2*t2
        rr3(1,2)= t2*t1
        rr3(1,3)= t3*t1
        rr3(2,1)= t1*t2
        rr3(2,2)=-t3*t3-t1*t1
        rr3(2,3)= t3*t2
        rr3(3,1)= t1*t3
        rr3(3,2)= t2*t3
        rr3(3,3)=-t2*t2-t1*t1
        !
        do    i=1,3
        do    j=1,3
            if    (tt .lt. 1.0d-10) then
                c1=1.0d0
                ! c2=1.0
                c2=0.5d0
            else
                c1 = ss/tt
                c2 = (1.0d0-cc)/tt**2
            endif
            rr(i,j) = rr1(i,j) + rr2(i,j)*c1 + rr3(i,j)*c2
        enddo
        enddo
        !
        return
    end subroutine Segment_FiniteRot

    subroutine Beam_InitTriadANDFormMass
        implicit none
        integer:: i
    
        do i = 1, m_nelmts
            ! InitTriad
            call m_elements(i)%InitTriad_D
            ! FormMass
            call m_elements(i)%FormMassMatrix
            call m_elements(i)%RotateMatrix(m_elements(i)%xll0,m_elements(i)%xmm0,m_elements(i)%xnn0)
            call m_elements(i)%RKR(m_elements(i)%m_masMat)
        enddo

        return
    end subroutine Beam_InitTriadANDFormMass
end module BeamStrucutre

program main
    use BeamStrucutre
    implicit none
    ! character*1 ::     trans
    ! real(8) :: alpha = 1., beta = 0.0
    ! real(8) :: A(12, 12), x(12), y(12)
    ! integer(8) :: i, j, m=12, lda=12, icx=1, icy=1
    ! external         DGEMV, dsymv
    ! do i = 1, 12
    !     do j = 1, 12
    !         A(i, j) = 0.
    !     enddo
    !     A(i, i) = 1.
    !     x(i) = dble(i)
    !     y(i) = 0.
    ! enddo
    ! trans(1:1) = "T"
    ! write(*, *)x, trans
    ! write(*, *) "======DGEMV======="
    ! call DGEMV(trans, m, m, alpha, A, lda, x, icx, beta, y, icy)
    ! write(*, *)y, trans
    ! write(*, *) "======dsymv======="
    ! trans(1:1) = "L"
    ! call dsymv(trans, m, alpha, A, lda, x, icx, beta, y, icy)
    ! write(*, *)y, trans
    
    ! test example: ISBN 9787040258417 Zeng Pan. P45
    ! character (LEN=20)::filename
    ! filename = 'Beam.dat'
    ! call Beam_initialise(filename, 5d-1, 1d0, 5d0, 0d0, 0d0, 1d0, 0d0)
    ! call Beam_Solve(1, 200, 1d-3)
    character (LEN=20)::filename
    filename = 'Dynamic.dat'
    call Beam_initialise(filename, 0.5d0, 0.25d0, 0.12d0, 0d0, 0d0, 1d0, 0d0)
    call Beam_Solve(10, 200, 1d-3)
end program