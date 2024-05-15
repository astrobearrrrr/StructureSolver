# StructureSolver
Structure Solver for structure
running programme
***********************************************************************
## How to run
```
gfortran -ffree-line-length-none -o StructureSolver StructureSolver.f90
.\StructureSolver.exe
```
```
gfortran -ffree-line-length-none -o BeamStructure BeamStrucutre.f90
.\BeamStructure.exe
```
***********************************************************************
## Book
- ISBN 9787040258417 有限元基础教程 曾攀
- ISBN 9781441929105 Nonlinear Analysis of Thin-Walled Structures: Statics, Dynamics, and Stability James F. Doyle
- ISBN 9787302388333 计算动力学 张雄
- ISBN 9780792312086 Static and Dynamic Analysis of Structures with An Emphasis on Mechanics and Computer Matrix Methods James F. Doyle
***********************************************************************
## test case
- static: ISBN 9787040258417 有限元基础教程 Zeng Pan 曾攀. P45
- dynamic: ISBN 9787576318555 计算固体力学 Dong Chunying 董春迎. P146 8.2
***********************************************************************
### static: ISBN 9787040258417 有限元基础教程 Zeng Pan 曾攀. P45
- Question
```
 y       25000N
 ^         |
 |         v
 |>--------|  n3
 n4       /|
        /  |  300mm
      /    |
    /      |
 |>--------|->20000N  ->x
    400mm  A
 n1        n2
 ```
node 1 and 3 is hinge joint, node 2 is sliding hinge which displace can happen in the x direction, node 3 is free.

$A = 100mm^{2}, E = 2.95*10^{5}N/mm^{2}$

- Results:
- displacement

  node 2: $d _{x}= 0.2712mm$

  node 3: $d _{x}= 0.0565mm$

  node 3: $d _{y}=-0.2225mm$

- Main program can be written as
```
program main
    use BeamStrucutre
    implicit none
    character (LEN=20)::filename
    filename = 'Beam.dat'
    call Beam_initialise(filename, 5d-1, 1d0, 5d0, 0d0, 0d0, 1d0, 0d0)
    call Beam_Solve(1, 200, 1d-3)
end program
```
***********************************************************************
### dynamic: ISBN 9787576318555 计算固体力学 Dong Chunying 董春迎. P146 8.2
- Question
```
    beam 1    beam 2
    0.5cm      0.5cm  
 |----------__________
 |          __________|---->F  ->x
 |----------
n1         n2         n3
```
$A_{1} = 2A, A_{2} = A, A = 2cm^{2}, E = 50 Pa, \rho = 5*10^{6}kg/m^{3}$

Initional: node 1-3 

Displacement: 
$\textbf{d} _{0} = 0$

Velocity: 
$\textbf{v} _{0} = 0$

Acceleration: 
$\textbf{a} _{0} = \textbf{M}^{-1}(\textbf{F} _{0}-\textbf{Cv} _{0}-\textbf{Kd} _{0})$

Force: 
$F=1$ at $t=0 \sim 0.5$, $F=0$ at $t \textgreater 0.5$ 

Two degrees of freedom: node 2 $x$-direction, node 3 $x$-direction
$$M = 1/4*[3 ,  0 ;  0 , 1]$$
$$K =   2*[3 , -1 ; -1 , 1]$$
$$C =     [0 ,  0 ;  0 , 0]$$

- Results of calculation using Newmark method.
 **MATLAB** code source from Spring course, **Computational Solid Mechanics, Chunying Dong**.
```
=========   Newmark method   =================
    Time point         Node 2       Node 3 
      0.1200         0.00000275   0.00287174
      0.2400         0.00002193   0.01145402
      0.3600         0.00009009   0.02564859
      0.4800         0.00026063   0.04529332
      0.6000         0.00060524   0.06872886
      0.7200         0.00120669   0.09281847
      0.8400         0.00214866   0.11585673
      0.9600         0.00350832   0.13759002
      1.0800         0.00535303   0.15778453
      1.2000         0.00773873   0.17622961
```

#### Before calculating the dynamic case, some changes needed to modify the BeamStructture.f90 file. 
1. change Beam_Solve
```
subroutine Beam_Solve(maxDynamic, maxNewtonRaphson, dtol)
...
    do i = 1, maxDynamic
        call Beam_InitDspVelAccATTimeT(dspO, velO, accO, dsp, vel, acc)
        dnorm=1.0d0
        do j = 1, maxNewtonRaphson
...
end subroutine Beam_Solve
```
change to
```
subroutine Beam_Solve(maxDynamic, maxNewtonRaphson, dtol)
...
do i = 1, maxDynamic
    call Beam_InitDspVelAccATTimeT(dspO, velO, accO, dsp, vel, acc)
    if(i>4)then
        lodExte=0.0d0
    endif
    dnorm=1.0d0
    do j = 1, maxNewtonRaphson
...
end subroutine Beam_Solve
```

2. change Beam_InitDspVelAcc
```
subroutine Beam_InitDspVelAcc(dsp, vel, acc)
    implicit none
    real(8), intent(inout) :: dsp(1:6, 1:m_npts), vel(1:6, 1:m_npts), acc(1:6, 1:m_npts)
    dsp(1:6, 1:m_npts) = 0.0d0
    vel(1:6, 1:m_npts) = 0.0d0
    acc(1:6, 1:m_npts) = 0.0d0
end subroutine Beam_InitDspVelAcc
```
change to
```
subroutine Beam_InitDspVelAcc(dsp, vel, acc)
    implicit none
    real(8), intent(inout) :: dsp(1:6, 1:m_npts), vel(1:6, 1:m_npts), acc(1:6, 1:m_npts)
    real(8) :: M(1:g_ndofs),temp
    integer :: i,j,node,dof
    dsp(1:6, 1:m_npts) = 0.0d0
    vel(1:6, 1:m_npts) = 0.0d0
    acc(1:6, 1:m_npts) = 0.0d0
    M=0.0d0
    do i=1,m_nelmts
        do j = 1,nElmtDofs
            M(m_elements(i)%m_localToGlobal(j))=m_elements(i)%m_masMat(j,j)+M(m_elements(i)%m_localToGlobal(j))
        enddo
    enddo
    do i=1,g_ndofs
        if(M(i).gt.1d-5)then
        M(i)=1/M(i)
        endif
    enddo
    do i = 1,g_ndofs
        node = (i+5)/6
        dof = i-(node-1)*6
        temp = M(i)*lodExte(i)
        acc(dof,node)=temp
    enddo
end subroutine Beam_InitDspVelAcc
```
3. change Beam_UpdateDspANDTride
```
subroutine Beam_UpdateDspANDTride(iter, dspn, dsp, dnorm)
...
        do i = 1, m_npts
            dspnn(1:6,i)= beta*dspn((i-1)*6+1:(i-1)*6+6)
        enddo
        dsp(1:6, 1:m_npts) = dsp(1:6, 1:m_npts) + dspnn(1:6, 1:m_npts)
...
```
change to
```
subroutine Beam_UpdateDspANDTride(iter, dspn, dsp, dnorm)
...
        do i = 1, m_npts
            dspnn(1:6,i)= beta*dspn((i-1)*6+1:(i-1)*6+6)
        enddo
        dsp(1:6, 1:m_npts) =  dspnn(1:6, 1:m_npts)
...
```
4. change Beam_UpdateMatrixANDLoad
```
subroutine Beam_UpdateMatrixANDLoad(iter,dspO,dsp,vel,acc)
...
    lodEffe = lodExte - lodInte
...
end subroutine
```
change to 
```
subroutine Beam_UpdateMatrixANDLoad(iter,dspO,dsp,vel,acc)
...
lodEffe = lodExte
...
end subroutine
```

5. change Segment_UpdateMatrix
```
subroutine Segment_UpdateMatrix(this)
...
    this%m_tanMat = this%m_stfMat + gamma * this%m_geoMat
...
end subroutine Segment_UpdateMatrix
```
change to 
```
...
    this%m_tanMat = this%m_stfMat
...
```

6. change Segment_UpdateLoad
```
subroutine Segment_UpdateLoad(this,b,dspO,dsp,vel,acc)
...
    do i = 1,6
        mss(i)   = (coeffs(0)*(dspO(i,node0)-dsp(i,node0))+coeffs(2)*vel(i,node0)+coeffs(3)*acc(i,node0))*this%m_masMat(i,i)   &
                    +(coeffs(1)*(dspO(i,node0)-dsp(i,node0))+coeffs(4)*vel(i,node0)+coeffs(5)*acc(i,node0))*dampM*this%m_masMat(i,i)
        mss(i+6) = (coeffs(0)*(dspO(i,node1)-dsp(i,node1))+coeffs(2)*vel(i,node1)+coeffs(3)*acc(i,node1))*this%m_masMat(i+6,i+6)   &
                    +(coeffs(1)*(dspO(i,node1)-dsp(i,node1))+coeffs(4)*vel(i,node1)+coeffs(5)*acc(i,node1))*dampM*this%m_masMat(i+6,i+6)
    enddo
    call this%LocToGlobal(mss, b)
    if (dampK .gt. 0.0d0) then
        do i = 1,6
            wk1(i)  = coeffs(1)*(dspO(i,node0)-dsp(i,node0)) + coeffs(4)*vel(i,node0) +coeffs(5)*acc(i,node0)
            wk1(i+6)= coeffs(1)*(dspO(i,node1)-dsp(i,node1)) + coeffs(4)*vel(i,node1) +coeffs(5)*acc(i,node1)
        enddo
        wk2 = dampK*matmul(this%m_stfMat, wk1)
        call this%LocToGlobal(wk2, b)
    endif
    return
end subroutine
```
change to
```
subroutine Segment_UpdateLoad(this,b,dspO,dsp,vel,acc)
...
    do i = 1,6
        mss(i)   = (coeffs(0)*(dsp(i,node0))+coeffs(2)*vel(i,node0)+coeffs(3)*acc(i,node0))*this%m_masMat(i,i)   &
                    +(coeffs(1)*(dsp(i,node0))+coeffs(4)*vel(i,node0)+coeffs(5)*acc(i,node0))*dampM*this%m_masMat(i,i)
        mss(i+6) = (coeffs(0)*(dsp(i,node1))+coeffs(2)*vel(i,node1)+coeffs(3)*acc(i,node1))*this%m_masMat(i+6,i+6)   &
                    +(coeffs(1)*(dsp(i,node1))+coeffs(4)*vel(i,node1)+coeffs(5)*acc(i,node1))*dampM*this%m_masMat(i+6,i+6)
    enddo
    call this%LocToGlobal(mss, b)
    if (dampK .gt. 0.0d0) then
        do i = 1,6
            wk1(i)  = coeffs(1)*(dsp(i,node0)) + coeffs(4)*vel(i,node0) +coeffs(5)*acc(i,node0)
            wk1(i+6)= coeffs(1)*(dsp(i,node1)) + coeffs(4)*vel(i,node1) +coeffs(5)*acc(i,node1)
        enddo
        wk2 = dampK*matmul(this%m_stfMat, wk1)
        call this%LocToGlobal(wk2, b)
    endif
    return
end subroutine
```
7. main program can be written as
```
program main
    use BeamStrucutre
    implicit none
    character (LEN=20)::filename
    filename = 'Dynamic.dat'
    call Beam_initialise(filename, 0.5d0, 0.25d0, 0.12d0, 0d0, 0d0, 1d0, 0d0)
    call Beam_Solve(10, 1, 1d-3)
end program
```
