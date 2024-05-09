! test example: ISBN 9787040258417 Zeng Pan. P45
! program algorithm: ISBN 9781441929105 James F. Doyle. P354
program main
    implicit none

    integer:: nND,nEL,nEQ,nMT
    integer, allocatable:: jBC(:,:),ele(:,:)
    real(8), allocatable:: vBC(:,:),xyzful0(:,:),xyzful(:,:),prop(:,:)
    real(8), allocatable:: mss(:),dspful(:,:),velful(:,:),accful(:,:)
    real(8):: deltat
    real(8):: dampK,dampM

    real(8), allocatable:: triad_nn(:,:,:),triad_ee(:,:,:)
    real(8), allocatable:: triad_n1(:,:,:),triad_n2(:,:,:)
!   ===============================================================================================
    real(8):: Newmarkdelta,Newmarkalpha
    real(8):: dtol
    integer:: iterMax
    integer:: i,iND
!   ===============================================================================================
    real(8), allocatable:: xyzful00(:,:),mssful(:,:)
    real(8), allocatable:: lodful(:,:),extful(:,:),grav(:,:)
    real(8), allocatable:: xyzfulnxt(:,:)
    real(8):: alphaf,alpham,alphap,g(3)
    integer:: step,numsubstep,isubstep
    real(8):: dt,time,timeSimTotl,subdeltat,Tref=1
    integer:: idFlow = 11, idat = 12
    
!   ===============================================================================================
    open(unit=idat, file = 'Beam.dat')
        rewind(idat)
        read(idat,*)
        read(idat,*)nND,nEL,nMT
        read(idat,*)
    close(idat)
!   ===============================================================================================
    allocate(jBC(nND,6),ele(nEL,5) )
    allocate(vBC(nND,6),xyzful0(nND,6),xyzful(nND,6),prop(nMT,10) )
    allocate(mss(nND*6),dspful(nND,6),velful(nND,6),accful(nND,6) )

    allocate(triad_nn(3,3,nND),triad_ee(3,3,nND) )
    allocate(triad_n1(3,3,nND),triad_n2(3,3,nND) )

    allocate(xyzful00(nND,6),mssful(nND,6) )
    allocate(lodful(nND,6),extful(nND,6),grav(nND,6) )
    allocate(xyzfulnxt(nND,6) )
!   ===============================================================================================
    open(unit=idFlow, file = 'inFlow.dat')
    do i = 1,5
        read(idFlow,*)
    enddo
    read(idFlow,*)     timeSimTotl
    close(idFlow)
    open(unit=idFlow, file = 'inFlow.dat')
    do i = 1,15
        read(idFlow,*)
    enddo
    read(idFlow,*)     dt
    close(idFlow)
    open(unit=idFlow, file = 'inFlow.dat')
    do i = 1,39
        read(idFlow,*)
    enddo
    read(idFlow,*)     numsubstep
    read(idFlow,*)     dampK,dampM
    read(idFlow,*)     Newmarkdelta,Newmarkalpha
    read(idFlow,*)     alphaf,alpham,alphap
    read(idFlow,*)     dtol,iterMax
    close(idFlow)
!   ===============================================================================================
    open(unit=idat, file = 'Beam.dat')
    rewind(idat)
    read(idat,*)
    read(idat,*)nND,nEL,nMT
    read(idat,*)
    call readdt(jBC(1:nND,1:6),ele(1:nEL,1:5),xyzful00(1:nND,1:6),prop(1:nMT,1:10), &
                nND,nEL,nEQ,nMT,idat)
    close(idat)
!   ===============================================================================================
    xyzful0(1:nND,1:6)=xyzful00(1:nND,1:6)

    xyzful(1:nND,1:6)=xyzful0(1:nND,1:6)     
    velful(1:nND,1:6)=0.0d0
                
    dspful(1:nND,1:6)=0.0d0
    accful(1:nND,1:6)=0.0d0
!   ===============================================================================================
    g(1:3)=0.0d0
    xyzfulnxt(1:nND,1:6)=xyzful(1:nND,1:6)
    extful(1:nND,1:6)=0.0d0
!   ===============================================================================================
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Assemble the mass matrices [M].
!   -------------------------------------------------------------------
    CALL formmass_D(  ele(1:nEL,1:5),xyzful0(1:nND,1),xyzful0(1:nND,2),xyzful0(1:nND,3), &
                     prop(1:nMT,1:10),mss(1:nND*6),nND,nEL,nEQ,nMT,alphaf )
    do  iND = 1, nND
        mssful(iND,1:6)=mss((iND-1)*6+1:(iND-1)*6+6)
        grav(iND,1:6) = mssful(iND,1:6)*[g(1),g(2),g(3),0.0d0,0.0d0,0.0d0]
    enddo
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Initialize all triads.
!   -------------------------------------------------------------------
    CALL init_triad_D( ele(1:nEL,1:5),xyzful0(1:nND,1),xyzful0(1:nND,2),xyzful0(1:nND,3), &
                      triad_nn(1:3,1:3,1:nND),triad_n1(1:3,1:3,1:nEL),triad_n2(1:3,1:3,1:nEL), &
                      triad_ee(1:3,1:3,1:nEL),nND,nEL ) 
!   ===============================================================================================
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Begin loop over time increments:
!   -------------------------------------------------------------------
    deltat = dt
    time=0.0d0
    step=0
    do while(time/Tref .le. timeSimTotl)
        time=time+dt
        step=step+1
        !using substep for bodies
        subdeltat=deltat/dble(numsubstep)
        do isubstep=1,numsubstep
            !displacement condition
            vBC(1:nND,1:6) = xyzfulnxt(1:nND,1:6)-xyzful(1:nND,1:6)
            !loading vector
            lodful(1:nND,1:6) = extful(1:nND,1:6)+grav(1:nND,1:6)
            !-------------------------------------------------------------------
            CALL StructureSolver(jBC(1:nND,1:6),vBC(1:nND,1:6),ele(1:nEL,1:5), &
                                prop(1:nMT,1:10),mss(1:nND*6), &
                                xyzful0(1:nND,1:6),xyzful(1:nND,1:6),dspful(1:nND,1:6), &
                                velful(1:nND,1:6),accful(1:nND,1:6),lodful(1:nND,1:6),  &
                                subdeltat,dampK,dampM,  &
                                triad_nn(1:3,1:3,1:nND),triad_ee(1:3,1:3,1:nEL),  &
                                triad_n1(1:3,1:3,1:nEL),triad_n2(1:3,1:3,1:nEL),  &
                                nND,nEL,nEQ,nMT,Newmarkdelta,Newmarkalpha,dtol,iterMax,alphaf)
        enddo
    enddo


end program main

subroutine StructureSolver(jBC,vBC,ele,prop,mss,xyzful0,xyzful,dspful,velful,accful,lodExteful,deltat,dampK,dampM,     &
                           triad_nn,triad_ee,triad_n1,triad_n2,nND,nEL,nEQ,nMT,Newmarkdelta,Newmarkalpha, &
                           dtol,iterMax,alphaf)
    implicit none
    
    integer:: nND,nEL,nEQ,nMT
    integer:: jBC(nND,6),ele(nEL,5)
    real(8):: vBC(nND,6),xyzful0(nND,6),xyzful(nND,6),prop(nMT,10)
    real(8):: mss(nEQ),dspful(nND,6),velful(nND,6),accful(nND,6),lodExteful(nND,6)
    real(8):: deltat
    real(8):: dampK,dampM

    real(8):: triad_nn(3,3,nND),triad_ee(3,3,nEL)
    real(8):: triad_n1(3,3,nEL),triad_n2(3,3,nEL)
!   ----------------------------------------------------------------------------------------------
    real(8):: dspO(nEQ),velO(nEQ),accO(nEQ) 
    real(8):: du(nEQ),ddu(nEQ),lodEffe(nEQ),lodInte(nEQ),lodExte(nEQ),dsp(nEQ),vel(nEQ),acc(nEQ)

    real(8):: wk1(nEQ),wk2(nEQ)
    real(8):: Newmarkdelta,Newmarkalpha
    real(8):: a0,a1,a2,a3,a4,a5,a6,a7
    real(8):: beta0,beta,gamma,zi,z0
    real(8):: dtol,dnorm
    real(8):: geoFRM(nEL)
    integer:: i,iter,iterMax,maxramp,iModify,i1,i2
    real(8):: alphaf
!   -----------------------------------------------------------------------------------------------
    a0 = 1.0d0/(Newmarkalpha*deltat*deltat)
    a1 = Newmarkdelta/(Newmarkalpha*deltat)
    a2 = 1.0d0/(Newmarkalpha*deltat)
    a3 = 1.0d0/(Newmarkalpha*2.0d0) - 1.0d0
    a4 = Newmarkdelta/Newmarkalpha - 1.0d0
    a5 = (Newmarkdelta/Newmarkalpha - 2.0d0)*0.5d0*deltat
    a6 = (1.0d0 - Newmarkdelta)*deltat
    a7 = Newmarkdelta*deltat

    beta0 =1.0d0
    gamma =1.0d0
    maxramp =0

    iModify = 1
!   -----------------------------------------------------------------------------------------------
    do  i=1,nND
        i1 = (i-1)*6+1
        i2 = (i-1)*6+6
        dsp(i1:i2)=dspful(i,1:6)
        vel(i1:i2)=velful(i,1:6)
        acc(i1:i2)=accful(i,1:6)
        lodExte(i1:i2)=lodExteful(i,1:6)
    enddo

    dspO(1:nEQ) = dsp(1:nEQ)    ! displacement at time t
    velO(1:nEQ) = vel(1:nEQ)    ! velocity at time t
    accO(1:nEQ) = acc(1:nEQ)    ! acceleration at time t
!   ------------------------------------------------------------------------------
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Begin loop over iterations:
!   -------------------------------------------------------------------
    iter=0
    dnorm=1.0d0
    do while(dnorm >= dtol .and. iter<= iterMax )
!       ------------------------------------------------------------------------------
!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       forces from body stresses
!       Assemble the element nodal force vector
!       -------------------------------------------------------------------
        call body_stress_D( lodInte,xyzful0(1:nND,1),xyzful0(1:nND,2),xyzful0(1:nND,3),xyzful(1:nND,1),xyzful(1:nND,2),xyzful(1:nND,3), &
                            ele,prop,triad_n1,triad_n2,triad_ee,    &
                            nND,nEL,nEQ,nMT,geoFRM &
                          )
!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Form the effective load vector {F}_(t+delta t)
!       ISBN 9787302388333 Xiong Zhang. P113-114
!       [C]=dampM*[M]+dampK*[K]
!       ISBN 9781441929105 James F. Doyle. P268
!       -------------------------------------------------------------------
        lodExte = (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 20000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -25000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
        do    i= 1, nEQ
            lodEffe(i)=lodExte(i)-lodInte(i)+(a0*(dspO(i)-dsp(i))+a2*vel(i)+a3*acc(i))*mss(i)   &
                                            +(a1*(dspO(i)-dsp(i))+a4*vel(i)+a5*acc(i))*dampM*mss(i)
        enddo

        if    (dampK .gt. 0.0d0) then
            do    i= 1, nEQ
                wk1(i)= a1*(dspO(i)-dsp(i)) +a4*vel(i) +a5*acc(i)
            enddo

            call matrixfree(jBC(1:nND,:),wk1,wk2,ele(1:nEL,1:5), &
                            xyzful0(1:nND,1),xyzful0(1:nND,2),xyzful0(1:nND,3), &
                            xyzful(1:nND,1),xyzful(1:nND,2),xyzful(1:nND,3), &
                            prop(1:nMT,:),nND,nEL,nEQ,nMT,geoFRM,alphaf,gamma,a0,a1,dampM,dampK,0,lodEffe,vBC,iter)

            do    i= 1, nEQ
                lodEffe(i) = lodEffe(i) + dampK*wk2(i)
            enddo
        endif
!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       apply force boundary conditions
!       multiplied with bigger number method
!       -------------------------------------------------------------------
        call matrixfree(jBC(1:nND,:),wk1,wk2,ele(1:nEL,1:5), &
                        xyzful0(1:nND,1),xyzful0(1:nND,2),xyzful0(1:nND,3), &
                        xyzful(1:nND,1),xyzful(1:nND,2),xyzful(1:nND,3), &
                        prop(1:nMT,:),nND,nEL,nEQ,nMT,geoFRM,alphaf,gamma,a0,a1,dampM,dampK,2,lodEffe,vBC,iter)

!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Conjugate Gradient Method
!       Iterative solution for displacement vector.
!       -------------------------------------------------------------------
        call cg(du,lodEffe,jBC(1:nND,1:6),ele(1:nEL,1:5), &
                xyzful0(1:nND,1),xyzful0(1:nND,2),xyzful0(1:nND,3), &
                xyzful(1:nND,1),xyzful(1:nND,2),xyzful(1:nND,3), &
                prop(1:nMT,1:10),nND,nEL,nEQ,nMT, &
                geoFRM,alphaf,gamma,a0,a1,dampM,dampK,lodEffe,vBC,iter &
                ) 

!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Update the displacements
!       -------------------------------------------------------------------  
        if    (iter <= maxramp) then
            zi=2.0d0**(iter)
            z0=2.0d0**(maxramp)
            beta=zi/z0*beta0
        else
            beta=1.0d0*beta0
        endif

        ddu(1:nEQ)= beta*du(1:nEQ)
        dsp(1:nEQ)= dsp(1:nEQ) + ddu(1:nEQ)
        
        do  i=1,nND
            dspful(i,1:6)=dsp((i-1)*6+1:(i-1)*6+6)
        enddo

        xyzful(1:nND,1:6)=xyzful0(1:nND,1:6)+dspful(1:nND,1:6)

!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Update geometry and triads
!       -------------------------------------------------------------------  
        call update_triad_D(ele,ddu,triad_nn,triad_n1,triad_n2,nND,nEL,nEQ)

        call make_triad_ee(ele,xyzful(1:nND,1),xyzful(1:nND,2),xyzful(1:nND,3),triad_ee,triad_n1,triad_n2,nND,nEL)

!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       test for convergence
!       -------------------------------------------------------------------
        dnorm=dabs(maxval((du(1:nEQ)*beta)**2))
        iter=iter+1
    enddo

    open(unit=22,file='5.dat',position='append')
    write(22,'(24E28.5)')dspful
    close(22)

    open(unit=1111,file='6.dat',position='append')
    write(1111,'(24E28.5)') xyzful
    close(1111)

    acc(1:nEQ)  = a0*(dsp(1:nEQ) - dspO(1:nEQ)) -a2*velO(1:nEQ) - a3*accO(1:nEQ)
    vel(1:nEQ)  = velO(1:nEQ) + a6*accO(1:nEQ) + a7*acc(1:nEQ)

    do  i=1,nND
        velful(i,1:6)=vel((i-1)*6+1:(i-1)*6+6)
        accful(i,1:6)=acc((i-1)*6+1:(i-1)*6+6)        
    enddo

endsubroutine


!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   READ structural DaTafile
subroutine readdt(jBC,ele,xyzful0,prop,nND,nEL,nEQ,nMT,idat)
    implicit none
    integer:: nND,nEL,nEQ,nMT,idat
    integer:: ele(nEL,5),jBC(nND,6)
    real(8):: xyzful0(nND,6),prop(nMT,10)
!   ---------------------------------------------------------------------------
    integer:: i,j,nbc,node,nmp
    character*50 endin
!   -----------------------------------------------------------------------------------------------
!   READ  node
    read(idat,*)  nND
    do    i= 1, nND
        read(idat,*) node,xyzful0(node,1),xyzful0(node,2),xyzful0(node,3)
    enddo
    read(idat,'(1a50)') endin
!   -----------------------------------------------------------------------------------------------
!   READ elem data
    ! element number, node left, node right, node right, element type, material property index
    read(idat,*) nEL
    do  i= 1, nEL
        read(idat,*) j,ele(j,1:5)
    enddo
    read(idat,'(1a50)') endin

!   -----------------------------------------------------------------------------------------------
!    READ  bcs  default is 0=free
    jBC(1:nND,1:6) = 0
    read(idat,*)  nbc
    do  i=1,nbc
        read(idat,*)node,jBC(node,1),jBC(node,2),jBC(node,3), &
                         jBC(node,4),jBC(node,5),jBC(node,6)

    enddo
    read(idat,'(1a50)') endin
!   -----------------------------------------------------------------------------------------------
!   READ element material properties
    ! can have nMT types of material
    ! property data will be overwritten if isKB = 0 or 1
    read(idat,*) nMT
    do    i= 1, nMT
        read(idat,*) nmp,prop(nmp,1:8)
    enddo
    read(idat,'(1a50)') endin
!   -----------------------------------------------------------------------------------------------

    nEQ=nND*6

    !write(*,'(3(A,1x,I8,2x))')'nND=',nND,'nEL=',nEL,'nEQ=',nEQ
    !write(*,'(3(A,1x,I8,2x))')'nMT=',nMT,'nBD=',nBD,'nSTF=',nSTF
    return
endsubroutine

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Conjugate Gradient Method
!    Iterative solution for displacement vector.
!
!    https://www.detailedpedia.com/wiki-Conjugate_gradient_method
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine cg(x,b,jBC,ele,xord0,yord0,zord0,xord,yord,zord,prop,nND,nEL,nEQ,nMT, &
              geoFRM,alphaf,gamma,Newmarka0,Newmarka1,dampM,dampK,lodEffe,vBC,Newmarkiter)
    implicit none
    integer:: nEQ
    real(8):: b(nEQ)

    real(8):: x(nEQ), r(nEQ), p(nEQ), Ap(nEQ)
    integer :: iter, max_iter
    double precision :: alpha, beta, rsold, rsnew

    integer:: nND,nEL,nMT
    integer:: ele(nEL,5)
    real(8):: prop(nMT,10)
    real(8):: xord0(nND), yord0(nND),zord0(nND)
    real(8):: xord(nND), yord(nND),zord(nND)

    integer:: jBC(nND,6)
    double precision :: start, finish, total_time

    real(8):: geoFRM(nEL),alphaf,gamma,Newmarka0,Newmarka1,dampM,dampK,lodEffe(nEQ),vBC(nND,6)
    integer:: Newmarkiter
    
    
    ! initialize the solution vector x to 0.
    x = 0.0
    call CPU_TIME(start)
    ! initialize the residual vector. r = b - matmul(A, x)
    call matrixfree(jBC(1:nND,:),x,Ap,ele(1:nEL,1:5),xord0(1:nND),yord0(1:nND),zord0(1:nND),xord(1:nND),yord(1:nND),zord(1:nND), &
                        prop(1:nMT,:),nND,nEL,nEQ,nMT,geoFRM,alphaf,gamma,Newmarka0,Newmarka1,dampM,dampK,1,lodEffe,vBC,Newmarkiter)    
    r=b-Ap

    p = r

    ! maximum iteration count.
    max_iter = 10000

    rsold = dot_product(r, r)

    do iter = 1, max_iter
        ! matrix free method, Calculate the matrix-vector product. Ap = matmul(A, p)
        call matrixfree(jBC(1:nND,:),p,Ap,ele(1:nEL,1:5),xord0(1:nND),yord0(1:nND),zord0(1:nND),xord(1:nND),yord(1:nND),zord(1:nND), &
                        prop(1:nMT,:),nND,nEL,nEQ,nMT,geoFRM,alphaf,gamma,Newmarka0,Newmarka1,dampM,dampK,1,lodEffe,vBC,Newmarkiter)
        alpha = rsold / dot_product(p, Ap)
        x = x + alpha * p
        r = r - alpha * Ap
        rsnew = dot_product(r, r)
        if(sqrt(rsnew) < 1e-10) exit
        beta = rsnew/rsold
        p = r + beta * p
        rsold = rsnew
    enddo

    write(*,*) 'iter is', iter

    open(unit = 111, file = 'x.dat', position = 'append')
    write(111,'(24E28.5)') x
    close(111)
    
    
    call CPU_TIME(finish)  ! Get the end time.
    total_time = finish - start  ! Calculate the difference, namely the program running time.
    write(*,*) 'Program ran for ', total_time, ' seconds.'
    return
endsubroutine

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Matrix Free Method
!    Implement matrix-vector multiplication for individual elements, without assembling the global matrix.
!
!    Element‐by‐Element(Hughes, Thomas J. R.(1983).doi:10.1061/(ASCE)0733-9399(1983)109:2(576))
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine matrixfree(jBC,x,b,ele,xord0,yord0,zord0,xord,yord,zord,prop,nND,nEL,nEQ,nMT, &
                      geoFRM,alphaf,gamma,Newmarka0,Newmarka1,dampM,dampK,isjBC,lodEffe,vBC,iter)
    implicit none
    integer:: nND,nEL,nEQ,nMT
    integer:: ele(nEL,5)
    real(8):: prop(nMT,10)
    real(8):: xord0(nND), yord0(nND),zord0(nND)
    real(8):: xord(nND), yord(nND),zord(nND)
    real(8):: et(18,18),ek(18,18), ek12(12,12),eg(18,18), eg12(12,12),em(18,18), em12(12,12)
    real(8):: s,geoFRM(nEL),alphaf,gamma,Newmarka0,Newmarka1,dampM,dampK,lodEffe(nEQ),vBC(nND,6)
    integer:: isjBC,iter

    integer:: i,j,n,i1,j1,k1,mat,nELt
    real(8):: e0,g0,a0,r0,b0,zix0,ziy0,ziz0,dx,dy,dz,xl0,xl,xll,xmm,xnn,xl9,xll0,xmm0,xnn0

    real(8):: b(nEQ),x(nEQ),x18(18),b18(18)
    integer:: idof(18),jBC(nND,6),iND,iEQ,iEQ2
    b=0.0
    ! form each element matrix, and assemble
    do    n=1,nEL
        i1  = ele(n,1)
        j1  = ele(n,2)
        k1  = ele(n,3)
        nELt= ele(n,4)
        mat = ele(n,5)

        do  i= 1, 6
            idof(i)    = (i1-1)*6 + i
            idof(i+6)  = (j1-1)*6 + i
        enddo
        !
        if    (nELt .eq. 1 .OR. nELt .eq. 2 .OR. nELt .eq. 3) then
            ! frame
            ! material props not change
            e0=prop(mat,1)
            g0=prop(mat,2)
            a0=prop(mat,3)
            r0=prop(mat,4)
            b0=prop(mat,5)
            zix0=prop(mat,6)
            ziy0=prop(mat,7)
            ziz0=prop(mat,8)
            dx= xord0(j1) - xord0(i1)
            dy= yord0(j1) - yord0(i1)
            dz= zord0(j1) - zord0(i1)
            xl0=dsqrt(dx*dx+dy*dy+dz*dz)
            xll0=dx/xl0
            xmm0=dy/xl0
            xnn0=dz/xl0
            !
            ! orientation
            dx= xord(j1) - xord(i1)
            dy= yord(j1) - yord(i1)
            dz= zord(j1) - zord(i1)
            xl=dsqrt(dx*dx+dy*dy+dz*dz)
            xll=dx/xl
            xmm=dy/xl
            xnn=dz/xl

            ! Calculate the linear stiffness matrix
            xl9= xl0
            call elmstfFRM_D(xl9,zix0,ziy0,ziz0,a0,e0,g0,ek12, nELt)
            ! use trans3d to transform from local to global stifness
            call trans3d_D(xll,xmm,xnn,ek12,b0)
            ! expand to [18x18]
            ek(1:18,1:18)=0.0
            ek(1:12,1:12)=ek12(1:12,1:12)

            ! Calculate the geometric stiffness matrix
            !read(igeoFRM,*) sxx
            s=geoFRM(n)
            xl9= xl0
            call elmgeomFRM_D(xl9,eg12,s)
            call trans3d_D(xll,xmm,xnn,eg12,b0)
            eg(1:18,1:18)=0.0d0
            eg(1:12,1:12)=eg12(1:12,1:12)

            ! Calculate the mass matrix
            xl9= xl0
            call elmmasFRM_D(r0,a0,xl9,zix0,em12,alphaf)
            call trans3d_D(xll0,xmm0,xnn0,em12,b0)
            em(1:18,1:18)=0.0d0
            em(1:12,1:12)=em12(1:12,1:12)

!           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           Form the tangent stiffness matrix
!           ISBN 9781441929105 James F. Doyle. P268
!           [C]=dampM*[M]+dampK*[K]
!           -------------------------------------------------------------------
!           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           ISBN 9781441929105 James F. Doyle. P193,196
!           K[T] = K[E] + gamma*K[G]
!           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            et = ek + gamma * eg
            et = et + Newmarka0 * em + Newmarka1 * dampM * em
            if(dampK .gt. 0.0d0) then
                et = et + Newmarka1 * dampK * ek
            endif

!           -------------------------------------------------------------------
!           apply boundary conditions
!           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
            if    (isjBC .eq. 1) then
                do  i=1, 3
                    do  j=1, 6
                        iND = ele(n,i)
                        if(jBC(iND,j).gt.0)then
                            iEQ = (i-1)*6+j
                            et(iEQ,iEQ)=et(iEQ,iEQ)*1.0d20
                        endif
                    enddo
                enddo
!               -------------------------------------------------------------------
!               Multiply the element stiffness matrix by the displacement vector, and assemble the force vector
!               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
                do i = 1, 12
                    x18(i)=x(idof(i))
                enddo
                b18=matmul(et,x18)
                do i = 1, 12
                    b(idof(i))=b(idof(i))+b18(i)
                enddo
            else  if  (isjBC .eq. 0) then
                do i = 1, 12
                    x18(i)=x(idof(i))
                enddo
                b18=matmul(ek,x18)
                do i = 1, 12
                    b(idof(i))=b(idof(i))+b18(i)
                enddo
            else if (isjBC .eq. 2) then
                do  i=1, 3
                    do  j=1, 6
                        iND = ele(n,i)
                        if(jBC(iND,j).gt.0)then
                            iEQ = (i-1)*6+j
                            iEQ2 = (iND-1)*6+j
                            lodEffe(iEQ2) = 0.0d0
                            if(iter==0)then
                                lodEffe(iEQ2)=lodEffe(iEQ2) + et(iEQ,iEQ)*vBC(iND,j)*1.0d20
                            else
                                lodEffe(iEQ2)=et(iEQ,iEQ)*0.0d0
                            endif
                        endif
                    enddo
                enddo
            endif
        !
        else
            write(*,*)'not this nELt:',nELt
            stop
        endif
    enddo

    return    
endsubroutine

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    set angle of initial triads
!    
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
subroutine init_triad_D(ele,xord0,yord0,zord0,triad_nn,triad_n1,triad_n2,triad_ee,nND,nEL)
    implicit none
    integer:: nND,nEL         
    integer:: ele(nEL,5)
    real(8):: xord0(nND),yord0(nND),zord0(nND)
    real(8):: triad_nn(3,3,nND),triad_ee(3,3,nEL)
    real(8):: triad_n1(3,3,nEL),triad_n2(3,3,nEL)
    !
    real(8):: dx,dy,dz,xl0
    real(8):: xll1,xmm1,xnn1
    real(8):: dd

    integer:: i,j,n,i1,j1,k1,nELt

    ! For each node, set the triad  to global system
    do  n=1,nND

        do i=1,3
        do j=1,3
            triad_nn(i,j,n)=0.0d0
        enddo
        enddo
        triad_nn(1,1,n)=1.0d0
        triad_nn(2,2,n)=1.0d0
        triad_nn(3,3,n)=1.0d0

    enddo
    !
    ! For each element, calculate the current orientation triad
    do  n=1,nEL

        i1  = ele(n,1)
        j1  = ele(n,2)
        k1  = ele(n,3)
        nELt= ele(n,4)
        !
        if    (nELt .eq. 1 .OR. nELt .eq. 2 .OR. nELt .eq. 3) then
            ! frame
            dx = xord0(j1) - xord0(i1)
            dy = yord0(j1) - yord0(i1)
            dz = zord0(j1) - zord0(i1)
            xl0=dsqrt(dx*dx+dy*dy+dz*dz)
            xll1=dx/xl0
            xmm1=dy/xl0
            xnn1=dz/xl0
            dd=dsqrt(xll1*xll1+xmm1*xmm1)
            triad_n1(1,1,n)=xll1
            triad_n1(2,1,n)=xmm1
            triad_n1(3,1,n)=xnn1
            !
            if    (dd .lt. 0.001d0) then
                triad_n1(1,2,n)=0.0d0
                triad_n1(2,2,n)=1.0d0
                triad_n1(3,2,n)=0.0d0
                triad_n1(1,3,n)=-xnn1
                triad_n1(2,3,n)=0.00d0
                triad_n1(3,3,n)=0.00d0
            else
                triad_n1(1,2,n)=-xmm1/dd
                triad_n1(2,2,n)=+xll1/dd
                triad_n1(3,2,n)=0.0d0

                triad_n1(1,3,n)=-xll1*xnn1/dd
                triad_n1(2,3,n)=-xmm1*xnn1/dd                                    
                triad_n1(3,3,n)= dd 
            endif
            ! all element triads have same initial orientation
            triad_n2(1:3,1:3,n)=triad_n1(1:3,1:3,n)
            triad_ee(1:3,1:3,n)=triad_n1(1:3,1:3,n)
        else
            write(*,*)'not this nELt:',nELt
            stop
        endif
    enddo

    return
endsubroutine


!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    update angle of  triads
!    triad_n1: the triad of node 1 in beam
!    triad_n2: the triad of node 2 in beam
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
subroutine update_triad_D(ele,ddut,triad_nn,triad_n1,triad_n2,nND,nEL,nEQ)
    implicit none
    integer:: nND,nEL,nEQ         
    integer:: ele(nEL,5)
    real(8):: ddut(nEQ)
    real(8):: triad_nn(3,3,nND),triad_n1(3,3,nEL),triad_n2(3,3,nEL)
   
    real(8):: ddutful(nND,6),rr(3,3)
    real(8):: dtx1,dty1,dtz1
    integer:: i,n,i1,j1,k1,node,jdof,ireduc,nELt

    do    i=1,nND*6
        node=(i+5)/6
        jdof=i-(node-1)*6
        ireduc=i
        ddutful(node,jdof)=ddut(ireduc)
    enddo
    !
    ! For each node, set the triad_nn = [R]triad_nn
    do    n=1,nND
        dtx1=ddutful(n,4)
        dty1=ddutful(n,5)
        dtz1=ddutful(n,6)
        call finite_rot(dtx1,dty1,dtz1,rr)
        triad_nn(1:3,1:3,n)=matmul(rr(1:3,1:3),triad_nn(1:3,1:3,n))
    enddo

    ! For each element, calculate the current orientation triad
    do  n=1,nEL
        i1  = ele(n,1)
        j1  = ele(n,2)
        k1  = ele(n,3)
        nELt= ele(n,4)
        ! n1 node
        dtx1=ddutful(i1,4)
        dty1=ddutful(i1,5)
        dtz1=ddutful(i1,6)
        call finite_rot(dtx1,dty1,dtz1,rr)
        triad_n1(1:3,1:3,n)=matmul(rr(1:3,1:3),triad_n1(1:3,1:3,n))
        ! n2 node
        dtx1=ddutful(j1,4)
        dty1=ddutful(j1,5)
        dtz1=ddutful(j1,6)
        call finite_rot(dtx1,dty1,dtz1,rr)
        triad_n2(1:3,1:3,n)=matmul(rr(1:3,1:3),triad_n2(1:3,1:3,n))
    enddo

    return
endsubroutine

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Get orientation of element
!    ISBN 9781441929105 James F. Doyle. P189-191
!    triad_n1: the triad of node 1 in beam
!    triad_n2: the triad of node 2 in beam
!    triad_ee: the triad of beam element(beam center)
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
subroutine make_triad_ee(ele,xord,yord,zord,triad_ee,triad_n1,triad_n2,nND,nEL)
    implicit none
    integer:: nND,nEL
    integer:: ele(nEL,5)
    real(8):: xord(nND), yord(nND), zord(nND)
    !
    real(8):: triad_aa(3,3)
    real(8):: triad_ee(3,3,nEL)
    real(8):: triad_n1(3,3,nEL),triad_n2(3,3,nEL)
    real(8):: rr(3,3),tx,ty,tz
    real(8):: triad_11(3,3),triad_22(3,3)
    !
    real(8):: dx,dy,dz,xl0
    
    real(8):: xll,xmm,xnn,dd,r2e1,r3e1
    integer:: i,j,k,n,i1,j1,k1,nELt
    !
    !
    ! For each element, calculate the current orientation triad
    do    n=1,nEL
        i1  = ele(n,1)
        j1  = ele(n,2)
        k1  = ele(n,3)
        nELt= ele(n,4)
        !
        if    (nELt .eq. 1 .OR. nELt .eq. 2 .OR. nELt .eq. 3) then
            ! frame
            dx = xord(j1) - xord(i1)
            dy = yord(j1) - yord(i1)
            dz = zord(j1) - zord(i1)
            xl0=dsqrt(dx*dx+dy*dy+dz*dz)
            xll=dx/xl0
            xmm=dy/xl0
            xnn=dz/xl0
            dd =dsqrt(xll*xll+xmm*xmm)
            do    i=1,3
            do    j=1,3
                triad_ee(i,j,n)=0.0d0
            enddo
            enddo
            triad_ee(1,1,n)=xll
            triad_ee(2,1,n)=xmm
            triad_ee(3,1,n)=xnn
            !
            ! get angle between two triads
            do    i=1,3
            do    j=1,3
                triad_11(i,j)=triad_n1(i,j,n)
                triad_22(i,j)=triad_n2(i,j,n)
            enddo
            enddo
            call get_angle_triad( triad_11,triad_22,tx,ty,tz)
            !
            ! rotate n1 to intermediate
            tx=tx/2.0d0
            ty=ty/2.0d0
            tz=tz/2.0d0
            call finite_rot(tx,ty,tz,rr)
            triad_aa(1:3,1:3)=matmul(rr(1:3,1:3),triad_n1(1:3,1:3,n))
            !
            ! vectors e2 e3
            r2e1 = 0.0d0
            r3e1 = 0.0d0
            do    k=1,3
                r2e1 = r2e1 + triad_aa(k,2)*triad_ee(k,1,n)
                r3e1 = r3e1 + triad_aa(k,3)*triad_ee(k,1,n)
            enddo
            do    j=1,3
                triad_ee(j,2,n)=triad_aa(j,2) - r2e1*(triad_aa(j,1)+triad_ee(j,1,n))/2.0d0
                triad_ee(j,3,n)=triad_aa(j,3) - r3e1*(triad_aa(j,1)+triad_ee(j,1,n))/2.0d0
            enddo
        else
            write(*,*)'not this nELt:',nELt
            stop
        endif
    enddo
    ! end of loop over elements
    return
endsubroutine

!
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Finite rotation
!    ISBN 9781441929105 James F. Doyle. P184
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
subroutine finite_rot(t1,t2,t3,rr)
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
endsubroutine


!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    GET ANGLE of between TRIADs
!    ISBN 9781441929105 James F. Doyle. P186
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
subroutine get_angle_triad(triad_n1,triad_n2,tx,ty,tz)
    implicit none
    real(8):: triad_n1(3,3),triad_n2(3,3)
    real(8):: rr(3,3)
    real(8):: tx,ty,tz, dtx,dty,dtz,c1,tt,sint
    integer:: i,j,k
    !
    ! get angle between two triads
    do    i=1,3
    do    j=1,3
        rr(i,j)=0.0d0
        do    k=1,3
            rr(i,j)=rr(i,j) + triad_n2(i,k)*triad_n1(j,k)
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
endsubroutine

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    triad: from global to local
!
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
subroutine global_to_local(triad,tx,ty,tz,tx2,ty2,tz2)
    implicit none
    real(8):: triad(3,3)
    real(8):: tx,ty,tz,tx2,ty2,tz2
    tx2 = triad(1,1)*tx+triad(2,1)*ty+triad(3,1)*tz
    ty2 = triad(1,2)*tx+triad(2,2)*ty+triad(3,2)*tz
    tz2 = triad(1,3)*tx+triad(2,3)*ty+triad(3,3)*tz
    return
endsubroutine

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Nodal loads due to body stresses
!   input: xord0,yord0,zord0,xord,yord,zord undeformed and deformed coordinates
!   ele elements [nodex, nodey]
!   prop material property
!   nMT number of material types
!   triad_n1,triad_n2,triad_ee
!   output: gforce global force
!   geoFRM
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine body_stress_D(gforce,xord0,yord0,zord0,xord,yord,zord,ele,prop,triad_n1,triad_n2,triad_ee, &
                            nND,nEL,nEQ,nMT,geoFRM &
                            )
    implicit none
    integer:: nMT,nEL,nND,nEQ
    integer:: ele(nEL,5)
    real(8):: xord0(nND), yord0(nND), zord0(nND)
    real(8):: xord(nND), yord(nND), zord(nND)
!
    real(8):: prop(nMT,10),geoFRM(nEL)
 
    real(8):: force(18),forceb(18)
    real(8):: gforce(nEQ)
    real(8):: ekb12(12,12)
!
    real(8):: triad_ee(3,3,nEL)
    real(8):: triad_n1(3,3,nEL),triad_n2(3,3,nEL)

    real(8):: triad_00(3,3),triad_11(3,3),triad_22(3,3)
    real(8):: rr(3,3)
    real(8):: ub(18),dl
!
    real(8):: dx0,dy0,dz0,du,dv,dw
    real(8):: dx,dy,dz,xl0
    real(8):: tx1,tx2,ty1,ty2,tz1,tz2,tx,ty,tz

    real(8):: e0,g0,a0,b0,r0,zix0,ziy0,ziz0,xl,fxx
    integer:: i,j,k,n,i1,j1,k1,mat,nELt

    ! pi=4.0*datan(1.0d0)

    gforce(1:nEQ)=0.0
  
    ! rewind(igeoFRM)
    ! rewind(igeoPLT)


    ! For each element, calculate the nodal forces
    do    n=1,nEL

        i1  = ele(n,1)
        j1  = ele(n,2)
        k1  = ele(n,3)        
        nELt= ele(n,4)
        mat = ele(n,5)
        ! Define the nodes and cell properties of each cell (e.g., beam cell triangle cell)
        if    (nELt .eq. 1 .OR. nELt .eq. 2 .OR. nELt .eq. 3) then
            ! frame
            e0  =prop(mat,1)
            g0  =prop(mat,2)
            a0  =prop(mat,3)
            r0  =prop(mat,4)
            b0  =prop(mat,5)
            zix0=prop(mat,6)
            ziy0=prop(mat,7)
            ziz0=prop(mat,8)

            dx0 = xord0(j1) - xord0(i1)
            dy0 = yord0(j1) - yord0(i1)
            dz0 = zord0(j1) - zord0(i1)
            xl0 = dsqrt(dx0*dx0+dy0*dy0+dz0*dz0)
            !Define the xyz initial directional distance and the initial beam length, xl0, for the two end nodes of the beam unit
            !
            ! orientation
            du = (xord(j1)-xord0(j1))-(xord(i1)-xord0(i1))
            dv = (yord(j1)-yord0(j1))-(yord(i1)-yord0(i1))
            dw = (zord(j1)-zord0(j1))-(zord(i1)-zord0(i1))
            ! Define the axial strain in the xyz direction (difference between displacements at two nodes)

            dx = dx0 + du
            dy = dy0 + dv
            dz = dz0 + dw
            xl =dsqrt(dx*dx+dy*dy+dz*dz)
            ! Define the the xyz directional distance of the two nodes of the beam element after deformation
            ! as well as the length, xl, of the element after deformation
            !
            dl = ( (2*dx0+du)*du +(2*dy0+dv)*dv +(2*dz0+dw)*dw )/ (xl+xl0)
            ! get twisting angles
            do    i=1,3
            do    j=1,3
                triad_00(i,j)=triad_ee(i,j,n)
                triad_11(i,j)=triad_ee(i,j,n)
                triad_22(i,j)=triad_n1(i,j,n)
            enddo
            enddo
            call get_angle_triad(triad_11,triad_22,tx,ty,tz)
            call global_to_local(triad_00,tx,ty,tz,tx1,ty1,tz1)
            !
            do    i=1,3
            do    j=1,3
                triad_11(i,j)=triad_ee(i,j,n)
                triad_22(i,j)=triad_n2(i,j,n)
            enddo
            enddo
            call get_angle_triad(triad_11,triad_22,tx,ty,tz)
            call global_to_local(triad_00,tx,ty,tz,tx2,ty2,tz2)

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
            ! get current stiffness in local coords. use L0

            call elmstfFRM_D(xl0,zix0,ziy0,ziz0,a0,e0,g0,ekb12,nELt )
            !
            ! compute axial force
            fxx=dl*e0*a0/xl0
            ! save local force for geo stiff
            ! write(igeoFRM,'(D25.15)') fxx
            geoFRM(n)=fxx
            ! nodal forces in local coords
            ! {F}=[k]{u}
            forceb(1:12) =matmul(ekb12(1:12,1:12),ub(1:12))
            forceb(13:18)=0.0
            !
            ! transform to global
            do  i=1,3
                do  j=1,3
                    rr(i,j)=triad_ee(i,j,n)
                enddo
            enddo
            do  i=1,18
                force(i)=0.0
            enddo
            do    i=1,3
                do    k=1,3
                    force(0+i) = force(0+i) + rr(i,k)*forceb(0+k)
                    force(3+i) = force(3+i) + rr(i,k)*forceb(3+k)
                    force(6+i) = force(6+i) + rr(i,k)*forceb(6+k)
                    force(9+i) = force(9+i) + rr(i,k)*forceb(9+k)
                enddo
            enddo
            !
            call assembFOR(nEQ,gforce,force,i1,j1,k1)
            !
        else
            write(*,*)'not this nELt:',nELt
            stop
        endif

    enddo
    !
    return
endsubroutine

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    ELeMent STiFfness for FRaMe
!    calculates the element stiffness matrices.
!
!    ISBN 9787040258417 Zeng Pan. P70
!    ISBN 9780792312086 James F. Doyle. P81
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine elmstfFRM_D(length, ix, iy, iz,area, emod, gmod, ek ,nELt)
    implicit none
    real(8):: area, length, Invlength, ix, iy, iz, emod, gmod
    real(8):: ek(12,12)
    integer:: nELt

    integer:: i,j
    real(8):: emlen,emlen2,emlen3
    !
    ! initialize all ek elements to zero
    ek(1:12,1:12)=0.0
    !
    ! STIFFNESS matrix in local coordinates
    ! emod is Modulus of elasticity
    !
    Invlength = 1/length
    emlen  = emod*Invlength
    emlen2 = emlen*Invlength
    emlen3 = emlen2*Invlength

    if (nELt .eq. 1) then
        ek(1,1)   =   area*emlen
    !
    elseif (nELt .eq. 2) then
        ek(1,1)   =   area*emlen
        ek(2,2)   =   area*emlen
    !
    elseif (nELt .eq. 3) then
        ek(1,1)   =   area*emlen
        ek(2,2)   =   12.0*emlen3*iz
        ek(3,3)   =   12.0*emlen3*iy
        ek(4,4)   =   gmod*ix/length
        ek(5,5)   =   4.0*emlen*iy
        ek(6,6)   =   4.0*emlen*iz
    !
        ek(2,6)   =   6.0*emlen2*iz
        ek(3,5)   =  -6.0*emlen2*iy
    else
        write(*,*)'not this nELt:',nELt
        stop
    endif
    !
    ek(7,7)   =   ek(1,1)
    ek(8,8)   =   ek(2,2)
    ek(9,9)   =   ek(3,3)
    ek(10,10) =   ek(4,4)
    ek(11,11) =   ek(5,5)
    ek(12,12) =   ek(6,6)
    !
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
    !
    ek(8,12)  =   -ek(2,6)
    ek(9,11)  =   -ek(3,5)
    !
    ! impose the symmetry
    do  i= 1, 12
        do  j= i, 12
            ek(j,i) = ek(i,j)
        enddo
    enddo
    !
    return
endsubroutine

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   ELeMent MASs matrix for the FRaMe
!   Lumped mass matrix
!   ISBN 9781441929105 James F. Doyle. P273
!   ISBN 9780792312086 James F. Doyle. P423
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine elmmasFRM_D(rho,area,length,zix,em,alphaf)
    implicit none
    real(8):: rho, area, length,zix
    real(8):: em(12,12)
    real(8):: alphaf,roal
    !
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
    return
endsubroutine

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   ELeMent GEOMetric stiffness matrix for a FRaMe
!   ISBN 9781441929105 James F. Doyle. P217,228,229
!   ISBN 9780792312086 James F. Doyle. P424
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! s tension force
    ! element initial length
    subroutine elmgeomFRM_D(length,eg,s)
    implicit none
    real(8):: length,eg(12,12),s
    real(8):: emlenz,alpha
    integer:: i,j
    !
    ! initialize all eg elements to zero
    eg(1:12,1:12)=0.0d0
    ! Stiffness matrix in local coordinates
    ! if (s .gt. 200.) s=200.
    ! emlenz  =   s/(30.0*length)
    ! alpha   =   s*1.0e-0
    ! beta    =   1.0
    alpha   =   (s/length)*1.0d-6
    emlenz  =   s/(30.0*length)
    ! emlenz  =   s/length
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
    !
    eg(7,7)   =   eg(1,1)
    eg(8,8)   =   eg(2,2)
    eg(9,9)   =   eg(3,3)
    eg(10,10) =   eg(4,4)
    eg(11,11) =   eg(5,5)
    eg(12,12) =   eg(6,6)
    !
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
    !
    eg(8,12)  =   -eg(2,6)
    eg(9,11)  =   -eg(3,5)
    !
    ! impose the symmetry
    do  i= 1, 12
    do  j= i, 12
        eg(j,i) = eg(i,j)
    enddo
    enddo
    !
    ! check diagonal terms
    do    i=1,12
        if (dabs(eg(i,i)) .lt. 1.0d-20) eg(i,i)=1.0d-20
    enddo
    !
    return
endsubroutine

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    ISBN 9780792312086 James F. Doyle. P83-85
!
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   TRANSformation of matrix in 3D. L->G
subroutine trans3d_D(l,m,n,ek,beta)
    implicit none
    real(8):: ek(12,12),rt(3,3),r(3,3),ktemp(12,12)
    real(8):: m,n,l,beta,pi,sb,cb,d,Invd
    integer:: i,j,k,j1,j2,ii,jj,in,jn

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
    !
    do  in=1,3
    do  jn=1,3
        rt(jn,in)=r(in,jn)
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
endsubroutine

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   ASSEMBle consistent FORce for pressurized flex plate
!   input: a local value; nEQ length of aa, number of dof (6 * elements);
!   output:: aa global value
subroutine assembFOR(nEQ, aa, a,i1, j1,k1)
    implicit none
    integer:: nEQ   
    real(8):: aa(nEQ),a(18)
    integer:: i1,j1,k1
    integer:: idof(18),i,ieqn1
    !
    ! Set idof to posible DoF number of each nodes
    do    i= 1, 6
        idof(i)    = (i1-1)*6 + i
        idof(i+6)  = (j1-1)*6 + i
        idof(i+12) = (k1-1)*6 + i
    enddo
    !
    ! Store the values for individual array in global array
    do  i= 1, 18
        ieqn1 = idof(i)
        aa(ieqn1) = aa(ieqn1) + a(i)
    enddo
    !
    return
endsubroutine

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   FORM MASS matrix: only lumped
!
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine formmass_D(ele,xord,yord,zord,prop,mss,nND,nEL,nEQ,nMT,alphaf)
    implicit none                  
    integer:: nND,nEL,nEQ,nMT
    integer:: ele(nEL,5)
    real(8):: prop(nMT,10)
    !
    real(8):: xord(nND), yord(nND), zord(nND)
    real(8):: mss(nEQ), em(18,18)
    real(8):: em12(12,12)
 
    integer:: i,i1,j1,k1,mat,nELt
    real(8):: a0,r0,b0,zix0,ziy0,ziz0,dx,dy,dz,xl,xll,xmm,xnn

    real(8):: alphaf

    ! zero array before assembling

    mss(1:nEQ)=0.0d0
    !
    ! form the element form matrix, and assemble
    do  i=1,nEL
        i1  = ele(i,1)
        j1  = ele(i,2)
        k1  = ele(i,3)        
        nELt= ele(i,4) ! element type, 2 line segment
        mat = ele(i,5)
        if    (nELt .eq. 1 .OR. nELt .eq. 2 .OR. nELt .eq. 3) then
            a0=prop(mat,3)
            r0=prop(mat,4)
            b0=prop(mat,5)
            zix0=prop(mat,6)
            ziy0=prop(mat,7)
            ziz0=prop(mat,8)
            dx= xord(j1) - xord(i1)
            dy= yord(j1) - yord(i1)
            dz= zord(j1) - zord(i1)
            xl=dsqrt(dx*dx+dy*dy+dz*dz)
            xll=dx/xl
            xmm=dy/xl
            xnn=dz/xl
            call elmmasFRM_D(r0,a0,xl,zix0,em12,alphaf)
            call trans3d_D(xll,xmm,xnn,em12,b0)
            em(1:18,1:18)=0.0d0
            em(1:12,1:12)=em12(1:12,1:12)
            call assembLUM(nEQ,mss,em,i1,j1,k1)
        else
            write(*,*)'not this nELt:',nELt
            stop
        endif
    enddo    
    return
endsubroutine

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   ASSEMBle LUMped element matrices
subroutine assembLUM(nEQ,aa,a,i1,j1,k1)
    implicit none
    integer:: nEQ
    real(8):: aa(nEQ),a(18,18)
    integer:: i1,j1,k1
    integer:: idof(18),i,ieqn1
    !
    ! Set idof to posible DoF number of each nodes
    do  i= 1, 6
        idof(i)    = (i1-1)*6 + i
        idof(i+6)  = (j1-1)*6 + i
        idof(i+12) = (k1-1)*6 + i
    enddo
    !
    ! Store the values for individual array in global array
    do  i= 1, 18
        ieqn1 = idof(i)
        aa(ieqn1) = aa(ieqn1) + a(i,i)
    enddo
    !
    return
endsubroutine