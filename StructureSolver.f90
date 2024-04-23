program StructureSolver

    integer, allocatable:: nND(:),nEL(:),nEQ(:),nMT(:),nBD(:),nSTF(:)
    integer, allocatable:: ele(:,:,:),nloc(:,:),nprof(:,:),nprof2(:,:),jBC(:,:,:)
    real(8), allocatable:: xyzful00(:,:,:),mssful(:,:,:),vBC(:,:,:),mss(:,:),prop(:,:,:)
    real(8), allocatable:: lodful(:,:,:),extful(:,:,:),extful1(:,:,:),extful2(:,:,:),grav(:,:,:),repful(:,:,:)
    real(8), allocatable:: xyzful0(:,:,:),xyzfulnxt(:,:,:),dspful(:,:,:),accful(:,:,:)
    real(8), allocatable:: xyzful(:,:,:),xyzfulIB(:,:,:),velful(:,:,:)
    real(8), allocatable:: du(:),lodEffe(:)
    integer :: idat, iFish, nFish, nND_max, nEL_max, nMT_max
    real(8):: triad_e0(3,3,1), triad_ee(3,3,1)
    idat = 12
    iFish = 1
    nFish = 1
    triad_e0 = 0.0
    triad_ee = 0.0

    allocate( nND(1:nFish),nEL(1:nFish),nMT(1:nFish),nEQ(1:nFish),nBD(1:nFish),nSTF(1:nFish))

    open(unit=idat, file = 'Beam.dat')
        rewind(idat)
        read(idat,*)
        read(idat,*)nND(iFish),nEL(iFish),nMT(iFish)
        read(idat,*)
    close(idat)

    nND_max=maxval(nND(:))
    nEL_max=maxval(nEL(:))
    nMT_max=maxval(nMT(:))

    allocate( ele(1:nFish,nEL_max,5),nloc(1:nFish,nND_max*6),nprof(1:nFish,nND_max*6),nprof2(1:nFish,nND_max*6),jBC(1:nFish,nND_max,6))
    allocate( xyzful00(1:nFish,nND_max,6),xyzful0(1:nFish,nND_max,6),mssful(1:nFish,nND_max,6),lodful(1:nFish,nND_max,6), &
              extful(1:nFish,nND_max,6),extful1(1:nFish,nND_max,6),extful2(1:nFish,nND_max,6),grav(1:nFish,nND_max,6),repful(1:nFish,nND_max,6),    &
              vBC(1:nFish,nND_max,6),mss(1:nFish,nND_max*6),prop(1:nFish,nMT_max,10))
    allocate( xyzful(1:nFish,nND_max,6),xyzfulIB(1:nFish,nND_max,6),xyzfulnxt(1:nFish,nND_max,6),dspful(1:nFish,nND_max,6),velful(1:nFish,nND_max,6),accful(1:nFish,nND_max,6))



    open(unit=idat, file = 'Beam.dat')
    rewind(idat)
    read(idat,*)
    read(idat,*)nND(iFish),nEL(iFish),nMT(iFish)
    read(idat,*)
    !write(*,*)'read FEMeshFile ',iFish
!   ===============================================================================================
    call readdt(jBC(iFish,1:nND(iFish),1:6),ele(iFish,1:nEL(iFish),1:5),xyzful00(iFish,1:nND(iFish),1:6),prop(iFish,1:nMT(iFish),1:10), &
                nND(iFish),nEL(iFish),nEQ(iFish),nMT(iFish),idat)
    close(idat)


    do  iND=1,nND(iFish)
        xyzful0(iFish,iND,1:6)=xyzful00(iFish,iND,1:6)
    enddo
    xyzful(iFish,1:nND(iFish),1:6)=xyzful0(iFish,1:nND(iFish),1:6)

    allocate( du(nEQ(iFish)),lodEffe(nEQ(iFish)))
    lodEffe = (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 20000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -25000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
    
    call CG(lodEffe,jBC(iFish,1:nND(iFish),1:6),ele(iFish,1:nEL(iFish),1:5),xyzful0(iFish,1:nND(iFish),1),xyzful0(iFish,1:nND(iFish),2),xyzful0(iFish,1:nND(iFish),3),xyzful(iFish,1:nND(iFish),1),xyzful(iFish,1:nND(iFish),2),xyzful(iFish,1:nND(iFish),3), &
            prop(iFish,1:nMT(iFish),1:10),nND(iFish),nEL(iFish),nEQ(iFish),nMT(iFish)) 


end program StructureSolver


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
!    ELeMent STiFfness for FRaMe
!    calculates the element stiffness matrices.
!
!    ISBN 9787040258417 Zeng Pan. P70
!    ISBN 9780792312086 James F. Doyle. P81
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine elmstfFRM_D(length, ix, iy, iz,area, emod, gmod, ek ,nELt)
    implicit none
    real(8):: area, length, ix, iy, iz, emod, gmod
    real(8):: ek(12,12)
    integer:: nELt

    integer:: i,j
    real(8):: emlen,emlen2,emlen3
!
!   initialize all ek elements to zero
    ek(1:12,1:12)=0.0
!
!   STIFFNESS matrix in local coordinates
!   emod is Modulus of elasticity
    emlen  = emod/length
    emlen2 = emlen/length
    emlen3 = emlen2/length

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
!
!   impose the symmetry
    do  i= 1, 12
        do  j= i, 12
            ek(j,i) = ek(i,j)
        enddo
    enddo
!
    return
endsubroutine

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    ISBN 9780792312086 James F. Doyle. P83
!
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   TRANSformation of matrix in 3D. L->G
subroutine trans3d_D(l,m,n,ek,beta)
    implicit none
    real(8):: ek(12,12),rt(3,3),r(3,3),ktemp(12,12)
    real(8):: m,n,l,beta,pi,sb,cb,d
    integer:: i,j,k,j1,j2,ii,jj,in,jn

    pi=4.0*datan(1.0d0)
!
    sb=dsin(beta*pi/180)
    cb=dcos(beta*pi/180)
    d=dsqrt(1.0-n**2)
!   if (abs(l).ge. 0.995 .and. abs(beta).le. 0.01) return
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
              r(2,1)  =  -m/d
              r(2,2)  =  l/d
              r(2,3)  =  0.0
              r(3,1)  =  -l*n/d
              r(3,2)  =  -m*n/d
              r(3,3)  =  d
           else
              r(2,1)  =  -(m*cb+l*n*sb)/d
              r(2,2)  =  (l*cb-m*n*sb)/d
              r(2,3)  =  d*sb
              r(3,1)  =  (m*sb-l*n*cb)/d
              r(3,2)  =  -(l*sb+m*n*cb)/d
              r(3,3)  =  d*cb
           endif
    endif
!
    do  in=1,3
    do  jn=1,3
        rt(jn,in)=r(in,jn)
    enddo
    enddo
!    take [Rtrans][K][R] using the nature of [R] to speed computation.
!    k is sectioned off into 3x3s then multiplied [rtrans][k][r]
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
!    Conjugate Gradient Method
!    Iterative solution for displacement vector.
!
!    https://www.detailedpedia.com/wiki-Conjugate_gradient_method
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine CG(b,jBC,ele,xord0,yord0,zord0,xord,yord,zord,prop,nND,nEL,nEQ,nMT)
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
    
    
    !Initialize the solution vector x to 0.
    x = 0.0
    call CPU_TIME(start)
    !Initialize the residual vector. r = b - matmul(A, x)
    call matrixfree(jBC(1:nND,:),x,Ap,ele(1:nEL,1:5),xord0(1:nND),yord0(1:nND),zord0(1:nND),xord(1:nND),yord(1:nND),zord(1:nND), &
                        prop(1:nMT,:),nND,nEL,nEQ,nMT)    
    r=b-Ap

    !!!!!!!!!!!!
    p = r

    !Maximum iteration count.
    max_iter = 10000

    rsold = dot_product(r, r)

    do iter = 1, max_iter
        call matrixfree(jBC(1:nND,:),p,Ap,ele(1:nEL,1:5),xord0(1:nND),yord0(1:nND),zord0(1:nND),xord(1:nND),yord(1:nND),zord(1:nND), &
                        prop(1:nMT,:),nND,nEL,nEQ,nMT)    ! matrix free method, Calculate the matrix-vector product. Ap = matmul(A, p)
        alpha = rsold / dot_product(p, Ap)
        x = x + alpha * p
        r = r - alpha * Ap
        rsnew = dot_product(r, r)
        if(sqrt(rsnew) < 1e-10) exit
        beta = rsnew/rsold
        p = r + beta * p
        rsold = rsnew
    enddo

    write(*,*) 'end', iter
    
    call CPU_TIME(finish)  ! Get the end time.
    total_time = finish - start  ! Calculate the difference, namely the program running time.
    write(*,*) 'Program ran for ', total_time, ' seconds.'
    
endsubroutine

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Matrix Free Method
!    Implement matrix-vector multiplication for individual elements, without assembling the global matrix.
!
!    Element‐by‐Element(Hughes, Thomas J. R.(1983).doi:10.1061/(ASCE)0733-9399(1983)109:2(576))
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine matrixfree(jBC,x,b,ele,xord0,yord0,zord0,xord,yord,zord,prop,nND,nEL,nEQ,nMT)
    implicit none
    integer:: nND,nEL,nEQ,nMT
    integer:: ele(nEL,5)
    real(8):: prop(nMT,10)
    real(8):: xord0(nND), yord0(nND),zord0(nND)
    real(8):: xord(nND), yord(nND),zord(nND)
    real(8):: ek(18,18), ek12(12,12)

    integer:: i,j,k,n,i1,j1,k1,mat,nELt
    real(8):: e0,g0,a0,r0,b0,zix0,ziy0,ziz0,dx,dy,dz,xl0,xl,xll,xmm,xnn,xl9

    real(8):: b(nEQ),x(nEQ),x18(18),b18(18)
    integer:: idof(18),jBC(nND,6),iND,iEQ
    b=0.0
!   form each element matrix, and assemble
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
!           frame
!           material props not change
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
!
!           orientation
            dx= xord(j1) - xord(i1)
            dy= yord(j1) - yord(i1)
            dz= zord(j1) - zord(i1)
            xl=dsqrt(dx*dx+dy*dy+dz*dz)
            xll=dx/xl
            xmm=dy/xl
            xnn=dz/xl

!           Calculate the stiffness matrix and assemble
            xl9= xl0
            call elmstfFRM_D(xl9,zix0,ziy0,ziz0,a0,e0,g0,ek12, nELt)
!           use trans3d to transform from local to global stifness
            call trans3d_D(xll,xmm,xnn,ek12,b0)

!
!           expand to [18x18]
            ek(1:18,1:18)=0.0

            do    j=1,12
                do    k=1,12
                    ek(j,k)=ek12(j,k)
                enddo
            enddo
!           -------------------------------------------------------------------
!           apply boundary conditions
!           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
            do  iND=1, 3
                do  j=1, 6
                    if(jBC(ele(n,iND),j).gt.0)then
                        iEQ = (iND-1)*6+j
                        ek(iEQ,iEQ)=ek(iEQ,iEQ)*1.0d20
                        ! if(iter==0)then
                        !     lodEffe(iEQ)=stfEffe(nloc(iFish,iEQ))*vBC(iND,j)
                        ! else
                        !    lodEffe(iEQ)=stfElas(nloc(iFish,iEQ))*0.0d0
                        ! endif
                    endif
                enddo
            enddo

!           -------------------------------------------------------------------
!           Multiply the element stiffness matrix by the displacement vector, and assemble the force vector
!           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
            do i = 1, 12
                x18(i)=x(idof(i))
            enddo
            b18=matmul(ek,x18)
            do i = 1, 12
                b(idof(i))=b(idof(i))+b18(i)
            enddo
!
        else
            write(*,*)'not this nELt:',nELt
            stop
        endif
    enddo

    return    
endsubroutine