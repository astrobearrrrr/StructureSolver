program StructureSolver

    integer, allocatable:: nND(:),nEL(:),nEQ(:),nMT(:),nBD(:),nSTF(:)
    integer, allocatable:: ele(:,:,:),nloc(:,:),nprof(:,:),nprof2(:,:),jBC(:,:,:)
    real(8), allocatable:: xyzful00(:,:,:),mssful(:,:,:),vBC(:,:,:),mss(:,:),prop(:,:,:)
    real(8), allocatable:: lodful(:,:,:),extful(:,:,:),extful1(:,:,:),extful2(:,:,:),grav(:,:,:),repful(:,:,:)
    real(8), allocatable:: xyzful0(:,:,:),xyzfulnxt(:,:,:),dspful(:,:,:),accful(:,:,:)
    real(8), allocatable:: xyzful(:,:,:),xyzfulIB(:,:,:),velful(:,:,:)
    real(8), allocatable:: stfEffe(:),stfElas(:),stfGeom(:)
    real(8), allocatable:: du(:),lodEffe(:)
    integer :: idat, iFish, nFish, nND_max, nEL_max, nMT_max, ierror
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
    call readdt(jBC(iFish,1:nND(iFish),1:6),ele(iFish,1:nEL(iFish),1:5),nloc(iFish,1:nND(iFish)*6),nprof(iFish,1:nND(iFish)*6), &
                nprof2(iFish,1:nND(iFish)*6),xyzful00(iFish,1:nND(iFish),1:6),prop(iFish,1:nMT(iFish),1:10), &
                nND(iFish),nEL(iFish),nEQ(iFish),nMT(iFish),nBD(iFish),nSTF(iFish),idat)
    close(idat)


    do  iND=1,nND(iFish)
        xyzful0(iFish,iND,1:6)=xyzful00(iFish,iND,1:6)
    enddo
    xyzful(iFish,1:nND(iFish),1:6)=xyzful0(iFish,1:nND(iFish),1:6)

    allocate( stfEffe(nSTF(iFish)),stfElas(nSTF(iFish)),stfGeom(nSTF(iFish)))

    call formstif_s(stfElas,ele(iFish,1:nEL(iFish),1:5),xyzful0(iFish,1:nND(iFish),1),xyzful0(iFish,1:nND(iFish),2),xyzful0(iFish,1:nND(iFish),3),xyzful(iFish,1:nND(iFish),1),xyzful(iFish,1:nND(iFish),2),xyzful(iFish,1:nND(iFish),3), &
                        prop(iFish,1:nMT(iFish),1:10),nprof(iFish,1:nND(iFish)*6),nloc(iFish,1:nND(iFish)*6),triad_e0,triad_ee,nND(iFish),nEL(iFish),nEQ(iFish),nMT(iFish),nSTF(iFish))


    ! open(unit=11, file = '1.dat')
    ! do i=1,20
    !     do j=1,nSTF(iFish)
    !         if(abs(stfElas(j)) .lt. 0.1)then
    !             stfElas(j)=0.0
    !         endif
    !     enddo
    !     write(11,'(24E28.5)') stfElas(i+i*(i-1)/2-i+1:i+i*(i-1)/2)/29.5/10000/100*6
    ! enddo
    ! close(11)

    allocate( du(nEQ(iFish)),lodEffe(nEQ(iFish)))
    lodEffe = (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 20000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -25000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
!       -------------------------------------------------------------------
!       apply boundary conditions
!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do  iND=1,nND(iFish)
        do  j=1  ,6
            if(jBC(iFish,iND,j)>0)then
                iEQ=(iND-1)*6+j
                stfElas(nloc(iFish,iEQ))=stfElas(nloc(iFish,iEQ))*1.0d20
                ! if(iter==0)then
                !     lodEffe(iEQ)=stfEffe(nloc(iFish,iEQ))*vBC(iND,j)
                ! else
                !    lodEffe(iEQ)=stfElas(nloc(iFish,iEQ))*0.0d0
                ! endif
            endif
        enddo
        enddo

!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       -------------------------------------------------------------------
        call CG(stfElas,nSTF(iFish),lodEffe,nEQ(iFish),nBD(iFish),du,nprof(iFish,:),nprof2(iFish,:),nloc(iFish,:))

end program StructureSolver


!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    copyright@ RuNanHua
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   READ structural DaTafile
subroutine readdt(jBC,ele,nloc,nprof,nprof2,xyzful0,prop,nND,nEL,nEQ,nMT,nBD,nSTF,idat)
    implicit none
    integer:: nND,nEL,nEQ,nMT,nBD,nSTF,idat
    integer:: ele(nEL,5),jBC(nND,6),nloc(nND*6),nprof(nND*6),nprof2(nND*6)
    real(8):: xyzful0(nND,6),prop(nMT,10)
!   ---------------------------------------------------------------------------
    integer:: i,j,nm1,nm2,material,nbc,node,nmp,ii,ind,ibandh,iend,ibandv,ji1
    character*50 title,endin
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
    read(idat,*) nMT
    do    i= 1, nMT
        read(idat,*) nmp,prop(nmp,1:8)
    enddo
    read(idat,'(1a50)') endin
!   -----------------------------------------------------------------------------------------------

    nEQ=nND*6

    call maxbnd(ele,nprof,nND,nEL,nEQ,nBD)

!   nprof2
    do i=1,nEQ
       ibandh=1
       iend=i+nBD-1
       if (iend .gt. nEQ) iend=nEQ
       do j=i+1,iend
            ibandv=nprof(j)
            ji1=j-i+1
            if  (ibandv .ge. ji1) then
                ibandh = ji1
            endif
       enddo
       nprof2(i)=ibandh
    enddo

    nloc(1)=1
    do    i=1,nEQ-1
        nloc(i+1) = nloc(i) + nprof(i)
    enddo
    nSTF=nloc(nEQ)+nprof(nEQ)
    !nprof  :  1  2  3  4  5  6  7  8  9 10 11 12  7  8  9 10 11 12  7  8  9 ...
    !nprof2 : 12 11 10  9  8  7 12 11 10  9  8  7 12 11 10  9  8  7 12 11 10 ...

    !write(*,'(3(A,1x,I8,2x))')'nND=',nND,'nEL=',nEL,'nEQ=',nEQ
    !write(*,'(3(A,1x,I8,2x))')'nMT=',nMT,'nBD=',nBD,'nSTF=',nSTF
    return
endsubroutine

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   MAX Band Width calculation
subroutine maxbnd(ele,nprof,nND,nEL,nEQ,nBD)
    implicit none
    integer:: nBD,nEQ,nEL,nND
    integer:: ele(nEL,5),nprof(nEQ)

    integer:: ipv(18)
    integer:: i,j,n,ihfbnd,idof,jdof,kdof,ieqn1,ieqn2,jband,ieq,ieq2

    nprof(1:nEQ)=0

    ihfbnd=0
    do  n=1,nEL
        idof=(ele(n,1)-1)*6
        jdof=(ele(n,2)-1)*6
        kdof=(ele(n,3)-1)*6
        do  i=1,6
            ipv(i   )=idof+i
            ipv(i+6 )=jdof+i
            ipv(i+12)=kdof+i
        enddo

        do  i=1,18
            ieqn1=ipv(i)
            do  j=i,18
                ieqn2=ipv(j)
                ihfbnd = max0(ihfbnd,iabs(ieqn1-ieqn2))
                jband=abs(ieqn1-ieqn2)+1
                ieq=max(ieqn1,ieqn2)
                if  (jband .gt. nprof(ieq)) then
                    nprof(ieq)=jband
                endif
                ieq2=min(ieqn1,ieqn2)
!               if  (jband .gt. nprof2(ieq2)) then
!                   nprof2(ieq2)=jband
!               endif
            enddo
        enddo

    enddo
    nBD=ihfbnd+1

    return
endsubroutine

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    FORM STIFfness matrix  [K]
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine formstif_s(stf,ele,xord0,yord0,zord0,xord,yord,zord,prop,nprof,nloc,triad_e0,triad_ee,nND,nEL,nEQ,nMT,nSTF)
    implicit none
    integer:: nND,nEL,nEQ,nMT,nSTF
    integer:: ele(nEL,5),nprof(nEQ), nloc(nEQ)
    real(8):: prop(nMT,10)
    real(8):: xord0(nND), yord0(nND),zord0(nND)
    real(8):: xord(nND), yord(nND),zord(nND)
    real(8):: stf(nSTF)
    real(8):: ek(18,18), ekb(18,18),ek12(12,12),ek9(9,9)

    real(8):: triad_e0(3,3,nEL),triad_ee(3,3,nEL),rr(3,3)

    real(8):: xyz12(3),xyz13(3),xyzb12(3),xyzb13(3)
    real(8):: x1,x2,x3,y1,y2,y3,z1,z2,z3
    real(8):: xb1,xb2,xb3,yb1,yb2,yb3

    integer:: i,j,k,ipress,n,i1,j1,k1,mat,nELt,ii,jj
    real(8):: e0,g0,a0,r0,b0,zix0,ziy0,ziz0,dx,dy,dz,xl0,xl,xll,xmm,xnn,xl9,area,t0,zip0,zia0,zib0,pl0,alpha,beta
    real(8):: temp

!   initialize [K]  to zero
    stf(1:nSTF)=0.0

!   form each element matrix, and assemble
    do    n=1,nEL
        i1  = ele(n,1)
        j1  = ele(n,2)
        k1  = ele(n,3)
        nELt= ele(n,4)
        mat = ele(n,5)
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

            call assembCOL(nSTF,nND,nEQ,stf,ek,i1,j1,k1,nloc)
!
        else
            write(*,*)'not this nELt:',nELt
            stop
        endif
    enddo

    return
endsubroutine

!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    ELeMent STiFfness for FRaMe
!    calculates the element stiffness matrices.
!    copyright@ RuNanHua
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
!    copyright@ RuNanHua
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
!    ASSEMBle element matrices in COLumn form
!    copyright@ RuNanHua
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine assembCOL(nSTF,nND,nEQ,aa,a,i1,j1,k1,nloc)
    implicit none
    integer:: nSTF,nND,nEQ
    real(8):: aa(nSTF),a(18,18)
    integer:: idof(18)
    integer:: i1,j1,k1,nloc(nEQ)
    integer:: i,j,imax,jmax,ieqn1,ieqn2,jband,iloc
    imax = 0
!
!   Set idof to posible DoF number of each nodes
    if    (j1 .gt. 0  .AND. k1 .gt. 0) then
        imax=18
        jmax=18
    elseif (j1 .eq. 0  .AND. k1 .eq. 0) then
        imax=6
        jmax=6
    endif

    do  i= 1, 6
        idof(i)    = (i1-1)*6 + i
        idof(i+6)  = (j1-1)*6 + i
        idof(i+12) = (k1-1)*6 + i
    enddo
!
!   Store the values for individual array in global array
    do  i= 1, imax
        ieqn1 = idof(i)
        do    j= i, jmax
            ieqn2 = idof(j)
            if    (ieqn1 .gt. ieqn2) then
                    jband= (ieqn1-ieqn2)+1
                    iloc = nloc(ieqn1)
                    aa(iloc +jband-1) = aa(iloc +jband-1) + a(i,j)
!                   aa(ieqn2,jband) = aa(ieqn2,jband) + a(i,j)
            else
                    jband= (ieqn2-ieqn1)+1
                    iloc = nloc(ieqn2)
                    aa(iloc +jband-1) = aa(iloc +jband-1) + a(i,j)
!                   aa(ieqn1,jband) = aa(ieqn1,jband) + a(i,j)
            endif
        enddo
    enddo
!
    return
endsubroutine


subroutine CG(a,maxstore,b,nEQ,nBD,wk,nprof,nprof2,nloc)
    implicit none
    integer:: maxstore,nEQ,nBD,imult
    real(8):: a(maxstore),b(nEQ),wk(nEQ)
    integer:: nprof(nEQ), nprof2(nEQ),nloc(nEQ)

    real(8):: AA(nEQ,nEQ),x(nEQ),r(nEQ), p(nEQ), Ap(nEQ), z(nEQ)
    integer :: i, j, n, iter, max_iter, k
    double precision :: alpha, beta, rsold, rsnew
    write(*,*)nloc
    AA=0.0
    do i = 1,nEQ
        do j = 1,i
            AA(i,i-j+1) = a(nloc(i)+j-1)
            AA(i-j+1,i) = AA(i,i-j+1)
        enddo
    enddo
    !初始化解向量x为0
    x = 0.0

    !初始化残差向量
    r = b - matmul(AA, x)
    !!!!!!!!!!!!
    p = r

    !最大迭代次数
    max_iter = 10000

    rsold = dot_product(r, r)

    do iter = 1, max_iter
      Ap = matmul(AA, p)
      alpha = rsold / dot_product(p, Ap)
      x = x + alpha * p
      r = r - alpha * Ap
      rsnew = dot_product(r, r)
      if(sqrt(rsnew) < 1e-10) exit
      p = r + (rsnew/rsold) * p
      rsold = rsnew
    end do
    write(*,*) 'end', iter

    open(unit=11, file = '3.dat')
    do  i = 1,nEQ
        write(11,'(24E28.5)') AA(i,:)/29.5/10000/100*6
    enddo
    write(11,'(A      )') 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
        write(11,'(24E28.5)') x
    close(11)

    
    end