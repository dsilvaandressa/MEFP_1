!  Trabalho1MEFP.f90 
!
!  FUNCTIONS:
!  Trabalho1MEFP - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: Trabalho1MEFP
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program Trabalho1MEFP

    implicit none

    integer :: nnos, nel, nt, no, k, i, k1, k2, beta, dir, j
    double precision :: tol, l0, lf
    double precision, allocatable, dimension(:,:,:) :: HEX
    double precision, allocatable, dimension(:,:) :: x, y, prop
    double precision, allocatable, dimension(:) :: E, s_svk, ue_svk, fn_svk, fint_svk
    integer, allocatable, dimension (:,:) :: inc

    open (15,FILE="entrada.txt",STATUS="OLD")
    read(15,*) !Num de nós, Num de elem, Num passos, tolerância
    read(15,*) nnos, nel, nt, tol
    read(15,*) !No, X1, X2, X3

    allocate (x(nnos,3), y(nnos,3), prop(nel,2))
    allocate (inc(nel,2))
    allocate (s_svk(nel), E(nel), ue_svk(nel), fn_svk(nel), fint_svk(6*nel))
    allocate (HEX(nel,nnos*3,nnos*3))

    !Leitura da entrada
    !x(1, 1-3) coordenadas do primeiro nó
    !inc(1, 1-2) nós do primeiro elemento
    !prop(1, 1-2) elasticidade e área respectivamente do primeiro elemento

    open (20,FILE="pos_final.txt",STATUS="OLD")
    read (20,*) !No, X1, X2, X3
    i=0
    do i=1,nnos
        read(15,*) no, x(no,1), x(no,2), x(no,3)
        read(20,*) no, y(no,1), y(no,2), y(no,3)
!        y(no,1) = x(no,1) + 1
!        y(no,2) = x(no,2) + 2
!        y(no,3) = x(no,3) - 1
    end do

    read(15,*) !EL, no1, no2
    i=0
    do i=1,nel
        read(15,*) k, inc(k,1), inc(k,2)
    end do

    read(15,*) !EL, Elasticidade, Área
    i=0
    do i=1,nel
        read(15,*) k, prop(k,1), prop(k,2)
    end do

    close(15)
    close(20)

    !Comprimento incial e final dos elementos e forças internas
    i=0
    do i=1,nel
		
		k1=inc(i,1)
        k2=inc(i,2)
		
        l0=((x(k1,1)-x(k2,1))**2 + (x(k1,2)-x(k2,2))**2 + (x(k1,3)-x(k2,3))**2)**0.5
        lf=((y(k1,1)-y(k2,1))**2 + (y(k1,2)-y(k2,2))**2 + (y(k1,3)-y(k2,3))**2)**0.5
        E(i)=((lf**2/l0**2)-1)/2
        s_svk(i) = E(i)*prop(i,1)
        ue_svk(i) = (E(i)**2)*prop(i,1)/2
        fn_svk(i) = s_svk(i)*prop(i,2)*lf/l0

        do k=1,6
            if (k<=3) then
                beta = 1
                dir = k
            else
                beta = 2
                dir = k-3
            end if
            !Cálculo das forças internas
            fint_svk(i+k-1) = prop(i,2)*s_svk(i)*(-1**beta)*(y(k2,dir)-y(k1,dir))/l0
        end do

    end do

    !Cálculo da Matriz Hessiana
     
    do k=1,nel
        
        k1=inc(k,1)
        k2=inc(k,2)
        
        l0=((x(k1,1)-x(k2,1))**2 + (x(k1,2)-x(k2,2))**2 + (x(k1,3)-x(k2,3))**2)**0.5
        lf=((y(k1,1)-y(k2,1))**2 + (y(k1,2)-y(k2,2))**2 + (y(k1,3)-y(k2,3))**2)**0.5
        
        do i=1,3*nnos
            do j=1,3*nnos
                HEX(k,i,j) = (prop(k,2)/l0)
            end do        
        end do        
    end do
    

    
    end program Trabalho1MEFP

