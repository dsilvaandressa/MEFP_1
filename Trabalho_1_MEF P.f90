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

    integer(8) :: nnos, nel, nt, no, k, i, k1, k2, dir, j, nnoscar, n, nnres, iter, passo, tipo
    integer(8), allocatable, dimension(:) :: ip
    double precision :: tol, l0, lf, aux, K_svk, err, normx
    double precision, allocatable, dimension(:,:) :: HEX, HEXg
    double precision, allocatable, dimension(:,:) :: x, y, prop, dy
    double precision, allocatable, dimension(:) :: E, s_svk, ue_svk, fn_svk, fint_svk, Fnodais, g, Rnodais, res, xaux, Fnodais_inc, Rnodais_inc
    integer, allocatable, dimension (:,:) :: inc
    double precision, allocatable, dimension(:,:,:) :: u
    character(len=5):: texto
    character(len=3):: ext, tostr
    
    open (15,FILE="entrada_cupula.txt",STATUS="OLD")
    read(15,*) !Num de nós, Num de elem, Num passos, tolerância, Num de forças aplicadas, Num de nós restritos, Tipo de Anlálise (1 - estat, 2 - dinam)
    read(15,*) nnos, nel, nt, tol, nnoscar, nnres, tipo
    read(15,*) !No, X1, X2, X3

    allocate (x(nnos,3), y(nnos,3), prop(nel,2), dy(nnos,3), u(nt,nnos,3))
    allocate (inc(nel,2))
    allocate (s_svk(nel), E(nel), ue_svk(nel), fn_svk(nel), fint_svk(3*nnos))
    allocate (HEX(6,6), HEXg(3*nnos,3*nnos))
    allocate (Fnodais(3*nnos), g(3*nnos), ip(3*nnos), xaux(3*nnos), Fnodais_inc(3*nnos))
    allocate (Rnodais(3*nnos), res(3*nnos), Rnodais_inc(3*nnos))
    
    HEXg=0
    HEX=0
    
    !Leitura da entrada
    !x(1, 1-3) coordenadas do primeiro nó
    !inc(1, 1-2) nós do primeiro elemento
    !prop(1, 1-2) elasticidade e área respectivamente do primeiro elemento

    i=0
    do i=1,nnos
        read(15,*) no, x(no,1), x(no,2), x(no,3)
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

    read(15,*) !No, Direção, Força
    i=0
    Fnodais = 0
    if (nnoscar>0) then
        do i=1,nnoscar
            read(15,*) k, dir, Fnodais(3*k-3+dir)
        end do
    else
    end if
    
    read(15,*) !No, R01 R02 R03
    i=0
    res=0
    do i=1,nnres
        read(15,*) no, res(no*3-2), res(no*3-1), res(no*3)    
    end do
    
    Rnodais = 0
    
    read(15,*) texto !No, dY01 dY02 dY03
    i=0
    do i=1,nnres
        read(15,*) no, Rnodais(no*3-2), Rnodais(no*3-1), Rnodais(no*3)      
    end do
    
        
    close(15)
   

    
    !Passos de carga
    
    !Tentativa inicial 
    
        y=x
        do i=1,nnos
            xaux(3*i-2) = x(i,1)
            xaux(3*i-1) = x(i,2)
            xaux(3*i) = x(i,3)
        end do
    
        normx = abs(norm2(xaux))
          
        Fnodais_inc = Fnodais/nt !Inicial
        Fnodais = 0
        Rnodais_inc = Rnodais/nt
    
    
    do passo=1,nt
    
        Fnodais = Fnodais + Fnodais_inc
        
        do i=1,nnos*3
            if (Rnodais(i)/=0) then
                if (mod(i,3)==0) then
                    y(int(i/3),3) = y(int(i/3),3) + Rnodais_inc(i)
                else
                    y(int(i/3)+1,mod(i,3)) = y(int(i/3)+1,mod(i,3)) + Rnodais_inc(i)
                end if
            end if
        end do
               
        i=0
        
    !Início Newton-Raphson
        iter=0
        do 
            !Cálculo da Matriz Hessiana
            fint_svk = 0
            HEX=0 
            HEXg=0
            do k=1,nel
        
                k1=inc(k,1)
                k2=inc(k,2)
        
                l0=((x(k1,1)-x(k2,1))*(x(k1,1)-x(k2,1)) + (x(k1,2)-x(k2,2))*(x(k1,2)-x(k2,2)) + (x(k1,3)-x(k2,3))*(x(k1,3)-x(k2,3)))**0.5
                lf=((y(k1,1)-y(k2,1))*(y(k1,1)-y(k2,1)) + (y(k1,2)-y(k2,2))*(y(k1,2)-y(k2,2)) + (y(k1,3)-y(k2,3))*(y(k1,3)-y(k2,3)))**0.5
                
                K_svk = prop(k,1)
                s_svk(k) = (((lf*lf)/(l0*l0))-1)*0.5*prop(k,1)
        
                aux = (prop(k,2)/l0)   
            
                !Forças internas 
                    fint_svk(k1*3-2) =  fint_svk(k1*3-2) + prop(k,2)*s_svk(k)*(-1)*(y(k2,1)-y(k1,1))/l0
                    fint_svk(k1*3-1) = fint_svk(k1*3-1) + prop(k,2)*s_svk(k)*(-1)*(y(k2,2)-y(k1,2))/l0
                    fint_svk(k1*3) = fint_svk(k1*3) + prop(k,2)*s_svk(k)*(-1)*(y(k2,3)-y(k1,3))/l0
            
                    fint_svk(k2*3-2) = fint_svk(k2*3-2) + prop(k,2)*s_svk(k)*(1)*(y(k2,1)-y(k1,1))/l0
                    fint_svk(k2*3-1) = fint_svk(k2*3-1) + prop(k,2)*s_svk(k)*(1)*(y(k2,2)-y(k1,2))/l0
                    fint_svk(k2*3) = fint_svk(k2*3) + prop(k,2)*s_svk(k)*(1)*(y(k2,3)-y(k1,3))/l0
                   
                
                do i=1,3
                        HEX(i,1) = aux*(K_svk*(y(k2,i)-y(k1,i))*(y(k2,1)-y(k1,1)))/(l0*l0) 
                        HEX(i,2) = aux*(K_svk*(y(k2,i)-y(k1,i))*(y(k2,2)-y(k1,2)))/(l0*l0)
                        HEX(i,3) = aux*(K_svk*(y(k2,i)-y(k1,i))*(y(k2,3)-y(k1,3)))/(l0*l0)
                        HEX(i,4) = -aux*(K_svk*(y(k2,i)-y(k1,i))*(y(k2,1)-y(k1,1)))/(l0*l0)
                        HEX(i,5) = -aux*(K_svk*(y(k2,i)-y(k1,i))*(y(k2,2)-y(k1,2)))/(l0*l0)
                        HEX(i,6) = -aux*(K_svk*(y(k2,i)-y(k1,i))*(y(k2,3)-y(k1,3)))/(l0*l0)                
                end do    
        
                do i=4,6
                        HEX(i,1) = -aux*(K_svk*(y(k2,i-3)-y(k1,i-3))*(y(k2,1)-y(k1,1)))/(l0*l0) 
                        HEX(i,2) = -aux*(K_svk*(y(k2,i-3)-y(k1,i-3))*(y(k2,2)-y(k1,2)))/(l0*l0)
                        HEX(i,3) = -aux*(K_svk*(y(k2,i-3)-y(k1,i-3))*(y(k2,3)-y(k1,3)))/(l0*l0)
                        HEX(i,4) = aux*(K_svk*(y(k2,i-3)-y(k1,i-3))*(y(k2,1)-y(k1,1)))/(l0*l0)
                        HEX(i,5) = aux*(K_svk*(y(k2,i-3)-y(k1,i-3))*(y(k2,2)-y(k1,2)))/(l0*l0)
                        HEX(i,6) = aux*(K_svk*(y(k2,i-3)-y(k1,i-3))*(y(k2,3)-y(k1,3)))/(l0*l0)                
                end do   
        
                HEX(1,1) = HEX(1,1) + aux*s_svk(k)
                HEX(2,2) = HEX(2,2) + aux*s_svk(k)
                HEX(3,3) = HEX(3,3) + aux*s_svk(k)
                HEX(4,4) = HEX(4,4) + aux*s_svk(k)
                HEX(5,5) = HEX(5,5) + aux*s_svk(k)
                HEX(6,6) = HEX(6,6) + aux*s_svk(k)
        
                HEX(1,4) = HEX(1,4) - aux*s_svk(k)
                HEX(2,5) = HEX(2,5) - aux*s_svk(k)
                HEX(3,6) = HEX(3,6) - aux*s_svk(k)
                HEX(4,1) = HEX(4,1) - aux*s_svk(k)
                HEX(5,2) = HEX(5,2) - aux*s_svk(k)
                HEX(6,3) = HEX(6,3) - aux*s_svk(k)
            
            
                do i=1,3
                    do j=1,3
                        HEXg(3*k1-3+i,3*k1-3+j) =  HEXg(3*k1-3+i,3*k1-3+j)  + HEX(i,j)
                        HEXg(3*k1-3+i,3*k2-3+j) =  HEXg(3*k1-3+i,3*k2-3+j)  + HEX(i,j+3)                
                        HEXg(3*k2-3+i,3*k1-3+j) =  HEXg(3*k2-3+i,3*k1-3+j)  + HEX(i+3,j)
                        HEXg(3*k2-3+i,3*k2-3+j) =  HEXg(3*k2-3+i,3*k2-3+j)  + HEX(i+3,j+3)                
                    end do
                end do
            
                !HEXg(3*k1-2:3*k1, 3*k1-2:3*k1) = HEXg(3*k1-2:3*k1, 3*k1-2:3*k1) + HEX(1:3,1:3)
                !HEXg(3*k1-2:3*k1, 3*k2-2:3*k2) = HEXg(3*k1-2:3*k1, 3*k2-2:3*k2) + HEX(1:3,4:6)
                !HEXg(3*k2-2:3*k2, 3*k1-2:3*k1) = HEXg(3*k2-2:3*k2, 3*k1-2:3*k1) + HEX(4:6,1:3)
                !HEXg(3*k2-2:3*k2, 3*k2-2:3*k2) = HEXg(3*k2-2:3*k2, 3*k2-2:3*k2) + HEX(4:6,4:6)
                !
            end do
            !Condições de contorno
            g = Fnodais - fint_svk
            do i=1,3*nnos
                if (res(i)==1) then
                    HEXg(i,:) = 0
                    HEXg(:,i) = 0
                    HEXg(i,i) = 1
                    g(i)=0
                else
                end if
            end do
        
            ! Solução sistema        
            n = 3*nnos
            ip=0
            i=0
            j=1
    
            call dgesv(n,j,HEXg,n,ip,g,n,i)
        
            
            do i=1,nnos
                dy(i,1) = g(3*i-2)
                dy(i,2) = g(3*i-1)
                dy(i,3) = g(3*i)
            end do
    
            y = y+dy
            err = abs(norm2(g))/normx
        
            iter = iter + 1
        
            print*, iter
        
            if (err<=tol) then 
                exit
            else if (iter>1000) then
                print*, "Não convergiu"
                exit
            end if        
        end do ! Newton-Raphson
    if (iter>1000) then
        pause
        exit
    end if
    
        u(passo,:,:) = y - x

    end do !passos de carga
    
      print*, u(1,13,3)
      
    ! Escrita de dados para visualização
    
          
      open (25,FILE="saida0.vtu",STATUS="REPLACE")
    
1     format(A,I0,A,I0,A)
2     format(F0.6,1x,F0.6,1x,F0.6)     
3     format(I0,1x,I0) 
4     format(I0) 
5     format(A) 
      
     write (25,5) "<?xml version=""1.0""?>"
     write (25,5) "<VTKFile type=""UnstructuredGrid"">"
     write (25,5) "  <UnstructuredGrid>"
     write (25,1) "  <Piece NumberOfPoints=""",nnos,"""  NumberOfCells=""",nel,""">"
     write (25,5) "    <Points>"
     write (25,5) "      <DataArray type=""Float64"" NumberOfComponents=""3"" format=""ascii"">"
     
     do i=1,nnos
         write (25,2) x(i,1), x(i,2), x(i,3)
     end do
     
     write (25,5) "      </DataArray>"
     write (25,5) "    </Points>"
     write (25,5) "    <Cells>"
     write (25,5) "      <DataArray type=""Int32"" Name=""connectivity"" format=""ascii"">"
     
     do i=1,nel
         write(25,3) inc(i,1)-1, inc(i,2)-1
     end do
          
     write (25,5) "      </DataArray>"
     write (25,5) "      <DataArray type=""Int32"" Name=""offsets"" format=""ascii"">"
     
     do i=1,nel
         write(25,4) 2*i
     end do
     
     write (25,5) "      </DataArray>"
     write (25,5) "      <DataArray type=""UInt8"" Name=""types"" format=""ascii"">"
     
     do i=1,nel
         write(25,4) 4
     end do
     
     write (25,5) "      </DataArray>"
     write (25,5) "    </Cells>"
     write (25,5) "    <PointData>"
     write (25,5) "      <DataArray type=""Float64"" NumberOfComponents=""3"" Name=""Displacement"" format=""ascii"">"
     
     do i=1,nnos
         write (25,2) u(nt,i,1), u(nt,i,2), u(nt,i,3)
     end do
     
     write (25,5) "      </DataArray>"
     write (25,5) "    </PointData>"
     write (25,5) "    <CellData>"
     write (25,5) "    </CellData>"
     write (25,5) "  </Piece>"
     write (25,5) "  </UnstructuredGrid>"
     write (25,5) "</VTKFile>"
     
     close(25)
     
     do j=1,nt
        
        ext = tostr(j)
     
         open (35,FILE="saida"//trim(ext)//".vtu",STATUS="REPLACE")
     
         write (35,5) "<?xml version=""1.0""?>"
         write (35,5) "<VTKFile type=""UnstructuredGrid"">"
         write (35,5) "  <UnstructuredGrid>"
         write (35,1) "  <Piece NumberOfPoints=""",nnos,"""  NumberOfCells=""",nel,""">"
         write (35,5) "    <Points>"
         write (35,5) "      <DataArray type=""Float64"" NumberOfComponents=""3"" format=""ascii"">"
     
         do i=1,nnos
             write (35,2) (x(i,1)+u(j,i,1)), (x(i,2)+u(j,i,2)), (x(i,3)+u(j,i,3))
         end do
     
         write (35,5) "      </DataArray>"
         write (35,5) "    </Points>"
         write (35,5) "    <Cells>"
         write (35,5) "      <DataArray type=""Int32"" Name=""connectivity"" format=""ascii"">"
     
         do i=1,nel
             write(35,3) inc(i,1)-1, inc(i,2)-1
         end do
          
         write (35,5) "      </DataArray>"
         write (35,5) "      <DataArray type=""Int32"" Name=""offsets"" format=""ascii"">"
     
         do i=1,nel
             write(35,4) 2*i
         end do
     
         write (35,5) "      </DataArray>"
         write (35,5) "      <DataArray type=""UInt8"" Name=""types"" format=""ascii"">"
     
         do i=1,nel
             write(35,4) 4
         end do
     
         write (35,5) "      </DataArray>"
         write (35,5) "    </Cells>"
         write (35,5) "    <PointData>"
         write (35,5) "      <DataArray type=""Float64"" NumberOfComponents=""3"" Name=""Displacement"" format=""ascii"">"
     
         do i=1,nnos
             write (35,2) u(j,i,1), u(j,i,2), u(j,i,3)
         end do
     
         write (35,5) "      </DataArray>"
         write (35,5) "    </PointData>"
         write (35,5) "    <CellData>"
         write (35,5) "    </CellData>"
         write (35,5) "  </Piece>"
         write (35,5) "  </UnstructuredGrid>"
         write (35,5) "</VTKFile>"
     
         close(35)
     
     end do
     
    end program Trabalho1MEFP

function tostr(i) result(ext)
    implicit none
    character(len=15) :: ext
    integer, intent(in) :: i
    character(range(i)+2) :: tmp
    write(tmp,'(i0)') i
    ext = trim(tmp)
end function