MODULE COMMON
  USE INPUT_VARS
  implicit none

  complex(8),parameter             :: zero=(0d0,0d0),one=(1d0,0d0),xi=(0d0,1d0)
  real(8),parameter                :: pi    = 3.14159265358979323846264338327950288419716939937510d0
  complex(8),allocatable           :: Hij(:,:,:)      ![Nlso][Nlso][Nk=1]
  complex(8),allocatable           :: Hloc(:,:)   ![Nlso][Nlso]
  integer                          :: Nlat,Nineq,Norb
  integer,dimension(:),allocatable :: lat2ineq,ineq2lat
  real(8),dimension(:),allocatable :: sb_field_sign



contains





  !----------------------------------------------------------------------------------------!
  ! purpose: read the real local self-energy from disk
  !----------------------------------------------------------------------------------------!
  subroutine read_sigma(Sreal)
    complex(8),intent(inout) :: Sreal(Nineq,Nspin,Nspin,Norb,Norb,Lreal)
    character(len=30)        :: suffix
    integer                  :: i,iineq,ispin,iorb,unit
    real(8)                  :: wr(Lreal),reS,imS
    !
    wr = linspace(wini,wfin,Lreal)
    !
    write(*,*)"write spin-orbital diagonal elements:"
    do ispin=1,Nspin
       do iorb=1,Norb
          write(suffix,"(A2,I1,A2,I1,A10)")"_l",iorb,"_s",ispin,"_realw.dat"
          open(newunit=unit,file=trim(Sfile)//trim(suffix),status='old')
          do iineq=1,Nineq
             do i=1,Lreal
                read(unit,*)wr(i),imS,reS
                Sreal(iineq,ispin,ispin,iorb,iorb,i) = dcmplx(reS,imS)
             enddo
          enddo
          close(unit)
       enddo
    enddo
    print*,"ciao"
  end subroutine read_sigma





  !-----------------------------------------------------------------------------
  ! Purpose:
  !-----------------------------------------------------------------------------
  function get_file_length(file) result(lines)
    integer           :: lines
    character(len=*)  :: file
    integer           :: ifile,ierr,pos
    logical           :: IOfile,bool,bool1,bool2
    character(len=256)::buffer
    inquire(file=trim(file),exist=IOfile)
    lines=0
    if(.not.IOfile)then
       write(*,*) 'Cannot read +'//trim(file)//'. Skip file_size'
       return
    endif
    open(99,file=trim(file))
    ierr=0
    do while(ierr==0)
       lines=lines+1
       read(99,*,iostat=ierr)buffer
       bool1=scan(buffer,"#").ne.0
       bool2=len_trim(buffer).eq.0       
       if(bool1 .OR. bool2)lines=lines-1
    enddo
    lines=lines-1
    rewind(99)
    close(99)
  end function get_file_length




  !-----------------------------------------------------------------------------
  ! Purpose:
  !-----------------------------------------------------------------------------
  function linspace(start,stop,num,istart,iend,mesh) result(array)
    real(8)          :: start,stop,step,array(num)
    integer          :: num,i
    logical,optional :: istart,iend
    logical          :: startpoint_,endpoint_
    real(8),optional :: mesh
    if(num<0)stop "linspace: N<0, abort."
    startpoint_=.true.;if(present(istart))startpoint_=istart
    endpoint_=.true.;if(present(iend))endpoint_=iend
    if(startpoint_.AND.endpoint_)then
       if(num<2)stop "linspace: N<2 with both start and end points"
       step = (stop-start)/real(num-1,8)
       forall(i=1:num)array(i)=start + real(i-1,8)*step
    elseif(startpoint_.AND.(.not.endpoint_))then
       step = (stop-start)/real(num,8)
       forall(i=1:num)array(i)=start + real(i-1,8)*step
    elseif(.not.startpoint_.AND.endpoint_)then
       step = (stop-start)/real(num,8)
       forall(i=1:num)array(i)=start + real(i,8)*step
    else
       step = (stop-start)/real(num+1,8)
       forall(i=1:num)array(i)=start + real(i,8)*step
    endif
    if(present(mesh))mesh=step
  end function linspace





  !-----------------------------------------------------------------------------
  ! Purpose:
  !-----------------------------------------------------------------------------
  function blocks_to_matrix(Vblocks,Nlat,Nso) result(Matrix)
    complex(8),dimension(Nlat,Nso,Nso)      :: Vblocks
    integer                                 :: Nlat,Nso
    complex(8),dimension(Nlat*Nso,Nlat*Nso) :: Matrix
    integer                                 :: i,j,ip
    Matrix=zero
    do ip=1,Nlat
       i = 1 + (ip-1)*Nso
       j = ip*Nso
       Matrix(i:j,i:j) =  Vblocks(ip,:,:)
    enddo
  end function blocks_to_matrix




  !-----------------------------------------------------------------------------
  ! Purpose:
  !-----------------------------------------------------------------------------
  function trace_matrix(M,dim) result(tr)
    integer                       :: dim
    complex(8),dimension(dim,dim) :: M
    complex(8) :: tr
    integer                       :: i
    tr=dcmplx(0d0,0d0)
    do i=1,dim
       tr=tr+M(i,i)
    enddo
  end function trace_matrix



  !-----------------------------------------------------------------------------
  ! Purpose:
  !-----------------------------------------------------------------------------
  function nn2so_reshape(Hnn,Nspin,Norb) result(Hso)
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: Hnn
    integer                                     :: Nspin,Norb
    complex(8),dimension(Nspin*Norb,Nspin*Norb) :: Hso
    integer                                     :: iorb,ispin,is
    integer                                     :: jorb,jspin,js
    Hso=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Hso(is,js) = Hnn(ispin,jspin,iorb,jorb)
             enddo
          enddo
       enddo
    enddo
  end function nn2so_reshape




  !-----------------------------------------------------------------------------
  ! Purpose:
  !-----------------------------------------------------------------------------
  function eye(n) result(A)
    integer, intent(in) :: n
    real(8)             :: A(n, n)
    integer             :: i
    A = 0d0
    do i = 1, n
       A(i,i) = 1d0
    end do
  end function eye



  !+------------------------------------------------------------------+
  !PURPOSE  : calculate the Heaviside  function
  !+------------------------------------------------------------------+
  elemental function heaviside(x)
    real(8),intent(in) :: x
    real(8)            :: heaviside
    if(x < 0.d0) then
       heaviside = 0.0d0
    elseif(x==0.d0)then
       heaviside = 0.50d0
    else
       heaviside = 1.0d0
    endif
  end function heaviside



  !+-------------------------------------------------------------------+
  !PURPOSE  : calculate the Fermi-Dirac distribution
  !+-------------------------------------------------------------------+
  elemental function fermi(x,beta)
    real(8),intent(in) :: x, beta 
    real(8)            :: fermi
    if(x*beta > 100.d0)then
       fermi=0.d0
       return
    endif
    fermi = 1.d0/(1.d0+exp(beta*x))
  end function fermi




  !-----------------------------------------------------------------------------
  ! Purpose:
  !-----------------------------------------------------------------------------
  subroutine Inv(Am)
    ! Inverts the general complex matrix Am
    complex(8), intent(inout) :: Am(:,:)   ! Matrix to be inverted
    complex(8), allocatable   :: Amt(:,:), work(:)
    integer                   :: n, nb
    integer                   :: lwork, info
    integer, allocatable      :: ipiv(:)
    real(8),external ::  ilaenv
    n = size(Am, 1)
    if(any(shape(Am)/=[n,n]))stop "error Inv: wrong shape(Am)"
    nb = ilaenv(1, 'ZGETRI', "UN", n, -1, -1, -1)  ! TODO: check UN param
    if (nb < 1) nb = max(1, n)
    lwork = n*nb
    allocate(Amt(n,n), ipiv(n), work(lwork))
    Amt = Am
    call zgetrf(n, n, Amt, n, ipiv, info)
    if (info /= 0) then
       print *, "zgetrf returned info =", info
       if (info < 0) then
          print *, "the", -info, "-th argument had an illegal value"
       else
          print *, "U(", info, ",", info, ") is exactly zero; The factorization"
          print *, "has been completed, but the factor U is exactly"
          print *, "singular, and division by zero will occur if it is used"
          print *, "to solve a system of equations."
       end if
       stop 'Z_mat_invert: zgetrf error'
    end if
    call zgetri(n, Amt, n, ipiv, work, lwork, info)
    if (info /= 0) then
       print *, "zgetri returned info =", info
       if (info < 0) then
          print *, "the", -info, "-th argument had an illegal value"
       else
          print *, "U(", info, ",", info, ") is exactly zero; the matrix is"
          print *, "singular and its inverse could not be computed."
       end if
       stop 'Z_mat_invert: zgetri error'
    end if
    Am = Amt
  end subroutine Inv




END MODULE COMMON
