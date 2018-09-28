MODULE TRANSPORT
  USE COMMON
  USE INPUT_VARS
  USE BUILD_H
  USE GREENS_FUNCTIONS
  implicit none
  private
  
  public :: eval_transport

contains


  !----------------------------------------------------------------------------------------!
  ! purpose: evaluate 
  !  - conductance (without vertex corrections) 
  !  - bias-driven current
  ! for a nanostructure on the real axis, given the non-local Green's function 
  ! and the L/R hybridization matrix, of size [Nlat*Nspin*Norb**2*Lreal]
  !----------------------------------------------------------------------------------------!
  subroutine eval_transport(Hij,Sreal)
    ! inputs: Hamiltonian and retarded self-energy
    complex(8),intent(in)                         :: Hij(:,:,:)           ![Nlso][Nlso][Nk=1]
    complex(8),intent(in)                         :: Sreal(:,:,:,:,:,:)   ![Nlat][Nspin][Nspin][Norb][Norb][Lreal]
    ! retarded Green's function
    complex(8),dimension(:,:,:,:,:,:),allocatable :: Greal(:,:,:,:,:,:)   ![Nlat][Nlat][Nspin][Nspin][Norb][Norb]
    ! hybridization function to environment
    complex(8),dimension(:,:),allocatable         :: self_lead            ![Nlso][Nlso]
    !
    complex(8),dimension(:,:),allocatable         :: GR,HR,GA,HL,Re,Le,Te ![Nlat*Norb][Nlat*Norb]
    integer,dimension(:,:),allocatable            :: rmask,lmask          ![Nlat*Norb][Nlat*Norb]
    !
    real(8),dimension(:),allocatable              :: wr
    !
    complex(8),dimension(:,:),allocatable         :: transe               ![Nspin][Lreal]
    real(8),dimension(:),allocatable              :: jcurr                ![Nspin]
    !
    real(8)                                       :: lbias,rbias
    !
    integer                                       :: ilat,jlat,ispin,jspin,iorb,jorb,io,jo,is,js,i,Nlso,Nlo
    integer                                       :: unit,unit_in,unit_out,eof,lfile
    character(len=30)                             :: suffix
    !
    Nlso=Nlat*Nspin*Norb
    Nlo=Nlat*Norb
    !
    allocate(wr(Lreal))
    wr = linspace(wini,wfin,Lreal)
    !
    ! allocate variables for matrix-matrix multiplication
    allocate(GR(Nlo,Nlo));GR=zero
    allocate(HR(Nlo,Nlo));HR=zero
    allocate(GA(Nlo,Nlo));GA=zero
    allocate(HL(Nlo,Nlo));HL=zero
    allocate(Re(Nlo,Nlo));Re=zero
    allocate(Le(Nlo,Nlo));Le=zero
    allocate(Te(Nlo,Nlo));Te=zero
    !
    ! -----------------------------------------------------------
    ! temporary patch:
    ! to be back-compatible with masks from Norb=1 calculations
    ! the following loop does not read the orbital indexes
    if(Norb==1)then
       ! set masks in latiice indexes
       allocate(lmask(Nlo,Nlo),rmask(Nlo,Nlo))
       lmask(:,:)=0
       rmask(:,:)=0
       lfile = get_file_length("lmask.in")
       open(newunit=unit,file='lmask.in',status='old')
       do i=1,lfile
          read(unit,*) ilat,jlat !does not read iorb & jorb
          ilat=ilat+1
          jlat=jlat+1
          lmask(ilat,jlat)=1
          write(6,*) ilat,jlat,lmask(ilat,jlat)
       enddo
       lfile = get_file_length("rmask.in")
       open(newunit=unit,file='rmask.in',status='old')
       do i=1,lfile
          read(unit,*) ilat,jlat
          ilat=ilat+1
          jlat=jlat+1
          rmask(ilat,jlat)=1
          write(6,*) ilat,jlat,rmask(ilat,jlat)
       enddo
    else
       ! set masks in latiice-orbital indexes
       allocate(lmask(Nlo,Nlo),rmask(Nlo,Nlo))
       lmask(:,:)=0
       rmask(:,:)=0
       lfile = get_file_length("lmask.in")
       open(newunit=unit,file='lmask.in',status='old')
       do i=1,lfile
          read(unit,*) ilat,iorb,jlat,jorb
          ilat=ilat+1
          iorb=iorb+1
          jlat=jlat+1
          jorb=jorb+1
          io = iorb +  (ilat-1)*Norb
          jo = jorb +  (jlat-1)*Norb
          lmask(io,jo)=1
          write(6,*) ilat,iorb,jlat,jorb,lmask(io,jo)
       enddo
       lfile = get_file_length("rmask.in")
       open(newunit=unit,file='rmask.in',status='old')
       do i=1,lfile
          read(unit,*) ilat,iorb,jlat,jorb
          ilat=ilat+1
          iorb=iorb+1
          jlat=jlat+1
          jorb=jorb+1
          io = iorb +  (ilat-1)*Norb
          jo = jorb +  (jlat-1)*Norb
          rmask(io,jo)=1
          write(6,*) ilat,iorb,jlat,jorb,rmask(io,jo)
       enddo
    endif
    !
    ! allocate spin-resolved transmission coefficient
    allocate(transe(Nspin,Lreal))
    !
    ! allocate self-energy of the leads for a give frequency 
    allocate(self_lead(Nlso,Nlso))
    !
    ! allocate non-local Green's function for a given frequency
    allocate(Greal(Nlat,Nlat,Nspin,Nspin,Norb,Norb))
    !
    do i=1,Lreal
       ! evaluate hybridization function for wr(i)
       call get_self_lead(wr(i),self_lead)
       ! evaluate non-local Green's function for wr(i)
       call get_gij_wr(Hij(:,:,1),Greal,Sreal(:,:,:,:,:,i),wr(i),self_lead)
       !
       do ispin=1,Nspin
          ! fill auxiliary matrix [Nlso][Nlso]
          do ilat=1,Nlat
             do jlat=1,Nlat
                do iorb=1,Norb
                   do jorb=1,Norb
                      io = iorb +  (ilat-1)*Norb
                      jo = jorb +  (jlat-1)*Norb
                      is = iorb + (ispin-1)*Norb + (ilat-1)*Nspin*Norb !== ilat
                      js = jorb + (ispin-1)*Norb + (jlat-1)*Nspin*Norb !== jlat
                      !
                      ! auxiliary retarded Green's function
                      GR(io,jo)=Greal(ilat,jlat,ispin,ispin,iorb,jorb)
                      !
                      ! set \Gamma matrix for L/R according to masks to select L-subset OR R-subset
                      ! R-subset
                      HR(io,jo)=zero
                      if(rmask(io,jo)==1) HR(io,jo) = cmplx(2.d0*dimag(self_lead(is,js)),0d0)
                      ! L-subset
                      HL(io,jo)=zero
                      if(lmask(io,jo)==1) HL(io,jo) = cmplx(2.d0*dimag(self_lead(is,js)),0d0)
                   enddo
                enddo
             enddo
          enddo
          ! advanced Green's function
          GA=conjg(transpose(GR))
          !
          ! get transmission function as T(ispin,i)=Tr[Gadvc*Hybl*Gret*Hybr]
          Re = matmul(GR,HR)
          Le = matmul(GA,HL)
          Te = matmul(Le,Re)
          transe(ispin,i) = trace_matrix(Te,Nlo)
       enddo
    enddo
    !
    ! write transport coefficient of disk
    do ispin=1,Nspin
       write(suffix,"(A2,I1,A10)")"_s",ispin,"_realw.dat"
       open(newunit=unit,file="Te"//trim(suffix))
       do i=1,Lreal
          write(unit,*)wr(i),dimag(transe(ispin,i)),dreal(transe(ispin,i))
       enddo
    enddo
    !
    deallocate(GR,HR,GA,HL)
    deallocate(rmask,lmask)
    deallocate(Re,Le)
    !
    if(jbias)then
       !
       ! evaluate spin-resolved current as:
       ! J = \int_{-\infty}^{\infty} de T(e)*(f_L(e)-f_R(e))
       allocate(jcurr(Nspin));jcurr=0.d0
       !
       open(newunit=unit_in,file='jbias.in',status='old')
       open(newunit=unit_out,file="jbias.ed")
       do
          read(unit_in,*,IOSTAT=EOF)lbias,rbias
          if(EOF<0)exit
          !
          ! write L/R bias voltages
          write(unit_out,'(2f16.9)',advance='no')lbias,rbias
          !
          jcurr=0.d0
          do ispin=1,Nspin
             do i=1,Lreal
                jcurr(ispin) = jcurr(ispin) + transe(ispin,i)* &
                     (fermi(wr(i)-lbias,beta)-fermi(wr(i)-rbias,beta))* &
                     abs(wfin-wini)/Lreal
             enddo
             !
             ! write spin-resolved current on disk
             write(unit_out,'(1f16.9)',advance='no')jcurr(ispin)
          enddo
          write(unit_out,*) ! newline
       enddo
       close(unit_in)
       close(unit_out)
       !
       deallocate(jcurr)
       !
    endif
    !
    deallocate(Te)
    !
    return
    !
  end subroutine eval_transport


END MODULE TRANSPORT
