MODULE GREENS_FUNCTIONS
  USE COMMON
  USE INPUT_VARS
  implicit none
  private

  public :: get_gij_wr
  public :: get_self_lead

contains



  !----------------------------------------------------------------------------------------!
  ! purpose: evaluate the non-local Green's function for a given frequency wr
  !----------------------------------------------------------------------------------------!
  subroutine get_gij_wr(Hij,Greal,Sreal,wr,self_lead)
    complex(8),dimension(:,:),intent(in)            :: Hij       ![Nlso][Nlso]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Greal     ![Nlat][Nlat][Nspin][Nspin][Norb][Norb]
    complex(8),dimension(:,:,:,:,:),intent(in)      :: Sreal     !      [Nlat][Nspin][Nspin][Norb][Norb]
    complex(8),dimension(:,:),intent(in),optional   :: self_lead ![Nlso][Nlso]
    complex(8),dimension(:,:,:),allocatable         :: zeta_real ![Nlat][Nspin*Norb][Nspin*Norb]
    real(8),intent(in)                              :: wr
    integer                                         :: Nso,Nlso,ilat
    !
    !
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    !Testing part:
    if(any(shape(Hij)/=[Nlso,Nlso]))stop "get_gij_wr error: wrong shape(Hij)"
    if(any(shape(Sreal)/=[Nlat,Nspin,Nspin,Norb,Norb]) )stop "get_gij_wr error: wrong shape(Sreal) "
    if(any(shape(Greal)/=[Nlat,Nlat,Nspin,Nspin,Norb,Norb]) )stop "get_gij_wr error: wrong shape(Greal) "  
    if(present(self_lead).AND.any(shape(self_lead)/=[Nlso,Nlso]) )stop "get_gij_wr error: wrong shape(self_lead) "
    !
    allocate(zeta_real(Nlat,Nso,Nso));zeta_real=zero
    !
    do ilat=1,Nlat
       zeta_real(ilat,:,:) = (wr+xi*eps+xmu)*eye(Nso) - nn2so_reshape(Sreal(ilat,:,:,:,:),NSpin,Norb)
    enddo
    !
    !pass each Z_site to the routines that invert (Z-Hk) for each k-point 
    Greal=zero
    if(present(self_lead))then
       call invert_gij(zeta_real,Hij,Greal,self_lead)
    else
       call invert_gij(zeta_real,Hij,Greal)
    endif
  end subroutine get_gij_wr



  !----------------------------------------------------------------------------------------!
  ! purpose: embed self_lead into Gij
  !----------------------------------------------------------------------------------------!
  subroutine invert_gij(zeta,Hk,Gij,self_lead)
    complex(8),dimension(:,:,:),intent(in)          :: zeta      ![Nlat][Nspin*Norb][Nspin*Norb]
    complex(8),dimension(:,:),intent(in)            :: Hk        ![Nlso][Nlso]
    complex(8),dimension(:,:,:,:,:,:),intent(inout) :: Gij       ![Nlat][Nlat][Nspin][Nspin][Norb][Norb]
    complex(8),dimension(:,:),intent(in),optional   :: self_lead ![Nlso][Nlso]
    !allocatable arrays
    complex(8),dimension(:,:),allocatable           :: Gmatrix   ![Nlso][Nlso]
    integer                                         :: Nso,Nlso
    integer                                         :: ilat,jlat,iorb,jorb,ispin,jspin,io,jo
    !
    Nso   = Nspin*Norb
    Nlso  = Nlat*Nspin*Norb
    !
    allocate(Gmatrix(Nlso,Nlso))
    Gij=zero
    Gmatrix  = blocks_to_matrix(zeta(:,:,:),Nlat,Nso) - Hk
    if(present(self_lead)) Gmatrix = Gmatrix - self_lead
    call inv(Gmatrix) 
    !store the diagonal blocks directly into the tmp output 
    do ilat=1,Nlat
       do jlat=1,Nlat
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      io = iorb + (ispin-1)*Norb + (ilat-1)*Nspin*Norb
                      jo = jorb + (jspin-1)*Norb + (jlat-1)*Nspin*Norb
                      Gij(ilat,jlat,ispin,jspin,iorb,jorb) = Gmatrix(io,jo)
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
  end subroutine invert_gij





  !----------------------------------------------------------------------------------------!
  ! purpose: define the hybridization matrix of size [Nlat][Nlat][Nspin][Norb][Norb][Lreal] 
  ! reading the parameters from an input file
  !----------------------------------------------------------------------------------------!
  subroutine get_self_lead(wr,self_lead)
    real(8),intent(in)                      :: wr
    complex(8),dimension(:,:),intent(inout) :: self_lead(:,:) ![Nlso][Nlso]
    complex(8),dimension(:,:),allocatable   :: lead_real ![Nlead][Nspin]
    integer                                 :: ilat,jlat,ispin,jspin,iorb,jorb,io,jo,i,Nlso
    integer                                 :: unit,l,lfile
    integer                                 :: ikind,ilead,Nlead
    real(8)                                 :: D,mu,V,epsk
    integer                                 :: k,kmax
    complex(8)                              :: ksum
    character(50)                           :: suffix
    !
    Nlso = size(self_lead,1)
    !
    ! initialize embedding hybridization function
    self_lead=zero

    ! determine Nleads & allocate lead matrix
    lfile = get_file_length("lead.in")
    open(newunit=unit,file='lead.in',status='old')
    read(unit,*)Nlead
    !
    allocate(lead_real(Nlead,Nspin))
    lead_real(:,:)=zero

    ! lead file setup lead by kind, half-bandwitdh (D) and chemical potential (mu)
    ! *** note: reading the file for each wr is slow and inefficient 
    do l=1,lfile-1 ! because Nlead was read separately above
       read(unit,*) ilead, ispin, D, mu, ikind
       ilead=ilead+1
       ispin=ispin+1
       if(ilead>Nlead)stop "set_hyb error: in input file 'lead.in' ilead > Nlead"
       if(ispin>Nspin)stop "set_hyb error: in input file 'lead.in' ispin > Nspin"
       !
       ! set the lead's Green's function, depending on ikind
       if(ikind==0)then
          ! flat DOS (analytic)
          !write(*,*) "flat DOS (analytic)"
          lead_real(ilead,ispin)=dcmplx( log(abs((D+wr+mu)/(D-wr-mu))) , -pi*heaviside(D-abs(wr+mu)) )/(2d0*D)
       elseif(ikind==1)then
          ! flat DOS (k-sum)
          !write(*,*) "flat DOS (k-sum)"
          kmax=10000
          ksum=zero
          do k=1,kmax
             epsk = -D + 2*D/kmax*(k-1)
             ksum = ksum + 1d0/( wr+xi*eps+mu - epsk)
          enddo
          lead_real(ilead,ispin)=ksum/kmax
       elseif(ikind==2)then
          ! broad-band limit
          !write(*,*) "broad-band limit (analytic)" 
          lead_real(ilead,ispin)=dcmplx(0d0,-1.d0*pi) ! to ensure DOS normalization
       elseif(ikind==3)then
          ! semicircular DOS (k-sum) 
          !write(*,*) "semicircular DOS (k-sum)"
          ksum=zero
          do k=1,kmax
             epsk = -D + 2*D/kmax*(k-1)
             ksum = ksum + (4d0/(pi*kmax))*sqrt(1d0-(epsk/D)**2)/( wr+xi*eps+mu - epsk)
          enddo
          lead_real(ilead,ispin)=ksum
       elseif(ikind==4)then
          ! readin hk DOS
          write(*,*) "readin hk DOS to be implemented and benchmarked w/ w2dynamics"
          stop
       else
          write(*,*) "set_hyb error: in input file 'lead.in' invalid ikind"
          stop
       endif
       !*** broken if wr is not an array: fix it!
       !! store lead(s) DOS on disk
       !!suffix="_ilead"//reg(txtfy(ilead))//"_s"//reg(txtfy(ispin))//"_realw.ed"
       !!call splot("lead"//trim(suffix),wr,lead_real(ilead,ispin,:))
    enddo
    close(unit)

    ! hybridization file determine lead-site connections 
    ! *** note: reading the file for each wr is slow and inefficient 
    lfile = get_file_length("vij.in")
    open(newunit=unit,file='vij.in',status='old')
    do i=1,lfile
       read(unit,*) ilat, iorb, jlat, jorb, ilead, V
       ilat=ilat+1
       iorb=iorb+1
       jlat=jlat+1
       jorb=jorb+1
       ilead=ilead+1
       if((iorb>Norb).or.(jorb>Norb))stop "set_hyb error: in input file 'vij.in' i/jorb > Norb"
       if((ilat>Nlat).or.(jlat>Nlat))stop "set_hyb error: in input file 'vij.in' i/jlat > Nlat"
       if(ilead>Nlead)stop "set_hyb error: in input file 'vij.in' ilead > Nlead"
       do ispin=1,Nspin
          ! get stride and set matrix element: no symmetrization
          io = iorb + (ispin-1)*Norb + (ilat-1)*Nspin*Norb !== ilat
          jo = jorb + (ispin-1)*Norb + (jlat-1)*Nspin*Norb !== jlat
          self_lead(io,jo) = self_lead(io,jo) + lead_real(ilead,ispin)*V**2
          !*** broken if wr is not an array: fix it!
          !! store self-energy of the lead(s) on disk
          !!suffix="_i"//reg(txtfy(ilat))//"_j"//reg(txtfy(jlat))//"_s"//reg(txtfy(ispin))//"_realw.ed"
          !!call splot("Hyb"//trim(suffix),wr,Hyb_real(io,jo,:))
       enddo
    enddo
    close(unit)
    deallocate(lead_real)
  end subroutine get_self_lead




END MODULE GREENS_FUNCTIONS










