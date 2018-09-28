program trans
  USE COMMON
  USE INPUT_VARS
  USE BUILD_H
  USE TRANSPORT
  implicit none
  !
  !
  integer                                       :: Nk,Nlso
  integer                                       :: ilat,ineq,ispin,iorb
  !
  !
  ! Green's functions
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Greal           ![Nlat][Nspin]{Nspin][Norb][Norb][Lreal]
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Sreal           ![Nlat][Nspin]{Nspin][Norb][Norb][Lreal]
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Sreal_ineq      ![Nineq][Nspin]{Nspin][Norb][Norb][Lreal]
  !
  !
  ! read input
  call read_input('inputTRANS.conf')
  !
  ! set input structure hamiltonian: allocates and sets Hij & nanoHloc
  !
  ! instead of building Hij, implement also a readin option
  ! format: wannier90-like (also compatible with build_Hij)
  ! in that case the input file would also be much "slimmer"
  ! note that many integers (e.g., Nlat, Nlso, ...) are set within build_Hij!
  if(read_Hij)then
     write(*,*)"error: readin option for Hij not implemented yet!"
  else
     call build_Hij([nfile,hijfile])
  endif
  !
  ! allocate and read self-energy
  allocate(Sreal(Nlat,Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Sreal_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lreal))
  call read_sigma(Sreal_ineq)
  do ilat=1,Nlat
     ineq = lat2ineq(ilat)
     Sreal(ilat,:,:,:,:,:) = Sreal_ineq(ineq,:,:,:,:,:)
  enddo
  !
  ! evaluate the linear response (zero-bias) transmission function 
  ! if jbias=T evaluate the corresponding bias-driven current
  call eval_transport(Hij,Sreal)
  !
end program trans
