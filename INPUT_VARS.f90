MODULE INPUT_VARS
  implicit none


  !input variables
  !=========================================================
  integer                                       :: Nspin      !spin degeneracy (max 2)
  real(8)                                       :: xmu        !chemical potential
  real(8)                                       :: beta       !inverse temperature
  real(8)                                       :: eps        !broadening
  real(8)                                       :: wini,wfin
  integer                                       :: Lreal
  logical                                       :: jbias
  character(len=40)                             :: Sfile
  character(len=40)                             :: Nfile
  character(len=40)                             :: hijfile
  logical                                       :: read_Hij
  
  namelist/variables/Nspin,xmu,beta,eps,wini,wfin,Lreal,jbias,Sfile,Nfile,HijFile,Read_Hij

contains



  !+-------------------------------------------------------------------+
  !PURPOSE  : READ THE INPUT FILE AND SETUP GLOBAL VARIABLES
  !+-------------------------------------------------------------------+
  subroutine read_input(INPUTunit)
    character(len=*)   :: INPUTunit
    integer            :: i,n_command_arg,pos
    character(len=128) :: arg_buffer,nml_value,nml_name
    logical            :: bool
    Nspin=1
    xmu=0d0
    beta=100d0
    eps=0.03d0
    wini=-15d0
    wfin=15d0
    Lreal=1000
    jbias=.false.
    Sfile="LSigma"
    Nfile="nano.in"
    Hijfile="hij.in"
    read_hij=.false.
    !
    n_command_arg=command_argument_count()
    inquire(file=trim(INPUTunit),EXIST=bool)
    if(bool)then
       open(10,file=trim(INPUTunit))          
       read(10,nml=variables)
       close(10)
    else
       write(*,*)"Using default variables!"
    endif
    if(n_command_arg/=0)then
       do i=1,n_command_arg
          call get_command_argument(i,arg_buffer)
          pos       = scan(arg_buffer,"=")
          nml_name  = arg_buffer(1:pos-1)
          nml_value = arg_buffer(pos+1:)
          print*,"name ",trim(nml_name)
          print*,"value ",trim(nml_value)
          select case(nml_name)
          case ("Nspin")
             read(nml_value,*)Nspin
             write(*,*)"Update variable "//trim(nml_name)//" with: ",nml_value
          case("xmu")
             read(nml_value,*)xmu
             write(*,*)"Update variable "//trim(nml_name)//" with: ",nml_value
          case("beta")
             read(nml_value,*)beta
             write(*,*)"Update variable "//trim(nml_name)//" with: ",nml_value
          case ("eps")
             read(nml_value,*)eps
             write(*,*)"Update variable "//trim(nml_name)//" with: ",nml_value
          case ("wini")
             read(nml_value,*)wini
             write(*,*)"Update variable "//trim(nml_name)//" with: ",nml_value
          case ("wfin")
             read(nml_value,*)wfin
             write(*,*)"Update variable "//trim(nml_name)//" with: ",nml_value
          case ("Lreal")
             read(nml_value,*)Lreal
             write(*,*)"Update variable "//trim(nml_name)//" with: ",nml_value
          case ("jbias")
             read(nml_value,*)jbias
             write(*,*)"Update variable "//trim(nml_name)//" with: ",nml_value
          case ("Nfile")
             read(nml_value,*)Nfile
             write(*,*)"Update variable "//trim(nml_name)//" with: ",nml_value
          case ("Hijfile")
             read(nml_value,*)Hijfile
             write(*,*)"Update variable "//trim(nml_name)//" with: ",nml_value
          case ("Sfile")
             read(nml_value,*)Sfile
             write(*,*)"Update variable "//trim(nml_name)//" with: ",nml_value
          case ("read_Hij")
             read(nml_value,*)read_Hij
             write(*,*)"Update variable "//trim(nml_name)//" with: ",nml_value
          case default
             print*,"No vars update from CL"
          end select
       enddo
    endif
    !
    write(*,nml=variables)
    !
    open(10,file="used."//trim(INPUTunit))
    write(10,nml=variables)
    close(10)
    !
  end subroutine read_input


END MODULE INPUT_VARS
