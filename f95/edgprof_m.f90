!> @addtogroup groupname4
!> @{
module edgprof_m
!> @}
  use edgprof_h
  use log_m
  use misc_m
  use const_numphys_h
  use const_kind_m
  use pcontrol_h

  implicit none
  private

! public subroutines
  public :: &
  edgprof_init,  & !< open file
  edgprof_readcon,  & !< read data from file
  edgprof_factors,  & !< factors for data
  edgprof_exp,  & !< exponential profile
  edgprof_expdouble,  & !< double exponential profile
  edgprof_eich,  & !< Eich profile
  edgprof_userdefined,  & !< userdefined profile
  edgprof_samples,  & !< samples profile
  edgprof_fn, &  !< general external function call
  edgprof_region, &  !< general external function call
  edgprof_initwrite, & !< open new output file
  edgprof_write, &  !< write out object
  edgprof_delete, & !< delete object
  edgprof_close, & !< close file
  edgprof_closewrite, & !< close write file
  edgprof_lambda_q !< Calculates_lambda_q

! private variables
  character(*), parameter :: m_name='edgprof_m' !< module name
  integer  :: status   !< error status
  integer, save  :: ninep=5     !< control file unit number
  integer, save  :: noutep=6      !< output file unit number
  character(len=80), save :: controlfile !< control file name
  character(len=80), save :: outputfile !< output file name
  integer  :: ilog      !< for namelist dump after error
  integer(ki4) :: i !< loop counter
  integer(ki4) :: j !< loop counter
  integer(ki4) :: k !< loop counter
  integer(ki4) :: l !< loop counter
  integer(ki4) :: ij !< loop counter

  contains
!---------------------------------------------------------------------
!> open file
subroutine edgprof_init(file,channel)

  !! arguments
  character(*), intent(in) :: file !< file name
  integer, intent(out),optional :: channel   !< input channel for object data structure
  !! local
  character(*), parameter :: s_name='edgprof_init' !< subroutine name
  !! logical :: unitused !< flag to test unit is available

  !! get file unit do i=99,1,-1 inquire(i,opened=unitused) if(.not.unitused)then ninep=i channel=i exit end if end do

  !! open file
  controlfile=trim(file)
  call log_value("Control data file",trim(controlfile))
  call misc_getfileunit(ninep)
  open(unit=ninep,file=controlfile,status='OLD',iostat=status)
  if(status/=0)then
     !! error opening file
     print '("Fatal error: Unable to open control file, ",a)',controlfile
     call log_error(m_name,s_name,1,error_fatal,'Cannot open control data file')
     stop
  end if
  if (present(channel)) channel=ninep

end subroutine edgprof_init
!---------------------------------------------------------------------
!> read data from file
subroutine edgprof_readcon(self,pnumerics,kin)

  !! arguments
  type(edgprof_t), intent(out) :: self !< type which data will be assigned to
  type(pnumerics_t), intent(inout) :: pnumerics !< powcal general numerical paramete
  integer, intent(in),optional :: kin   !< input channel for object data structure

  !! local
  character(*), parameter :: s_name='edgprof_readcon' !< subroutine name
  character(len=80) :: profile_formula(4) !< formula to be used
  integer(ki4), parameter :: MAX_NUMBER_OF_POSITIONS=500 !< maximum number of positions allowed
  integer(ki4), parameter :: MAX_NUMBER_OF_PARAMETERS=10 !< maximum number of parameters allowed
  real(kr8) :: power_split !< power split ion to electron direction
  real(kr8), dimension(4) :: decay_length !< power decay length at midplane (m)
  real(kr8) :: power_loss !< power crossing the separatrix (W)
  real(kr8),dimension(4) :: diffusion_length !< (m)
  real(kr8) :: q_parallel0 !< parallel reference flux (W/sqm)
  real(kr8) :: ratio_of_q_parallel0_near !< near-SOL parallel relative reference flux (dimensionless)
  real(kr8) :: decay_length_near !< near-SOL power decay length at midplane (m)

  real(kr8), dimension(MAX_NUMBER_OF_POSITIONS) :: positions  !< local variable
  real(kr8), dimension(MAX_NUMBER_OF_POSITIONS) :: deposition_profile  !< local variable
  real(kr8):: c_prop_1  !< Constant of proportionality for lambda_q custom-1
  real(kr8):: c_prop_2  !< Constant of proportionality for lambda_q custom-2
  real(kr8):: c_prop_3  !< Constant of proportionality for lambda_q custom-3
  real(kr8):: c_prop_4  !< Constant of proportionality for lambda_q custom-4
  real(kr8):: power_loss_exp !< Exponent for power loss in lambda_q custom
  real(kr8):: k_0_exp !< Exponent for k_0 in lambda_q custom
  real(kr8):: q_exp !< Exponent for q in lambda_q custom
  real(kr8):: a_exp !< Exponent for a in lambda_q custom
  real(kr8):: R_exp !< Exponent for R in lambda_q custom
  real(kr8):: L_exp !< Exponent for L in lambda_q custom
  real(kr8):: ratio_B_exp !< Exponent for ratio_B in lambda_q custom
  real(kr8):: P_V_ratio_exp !< Exponent for P_V_ratio in lambda_q custom
  real(kr8):: a_R0_elongation_ratio_exp !< Exponent for a_R0_elogation_ratio in lambda_q custom
  integer(ki4) :: number_of_positions  !< local variable
  character(len=80) :: coord_of_positions !< whether in distance or flux from SOL
  character(len=80), dimension(4) :: lambda_q_choice !< Choice of calculation of lambda_q

  real(kr8), dimension(MAX_NUMBER_OF_PARAMETERS) :: general_real_parameters  !< local variable
  real(kr8) :: pi=3.1415926535 !< Pi
  real(kr8) :: n_u!< Plasma upstream density
  real(kr8) :: e_charge=1.602E-19 !< Electron charge
  real(kr8) :: k_0 !< Electron conductivity
  real(kr8) :: q !< Safety factor
  real(kr8) :: R !< 
  real(kr8) :: a !< 
  real(kr8) :: L !< Length of the scrape off layer
  real(kr8) :: ratio_B !< 
  real(kr8) :: P_V_ratio !<
  real(kr8) :: a_R0_ratio !<
  real(kr8) :: elongation !<
  real(kr8) :: x_perpendicular  !< 
  integer(ki4), dimension(MAX_NUMBER_OF_PARAMETERS) :: general_integer_parameters  !< local variable
  integer(ki4) :: number_of_real_parameters  !< local variable
  integer(ki4) :: number_of_integer_parameters  !< local variable
  
  !! edgprof parameters
  namelist /edgprofparameters/ &
 &profile_formula, &
 &power_split, decay_length, power_loss, &
 &q_parallel0, &
 &ratio_of_q_parallel0_near, &
 &diffusion_length, &
 &decay_length_near, &
 &positions, deposition_profile, coord_of_positions, number_of_positions,&
 &general_real_parameters, number_of_real_parameters, &
 &general_integer_parameters, number_of_integer_parameters, lambda_q_choice
 
   namelist /lambdaqparameters/ &
  &c_prop_1 , c_prop_2, c_prop_3, &
  &c_prop_4, power_loss_exp, n_u, &
  &k_0, k_0_exp, q, q_exp, &
  &R, R_exp, a, a_exp, &
  &x_perpendicular, L, L_exp,  &
  &ratio_B, ratio_B_exp, P_V_ratio_exp, &
  &P_V_ratio, a_R0_ratio, a_R0_elongation_ratio_exp, &
  &elongation

  !! Set default lambdaqparameters
  c_prop_1 = ((7.0d0/8.0d0)**(2.0d0/9.0d0))*((8.0d0*(pi**2.0d0)/7.0d0)**(7.0d0/9.0d0))  
  c_prop_2 = ((8.0d0*(pi**2.0d0)/7.0d0)**(7.0d0/9.0d0))  
  c_prop_3 = (8**(5.0d0/9.0d0))*(pi**(4.0d0/3.0d0))/(7.0d0**(5.0d0/9.0d0))
  c_prop_4 = 10.0d0
  power_loss_exp=-5.0d0/9.0d0
  n_u=1E+19
  k_0=2000.0d0
  k_0_exp=-2.0d0/9.0d0
  q=3.0d0
  q_exp=4.0/9.0d0
  R=2.5d0
  R_exp=5.0d0/9.0d0
  a=1.5d0
  a_exp=5.0d0/9.0d0
  x_perpendicular=2.0d0
  L=50.0d0
  L_exp=4.0d0/9.0d0
  ratio_B=0.3d0
  ratio_B_exp=-2.0d0/9.0d0
  P_V_ratio_exp=-0.38d0
  P_V_ratio = 4100.0d0
  a_R0_ratio=0.33d0
  a_R0_elongation_ratio_exp=1.3d0
  elongation=1.3d0


  !! set default edgprof parameters
  power_split=0.5_kr8
  decay_length=.012_kr8
  diffusion_length=0
  power_loss=7.5e+06_kr8
  !! ignore q_parallel0 if negative
  q_parallel0=-1_kr8
  ratio_of_q_parallel0_near=0

  lambda_q_choice='default'
  profile_formula='unset'
  decay_length_near=0
  positions=0
  deposition_profile=0
  coord_of_positions='radius'
  number_of_positions=0
  general_real_parameters=0
  general_integer_parameters=0
  number_of_real_parameters=0
  number_of_integer_parameters=0

  if(present(kin)) then
     !! assume unit already open and reading infile
     ninep=kin
  end if

  if (pnumerics%ledgprof) then
     !!read edgprof parameters
     read(ninep,nml=edgprofparameters,iostat=status)
  Do i=1,4
  formula_chosen_lambda_pre_read: select case (lambda_q_choice(i))
  case('custom-1') 
    read(ninep,nml=lambdaqparameters,iostat=status)  
  case('custom-2') 
    read(ninep,nml=lambdaqparameters,iostat=status)  
  case('custom-3') 
    a_exp=7.0d0/9.0d0
    L_exp=2.0d0/9.0d0
    read(ninep,nml=lambdaqparameters,iostat=status)
  case('custom-4') 
    read(ninep,nml=lambdaqparameters,iostat=status)
  end select formula_chosen_lambda_pre_read
  end do
  
     if(status/=0) then
        print '("Fatal error reading edgprofparameters or lambda q parameters")'
        call log_getunit(ilog)
        write(ilog,nml=edgprofparameters)
        call log_error(m_name,s_name,1,error_fatal,'Error reading edgprofparameters or lambdaqparameters')
     end if
  else
     !! set edgprof parameters using values from powcalparameters
     power_split=pnumerics%f
     decay_length=pnumerics%lmid
     diffusion_length=pnumerics%sigma
     power_loss=pnumerics%ploss
     q_parallel0=pnumerics%qpara0
     DO i=1,4
		if (diffusion_length(i)>1.e-6_kr8) then
			profile_formula='eich'
		else
			profile_formula='unset'
		end if
     END DO
  end if

  call lowor(profile_formula,1,len_trim(profile_formula))

  !!Calculates user choice of lambda_q
  DO i=1,4
  formula_chosen_lambda: select case (lambda_q_choice(i))
  case('default')
  
  case('custom-1') 
    decay_length(i)=c_prop_1*(power_loss**power_loss_exp)*&
                 (k_0**k_0_exp)*(q**q_exp)*R*(a**a_exp)*&
                 ((e_charge*n_u*x_perpendicular)**(7.0d0/9.0d0))
                 
  case('custom-2') 
    decay_length(i)=c_prop_2*(power_loss**power_loss_exp)*&
                 (k_0**k_0_exp)*(R**R_exp)*(a**a_exp)*&
                 ((e_charge*n_u*x_perpendicular)**(7.0d0/9.0d0))*(L**L_exp)
                 
  case('custom-3') 
    decay_length(i)=c_prop_3*(power_loss**power_loss_exp)*&
                 (k_0**k_0_exp)*(q**q_exp)*(R**R_exp)*(a**a_exp)*&
                 ((e_charge*n_u*x_perpendicular)**(7.0d0/9.0d0))*(L**L_exp)*&
                 (ratio_B**ratio_B_exp)
  case('custom-4') 
    decay_length(i)=c_prop_4*((P_V_ratio)**P_V_ratio_exp)*((a_R0_ratio/elongation)**a_R0_elongation_ratio_exp)
  end select formula_chosen_lambda
  end do
  !! check for valid data
DO j=1,4
  formula_chosen: select case (profile_formula(j))
  case('unset','exp')
     if(power_split<0.OR.power_split>1) &
 &   call log_error(m_name,s_name,11,error_fatal,'power_split must be >=0 and <=1')
 do i=1,4
     if(decay_length(i)<=0) &
 &   call log_error(m_name,s_name,12,error_fatal,'decay_length (m) must be >0')
 end do
     if(q_parallel0<=0.AND.power_loss<=0) &
 &   call log_error(m_name,s_name,14,error_fatal,'power_loss (W) must be >0')

  case('expdouble')
     if(power_split<0.OR.power_split>1) &
 &   call log_error(m_name,s_name,21,error_fatal,'power_split must be >=0 and <=1')
 do i=1,4
     if(decay_length(i)<=0) &
 &   call log_error(m_name,s_name,22,error_fatal,'decay_length (m) must be >0')
 end do
     if(q_parallel0<=0.AND.power_loss<=0) &
 &   call log_error(m_name,s_name,24,error_fatal,'power_loss (W) must be >0')
     if(ratio_of_q_parallel0_near<=0.AND.power_loss<=0) &
 &   call log_error(m_name,s_name,25,error_fatal,'ratio_of_q_parallel0_near (W) must be >=0')
     if(decay_length_near<0) &
 &   call log_error(m_name,s_name,22,error_fatal,'decay_length_near (m) must be >=0')

  case('eich')
     if(power_split<0.OR.power_split>1) &
 &   call log_error(m_name,s_name,31,error_fatal,'power_split must be >=0 and <=1')
 do i=1,4
     if(decay_length(i)<=0) &
 &   call log_error(m_name,s_name,32,error_fatal,'decay_length (m) must be >0')
     if(diffusion_length(i)<0) &
 &   call log_error(m_name,s_name,33,error_fatal,'diffusion_length (m) must be >=0')
  end do
     if(q_parallel0<=0.AND.power_loss<=0) &
 &   call log_error(m_name,s_name,34,error_fatal,'power_loss (W) must be >0')

  case('samples')
     if(coord_of_positions/='radius'.AND.coord_of_positions/='flux') &
 &   call log_error(m_name,s_name,39,error_fatal,'position type not recognised')
     if(q_parallel0<=0.AND.power_loss<=0) &
 &   call log_error(m_name,s_name,40,error_fatal,'power_loss (W) must be >0')
     if(number_of_positions<=0) &
 &   call log_error(m_name,s_name,41,error_fatal,'number of positions must be >0')
     if(number_of_positions>MAX_NUMBER_OF_POSITIONS) then
        call log_value("max number of positions",MAX_NUMBER_OF_POSITIONS)
        call log_error(m_name,s_name,42,error_fatal,'too many positions: increase MAX_NUMBER_OF_POSITIONS')
     end if

  case('userdefined')
     if(q_parallel0<=0.AND.power_loss<=0) &
 &   call log_error(m_name,s_name,43,error_fatal,'power_loss (W) must be >0')
     if(number_of_real_parameters<0) &
 &   call log_error(m_name,s_name,44,error_fatal,'number of real parameters must be >=0')
     if(number_of_real_parameters>MAX_NUMBER_OF_PARAMETERS) then
        call log_value("max number of real parameters",MAX_NUMBER_OF_PARAMETERS)
        call log_error(m_name,s_name,45,error_fatal,'too many parameters: increase MAX_NUMBER_OF_PARAMETERS')
     end if
     if(number_of_integer_parameters<0) &
 &   call log_error(m_name,s_name,46,error_fatal,'number of integer parameters must be >=0')
     if(number_of_integer_parameters>MAX_NUMBER_OF_PARAMETERS) then
        call log_value("max number of integer parameters",MAX_NUMBER_OF_PARAMETERS)
        call log_error(m_name,s_name,47,error_fatal,'too many parameters: increase MAX_NUMBER_OF_PARAMETERS')
     end if
     if(number_of_integer_parameters==0.AND.number_of_real_parameters==0) &
 &   call log_error(m_name,s_name,48,error_fatal,'no parameters set')

  end select formula_chosen
end do
  !! store values
  DO i=1,4
  self%formula(i)=TRIM(profile_formula(i))
  END DO
  self%f=power_split
  do i=1,4
  self%lmid(i)=decay_length(i)
  end do
  self%ploss=power_loss
  self%sigma=diffusion_length
  self%qpara0=q_parallel0
  self%lmidnr=decay_length_near
  self%rqpara0=ratio_of_q_parallel0_near
  self%postype=coord_of_positions

  !! allocate arrays and assign

  self%npos=number_of_positions
  self%nrpams=number_of_real_parameters
  self%nipams=number_of_integer_parameters
DO i=1,4
  formula_allocate: select case (profile_formula(i))
  case('samples')
     allocate(self%pos(number_of_positions), self%prof(number_of_positions), stat=status)
     call log_alloc_check(m_name,s_name,60,status)
     self%pos=positions(:number_of_positions)
     self%prof=deposition_profile(:number_of_positions)
     ! check positions
     do j=1,number_of_positions-1
        if (self%pos(j)<self%pos(j+1)) then
           cycle
        else
           call log_error(m_name,s_name,61,error_fatal,'positions are not monotone increasing')
        end if
     end do

  case('userdefined')
     if (number_of_real_parameters>0) allocate(self%rpar(number_of_real_parameters), stat=status)
     call log_alloc_check(m_name,s_name,65,status)
     if (number_of_integer_parameters>0) allocate(self%npar(number_of_integer_parameters), stat=status)
     call log_alloc_check(m_name,s_name,66,status)
     self%rpar=general_real_parameters(:number_of_real_parameters)
     self%npar=general_integer_parameters(:number_of_integer_parameters)
  case default
  end select formula_allocate
end do
end  subroutine edgprof_readcon

!---------------------------------------------------------------------
!> factors for power deposition
subroutine edgprof_factors(self,rbdry,bpbdry,btotbdry,psign)
  !! arguments
  type(edgprof_t), intent(inout) :: self !< type containing profile parameters
  real(kr8), intent(in) :: rbdry(2) !< boundary value
  real(kr8), intent(in) :: bpbdry !< boundary value
  real(kr8), intent(in) :: btotbdry!< boundary value
  real(kr8), intent(in) :: psign !< adjust sign of exponential factor

  !! local variables
  character(*), parameter :: s_name='edgprof_factors' !< subroutine name
  real(kr8) :: zrbfac(4) !< power factor in \f$ B \f$ only
  real(kr8),dimension(4) :: zrblfac !< power factor scaled by lmid
  integer(ki4) :: i   !<  local variable

  self%fpfacnr=0.0
  self%rblfacnr=0.0
  
  ! power normalisation factor
  zrbfac=1/(2*const_pid*rbdry(2)*bpbdry)
  ! default diffusion factor
  self%slfac=0
  DO i=1,4
  formula_chosen: select case (self%formula(i))
  case('unset','exp')
     zrblfac(i)=zrbfac(i)/self%lmid(i)
     self%rblfac(i)=2*const_pid*zrblfac(i)*((-1.)*psign)
     if (self%qpara0>0) then
        self%fpfac(i)=self%f*self%qpara0/btotbdry
     else
        self%fpfac(i)=self%f*self%ploss*zrblfac(i)
     end if
  case('expdouble')
     zrblfac(i)=zrbfac(i)/(self%lmid(i)+self%rqpara0*self%lmidnr(i))
     self%rblfac(i)=2*const_pid*zrbfac(i)*((-1.)*psign)/self%lmid(i)
     self%rblfacnr(i)=2*const_pid*zrbfac(i)*((-1.)*psign)/self%lmidnr(i)
     self%fpfac(i)=self%f*self%ploss*zrblfac(i)
     self%fpfacnr(i)=self%rqpara0*self%fpfac(i)

  case('eich')
     zrblfac=zrbfac/self%lmid
     self%slfac=self%sigma/(2*self%lmid)
     self%rblfac=2*const_pid*zrblfac*((-1.)*psign)
     self%fpfac=(self%f/2)*self%ploss*zrblfac

  case('samples')
     position_type: select case (self%postype)
     case('radius')
        self%rblfac=2*const_pid*zrbfac*((-1.)*psign)
     case default
        self%rblfac=1
     end select position_type
     self%fpfac=self%f*self%ploss*zrbfac

  case('userdefined')
     usrposition_type: select case (self%postype)
     case('radius')
        self%rblfac=2*const_pid*zrbfac*((-1.)*psign)
     case default
        self%rblfac=1
     end select usrposition_type
     self%fpfac=self%f*self%ploss*zrbfac/self%fint

  end select formula_chosen
end do
  !dbg write(*,*) "psign",psign,zrblfac,self%slfac,self%rblfac,self%fpfac

end subroutine edgprof_factors
!---------------------------------------------------------------------
!> exponential profile
function edgprof_exp(self,psid,j)

  !! arguments
  type(edgprof_t), intent(in) :: self !< type containing profile parameters
  real(kr8) :: edgprof_exp !< local variable
  real(kr8), intent(in) :: psid !< position in \f$ \psi \f$
  integer (ki4) :: j

  !! local variables
  character(*), parameter :: s_name='edgprof_exp' !< subroutine name
  real(kr8) :: pow !< local variable

  pow=self%fpfac(j)*exp(self%rblfac(j)*psid)

  !! return profile
  edgprof_exp=pow

end function edgprof_exp

!---------------------------------------------------------------------
!> double exponential profile
function edgprof_expdouble(self,psid,j)

  !! arguments
  type(edgprof_t), intent(in) :: self !< type containing profile parameters
  real(kr8) :: edgprof_expdouble !< local variable
  real(kr8), intent(in) :: psid !< position in \f$ \psi \f$
  integer (ki4) :: j
  
  !! local variables
  character(*), parameter :: s_name='edgprof_expdouble' !< subroutine name
  real(kr8) :: pow !< local variable
  
				  pow=self%fpfac(j)*exp(self%rblfac(j)*psid) &
 &+self%fpfacnr(j)*exp(self%rblfacnr(j)*psid)


  !! return profile
  edgprof_expdouble=pow

end function edgprof_expdouble

!---------------------------------------------------------------------
!> Eich profile
function edgprof_eich(self,psid,j)

  !! arguments
  type(edgprof_t), intent(in) :: self !< type containing profile parameters
  real(kr8) :: edgprof_eich !< local variable
  real(kr8), intent(in) :: psid !< position in \f$ \psi \f$
  integer (ki4) :: j
    
  !! local variables
  character(*), parameter :: s_name='edgprof_eich' !< subroutine name
  real(kr8) :: pow !< local variable

  !! calculate profile

  pow=self%fpfac(j)*&
  exp(self%slfac(j)**2+self%rblfac(j)*psid)*&
  erfc(self%slfac(j)+self%rblfac(j)*psid/(2*self%slfac(j)))

  !! return profile
  edgprof_eich=pow

end function edgprof_eich

!---------------------------------------------------------------------
!> userdefined profile
function edgprof_userdefined(self,psid,j)

  !! arguments
  type(edgprof_t), intent(in) :: self !< type containing profile parameters
  real(kr8) :: edgprof_userdefined !< local variable
  real(kr8), intent(in) :: psid !< position in \f$ \psi \f$
  integer (ki4) :: j

  !! local variables
  character(*), parameter :: s_name='edgprof_userdefined' !< subroutine name
  real(kr8) :: pow !< local variable
  real(kr8) :: zpos !< position
  real(kr8) :: zfint !< integral

  !! user defines profile here

  ! convert to r if necessary
  zpos=self%rblfac(j)*psid
  ! must define f(zpos)
  ! integral of f over r to normalise
  zfint=1
  pow=self%fpfac(j)/zfint

  !! return profile
  edgprof_userdefined=pow

end function edgprof_userdefined

!---------------------------------------------------------------------
!> profile from samples
function edgprof_samples(self,psid,j)

  !! arguments
  type(edgprof_t), intent(in) :: self !< type containing profile parameters
  real(kr8) :: edgprof_samples !< local variable
  real(kr8), intent(in) :: psid !< position in \f$ \psi \f$
  integer (ki4) :: j
    
  !! local variables
  character(*), parameter :: s_name='edgprof_samples' !< subroutine name
  real(kr8) :: pow !< local variable
  real(kr8) :: zpos !< position
  real(kr8) :: zpprof !< posn \f$ x \f$ profile
  real(kr8) :: zrh !< reciprocal mesh spacing
  integer(ki4) :: ip   !<  first integer coordinate of point
  integer(ki4) :: iflag   !<  warning flag

  !! calculate profile
  ! convert to r if necessary
  zpos=self%rblfac(j)*psid
  call interv(self%pos,self%npos,zpos,ip,iflag)
  if (iflag/=0) then
     call log_error(m_name,s_name,1,error_warning,'Point not in range of positions')
     iflag=0
  end if
  zrh=1/(self%pos(ip+1)-self%pos(ip))
  zpprof=(self%pos(ip+1)-zpos)*self%prof(ip) + (zpos-self%pos(ip))*self%prof(ip+1)
  pow=self%fpfac(1)*zrh*zpprof
  !! return profile
  edgprof_samples=pow

end function edgprof_samples


!---------------------------------------------------------------------
!> Determines the region needed for edgprof
function edgprof_region(R,Z,cenz,rxpt,psi,psixpt)

  !! arguments
  integer(ki4) :: edgprof_region !< local variable
  real(kr8) :: R !< position \f$ R \f$
  real(kr8) :: Z !< position  \f$ Z \f$
  real(kr8) :: cenz !< position centre of plasma \f$ Z \f$
  real(kr8) :: rxpt(2)   !<  position x point \f$ R \f$
  real(kr8), intent(in) :: psi !<  \f$ \psi \f$
  real(kr8) :: psixpt(2) !<  flux xpoints
  
  !! local variables
  character(*), parameter :: s_name='edgprof_region' !< subroutine name
  integer(ki4) :: psixplow !< lowest x-point psi
  integer(ki4) :: pow !< local variable
  
  psixplow=MINVAL(psixpt)
  
  		IF(Z<=cenz) THEN
			IF(R<=rxpt(1)) THEN
				IF(psi<=psixplow) THEN
					pow=1
				ELSE 
					pow=2
				END IF
			ELSE
				IF(psi<=psixplow) THEN
					pow=4
				ELSE 
					pow=3
				END IF			
			END IF
		ELSE
			IF(R<=rxpt(2)) THEN
				IF(psi<=psixplow) THEN
					pow=1
				ELSE 
					pow=2
				END IF
			ELSE
				IF(psi<=psixplow) THEN
					pow=4
				ELSE 
					pow=3
				END IF			
			END IF
		END IF

  !!return region
  edgprof_region=pow
 end function edgprof_region
!---------------------------------------------------------------------
!> profile from samples
function edgprof_fn(self,psi,psid,R,Z,cenz,rxpt,psixpt)

  !! arguments
  type(edgprof_t), intent(in) :: self !< type containing profile parameters
  real(kr8) :: edgprof_fn !< local variable
  real(kr8), intent(in) :: psid !< position in \f$ \psi \f$
  real(kr8), intent(in) :: psi !<  \f$ \psi \f$
  real(kr8) :: R !< position \f$ R \f$
  real(kr8) :: Z !< position  \f$ Z \f$
  real(kr8) :: cenz !< position centre of plasma \f$ Z \f$
  real(kr8) :: rxpt(2)   !<  position x point \f$ R \f$
  real(kr8) :: psixpt(2) !<  flux xpoints
  integer(ki4) :: i   !<  local variable
  
  !! local variables
  character(*), parameter :: s_name='edgprof_fn' !< subroutine name
  real(kr8) :: pow !< local variable

  i=edgprof_region(R,Z,cenz,rxpt,psi,psixpt)
  !! select profile
  formula_chosen: select case (self%formula(i))
  case('unset','exp')
	pow=edgprof_exp(self,psid,i)
  case('expdouble')
	pow=edgprof_expdouble(self,psid,i)
  case('eich')
    pow=edgprof_eich(self,psid,i)
  case('userdefined')
	pow=edgprof_userdefined(self,psid,i)
  case('samples')
    pow=edgprof_samples(self,psid,i)
  end select formula_chosen
	
  !! return profile
  edgprof_fn=pow

end function edgprof_fn
!---------------------------------------------------------------------
!> open new file
subroutine edgprof_initwrite(fileroot,channel)

  !! arguments
  character(*), intent(in) :: fileroot !< file root
  integer, intent(out),optional :: channel   !< output channel for object data structure
  !! local
  character(*), parameter :: s_name='edgprof_initwrite' !< subroutine name
  !! logical :: unitused !< flag to test unit is available
  character(len=80) :: outputfile !< output file name

  !! get file unit do i=99,1,-1 inquire(i,opened=unitused) if(.not.unitused)then channel=i exit end if end do noutep=i

  !! open file
  outputfile=trim(fileroot)//"_edgprof.out"
  call log_value("Control data file",trim(outputfile))
  call misc_getfileunit(noutep)
  open(unit=noutep,file=outputfile,status='NEW',iostat=status)
  if(status/=0)then
     open(unit=noutep,file=outputfile,status='REPLACE',iostat=status)
  end if
  if(status/=0)then
     !! error opening file
     print '("Fatal error: Unable to open output file, ",a)',outputfile
     call log_error(m_name,s_name,1,error_fatal,'Cannot open output data file')
     stop
  end if
  if (present(channel)) channel=noutep

end subroutine edgprof_initwrite

!---------------------------------------------------------------------
!> write edgprof data
subroutine edgprof_write(self,kout)

  !! arguments
  type(edgprof_t), intent(in) :: self   !< edgprof data structure
  integer, intent(in), optional :: kout   !< output channel for edgprof data structure

  !! local
  character(*), parameter :: s_name='edgprof_write' !< subroutine name
  integer :: iout   !< output channel for edgprof data structure

  !! sort out unit
  if(present(kout)) then
     iout=kout
  else
     iout=noutep
  end if

  write(iout,*,iostat=status) 'rblfac'
  call log_write_check(m_name,s_name,1,status)
  write(iout,*,iostat=status) self%rblfac
  call log_write_check(m_name,s_name,2,status)
  write(iout,*,iostat=status) 'fpfac'
  call log_write_check(m_name,s_name,3,status)
  write(iout,*,iostat=status) self%fpfac
  call log_write_check(m_name,s_name,4,status)
  write(iout,*,iostat=status) 'rblfacnr'
  call log_write_check(m_name,s_name,5,status)
  write(iout,*,iostat=status) self%rblfacnr
  call log_write_check(m_name,s_name,6,status)
  write(iout,*,iostat=status) 'fpfacnr'
  call log_write_check(m_name,s_name,7,status)
  write(iout,*,iostat=status) self%fpfacnr
  call log_write_check(m_name,s_name,8,status)
  write(iout,*,iostat=status) 'slfac'
  call log_write_check(m_name,s_name,9,status)
  write(iout,*,iostat=status) self%slfac
  call log_write_check(m_name,s_name,10,status)
  !
  write(iout,*,iostat=status) 'profile_formula'
  call log_write_check(m_name,s_name,18,status)
  write(iout,*,iostat=status) self%formula
  call log_write_check(m_name,s_name,19,status)
  write(iout,*,iostat=status) 'f'
  call log_write_check(m_name,s_name,20,status)
  write(iout,*,iostat=status) self%f
  call log_write_check(m_name,s_name,21,status)
  write(iout,*,iostat=status) 'lmid'
  call log_write_check(m_name,s_name,22,status)
  write(iout,*,iostat=status) self%lmid
  call log_write_check(m_name,s_name,23,status)
  write(iout,*,iostat=status) 'ploss'
  call log_write_check(m_name,s_name,24,status)
  write(iout,*,iostat=status) self%ploss
  call log_write_check(m_name,s_name,25,status)
  write(iout,*,iostat=status) 'sigma'
  call log_write_check(m_name,s_name,26,status)
  write(iout,*,iostat=status) self%sigma
  call log_write_check(m_name,s_name,27,status)
  write(iout,*,iostat=status) 'qpara0'
  call log_write_check(m_name,s_name,28,status)
  write(iout,*,iostat=status) self%qpara0
  call log_write_check(m_name,s_name,29,status)
  write(iout,*,iostat=status) 'lmidnr'
  call log_write_check(m_name,s_name,30,status)
  write(iout,*,iostat=status) self%lmidnr
  call log_write_check(m_name,s_name,31,status)
  write(iout,*,iostat=status) 'rqpara0'
  call log_write_check(m_name,s_name,32,status)
  write(iout,*,iostat=status) self%rqpara0
  call log_write_check(m_name,s_name,33,status)
  write(iout,*,iostat=status) 'profile_type'
  call log_write_check(m_name,s_name,34,status)
  write(iout,*,iostat=status) self%postype
  call log_write_check(m_name,s_name,35,status)
  write(iout,*,iostat=status) 'npos'
  call log_write_check(m_name,s_name,36,status)
  write(iout,*,iostat=status) self%npos
  call log_write_check(m_name,s_name,37,status)
  if (self%npos>0) then
     write(iout,*,iostat=status) 'profile_positions'
     call log_write_check(m_name,s_name,38,status)
     write(iout,*,iostat=status) self%pos
     call log_write_check(m_name,s_name,39,status)
     write(iout,*,iostat=status) 'profile_values'
     call log_write_check(m_name,s_name,40,status)
     write(iout,*,iostat=status) self%prof
     call log_write_check(m_name,s_name,41,status)
  end if
  write(iout,*,iostat=status) 'nrpams'
  call log_write_check(m_name,s_name,46,status)
  write(iout,*,iostat=status) self%nrpams
  call log_write_check(m_name,s_name,47,status)
  if (self%nrpams>0) then
     write(iout,*,iostat=status) 'real_parameters'
     call log_write_check(m_name,s_name,48,status)
     write(iout,*,iostat=status) self%rpar
     call log_write_check(m_name,s_name,49,status)
  end if
  write(iout,*,iostat=status) 'nipams'
  call log_write_check(m_name,s_name,50,status)
  write(iout,*,iostat=status) self%nipams
  call log_write_check(m_name,s_name,51,status)
  if (self%nipams>0) then
     write(iout,*,iostat=status) 'integer_parameters'
     call log_write_check(m_name,s_name,52,status)
     write(iout,*,iostat=status) self%npar
     call log_write_check(m_name,s_name,53,status)
  end if

end subroutine edgprof_write

!---------------------------------------------------------------------
!> close write file
subroutine edgprof_closewrite

  !! local
  character(*), parameter :: s_name='edgprof_closewrite' !< subroutine name

  !! close file
  close(unit=noutep,iostat=status)
  if(status/=0)then
     !! error closing file
     print '("Fatal error: Unable to close output file, ",a)',outputfile
     call log_error(m_name,s_name,1,error_fatal,'Cannot close output data file')
     stop
  end if

end subroutine edgprof_closewrite

!---------------------------------------------------------------------
!> delete object
subroutine edgprof_delete(self)

  !! arguments
  type(edgprof_t), intent(inout) :: self !< type containing profile parameters
  !! local
  character(*), parameter :: s_name='edgprof_delete' !< subroutine name
   integer(ki4) :: i  !< local variable
   
DO i=1,4
  formula_deallocate: select case (self%formula(i))
  case('samples')
     deallocate(self%pos)
     deallocate(self%prof)

  case('userdefined')
     if (self%nrpams>0) deallocate(self%rpar)
     if (self%nipams>0) deallocate(self%npar)
  case default
  end select formula_deallocate
END DO
end subroutine edgprof_delete
!---------------------------------------------------------------------
!> close file
subroutine edgprof_close

  !! local
  character(*), parameter :: s_name='edgprof_close' !< subroutine name

  !! close file
  close(unit=ninep,iostat=status)
  if(status/=0)then
     !! error closing file
     print '("Fatal error: Unable to close control file, ",a)',controlfile
     call log_error(m_name,s_name,1,error_fatal,'Cannot close control data file')
     stop
  end if

end subroutine edgprof_close
!---------------------------------------------------------------------
!> Calculates custom lambda q
subroutine edgprof_lambda_q

  !! local
  character(*), parameter :: s_name='edgprof_lambda_q' !< subroutine name


end subroutine edgprof_lambda_q

end module edgprof_m
