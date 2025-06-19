! ==============================================================================
! File:        p2_double_corrector.f90
! Author:      Giuseppe Francesco Conte (personal email: giuseppefconte00@gmail.com)
!
! Description: This FORTRAN code implements the evolution of the Schwinger–Dyson
!              equations for the Hamiltonian dynamics of the specialized p=2-spin model.
!
!              The model evolves under conservative dynamics with spherical
!              constraints. This implementation numerically integrates the 
!              correlation and response functions over time using a predictor–
!              corrector scheme.
!
!              Respect to the p>=3 case, this contains also another equation in order to 
!              correctly evolve the correltion function between different replicas, arising
!              in the low temperature regime (T'<T_c).
!
! Dependencies: None (standard Fortran 90)
!
! Compilation:  
!             - For small matrices:
!                  gfortran -g -o p2_double_corrector p2_double_corrector.f90
!             - For large matrices (better RAM handling):
!                  gfortran -g -mcmodel=medium -fPIC -o p2_double_corrector p2_double_corrector.f90
!
! Usage:   ./p2_double_corrector mass diss itemp Jf
!
! Notes:
!   - Time is discretized with step size h.
!   - Initial conditions assume T'.
!   - Lagrange multipliers and Green's functions are evolved dynamically.
!
! Acknowledgments:
!   I thank Lorenzo Correale (SISSA) for helpful discussions and insightful
!   suggestions regarding both theoretical and numerical aspects of the problem.
!
!   For a complete explanation of the theoretical background and numerical 
!   details, please refer to my Master's thesis.
!
! Last update: 2025-06-18
! ==============================================================================

program matrices
  !function to find Nan
  use ieee_arithmetic
  implicit none

  ! Dimension of arrays (n-dimensional vectors) and matrix (n x n)
  integer, parameter :: n = 800
  integer ntot,number
  parameter(ntot=(n*(n+1))/2+n)

  ! plot parameters
  integer jw1,jw2,jw3,jw4,jw5,jw6,jw7,jwn
  integer tw1,tw2,tw3,tw4,tw5,tw6,tw7
  double precision mtw1,mtw2,mtw3,mtw4,mtw5,mtw6,mtw7

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! declaring matrix, arrays, integers and float
  ! RESPECT TO THE p=3 case, we also have the contribution 
  ! from the correlation fucntion between different replicas: 
  ! this is encoded by q_pred, pi_q_pred. However, in the eqaution 
  ! appears only the contribution from the FIRST COLUMN of the
  ! different-replcias correlation function: for this reason we only
  ! need this contrubution
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  real(8), dimension(n+1, n+1) :: r_pred, c_pred, pi_c_pred, pi_r_pred,q_pred,pi_q_pred
  real(8), dimension(n+1, n+1) :: r_corr, c_corr, pi_c_corr, pi_r_corr,q_corr,pi_q_corr
  real(8), dimension(n+1)    :: mu_pred,ci,F3_pred,ci_corr,F3_corr,F3,F1_corr,F1_pred,F1,F2_pred,F2_corr,F2,F4,F4_pred,F4_corr! Vettore mu_pred
  double precision plot_r(1:ntot+1),plot_c(1:ntot+1)
  ! loop varibales
  integer :: i,j,k,l
  ! integers parameters
  integer :: p,pminus
  ! real valued varibales
  double precision J0,J1,h,itemp,mass,coeff0,zf,Jf,diss,q_in,z_in
  ! prediction quanitites
  double precision I_mu_pred,I1_mu_pred,I2_mu_pred
  double precision I1_pred,I3_pred,I5_pred,I7_pred
  ! correction quantities
  double precision I_mu_corr,I1_mu_corr,I2_mu_corr
  double precision I1_corr,I3_corr,I5_corr,I7_corr
 
  character(len=50) :: filename
  character(len=50) :: filen
  character(len=32) :: arg1,arg2,arg3,arg4
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! The following block is used to start the simulation with the following input parameters
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  call getarg(1,arg1)
  read(arg1,*) mass
  call getarg(2,arg2)
  read(arg2,*) diss
  call getarg(3,arg3)
  read(arg3,*) itemp
  call getarg(4,arg4)
  read(arg4,*) Jf
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! we specialize this code to the p=2 case, so we do not selct the p-value in agument of the esecutable
  ! we fix the parameters of the simulation below
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  p=2
  pminus=p-1
  J0=1.d0
  J1=J0
  ! choiche of the discretization step
  h=0.01
  coeff0=h/mass
  z_in=0.d0
  zf=0.d0
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! we define the intial value of the INITIAL overlap.
  ! This depend from the temperature of the inital state configuration
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  q_in=0.d0
  if (itemp .ge. 1) then
    q_in = 0.d0
  elseif (itemp .lt. 1) then
    q_in = 1 - (itemp)/J0
  end if

  ! We define the initial value of the Lagrange multiplier.
  ! By inserting the value of q_in into the following expression
  ! for z_in, we obtain the same analytical result that can be
  ! derived from an asymptotic evaluation in the non-quenched
  ! dynamical regime (i.e., equilibrium dynamics).
  ! In other words, z_f = z_in at equilibrium.

  z_in = itemp + (J1*J0)/(itemp)*(1-(q_in)**2)

  ! we define the aymptotica value for the lagrange multilier 
  ! in the following lines as the varible zf. Those expression
  ! can be recovered analytically
  
  if (itemp .ge. 1) then
    zf = itemp + (J1**2)/(itemp)
  elseif (itemp .lt. 1) then
    zf = 2.d0*J1
  end if

  ! Here we fix the paramter of the printing cicles
  tw1 = 1
  tw2 = n/8 !n/8
  tw3 = n/4 !n/4
  tw4 = 3*n/8 !3*n/8
  tw5 = n/2 !n/2
  tw6 = 5*n/8 !5*n/8
  tw7 = 3*n/4 !3*n/4
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  jw1 = 1
  jw2 = n/8
  jw3 = n/4
  jw4 = 3*n/8
  jw5 = n/2
  jw6 = 5*n/8
  jw7 = 3*n/4
  jwn = n  
  
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Here we open the log file and we write the parameters setted to obtain THE LAST output file
  write(filen,'(a5,f5.3,a2,f4.2,a5,f4.2)') 'itemp',itemp,"_J",Jf,'_diss',diss !'_mass',mass,'_diss',diss,"_x",x(a5,f4.2,a5,f4.1,a5,f4.2,a2,f4.2)
  filename = trim("pspin_param_"//filen)//".log"
  open(18,file=filename,status='unknown',recl=256)

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! We create and open each of the input files 
  filename = trim("tw0_"//filen)//".dat"
  open(unit=10, file=filename,status="unknown")
  filename = trim("tw1_"//filen)//".dat"
  open(unit=11, file=filename,status="unknown")
  filename = trim("tw2_"//filen)//".dat"
  open(unit=12, file=filename,status="unknown")
  filename = trim("tw3_"//filen)//".dat"
  open(unit=13, file=filename,status="unknown")
  filename = trim("tw4_"//filen)//".dat"
  open(unit=14, file=filename,status="unknown")
  filename = trim("tw5_"//filen)//".dat"
  open(unit=15, file=filename,status="unknown")
  filename = trim("tw6_"//filen)//".dat"
  open(unit=16, file=filename,status="unknown")
  filename = trim("tw7_"//filen)//".dat"
  open(unit=17, file=filename,status="unknown")
  !
  !   writing at t fixed moving tw
  !
  filename = trim("jw1_"//filen)//".dat"
  open(unit=21, file=filename,status="unknown")
  filename = trim("jw2_"//filen)//".dat"
  open(unit=22, file=filename,status="unknown")
  filename = trim("jw3_"//filen)//".dat"
  open(unit=23, file=filename,status="unknown")
  filename = trim("jw4_"//filen)//".dat"
  open(unit=24, file=filename,status="unknown")
  filename = trim("jw5_"//filen)//".dat"
  open(unit=25, file=filename,status="unknown")
  filename = trim("jw6_"//filen)//".dat"
  open(unit=26, file=filename,status="unknown")
  filename = trim("jw7_"//filen)//".dat"
  open(unit=27, file=filename,status="unknown")
  filename = trim("jwn_"//filen)//".dat"
  open(unit=28, file=filename,status="unknown")
  !
  !   writing at chi(t,tw) and C(t,tw) in the same file
  !
  filename = trim("chiCw1_"//filen)//".dat"
  open(unit=31, file=filename,status="unknown")
  filename = trim("chiCw2_"//filen)//".dat"
  open(unit=32, file=filename,status="unknown")
  filename = trim("chiCw3_"//filen)//".dat"
  open(unit=33, file=filename,status="unknown")
  filename =trim("chiCw4_"//filen)//".dat"
  open(unit=34, file=filename,status="unknown")
  filename =trim("chiCw5_"//filen)//".dat"
  open(unit=35, file=filename,status="unknown")
  filename = trim("chiCw6_"//filen)//".dat"
  open(unit=36, file=filename,status="unknown")
  filename = trim("chiCw7_"//filen)//".dat"
  open(unit=37, file=filename,status="unknown")
  filename = trim("chiCwn_"//filen)//".dat"
  open(unit=38, file=filename,status="unknown")
  filename = trim("lagrange_"//filen)//".dat"
  open(unit=40, file=filename,status="unknown")
  filename = trim("q-evolution_"//filen)//".dat"
  open(unit=41, file=filename,status="unknown")
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Matrix and vectors initializations to ZERO
  r_pred= 0.0d0
  pi_c_pred = 0.0d0
  pi_r_pred= 0.0d0
  r_corr= 0.0d0
  pi_c_corr = 0.0d0
  pi_r_corr= 0.0d0
  mu_pred = 0.0d0
  c_corr = 0.0d0
  c_pred = 0.d0
  ! we initialize also the matrix designed to the evolution od the correlation function between different replicas
  q_pred=0.d0
  pi_q_pred=0.d0
  q_corr=0.d0
  pi_q_corr=0.d0
  ! We set the equal-time conditions.
  ! The corresponding condition for the correlation matrix between different replicas
  ! has, as initial value, the initial static overlap.
  do i = 1, n
     c_pred(i, i) = 1.0d0
     mu_pred(i)=0
     ci(i)=0.d0
     F3_pred(i)=0.d0
     F3_corr(i)=0.d0
     F3(i)=0.d0
     q_pred(i,i)=q_in
     q_corr(i,i)=q_in
  end do

  ! we set the inital value of the Lagrange multiplier 
  mu_pred(1)=z_in
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! We start the dynamics   
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  do i=1,n
   
    ! Integral for the lagrange multiplier
    I_mu_pred=0.d0
    do k=1,i
          !I1_mu_pred=(pminus)*p*0.5d0*(J0**2)*(c_pred(i,k))**(p-2)*((p-2)*pi_c_pred(i,k)*r_pred(i,k)+pi_r_pred(i,k)*c_pred(i,k))
          !I2_mu_pred=(pminus)*p*0.5d0*(J0**2)*pi_c_pred(i,k)*((c_pred(i,k))**(p-2))*r_pred(i,k)
          !I3_mu_pred=3.d0*p*0.5d0*(J0**2)*(pminus)*r_pred(i,k)*(c_pred(i,k))**(p-2)*pi_c_pred(i,k)
          !I4_mu_pred=3.d0*p*0.5d0*(J0**2)*((c_pred(i,k))**(pminus))*pi_r_pred(i,k)
          !I_mu_pred= I_mu_pred + (I1_mu_pred + I2_mu_pred + I3_mu_pred + I4_mu_pred)
          I1_mu_pred=4.d0*(J0**2)*c_pred(i,k)*pi_r_pred(i,k)
          I2_mu_pred=4.d0*(J0**2)*r_pred(i,k)*pi_c_pred(i,k)
          I_mu_pred = I_mu_pred + (I1_mu_pred + I2_mu_pred)
    enddo
    I_mu_pred = I_mu_pred*h
    ci(i) = (4.d0*J1*J0)/(itemp)*(pi_c_pred(i,1)*c_pred(i,1) - pi_q_pred(i,1)*q_pred(i,1)) ! thise line contains also the contribution due to the
                                                                                          ! correlation fucntion between different replicas
    F3_pred(i) = I_mu_pred + ci(i)
    mu_pred(i+1)= mu_pred(i) + h*F3_pred(i)
    !print*,'i',i,'mu_pred(i)',mu_pred(i)
    write(40,*) h*i,mu_pred(i),zf,mu_pred(i)-zf,q_pred(i,1)
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! start the evolution all over the columns
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    do j=1,i
      ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !!! I1 = int_{t'}^{t}dt'' r(t,t'')c(t,t'')^(p-2)r(t'',t')
      ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
      I1_pred = 0.d0
      do k=j,i
        I1_pred = I1_pred + r_pred(i,k)*r_pred(k,j)
      enddo
      I1_pred = (J0**2)*p*pminus*0.5*I1_pred
      I1_pred = I1_pred*h
      F1_pred(j) = I1_pred
      !print*,'F1_pred(j)',F1_pred(j)

      ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !  I3 = int_0^t' dt'' c(t,t'')**pminus r(t',t'')	
      ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
      I3_pred = 0.d0
      do k =1,j
        I3_pred=I3_pred+c_pred(i,k)*r_pred(j,k)
      enddo
      I3_pred = p*(J0**2)*0.5*I3_pred
      I3_pred = h*I3_pred

      ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! I5 = int_0^t dt'' r(t,t'')c(t,t'')**(p-2)*c(t'',t')	
      ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
      I5_pred = 0.d0
      do k=1,i
        I5_pred=I5_pred+r_pred(i,k)*c_pred(k,j)
      enddo
      I5_pred = (J0**2)*p*pminus*0.5*I5_pred
      I5_pred = h*I5_pred
      F2_pred(j) = I3_pred + I5_pred + ((J1*J0)/itemp) * (c_pred(i,1)*c_pred(j,1) - q_pred(i,1)*q_pred(j,1))! the last terms is used to propagate correctlu the information of 
                                                                                                            ! the inital state configuation, sampled at
                                                                                                            ! temperature T', in order to taking into account the inital contribution 
                                                                                                            ! from both condesed and paramagnetic intial conditions
    
      ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
      ! Here appers the lines used to compute the integral aassociated to the adjuctive equation that appers
      ! in the low temperature regime, in order evelve correctly the correlation function between different replicas 
      ! I7 = int_0^t dt'' r(t,t'')q(t'',t)
      ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

      I7_pred=0.d0
      do k=1,i
            I7_pred = I7_pred + r_pred(i,k)*q_pred(k,1)
      enddo
      I7_pred=I7_pred*h
      I7_pred = I7_pred*J1*J1
      !if (i==j) then
      !print*,i,I7_pred
      !endif
      F4_pred(i) = I7_pred + ((J0*J1)/itemp)*(q_in*c_pred(i,1)+(1-(2*q_in))*q_pred(i,1))
                
      !        
      ! Updating the valued computed before in order to evolve correctly the dynamics
      !  
      if (i==j) then
      pi_r_pred(i+1,j) = pi_r_pred(i,j) + h*(-mu_pred(i)*r_pred(i,j)+F1_pred(j)) + 1
      else
      pi_r_pred(i+1,j) = pi_r_pred(i,j) + h*(-mu_pred(i)*r_pred(i,j)+F1_pred(j))
      endif
      pi_c_pred(i+1,j) = pi_c_pred(i,j) + h*(-mu_pred(i)*c_pred(i,j) + F2_pred(j))

      r_pred(i+1,j) = r_pred(i,j) + coeff0*pi_r_pred(i,j)
      c_pred(i+1,j) = c_pred(i,j) + coeff0*pi_c_pred(i,j)
      c_pred(j,i)=c_pred(i,j)
      q_pred(i,i)=q_in
      c_pred(i,i)=1                  
      !print*,c_pred(i+1,j),c_pred(i,j),r_pred(i+1,j),r_pred(i,j)
  
    enddo
    ! Here we compute the evolution of the correlation function between different replicas (ONLY THE FIRS COLUMN)
    pi_q_pred(i+1,1) = pi_q_pred(i,1) + h*(-mu_pred(i)*q_pred(i,1) + F4_pred(i))
    q_pred(i+1,1) = q_pred(i,1) + coeff0*pi_q_pred(i,1)

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% start FIRST   correction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    I_mu_corr=0.d0
    do k=1,i+1
          !I1_mu_corr=(pminus)*p*0.5d0*(J0**2)*(c_pred(i+1,k))**(p-2)*((p-2)*pi_c_pred(i+1,k)*r_pred(i+1,k)+pi_r_pred(i+1,k)*c_pred(i+1,k))
          !I2_mu_corr=(pminus)*p*0.5d0*(J0**2)*pi_c_pred(i+1,k)*((c_pred(i+1,k))**(p-2))*r_pred(i+1,k)
          !I3_mu_corr=3.d0*p*0.5d0*(J0**2)*(pminus)*r_pred(i+1,k)*(c_pred(i+1,k))**(p-2)*pi_c_pred(i+1,k)
          !I4_mu_corr=3.d0*p*0.5d0*(J0**2)*((c_pred(i+1,k))**(pminus))*pi_r_pred(i+1,k)
          !I_mu_corr= I_mu_corr + (I1_mu_corr + I2_mu_corr + I3_mu_corr + I4_mu_corr)
          I1_mu_corr = 4.d0*(J0**2)*c_pred(i+1,k)*pi_r_pred(i+1,k)
          I2_mu_corr = 4.d0*(J0**2)*r_pred(i+1,k)*pi_c_pred(i+1,k)
          I_mu_corr = I_mu_corr + (I1_mu_corr+I2_mu_corr)
    enddo
    I_mu_corr = I_mu_corr*h
    ci_corr(i) = ((4.d0*(J0*J1))/itemp)*(pi_c_pred(i+1,1)*c_pred(i+1,1)- pi_q_pred(i+1,1)*q_pred(i+1,1))
    F3_corr(i) = I_mu_corr + ci_corr(i)
    F3(i)=0.5d0*(F3_corr(i)+F3_pred(i))
    mu_pred(i+1)= mu_pred(i) + h*F3(i)    
    !print*,'i',i,'F3(i)',F3(i),'F3_pred(i)',F3_pred(i),'F3_corr(i)',F3_corr(i)
    do j=1,i
      I1_corr = 0.d0
      do k=j,i+1
      I1_corr = I1_corr + r_pred(i+1,k)*((c_pred(i+1,k))**(p-2))*r_pred(k,j)
      enddo
      I1_corr = (J0**2)*p*pminus*0.5*I1_corr
      I1_corr = I1_corr*h
      F1_corr(j) = I1_corr
      F1(j)=0.5d0*(F1_corr(j)+F1_pred(j))
            
      ! I3 = int_0^t' dt'' c(t,t'')**pminus r(t',t'')
      I3_corr = 0.d0
      do k =1,j
        I3_corr=I3_corr+(c_pred(i+1,k)**(pminus))*r_pred(j,k)
      enddo
      I3_corr = p*(J0**2)*0.5*I3_corr
      I3_corr = h*I3_corr

      ! I5 = int_0^t dt'' r(t,t'')c(t,t'')**(p-2)*c(t'',t')
      I5_corr = 0.d0
      do k=1,i+1
        I5_corr=I5_corr+r_pred(i+1,k)*(c_pred(i+1,k)**(p-2))*c_pred(k,j)
      enddo
      I5_corr = (J0**2)*p*pminus*0.5*I5_corr
      I5_corr = h*I5_corr
      F2_corr(j)= I3_corr + I5_corr + (J1*J0)/itemp*(c_pred(i+1,1)*c_pred(j,1) - q_pred(i+1,1)*q_pred(j,1))  
      F2(j) = 0.5d0*(F2_pred(j) + F2_corr(j))
          
      ! I7 = int_0^t dt'' r(t,t'')q(t'',t0)
      I7_corr=0.d0
      do k=1,i+1
        I7_corr = I7_corr + r_pred(i+1,k)*q_pred(k,1)
      enddo
      I7_corr=I7_corr*h
      F4_corr(i) = I7_corr + ((J0**2)/itemp)*(q_in*c_pred(i+1,1)+(1-2*q_in)*q_pred(i+1,1))
      F4(i) = 0.5d0*(F4_pred(i)+F4_corr(i))
      
      if (i==j) then
      pi_r_pred(i+1,j) = pi_r_pred(i,j) + h*(0.5d0*(-mu_pred(i)*r_pred(i,j) -mu_pred(i+1)*r_pred(i+1,j)) + F1(j)) + 1
      else
      pi_r_pred(i+1,j) = pi_r_pred(i,j) + h*(0.5d0*(-mu_pred(i)*r_pred(i,j) -mu_pred(i+1)*r_pred(i+1,j)) + F1(j))
      endif
      pi_c_pred(i+1,j) = pi_c_pred(i,j) + h*(0.5d0*(-mu_pred(i)*c_pred(i,j) -mu_pred(i+1)*c_pred(i+1,j)) + F2(j))

      r_pred(i+1,j) = r_pred(i,j) + coeff0*0.5d0*(pi_r_pred(i,j)+pi_r_pred(i+1,j))
      c_pred(i+1,j) = c_pred(i,j) + coeff0*0.5d0*(pi_c_pred(i,j)+pi_c_pred(i+1,j))

      c_pred(j,i)=c_pred(i,j)
      q_pred(j,i)=q_pred(i,j)
      c_pred(i,i)=1
      q_pred(i,i)=q_in
      !print*,'j',j,'F1(j)',F1(j),'F1_pred(j)',F1_pred(j),'F1_corr(j)',F1_corr(j)
    enddo ! end su j
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! correction step over q_pred
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pi_q_pred(i+1,1) = pi_q_pred(i,1) + h*(0.5d0*(-mu_pred(i)*q_pred(i,1) -mu_pred(i+1)*q_pred(i+1,1)) + F4(i))
    q_pred(i+1,1) = q_pred(i,1) + coeff0*0.5d0*(pi_q_pred(i,1) + pi_q_pred(i+1,1))
      
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SECOND correction step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    I_mu_corr=0.d0
    do k=1,i+1
          !I1_mu_corr=(pminus)*p*0.5d0*(J0**2)*(c_pred(i+1,k))**(p-2)*((p-2)*pi_c_pred(i+1,k)*r_pred(i+1,k)+pi_r_pred(i+1,k)*c_pred(i+1,k))
          !I2_mu_corr=(pminus)*p*0.5d0*(J0**2)*pi_c_pred(i+1,k)*((c_pred(i+1,k))**(p-2))*r_pred(i+1,k)
          !I3_mu_corr=3.d0*p*0.5d0*(J0**2)*(pminus)*r_pred(i+1,k)*(c_pred(i+1,k))**(p-2)*pi_c_pred(i+1,k)
          !I4_mu_corr=3.d0*p*0.5d0*(J0**2)*((c_pred(i+1,k))**(pminus))*pi_r_pred(i+1,k)
          !I_mu_corr= I_mu_corr + (I1_mu_corr + I2_mu_corr + I3_mu_corr + I4_mu_corr)
          I1_mu_corr = 4.d0*(J0**2)*c_pred(i+1,k)*pi_r_pred(i+1,k)
          I2_mu_corr = 4.d0*(J0**2)*r_pred(i+1,k)*pi_c_pred(i+1,k)
          I_mu_corr = I_mu_corr + (I1_mu_corr+I2_mu_corr)
    enddo
    I_mu_corr = I_mu_corr*h
    ci_corr(i) = ((4.d0*(J0*J1))/itemp)*(pi_c_pred(i+1,1)*c_pred(i+1,1)- pi_q_pred(i+1,1)*q_pred(i+1,1))
    F3_corr(i) = I_mu_corr + ci_corr(i)
    F3(i)=0.5d0*(F3_corr(i)+F3_pred(i))
    mu_pred(i+1)= mu_pred(i) + h*F3(i)
    !print*,'i',i,'F3(i)',F3(i),'F3_pred(i)',F3_pred(i),'F3_corr(i)',F3_corr(i)
    

    do j=1,i
    
      I1_corr = 0.d0
      do k=j,i+1
        I1_corr = I1_corr + r_pred(i+1,k)*((c_pred(i+1,k))**(p-2))*r_pred(k,j)
      enddo
      I1_corr = (J0**2)*p*pminus*0.5*I1_corr
      I1_corr = I1_corr*h
      F1_corr(j) = I1_corr
      F1(j)=0.5d0*(F1_corr(j)+F1_pred(j))
          
      !  I3 = int_0^t' dt'' c(t,t'')**pminus r(t',t'')
      I3_corr = 0.d0
      do k =1,j
        I3_corr=I3_corr+(c_pred(i+1,k)**(pminus))*r_pred(j,k)
      enddo
      I3_corr = p*(J0**2)*0.5*I3_corr
      I3_corr = h*I3_corr

      ! I5 = int_0^t dt'' r(t,t'')c(t,t'')**(p-2)*c(t'',t')
      I5_corr = 0.d0
      do k=1,i+1
          I5_corr=I5_corr+r_pred(i+1,k)*(c_pred(i+1,k)**(p-2))*c_pred(k,j)
      enddo
      I5_corr = (J0**2)*p*pminus*0.5*I5_corr
      I5_corr = h*I5_corr
      F2_corr(j)= I3_corr + I5_corr + ((J1*J0)/itemp)*(c_pred(i+1,1)*c_pred(j,1) - q_pred(i+1,1)*q_pred(j,1))  
      F2(j) = 0.5d0*(F2_pred(j) + F2_corr(j))
          
      ! I7 = int_0^t dt'' r(t,t'')q(t'',t0)
      I7_corr=0.d0
      do k=1,i+1
            I7_corr = I7_corr + r_pred(i+1,k)*q_pred(k,1)
      enddo
      I7_corr=I7_corr*h
      F4_corr(i) = I7_corr + ((J0**2)/itemp)*(q_in*c_pred(i+1,1)+(1.d0-(2.d0*q_in))*q_pred(i+1,1))
      F4(i) = 0.5d0*(F4_pred(i)+F4_corr(i))  
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
      !update     
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if (i==j) then
      pi_r_pred(i+1,j) = pi_r_pred(i,j) + h*(0.5d0*(-mu_pred(i)*r_pred(i,j) -mu_pred(i+1)*r_pred(i+1,j)) + F1(j)) + 1
      else
      pi_r_pred(i+1,j) = pi_r_pred(i,j) + h*(0.5d0*(-mu_pred(i)*r_pred(i,j) -mu_pred(i+1)*r_pred(i+1,j)) + F1(j))
      endif
      pi_c_pred(i+1,j) = pi_c_pred(i,j) + h*(0.5d0*(-mu_pred(i)*c_pred(i,j) -mu_pred(i+1)*c_pred(i+1,j)) + F2(j))

      r_pred(i+1,j) = r_pred(i,j) + coeff0*0.5d0*(pi_r_pred(i,j)+pi_r_pred(i+1,j))
      c_pred(i+1,j) = c_pred(i,j) + coeff0*0.5d0*(pi_c_pred(i,j)+pi_c_pred(i+1,j))

      c_pred(j,i)=c_pred(i,j)
      q_pred(j,i)=q_pred(i,j)
      c_pred(i,i)=1
      q_pred(i,i)=q_in
    enddo ! end over j
    pi_q_pred(i+1,1) = pi_q_pred(i,1) + h*(0.5d0*(-mu_pred(i)*q_pred(i,1) -mu_pred(i+1)*q_pred(i+1,1) ) + F4(i))
    q_pred(i+1,1) = q_pred(i,1) + coeff0*0.5d0*(pi_q_pred(i,1) + pi_q_pred(i+1,1))
    !print*,'q',q_pred(i,1)
    write(41,*) h*i,q_pred(i,1),q_in
  enddo
  ! end evolution
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! print the 10 x 10 block the evolution matrices AFTER the Predictor-Corrector scheme
    print *, "===== MATRICE r_pred (parte inferiore) ====="
    do i = 1, 10
    write(*,'(1X, "( ", *(F10.5,1X))') (r_pred(i,j), j=1,i)
    end do
    
    print *, "===== MATRICE pi_r_pred (parte inferiore) ====="
    do i = 1, 10
      write(*,'(1X, "( ", *(F10.5,1X))') (pi_r_pred(i,j), j=1,i)
    end do

    print *, "===== MATRICE c_pred (parte inferiore) ====="
    do i = 1, 10
      write(*,'(1X, "( ", *(F10.5,1X))') (c_pred(i,j), j=1,i)
    end do
    print *, "===== MATRICE pi_c_pred (parte inferiore) ====="
    do i = 1, 10
      write(*,'(1X, "( ", *(F10.5,1X))') (pi_c_pred(i,j), j=1,i)
    end do
  
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! now we VECTORIZE the matrix in order to use the same writing algoritm used in tartaglia-nessi codes

  do i=1,n
    do j=1,i
      plot_c(i*(i+1)/2+j)=c_pred(i,j)
      plot_r(i*(i+1)/2+j)=r_pred(i,j)
    enddo
   
    if(i.eq.jw1) then
    number = (i*(i+1)/2)
    do j=1,i
      write(21,*) (i-j)*h,plot_c(number+j),plot_r(number+j),mu_pred(i)
    enddo
    endif
    !       ----------  i=jw2 --------------
    if(i.eq.jw2) then
    number = (i*(i+1)/2)
      do j=1,i
          write(22,*) (i-j)*h,plot_c(number+j),plot_r(number+j),mu_pred(i)
      enddo
    endif
    !       ----------  i=jw3 --------------
    if(i.eq.jw3) then
      number = (i*(i+1)/2)
      do j=1,i
        write(23,*) (i-j)*h,plot_c(number+j),plot_r(number+j),mu_pred(i)
      enddo
    endif
    !       ----------  i=jw4 --------------
    if(i.eq.jw4) then
      number = (i*(i+1)/2)
      do j=1,i
        write(24,*) (i-j)*h,plot_c(number+j),plot_r(number+j),mu_pred(i)
      enddo
    endif
    !       ----------  i=jw5 --------------
    if(i.eq.jw5) then
      number = (i*(i+1)/2)
      do j=1,i
        write(25,*) (i-j)*h,plot_c(number+j),plot_r(number+j),mu_pred(i)
      enddo
    endif
    !       ----------  i=jw6 --------------
    if(i.eq.jw6) then
      number = (i*(i+1)/2)
      do j=1,i
        write(26,*) (i-j)*h,plot_c(number+j),plot_r(number+j),mu_pred(i)
      enddo
    endif
    !       ----------  i=jw7 --------------
    if(i.eq.jw7) then
      number = (i*(i+1)/2)
      do j=1,i
        write(27,*) (i-j)*h,plot_c(number+j),plot_r(number+j),mu_pred(i)
      enddo
    endif
    !       ----------  i=jwn --------------
    if(i.eq.jwn) then
      number = (i*(i+1)/2)
      do j=1,i
        write(28,*) (i-j)*h,plot_c(number+j),plot_r(number+j),mu_pred(i)
      enddo
    endif
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !        OUTPUT writing results at t' fixed
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    number = (i*(i+1)/2)
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !        Writing c(tau+tw,t_w) and c(tau+tw,t_w vs tau for some tw's
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !
      if(i.ge.0) then
          write(10,*) i*h,plot_c(number),plot_r(number)
      endif
      if(i.ge.tw1) then
          write(11,*) (i-tw1)*h, plot_c(number+tw1),plot_r(number+tw1)
      endif
      if(i.ge.tw2) then
          write(12,*) (i-tw2)*h, plot_c(number+tw2),plot_r(number+tw2)
      endif
      if(i.ge.tw3) then
          write(13,*) (i-tw3)*h, plot_c(number+tw3),plot_r(number+tw3)
      endif
      if(i.ge.tw4) then
          write(14,*) (i-tw4)*h, plot_c(number+tw4),plot_r(number+tw4)
      endif
      if(i.ge.tw5) then
          write(15,*) (i-tw5)*h, plot_c(number+tw5),plot_r(number+tw5)
      endif
      if(i.ge.tw6) then
          write(16,*) (i-tw6)*h, plot_c(number+tw6),plot_r(number+tw6)
      endif
      if(i.ge.tw7) then
          write(17,*) (i-tw7)*h, plot_c(number+tw7),plot_r(number+tw7)
      endif
    !
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !                   end printing dynamics
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
  enddo
 
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !  We calculate the integrated response ad differetn waiting time 
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !
  ! -----------------tw1-------------------------
  !
        mtw1=0.
      write(31,*) 0.,1.,0.
      do k=1,n-tw1
        mtw1=.5*plot_r(tw1*(tw1+1)/2+tw1)
        do l=tw1+1,tw1+k-1
          mtw1=mtw1+plot_r((tw1+k)*(tw1+k+1)/2+l)
        enddo
        write(31,*) (tw1+k)*h,plot_c((tw1+k)*(tw1+k+1)/2+tw1), h*mtw1,(1/itemp)*(1-plot_c((tw1+k)*(tw1+k+1)/2+tw1))
      enddo
  !
  ! -----------------tw2-------------------------
  !
      mtw2=0.
      write(32,*) 0.,1.,0.
      do k=1,n-tw2
        mtw2=.5*plot_r(tw2*(tw2+1)/2+tw2)
        do l=tw2+1,tw2+k-1
          mtw2=mtw2+plot_r((tw2+k)*(tw2+k+1)/2+l)
        enddo
        write(32,*) (tw2+k)*h,plot_c((tw2+k)*(tw2+k+1)/2+tw2),h*mtw2,(1/itemp)*(1-plot_c((tw2+k)*(tw2+k+1)/2+tw2))
      enddo
  !
  ! -----------------tw3-------------------------
  !
      mtw3=0.
      write(33,*) 0.,1.,0.
      do k=1,n-tw3
        mtw3=.5*plot_r(tw3*(tw3+1)/2+tw3)
        do l=tw3+1,tw3+k-1
          mtw3=mtw3+plot_r((tw3+k)*(tw3+k+1)/2+l)
        enddo
        write(33,*) (tw3+k)*h,plot_c((tw3+k)*(tw3+k+1)/2+tw3),h*mtw3,(1/itemp)*(1-plot_c((tw3+k)*(tw3+k+1)/2+tw3))
      enddo
  !
  ! -----------------tw4-------------------------
  !
      mtw4=0.
      write(34,*) 0.,1.,0.
      do k=1,n-tw4
        mtw4=.5*plot_r(tw4*(tw4+1)/2+tw4)
        do l=tw4+1,tw4+k-1
          mtw4=mtw4+plot_r((tw4+k)*(tw4+k+1)/2+l)
        enddo
        write(34,*) (tw4+k)*h,plot_c((tw4+k)*(tw4+k+1)/2+tw4),h*mtw4,(1/itemp)*(1-plot_c((tw4+k)*(tw4+k+1)/2+tw4))
      enddo
  !
  ! -----------------tw5-------------------------
  !
      mtw5=0.
      write(35,*) 0.,1.,0.
      do k=1,n-tw5
        mtw5=.5*plot_r(tw5*(tw5+1)/2+tw5)
        do l=tw5+1,tw5+k-1
          mtw5=mtw5+plot_r((tw5+k)*(tw5+k+1)/2+l)
        enddo
        write(35,*) (tw5+k)*h,plot_c((tw5+k)*(tw5+k+1)/2+tw5),h*mtw5,(1/itemp)*(1-plot_c((tw5+k)*(tw5+k+1)/2+tw5))
      enddo
  !
  ! -----------------tw6-------------------------
  !
      mtw6=0.
      write(36,*) 0.,1.,0.
      do k=1,n-tw6
        mtw6=.5*plot_r(tw6*(tw6+1)/2+tw6)
        do l=tw6+1,tw6+k-1
          mtw6=mtw6+plot_r((tw6+k)*(tw6+k+1)/2+l)
        enddo
        write(36,*) (tw6+k)*h,plot_c((tw6+k)*(tw6+k+1)/2+tw6),h*mtw6,(1/itemp)*(1-plot_c((tw6+k)*(tw6+k+1)/2+tw6))
      enddo
  !
  ! -----------------tw7-------------------------
  !
      mtw7=0.
      write(37,*) 0.,1.,0.
      do k=1,n-tw7
        mtw7=.5*plot_r(tw7*(tw7+1)/2+tw7)
        do l=tw7+1,tw7+k-1
          mtw7=mtw7+plot_r((tw7+k)*(tw7+k+1)/2+l)
        enddo
        write(37,*) (tw7+k)*h,plot_c((tw7+k)*(tw7+k+1)/2+tw7),h*mtw7,(1/itemp)*(1-plot_c((tw7+k)*(tw7+k+1)/2+tw7))
      enddo
 
end program matrices