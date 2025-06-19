program matrices
  implicit none

  ! Dichiarazione delle dimensioni
  integer, parameter :: n = 31000
  integer ntot,number
  parameter(ntot=(n*(n+1))/2+n)

!dichiarazione paramentri plot
  integer jw1,jw2,jw3,jw4,jw5,jw6,jw7,jwn
  integer tw1,tw2,tw3,tw4,tw5,tw6,tw7
  double precision mtw1,mtw2,mtw3,mtw4,mtw5,mtw6,mtw7

  ! Dichiarazione delle matrici in doppia precisione
  real(8), dimension(n+1, n+1) :: r_pred, c_pred, pi_c_pred, pi_r_pred, q_pred,pi_q_pred,F1_pred,F2_pred
  !real(8), dimension(n+1, n+1) :: r_corr, c_corr, pi_c_corr, pi_r_corr,q_corr,pi_q_corr,F1_corr,F1_pred,F1,F2_pred,F2_corr,F2
  real(8), dimension(n+1)    :: mu_pred,mu_corr,ci,F3_pred,ci_corr,F3_corr,F3 ! Vettore mu_pred
  double precision plot_r(0:ntot),plot_c(0:ntot)
  double precision I_mu_pred_vec(0:ntot), I_mu_corr_vec(0:ntot), F4_pred(0:ntot),F4_corr(0:ntot),F4(0:ntot)
  ! Variabili di loop
  integer :: i,j,k,l
  !parametri interi
  integer :: p,pminus
  double precision J0,h,itemp,mass,coeff0,zf,Jf,diss,z_in,q_in
  double precision I_mu_pred,I1_mu_pred,I2_mu_pred,I3_mu_pred,I4_mu_pred
  double precision I1_pred,I3_pred,I5_pred,I7_pred
  !correction
  double precision I_mu_corr,I1_mu_corr,I2_mu_corr,I3_mu_corr,I4_mu_corr
  double precision I1_corr,I3_corr,I5_corr,I7_corr
  
  character(len=50) :: filename
  character(len=50) :: filen 
  character(len=32) :: arg1,arg2,arg3,arg4,arg5
  
  !inizializzo i parametri della simulazione
  
  
      	    call getarg(1,arg1)
    	    read(arg1,*) mass
            call getarg(2,arg2)
    	    read(arg2,*) diss
	    call getarg(3,arg3)
    	    read(arg3,*) itemp
	    call getarg(4,arg4)
    	    read(arg4,*) Jf
    	    call getarg(5,arg5)
    	    read(arg5,*) p
  
  
  !p=3
  pminus=p-1
  !mass=1.d0
  J0=1.d0
  !itemp=0.8
  !scelta passo discretizzazione
  h=0.001d0
  coeff0=h/mass
  q_in=0.d0
  z_in=0.d0
  zf=0.d0

 
 
   
          tw1 = 1 
          tw2 = n/8 !n/8
          tw3 = n/4 !n/4
          tw4 = 3*n/8 !3*n/8
          tw5 = n/2 !n/2
          tw6 = 5*n/8 !5*n/8
          tw7 = 3*n/4 !3*n/4
!
          jw1 = 1
          jw2 = n/8
          jw3 = n/4
          jw4 = 3*n/8
          jw5 = n/2
          jw6 = 5*n/8
          jw7 = 3*n/4
          jwn = n  
          
        write(filen,'(a5,f5.3,a2,f4.2,a5,f4.2)') 'itemp',itemp,"_J",Jf,'_diss',diss
!'_mass',mass,'_diss',diss,"_x",x(a5,f4.2,a5,f4.1,a5,f4.2,a2,f4.2)
  	filename = trim("pspin_param_"//filen)//".log"
  	open(18,file=filename,status='unknown',recl=256)
          
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
!   	writing at t fixed moving tw
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
!   	writing at chi(t,tw) and C(t,tw) in the same file
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
        filename = trim("integrals_"//filen)//".dat"
        open(unit=41, file=filename,status="unknown")
        filename = trim("evolution_q_"//filen)//".dat"
        open(unit=42, file=filename,status="unknown")        
        filename = trim("error_"//filen)//".dat"
        open(unit=43, file=filename,status="unknown")
        filename = trim("more_integrals_"//filen)//".dat"
        open(unit=44, file=filename,status="unknown")        
        
  ! Inizializza le matrici a zero
  r_pred= 0.0d0
  pi_c_pred = 0.0d0
  pi_r_pred= 0.0d0
  !r_corr= 0.0d0
  !pi_c_corr = 0.0d0
  !pi_r_corr= 0.0d0
  ! Inizializza mu_pred a zero
  mu_pred = 0.0d0
  !parte per temp basse
  !pi_q_corr=0.d0
  pi_q_pred=0.d0

  ! C a tempoi uguali è l'identità
  !c_corr = 0.0d0
  c_pred = 0.d0
  
  
  if(p==2) then
  	if (itemp .ge. J0) then
  	q_in = 0.d0
    	zf = itemp + (J0**2)/(itemp)
  	elseif (itemp .lt. J0) then
  	q_in = 1.0 - (itemp)/J0
    	zf = 2.d0*J0
  	end if
  	z_in = itemp + (Jf*J0*p)/(2*itemp)*(1-(q_in)**2)
  else
  	q_in=0.d0
  	z_in = itemp + (Jf*J0*p)/(2*itemp)
  	zf=z_in
  endif	
  
  do i = 1, n
     c_pred(i, i) = 1.0d0
     q_pred(i,i) = q_in
     !q_corr(i,i) = q_in
     mu_pred(i)=0.d0
     !mu_corr(i)=0.d0
     ci(i)=0.d0
     F3_pred(i)=0.d0
     !F3_corr(i)=0.d0
     F3(i)=0.d0
     I_mu_pred_vec(i)=0.d0
  end do
  
  pi_r_pred(1,1) = 1.d0
  r_pred(1,1) = h/mass
  pi_c_pred(1,1) = -h*itemp/mass
  print*, q_pred(1,1)
  mu_pred(1)=z_in
  
  do i=1,n
   
  !integrale per il moltiplicatore
  I_mu_pred=0.d0
  do k=1,i
        !I1_mu_pred=(pminus)*p*0.5d0*(J0**2)*(c_pred(i,k))**(p-2)*((p-2)*pi_c_pred(i,k)*r_pred(i,k)+pi_r_pred(i,k)*c_pred(i,k)) 
        !I2_mu_pred=(pminus)*p*0.5d0*(J0**2)*pi_c_pred(i,k)*((c_pred(i,k))**(p-2))*r_pred(i,k)
        !I3_mu_pred=3.d0*p*0.5d0*(J0**2)*(pminus)*r_pred(i,k)*(c_pred(i,k))**(p-2)*pi_c_pred(i,k)
        !I4_mu_pred=3.d0*p*0.5d0*(J0**2)*((c_pred(i,k))**(pminus))*pi_r_pred(i,k)
        !I_mu_pred= I_mu_pred + (I1_mu_pred + I2_mu_pred + I3_mu_pred + I4_mu_pred)
        I_mu_pred = I_mu_pred + 0.5d0*(J0**2)*p*(p+2)*(c_pred(i,k))**(p-2)*((pminus*pi_c_pred(i,k)*r_pred(i,k)) +(c_pred(i,k)*pi_r_pred(i,k)))
  enddo
  I_mu_pred = I_mu_pred*h
  I_mu_pred_vec(i) = I_mu_pred
  ci(i) = ((0.5*p*(p+2)*(J0**2))/itemp)*(pi_c_pred(i,1)*(c_pred(i,1))**(pminus)) - 4.d0*(J0**2)/itemp*(pi_q_pred(i,1)*q_pred(i,1))
  F3_pred(i) = I_mu_pred_vec(i) + ci(i)
  
  if (i==1) then
  mu_pred(i+1)= mu_pred(i) + h*F3_pred(i)
  elseif (i==2) then
  mu_pred(i+1)= mu_pred(i) + h*F3_pred(i)
  else
  mu_pred(i+1)=mu_pred(i-1) + 2.d0*h*F3_pred(i)
  endif
  
  print*,'i',i,'mu_pred(i)',mu_pred(i)
  write(40,*) h*i,mu_pred(i),zf,mu_pred(i)-zf 
  write(41,*) h*i,I_mu_pred_vec(i),ci(i),F3_pred(i),mu_pred(i+1)/mu_pred(i)
  if(i>2) then
    	write(43,*) h*i,mu_pred(i)-mu_pred(i-1)+0.5d0*h*(F3_pred(i)+F3_pred(i-1))
  else
  	write(43,*) 0,0.d0
  endif
      ! I7 = int_0^t dt'' r(t,t'')q(t'',t0)
        I7_pred=0.d0
        do k=1,i
                I7_pred = I7_pred + r_pred(i,k)*q_pred(k,1)
        enddo
        I7_pred=I7_pred*h
        F4_pred(i) = I7_pred + ((J0**2)/itemp)*(q_in*c_pred(i,1)+(1-2*q_in)*q_pred(i,1))
        
          
	  if (i==1) then
		pi_q_pred(i+1,1) = pi_q_pred(i,1) + h*(-mu_pred(i)*q_pred(i,1) + F4_pred(i))
        	q_pred(i+1,1) = q_pred(i,1) + coeff0*pi_q_pred(i,1)
	  elseif (i==2) then
	  	pi_q_pred(i+1,1) = pi_q_pred(i,1) + h*(-mu_pred(i)*q_pred(i,1) + F4_pred(i))
        	q_pred(i+1,1) = q_pred(i,1) + coeff0*pi_q_pred(i,1)
	  else
	        pi_q_pred(i+1,1) = pi_q_pred(i-1,1) + 2.d0*h*(-mu_pred(i)*q_pred(i,1) + F4_pred(i))
        	q_pred(i+1,1) = q_pred(i-1,1) + 2.d0*coeff0*pi_q_pred(i,1)
	  endif
        
  write(42,*) i*h,q_pred(i+1,1),q_pred(i+1,1)-q_in
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  do j=1,i
 
!!! I1 = int_{t'}^{t}dt'' r(t,t'')c(t,t'')^(p-2)r(t'',t')
 
  	I1_pred = 0.d0
	do k=j,i
		I1_pred = I1_pred + r_pred(i,k)*((c_pred(i,k))**(p-2))*r_pred(k,j)
	enddo
	I1_pred = (J0**2)*p*pminus*0.5*I1_pred
	I1_pred = I1_pred*h
	F1_pred(i,j) = I1_pred
	!print*,'F1_pred(j)',F1_pred(j)
	
!  I3 = int_0^t' dt'' c(t,t'')**pminus r(t',t'')	
	I3_pred = 0.d0
	do k =1,j
		I3_pred=I3_pred+(c_pred(i,k)**(pminus))*r_pred(j,k)
	enddo
	I3_pred = p*(J0**2)*0.5*I3_pred
	I3_pred = h*I3_pred
	
! I5 = int_0^t dt'' r(t,t'')c(t,t'')**(p-2)*c(t'',t')	
	I5_pred = 0.d0
        do k=1,i
              if(k.le.j) then
              	I5_pred=I5_pred+r_pred(i,k)*(c_pred(i,k)**(p-2))*c_pred(j,k)
              else
              	I5_pred=I5_pred+r_pred(i,k)*(c_pred(i,k)**(p-2))*c_pred(k,j)
              endif
        enddo
        I5_pred = (J0**2)*p*pminus*0.5*I5_pred
	I5_pred = h*I5_pred		
	F2_pred(i,j)= I3_pred + I5_pred + ((0.5*p*(J0**2))/itemp)*((c_pred(i,1)**(pminus))*c_pred(j,1) - q_pred(i,1)*q_pred(j,1)) !ulrimo termine è legto ale condizioni iniziali nell'integrale di C
	
	if(j == tw4) then
		write(44,*) i*h,I3_pred,I5_pred,F2_pred(i,j)
	endif	
  

	    if (i == j) then
	      ! impulso iniziale
	      pi_r_pred(i+1,j) = pi_r_pred(i,j) + h*(-mu_pred(i)*r_pred(i,j) + F1_pred(i,j)) + 1
	      pi_c_pred(i+1,j) = pi_c_pred(i,j) + h*(-mu_pred(i)*c_pred(i,j) + F2_pred(i,j))
	      r_pred(i+1,j) = r_pred(i,j) + coeff0*pi_r_pred(i,j)
	      c_pred(i+1,j) = c_pred(i,j) + coeff0*pi_c_pred(i,j)
	      pi_c_pred(j,i+1) = pi_c_pred(i+1,j)
	      c_pred(j,i+1) = c_pred(i+1,j)
	      c_pred(i,i) = 1.d0
	    elseif (i == j+1) then
	      ! secondo step esplicito
	      pi_r_pred(i+1,j) = pi_r_pred(i,j) + h*(-mu_pred(i)*r_pred(i,j) + F1_pred(i,j))
	      pi_c_pred(i+1,j) = pi_c_pred(i,j) + h*(-mu_pred(i)*c_pred(i,j) + F2_pred(i,j))
	      r_pred(i+1,j) = r_pred(i,j) + coeff0*pi_r_pred(i,j)
	      c_pred(i+1,j) = c_pred(i,j) + coeff0*pi_c_pred(i,j)
	      pi_c_pred(j,i+1) = pi_c_pred(i+1,j)
	      c_pred(j,i+1) = c_pred(i+1,j)
	      c_pred(i,i) = 1.d0
	    else
	      ! schema a due passi
	      pi_r_pred(i+1,j) = pi_r_pred(i-1,j) + 2.d0*h*(-mu_pred(i)*r_pred(i,j) + F1_pred(i,j))
	      pi_c_pred(i+1,j) = pi_c_pred(i-1,j) + 2.d0*h*(-mu_pred(i)*c_pred(i,j) + F2_pred(i,j))
	      r_pred(i+1,j) = r_pred(i-1,j) + 2.d0*coeff0*pi_r_pred(i,j)
	      c_pred(i+1,j) = c_pred(i-1,j) + 2.d0*coeff0*pi_c_pred(i,j)
	      pi_c_pred(j,i+1) = pi_c_pred(i+1,j)
	      c_pred(j,i+1) = c_pred(i+1,j)
	      c_pred(i,i) = 1.d0
	    endif
	    !print*,'pi_r_pred(',i,j,')', pi_r_pred(i,j), 'pi_c_pred(',i,j,')', pi_c_pred(i,j)
	    !print*,'r_pred(',i,j,')', r_pred(i,j), 'c_pred(',i,j,')', c_pred(i,j)
	    ! comuni a tutti i casi
	    !r_pred(i+1,j) = r_pred(i,j) + coeff0*pi_r_pred(i,j)
	    !c_pred(i+1,j) = c_pred(i,j) + coeff0*pi_c_pred(i,j)
	    !c_pred(j,i+1) = c_pred(i+1,j)

   if (c_pred(i,j) > 1 ) then
   	print*,' errore'
   	!stop
   endif
  enddo
	
 enddo
 
 	print *, "===== MATRICE r_pred (parte inferiore) ====="
	do i = 1, n
	  write(*,'(1X, "( ", *(F10.5,1X))') (r_pred(i,j), j=1,n)
	end do
	
	print *, "===== MATRICE pi_r_pred (parte inferiore) ====="
	do i = 1, n
	  write(*,'(1X, "( ", *(F10.5,1X))') (pi_r_pred(i,j), j=1,n)
	end do

	print *, "===== MATRICE c_pred (parte inferiore) ====="
	do i = 1, n
	  write(*,'(1X, "( ", *(F10.5,1X))') (c_pred(i,j), j=1,n)
	end do
 	print *, "===== MATRICE pi_c_pred (parte inferiore) ====="
	do i = 1, n
	  write(*,'(1X, "( ", *(F10.5,1X))') (pi_c_pred(i,j), j=1,n)
	end do
 
 
  do i=1,n
   do j=1,i
   
    plot_c(i*(i+1)/2+j)=c_pred(i,j)
    plot_r(i*(i+1)/2+j)=r_pred(i,j)
    enddo
    
             if(i.eq.jw1) then
           number = (i*(i+1)/2)
           do j=0,i
            write(21,*) (i-j)*h,plot_c(number+j),plot_r(number+j),mu_pred(i)
           enddo
         endif
!       ----------  i=jw2 --------------
         if(i.eq.jw2) then
           number = (i*(i+1)/2)
           do j=0,i
            write(22,*) (i-j)*h,plot_c(number+j),plot_r(number+j),mu_pred(i)
           enddo
         endif
!       ----------  i=jw3 --------------
         if(i.eq.jw3) then
           number = (i*(i+1)/2)
           do j=0,i
            write(23,*) (i-j)*h,plot_c(number+j),plot_r(number+j),mu_pred(i)
           enddo
         endif
!       ----------  i=jw4 --------------
         if(i.eq.jw4) then
           number = (i*(i+1)/2)
           do j=0,i
            write(24,*) (i-j)*h,plot_c(number+j),plot_r(number+j),mu_pred(i)
           enddo
         endif
!       ----------  i=jw5 --------------
         if(i.eq.jw5) then
           number = (i*(i+1)/2)
           do j=0,i
            write(25,*) (i-j)*h,plot_c(number+j),plot_r(number+j),mu_pred(i)
           enddo
         endif
!       ----------  i=jw6 --------------
         if(i.eq.jw6) then
           number = (i*(i+1)/2)
           do j=0,i
            write(26,*) (i-j)*h,plot_c(number+j),plot_r(number+j),mu_pred(i)
           enddo
         endif
!       ----------  i=jw7 --------------
         if(i.eq.jw7) then
           number = (i*(i+1)/2)
           do j=0,i
            write(27,*) (i-j)*h,plot_c(number+j),plot_r(number+j),mu_pred(i)
           enddo
         endif
!       ----------  i=jwn --------------
         if(i.eq.jwn) then
           number = (i*(i+1)/2)
           do j=0,i
            write(28,*) (i-j)*h,plot_c(number+j),plot_r(number+j),mu_pred(i)
           enddo
         endif
!
!     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!        OUTPUT Ecriture des resultats
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
         number = (i*(i+1)/2)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!        Writing c(tau+tw,t_w) and c(tau+tw,t_w vs tau for some tw's
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
         if(i.ge.1) then
            write(10,*) i*h,plot_c(number+1),plot_r(number+1)
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
!                   FIN DE LA DINAMICA
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  enddo
  
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!                   CALCULO DE LAS MAGN
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! -----------------tw1-------------------------
!
      mtw1=0.
      write(31,*) 0.,1.,0.
      do k=0,n-tw1
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
      do k=0,n-tw2
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
      do k=0,n-tw3
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
      do k=0,n-tw4
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
      do k=0,n-tw5
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
      do k=0,n-tw6
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
      do k=0,n-tw7
        mtw7=.5*plot_r(tw7*(tw7+1)/2+tw7)
        do l=tw7+1,tw7+k-1
          mtw7=mtw7+plot_r((tw7+k)*(tw7+k+1)/2+l)
        enddo
        write(37,*) (tw7+k)*h,plot_c((tw7+k)*(tw7+k+1)/2+tw7),h*mtw7,(1/itemp)*(1-plot_c((tw7+k)*(tw7+k+1)/2+tw7))
      enddo
 

end program matrices

