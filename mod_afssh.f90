Module mod_afssh
!! Hammes-Schiffer, Tully, JCP 101, 4657 (1994)
implicit none
real*8, parameter :: clight=2.99792458D10,av=6.0221367D23,hbar=1.05457266D-34
real*8, parameter :: kb=1.3806503d-23
real*8, parameter :: pascal_to_atm=9.86923267d-6,kcal_to_J=6.95d-21
real*8, parameter :: amu2kg=1.66053892d-27
real*8, parameter :: au2kg=9.10938291d-31,au2J=4.35974417d-18,au2m=5.2917721092d-11,au2s=2.418884326505d-17
real*8, parameter :: q_el=2.307075d-28
complex*16,parameter :: iota = (0.d0,1.d0)
real*8 pi,wave_to_J

!! Potential
integer nquant
real*8 g_coup,epsilon
real*8 V_exothermicity,Vc,omg_B,gamma_B,temperature
real*8 beta,gamma_D,lambda_B,V_reorg,V_barrier
real*8 omg_c,omg_scaled
real*8 s01,s02,x_cr
real*8,allocatable :: mass(:),omg(:),ck(:)

!! Input/Output
real*8 cnt_frust,cnt_collapse,cnt_init,cnt_term
real*8,allocatable :: pop(:,:),pop_surf(:,:),pop_amp(:,:)
complex*16,allocatable :: rho(:,:,:)

!! Classical
integer nclass,idistribution
real*8,allocatable :: x(:),v(:),acc(:)
real*8,allocatable :: x_old(:),v_old(:),acc_old(:),x_hop(:)
real*8 tim_hop
integer iforward,flag_terminate,flag_frust,flag_reactant,flag_hop,flag_ortho
complex*16,allocatable :: delr(:,:,:),delp(:,:,:),delacc(:,:,:)
complex*16,allocatable :: delr_old(:,:,:),delp_old(:,:,:)

!! Quantized vibration
integer ncl_site
integer nb_vib,n_dvr
real*8,allocatable ::si_sho(:,:,:),sho_overlap(:,:,:,:),q_exp(:,:,:)
real*8,allocatable ::qsq_exp(:,:,:)
real*8,allocatable ::en_sho(:,:),fc_init(:)
real*8,allocatable ::Hamil_diab_0(:,:)

!! Quantum
integer state,nbasis,state_tentative
integer state_old
real*8,allocatable :: si_adiab(:,:),V_k(:),d_ij(:,:,:),vdotd(:,:),V_k_old(:)
real*8,allocatable :: Hamil_site(:,:),Hamil_diab(:,:),delH_dels(:,:,:),delH_dels_ad(:,:,:)
real*8,allocatable :: pot(:,:),force(:,:,:),force_old(:,:,:),delF(:,:,:),Hamil_diab_old(:,:)
complex*16,allocatable :: ci(:),ci_old(:),sigma(:,:)
real*8,allocatable :: si_adiab_prev(:,:)
complex*16,allocatable :: mat(:,:),mat_adiab(:,:)
real*8,allocatable :: hop_prob(:),W_overlap(:,:),hop_prob_net(:)

!! Evolution
integer n_traj,nsteps,nsteps_eq,nstep_write,iwrite,nstep_avg
real*8 dtc,total_time,curr_time,traj_num,tim_eq
real*8 energy,pot_en,KE_en,energy_cutoff,energy_old
real*8 ensq_avg,en_avg
integer nst_av
integer ihop,icollapse,iterminate,ifriction,iaverage
real*8 prob_tun

!! Parallelization
integer iparallel,iflow,iproc,iwait,ifolder

!! Misc
integer nold,cnt_rate
real*8 tim_tot,tim_ev_cl,tim_diag,tim_cl,tim_rattle,tim_pbc,tim_LJ_tot,tim_solv_solv,tim_check,tim_check2
real*8 tim_T_jk
integer,allocatable:: seed(:)
real*8,allocatable:: work(:)
complex*16,allocatable:: cwork(:)
integer,allocatable:: iwork(:),isuppz(:)

contains
!---------------------------------------------------------- 
!---------------------------------------------------------- 

subroutine setup
  implicit none
  character st_ch
  integer i,size_seed,seed2(2)
  real*8 rnd,c_0,c_e,kt

  pi=dacos(-1.d0)
  wave_to_J=2*pi*clight*hbar
  nold=0

  open(10,file="AFSSH.inp")
  read(10,*) iflow
  read(10,*) iproc
  read(10,*) iparallel
  read(10,*) iwait
  read(10,*) N_traj
  read(10,*) dtc
  read(10,*) total_time
  read(10,*) iwrite
  read(10,*) nstep_write
  read(10,*) nstep_avg
  read(10,*) idistribution
  read(10,*) flag_frust
  read(10,*) flag_ortho
  read(10,*) energy_cutoff
  read(10,*) nclass
  read(10,*) nquant
  read(10,*) n_dvr
  read(10,*) nb_vib
  read(10,*) V_exothermicity
  read(10,*) Vc
  read(10,*) omg_B
  read(10,*) lambda_B
  read(10,*) gamma_B
  read(10,*) temperature
  read(10,*) iforward
  read(10,*) icollapse
  read(10,*) ifriction
  read(10,*) seed2
  read(10,*) st_ch
  close(10)
  !----------------------------------------------------------

  if(st_ch.ne.'x') then
    write(6,*) "problem in reading input file"
    stop
  endif
  !---------------------------------------------------------- 

  !nquant=nquant*nb_vib
  nbasis=nquant
  ncl_site=nclass

  energy_cutoff=energy_cutoff*wave_to_J
  !temperature=temperature*wave_to_J/kb
  kt=kb*temperature
  Vc=Vc*wave_to_J
  V_exothermicity=V_exothermicity*wave_to_J
  omg_B=omg_B*2*pi*clight
  gamma_B=gamma_B*2*pi*clight
  !V_barrier=V_barrier*wave_to_J
  !gamma_B=gamma_B*omg_B
  lambda_B=lambda_B*wave_to_J
  !gamma_B=omg_B**2/gamma_D
  !V_reorg=4*lambda_D

!write(6,*) omg_B/(2*pi*clight),gamma_B/(2*pi*clight)
!stop
  nsteps=nint(total_time/dtc)+1
  beta=1.d0/(kb*temperature)
  !lambda_D=V_reorg/4.d0

  !-----------------------------------------------------------------  
  i=nsteps/nstep_avg+1
  allocate(pop(nquant,i),pop_surf(nquant,i),pop_amp(nquant,i))
  allocate(rho(nquant,nquant,i))
  allocate(x(nclass),v(nclass),acc(nclass))
  allocate(x_old(nclass),v_old(nclass),acc_old(nclass),x_hop(nclass))
  allocate(mass(nclass),omg(nclass),ck(nclass))
  allocate(delr(nquant,nquant,nclass),delp(nquant,nquant,nclass),delacc(nquant,nquant,nclass))
  allocate(delr_old(nquant,nquant,nclass),delp_old(nquant,nquant,nclass))
  allocate(si_adiab(nbasis,nquant),ci(nquant),V_k(nquant),V_k_old(nquant),sigma(nquant,nquant))
  allocate(Hamil_site(nbasis,nbasis),Hamil_diab(nbasis,nbasis),delH_dels(nbasis,nbasis,nclass),delH_dels_ad(nquant,nquant,nclass))
  allocate(pot(nquant,nquant),force(nquant,nquant,nclass),force_old(nquant,nquant,nclass),delf(nquant,nquant,nclass))
  allocate(Hamil_diab_old(nbasis,nbasis))
  allocate(mat(nbasis,nbasis),mat_adiab(nquant,nquant))
  allocate(d_ij(nquant,nquant,nclass),vdotd(nquant,nquant),hop_prob(nquant),W_overlap(nquant,nquant))
  allocate(hop_prob_net(nquant))
  allocate(ci_old(nquant),si_adiab_prev(nbasis,nquant))
  allocate(si_sho(n_dvr,nb_vib,2))
  allocate(sho_overlap(nb_vib,nb_vib,2,2))
  allocate(q_exp(nb_vib,nb_vib,2),qsq_exp(nb_vib,nb_vib,2))
  allocate(en_sho(nb_vib,2),fc_init(nb_vib))
  allocate(Hamil_diab_0(nbasis,nbasis))

  call random_seed(size=size_seed)
  allocate(seed(size_seed))
  do i=1,size_seed/2
    seed(i)=seed2(1)*(2*i+i*i-7)
  enddo
  do i=size_seed/2+1,size_seed
    seed(i)=seed2(2)*(i/2+34-i**3)
  enddo
  call random_seed(put=seed)
  call system_clock(count_rate=cnt_rate)
  !-----------------------------------------------------------------  

  if(iflow==2) then
    open(10,file="ifolder.inp")
    read(10,*) ifolder
    close(10)
    N_traj=N_traj/iparallel
    if(ifolder>1) then
      do i=1,(ifolder-1)*N_traj
        seed=seed+1
        call init_cond
      enddo
      call random_seed(put=seed)
    endif
  else
    ifolder=1
  endif

  tim_tot=0.d0

end subroutine setup
!---------------------------------------------------------- 

subroutine main
  implicit none
  integer i,j,k,n
  real*8 t1,t2

  call files(0)

  call cpu_time(t1)

  call setup_parameters
  call initialize_averages

!call check_acceleration

  do i=1,N_traj
    traj_num=i
    call init_cond
    call evolve(nsteps)
    call average_end
  enddo
  call write_average

  call cpu_time(t2);tim_tot=tim_tot+t2-t1
  call files(1)

end subroutine main
!---------------------------------------------------------- 

subroutine files(flag)
  implicit none
  integer,intent(in)::flag

  if(flag==0) then
    open(10,file="output")
    open(11,file="output_cl")
    open(12,file="output_qm")
    open(13,file="output_hop")
    open(14,file="output_overlap")
    open(15,file="output_dec")

    open(100,file="pop.out")
    open(101,file="cnts.out")
  else
    write(10,*)
    write(10,*)"Total time=",tim_tot
    close(10);close(11);close(12);close(13);close(14);close(15)
    close(100);close(101)
  endif

end subroutine files
!-----------------------------------------------------------------  

subroutine initialize_averages
  implicit none

  cnt_frust=0.d0
  cnt_collapse=0.d0
  cnt_init=0.d0
  cnt_term=0.d0
  pop=0.d0
  rho=0.d0
  pop_surf=0.d0
  pop_amp=0.d0

end subroutine initialize_averages
!-----------------------------------------------------------------  

subroutine init_cond
  implicit none
  integer i
  real*8 sig_x,sig_p,rnd,ak,su
  real*8 energy_0

  do i=1,nclass
    !ak=2/(hbar*omg(i))*dtanh(beta*hbar*omg(i)/2.d0) !! Wigner
    ak=beta    !! Classical
    sig_x=1.d0/dsqrt(ak*mass(i)*omg(i)**2)
    sig_p=dsqrt(mass(i)/ak)
    call gaussian_random_number(rnd)
    x(i)=rnd*sig_x-g_coup/(mass(1)*omg(1)**2)
    call gaussian_random_number(rnd)
    v(i)=(1.d0/mass(i)*(rnd*sig_p))
  enddo

  state=1

  !! quantum state initialized on diabat 1
  ci=0.d0
  ci(1)=1.d0

  iaverage=1
  iterminate=0
  flag_terminate=0

  curr_time=0.d0
  call evaluate_variables(0)
  call evaluate_variables(1)

  !! to compute the standard deviation of the energy of the trajectory
  en_avg=0.d0;ensq_avg=0.d0
  nst_av=0

end subroutine init_cond
!-----------------------------------------------------------------  

subroutine evolve(nsteps)
  implicit none
  integer,intent(in) :: nsteps
  integer i,j,nstep_sm,iflag_coll,i_do_something
  real*8 t1,t2
  integer iterm

  !call cpu_time(t1)

  call write_output(1,1)
  iterm=0
  do i=1,nsteps
    call write_output(i,0)
    call average(i)
    call save_old_state
    call evolve_classical(dtc)
    call evolve_quantum_small_dtq
    curr_time=curr_time+dtc
  enddo
  call write_output(1,1)

  !call cpu_time(t2)
  !tim_evolve=tim_evolve+t2-t1

end subroutine evolve
!-----------------------------------------------------------------  

subroutine average(i)
  implicit none
  integer,intent(in) :: i
  integer j,i1,j1,k,kp
  complex*16 ci_diab(nquant),rho_ad(nquant,nquant)
  real*8 r_avg,U(nquant,nquant),U_exc(nquant,nquant)
  integer if_reactant
  real*8 t1,t2

  !call cpu_time(t1)

  if(iwrite==1) then
    en_avg=en_avg+energy
    ensq_avg=ensq_avg+energy*energy
    nst_av=nst_av+1
  endif

  if(iaverage==1.and.(mod(i,nstep_avg)==1.or.nstep_avg==1)) then
    if(nstep_avg==1) then
      j=i
    else
      j=i/nstep_avg+1
    endif

    rho_ad=0.d0
    do i1=1,nquant
      do j1=1,nquant
        rho_ad(i1,j1)=ci(i1)*dconjg(ci(j1))
      enddo
    enddo
    rho(:,:,j)=rho(:,:,j)+rho_ad!matmul(U,matmul(rho_ad,transpose(U)))

  endif

  !call cpu_time(t2)
  !tim_coll=tim_coll+t2-t1

end subroutine average
!-----------------------------------------------------------------  

subroutine average_end
  implicit none

end subroutine average_end
!-----------------------------------------------------------------  

subroutine save_old_state
  implicit none

  x_old=x
  v_old=v
  acc_old=acc
  ci_old=ci
  state_old=state
  !ci2_old=ci2
  si_adiab_prev=si_adiab
  V_k_old=V_k
  force_old=force
  energy_old=energy
  delr_old=delr
  delp_old=delp
  Hamil_diab_old=Hamil_diab

end subroutine save_old_state
!-----------------------------------------------------------------  

subroutine revert_state
  implicit none

  x=x_old
  v=v_old
  state=state_old
  ci=ci_old
  delr=delr_old
  delp=delp_old
  force=force_old
  !ci2=ci2_old
  call evaluate_variables(0)
  call evaluate_variables(1)

end subroutine revert_state
!-----------------------------------------------------------------  

subroutine evolve_quantum_small_dtq
  implicit none
  complex*16,dimension(1:nquant):: k1,k2,k3,k4
  real*8 mat(nbasis,nbasis)

  mat=Hamil_diab_old
  k1=-iota/hbar*matmul(mat,ci)

  mat=(Hamil_diab+Hamil_diab_old)/2.d0
  k2=-iota/hbar*matmul(mat,ci+dtc/2.d0*k1)
  k3=-iota/hbar*matmul(mat,ci+dtc/2.d0*k2)

  mat=Hamil_diab
  k4=-iota/hbar*matmul(mat,ci+dtc*k3)

  ci=ci+dtc/6.d0*(k1+2*k2+2*k3+k4)

end subroutine evolve_quantum_small_dtq
!-----------------------------------------------------------------  

subroutine evolve_classical(dt)
  !! Velocity Verlet
  implicit none
  integer i
  real*8,intent(in) :: dt
  real*8 gama_dt,c0,c1,c2
  real*8 delta_r(nclass),delta_v(nclass),acc_sav(nclass)
  real*8 t1,t2

  !call cpu_time(t1)

  if(ifriction==0) then
    !! Step 1
    x=x+v*dt+0.5*acc*dt*dt
    v=v+0.5*acc*dt
    acc_old=acc
    call evaluate_variables(0)
    v=v+0.5*dt*acc
    call evaluate_variables(1)
  endif

  if(ifriction==1) then
    gama_dt=gamma_B*dt
    c0=dexp(-gama_dt)
    c1=1.d0/gama_dt*(1.d0-c0)
    c2=1.d0/gama_dt*(1.d0-c1)
     call stochastic_force(delta_r,delta_v,dt)
     x=x+c1*dt*v+c2*dt*dt*acc+delta_r
     acc_old=acc
     call evaluate_variables(0)
     v=c0*v+(c1-c2)*dt*acc_old+c2*dt*acc+delta_v
     call evaluate_variables(1)
  endif

  !call cpu_time(t2);tim_ev_cl=tim_ev_cl+t2-t1

end subroutine evolve_classical
!-----------------------------------------------------------------  

subroutine traj_terminate(iterm)
  implicit none
  integer,intent(out) :: iterm

  iterm=0

end subroutine traj_terminate
!-----------------------------------------------------------------  

subroutine write_output(n,nflag)
  !! nflag=0: Writes various variables as a function of time
  !! nflag=1: writes minimal useful information at the start and end of trajectory
  implicit none
  integer,intent(in)::nflag,n
  integer i
  real*8 t1,t2
  real*8 phase

  !call cpu_time(t1)

  if(nflag==0) then
    if(iwrite==1) then
      if(mod(n,nstep_write)==1.or.nstep_write==1) then
        write(10,'(4es17.7,i5)')curr_time*1.d15,energy/wave_to_J,sum(cdabs(ci)**2),temperature,state
        write(11,'(es15.5$)')curr_time*1.d15
        write(12,'(5f15.5)')curr_time*1.d15,cdabs(ci(1:2))**2,datan2(dimag(ci(1:2)),real(ci(1:2)))*180/pi
        !do i=1,nclass
        !  write(11,'(2es15.5$)')x(i)*1.d10,v(i)
        !enddo
        !write(11,*)
      endif
    endif
  endif

  if(nflag==1) then
    if(iwrite==0)then
      write(10,'(5es15.5)')traj_num,energy/wave_to_J,sum(cdabs(ci)**2),temperature
      write(11,*) traj_num
      write(11,'(es15.5$)')curr_time*1.d15
      do i=1,nclass
        write(11,'(2es15.5$)')x(i)*1.d10,v(i)
      enddo
      write(11,*)
      write(11,*)
    endif
    if(iwrite==1) then
      write(10,*)"traj num=",traj_num
      write(10,*)"standard deviation=",dsqrt((ensq_avg-en_avg**2/dfloat(nst_av))/dfloat(nst_av))/wave_to_J
      write(10,*)"ci**2=",sum(cdabs(ci)**2)
      write(10,*);write(10,*)
      write(11,*);write(11,*)
      write(12,*);write(12,*)
    endif
  endif

  !call cpu_time(t2)
  !tim_wr_out=tim_wr_out+t2-t1

end subroutine write_output
!-----------------------------------------------------------------  

subroutine write_average
  !! Writes the final useful output
  implicit none
  integer i,j,i1,k
  real*8 nf,pop_el(2)

  nf=dfloat(n_traj)
  cnt_frust=cnt_frust/nf
  cnt_collapse=cnt_collapse/nf

  pop=pop/nf
  rho=rho/nf
  pop_surf=pop_surf/nf
  pop_amp=pop_amp/nf

  do i=1,nsteps/nstep_avg
    write(100,'(21f15.7)')(i-1)*nstep_avg*dtc*1.d15,rho(1,:,i),rho(2,:,i)
  enddo

  write(101,*) Vc/wave_to_J,cnt_frust,cnt_collapse

end subroutine write_average
!-----------------------------------------------------------------  

subroutine evaluate_variables(flag)
  implicit none
  integer,intent(in):: flag
  integer i,j

  if(flag==0) then
    !! position dependant variables only
    call tise
  endif

  if(flag==1) then
    KE_en=0.d0
    do i=1,nclass
      KE_en=KE_en+0.5*mass(i)*v(i)*v(i)
    enddo

    energy=pot_en+KE_en
    !temperature=2*KE_en/(nclass*kb)

    !vdotd=0.d0
    !do i=1,nclass
    !  vdotd=vdotd+v(i)*d_ij(:,:,i)
    !enddo
    !call compute_vdotd
    
  endif

end subroutine evaluate_variables
!-----------------------------------------------------------------  

subroutine tise
  !! time independent schrodinger equation
  !! Output - pot_en,acc
  !! Output - V_k,d_ij
  implicit none
  integer i,j,k
  real*8 Hamil(nbasis,nbasis),ens(nbasis),vect(nbasis,nquant)
  real*8 pot_cl,acc_cl(nclass),acc_qm(nclass),dpotcl_dx(nclass)
  real*8 pot_qm,dpotqm_dx(nclass)
  real*8 si_adiab_old(nquant,nbasis)
  real*8 t1,t2

  !call cpu_time(t1)

  call compute_potential(Hamil,delH_dels)
  Hamil_diab=Hamil

  do i=1,nclass
    dpotqm_dx(i)=real(sum(conjg(ci)*matmul(delH_dels(:,:,i),ci)))
  enddo
  pot_qm=real(sum(conjg(ci)*matmul(Hamil(:,:),ci)))
!write(6,*)curr_time*1.d15,sum(conjg(ci)*matmul(Hamil(:,:),ci)),sum(conjg(ci)*matmul(delH_dels(:,:,i),ci))
  !call cpu_time(t2);tim_diag=tim_diag+(t2-t1)

  !call cpu_time(t1)

  call potential_classical(pot_cl,dpotcl_dx)
  acc_qm=-1.d0/mass*dpotqm_dx
  acc_cl=-1.d0/mass*dpotcl_dx

  pot_en=pot_cl+pot_qm
  V_k=pot_cl+pot_qm
  acc=acc_cl+acc_qm

  !call cpu_time(t2);tim_cl=tim_cl+(t2-t1)

end subroutine tise
!-----------------------------------------------------------------  

subroutine setup_parameters
  implicit none
  integer i
  real*8 si_diab(nbasis,2),Vb
  real*8 c_0,c_e
  real*8 omg_max,delw


  mass=1836.d0*au2kg
  omg_max=3*omg_B
  omg_c=2*omg_B
  delw=omg_max/real(nclass)

  g_coup=dsqrt(lambda_B*mass(1)*omg_B**2/2.d0)

  g_coup=2.61d-3*au2J/au2m
  lambda_B=2*g_coup**2/(mass(1)*omg_B**2)

  Hamil_site=0.d0
  Hamil_site(1,1)=410.d0
  Hamil_site(2,2)=530.d0

  Hamil_site(1,2)=-87.7d0; Hamil_site(2,1)=Hamil_site(1,2)

  Hamil_site=Hamil_site*wave_to_J

  omg=omg_B

end subroutine setup_parameters
!-----------------------------------------------------------------  

subroutine compute_potential(H_diab,delV_dels)
  implicit none
  real*8,intent(out) :: H_diab(nquant,nquant),delV_dels(nquant,nquant,nclass)
  integer i,j,k1
  real*8 h1,dh1(ncl_site)

  H_diab=Hamil_site
  delv_dels=0.d0

  H_diab(1,1)=H_diab(1,1)+0.5*mass(1)*omg_B**2*x(1)**2+g_coup*x(1)
  H_diab(2,2)=H_diab(2,2)+0.5*mass(1)*omg_B**2*x(1)**2-g_coup*x(1)

  delv_dels(1,1,1)=delv_dels(1,1,1)+mass(1)*omg_B**2*x(1)+g_coup
  delv_dels(2,2,1)=delv_dels(2,2,1)+mass(1)*omg_B**2*x(1)-g_coup

end subroutine compute_potential
!-----------------------------------------------------------------  

subroutine potential_classical(pot_cl,acc_cl)
  implicit none
  real*8,intent(out) :: pot_cl,acc_cl(nclass)
  integer i
  real*8 q1,q3

  pot_cl=0.d0
  acc_cl=0.d0

end subroutine potential_classical
!-----------------------------------------------------------------  

subroutine check_acceleration
  !! A test subroutine that compares analytical accelerations with numerical
  !accelerations
  implicit none
  integer i,nflag
  real*8 delx,en_old,acc_sav(nclass)
  real*8 q0,rnd

  delx=1.d-17
  state=1

  do i=1,nclass
    call random_number(rnd)
    x(i)=(rnd*2-1.d0)*1.d-10
  enddo

  call init_cond

  call evaluate_variables(0)
  en_old=pot_en;acc_sav=acc

  write(6,*) "delx=",delx
  write(6,*)

  do i=1,nclass
      x(i)=x(i)+delx
      call evaluate_variables(0)
      acc(i)=-(pot_en-en_old)/delx/mass(i)
      write(6,*)"Analytical acceleration =",acc_sav(i)
      write(6,*)"Numerical acceleration  =",acc(i)
      write(6,*)"Error =",(acc(i)-acc_sav(i))/acc(i)*100.d0
      write(6,*)
      x(i)=x(i)-delx
  enddo

  stop

end subroutine check_acceleration
!---------------------------------------------------------- 

subroutine stochastic_force(delr,delv,dt)
  !! stoachastic forces for langevin equation
  !! Not used for the Holstein model results 
  implicit none
  real*8,intent(in)::dt
  real*8,intent(out) :: delr(nclass),delv(nclass)!f(nclass)
  integer i
  real*8 rnd1,rnd2,sig_r,sig_v,sig_rv,gdt

  gdt=gamma_B*dt

  do i=1,nclass

    sig_r=dt*dsqrt(kb*temperature/mass(i) *1.d0/gdt*(2-1.d0/gdt*(3-4*dexp(-gdt)+dexp(-2*gdt))))
    sig_v=dsqrt(kb*temperature/mass(i)*(1-dexp(-2*gdt)))
    sig_rv=(dt*kb*temperature/mass(i)* 1.d0/gdt *(1-dexp(-gdt))**2)/(sig_r*sig_v)  !! correlation coeffecient

    call gaussian_random_number(rnd1)
    call gaussian_random_number(rnd2)
    delr(i)=sig_r*rnd1
    delv(i)=sig_v*(sig_rv*rnd1+dsqrt(1-sig_rv**2)*rnd2)
  enddo

!  delr=delr-sum(delr)/dfloat(nclass)
!  delv=delv-sum(delv)/dfloat(nclass)

end subroutine stochastic_force
!-----------------------------------------------------------------  

subroutine gaussian_random_number(rnd)
  !! generates gaussian distribution with center 0, sigma 1
  !! q0+sig*rnd gives center=q0, sigma=sig
  implicit none
  real*8,intent(out)::rnd
  real*8 rnd1,rnd2,pi

  pi=dacos(-1.d0)

  call random_number(rnd1)
  call random_number(rnd2)
  rnd = dsqrt(-2*log(rnd1))*dcos(2*pi*rnd2)

end subroutine gaussian_random_number
!---------------------------------------------------------- 

End Module mod_afssh
