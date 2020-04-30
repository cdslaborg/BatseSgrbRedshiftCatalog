! Amir Shahmoradi, October 17, 2012, 5:58 PM, IFS, UTEXAS.

SUBROUTINE AMH(npar,iter,ireport,idelay,X,COV,SR,func,mcmc_unit,mcmc_timing_unit,LABEL,idum_in)

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: npar,iter,ireport,idelay,mcmc_unit,mcmc_timing_unit
  INTEGER, OPTIONAL :: idum_in
  INTEGER :: i,j,ii,jj,k,idum
  DOUBLE PRECISION, EXTERNAL :: func              
  DOUBLE PRECISION, DIMENSION(npar), INTENT(INOUT) :: X    
  DOUBLE PRECISION, DIMENSION(npar,npar), INTENT(INOUT) :: COV
  DOUBLE PRECISION, DIMENSION(npar,2), INTENT(IN) :: SR    
  DOUBLE PRECISION, DIMENSION(npar) :: PROPOSAL        
  DOUBLE PRECISION, DIMENSION(0:idelay,npar) :: X0      
  DOUBLE PRECISION, DIMENSION(1,npar) :: MEAN,NEWMEAN      
  DOUBLE PRECISION, DIMENSION(npar,npar) :: SUMX0X0      
  DOUBLE PRECISION :: accr,accrdum              
  DOUBLE PRECISION :: ran2                  
  DOUBLE PRECISION :: sumaccr                  
  DOUBLE PRECISION :: sumaccr_ireport              
  DOUBLE PRECISION :: loglike_current,loglike_proposal    
  DOUBLE PRECISION :: tbegin0,tbegin,tend            
  DOUBLE PRECISION :: scalefactor                
  DOUBLE PRECISION :: dummy                  
  CHARACTER(LEN=15) :: string
  CHARACTER(LEN=8) :: labels
  CHARACTER(LEN=30), INTENT(IN), DIMENSION(npar) :: LABEL

  write(mcmc_timing_unit,'(5A25)') 'iteration','avg_progressive_accr','avg_time','time_taken_sofar','time_remained'
  labels='(   A30)'
  write(labels(2:4),'(I3)') npar+4
  write(mcmc_unit,labels) 'accepted_iter','iteration','Log(Likelihood)','avg(accr)',(TRIM(LABEL(j)),j=1,npar)
  CALL CPU_TIME(tbegin)
  CALL CPU_TIME(tbegin0)
  string='(2I30,   E30.8)'
  write(string(7:9),'(I3)') npar+2
    if (present(idum_in)) then
    idum = idum_in
    write (*,*) 'Random Number Generator (RNG) seed was set by the user to : ', idum
  else
    call random_seed()
    call random_number(dummy)
    idum = -999-NINT(98999*dummy)
    write(*,*) 'Random Number Generator Seed, idum, is: ',idum
  end if
  if (idum>-999) then
    write(*,*) 'The seed for random number generator is not appropriate...'
    write(*,*) 'idum should be a large negative number. It is: ',idum
    STOP
  end if    
  i=0
  j=0  
  jj=0
  loglike_current=func(X)
  accr=0.d0/0.d0
  write(mcmc_unit,string) i,j,loglike_current,accr,(X(k),k=1,npar)
  sumaccr=0.0d0
  sumaccr_ireport=0.0d0
  X0(i,1:npar)=X
  InitialMarkovChain: do      
    j=j+1
    jj=jj+1
    do  ! Check for the support Region consistency:
      call MVNRND(idum,npar,X0(i,1:npar),COV,PROPOSAL)
      if (any(PROPOSAL<=SR(1:npar,1)).or.any(PROPOSAL>=SR(1:npar,2))) then
        write(*,*) 'MH PROPOSAL not in the Support Region...MH will cycle.'
        write(*,*) (SR(ii,1),ii=1,npar)
        write(*,*) (PROPOSAL(ii),ii=1,npar)
        write(*,*) (SR(ii,2),ii=1,npar)
        CYCLE
      end if
      EXIT
    end do
    loglike_proposal=func(PROPOSAL)
    accrdum=dexp(loglike_proposal-loglike_current)
    if (accrdum<0.d0) then 
      write(*,*) 'Model is wrong, acceptance ratio < 0: ',accrdum
      write(*,*) 'Press Enter to Stop the program...'
      STOP
    end if
    if (loglike_proposal==-huge(loglike_proposal)) then
      accr=0.d0
      write(*,*) 'likelihood of data given model parameters is zero!'
    else
      accr=dmin1(1.,accrdum)
    end if
    sumaccr=sumaccr+accr
    sumaccr_ireport=sumaccr_ireport+accr
    if (ran2(idum)<accr) then  ! Update the stochastic variable and if the condition is met, accept the proposal as the new point.
      i=i+1
      X0(i,1:npar)=PROPOSAL
      loglike_current=loglike_proposal
      write(mcmc_unit,string) i,j,loglike_current,sumaccr/dble(j),(X0(i,k),k=1,npar)
    end if
    if (jj==ireport) then
      CALL CPU_TIME(tend)
      write(*,'(A,I5,A13,1F15.6)') 'AVG(accr) for the last ',ireport,' iterations: ',sumaccr_ireport/dble(ireport)
      write(*,'(A,I5,A,1F15.6,A)') 'Execution time for the last ',ireport,' iterations: ',tend-tbegin,' seconds.'
      write(*,'(A,1E15.6,A)') 'Elapsed time since the beginning of Markov Chain: ',tend-tbegin0,' seconds.'
      write(mcmc_timing_unit,'(1I25,4E25.8)') j,sumaccr_ireport/dble(ireport),tend-tbegin,tend-tbegin0,dble(iter-i)*(tend-tbegin0)/dble(i)
      tbegin=tend
      jj=0
      sumaccr_ireport=0.0d0
    end if
    if (mod(j,ireport)==0 .or. j==1) then
      write(*,'(I7,A,I7,A22,1F15.6)') i,' sampled. AVG(accr) for ',j,' iteration so far is: ',sumaccr/dble(j)
    end if
    if (i>=idelay) exit
    cycle
  end do InitialMarkovChain
  write(*,*) 'Average acceptance Ratio for the entire Chain was: ',sumaccr/dble(j)
  write(*,*) 'Now starting the Adaptive MH MCMC sampling...'

  scalefactor=0.8*2.38**2/dble(npar)  ! Haar scale factor
  MEAN=0.d0
  SUMX0X0=0.d0
  do ii=1,idelay
    MEAN(1,1:npar)=MEAN(1,1:npar)+X0(ii,1:npar)
    SUMX0X0=SUMX0X0+matmul(transpose(X0(ii:ii,1:npar)),X0(ii:ii,1:npar))
  end do
  MEAN(1,1:npar)=MEAN(1,1:npar)/dble(i)
  COV=(SUMX0X0-dble(i)*matmul(transpose(MEAN(1:1,1:npar)),MEAN(1:1,1:npar)))*scalefactor/dble(i-1)

  X=X0(i,1:npar)
  MarkovChain: do
    j=j+1
    jj=jj+1
    do  ! Check for the support Region consistency:
      call MVNRND(idum,npar,X,COV,PROPOSAL)
      if (any(PROPOSAL<=SR(1:npar,1)).or.any(PROPOSAL>=SR(1:npar,2))) then
        write(*,*) 'MH PROPOSAL not in the Support Region...MH will cycle.'
        write(*,*) (SR(ii,1),ii=1,npar)
        write(*,*) (PROPOSAL(ii),ii=1,npar)
        write(*,*) (SR(ii,2),ii=1,npar)
        CYCLE
      end if
      EXIT
    end do
    loglike_proposal=func(PROPOSAL)
    accrdum=dexp(loglike_proposal-loglike_current)
    if (accrdum<0.d0) then 
      write(*,*) 'Model is wrong, acceptance ratio < 0: ',accrdum
      STOP
    end if
    if (loglike_proposal==-huge(loglike_proposal)) then
      accr=0.d0
    else
      accr=dmin1(1.,accrdum)
    end if
    sumaccr=sumaccr+accr
    sumaccr_ireport=sumaccr_ireport+accr
    if (ran2(idum)<accr) then
      i=i+1
      X=PROPOSAL
      loglike_current=loglike_proposal
      write(mcmc_unit,string) i,j,loglike_current,sumaccr/dble(j),(X(k),k=1,npar)
      MEAN(1,1:npar)=(MEAN(1,1:npar)*dble(i)+X)/dble(i+1)
      X0(0,1:npar)=X-MEAN(1,1:npar)
      COV=COV*dble(i-1)/dble(i)+matmul(transpose(X0(0:0,1:npar)),X0(0:0,1:npar))*scalefactor/dble(i-1)
    end if
    if (jj==ireport) then
      CALL CPU_TIME(tend)
      write(*,'(A,I5,A13,1F15.6)') 'AVG(accr) for the last ',ireport,' iterations: ',sumaccr_ireport/dble(ireport)
      write(*,'(A,I5,A,1F15.6,A)') 'Execution time for the last ',ireport,' iterations: ',tend-tbegin,' seconds.'
      write(*,'(A,1E15.6,A)') 'Elapsed time since the beginning of Markov Chain: ',tend-tbegin0,' seconds.'
      write(mcmc_timing_unit,'(1I25,4E25.9)') j,sumaccr_ireport/dble(ireport),tend-tbegin,tend-tbegin0,dble(iter-i)*(tend-tbegin0)/dble(i)
      tbegin=tend
      jj=0
      sumaccr_ireport=0.0d0
    end if
    if (mod(j,ireport)==0 .or. j==1) then
      write(*,'(I7,A,I7,A22,1F15.6)') i,' sampled. AVG(accr) for ',j,' iteration so far is: ',sumaccr/dble(j)
    end if
    if (i>iter) exit
    cycle
  end do MarkovChain
  write(*,*) 'Average acceptance Ratio for the entire Chain was: ',sumaccr/dble(j)
  END SUBROUTINE AMH