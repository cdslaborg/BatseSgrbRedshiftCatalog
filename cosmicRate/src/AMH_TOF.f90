!   This code is a modified version of the code AMH_TO.f90, which takes in also the name of the output MCMC and MCMC_timing files.
! Amir Shahmoradi, Sunday 1:41 AM, Dec 11, 2013, IFS/ICMB, The University of Texas at Austin	 
! 
! This code is an optimized version of AMH_T.f90 subroutine. I cut the sampling time to half of AMH_T.f90 by eliminating the redundant twice
! calculation of the likelihood funtion in the subroutine.
! Amir Shahmoradi, December 10, 2012, 2:02 PM, IFS, UTEXAS.
! NOTE: In this version I'm also writing the time it takes for each new proposal to an output file. This is OK for heavy likelihood calculations. However, if the
! Markov sampling is very fast (that is the likelihood calculation is fast) then this might slow down the entire code.
! Amir Shahmoradi, November 5, 2012, 2:09 PM, IFS, UTEXAS.
! NOTE: This code assumes the maximum number of model parameters npar<100, otherwise there MIGHT be format writing issues in the output files.
! UPDATE: This version has automatic RNG seed generator, so that it automatically avoids generating the same chain of random numbers on every run.
! 	Also, the labels of the parameters are taken in the input and written in the output. It is assumed that no label is longer than 30 characters.
! Amir Shahmoradi, November 5, 2012, 1:22 PM, IFS, UTEXAS.
! Given an external function f(X), where X is the array of length npar of the variables,
! this subroutine tries to sample f(X) using an Adaptive Metropolis-Hastings Markov Chain Monte Carlo (MCMC) algorithm.
! An example AMH algorithm has been discussed by Haario (2001).
! There are 4 input integer parameters:
! 	1. npar: the number of stochastic variables (dimension of the space to be explored).
! 	2. iter: the total number of iterations for the Markov Chain.
! 	3. ireport: every ireoprt loops of MCMC chain, the average of the acceptance ratio AVG(accr) will reported on the screen. Also reported
! 		are the average time taken to obtain ireoprt sample, and the time since the beginning of the chain.
! 	4. iupdate: every iupdate, update the covariance matrix.
! 	5. idelay: integer number, 0 or 1. if idelay=0, iupdate remains constant throughout sampling. If idelay=1, double iupdate for
! 		 the next covariance matrix update.
! The algorithm uses an npar dimensional Guassian as the proposal distribution.
! The algorithm gets three input arrays as the initializers:
! 	1. X: the npar-long 1D array of the initial starting point for the chain.
! 	2. COV: the initial npar*npar 2D covariance matrix of the proposal distribution.
! The user must provide an external function "func" which is the function to be sampled.
! Note that func could be in fact the logarithm of the original function to be sampled.
! The resulting sample from the Markov chain will be stored in two files:
! 	1. "MHMCMC_SAMPLE.txt" in real time. This file only has the sample after burn iterations.
! 	1. "MHMCMC_RAW_SAMPLE.txt" in real time. This file includes all of the accepted proposals. 
! USES subroutines & functions: samcovmat.f90, mvndevran2.f90, ran2.f90
! Amir Shahmoradi, April 13, 2012, 7:39 PM, IFS, UTEXAS.
! VERSION 2. This version of the code, also takes in as the input, the boundaries (the support region) of the distribution to be sampled.
! Amir Shahmoradi, April 19, 2012, 2:10 AM, IFS, UTEXAS.
! Amir Shahmoradi, July 17, 2012, 1:08 AM, IFS, UTEXAS.
! Version.3, Amir Shahmoradi, October 17, 2012, 5:58 PM, IFS, UTEXAS.

SUBROUTINE AMH(npar,iter,ireport,idelay,X,COV,SR,func,mcmc_unit,mcmc_timing_unit,LABEL,idum_in)

	IMPLICIT NONE
	INTEGER, INTENT(IN) :: npar,iter,ireport,idelay,mcmc_unit,mcmc_timing_unit	!iupdate
	! mcmc_unit,mcmc_timing_unit are the output file unit numbers.
	INTEGER, OPTIONAL :: idum_in						! If idum_in is present, fix the RNG seed to idum_in.
	INTEGER :: i,j,ii,jj,k,idum
	DOUBLE PRECISION, EXTERNAL :: func								! function to be sampled
	DOUBLE PRECISION, DIMENSION(npar), INTENT(INOUT) :: X			! X is the npar-long array of the stochastic variable.
	DOUBLE PRECISION, DIMENSION(npar,npar), INTENT(INOUT) :: COV	! covariance matrix of the proposal distribution.
	DOUBLE PRECISION, DIMENSION(npar,2), INTENT(IN) :: SR			! Suport Region of the target distribution.
	DOUBLE PRECISION, DIMENSION(npar) :: PROPOSAL					! npar-long random deviate array drawn from the proposal distribution.
	DOUBLE PRECISION, DIMENSION(0:idelay,npar) :: X0				! npar-long random deviate array drawn initialy from the proposal distribution.
	DOUBLE PRECISION, DIMENSION(1,npar) :: MEAN,NEWMEAN				! npar-long MEAN vector of X0 vectors.
	DOUBLE PRECISION, DIMENSION(npar,npar) :: SUMX0X0				! npar*npar matrix of X0 cross its transpose.
	DOUBLE PRECISION :: accr,accrdum								! accr stands for the ACCeptance Ratio.
	DOUBLE PRECISION :: ran2										! the random-number-generator function
	DOUBLE PRECISION :: sumaccr										! This to figure out what the average accpetance ratio for the entire chain is.
	DOUBLE PRECISION :: sumaccr_ireport								! This to figure out what the average accpetance ratio for ireport length of chain is.
	DOUBLE PRECISION :: loglike_current,loglike_proposal			! This is the loglikelihood (to be written in an output file) when proposal is accpted.
	DOUBLE PRECISION :: tbegin0,tbegin,tend							! for measuring the execution time.
	DOUBLE PRECISION :: scalefactor									! Gelman et al. (1996)
	DOUBLE PRECISION :: dummy										! used for random seed generation.
	CHARACTER(LEN=15) :: string
	CHARACTER(LEN=8) :: labels
	CHARACTER(LEN=30), INTENT(IN), DIMENSION(npar) :: LABEL						! Parameter labels (names) to be written to the output file header.
	!CHARACTER(LEN=300) :: file_mcmc_out,file_mcmc_timing_out

	!open(unit=mcmc_unit,file=file_mcmc_out,status='replace')			! ,status='new'
	!open(unit=mcmc_timing_unit,file=file_mcmc_timing_out,status='replace')	! ,status='new'
	write(mcmc_timing_unit,'(5A25)') 'iteration','avg_progressive_accr','avg_time','time_taken_sofar','time_remained'
	labels='(   A30)'
	write(labels(2:4),'(I3)') npar+4
	!write(*,*) labels
	!read(*,*)
	write(mcmc_unit,labels) 'accepted_iter','iteration','Log(Likelihood)','avg(accr)',(TRIM(LABEL(j)),j=1,npar)
	CALL CPU_TIME(tbegin)
	CALL CPU_TIME(tbegin0)
	string='(2I30,   E30.8)'
	write(string(7:9),'(I3)') npar+2
    if (present(idum_in)) then
		!write (*,*) 'Amir',present(idum_in),idum_in
		idum = idum_in
		write (*,*) 'Random Number Generator (RNG) seed was set by the user to : ', idum
	else
		call random_seed()
		call random_number(dummy)
		idum = -999-NINT(98999*dummy)
		write(*,*) 'Random Number Generator Seed, idum, is: ',idum
	end if
    !write(*,*) dummy,idum
    !read(*,*)
	if (idum>-999) then
		write(*,*) 'The seed for random number generator is not appropriate...'
		write(*,*) 'idum should be a large negative number. It is: ',idum
		STOP
	end if		
	i=0	! Marcov Chain acceptance counter
	j=0	! Marcov Chain counter
	jj=0	! Marcov Chain ireport counter
	loglike_current=func(X)
	accr=0.d0/0.d0	! result is NaN
	write(mcmc_unit,string) i,j,loglike_current,accr,(X(k),k=1,npar)
	sumaccr=0.0d0
	sumaccr_ireport=0.0d0
	X0(i,1:npar)=X
	InitialMarkovChain: do			
		j=j+1
		jj=jj+1
		do	! Check for the support Region consistency:
			call MVNRND(idum,npar,X0(i,1:npar),COV,PROPOSAL)	! Draw a proposal as the next poential step in the chain.
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
			!write(*,*) 'loglike==-huge(loglike)'
			!read(*,*)
			accr=0.d0
			write(*,*) 'likelihood of data given model parameters is zero!'
		else
			accr=dmin1(1.,accrdum)
		end if
		sumaccr=sumaccr+accr
		sumaccr_ireport=sumaccr_ireport+accr
		if (ran2(idum)<accr) then	! Update the stochastic variable and if the condition is met, accept the proposal as the new point.
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
!	----------------------------------------------------------------------------------------------------------------------------
	! First calculate the mean vector of X0 and also SUMX0X0 matrix, then the covariance matrix COV:
		scalefactor=0.8*2.38**2/dble(npar)	! Haar (2001) scale factor.=
		MEAN=0.d0
		SUMX0X0=0.d0
		do ii=1,idelay
			MEAN(1,1:npar)=MEAN(1,1:npar)+X0(ii,1:npar)
			SUMX0X0=SUMX0X0+matmul(transpose(X0(ii:ii,1:npar)),X0(ii:ii,1:npar))
		end do
		MEAN(1,1:npar)=MEAN(1,1:npar)/dble(i)	! i==idelay
		COV=(SUMX0X0-dble(i)*matmul(transpose(MEAN(1:1,1:npar)),MEAN(1:1,1:npar)))*scalefactor/dble(i-1)	! i==idelay
!	----------------------------------------------------------------------------------------------------------------------------
	X=X0(i,1:npar)
	MarkovChain: do
		j=j+1
		jj=jj+1
		do	! Check for the support Region consistency:
			call MVNRND(idum,npar,X,COV,PROPOSAL)	! Draw a proposal as the next poential step in the chain.
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
			!write(*,*) 'loglike==-huge(loglike)'
			!read(*,*)
			accr=0.d0
		else
			accr=dmin1(1.,accrdum)
		end if
		sumaccr=sumaccr+accr
		sumaccr_ireport=sumaccr_ireport+accr
		if (ran2(idum)<accr) then	! Update the stochastic variable.
			i=i+1
			X=PROPOSAL
			loglike_current=loglike_proposal
			write(mcmc_unit,string) i,j,loglike_current,sumaccr/dble(j),(X(k),k=1,npar)
			! Update MEAN of X:
				!NEWMEAN(1,1:npar)=(MEAN(1,1:npar)*dble(i)+X)/dble(i+1)
				MEAN(1,1:npar)=(MEAN(1,1:npar)*dble(i)+X)/dble(i+1)
			! Update the covariance matrix:
				X0(0,1:npar)=X-MEAN(1,1:npar)
				COV=COV*dble(i-1)/dble(i)+matmul(transpose(X0(0:0,1:npar)),X0(0:0,1:npar))*scalefactor/dble(i-1)
			!MEAN=NEWMEAN
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
	!close(mcmc_unit)
	write(*,*) 'Average acceptance Ratio for the entire Chain was: ',sumaccr/dble(j)
	END SUBROUTINE AMH