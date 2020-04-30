! This file contains all the modules containing common variables and parameters used by different functions and subroutines of the program.
! Amir Shahmoradi, Monday Dec 9, 2013, 9:04, IFS, The University of Texas at Austin.

MODULE MODELparameters

	IMPLICIT NONE
	INTEGER, PARAMETER :: nvar=4							! The number of GRB observational (i.e. model) variables.
	DOUBLE PRECISION, DIMENSION(nvar,nvar) :: MVNCOV		! The 3*3 covariance matrix of the MVN distribution
	DOUBLE PRECISION, DIMENSION(nvar,nvar) :: INVMVNCOV		! Inverse of the 3*3 covariance matrix of the MVN distribution
	DOUBLE PRECISION, DIMENSION(nvar) :: X					! The vector of rest frame GRB variables.
	! For the parameters below, the array elements correspond to those of LGRBs (1) and SGRBs (2) respectively.
	DOUBLE PRECISION :: conlisomeanD	        ! mean of the conditional dist. of Log(Liso) given Log(duration_z).
	DOUBLE PRECISION :: conepkzmeanD	        ! mean of the conditional dist. of Log(Epks) given Log(Epkz).
	DOUBLE PRECISION :: conlisosigmaD	        ! Stdev of the conditional dist. of Log(Liso) given Log(duration_z).
	DOUBLE PRECISION :: conepkzsigmaD	        ! Stdev of the conditional dist. of Log(Epkz) given Log(duration_z).
	DOUBLE PRECISION :: conepkzmeanLD	        ! mean of the conditional dist. of Log(Liso) given Log(Liso)_given Log(duration_z).
	DOUBLE PRECISION :: conepkzsigmaLD	        ! Stdev of the conditional dist. of Log(Epkz) given Log(Liso)_given Log(duration_z).
	DOUBLE PRECISION :: aELD,bELD		        ! terms used in the conditional dist. of Log(Epkz) given Log(Liso) given Log(durationz),
	DOUBLE PRECISION :: aLD,bLD			        ! terms used in the conditional dist. of Log(Liso) given Log(durationz),
	DOUBLE PRECISION :: aED,bED			        ! terms used in the conditional dist. of Log(Epkz) given Log(durationz),
	DOUBLE PRECISION :: dummy2			        ! a term that is in common with Loglike function and the model integration function.
	DOUBLE PRECISION :: sqrt2lisosigma	        ! = sqrt2*lisosigma, used in model integration, defined for better code performance.
	DOUBLE PRECISION :: sqrt2t90zsigma	        ! = sqrt2*t90zsigma, used in model integration, defined for better code performance.
	DOUBLE PRECISION :: sqrt2Pit90zsigma		! = sqrt2Pi*t90zsigma, used in model integration, defined for better code performance.
	DOUBLE PRECISION :: sqrt2conlisosigmaD		! = sqrt2*conlisosigmaD, used in model integration, defined for better code performance.
	DOUBLE PRECISION :: sqrt2conepkzsigmaLD		! = sqrt2*conepkzsigmaLD, used in model integration, defined for better code performance.
	DOUBLE PRECISION :: sqrt2PiconlisosigmaD	! = sqrt2Pi*conlisosigmaD, used in model integration, defined for better code performance.
	DOUBLE PRECISION :: sqrt2PiconepkzsigmaLD	! = sqrt2Pi*conepkzsigmaLD, used in model integration, defined for better code performance.
												! It is in fact, the normalization factor of the conditional distribution of Liso on durationz.
	DOUBLE PRECISION :: loglisomean		        ! mean of the marginal dist. of Log(Liso).
	DOUBLE PRECISION :: logepkzmean		        ! mean of the marginal dist. of Log(Epkz).
	DOUBLE PRECISION :: logeisomean		        ! mean of the marginal dist. of Log(Eiso).
	DOUBLE PRECISION :: logt90zmean		        ! mean of the marginal dist. of Log(T90z).
	DOUBLE PRECISION, DIMENSION(9:14) :: RHO	! The array of the correlation coefficients obtained from the inverse Fisher transformation.
	DOUBLE PRECISION :: modelint				! integral of the model over the redshift range [zmin,zmax].
	DOUBLE PRECISION :: rhoLE_given_Durz
	! The parameter rhoLE_given_Durz is the partial correlation of Liso & Epkz conditional on rest-frame duration:
	! It is used in the integration of the model where the joint distribution of Lis-Epkz given rest-frame duration
	! has to be integrated. The relation can be found at http://en.wikipedia.org/wiki/Partial_correlation
	! The integration on duration is done according to qromb.f90 algorithm between the above two limits, and according to
	! qromo.f90 (with midexp.f90 & midexpMirror.f90) below and beyond this range.
	! DOUBLE PRECISION :: LGRBprop	! The proportion (fraction, or normalization) of the Gaussian component for LGRBs MVN dist.
									! the SGRBprop is by default SGRBprop=(1-LGRBprop) in the case of two component Gaussian Mixture.
END MODULE MODELparameters

MODULE Zparameters
	IMPLICIT NONE
	!INTEGER, PARAMETER :: nz=3	! 3^(nz-1) gives the total number of steps in redshift intergations.
	!DOUBLE PRECISION, PARAMETER :: zmin=0.1d0,zmax=6.5d0		! Refer to Butler et al (2010) for desciption of parameters.
	DOUBLE PRECISION, PARAMETER :: zmin=0.5d0,zmax=5.0d0		! Refer to Butler et al (2010) for desciption of parameters.
	!DOUBLE PRECISION, PARAMETER :: z0=0.993d0,z1=3.8d0			! Refer to Butler et al (2010) for desciption of parameters.
	!DOUBLE PRECISION, PARAMETER :: g0=3.3d0,g1=0.0549d0,g2=-4.46d0
	DOUBLE PRECISION, PARAMETER :: z0=0.97d0,z1=4.5d0			! Refer to Butler et al (2010) for desciption of parameters.
	DOUBLE PRECISION, PARAMETER :: g0=3.4d0,g1=-0.3d0,g2=-7.8d0
	DOUBLE PRECISION :: gamma1									! (1.+z0)**(g0-g1)
	DOUBLE PRECISION :: gamma2									! gamma1*(1.+z1)**(g1-g2)
																! Refer to the 2nd Blue Notebook for a description of the
																! comoving rate equation and the factors gamma1 and gamma2.
																! the normalization of the model, also when marginalizing the model over z.
	DOUBLE PRECISION :: zplus1,logzplus1,lumdisterm,dvdzOzplus1,logzplus1_dur
END MODULE Zparameters

MODULE COSMOparameters
	IMPLICIT NONE
	DOUBLE PRECISION, PARAMETER :: omega_DE=0.7d0			! Dark Energy density.
	DOUBLE PRECISION, PARAMETER :: omega_DM=0.3d0			! Dark Matter density.
END MODULE COSMOparameters

MODULE constants
	IMPLICIT NONE
	DOUBLE PRECISION, PARAMETER :: c=3.0d5,H=7.1d1			! C is the speed of light (Km/s), H is the Hubble constant in units of km/s/MPc.
	DOUBLE PRECISION, PARAMETER :: Mpc2cm=3.09d24			! 1 Mega Parsec = Mpc2cm centimeters.
	DOUBLE PRECISION, PARAMETER :: pi=dacos(-1.0d0)			! The irrational number pi.
	DOUBLE PRECISION, PARAMETER :: log10Mpc2cmSq4pi=dlog10(4.*pi)+2*dlog10(Mpc2cm)
															! Log10(Mpc2cm centimeters.
	DOUBLE PRECISION, PARAMETER :: sqrt2=dsqrt(2.0d0)		! Square root of 2.
	DOUBLE PRECISION, PARAMETER :: CoverH=C/H				! the speed of light in units of km/s divided by the Hubble constant.
	DOUBLE PRECISION, PARAMETER :: log10e=0.434294481903259	! Log10 of Napier constant.
	DOUBLE PRECISION, PARAMETER :: sqrtPiO2=1.2533141373155	! Square root of Pi/2. This is needed in the integration of
															! Epk and ELR density functions.
	DOUBLE PRECISION, PARAMETER :: sqrt2Pi=2.50662827463100	! Square root of 2*Pi used in integration of the model.
	!DOUBLE PRECISION, PARAMETER :: sqrt2=1.4142135623730951	! Square root of 2. This is needed in the integration of Epk and ELR density functions.
	!DOUBLE PRECISION, PARAMETER :: erf3=0.9999779095030014	! I assume this to be 1.
END MODULE constants

MODULE detection
	IMPLICIT NONE
	!DOUBLE PRECISION, PARAMETER :: signif=3.8d0 ! x0 & x1 assumed to be SIGNIF(icance)*stdevthresh far from meanthresh.
	DOUBLE PRECISION, PARAMETER :: signif=3.0d0 ! x0 & x1 assumed to be SIGNIF(icance)*stdevthresh far from meanthresh.
	DOUBLE PRECISION :: meanthresh	! MEAN value of BATSE efficiency for any GRB duration.
	DOUBLE PRECISION :: stdevthresh	! Standard deviation of BATSE efficiency. 
	DOUBLE PRECISION, PARAMETER :: Serfc=0.282313526464596		! The scale of the change in BATSE efficiency for different GRB durations.
	DOUBLE PRECISION, PARAMETER :: meandur=-0.483553339256463	! MEAN duration in the Error function used to model the connection between the peak fluxes in 64 and 1024 ms.
	DOUBLE PRECISION, PARAMETER :: scaledur=1.05146989846948	! scale of the duration in the Error function used to model the connection between the peak fluxes in 64 and 1024 ms.
	DOUBLE PRECISION :: meandurz	! meandur at redshift z
	DOUBLE PRECISION :: logdurzmax,logdurzmin	! rest-frame durations above and below which meanthresh_at_dur becomes practically independent of duration.
	! Below are observer-frame durations above and below which meanthresh_at_dur becomes practically independent of duration.
	DOUBLE PRECISION, PARAMETER :: logdurmax=2.5d0,logdurmin=-3.2d0
	!DOUBLE PRECISION, PARAMETER :: logdurmax=3.0d0,logdurmin=-4.0d0
	DOUBLE PRECISION :: eff_min_pph	! The EFFective Log(P1024_ph) below which the trigger efficiency becomes practically zero.
	DOUBLE PRECISION :: eff_max_pph	! The EFFective Log(P1024_ph) above which the trigger efficiency becomes practically 100%.
	DOUBLE PRECISION :: eff_min_lpb	! The EFFective Log(Pbol) below which the trigger efficiency becomes practically zero.
	DOUBLE PRECISION :: glb_max_lpb	! The GLoBal Log(Pbol) above which the trigger efficiency becomes practically 100%.
	! not needed: DOUBLE PRECISION :: x1			! The Log(P1024_ph) at which the trigger efficiency becomes practically maximum at the given duration.
									! This parameter depends on the duration (T90 or FPR64) of the burst.
	DOUBLE PRECISION :: sqrt2stdevthresh	! = sqrt2*stdevthresh, used in model integration, defined for better efficiency of the code.
	DOUBLE PRECISION, PARAMETER :: eminbol=20.d0,emaxbol=2000.d0	! Entire (bolometric) energy window of BATSE.
	DOUBLE PRECISION, PARAMETER :: emindet=50.d0,emaxdet=300.d0		! trigger energy window of BATSE.
	!DOUBLE PRECISION, PARAMETER :: logp1024maxlogepk=1.9744802198112485d0
	DOUBLE PRECISION, PARAMETER :: logp1024min=4.92d0,logp1024max=6.318167895318538d0
	!DOUBLE PRECISION, PARAMETER :: logepkmin=-2.915056638230699d0,logepkmax=5.4093868613659435d0
	DOUBLE PRECISION, PARAMETER :: logepkmin=-0.5d0,logepkmax=5.d0
	DOUBLE PRECISION, PARAMETER :: logphminmaxdiff=logp1024max-logp1024min
	DOUBLE PRECISION, PARAMETER :: lpb_correction=logphminmaxdiff-logp1024min+2.d0*Serfc
	DOUBLE PRECISION :: logepkzmin,logepkzmax
	DOUBLE PRECISION :: log10pbol	! used in the integration of model.
	DOUBLE PRECISION :: logp1024ph_eff	! effective peak flux Log(P1024ph) at the given duration.
	DOUBLE PRECISION :: min_lpb_at_dur	! MINimum Log(PBol) AT a Given observed DURation below which trigger efficiency becomes practially zero.
	DOUBLE PRECISION :: pph_correction	! The erfc term that corrects Log(P1024ph) to the effective Log(P1024ph).
END MODULE detection

MODULE GRBworld
	IMPLICIT NONE
	INTEGER, PARAMETER :: npar=16							! number of GRB world model parameters.
	DOUBLE PRECISION, DIMENSION(npar) :: STCV				! The SToChastic Variable of MCMC that contains GRB world model parameters.
	DOUBLE PRECISION, DIMENSION(npar) :: STDEV			! The STandard DEViation vector of the stochastic variable vector.
	DOUBLE PRECISION, DIMENSION(npar,npar) :: COV		! Covariance Matrix of the GRBworld model parameters.
	DOUBLE PRECISION, DIMENSION(npar,npar) :: COR		! Correlation Matrix of the GRBworld model parameters.
	DOUBLE PRECISION, DIMENSION(npar,2) :: SR			! Support Region of the target distribution.
END MODULE GRBworld

MODULE OBSGRBDATA
	IMPLICIT NONE
	INTEGER, PARAMETER :: ndata=565
	INTEGER :: idata	! idata is the index of the array of data which is required as global when doing likelihood integrations.
	TYPE :: GRBDATA
		INTEGER :: trigger
		DOUBLE PRECISION :: logp1024ph,logpbol,logepk,logsbol,logt90,logprob,efflogp1024ph	! LGRBprob,SGRBprob
		! logp1024ph here in the input is infact the effective 1-sec peak flux of BATSE GRBs based on their Log(T90).
		! efflogp1024ph in the input is the Effective Triggering peak flux required for detection of the burst.
	END TYPE GRBDATA
	!TYPE(GRBDATA), DIMENSION(ndata) :: GRB
END MODULE OBSGRBDATA

MODULE Bandmodel
	IMPLICIT NONE
	DOUBLE PRECISION, PARAMETER :: alpha=-1.1,beta=-2.3		! These are the photon indices of the Band model of GRBs.
	!DOUBLE PRECISION, PARAMETER :: bolemin=20.d0,bolemax=2.d3	! Bolometric BATSE energy range
	!DOUBLE PRECISION, PARAMETER :: batsemin=50.d0,batsemax=3.d2	! Detection BATSE energy range
END MODULE Bandmodel