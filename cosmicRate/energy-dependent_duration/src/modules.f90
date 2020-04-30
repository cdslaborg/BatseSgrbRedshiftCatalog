! Amir Shahmoradi, Monday Dec 9, 2013, 9:04, IFS, The University of Texas at Austin.

MODULE MODELparameters
  IMPLICIT NONE
  INTEGER, PARAMETER :: nvar=4              
  DOUBLE PRECISION, DIMENSION(nvar,nvar) :: MVNCOV    
  DOUBLE PRECISION, DIMENSION(nvar,nvar) :: INVMVNCOV
  DOUBLE PRECISION, DIMENSION(nvar) :: X          
  DOUBLE PRECISION :: conlisomeanD          
  DOUBLE PRECISION :: conepkzmeanD          
  DOUBLE PRECISION :: conlisosigmaD          
  DOUBLE PRECISION :: conepkzsigmaD          
  DOUBLE PRECISION :: conepkzmeanLD          
  DOUBLE PRECISION :: conepkzsigmaLD          
  DOUBLE PRECISION :: aELD,bELD            
  DOUBLE PRECISION :: aLD,bLD              
  DOUBLE PRECISION :: aED,bED              
  DOUBLE PRECISION :: dummy2              
  DOUBLE PRECISION :: sqrt2lisosigma          
  DOUBLE PRECISION :: sqrt2t90zsigma          
  DOUBLE PRECISION :: sqrt2Pit90zsigma    
  DOUBLE PRECISION :: sqrt2conlisosigmaD    
  DOUBLE PRECISION :: sqrt2conepkzsigmaLD    
  DOUBLE PRECISION :: sqrt2PiconlisosigmaD  
  DOUBLE PRECISION :: sqrt2PiconepkzsigmaLD  
                        
  DOUBLE PRECISION :: loglisomean            
  DOUBLE PRECISION :: logepkzmean            
  DOUBLE PRECISION :: logeisomean            
  DOUBLE PRECISION :: logt90zmean            
  DOUBLE PRECISION, DIMENSION(9:14) :: RHO  
  DOUBLE PRECISION :: modelint        
  DOUBLE PRECISION :: rhoLE_given_Durz
END MODULE MODELparameters

MODULE Zparameters
  IMPLICIT NONE
  DOUBLE PRECISION, PARAMETER :: zmin=0.5d0,zmax=5.0d0
  DOUBLE PRECISION, PARAMETER :: z0=0.97d0,z1=4.5d0  
  DOUBLE PRECISION, PARAMETER :: g0=3.4d0,g1=-0.3d0,g2=-7.8d0
  DOUBLE PRECISION :: gamma1
  DOUBLE PRECISION :: gamma2
  DOUBLE PRECISION :: zplus1,logzplus1,lumdisterm,dvdzOzplus1,logzplus1_dur
END MODULE Zparameters

MODULE COSMOparameters
  IMPLICIT NONE
  DOUBLE PRECISION, PARAMETER :: omega_DE=0.7d0
  DOUBLE PRECISION, PARAMETER :: omega_DM=0.3d0
END MODULE COSMOparameters

MODULE constants
  IMPLICIT NONE
  DOUBLE PRECISION, PARAMETER :: c=3.0d5,H=7.1d1
  DOUBLE PRECISION, PARAMETER :: Mpc2cm=3.09d24 
  DOUBLE PRECISION, PARAMETER :: pi=dacos(-1.0d0)
  DOUBLE PRECISION, PARAMETER :: log10Mpc2cmSq4pi=dlog10(4.*pi)+2*dlog10(Mpc2cm)
  DOUBLE PRECISION, PARAMETER :: sqrt2=dsqrt(2.0d0)
  DOUBLE PRECISION, PARAMETER :: CoverH=C/H        
  DOUBLE PRECISION, PARAMETER :: log10e=0.434294481903259  
  DOUBLE PRECISION, PARAMETER :: sqrtPiO2=1.2533141373155  
  DOUBLE PRECISION, PARAMETER :: sqrt2Pi=2.50662827463100  
END MODULE constants

MODULE detection
  IMPLICIT NONE
  DOUBLE PRECISION, PARAMETER :: signif=3.0d0
  DOUBLE PRECISION :: meanthresh  
  DOUBLE PRECISION :: stdevthresh 
  DOUBLE PRECISION, PARAMETER :: Serfc=0.282313526464596    
  DOUBLE PRECISION, PARAMETER :: meandur=-0.483553339256463 
  DOUBLE PRECISION, PARAMETER :: scaledur=1.05146989846948  
  DOUBLE PRECISION :: meandurz
  DOUBLE PRECISION :: logdurzmax,logdurzmin  
  DOUBLE PRECISION, PARAMETER :: logdurmax=2.5d0,logdurmin=-3.2d0
  DOUBLE PRECISION :: eff_min_pph
  DOUBLE PRECISION :: eff_max_pph
  DOUBLE PRECISION :: eff_min_lpb
  DOUBLE PRECISION :: glb_max_lpb
  DOUBLE PRECISION :: sqrt2stdevthresh
  DOUBLE PRECISION, PARAMETER :: eminbol=20.d0,emaxbol=2000.d0
  DOUBLE PRECISION, PARAMETER :: emindet=50.d0,emaxdet=300.d0
  DOUBLE PRECISION, PARAMETER :: logp1024min=4.92d0,logp1024max=6.318167895318538d0
  DOUBLE PRECISION, PARAMETER :: logepkmin=-0.5d0,logepkmax=5.d0
  DOUBLE PRECISION, PARAMETER :: logphminmaxdiff=logp1024max-logp1024min
  DOUBLE PRECISION, PARAMETER :: lpb_correction=logphminmaxdiff-logp1024min+2.d0*Serfc
  DOUBLE PRECISION :: logepkzmin,logepkzmax
  DOUBLE PRECISION :: log10pbol     
  DOUBLE PRECISION :: logp1024ph_eff
  DOUBLE PRECISION :: min_lpb_at_dur
  DOUBLE PRECISION :: pph_correction
END MODULE detection

MODULE GRBworld
  IMPLICIT NONE
  INTEGER, PARAMETER :: npar=16                
  DOUBLE PRECISION, DIMENSION(npar) :: STCV    
  DOUBLE PRECISION, DIMENSION(npar) :: STDEV   
  DOUBLE PRECISION, DIMENSION(npar,npar) :: COV
  DOUBLE PRECISION, DIMENSION(npar,npar) :: COR
  DOUBLE PRECISION, DIMENSION(npar,2) :: SR    
END MODULE GRBworld

MODULE OBSGRBDATA
  IMPLICIT NONE
  INTEGER, PARAMETER :: ndata=565
  INTEGER :: idata
  TYPE :: GRBDATA
    INTEGER :: trigger
    DOUBLE PRECISION :: logp1024ph,logpbol,logepk,logsbol,logt90,logprob,efflogp1024ph
  END TYPE GRBDATA
  TYPE(GRBDATA), DIMENSION(ndata) :: GRB
END MODULE OBSGRBDATA

MODULE Bandmodel
  IMPLICIT NONE
  DOUBLE PRECISION, PARAMETER :: alpha=-1.1,beta=-2.3
END MODULE Bandmodel