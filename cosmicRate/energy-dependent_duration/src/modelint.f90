FUNCTION modelintz(z)
  USE MODELparameters, ONLY: sqrt2lisosigma,loglisomean
  USE COSMOparameters
  USE Zparameters
  USE constants, ONLY: log10Mpc2cmSq4pi
  USE OBSGRBDATA
  USE detection
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) :: z
  DOUBLE PRECISION :: lumdisMPc,lumdisMPcSq,erfcc
  DOUBLE PRECISION :: modelintz,ldisWickram,delayed_rate_Belz_Li    ! function
  DOUBLE PRECISION, EXTERNAL :: probatdurz  ! function
  zplus1        = z+1.0d0
  logzplus1     = dlog10(zplus1)
  logzplus1_dur = logzplus1*0.66
  lumdisMPc     = ldisWickram(z)
  lumdisMPcSq   = lumdisMPc*lumdisMPc
  lumdisterm    = log10Mpc2cmSq4pi+dlog10(lumdisMPcSq)
  meandurz   = meandur-logzplus1_dur
  logdurzmin = logdurmin-logzplus1_dur
  logdurzmax = logdurmax-logzplus1_dur
  logepkzmin = logepkmin+logzplus1
  logepkzmax = logepkmax+logzplus1
  call qrombDur(probatdurz,logdurzmin,logdurzmax,modelintz)
  modelintz=modelintz+0.5d0*erfcc((glb_max_lpb+lumdisterm-loglisomean)/sqrt2lisosigma)
  modelintz=modelintz*delayed_rate_Belz_Li(z)
END FUNCTION modelintz

FUNCTION probatdurz(logdurz)
  USE MODELparameters
  USE COSMOparameters
  USE Zparameters
  USE detection
  IMPLICIT NONE
  DOUBLE PRECISION :: probatdurz,logdurz
  DOUBLE PRECISION :: erfcc  ! function
  DOUBLE PRECISION, EXTERNAL :: ProbatLisoGivenDurz  ! function
  pph_correction=Serfc*erfcc((logdurz-meandurz)/scaledur)
  min_lpb_at_dur=eff_min_lpb+pph_correction
  conlisomeanD=aLD+bLD*logdurz
  conepkzmeanD=aED+bED*logdurz
  aELD=conepkzmeanD-conlisomeanD*bELD
  call qrombPbol(ProbatLisoGivenDurz,min_lpb_at_dur+lumdisterm,glb_max_lpb+lumdisterm,probatdurz)
  probatdurz= probatdurz &
            / (sqrt2Pit90zsigma*dexp(((logdurz-logt90zmean)/sqrt2t90zsigma)**2))
END FUNCTION probatdurz

FUNCTION ProbatLisoGivenDurz(logliso)
  USE MODELparameters
  USE Zparameters, ONLY: lumdisterm
  USE detection
  IMPLICIT NONE
  DOUBLE PRECISION :: logliso,ProbatLisoGivenDurz,efficiency
  DOUBLE PRECISION, EXTERNAL :: EpkzProbGivenLisoDurz
  conepkzmeanLD=aELD+bELD*logliso
  log10pbol=logliso-lumdisterm
  call qrombEpk(EpkzProbGivenLisoDurz,logepkzmin,logepkzmax,ProbatLisoGivenDurz)
  ProbatLisoGivenDurz=ProbatLisoGivenDurz/&
  (sqrt2PiconlisosigmaD*dexp(((logliso-conlisomeanD)/sqrt2conlisosigmaD)**2))
END FUNCTION ProbatLisoGivenDurz

FUNCTION EpkzProbGivenLisoDurz(logepkz)
  USE MODELparameters
  USE Zparameters, ONLY: logzplus1
  USE detection
  IMPLICIT NONE
  DOUBLE PRECISION :: EpkzProbGivenLisoDurz,erfcc,logepkz,efficiency,PbolEpk2P1024ph
  efficiency = 0.5d0+0.5d0 &
             * (1.d0-erfcc((PbolEpk2P1024ph(logepkz-logzplus1,log10pbol)-pph_correction-meanthresh)/sqrt2stdevthresh))
  EpkzProbGivenLisoDurz = efficiency &
                        / (sqrt2PiconepkzsigmaLD*dexp(((logepkz-conepkzmeanLD)/sqrt2conepkzsigmaLD)**2))
END FUNCTION EpkzProbGivenLisoDurz