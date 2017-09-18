!======================================================================!
! MODULE COBWEB                                                        !
!======================================================================!
! Authors: G. Tommei, F. Spoto, A. Del Vigna, A. Milani                !
! Date of creation: Nov. 2004                                          !
! Last update: Feb 2017 (F. Spoto, A. Del Vigna)                       !
!======================================================================!
! Subroutines contained:                                               !
!    - spider                                                          !
!    - mov                                                             !
!    - sample_web_grid                                                 !
!    - immediate_imp                                                   !
!    - mc_res_att                                                      !
!======================================================================!
MODULE cobweb
  USE output_control
  USE fund_const
  USE astrometric_observations
  USE orbit_elements
  USE least_squares
  USE name_rules

  IMPLICIT NONE

  PRIVATE
  !*******************************!
  !  PUBLIC PARAMETERS/VARIABLES  !
  !*******************************!
  INTEGER, PARAMETER, PUBLIC :: nmax=300                  ! Maximum number of points along one axis
  INTEGER, PARAMETER, PUBLIC :: nmax2=nmax*nmax           ! Maximum number of points in the grid
  
  INTEGER,            PUBLIC :: ndir                      ! Number of directions (rows)
  INTEGER,            PUBLIC :: np                        ! Number of points per direction (columns)
  INTEGER,            PUBLIC :: ndir2                     ! Number of directions in the densified grid (row)
  INTEGER,            PUBLIC :: np2                       ! Number of points per direction in the densified grid (columns)
  LOGICAL,            PUBLIC :: cobflag                   ! Flag: true if the MOV is available (used in fitobs)
  DOUBLE PRECISION,   PUBLIC :: sigx                      ! Max of chi for impact search (now set to 5.d0 in fitobs.def)
  DOUBLE PRECISION,   PUBLIC :: hmax                      ! Max magnitude for interesting object (now 34.5)
  LOGICAL,            PUBLIC :: use_nominal               ! If true, do cobweb; if false, grid is used
  LOGICAL,            PUBLIC :: select_risk               ! True when an orbit with cov. is available and geoc_chi > 1
  DOUBLE PRECISION,   PUBLIC :: xogeo(3)                  ! Geocentric observer positiob (Equatorial)
  DOUBLE PRECISION,   PUBLIC :: vogeo(3)                  ! Geocentric observer velocity (Equatorial)
  DOUBLE PRECISION,   PUBLIC :: eearth(0:nmax,0:nmax)     ! Two-body energy w.r. to Earth
  DOUBLE PRECISION,   PUBLIC :: esun(0:nmax,0:nmax)       ! Two-body energy w.r. to Sun
  DOUBLE PRECISION,   PUBLIC :: elong                     ! Elongation from the Sun (output of the AR computation)
  !*******************!
  !  MOV computation  !
  !*******************!
  TYPE(orbit_elem),   PUBLIC :: elcob(0:nmax,0:nmax)      ! Corrected elements on MOV (output of mov)
  TYPE(orb_uncert),   PUBLIC :: unc6(0:nmax,0:nmax)       ! 6-dim uncertainty matrices (output of mov)
  DOUBLE PRECISION,   PUBLIC :: rms_points(0:nmax,0:nmax) ! RMS after correction of attributable
  LOGICAL         ,   PUBLIC :: succ_cob(0:nmax,0:nmax)   ! Success flag for the MOV computation
  DOUBLE PRECISION,   PUBLIC :: chi_cob(0:nmax,0:nmax)    ! Chi after the MOV computation
  DOUBLE PRECISION,   PUBLIC :: jacr(0:nmax,0:nmax)       ! Jacobian det. for the map from the sampling space to the AR
  DOUBLE PRECISION,   PUBLIC :: jacarho(0:nmax,0:nmax)    ! Determinant of C_rho/Gamma_rho
  DOUBLE PRECISION,   PUBLIC :: jacmu(0:nmax,0:nmax)      ! Determinant of the map from the AR to the MOV
  !*******************************!
  !  Imminent impactors analysis  !
  !*******************************!
  DOUBLE PRECISION,   PUBLIC :: tcla(0:nmax,0:nmax)       ! Time of close approach of the (i,j) solution
  DOUBLE PRECISION,   PUBLIC :: imppro(0:nmax,0:nmax)     ! Impact probability 
  ! Non public variables (only public in this module)
  DOUBLE PRECISION           :: csi(0:nmax,0:nmax)        ! Xi coordinate on the MTP
  DOUBLE PRECISION           :: zeta(0:nmax,0:nmax)       ! Zeta coordinate on the MTP
  DOUBLE PRECISION           :: vel(0:nmax,0:nmax)        ! Velocity at closest point

  !*******************!
  !  PUBLIC ROUTINES  !
  !*******************!
  PUBLIC :: spider, mov, sample_web_grid, immediate_imp, mc_res_att, compute_score

CONTAINS

  !============================================================!                  
  ! SAMPLE_WEB_GRID                                            !
  !============================================================!
  ! Sampling of admissible region, either by cobweb or by grid !
  !============================================================!
  SUBROUTINE sample_web_grid(astnam,ini0,cov0,conditional,elc,unc, &
       & csinor,m,obs,obsw,nd,artyp,geoc_chi,acce_chi,chi)
    USE attributable 
    USE arc_control 
    USE triangles
    USE station_coordinates, ONLY: statcode
    USE planet_masses,       ONLY: gmearth
    !=======================================================================================================
    CHARACTER*(name_len), INTENT(IN)    :: astnam      ! Asteroid name
    LOGICAL,              INTENT(IN)    :: ini0        ! Availability of orbital elements
    LOGICAL,              INTENT(IN)    :: cov0        ! Availability of covariance
    LOGICAL,              INTENT(IN)    :: conditional ! Conditional probability on the MOV
    TYPE(orbit_elem),     INTENT(INOUT) :: elc         ! Nominal elements
    TYPE(orb_uncert),     INTENT(IN)    :: unc         ! Uncertainty matrices
    DOUBLE PRECISION,     INTENT(IN)    :: csinor      ! Residuals norm
    INTEGER,              INTENT(IN)    :: m           ! Number of obervations
    TYPE(ast_obs),        INTENT(IN)    :: obs(m)      ! Observations
    TYPE(ast_wbsr),       INTENT(INOUT) :: obsw(m)     ! Weights,residuals
    INTEGER,              INTENT(IN)    :: nd          ! Dimension of the par space
    INTEGER,              INTENT(OUT)   :: artyp       ! Arc type
    DOUBLE PRECISION,     INTENT(OUT)   :: geoc_chi    ! Normalized value of the geodetic curvature
    DOUBLE PRECISION,     INTENT(OUT)   :: acce_chi    ! Normalized value for the along-track acceleration
    DOUBLE PRECISION,     INTENT(OUT)   :: chi         ! Chi-value of geodetic curvature and acceleration
    !===== Arc type computation ============================================================================
    TYPE(attrib)                :: attrc                      ! Attributable
    LOGICAL                     :: error                      ! Error in the attributable computation
    DOUBLE PRECISION            :: trou                       ! Rounded time for the attributable (at the integer MJD)
    INTEGER                     :: nigarc                     ! Number of nights for the arc type
    INTEGER                     :: fail_arty                  ! Arc type computation failure
    CHARACTER(LEN=3)            :: staz                       ! Observatory code (from attributable computation)
    DOUBLE PRECISION, PARAMETER :: sphx=2.d0                  ! Maximum arc span (in degrees)
    !===== Admissible Region computation ===================================================================
    DOUBLE PRECISION, PARAMETER :: a_max=100.d0               ! Maximum value for the semimajor axis
    DOUBLE PRECISION            :: refs(3,3)                  ! Matrix of the reference system
    DOUBLE PRECISION            :: E_bound                    ! Bound for E_Sun
    DOUBLE PRECISION            :: c(0:5)                     ! Coefficients 
    INTEGER                     :: nrootsf                    ! Number of roots of the 6 deg. poly. (V(r), comments below)
    DOUBLE PRECISION            :: roots(3)                   ! Roots of the polynomial V(r)
    INTEGER                     :: nrootsd                    ! Number of roots of the pol. derivative V'(r)
    DOUBLE PRECISION            :: rootsd(8)                  ! Roots of the pol. derivative V'(r)
    DOUBLE PRECISION            :: xo(3)                      ! Heliocentric observer position, equatorial     
    DOUBLE PRECISION            :: vo(3)                      ! Heliocentric observer velocity, equatorial
    DOUBLE PRECISION            :: r1                         ! Lower bound for rho
    DOUBLE PRECISION            :: r2                         ! Upper bound for rho (first pos. root of V(r))
    INTEGER                     :: ind                        ! Index in rootsd for a zero of V'(r) between r1 and r2
    DOUBLE PRECISION            :: rder0                      ! Zero of V'(r) between r1 and r2
    INTEGER                     :: ntot                       ! Total number of points for the boundary computation
    DOUBLE PRECISION            :: frv(npox)                  ! Metric value on the boundary of the AR
    DOUBLE PRECISION            :: rhodot(npox)               ! Values of rho_dot in the boundary sampling
    DOUBLE PRECISION            :: rdw                        ! Scaling factor for range in rhodot
    !===== Cobweb/Grid computation =========================================================================
    DOUBLE PRECISION            :: rms                        ! RMS for cobweb/grid computation
    DOUBLE PRECISION            :: rmax                       ! Maximum value for rho
    DOUBLE PRECISION            :: rmin                       ! Minimum value for rho
    DOUBLE PRECISION            :: rdotmax                    ! Maximum value for rho_dot
    DOUBLE PRECISION            :: rdotmin                    ! Minimum value for rho_dot
    LOGICAL                     :: log_rho                    ! Logarithmic grid in rho
    DOUBLE PRECISION            :: rho                        ! Rho values in the grid
    INTEGER                     :: npt                        ! Number of points used in the cobweb computation
    LOGICAL                     :: found_cob(0:nmax,0:nmax)   ! Flag which is true if e < 1 and rho > 0
    DOUBLE PRECISION            :: theta(nmax)                ! Array of angles (for cobweb)
    DOUBLE PRECISION            :: r_cob(nmax)                ! Array of distances (for cobweb)
    DOUBLE PRECISION            :: jacthr(0:nmax,0:nmax)      ! Jacobian of (R,theta) -> (rho,rho_dot)
    DOUBLE PRECISION            :: r_rdot(0:nmax,0:nmax,2)    ! Coordinates of the cobweb/grid points
    DOUBLE PRECISION            :: ecc(0:nmax,0:nmax)         ! Eccentricity of each point
    DOUBLE PRECISION            :: q(0:nmax,0:nmax)           ! Pericenter of each point
    DOUBLE PRECISION            :: qg(0:nmax,0:nmax)          ! Apocenter of each point
    DOUBLE PRECISION            :: enne(0:nmax,0:nmax)        ! Mean motion of each point
    INTEGER                     :: succ(0:nmax,0:nmax)        ! Flag for success 
    DOUBLE PRECISION            :: score_nea                  ! Score point: NEA
    DOUBLE PRECISION            :: score_mba                  ! Score point: MBA 
    DOUBLE PRECISION            :: score_distant              ! Score point: TNO
    DOUBLE PRECISION            :: lim                        ! Limit value for the energy 
    !===== MOV computation =================================================================================
    DOUBLE PRECISION            :: delnor_cob(0:nmax,0:nmax)  ! Diff. corrections norm
    !===== File variables ==================================================================================
    CHARACTER(LEN=100)          :: file                       ! File name
    INTEGER                     :: le                         ! Lenght of file
    INTEGER                     :: iunfla                     ! Unit for the att file
    INTEGER                     :: iunpol                     ! Unit for the pol file (AR computation output)
    INTEGER                     :: iuncob                     ! Unit for the web file (info on the sampling points)
    INTEGER                     :: iuncom                     ! Unit for the com file (info on the sampling points)
    INTEGER                     :: iunele                     ! Unit for the kepele file (KEP el. for the sampling points)
    INTEGER                     :: iunscore                   ! Unit for the score file
    !===== Functions =======================================================================================
    DOUBLE PRECISION            :: prscal                     ! Scalar product
    DOUBLE PRECISION            :: vsize                      ! Norm of a vector
    !===== Other variables =================================================================================
    DOUBLE PRECISION            :: crho(2,2)                  ! Submatrix of C for (rho,rho_dot)
    DOUBLE PRECISION            :: grho(2,2)                  ! Submatrix of G for (rho,rho_dot)
    DOUBLE PRECISION            :: att(4)                     ! Attributable for the energy computation
    TYPE(orbit_elem)            :: elk                        ! Keplerian elements of the sampling points
    TYPE(orbit_elem)            :: elcom                      ! Cometarian elements of the sampling points
    INTEGER                     :: fail_flag                  ! Failure of the coordinate change
    CHARACTER(LEN=20)           :: mulname                    ! Name of the sampling point
    INTEGER                     :: i, j, k, ii, jj, iii       ! Loop indexes
    !=======================================================================================================
    IF(nd.NE.6)THEN
       WRITE(*,*) 'sample_web_grid: not ready for non-grav, nd = ',nd
       STOP
    END IF
    !******************!
    !  Inizialization  !
    !******************!
    ! In the grid case, r_cob and theta are not computed
    r_cob = 0
    theta = 0

    !****************************!
    !  COMPUTE THE ATTRIBUTABLE  !
    !****************************!
    CALL tee(iun_log,'------------------------=')
    CALL tee(iun_log,'| COMPUTE ATTRIBUTABLE |=')
    CALL tee(iun_log,'------------------------=')
    CALL attri_comp(m,obs,obsw,attrc,error)
    IF(error)THEN
       WRITE(*,*) 'sample_web_grid: error in attri_comp'
       STOP
    END IF
    trou=NINT(attrc%tdtobs)
    ! Write the attributable computation output
    IF(attrc%sph*radeg.GT.sphx)THEN
       WRITE(*,*) 'sample_web_grid: arc too wide = ', attrc%sph*radeg
    ELSE
       CALL wri_attri(0,0,astnam,attrc,trou)
       CALL wri_attri(iun_log,iun_log,astnam,attrc,trou)
    ENDIF

    !************************!
    !  COMPUTE THE ARC TYPE  !
    !************************!
    CALL tee(iun_log,'--------------------=')
    CALL tee(iun_log,'| COMPUTE ARC TYPE |=')
    CALL tee(iun_log,'--------------------=')
    artyp=arc_type(obs,obsw,m,geoc_chi,acce_chi,chi,nigarc,fail_arty)
    WRITE(*,176) artyp
    WRITE(*,177) m,nigarc,fail_arty,geoc_chi,acce_chi,chi
    WRITE(iun_log,176) artyp
    WRITE(iun_log,177) m,nigarc,fail_arty,geoc_chi,acce_chi,chi
176 FORMAT(' The arc type is ',I3)
177 FORMAT(' nobs = ',I4,', nights = ',I3,', errcode = ',I4/ &
         & ' geoc_chi = ',1P,D9.2,', acce_chi = ',D9.2,', chi = ',D9.2)
    
    !*********************************!
    !  COMPUTE THE ADMISSIBLE REGION  !
    !*********************************!
    !------------------------------------------------------------------!
    ! The most important condition defining the AR is the              !
    ! non-positivity of the two body energy of the body w.r.t the      !
    ! Sun. This condition gives solutions for r_dot if and only if the !
    ! following inequality holds:                                      !
    !                                                                  !
    !                          V(r) <= 4k^4                            !
    !                                                                  !
    ! where V(r) is a 6 degree polynomial. We can have only three      !
    ! possibilities:                                                   !
    !                                                                  !
    ! (1) V(r) has 4 simple real roots -> the AR has 2 connected       !
    !     components;                                                  !
    !                                                                  !
    ! (2) V(r) has 3 distinct real roots, 2 simple and 1 with even     !
    !     multiplicity -> the second c.c. reduces to a point;          !
    !                                                                  !
    ! (3) V(r) has 2 distinct real roots, 1 simple and 1 with odd      !
    !     multiplicity -> the AR has 1 c.c.                            !
    !                                                                  !
    ! Moreover, V'(r) cannot have more than 3 distinct roots. If they  !
    ! are exactly three, then there cannot be any root with            ! 
    ! multiplicity 2.                                                  !
    !------------------------------------------------------------------!
    !***************************************************!
    !  Open .att file for output of the AR computation  !
    !***************************************************!
    file = astnam//'.att'
    CALL rmsp(file,le)
    CALL filopn(iunfla,file(1:le),'unknown')
    staz = attrc%obscod
    CALL tee(iun_log,'-----------------------------=')
    CALL tee(iun_log,'| COMPUTE ADMISSIBLE REGION |=')
    CALL tee(iun_log,'-----------------------------=')
    CALL admis_reg(astnam,iunfla,attrc%tdtobs,attrc%angles,staz,nrootsf, &
         & roots,c,refs,xo,vo,a_max,E_bound,nrootsd,rootsd,attrc%apm,xogeo,vogeo,elong)
    CALL filclo(iunfla,' ')
    ! Output the roots
    WRITE(*,*) 'Roots = ', roots(1:nrootsf) 
    WRITE(iun_log,*) 'Roots = ', roots(1:nrootsf) 
    !****************************************!
    !  Find the range of values for rho_dot  !
    !****************************************!
    r2 = roots(1) ! First positive root of V(r)
    r1 = 1.d-4    ! Arbitrary lower boundary for the AR
    ! Search for a zero (rder0) of the derivative V'(r) between r1 and r2 
    rder0 = -1.d0
    ind   = 0
    IF(nrootsd.GE.1) THEN
       IF(nrootsd.EQ.1) THEN
          ! Do nothing
       ELSEIF(nrootsd.EQ.2) THEN
          ! Do nothing
       ELSEIF(nrootsd.EQ.3) THEN
          DO iii=1,nrootsd
             IF(rootsd(iii).GT.r1 .AND. rootsd(iii).LT.r2) THEN
                ind=iii
             ENDIF
          ENDDO
          IF(ind.GT.2)THEN
             IF(rootsd(ind-2).GT.r1) THEN
                rder0=rootsd(ind-2)
                WRITE(*,*) 'r1 = ',r1,', rder0 = ',rder0,', r2 = ',r2
             ENDIF
          ENDIF
          IF(rder0.GE.r2) THEN
             WRITE(*,*) 'sample_web_grid: ERROR! rder0 > r2'
          ENDIF
       ELSE
          ! V'(r) cannot have more than 3 distinct roots
          WRITE(*,*) 'sample_web_grid: roots of the der. of V(r) (nrootsd) = ',nrootsd
          WRITE(*,*) '                 It cannot be possible!'
       ENDIF
    ENDIF
    nfunc = 2 ! Use logaritmic scale for the metric used in the boundary sampling
    CALL sample_bound_ne1(E_bound,attrc%eta**2,r1,rder0,r2,8,5, &
         & ntot,c,frv,rhodot,npox,rdw)
    rdotmin = MINVAL(rhodot(1:ntot))
    rdotmax = MAXVAL(rhodot(1:ntot))
    !*************************************************!
    !  Open .pol file, containing the data of the AR  !
    !*************************************************!
    file=astnam//'.pol'
    CALL rmsp(file,le)
    CALL filopn(iunpol,file(1:le),'unknown')
    ! Write header
    WRITE(iunpol,'(A)') '% N. roots      Roots(1)              Roots(2)              Roots(3)              &
         &c(0)                 c(1)                c(2)                 c(3)                 c(4)          &
         &       c(5)                rdotmin              rdotmax'

    WRITE(iunpol,178) nrootsf, roots, c, rdotmin, rdotmax
178 FORMAT(6X,I1,2X,11(1X,F20.16)) 
    CALL filclo(iunpol,' ')

    !***************************!
    !  COBWEB/GRID COMPUTATION  !
    !***************************!
    !***************************!
    !  Selection of the method  !
    !***************************!
    !------------------------------------------------------------------!
    ! The sampling of the AR depends on the availability of a nominal  !
    ! solution. If there is a nominal solution, with its covariance,   !
    ! and if the SNR of the geodetic curvature is > 3, then we trust   !
    ! the nominal solution and we do a cobweb sampling.                !
    !------------------------------------------------------------------!
    use_nominal = (ini0 .AND. cov0 .AND. ABS(geoc_chi).GT.3.d0) ! Flag for the risk type selection (used in immediate_imp)
    select_risk = (ini0 .AND. cov0 .AND. ABS(geoc_chi).GT.1.d0) ! True when an orbit with cov. is available and geoc_chi > 1
    DO k=1,2
       IF(k.EQ.1) THEN
          IF(use_nominal)THEN
             ! Spider 100x100
             np = np2
             ndir = ndir2
             !*******************!
             !  Cobweb sampling  !
             !*******************!
             IF(elc%coo.ne.'ATT')THEN
                WRITE(*,*) ' sample_web_grid: ERROR! Elements must be in ATT, they are ', elc%coo
                STOP
             ENDIF
             CALL tee(iun_log,'------------------=')
             CALL tee(iun_log,'| COMPUTE COBWEB |=')
             CALL tee(iun_log,'------------------=')
             CALL spider(elc,unc,csinor,npt,r_rdot,theta,r_cob,ecc,q,qg,enne,found_cob,jacthr)
             rms=csinor
          ELSE
             log_rho = .FALSE.       ! Sampling of the AR
             succ = 0                ! Succ flag for writing
             rms = -1.d0             ! RMS fixed value
             lim = -gms/(2*a_max)    ! Energy curve level
             !*****************!
             !  Grid sampling  !
             !*****************!
             CALL tee(iun_log,'----------------=')
             CALL tee(iun_log,'| COMPUTE GRID |=')
             CALL tee(iun_log,'----------------=')
             !--------------------------------------------------------!
             ! Case with two connnected components handled (roots(3)) !
             !--------------------------------------------------------!
             IF(nrootsf.EQ.3) THEN
                rmax = roots(3)
                ndir = ndir2
                np = np2
             ELSE
                rmax=roots(1)
             END IF
             ! Grid starting point
             rmin=r1
             IF(rmax.LT.10**(0.5d0)) log_rho = .TRUE.
             CALL sample_grid(rmin,rmax,rdotmin,rdotmax,log_rho,r_rdot)
          END IF
          ! Excluding orbits with e greater than or equal to 1
          elc=undefined_orbit_elem
          elc%coo='ATT'                                   ! ATT coordinates
          elc%t=attrc%tdtobs                              ! Time of the attributable
          CALL statcode(attrc%obscod,elc%obscode)         ! Observatory code
          elc%coord(1:4)=attrc%angles                     ! Attributable for the first 4 coordinates
          
          !***************************!
          !  ENERGY W.R.T. THE SUN    !
          !***************************!
          CALL energy_sun(r_rdot(0:ndir,0:np,1:2),c,ndir,np,esun(0:ndir,0:np))
          DO i=1,ndir
             DO j=1,np
                IF(esun(i,j).LT.lim)THEN
                   found_cob(i,j)=.TRUE.
                ELSE
                   found_cob(i,j)=.FALSE.
                   succ(i,j) = 2
                ENDIF
             ENDDO
          ENDDO
          elc%coord(5:6)=(/0.d0, 0.d0/) ! Restore coordinates (0,0) for (rho,rho_dot)
          
          !**************************!
          !  MANIFOLD OF VARIATIONS  !
          !**************************!
          CALL mov(m,obs,obsw,elc,rms,r_rdot,found_cob,nd,elcob,unc6, &
               & rms_points,delnor_cob,succ_cob,chi_cob)
          cobflag=.TRUE.
          
          !*************************!
          !  JACOBIAN DETERMINANTS  !
          !*************************!
          CALL compute_determinants(use_nominal,log_rho,r_rdot,conditional,succ_cob,&
               &unc6,unc,jacthr,jacarho)
          
          !***************************!
          !  ENERGY W.R.T. THE EARTH  !
          !***************************!
          CALL energy_earth(att,xogeo,vogeo,r_rdot(0:ndir,0:np,1:2),ndir,np,eearth(0:ndir,0:np))
          
          !****************!
          !  WRITE OUTPUT  !
          !****************!
          CALL write_output_files(k,use_nominal,r2,astnam,succ_cob,rms,succ,r_rdot,r_cob,theta,rms_points,csinor,chi_cob, &
               & jacmu,delnor_cob,elcob,elc,esun,eearth)
          
          !*****************!
          !  COMPUTE_SCORE  !
          !*****************!
          CALL compute_score(score_distant, score_mba, score_nea)
!**************************************************************************************************************
       ELSEIF(k.EQ.2.AND..NOT.use_nominal) THEN
          log_rho = .FALSE.      ! Sampling of the AR
          succ = 0               ! Succ flag for writing
          ndir = ndir2           ! New number of points in rho
          np = np2               ! New number of points in rho_dot
! Initialize AR sampling boundaries
          rmin    = 1000.d0
          rmax    = 0.d0
          rdotmin = 1000.d0
          rdotmax = -1000.d0
          CALL tee(iun_log,'--------------------------=')
          CALL tee(iun_log,'| COMPUTE DENSIFIED GRID |=')
          CALL tee(iun_log,'--------------------------=')
          DO i=1,ndir
             DO j=1,np
                CALL coo_cha(elcob(i,j),'KEP',elk,fail_flag)
                CALL coo_cha(elcob(i,j),'COM',elcom,fail_flag)
                IF(chi_cob(i,j).LT.5.d0.AND.elk%h_mag.LT.34.5.AND.succ_cob(i,j)) THEN
                   IF(r_rdot(i,j,1).LT.rmin) THEN
                      rmin=r_rdot(i,j,1)
                   END IF
                   IF(r_rdot(i,j,1).GT.rmax) THEN
                      rmax=r_rdot(i,j,1)
                   END IF
                   IF(r_rdot(i,j,2).LT.rdotmin) THEN
                      rdotmin=r_rdot(i,j,2)
                   END IF
                   IF(r_rdot(i,j,2).GT.rdotmax) THEN
                      rdotmax=r_rdot(i,j,2)
                   END IF
                END IF
             END DO
          END DO
          IF(rmin.EQ.1000.d0.OR.rmax.EQ.0.d0.OR.rdotmin.EQ.1000.d0.OR.rdotmax.EQ.-1000.d0)THEN
             WRITE(*,*) 'cobweb: sample_web_grid, no good points in the previous iteration'
             EXIT
          END IF
          ! If rmin and rmax are equal (vertical line)
          IF(rmin.EQ.rmax)THEN
             rmin = rmin - 0.05
             rmax = rmax + 0.05
          END IF
          ! If rdotmin and rdotmax are equal (horizontal line)
          IF(rdotmin.EQ.rdotmax)THEN
             rdotmin = rmin - 0.005
             rdotmax = rmax + 0.005
          END IF
          IF(score_nea.GT.0.5d0) log_rho = .TRUE.
          CALL sample_grid(rmin,rmax,rdotmin,rdotmax,log_rho,r_rdot)
          elc=undefined_orbit_elem
          elc%t=attrc%tdtobs
          elc%coo='ATT'
          CALL statcode(attrc%obscod,elc%obscode)
          elc%coord(1:4)=attrc%angles

          !***************************!
          !  ENERGY W.R.T. THE SUN    !
          !***************************!
          CALL energy_sun(r_rdot(0:ndir,0:np,1:2),c,ndir,np,esun(0:ndir,0:np))

          ! check that belongs to attributable region: if not, cut from MOV           
          DO i=1,ndir
             DO j=1,np
                 IF(esun(i,j).LT.lim)THEN
                   found_cob(i,j)=.true.
                ELSE
                   found_cob(i,j)=.false.
                   succ(i,j) = 2
                ENDIF
             ENDDO
          ENDDO
          elc%coord(5:6)=(/0.d0, 0.d0/)

          !**************************!
          !  MANIFOLD OF VARIATIONS  !
          !**************************!
          CALL mov(m,obs,obsw,elc,rms,r_rdot,found_cob,nd,elcob,unc6,rms_points,delnor_cob,succ_cob,chi_cob)
          cobflag=.true.

          !*************************!
          !  JACOBIAN DETERMINANTS  !
          !*************************!
          CALL compute_determinants(use_nominal,log_rho,r_rdot,conditional,succ_cob,&
               &unc6,unc,jacthr,jacarho)

          !***************************!
          !  ENERGY W.R.T. THE EARTH  !
          !***************************!
          CALL energy_earth(attrc%angles,xogeo,vogeo,r_rdot(0:ndir,0:np,1:2),ndir,np,eearth(0:ndir,0:np))
    
          !********************!
          !  WRITE THE OUTPUT  !
          !********************!
          CALL write_output_files(k,use_nominal,r2,astnam,succ_cob,rms,succ,r_rdot,r_cob,theta,rms_points,csinor,chi_cob, &
               & jacmu,delnor_cob,elcob,elc,esun,eearth)
       END IF
    END DO

    !*****************!
    !  COMPUTE_SCORE  !
    !*****************!
    CALL compute_score(score_distant, score_mba, score_nea)
    file = astnam//'.score'
    CALL rmsp(file,le)
    CALL filopn(iunscore,file(1:le),'unknown')
    ! Write header score file
    WRITE(iunscore,'(A)') '% Object score'
    WRITE(iunscore, 595) 'Score Near-Earth: ', score_nea
    WRITE(iunscore, 595) 'Score Main Belt : ', score_mba
    WRITE(iunscore, 595) 'Score Distant   : ', score_distant
595 FORMAT(A, 3X, F4.2)
    CALL filclo(iunscore,' ')
  END SUBROUTINE sample_web_grid

  !============================================================!                  
  ! SAMPLE_WEB                                                 !
  !============================================================!
  ! Sampling of admissible region, either by cobweb or by grid !
  !============================================================!
  SUBROUTINE sample_grid(rmin, rmax, rdotmin, rdotmax, log_grid, r_rdot)
    !==========================================================!
    ! Begin interface
    DOUBLE PRECISION, INTENT(IN)  :: rmin                      ! Starting point for the rho sampling
    DOUBLE PRECISION, INTENT(IN)  :: rmax                      ! Ending point for the rho sampling
    DOUBLE PRECISION, INTENT(IN)  :: rdotmin                   ! Starting point for the rhodot sampling
    DOUBLE PRECISION, INTENT(IN)  :: rdotmax                   ! Ending point for the rhodot sampling
    LOGICAL,          INTENT(IN)  :: log_grid                  ! Logarithmic grid in rho
    DOUBLE PRECISION, INTENT(OUT) :: r_rdot(0:nmax,0:nmax,2)   ! Rho and rhodot in the grid
    ! End interface
    !==========================================================!
    DOUBLE PRECISION :: rms                                    ! RMS value
    DOUBLE PRECISION :: dr                                     ! Sampling step in rho
    DOUBLE PRECISION :: drdot                                  ! Sampling in rhodot
    INTEGER          :: i,j                                    ! Loop indices
    !=======================================================================================================
    ! Initializing
    ! Succ flag (integer): 0 if differential corrections fail, 2 if not in the admissible region, 1 if good
    r_rdot = 0.d0
    !*****************!
    !  Grid sampling  !
    !*****************!
    !---------------------------------------------------------------!
    ! ndir and np are used as number of values of rho and rho_dot
    ! Uniform grid in rho and rhodot
    ! Old version: dr=(log10(rmax)+4.d0)/ndir, r_rdot(i,1:np,1)=10**(dr*i-4.d0) 
    !---------------------------------------------------------------!
    IF(log_grid)THEN
       dr=(log10(rmax)-log10(rmin))/ndir
    ELSE
       dr=(rmax-rmin)/(ndir-1)
    END IF
    drdot=(rdotmax-rdotmin)/(np-1)
    DO i=1,ndir
       IF(log_grid) THEN
          r_rdot(i,1:np,1)=10**(dr*i+log10(rmin))
       ELSE
          r_rdot(i,1:np,1)=dr*(i-1)+rmin
       END IF
    ENDDO
    DO j=1,np
       r_rdot(1:ndir,j,2)=rdotmin+drdot*(j-1)
    ENDDO
  END SUBROUTINE sample_grid


  !============================================================!                  
  ! COMPUTE_DETERMINANTS                                       !
  !============================================================!
  ! Jacobian determinant                                       !
  !============================================================!
  SUBROUTINE compute_determinants(nominal,log_grid,r_rdot,conditional,succ_cobw,unc_6,unc,jacthr,jac_rho)
    !========================================================!
    ! Begin interface                                        !
    LOGICAL,          INTENT(IN)        :: nominal                    ! Use nominal
    LOGICAL,          INTENT(IN)        :: log_grid                   ! Uniform sampling of the grid (in log or in rho)
    DOUBLE PRECISION, INTENT(IN)        :: r_rdot(0:nmax,0:nmax,2)    ! Coordinates of the cobweb/grid points
    LOGICAL,          INTENT(IN)        :: conditional                ! Conditional probability on the MOV
    LOGICAL,          INTENT(IN)        :: succ_cobw(0:nmax,0:nmax)   ! Success flag for the MOV computation
    TYPE(orb_uncert), INTENT(IN)        :: unc_6(0:nmax,0:nmax)       ! 6-dim uncertainty matrices (output of mov)
    TYPE(orb_uncert), INTENT(IN)        :: unc                        ! Uncertainty matrices
    DOUBLE PRECISION, INTENT(IN)        :: jacthr(0:nmax,0:nmax)      ! Jacobian of (R,theta) -> (rho,rho_dot)
    DOUBLE PRECISION, INTENT(OUT)       :: jac_rho(0:nmax,0:nmax)     ! Determinant of C_rho/Gamma_rho
    ! End interface
    !===================================================!
    DOUBLE PRECISION                    :: crho(2,2)                  ! Submatrix of C for (rho,rho_dot)
    DOUBLE PRECISION                    :: grho(2,2)                  ! Submatrix of G for (rho,rho_dot)
    INTEGER                             :: i, j                       ! Loop indices
    !===================================================!
    !*************************!
    !  JACOBIAN DETERMINANTS  !
    !*************************!
    !-------------------------------------------------------------------!
    ! To compute the probability density function on the sampling       !
    ! space induced by the astrometric errors we have to consider some  !
    ! maps:                                                             !
    !                                                                   !
    ! 1) from the sampling space to the AR. The sampling space can be   !
    !    (R,theta) or (Log(rho),rho_dot): in the first case we use      !
    !    jacthr, in the second one we use log(10)*rho                   !
    !                                                                   !
    ! 2) from the AR to the MOV                                         !
    !                                                                   !
    ! 3) from the MOV to the possible residuals space                   !
    !-------------------------------------------------------------------!

    !*************!
    !  First map  !
    !*************!
    IF(nominal)THEN
       ! Cobweb sampling, jacobian computed in spider (jacthr)
       jacr(1:ndir,1:np) = jacthr(1:ndir,1:np)
       jacr(0,0)         = jacthr(0,0)
    ELSE
       ! Semilogarithmic grid used (no nominal available)
       IF(log_grid) THEN
          jacr(1:ndir,1:np) = log(10.d0)*r_rdot(1:ndir,1:np,1)
       ELSE
          jacr(1:ndir,1:np) = 1.d0
       END IF
    ENDIF
    !**************!
    !  Second map  !
    !**************!
    !-----------------------------------------------------------------!
    ! Second order equations implemented                              !
    !-----------------------------------------------------------------!
    IF(conditional)THEN
       ! Use the conditional normal matrix from the MOV computation
       DO i=1,ndir
          DO j=1,np
             IF(succ_cobw(i,j))THEN
                crho         = unc_6(i,j)%c(5:6,5:6)
                jac_rho(i,j) = SQRT(crho(1,1)*crho(2,2)-crho(1,2)*crho(2,1))
             ELSE
                jac_rho(i,j) = 0.d0
             ENDIF
          END DO
       END DO
       IF(nominal)THEN
          crho               = unc%c(5:6,5:6)
          jac_rho(0,0)       = SQRT(crho(1,1)*crho(2,2)-crho(1,2)*crho(2,1))
       ENDIF
    ELSE
       ! hypothesis 2: use marginal normal matrix from MOV computation
       DO i=1,ndir
          DO j=1,np
             IF(succ_cobw(i,j))THEN
                grho           = unc6(i,j)%g(5:6,5:6)
                jac_rho(i,j)   = 1.d0/SQRT(grho(1,1)*grho(2,2)-grho(1,2)*grho(2,1))
             ELSE
                jac_rho(i,j)   = 0.d0
             ENDIF
          ENDDO
       ENDDO
       IF(nominal)THEN
          grho                 = unc%g(5:6,5:6)
          jac_rho(0,0)         = 1.d0/SQRT(grho(1,1)*grho(2,2)-grho(1,2)*grho(2,1))
       ENDIF
    ENDIF
  END SUBROUTINE compute_determinants

  !============================================================!                  
  ! WRITE_OUTPUT_FILES                                         !
  !============================================================!
  ! Write output files: web and element files                  !
  !============================================================!
  SUBROUTINE write_output_files(grid,use_nominal,root_2,astnam,succ_cobw,rms,succ,r_rdot,r_cob,theta,rms_points,csinor,chi_cobw, &
       & jac_mu,delnor_cob,el_cob,elc,e_sun,e_earth)
    !===========================================================!
    ! Begin interface                     
    INTEGER,              INTENT(IN)    :: grid                      ! Grid number (first or second)
    LOGICAL,              INTENT(IN)    :: use_nominal               ! Flag to write .web file even with grid.eq.1
    DOUBLE PRECISION,     INTENT(IN)    :: root_2                    ! Second root
    CHARACTER*(name_len), INTENT(IN)    :: astnam                    ! Asteroid name
    LOGICAL,              INTENT(IN)    :: succ_cobw(0:nmax,0:nmax)  ! Success flag for the MOV computation
    DOUBLE PRECISION,     INTENT(IN)    :: rms                       ! RMS
    INTEGER,              INTENT(INOUT) :: succ(0:nmax,0:nmax)       ! Flag for success
    DOUBLE PRECISION,     INTENT(IN)    :: r_rdot(0:nmax,0:nmax,2)   ! Coordinates of the cobweb/grid points
    DOUBLE PRECISION,     INTENT(IN)    :: r_cob(nmax)               ! Array of distances (for cobweb)
    DOUBLE PRECISION,     INTENT(IN)    :: theta(nmax)               ! Array of distances (for cobweb)
    DOUBLE PRECISION,     INTENT(IN)    :: rms_points(0:nmax,0:nmax) ! RMS of corrected element
    DOUBLE PRECISION,     INTENT(IN)    :: csinor                    ! Residuals norm
    DOUBLE PRECISION,     INTENT(IN)    :: chi_cobw(0:nmax,0:nmax)   ! Chi after the MOV computation
    DOUBLE PRECISION,     INTENT(IN)    :: jac_mu(0:nmax,0:nmax)     ! Determinant of the map from the AR to the MOV
    DOUBLE PRECISION,     INTENT(IN)    :: delnor_cob(0:nmax,0:nmax) ! Diff. corrections norm
    TYPE(orbit_elem),     INTENT(IN)    :: el_cob(0:nmax,0:nmax)     ! Corrected elements on MOV (output of mov)
    TYPE(orbit_elem),     INTENT(IN)    :: elc                       ! Nominal elements
    DOUBLE PRECISION,     INTENT(IN)    :: e_earth(0:nmax,0:nmax)    ! Two-body energy w.r. to Earth
    DOUBLE PRECISION,     INTENT(IN)    :: e_sun(0:nmax,0:nmax)      ! Two-body energy w.r. to Sun
    ! End interface
    !===========================================================!
    INTEGER                     :: i,j                        ! Loop indices
    CHARACTER(LEN=100)          :: file_web                   ! Cobweb file
    CHARACTER(LEN=100)          :: file_kep                   ! Keplerian file
    CHARACTER(LEN=100)          :: file_com                   ! Cometarian file
    INTEGER                     :: iuncob                     ! Unit for web file
    INTEGER                     :: iunele                     ! Unit for keplerian file
    INTEGER                     :: iuncom                     ! Unit for cometarian file
    INTEGER                     :: le                         ! File lenght
    TYPE(orbit_elem)            :: elk                        ! Keplerian elements of the sampling points
    TYPE(orbit_elem)            :: elcom                      ! Cometarian elements of the sampling points
    INTEGER                     :: fail_flag                  ! Failure of the coordinate change
    CHARACTER(LEN=20)           :: mulname                    ! Name of the sampling point 
    INTEGER                     :: region                     ! Number of components
    !========== WRITE OUTPUT =========================================
    IF(grid.EQ.2.OR.use_nominal) THEN
       file_web=astnam//'.web'
       file_kep=astnam//'.kepele'
       file_com=astnam//'.comele'
    ELSEIF(grid.EQ.1)THEN
       file_web=astnam//'_1.web'
       file_kep=astnam//'_1.kepele'
       file_com=astnam//'_1.comele'          
    END IF
    CALL rmsp(file_web,le)
    CALL filopn(iuncob,file_web(1:le),'unknown')
    CALL rmsp(file_kep,le)
    CALL filopn(iunele,file_kep(1:le),'unknown')
    CALL wro1lh_matlab(iunele,'ECLM','J2000','KEP')
    CALL rmsp(file_com,le)
    CALL filopn(iuncom,file_com(1:le),'unknown')
    CALL wro1lh_matlab(iuncom,'ECLM','J2000','COM')
    WRITE(iuncob,"(A)")'%  i    j   succ    rho         rho_dot        r          theta       RMS           &
         &chi      Jac(grid-AR)  Jac(AR-MOV)    delnor    h_mag    E_Sun        E_earth'
    IF(rms.GE.0.d0)THEN
       ! if rms of nominal not available, the nominal does not exist and
       ! should not be written
       WRITE(iuncob,130)0,0,1,r_rdot(0,0,1:2),0.d0,0.d0,csinor,chi_cobw(0,0),&
            & jac_mu(0,0),delnor_cob(0,0),elc%h_mag,e_sun(0,0),e_earth(0,0)
    ENDIF
    DO i=1,ndir
       DO j=1,np
          IF(succ_cobw(i,j)) succ(i,j) = 1
130       FORMAT(I4,1x,I4,4x,I1,1X,4(F12.7,1x),1p,d12.5,1x,d12.5,1x,d12.5,0p,2x,1p,d10.3,0p,1x,F6.3,1X,&
               &1P,D12.5,0P,1X,1P,D12.5)
          WRITE(iuncob,130)i,j,succ(i,j),r_rdot(i,j,1:2),r_cob(j),theta(i),rms_points(i,j),chi_cobw(i,j),jac_mu(i,j),&
               &delnor_cob(i,j),el_cob(i,j)%h_mag,e_sun(i,j),e_earth(i,j)
          IF(succ_cobw(i,j))THEN
             IF(root_2.GT.0.AND.r_rdot(i,j,1).GE.root_2)THEN
                region = 2
             ELSE
                region = 1
             END IF
             CALL coo_cha(el_cob(i,j),'KEP',elk,fail_flag)
             IF(fail_flag.gt.0)THEN
                WRITE(*,*)' cobweb.f90: failed conversion to kep', elk
             ENDIF
             CALL coo_cha(el_cob(i,j),'COM',elcom,fail_flag)
             IF(fail_flag.gt.0)THEN
                WRITE(*,*)' cobweb.f90: failed conversion to com', elcom
             ENDIF
             WRITE(mulname,591)i,j
591          FORMAT(I3,1X,I3)
             CALL wro1lr_matlab(iunele,mulname,elk%coord,elk%coo,elk%t,elk%h_mag,CHI=chi_cobw(i,j),&
                  & RMS=rms_points(i,j),SUCC=region)
             IF(succ_cobw(i,j))THEN
                CALL wro1lr_matlab(iuncom,mulname,elcom%coord,elcom%coo,elcom%t,elcom%h_mag,CHI=chi_cobw(i,j),&
                     &RMS=rms_points(i,j),SUCC=region)
             END IF
          END IF
       ENDDO
    ENDDO
    CALL filclo(iuncob,' ')
    CALL filclo(iunele,' ')
    CALL filclo(iuncom,' ')
  END SUBROUTINE write_output_files

  !====================================================================!
  ! IMMEDIATE_IMP                                                      !
  !====================================================================!
  ! Search for immediate impactors, given cobweb or grid               !
  !====================================================================!
  SUBROUTINE immediate_imp(astnam,tmcla,dtclo,csinor,artyp,geoc_chi,acce_chi,chi,m,obs)
    USE propag_state, ONLY: pro_ele 
    USE tp_trace,     ONLY: wri_tppoint, dtpde, tpplane, mtp_store, njc 
    INCLUDE 'vers.h90'
    !=======================================================================================================
    CHARACTER*(name_len), INTENT(IN) :: astnam   ! Asteroid name 
    DOUBLE PRECISION,     INTENT(IN) :: tmcla    ! End time for the scan
    DOUBLE PRECISION,     INTENT(IN) :: dtclo    ! Time (days) for the scan
    DOUBLE PRECISION,     INTENT(IN) :: csinor   ! Residuals norm 
    INTEGER,              INTENT(IN) :: artyp    ! Arc type
    DOUBLE PRECISION,     INTENT(IN) :: geoc_chi ! Normalized value of the geodetic curvature
    DOUBLE PRECISION,     INTENT(IN) :: acce_chi ! Normalized value for the along-track acceleration
    DOUBLE PRECISION,     INTENT(IN) :: chi      ! Chi-value of geodetic curvature and acceleration
    INTEGER,              INTENT(IN) :: m        ! Number of observations
    TYPE(ast_obs),        INTENT(IN) :: obs(m)   ! Observations
    !===== Impact Probability computation ==================================================================
    INTEGER            :: nposs       ! Number of points with chi <= 5, H <= H_max
    INTEGER            :: nimp        ! Number of impacting points
    LOGICAL            :: significant ! Significant flag
    TYPE(orbit_elem)   :: el1         ! Propagation of a sampling orbit
    DOUBLE PRECISION   :: pronorm     ! Normalization factor for the Impact Probability
    DOUBLE PRECISION   :: iptot       ! Total IP
    DOUBLE PRECISION   :: h_magmin    ! Minimum H among the sampling points
    DOUBLE PRECISION   :: h_magmax    ! Maximum H among the sampling points
    DOUBLE PRECISION   :: pearthb     ! Probability to be an Earth-bounded object
    DOUBLE PRECISION   :: t_minimp    ! Minimum impact time among the impacting points
    DOUBLE PRECISION   :: t_maximp    ! Maximum impact among the impacting points
    DOUBLE PRECISION   :: hour_min    ! Hour of t_minimp
    INTEGER            :: iday_min    ! Day of t_minimp
    INTEGER            :: imonth_min  ! Month of t_minimp
    INTEGER            :: iyear_min   ! Year of t_minimp
    DOUBLE PRECISION   :: hour_max    ! Hour of t_maximp
    INTEGER            :: iday_max    ! Day of t_maximp
    INTEGER            :: imonth_max  ! Month of t_maximp
    INTEGER            :: iyear_max   ! Year of t_maximp
    !===== Other variables =================================================================================
    CHARACTER(LEN=8)   :: date        ! Date of computation
    CHARACTER(LEN=10)  :: time        ! Time of computation
    INTEGER            :: i,j         ! Loop indexes 
    DOUBLE PRECISION   :: dt          ! Arc length in hours
    CHARACTER(LEN=100) :: file        ! File name
    INTEGER            :: le          ! Length of file
    INTEGER            :: iunimp      ! File unit for the .mtp file (MTP analysis)
    INTEGER            :: iunrisk     ! File unit for the risk file
    CHARACTER(LEN=20)  :: mulname     ! Name of the Virtual Asteroid
    INTEGER            :: dataflag    ! Data policy flag (values from 0 to 4)
    !=======================================================================================================
    ! Select the MTP plane (not the TP)
    tpplane=.FALSE.
    ! Clean up the previous close app. records
    REWIND(iuncla) 
    !**************************************!
    !  File .mtp for the MTP data storage  !
    !**************************************!
    file=astnam//'.mtp'
    CALL rmsp(file,le)
    CALL filopn(iunimp,file(1:le),'unknown')

    !***********************************!
    !  NORMALIZATION FACTOR FOR THE IP  !
    !***********************************!
    nposs=0      ! Counter of "good" points
    pronorm=0.d0 ! Normalization factor
    !*******************************!
    !  Loop on the sampling points  !
    !*******************************!
    DO i=0,ndir
       DO j=0,np
          IF(i.EQ.0 .AND. j.NE.0) CYCLE                          ! Points (0,j) (with j non-zero) does not exist
          IF(i.NE.0 .AND. j.EQ.0) CYCLE                          ! Points (i,0) (with i non-zero) does not exist
          IF(i.EQ.0 .AND. j.EQ.0 .AND. (.NOT.use_nominal)) CYCLE ! Skip the nominal if it is not available
          IF(.NOT.succ_cob(i,j) .OR. chi_cob(i,j).GT.5.d0) CYCLE ! Consider only the successful points with chi <= 5
          IF(elcob(i,j)%h_mag.LE.hmax)THEN                       ! Consider the points below the shooting star limit
             nposs=nposs+1
             pronorm=pronorm+EXP(-chi_cob(i,j)**2/2.d0)*jacr(i,j)*jacmu(i,j)!*jacarho(i,j)
          ENDIF
       ENDDO
    ENDDO

    !**********************************!
    !  COMPUTE THE IMPACT PROBABILITY  !
    !**********************************!
    !******************!
    !  Inizialization  !
    !******************!
    nimp=0
    imppro=0.d0
    iptot=0.d0
    pearthb=0.d0
    h_magmin=100
    h_magmax=0
    t_minimp=100000
    t_maximp=0
    DO i=0,ndir
       DO j=0,np
          IF(i.EQ.0 .AND. j.NE.0) CYCLE                          ! Points (0,j) (with j non-zero) does not exist
          IF(i.NE.0 .AND. j.EQ.0) CYCLE                          ! Points (i,0) (with i non-zero) does not exist
          IF(i.EQ.0 .AND. j.EQ.0 .AND. (.NOT.use_nominal)) CYCLE ! Skip the nominal if it is not available
          IF(.NOT.succ_cob(i,j) .OR. chi_cob(i,j).GT.5.d0) CYCLE ! Consider only the successful points with chi <= 5
          IF(elcob(i,j)%h_mag.GT.hmax) CYCLE                     ! Discard the objects with H > H_max
          ! Take into account the VA with negative two-body energy with the Earth
          IF(eearth(i,j).LE.0.d0)THEN
             pearthb=pearthb+EXP(-chi_cob(i,j)**2/2.d0)*jacr(i,j)*jacmu(i,j)/pronorm!*jacarho(i,j)
          ELSE
             WRITE(mulname,590) i,j
590          FORMAT('cob_',I2,'_',I2)
             CALL rmsp(mulname,le)
             WRITE(iuncla,*) mulname
             CALL cov_not_av ! Propagate without covariance
             njc=0           ! Number of minima during a close approach
             CALL pro_ele(elcob(i,j),tmcla,el1) 
             WRITE(*,*) 'cob point ',i,'  ',j,'  propagated'
             WRITE(iun_log,*) 'cob point ',i,'  ',j,'  propagated'
             !****************!
             !  MTP analysis  !
             !****************!
             IF(njc.GT.0)THEN
                !-------------------------------------------------------!
                ! WARNING: we should handle double minima with a DO on  !
                ! jc and selecting the minimum mtp_store(jc)%r          !
                !-------------------------------------------------------!
                IF(njc.GT.1)THEN
                   WRITE(iun_log,*) '   multiple minima, njc = ',njc
                END IF
                IF(mtp_store(1)%iplam.EQ.3)THEN
                   csi(i,j)=mtp_store(1)%tp_coord(1)/reau
                   zeta(i,j)=mtp_store(1)%tp_coord(2)/reau
                   vel(i,j)=mtp_store(1)%tp_coord(3)/gk
                   tcla(i,j)=mtp_store(1)%tcla
                   IF(csi(i,j)**2+zeta(i,j)**2.LT.1.d0)THEN
                      nimp=nimp+1 ! Counter for impacting points
                      !**********************************!
                      !  Compute the Impact Probability  !
                      !**********************************!
                      imppro(i,j)=EXP(-chi_cob(i,j)**2/2.d0)*jacr(i,j)*jacmu(i,j)/pronorm!*jacarho(i,j)
                      iptot=iptot+imppro(i,j)
                      h_magmax=MAX(h_magmax,elcob(i,j)%h_mag)
                      h_magmin=MIN(h_magmin,elcob(i,j)%h_mag)
                      t_minimp=MIN(t_minimp,tcla(i,j));
                      t_maximp=MAX(t_minimp,tcla(i,j));
                   ENDIF
                   WRITE(iunimp,591) i, j, csi(i,j), zeta(i,j), vel(i,j), tcla(i,j), imppro(i,j), 1
591                FORMAT(I4,1X,I4,1X,3(F12.6,1X),1X,F8.2,1X,F10.8,1X,I1)         
                ENDIF
             ENDIF
          ENDIF
       ENDDO
    ENDDO
    CALL filclo(iunimp,' ')

    !***********************!
    !  WRITE THE RISK FILE  !
    !***********************!
    !-------------------------------------------------------------!
    ! We consider significant a case that has more than 3         !
    ! observations, and an arc of more than 0.02 hours (about 28  !
    ! minutes), or a case which has been selected through the     !
    ! select_risk flag, which is true when geoc_chi > 1.          !
    !-------------------------------------------------------------!  
    dt=obs(m)%time_tdt-obs(1)%time_tdt                       ! Arc duration
    significant=(m.GT.2 .AND. dt.GT.0.02) .OR. (select_risk) ! Significant flag
    !------------------------------------------------------------------!
    ! If significant is false we open a .nosig file. Otherwise, if     !
    ! nimp = 0 (there are no impacting points) we open a .norisk file, !
    ! while if nimp > 0 (and thus IP > 0) we open a .risk file.        !
    !------------------------------------------------------------------!
    IF(.NOT.significant)THEN
       file=astnam//'.nosig'
    ELSEIF(nimp.GT.0)THEN
       file=astnam//'.risk'
       !-------------------------------!
       ! lowrisk file not used anymore !
       !-------------------------------!
       !IF(iptot.LE.1.d-2)THEN
       !   file=astnam//'.lowrisk'
       !ENDIF
    ELSE
       file=astnam//'.norisk' 
    ENDIF
    !-----------------------------------------------------------------------!
    ! Data policy flag assignment. We give the following data policy flag:  !
    !                                                                       !
    ! - 0 for IP <= 1E-6                                                    !
    !                                                                       !
    ! - 1 for 1E-6 < IP <= 1E-3                                             !
    !                                                                       !
    ! - 2 for 1E-3 < IP <= 1E-2                                             !
    !                                                                       !
    ! - 3 for IP > 1E-2, and when a nominal does not exists or exists       !
    !   but curvature chi <= 10                                             !
    !                                                                       !
    ! - 4 for IP > 1E-2, when a nominal exists and has curvature chi > 10   !
    !-----------------------------------------------------------------------!
    IF(iptot.LE.1.d-6)THEN
       dataflag=0
    ELSEIF(iptot.LE.1.d-3)THEN
       dataflag=1
    ELSEIF(iptot.LE.1.d-2)THEN
       dataflag=2
    ELSEIF(use_nominal .AND. chi.GT.10)THEN
       dataflag=4
    ELSE
       dataflag=3
    ENDIF
    !**************!
    !  File .risk  !
    !**************!
    CALL rmsp(file,le)
    CALL filopn(iunrisk,file(1:le),'unknown')
    IF(use_nominal)THEN
       !***************!
       !  Cobweb case  !
       !***************!
       ! Header
       WRITE(iunrisk,201) 'NEOCP      n.obs     dt   arc     ',    &
            &     'significance of curvature     rms    ndir*np  nposs  nimp      IP     data policy'
       WRITE(iunrisk,201) 'Name                days  type    ',    &
            &     'geod     accel    overall                                                 flag   ' 
       WRITE(iunrisk,201) '----------------------------------',    &
            &     '---------------------------------------------------------------------------------'   
201    FORMAT(A34,A81) 
       ! Content (IP in exp. format if < 1E-3, otherwise in floating point format)
       IF(iptot.GT.0.d0 .AND. iptot.LT.1.d-3)THEN
          WRITE(iunrisk,555) astnam,m,dt,artyp,geoc_chi,acce_chi,chi,csinor, &
               &      ndir*np,nposs,nimp,iptot,dataflag
       ELSE
          WRITE(iunrisk,556) astnam,m,dt,artyp,geoc_chi,acce_chi,chi,csinor, &
               &      ndir*np,nposs,nimp,iptot,dataflag
       END IF
555    FORMAT(A9,2X,I3,5X,F5.2,2X,I2,4X,F7.2,2X,F7.2,4X,F7.2,3X,F6.3,3X,I6,2X,I5,2X,I5,2X,1P,E10.3,5X,I1) 
556    FORMAT(A9,2X,I3,5X,F5.2,2X,I2,4X,F7.2,2X,F7.2,4X,F7.2,3X,F6.3,3X,I6,2X,I5,2X,I5,5X,F5.3,8X,I1)   
    ELSE
       !********************************!
       !  Grid case (no RMS avaliable)  !
       !********************************!
       ! Header
       WRITE(iunrisk,200) 'NEOCP      n.obs     dt   arc     ',    &
            &     'significance of curvature     rms    ndir*np  nposs  nimp      IP     data policy'
       WRITE(iunrisk,200) 'Name                days  type    ',    &
            &     'geod     accel    overall                                                 flag   ' 
       WRITE(iunrisk,200) '----------------------------------',    &
            &     '---------------------------------------------------------------------------------'   
200    FORMAT(A34,A81) 
       ! Content (IP in exp. format if < 1E-3, otherwise in floating point format)
       IF(iptot.GT.0.d0 .AND. iptot.LT.1.d-3)THEN
          WRITE(iunrisk,557) astnam,m,obs(m)%time_tdt-obs(1)%time_tdt,artyp,geoc_chi,acce_chi,chi,        &
               &      ndir*np,nposs,nimp,iptot,dataflag
       ELSE
          WRITE(iunrisk,558) astnam,m,obs(m)%time_tdt-obs(1)%time_tdt,artyp,geoc_chi,acce_chi,chi,        &
               &      ndir*np,nposs,nimp,iptot,dataflag
       END IF
557    FORMAT(A9,2X,I3,5X,F5.2,2X,I2,4X,F7.2,2X,F7.2,4X,F7.2,3X,1X,'  -  ',3X,I6,2X,I5,2X,I5,2x,1p,E10.3,5X,I1)
558    FORMAT(A9,2X,I3,5X,F5.2,2X,I2,4X,F7.2,2X,F7.2,4X,F7.2,3X,1X,'  -  ',3X,I6,2X,I5,2X,I5,5x,F5.3,8X,I1)
    ENDIF
    !***************!
    !  Write notes  !
    !***************!
    ! Elongation from the Sun
    WRITE(iunrisk,597) elong*degrad
    ! Maximum and minimum magnitude of impactors
    ! Maximum and minimum time of impact
    IF(nimp.GT.0)THEN
       CALL mjddat(t_minimp,iday_min,imonth_min,iyear_min,hour_min)
       CALL mjddat(t_maximp,iday_max,imonth_max,iyear_max,hour_max)
       WRITE(iunrisk,595) h_magmin,h_magmax,iyear_min,imonth_min,iday_min,hour_min,iyear_max,imonth_max,iday_max,hour_max
    ENDIF
    ! Probability to be an Earth-bound object
    WRITE(iunrisk,596) pearthb
    ! Shooting star limit
    WRITE(iunrisk,598) hmax
    ! OrbFit version, date of computation
    CALL date_and_time(date,time)
    WRITE(iunrisk,121) pvers, date, time
    ! Reason for which the case is not significant
    IF(.NOT.significant) THEN
       IF(m.LE.2) THEN
          WRITE(iunrisk,122) m
       ELSEIF(dt.LE.0.02)THEN
          WRITE(iunrisk,123) dt
       ENDIF
    ENDIF
597 FORMAT(/'Elongation from the Sun = ',F5.1)
595 FORMAT(/'Absolute magnitude of impactors between ',F5.2,' and ',F5.2/ &
         &  'Time of impact between ',I4,'/',I2,'/',I2,' at hh ',F6.3,' and ',I4,'/',I2,'/',I2,' at hh ',F6.3)
596 FORMAT(/'Probability to be an Earth-bounded object = ', F5.3)
598 FORMAT(/'Shooting star limit at ',F5.2, ' magnitudes')
121 FORMAT(/'OrbFit software version ',A15,' Date of computation = ',A10,1X,A10,' CET')
122 FORMAT(/'This case is not significant because the number of observations is ',I3)
123 FORMAT(/'This case is not significant because dt is ',F5.2)
    CALL filclo(iunrisk,' ')
  END SUBROUTINE immediate_imp
  
  
  !====================================================================!
  ! SPIDER                                                             !
  !====================================================================!
  ! Computation of the points in the (R,theta) space of polar elliptic !
  ! coordinates                                                        !
  !====================================================================!
  SUBROUTINE spider(el,unc,csinor,npt,r_rdot,theta,r,ecc,q,qg,enne,found,jacthr)
    !=======================================================================================================
    TYPE(orbit_elem), INTENT(IN)            :: el                      ! Nominal elements (from the main program)
    TYPE(orb_uncert), INTENT(IN)            :: unc                     ! Uncertainty matrices
    DOUBLE PRECISION, INTENT(IN)            :: csinor                  ! Residuals norm
    INTEGER,          INTENT(OUT)           :: npt                     ! Total number of points with rho > 0 and e < 1
    DOUBLE PRECISION, INTENT(OUT)           :: r_rdot(0:nmax,0:nmax,2) ! (rho,rho_dot) of the sampling points
    DOUBLE PRECISION, INTENT(OUT)           :: theta(nmax)             ! Array of angles
    DOUBLE PRECISION, INTENT(OUT)           :: r(nmax)                 ! Array of distances (polar coordinates)
    DOUBLE PRECISION, INTENT(OUT)           :: ecc(0:nmax,0:nmax)      ! Eccentricity for each point
    DOUBLE PRECISION, INTENT(OUT)           :: q(0:nmax,0:nmax)        ! Pericenter for each point
    DOUBLE PRECISION, INTENT(OUT)           :: qg(0:nmax,0:nmax)       ! Apocenter for each point
    DOUBLE PRECISION, INTENT(OUT)           :: enne(0:nmax,0:nmax)     ! Mean motion for each point
    LOGICAL,          INTENT(OUT)           :: found(0:nmax,0:nmax)    ! True if e<1, r>0
    DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: jacthr(0:nmax,0:nmax)   ! Jacobian det. for (R,theta)->(rho,rho_dot)
    !=======================================================================================================
    TYPE(orbit_elem) :: eltmp(nmax,nmax) ! Temporary array of elements
    DOUBLE PRECISION :: gamma_rho(2,2)   ! Covariance matrix for (rho,rho_dot)
    DOUBLE PRECISION :: eigenval(2)      ! Eigenvalues of gamma_rho
    DOUBLE PRECISION :: eigenvec(2,2)    ! Matrix of eigenvectors of gamma_rho
    DOUBLE PRECISION :: fv1(2),fv2(2)    ! Temporary storage arrays for rs
    DOUBLE PRECISION :: v1(2),v2(2)      ! Eigenvectors of gamma_rho
    DOUBLE PRECISION :: a,b              ! Ellipse's semiaxis
    INTEGER          :: ierr             ! Error flag
    INTEGER          :: i,j              ! Loop indexes
    !=======================================================================================================
    !******************!
    !  Inizialization  !
    !******************!
    r_rdot=0.d0
    theta=0.d0
    r=0.d0
    ecc=0.d0
    q=0.d0
    qg=0.d0
    enne=0.d0
    found=.FALSE.
    IF(PRESENT(jacthr)) jacthr=0.d0
    !***************************************!
    !  Covariance matrix for (rho,rho_dot)  !
    !***************************************!
    gamma_rho(1,1)=unc%g(5,5)
    gamma_rho(1,2)=unc%g(5,6)
    gamma_rho(2,1)=unc%g(6,5)
    gamma_rho(2,2)=unc%g(6,6)
    ! Eigenvalues and eigenvectors of gamma_rho  
    CALL rs(2,2,gamma_rho,eigenval,1,eigenvec,fv1,fv2,ierr)
    IF(ierr.NE.0)THEN
       WRITE(*,*) 'spider: failed rs, ierr = ',ierr
       STOP
    END IF
    v1=eigenvec(:,1)
    v2=eigenvec(:,2)
    ! Ellipse's semiaxis
    a=SQRT(eigenval(1))
    b=SQRT(eigenval(2))  
    !***************************************************************!
    !  Regular grid in (R,theta) plane: elliptic polar coordinates  !
    !***************************************************************!
    !------------------------------------------------------------------!
    ! Create a regular grid of points in the space of polar elliptic   !
    ! coordinates (R,theta), where R is in [0,sigx] (sigx is a         !
    ! parameter defining the maximum value of RMS we consider reliable !
    ! in our analysis), and theta is in [0,2*pi). sigx is currently    !
    ! 5.d0, and it is set into fitobs.def.                             !
    !------------------------------------------------------------------!
    DO i=1,ndir
       theta(i)=(dpig/ndir)*i  
    END DO
    DO j=1,np
       r(j)=(sigx/np)*j 
    END DO
    ! Counter for the points for which rho > 0 and e < 1
    npt=0
    DO i=1,ndir
       DO j=1,np
          !***************************************!
          !  Map from (R,theta) to (rho,rho_dot)  !
          !***************************************!
          r_rdot(i,j,1)=r(j)*(a*COS(theta(i))*v1(1)-b*SIN(theta(i))*v1(2))+el%coord(5)
          r_rdot(i,j,2)=r(j)*(a*COS(theta(i))*v1(2)+b*SIN(theta(i))*v1(1))+el%coord(6)
          IF(PRESENT(jacthr))THEN
             jacthr(i,j)=r(j)*a*b
          ENDIF
       END DO
    END DO
    !******************************************************!
    !  Compute orbital quantities for each sampling point  !
    !******************************************************!
    ! Check for the coordinates (paranoia check)
    IF(el%coo.NE.'ATT')THEN
       WRITE(*,*) 'spider: elements must be in ATT cordinates, they are in ', el%coo
       STOP
    END IF
    ! Definition of the nominal solution, with indexes (i,j)=(0,0)
    r_rdot(0,0,1:2)=el%coord(5:6)
    DO i=0,ndir
       DO j=0,np
          ! Assign the (i,j)-th point of the sampling to (rho,rho_dot)
          eltmp(i,j)=el
          eltmp(i,j)%coord(5)=r_rdot(i,j,1)
          eltmp(i,j)%coord(6)=r_rdot(i,j,2)
          ! Compute eccentricity, pericenter, apocenter and mean motion
          CALL ecc_peri(eltmp(i,j),ecc(i,j),q(i,j),qg(i,j),enne(i,j))
          IF(r_rdot(i,j,1).GT.0.d0 .AND. ecc(i,j).LT.1.d0)THEN
             npt=npt+1
             found(i,j)=.TRUE.
          ELSE
             found(i,j)=.FALSE.
          ENDIF
       END DO
    END DO
  END SUBROUTINE spider
  

  !====================================================================!
  ! MOV                                                                !
  !====================================================================!
  ! Computation of the Manifold of Variations                          !
  !====================================================================!
  SUBROUTINE mov(m,obs,obsw,el,rms,r_rdot,found,nd,elc,uncert,rms_points,delnor,succ,chi)
    INCLUDE 'parobx.h90' 
    !=======================================================================================================
    INTEGER,          INTENT(IN)            :: m                         ! Number of obervations
    TYPE(ast_obs),    INTENT(IN)            :: obs(m)                    ! Observations
    TYPE(ast_wbsr),   INTENT(INOUT)         :: obsw(m)                   ! Weights and observations
    TYPE(orbit_elem), INTENT(IN)            :: el                        ! Nominal elements
    DOUBLE PRECISION, INTENT(IN)            :: rms                       ! Residuals norm
    DOUBLE PRECISION, INTENT(IN)            :: r_rdot(0:nmax,0:nmax,2)   ! Sampling in the (rho,rho_dot) space
    LOGICAL,          INTENT(IN)            :: found(0:nmax,0:nmax)      ! True if e < 1 and rho > 0 
    INTEGER,          INTENT(IN)            :: nd                        ! Dimension of the parameters spaces
    TYPE(orbit_elem), INTENT(OUT)           :: elc(0:nmax,0:nmax)        ! Corrected elements
    TYPE(orb_uncert), INTENT(OUT)           :: uncert(0:nmax,0:nmax)     ! 6-dim uncertainty matrices
    DOUBLE PRECISION, INTENT(INOUT)         :: rms_points(0:nmax,0:nmax) ! RMS of corrected element
    DOUBLE PRECISION, INTENT(OUT)           :: delnor(0:nmax,0:nmax)     ! Differential corrections norm
    LOGICAL,          INTENT(OUT)           :: succ(0:nmax,0:nmax)       ! Success flag for each point
    DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: chi(0:nmax,0:nmax)        ! Chi function
    !===== Computation of the Jacobian =====================================================================
    TYPE(orb_uncert) :: uncert4(0:nmax,0:nmax) ! Uncertainty matrices of the 4-dim fit
    INTEGER          :: scal_obs               ! Number of scalar observations
    DOUBLE PRECISION :: des_mat(nob2x,nd)      ! Design matrix
    DOUBLE PRECISION :: b_rho(nob2x,2)
    DOUBLE PRECISION :: b_A(nob2x,4)
    DOUBLE PRECISION :: b_A_rho(4,2)
    DOUBLE PRECISION :: dAdrho(4,2)
    DOUBLE PRECISION :: M_A(2,2)
    INTEGER          :: izer
    DOUBLE PRECISION :: aval(1:4)
    DOUBLE PRECISION :: G(4,4), det
    !=======================================================================================================
    TYPE(orbit_elem) :: eltemp(0:nmax,0:nmax)  ! Temporary orbit elements
    INTEGER          :: nused(0:nmax,0:nmax)   ! Number of observations used in the differential corrections 
    DOUBLE PRECISION :: rmsh(0:nmax,0:nmax)    ! Weighted RMS of the residuals 
    INTEGER          :: i,j,ii                    ! Loop indexes
    LOGICAL          :: succ4                  ! Success of fourdim_fit
    DOUBLE PRECISION :: rmsmin                 ! Minimum RMS among the succesful sampling points
    INTEGER          :: verb_dif_old           ! Temporary verbosity level  
    !=======================================================================================================
    !******************!
    !  Inizialization  !
    !******************!
    succ = .FALSE.
    ! Anti-verbose
    verb_dif_old = verb_dif
    verb_dif     = 1
    !***************************************************************!
    !  Assign the nominal elements to the point with indexes (0,0)  !
    !***************************************************************!
    IF(rms.GE.0.d0)THEN
       elc(0,0)        = el
       rms_points(0,0) = rms
       chi(0,0)        = 0.d0
       succ(0,0)       = .TRUE.    
    ENDIF
    !******************************************************!
    !  4-dim differential correction along each direction  !
    !******************************************************!
    ! Initialization to determine the minimum value of the RMS
    rmsmin=100.d0
    jacmu=0.d0
    DO i=1,ndir
       eltemp(i,0)=el ! For j=0, give the nominal elements
       DO j=1,np
          IF(found(i,j))THEN
             IF(rms.GE.0.d0)THEN
                eltemp(i,j) = eltemp(i,j-1) ! Give the corrected elements of the previous step
             ELSE
                eltemp(i,j) = el            ! Give the nominal elements
             ENDIF
             eltemp(i,j)%coord(5) = r_rdot(i,j,1)
             eltemp(i,j)%coord(6) = r_rdot(i,j,2)
             CALL fourdim_fit(m,obs,obsw,eltemp(i,j),elc(i,j),nd,uncert4(i,j),uncert(i,j),&
                  & rms_points(i,j),delnor(i,j),rmsh(i,j),nused(i,j),succ4,des_mat,scal_obs)
             succ(i,j) = uncert(i,j)%succ.AND.succ4
             !***************************************
             !WRITE(111,*) 'i=',i,' j=',j
             !WRITE(111,*) 'DESIGN MATRIX after fourdim_fit='
             !DO ii=1,scal_obs
             !   WRITE(111,*) des_mat(ii,1:6)
             !END DO
             !*************************************
             !***********************************************!
             !  SECOND ORDER COMPUTATION OF THE DETERMINANT  !
             !***********************************************!
             IF(succ(i,j))THEN
                b_rho(1:scal_obs,:)=des_mat(1:scal_obs,5:6)
                b_A(1:scal_obs,:)=des_mat(1:scal_obs,1:4)
                
                !WRITE(111,*) 'C_A after fourdim_fit='
                !DO ii=1,4
                !   WRITE(111,*) uncert4(i,j)%c(ii,1:4)
                !END DO
                CALL qr_inv(uncert4(i,j)%c(1:4,1:4),G,4,izer,det,aval(1:4))
                !WRITE(*,*) 'G_A after fourdim_fit='
                !DO ii=1,4
                !   WRITE(*,*) uncert4(i,j)%g(ii,1:4)
                !END DO
                   
                b_A_rho=MATMUL(TRANSPOSE(b_A(1:scal_obs,:)),b_rho(1:scal_obs,:))
             !   !WRITE(111,*) 'B_A_rho ='
             !   !DO ii=1,4
             !   !   WRITE(111,*) G(ii,1:2)
             !   !END DO
             !   
                dAdrho=MATMUL(uncert4(i,j)%g(1:4,1:4),b_A_rho)
             !   !WRITE(111,*) 'dAdrho ='
             !   !DO ii=1,4
             !   !   WRITE(111,*) dAdrho(ii,1:2)
             !   !END DO
             !   
                M_A=MATMUL(TRANSPOSE(dAdrho),dAdrho)
             !   !WRITE(111,*) 'M_A ='
             !   !DO ii=1,2
             !   !   WRITE(111,*) M_A(ii,1:2)
             !   !END DO
                jacmu(i,j) = SQRT((M_A(1,1)+1.d0)*(M_A(2,2)+1.d0)-M_A(1,2)*M_A(2,1))
             !   WRITE(111,*) 'I= ',i,' J= ',j,'     DET_MU=', jacmu(i,j)
             !   !WRITE(111,*) '********************************************************'
                rmsmin=MIN(rms_points(i,j),rmsmin)
             ENDIF
          ELSE
             succ(i,j)=.FALSE.
          ENDIF
       END DO
    END DO
    ! Compare the minimum RMS with the RMS of the nominal, if given
    IF(rms.LT.0.d0)THEN
       ! no nominal, use minimum
       rms_points(0,0)=rmsmin  
    ELSE
       IF(rmsmin.LT.rms)THEN
          WRITE(*,*) 'mov: found RMS =',rmsmin,', lower than nominal RMS = ', rms
       ENDIF
    ENDIF
    !***********************************!
    !  Computation of the chi function  !
    !***********************************!
    DO i=1,ndir
       DO j=1,np
          chi(i,j)=SQRT(ABS(rms_points(i,j)**2-rms_points(0,0)**2)*nused(i,j))
       END DO
    END DO
    ! Restore the beginning verbosity level
    verb_dif=verb_dif_old
  END SUBROUTINE mov

  !**********************
  ! COMPUTE_SCORE
  !**********************
  ! Compute the score of the object to be a main belt or a NEA
  SUBROUTINE compute_score(score_distant, score_mba, score_nea)
    DOUBLE PRECISION, INTENT(OUT) :: score_distant, score_mba, score_nea  ! score of the object (Distant, Main Belt and NEA)
    !*******************************************************************************************
    DOUBLE PRECISION :: nposs, n_nea, n_mba, n_dis  ! number of possible points
    INTEGER          :: i, j                        ! loop indices
    TYPE(orbit_elem) :: elcom                       ! cometarian orbital elements
    INTEGER          :: fail_flag                   ! fail flag in changing coordinates
    DOUBLE PRECISION :: qmin, qgmax                 ! Limit values for perihelion and aphelion
    DOUBLE PRECISION :: q, qg                       ! Perihelion and aphelion
    !*******************************************************************************************
    nposs=0
    n_nea=0
    n_mba=0
    n_dis=0
    qmin=1.4d0
    qgmax=7.d0
    ! loop on cobweb points
    DO i=0,ndir
       DO j=0,np
          ! this is a nasty method to write separately the nominal, if it exists
          IF(i.eq.0.and.j.ne.0) CYCLE
          IF(i.ne.0.and.j.eq.0) CYCLE
          IF(i.eq.0.and.j.eq.0.and.(.not.use_nominal)) CYCLE
          IF(.not.succ_cob(i,j).or.chi_cob(i,j).gt.5.d0)THEN
             CYCLE
          ENDIF
          IF(succ_cob(i,j).AND.elcob(i,j)%h_mag.le.hmax)THEN
             nposs=nposs+EXP(-chi_cob(i,j)**2/2)
          ENDIF
          CALL coo_cha(elcob(i,j),'COM',elcom,fail_flag)
          ! Perihelion and aphelion
          q  = elcom%coord(1)
          qg = q*(1+elcom%coord(2))/(1-elcom%coord(2))
          !NEA or MBA
          IF(chi_cob(i,j).LT.5.d0)THEN
             IF(qg.LT.qgmax)THEN
                !NEA
                IF(q.LT.qmin)THEN
                   n_nea=n_nea+EXP(-chi_cob(i,j)**2/2)
                !MBA
                ELSE
                   n_mba=n_mba+EXP(-chi_cob(i,j)**2/2)
                END IF
             ELSE
                ! Distant
                n_dis=n_dis+EXP(-chi_cob(i,j)**2/2)
             END IF
          END IF
       ENDDO
    ENDDO
    score_nea     = n_nea/nposs
    score_mba     = n_mba/nposs
    score_distant = n_dis/nposs
  END SUBROUTINE compute_score
  
  !====================================================================!                  
  ! MC_RES_ATT                                                         !
  !====================================================================!
  ! Search for immediate impactors with MC method on residuals         !
  ! assuming a nominal solution is available                           !
  !====================================================================!
  SUBROUTINE mc_res_att(astnac,nmc,tmcla,m,obs,obsw,el0,obscodn,nd,nimp,isucc)
    USE tp_trace,     ONLY: tpplane, mtp_store, njc
    USE close_app,    ONLY: fix_mole, kill_propag
    USE propag_state, ONLY: pro_ele
    INCLUDE 'parobx.h90'
    !=======================================================================================================
    CHARACTER*(name_len), INTENT(IN)  :: astnac     ! Asteroid name
    INTEGER,              INTENT(IN)  :: nmc        ! Number of MC sample points
    DOUBLE PRECISION,     INTENT(IN)  :: tmcla      ! Final time of propagation
    INTEGER,              INTENT(IN)  :: obscodn    ! Observatory code
    INTEGER,              INTENT(IN)  :: m          ! Number of observations
    TYPE(ast_obs),        INTENT(IN)  :: obs(nobx)  ! Observations
    TYPE(ast_wbsr),       INTENT(IN)  :: obsw(nobx) ! Residuals and weigths
    TYPE(orbit_elem),     INTENT(IN)  :: el0        ! Nominal elements
    INTEGER,              INTENT(IN)  :: nd         ! Dimension of the par space
    INTEGER,              INTENT(OUT) :: nimp       ! Number of impacts
    INTEGER,              INTENT(OUT) :: isucc      ! Number of convergent differential corrections
    !=======================================================================================================
    CHARACTER(LEN=30) :: file        ! File name
    !===== MC computation ==================================================================================
    INTEGER           :: k           ! Index for the MC points
    INTEGER           :: kk          ! For one istance of noise
    DOUBLE PRECISION  :: rvnorm      ! Random noise
    TYPE(ast_obs)     :: obs1(nobx)  ! Perturbed observations
    TYPE(ast_wbsr)    :: obsw1(nobx) ! Perturbed weights
    INTEGER           :: icor(6)     ! Parameters to be corrected
    TYPE(orbit_elem)  :: elk         ! Keplerian elements (after diff_cor)
    TYPE(orb_uncert)  :: unck        ! Uncertainty matrices for elk
    TYPE(orbit_elem)  :: elatt       ! Attributable orbital elements
    TYPE(orbit_elem)  :: el1         ! Propagated elements
    LOGICAL           :: succ        ! Successful differential corrections
    DOUBLE PRECISION  :: csinok      ! Norm of the residuals
    DOUBLE PRECISION  :: delnok      ! Norm of the last correction
    DOUBLE PRECISION  :: csi         ! Xi on the TP
    DOUBLE PRECISION  :: zeta        ! Zeta on the TP
    !===== Other variables =================================================================================
    INTEGER           :: fail_flag   ! Failure flag for the coordinate change
    INTEGER           :: le          ! File length
    INTEGER           :: iuncob      ! Output file for the MC points
    !=======================================================================================================
    IF(nd.NE.6)THEN
       WRITE(*,*) 'mc_res_att: not ready for non-grav, nd = ',nd
       STOP
    END IF
    !******************!
    !  Inizialization  !
    !******************!
    verb_dif=1
    tpplane=.FALSE. ! Use MTP plane, not TP
    !***************************!
    !  File .mc for the output  !
    !***************************!
    file=astnac//'.mc'
    CALL rmsp(file,le)
    CALL filopn(iuncob,file(1:le),'unknown')
    
    !***************************************!
    !  LOOP ON INSTANCES OF GAUSSIAN NOISE  !
    !***************************************!
    nimp=0  ! Impact counter
    isucc=0 ! Counter of nominal solutions found
    DO k=1,nmc
       !********************************!
       !  Create one instance of noise  !
       !********************************!
       DO kk=1,m 
          IF(obs(kk)%type.EQ.'O')THEN
             obs1(kk)=obs(kk)
             obs1(kk)%coord(1)=obs1(kk)%coord(1)+rvnorm()*obsw(kk)%rms_coord(1)/COS(obs1(kk)%coord(2))
             obs1(kk)%coord(2)=obs1(kk)%coord(2)+rvnorm()*obsw(kk)%rms_coord(2)
             obsw1(kk)=obsw(kk)
          ENDIF
       ENDDO
       !****************************!
       !  Compute nominal solution  !
       !****************************!
       icor=1    
       CALL diff_cor(m,obs1,obsw1,el0,icor,-1,elk,unck,csinok,delnok,succ,nd)
       ! Check convergence
       IF(succ)THEN
          isucc=isucc+1
          ! First output the data point on the (rho,rho_dot) plane
          CALL coo_cha(elk,'ATT',elatt,fail_flag,OBSCODE=obscodn)
          WRITE(iuncob,130) k,isucc,elatt%coord(5:6),0.d0,0.d0,csinok,0.d0,0.d0,0.d0,&
               & 1,delnok,elatt%h_mag,0.d0
130       FORMAT(I4,1X,I4,1X,4(F12.7,1X),1P,D12.5,1X,D12.5,1X,D12.5,1X,D12.5,0P,1X,I1,1X,1P,D10.3,0P,1X,F6.3,1X,1P,D12.5)
          !******************************!
          !  Propagate until final time  !
          !******************************!
          CALL cov_not_av
          njc=0
          CALL pro_ele(elk,tmcla,el1)
          kill_propag=.FALSE.
          !******************************************!
          !  Read close approach data, if available  !
          !******************************************!
          IF(njc.GT.0)THEN
             IF(mtp_store(1)%iplam.EQ.3)THEN
                csi=mtp_store(1)%tp_coord(1)/reau
                zeta=mtp_store(1)%tp_coord(2)/reau
                ! Impactor?
                IF(csi**2+zeta**2.LT.1.d0)THEN
                   nimp=nimp+1
                ENDIF
             ENDIF
          ENDIF
       ENDIF
       WRITE(*,*) k, isucc, nimp, csinok
    ENDDO
    CALL filclo(iuncob,' ')
    verb_dif=20
  END SUBROUTINE mc_res_att

END MODULE cobweb
