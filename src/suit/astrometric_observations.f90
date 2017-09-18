! MODULE astrometric_observations
! Handling of astrometric observations

! TYPE ast_obs
!     astrometric observations (optical or radar)
! TYPE ast_wbsr
!     weighting, bias, selection flag and residuals of astrometric observation
! CONTAINS

! SUBROUTINE write_rwo(file,obs,obsw,n,error_model,rms,rmsmag,sort,unit0)
! output observations in .rwo format
! SUBROUTINE write_rwo_rad(unit, obs, obsw)
! SUBROUTINE write_rwo_opt(unit,obs,obsw,error,nrec)
! SUBROUTINE write_rwo_pos(unit,obs,obsw)
! SUBROUTINE write_rwo_header(unit, error_model, nrec, rms, rmsmag)

! SUBROUTINE read_rwo(file,obs,obsw,n,error_model,rms,rmsmag,unit0,eof0)
! input observations in .rwo format
! SUBROUTINE read_rwo_rec(unit,nr,obs1,obsw1,error, eof)
! SUBROUTINE read_rwo_header(unit, error_model, error, rms, rmsmag)

! SUBROUTINE mpc_obs_input(obs_found,obs,nobs,unit,filnam,eof,ons_mode)
! SUBROUTINE mpcrec_transform(mpcrec,obs,complete,error,skip_rad,ons_mode)
!     transformation of a MPC record into an observation (TYPE ast_obs)
! SUBROUTINE mpcrec_add_obspos(obs,error,mpcrec1,mpcrec2)

! SUBROUTINE jpl_radarobs(obs_found,file,obs,nobs)
! SUBROUTINE jplradar_transform(rec,obs,error)
! CHARACTER(LEN=11) FUNCTION radar_resid_string(value,defined)
! INTEGER FUNCTION radar_station(stastr)

! SUBROUTINE input_obs(obsdir,astna0,precob,error_model,obs0,obs,obsw,m,iun20,change,rms,rmsmag)
!     input of observations: reads sequentially the .rwo, .obs, .rad file and combines
!     the data according to the precedence rule specified by precob

! SUBROUTINE addobs_mpc(obs,obsw,m,obst,obswt,mt,mnew,change)
! SUBROUTINE addobs_rwo(obs,obsw,m,obst,obswt,mt,change)
! INTEGER FUNCTION find_obs(obst,obs,m,double)
! LOGICAL FUNCTION changed_obs(obst,obswt,obs,obsw)

! SUBROUTINE aster_radius(objid,obstype,m)

! SUBROUTINE rdanga(string,angle,acc,error)

! HEADERS AND MODULES
! astrometric_observations.o: \
!	../include/parcmc.h90 \
!	../include/parobx.h90 \
!	../include/sysdep.h90 \
!	fund_const.o \
!	output_control.o \
!	reference_systems.o \
!	station_coordinates.o 


MODULE astrometric_observations
USE fund_const
USE output_control
USE name_rules
IMPLICIT NONE
PRIVATE
! SUBROUTINEs
PUBLIC :: input_obs,mpc_obs_input, addobs_rwo, addobs_mpc, aster_radius,jpl_radarobs,find_obs
PUBLIC :: read_rwo,write_rwo,read_rwo_rec,read_rwo_header, write_rwo_opt,write_rwo_header
! OPERATOR OVERLOADING
PUBLIC drop_geopos, add_geopos
! CHANGE OF RWO FORMAT

! common data for sub-module observ_rms_bias
! to read cbm10.sta
   INTEGER, PARAMETER :: nstarmsx=100
   DOUBLE PRECISION, SAVE  :: rmsstabin(2,nstarmsx) ! rms(a*cosd), rmsd, in arcsec
   CHARACTER*3, SAVE :: obscodrms(nstarmsx) ! star catalog code (' '=not known)
   INTEGER, SAVE :: indstarms(nstarmsx) ! lexicographic sort index
   INTEGER, SAVE :: nstarms ! actual number
! to read cbm10.rules, fcct14.rules
   INTEGER, PARAMETER :: nrulesx=100
   DOUBLE PRECISION, SAVE  :: rmsbin(2,nrulesx) ! rms(a*cosd), rmsd, in arcsec
   CHARACTER*4, SAVE :: obscatcod(nrulesx) ! star catalog code (' '=not known)
   INTEGER, SAVE :: indrules(nrulesx) ! lexicographic sort index
   INTEGER, SAVE :: nrules ! actual number
! public routines
PUBLIC observ_rms

CHARACTER*(6), PUBLIC :: change_format ! INPUT,OUTPUT, ' ' 

LOGICAL, PUBLIC :: ons_name  ! if true, the identifier in an .obs file is not a designation

!INTEGER, PARAMETER :: len_objdes = 9   ! length of the object designation
!PUBLIC len_objdes

! Generic astrometric observation
TYPE ast_obs
   CHARACTER(LEN=name_len)         :: objdes     ! object designation
   CHARACTER(LEN=1)                :: type       ! observation type:
                                                 !    O = optical
                                                 !    R = radar range
                                                 !    V = radar range rate
                                                 !    S = satellite/space
   CHARACTER(LEN=1)                :: tech       ! observation technology (column 15 of MPC record):
                                                 !    P = Photographic
                                                 !    E = Encoder
                                                 !    C = CCD
                                                 !    T = Transit or meridian circle
                                                 !    M = Micrometer
                                                 !    c = corrected for center of mass (radar only)
                                                 !    s = surface bounce (radar only)
   CHARACTER(LEN=1)                :: par_unit   ! parallax unit for satellite observations:
                                                 !    1 = km
                                                 !    2 = au
   CHARACTER(LEN=1)                :: note       ! note on observation circumstances (column 14 of MPC record)
   CHARACTER(LEN=2)                :: catcodmpc  ! MPC catalog flag:
                                                 ! a = USNO-A1.0
                                                 ! b = USNO-SA1.0
                                                 ! c = USNO-A2.0
                                                 ! d = USNO-SA2.0
                                                 ! e = UCAC-1
                                                 ! f = Tycho-1
                                                 ! g = Tycho-2
                                                 ! h = GSC-1.0
                                                 ! i = GSC-1.1
                                                 ! j = GSC-1.2
                                                 ! k = GSC-2.2
                                                 ! l = ACT
                                                 ! m = GSC-ACT
                                                 ! n = TRC
                                                 ! o = USNO-B1.0
                                                 ! p = PPM
                                                 ! q = UCAC-2-beta
                                                 ! r = UCAC-2
                                                 ! s = USNO-B2.0
                                                 ! t = UCAC-3-beta
                                                 ! u = UCAC-3
                                                 ! v = NOMAD
                                                 ! w = CMC
                                                 ! x = Hip-2
                                                 ! z = GSC (unspecified version)
                                                 ! ' ' = unknown
 
   CHARACTER(LEN=7)                :: obscod_s   ! observer's MPC code (string, official); for radar, TX(A3)//' '//RX(A3)
   INTEGER                         :: obscod_i   ! observer's MPC code (integer, internal representation); for radar,
                                                 !  encoding of both transmitter ans receiver (see jplradar_transform)
   DOUBLE PRECISION                :: time_utc   ! observation UTC time (MJD)
   DOUBLE PRECISION                :: time_tdt   ! observation TDT time (MJD)
   DOUBLE PRECISION, DIMENSION(2)  :: coord      ! observation coordinates:
                                                 !    RA,DEC (optical & satellite) (rad)
                                                 !    range,range rate (radar)
   DOUBLE PRECISION                :: mag        ! object apparent magnitude
   CHARACTER(LEN=6)                :: mag_str    ! magnitude original string in the record
   CHARACTER(LEN=1)                :: mag_band   ! magnitude color band
   LOGICAL                         :: mag_def    ! TRUE if magnitude is defined
   DOUBLE PRECISION                :: acc_time   ! accuracy of time (as given in input format) (d)
   DOUBLE PRECISION, DIMENSION(2)  :: acc_coord  ! accuracy of coordinates (as given in input format)
   DOUBLE PRECISION                :: acc_mag    ! accuracy of magnitude (as given in input format)
! Observer's position and velocity are given with respect to the mean ecliptic and equinox J2000
! WARNING: in case of radar observations:
!         obspos = body-fixed position of transmitting station
!         obsvel = body-fixed position of receiving station
   DOUBLE PRECISION, DIMENSION(3)  :: obspos     ! observer's geocentric position (au)
   DOUBLE PRECISION, DIMENSION(3)  :: obsvel     ! observer's geocentric velocity (au/d)
! The following array is used only for roving observatory observations (obs%type='O' and obs%obscod_s='247')
! and is ALLOCATEd and initialized by SUBROUTINE mpcrec_add_obspos or SUBROUTINE read_rwo
!     geopos(1) = E longitude of observing site (rad)
!     geopos(2) = N latitude of observing site (rad)
!     geopos(3) = Altitude of observing site (m)
   DOUBLE PRECISION, DIMENSION(:), POINTER :: geopos => NULL()
END TYPE ast_obs
! Generic astrometric observation
TYPE ast_obs_direct
   CHARACTER(LEN=name_len)       :: objdes     ! object designation
   CHARACTER(LEN=1)                :: type       ! observation type:
                                                 !    O = optical
                                                 !    R = radar range
                                                 !    V = radar range rate
                                                 !    S = satellite/space
   CHARACTER(LEN=1)                :: tech       ! observation technology (column 15 of MPC record):
                                                 !    P = Photographic
                                                 !    E = Encoder
                                                 !    C = CCD
                                                 !    T = Transit or meridian circle
                                                 !    M = Micrometer
                                                 !    c = corrected for center of mass (radar only)
                                                 !    s = surface bounce (radar only)
   CHARACTER(LEN=1)                :: par_unit   ! parallax unit for satellite observations:
                                                 !    1 = km
                                                 !    2 = au
   CHARACTER(LEN=1)                :: note       ! note on observation circumstances (column 14 of MPC record)
   CHARACTER(LEN=2)                :: catcodmpc  ! MPC catalog flag:
                                                 ! a = USNO-A1.0
                                                 ! b = USNO-SA1.0
                                                 ! c = USNO-A2.0
                                                 ! d = USNO-SA2.0
                                                 ! e = UCAC-1
                                                 ! f = Tycho-1
                                                 ! g = Tycho-2
                                                 ! h = GSC-1.0
                                                 ! i = GSC-1.1
                                                 ! j = GSC-1.2
                                                 ! k = GSC-2.2
                                                 ! l = ACT
                                                 ! m = GSC-ACT
                                                 ! n = TRC
                                                 ! o = USNO-B1.0
                                                 ! p = PPM
                                                 ! q = UCAC-2-beta
                                                 ! r = UCAC-2
                                                 ! s = USNO-B2.0
                                                 ! t = UCAC-3-beta
                                                 ! u = UCAC-3
                                                 ! v = NOMAD
                                                 ! w = CMC
                                                 ! x = Hip-2
                                                 ! z = GSC (unspecified version)
                                                 ! ' ' = unknown
   CHARACTER(LEN=7)                :: obscod_s   ! observer's MPC code (string, official); for radar, TX(A3)//' '//RX(A3)
   INTEGER                         :: obscod_i   ! observer's MPC code (integer, internal representation); for radar,
                                                 !  encoding of both transmitter ans receiver (see jplradar_transform)
   DOUBLE PRECISION                :: time_utc   ! observation UTC time (MJD)
   DOUBLE PRECISION                :: time_tdt   ! observation TDT time (MJD)
   DOUBLE PRECISION, DIMENSION(2)  :: coord      ! observation coordinates:
                                                 !    RA,DEC (optical & satellite) (rad)
                                                 !    range,range rate (radar)
   DOUBLE PRECISION                :: mag        ! object apparent magnitude
   CHARACTER(LEN=6)                :: mag_str    ! magnitude original string in the record
   CHARACTER(LEN=1)                :: mag_band   ! magnitude color band
   LOGICAL                         :: mag_def    ! TRUE if magnitude is defined
   DOUBLE PRECISION                :: acc_time   ! accuracy of time (as given in input format) (d)
   DOUBLE PRECISION, DIMENSION(2)  :: acc_coord  ! accuracy of coordinates (as given in input format)
   DOUBLE PRECISION                :: acc_mag    ! accuracy of magnitude (as given in input format)
! Observer's position and velocity are given with respect to the mean ecliptic and equinox J2000
! WARNING: in case of radar observations:
!         obspos = body-fixed position of transmitting station
!         obsvel = body-fixed position of receiving station
   DOUBLE PRECISION, DIMENSION(3)  :: obspos     ! observer's geocentric position (au)
   DOUBLE PRECISION, DIMENSION(3)  :: obsvel     ! observer's geocentric velocity (au/d)
END TYPE ast_obs_direct

!INTERFACE OPERATOR(=)
!  MODULE PROCEDURE drop_geopos
!END INTERFACE

DOUBLE PRECISION, DIMENSION(3), PARAMETER :: zero_3d_vect = (/ 0.d0, 0.d0, 0.d0 /)

! Undefined astrometric observation (useful for initializing variables before assigning values)
TYPE(ast_obs), PARAMETER :: undefined_ast_obs = AST_OBS( &
&   ' ' ,                &  ! object designation
&   '?' ,                &  ! observation type
&   '?' ,                &  ! observation technology
&   ' ' ,                &  ! parallax unit (only for satellite observations)
&   ' ' ,                &  ! note on observation circumstances (column 14 of MPC record)
&   '  ',                &  ! astrometric catalog code
&   '       ' ,          &  ! observer's MPC code (official)
&   -999999 ,            &  ! observer's MPC code (internal)
&   1.d99 ,              &  ! observation UTC time
&   1.d99 ,              &  ! observation TDT time (MJD)
&   (/ 0.d0, 0.d0 /) ,   &  ! observation coordinates
&   1.d99 ,              &  ! object apparent magnitude
&   '      ' ,           &  ! magnitude original string in the record
&   '?' ,                &  ! magnitude color band
&   .false. ,            &  ! magnitude is not defined
&   1.d99 ,              &  ! accuracy of time (as given in input format) (d)
&   (/ 1.d99, 1.d99 /) , &  ! accuracy of coordinates (as given in input format)
&   1.d99 ,              &  ! accuracy of magnitude (as given in input format)
&   zero_3d_vect ,       &  ! observer's geocentric position (au)
&   zero_3d_vect ,       &  ! observer's geocentric velocity (au/d)
&   NULL()               )  ! observer's geodetic coordinates (longitude, latitude, altitude)

! Weighting, bias, selection flag and residuals of astrometric observation
TYPE ast_wbsr
   INTEGER                        :: sel_coord  ! selection flag for coordinates:
                                                !   0 = not used
                                                !   1 = used in LS fit
                                                !   2 = used in initial orbit determination + LS fit
   DOUBLE PRECISION, DIMENSION(2) :: rms_coord  ! RMS of coordinates to be used in weighting
   DOUBLE PRECISION, DIMENSION(2) :: bias_coord ! biases of coordinates
   LOGICAL,          DIMENSION(2) :: force_w    ! TRUE if weight of coordinates has been forced
   DOUBLE PRECISION, DIMENSION(2) :: res_coord  ! residuals of coordinates
   LOGICAL                        :: resc_def   ! TRUE if residuals of coordinates are defined
   DOUBLE PRECISION               :: chi        ! 2D chi (SQRT(chisquare)) of residuals
   DOUBLE PRECISION               :: rms_mag    ! RMS of magnitude to be used in weighting
   INTEGER                        :: sel_mag    ! selection flag for magnitude:
                                                !   0 = not used
                                                !   1 = used in LS fit
   DOUBLE PRECISION               :: res_mag    ! residual of magnitude
   LOGICAL                        :: resm_def   ! TRUE if residual of magnitude is defined
END TYPE ast_wbsr

! Undefined weight, bias, selection flag and residuals of astrometric observation
! (useful for initializing variables before assigning values)
TYPE(ast_wbsr), PARAMETER :: undefined_ast_wbsr = AST_WBSR( &
   1 ,                      &  ! selection flag for coordinates
   (/ 9.d9, 9.d9 /) ,       &  ! RMS of coordinates to be used in weighting
   (/ 0.d0, 0.d0 /) ,       &  ! biases of coordinates
   (/ .false., .false. /) , &  ! TRUE if weight of coordinates has been forced
   (/ 0.d0, 0.d0 /) ,       &  ! residuals of coordinates
   .false. ,                &  ! TRUE if residuals of coordinates are defined
   0.d0 ,                   &  ! 2D chi (SQRT(chisquare)) of residuals
   9.d9 ,                   &  ! RMS of magnitude to be used in weighting
   0 ,                      &  ! selection flag for magnitude
   0.d0 ,                   &  ! residual of magnitude
   .false.                  )  ! TRUE if residual of magnitude is define

INTEGER, PARAMETER :: rwo_version = 2           !  current version of RWO file format

INTEGER, PARAMETER, PUBLIC :: nradobsx=10000 ! max total number of radar observations
CHARACTER*100, PUBLIC, DIMENSION(nradobsx) :: radar_rec  ! JPL formatted records 

! LIST OF PUBLIC ENTITIES
! Derived TYPEs
PUBLIC :: ast_obs, ast_wbsr, ast_obs_direct
! PARAMETERs
PUBLIC :: undefined_ast_obs, undefined_ast_wbsr
! DATA: only asteroid radius (it is considered part of the observation data)
DOUBLE PRECISION radius
PUBLIC :: radius

CONTAINS

TYPE(ast_obs_direct) FUNCTION drop_geopos(obs)
TYPE(ast_obs), INTENT(IN) :: obs
  drop_geopos%objdes=obs%objdes
  drop_geopos%type=obs%type
  drop_geopos%tech=obs%tech
  drop_geopos%par_unit=obs%par_unit
  drop_geopos%note=obs%note
  drop_geopos%obscod_s=obs%obscod_s
  drop_geopos%obscod_i=obs%obscod_i
  drop_geopos%time_utc=obs%time_utc
  drop_geopos%time_tdt=obs%time_tdt
  drop_geopos%coord=obs%coord
  drop_geopos%mag=obs%mag
  drop_geopos%mag_str=obs%mag_str
  drop_geopos%mag_band=obs%mag_band
  drop_geopos%mag_def=obs%mag_def
  drop_geopos%acc_time=obs%acc_time
  drop_geopos%acc_coord=obs%acc_coord
  drop_geopos%acc_mag=obs%acc_mag
  drop_geopos%obspos=obs%obspos
  drop_geopos%obsvel=obs%obsvel
END FUNCTION drop_geopos
TYPE(ast_obs) FUNCTION add_geopos(obs)
TYPE(ast_obs_direct), INTENT(IN) :: obs
  add_geopos%objdes=obs%objdes
  add_geopos%type=obs%type
  add_geopos%tech=obs%tech
  add_geopos%par_unit=obs%par_unit
  add_geopos%note=obs%note
  add_geopos%obscod_s=obs%obscod_s
  add_geopos%obscod_i=obs%obscod_i
  add_geopos%time_utc=obs%time_utc
  add_geopos%time_tdt=obs%time_tdt
  add_geopos%coord=obs%coord
  add_geopos%mag=obs%mag
  add_geopos%mag_str=obs%mag_str
  add_geopos%mag_band=obs%mag_band
  add_geopos%mag_def=obs%mag_def
  add_geopos%acc_time=obs%acc_time
  add_geopos%acc_coord=obs%acc_coord
  add_geopos%acc_mag=obs%acc_mag
  add_geopos%obspos=obs%obspos
  add_geopos%obsvel=obs%obsvel
  add_geopos%geopos => NULL()
END FUNCTION add_geopos


! SUBMODULE iorwo
! ==========================================
! Write a .rwo file- structured
! allows to select sorting or not sorting by observation time
! however, radar observations are always after optical
! ==========================================
SUBROUTINE write_rwo(file,obs,obsw,n,error_model,rms,rmsmag,sort,unit0)

CHARACTER(LEN=*),             INTENT(IN)           :: file        ! output file
INTEGER,                      INTENT(IN)           :: n           ! number of observations to write
TYPE(ast_obs),  DIMENSION(n), INTENT(IN)           :: obs         ! set of n observations
TYPE(ast_wbsr), DIMENSION(n), INTENT(IN)           :: obsw        ! set of n weights and residuals
CHARACTER(LEN=*),             INTENT(IN)           :: error_model ! error model file name 
DOUBLE PRECISION,             INTENT(IN), OPTIONAL :: rms         ! RMS of astrometric fit
DOUBLE PRECISION,             INTENT(IN), OPTIONAL :: rmsmag      ! RMS of photometric fit
LOGICAL,                      INTENT(IN), OPTIONAL :: sort        ! sort times?
                                ! note default is to sort
INTEGER,                      INTENT(IN), OPTIONAL :: unit0 ! output unit if file is already opened
INCLUDE 'parcmc.h90'     ! comment character
INTEGER, DIMENSION(n) :: time_ord
INTEGER :: unit,i,k,lf,nrec
DOUBLE PRECISION :: rms0, rmsmag0
LOGICAL :: error,radar
INTEGER, EXTERNAL :: lench
EXTERNAL :: filopn,filclo,mjddat,real2string,sessag
DATA nrec /0/
! ===============================================
! check operating mode
IF(PRESENT(unit0))THEN
   unit=unit0
ELSE! open output file
   CALL filopn(unit,file,'UNKNOWN')
   nrec=0
   IF(PRESENT(rms)) THEN
      rms0=rms
   ELSE
      rms0=-1.
   ENDIF
   IF(PRESENT(rmsmag)) THEN
      rmsmag0=rmsmag
   ELSE
      rmsmag0=-1.
   ENDIF
! write header
   CALL write_rwo_header(unit, error_model, nrec, rms0, rmsmag0)
ENDIF
IF(PRESENT(sort))THEN
   IF(sort)CALL heapsort(obs(1:n)%time_tdt,n,time_ord)
ELSE
   CALL heapsort(obs(1:n)%time_tdt,n,time_ord) !default is to sort
ENDIF
! First pass: only optical observations (including roving observatory and satellite) are processed
radar=.false.
DO k=1,n
   IF(PRESENT(sort))THEN
      IF(sort)THEN
         i=time_ord(k)
      ELSE
         i=k
      ENDIF
   ELSE
      i=time_ord(k)
   ENDIF
   SELECT CASE (obs(i)%type)
      CASE ('O', 'S')                  ! OPTICAL/SATELLITE OBSERVATION
         CALL write_rwo_opt(unit,obs(i),obsw(i),error,nrec)
         IF(error)THEN
            lf=lench(file)
            WRITE(ierrou,121) file(1:lf)
121         FORMAT('write_rwo: error in DAY field: file ',A)
! numerr has already been incremeted inside write_rwo_opt
         ENDIF
! Second record for satellite observations
         IF(obs(i)%type == 'S'.or.obs(i)%obscod_s == '247') THEN
             CALL write_rwo_pos(unit,obs(i),obsw(i))
             nrec=nrec+1
         ENDIF
      CASE ('R', 'V')             ! RADAR OBSERVATION (delayed to next pass)
         radar=.TRUE.
      CASE DEFAULT
         lf=lench(file)
         WRITE(*,130) file(1:lf),obs(i)%type
         WRITE(ierrou,130) file(1:lf),obs(i)%type
130      FORMAT('write_rwo: cannot write in file ',A,' observation of unknown type ',A)
         numerr=numerr+1
   END SELECT
END DO

! Second pass: radar observations
IF(radar) THEN
   WRITE(unit,107) comcha,comcha
107 FORMAT(A1,' Object   Obser ====== Date =======  ============ Radar range/range rate (km or km/d) ============= ',  &
              'Station    Residual'/  &
           A1,' Design   K T N YYYY MM DD hh:mm:ss        Measure     Accuracy    rms    F      Bias       ',          &
              'Resid   TRX RCX     Chi   S')
   DO k=1,n
      IF(PRESENT(sort))THEN
         IF(sort)THEN
            i=time_ord(k)
         ELSE
            i=k
         ENDIF
      ELSE
         i=time_ord(k)
      ENDIF
      IF(obs(i)%type.eq.'R'.or.obs(i)%type.eq.'V')THEN
         CALL write_rwo_rad(unit,obs(i), obsw(i))
         nrec=nrec+1
      ENDIF
   END DO
END IF
IF(.not.PRESENT(unit0))CALL filclo(unit,' ')
END SUBROUTINE write_rwo

! Writes one record of radar observations
SUBROUTINE write_rwo_rad(unit, obs, obsw)
  INTEGER,        INTENT(IN) :: unit 
  TYPE(ast_obs),  INTENT(IN) :: obs
  TYPE(ast_wbsr), INTENT(IN) :: obsw
  INTEGER                    :: iday,month,year,iss,hh,mm 
  DOUBLE PRECISION           :: day,hour,rss
  CHARACTER(LEN=9)  :: chi_str
  CHARACTER(LEN=11) :: resid_rad,bias_rad
  CHARACTER(LEN=36) :: tmp2
  SELECT CASE (obs%type)
  CASE ('R', 'V')             ! RADAR OBSERVATION
     CALL mjddat(obs%time_utc,iday,month,year,hour)
     rss=hour*3600
     iss=NINT(rss)
     IF(iss<0)THEN
        WRITE(*,*) '**** write_rwo_rad: negative hour ****', hour, rss, iss
        STOP
     ELSEIF(ABS(rss-iss)>1.D-3)THEN
        WRITE(*,*) '**** write_rwo_rad: problem in rounding time ****',hour,rss,iss
        STOP
     ENDIF
     hh=iss/3600
     iss=iss-hh*3600
     mm=iss/60
     iss=iss-mm*60
     WRITE(tmp2,140) obs%objdes,obs%type,obs%tech,obs%note,  &
             &            year,month,iday,hh,mm,iss
140  FORMAT(1X,A9,1X,A1,1X,A1,1X,A1,1X,I4.4,1X,I2.2,1X,I2.2,1X,I2.2,':',I2.2,':',I2.2)
  END SELECT
  SELECT CASE (obs%type)
  CASE ('R')                  ! RADAR RANGE OBSERVATION
     bias_rad=radar_resid_string(obsw%bias_coord(1)*aukm,.TRUE.)
     resid_rad=radar_resid_string(obsw%res_coord(1)*aukm,obsw%resc_def)
     IF(obsw%resc_def) THEN
        IF(abs(obsw%chi).lt.9.9e5)THEN
           WRITE(chi_str,128) obsw%chi
        ELSE
           WRITE(chi_str,1128) obsw%chi
        ENDIF
     ELSE
        chi_str=' '
     END IF
     WRITE(unit,202) tmp2,obs%coord(1)*aukm,obs%acc_coord(1)*aukm,    &
          &                obsw%rms_coord(1)*aukm,obsw%force_w(1),   &
          &               bias_rad,resid_rad,obs%obscod_s,chi_str,obsw%sel_coord
  CASE ('V')                  ! RADAR RANGE RATE OBSERVATION
     bias_rad=radar_resid_string(obsw%bias_coord(2)*aukm,.TRUE.)
     resid_rad=radar_resid_string(obsw%res_coord(2)*aukm,obsw%resc_def)
     IF(obsw%resc_def) THEN
        IF(abs(obsw%chi).lt.9.9e5)THEN
           WRITE(chi_str,128) obsw%chi
128        FORMAT(F9.2)
        ELSE
           WRITE(chi_str,1128) obsw%chi
1128       FORMAT(1P,D9.2)
        ENDIF
     ELSE
        chi_str=' '
     END IF
     WRITE(unit,202) tmp2,obs%coord(2)*aukm,obs%acc_coord(2)*aukm,     &
          &                obsw%rms_coord(2)*aukm,obsw%force_w(2),   &
          &               bias_rad,resid_rad,obs%obscod_s,chi_str,obsw%sel_coord
202  FORMAT(A,F16.5,1X,F9.5,1X,F9.5,1X,L1,1X,A,1X,A,1X,A7,1X,A,1X,I1)
  END SELECT
END SUBROUTINE write_rwo_rad

! writes one record of a .rwo file with an optical observation
! If satellite, roving only the first line
SUBROUTINE write_rwo_opt(unit,obs,obsw,error,nrec)
  INTEGER,        INTENT(IN) :: unit 
  TYPE(ast_obs),  INTENT(IN) :: obs
  TYPE(ast_wbsr), INTENT(IN) :: obsw
  LOGICAL, INTENT(OUT) :: error
  INTEGER,     INTENT(INOUT) :: nrec
  INTEGER     :: iday,month,year,nd2,ra_h,ra_min,dec_deg,dec_min,hh,mm,iss 
  DOUBLE PRECISION  :: day,hour,ra_sec,resnor,dec_sec,rss
  CHARACTER(LEN=1)  :: ra_sign,dec_sign
  CHARACTER(LEN=5)  :: dec_sec_str
  CHARACTER(LEN=6)  :: ra_sec_str
  CHARACTER(LEN=8)  :: rmsa_str,biasa_str,rmsd_str,biasd_str
  CHARACTER(LEN=9)  :: resida_str,residd_str,chi_str
  CHARACTER(LEN=13) :: sday
  CHARACTER(LEN=19) :: mag_str
  CHARACTER(LEN=49) :: tmp1
  INTEGER, EXTERNAL :: lench
! check obs.type
  IF(obs%type.ne.'O'.and.obs%type.ne.'S')THEN
     error=.true.
     WRITE(ierrou,223) obs%type, nrec+1
223 FORMAT('write_rwo: attempt to write as optical type ',A1,'  : record',I5)
     numerr=numerr+1
  ENDIF
! Conversion of time
  CALL mjddat(obs%time_utc,iday,month,year,hour)
! Compose common part of output string (designation, observation type, time)
  day=iday+hour/24.d0
  nd2=-NINT(LOG10(obs%acc_time))
  CALL real2string(day,2,nd2,sday,error)
  IF(error)THEN
     WRITE(ierrou,121) nrec+1
121 FORMAT('write_rwo: error in DAY field: record',I5)
     numerr=numerr+1
  ENDIF
  WRITE(tmp1,120) obs%objdes,obs%type,obs%tech,obs%note,year,month,sday,obs%acc_time
120 FORMAT(1X,A9,1X,A1,1X,A1,1X,A1,1X,I4.4,1X,I2.2,1X,A,1X,1P,E10.3,0P)
! Conversion of RA
  CALL sessag(obs%coord(1)*hrad,ra_sign,ra_h,ra_min,ra_sec)
  IF(ra_sign == '-')THEN
     WRITE(*,*) '**** write_rwo_opt: negative right ascension **** ',obs%objdes,' ',ra_sign,ra_h,ra_min,ra_sec
     WRITE(ierrou,*)   '**** write_rwo_opt: negative right ascension **** ',obs%objdes,' ',ra_sign,ra_h,ra_min,ra_sec
     numerr=numerr+1 
  ENDIF
  WRITE(ra_sec_str,122) ra_sec
122 FORMAT(F6.3)
  IF(ra_sec_str(1:1) == ' ') ra_sec_str(1:1)='0'
  WRITE(rmsa_str,125)  obsw%rms_coord(1)*secrad*COS(obs%coord(2))
  WRITE(biasa_str,125) obsw%bias_coord(1)*secrad*COS(obs%coord(2))
!  WRITE(biasa_str,125) obsw%bias_coord(1)*secrad*COS(obs%coord(2)) wrong correction for cos(dec)
125 FORMAT(F8.3)
  IF(obsw%resc_def) THEN
     resnor=obsw%res_coord(1)*secrad*COS(obs%coord(2))
     IF(ABS(resnor) > 999.d0)THEN
        WRITE(resida_str,123) resnor
123     FORMAT(1P,E9.1,0P)
     ELSE
        WRITE(resida_str,124) resnor
124     FORMAT(F9.3)
     END IF
  ELSE
     resida_str=' '
  END IF
! Conversion of DEC
  CALL sessag(obs%coord(2)*degrad,dec_sign,dec_deg,dec_min,dec_sec)
  WRITE(dec_sec_str,126) dec_sec
  IF(dec_sec_str(1:1) == ' ') dec_sec_str(1:1)='0'
  WRITE(rmsd_str,125)  obsw%rms_coord(2)*secrad
  WRITE(biasd_str,125) obsw%bias_coord(2)*secrad
  IF(obsw%resc_def) THEN
     resnor=obsw%res_coord(2)*secrad
     IF(ABS(resnor) > 999.d0)THEN
        WRITE(residd_str,123) resnor
     ELSE
        WRITE(residd_str,124) resnor
     END IF
  ELSE
     residd_str=' '
  END IF
! Chi
  IF(obsw%resc_def) THEN
     IF(abs(obsw%chi).lt.9.9e5)THEN
        WRITE(chi_str,128) obsw%chi
128     FORMAT(F9.2)
     ELSE
        WRITE(chi_str,1128) obsw%chi
1128    FORMAT(1P,D9.2)
     ENDIF
  ELSE
     chi_str=' '
  END IF
! Conversion of MAG
  mag_str=' '
  IF(obs%mag_def) THEN
     mag_str(1:6)=obs%mag_str
! This is not necessary since accuracy can be determined from mag_str
    ! IF(obs%acc_mag >= 0.d0) WRITE(mag_str(8:12),126) obs%acc_mag
     IF(obsw%rms_mag >= 0.d0) WRITE(mag_str(8:12),126) obsw%rms_mag
126  FORMAT(F5.2)
     IF(obsw%resm_def) WRITE(mag_str(14:19),127) obsw%res_mag
127  FORMAT(F6.2)
  END IF
  IF(lench(obs%obscod_s)>3) THEN
!     lf=lench(file)
     WRITE(ierrou,131) obs%obscod_s !,nrec+1,file(1:lf)
131  FORMAT('write_rwo: abnormal observatory name "',A)
     numerr=numerr+1
  END IF
  WRITE(unit,201) tmp1,ra_h,ra_min,ra_sec_str,obs%acc_coord(1)*secrad*COS(obs%coord(2)),   &
&               rmsa_str,obsw%force_w(1),biasa_str,resida_str,            &
&               dec_sign,dec_deg,dec_min,dec_sec_str,obs%acc_coord(2)*secrad, &
&               rmsd_str,obsw%force_w(2),biasd_str,residd_str,            &
&               mag_str,obs%catcodmpc,obs%obscod_s,chi_str,obsw%sel_coord,obsw%sel_mag

201 FORMAT(A,1X,I2.2,1X,I2.2,1X,A6,1X,1P,E10.3,0P,1X,A8,1X,L1,1X,A8,A9, &
          1X,A1,I2.2,1X,I2.2,1X,A5,1X,1P,E10.3,0P,1X,A8,1X,L1,1X,A8,A9, &
          1X,A,2X,A2,1X,A3,1X,A,1X,I1,1X,I1)
  nrec=nrec+1
END SUBROUTINE write_rwo_opt

! Writes second record with position for satellite/roving observations
SUBROUTINE write_rwo_pos(unit,obs,obsw)
  INTEGER,        INTENT(IN) :: unit 
  TYPE(ast_obs),  INTENT(IN) :: obs
  TYPE(ast_wbsr), INTENT(IN) :: obsw
  DOUBLE PRECISION, DIMENSION(3)   :: equpos
  CHARACTER(LEN=33) :: tmp3
  CHARACTER(LEN=49) :: tmp1
  DOUBLE PRECISION  :: conv
  CHARACTER(LEN=13) :: sday
  INTEGER           :: iday,month,year,nd2
  DOUBLE PRECISION  :: day,hour
  LOGICAL           :: error
! Conversion of time
  CALL mjddat(obs%time_utc,iday,month,year,hour)
! Compose common part of output string (designation, observation type, time)
  day=iday+hour/24.d0
  nd2=-NINT(LOG10(obs%acc_time))
  CALL real2string(day,2,nd2,sday,error)
  WRITE(tmp1,120) obs%objdes,obs%type,obs%tech,obs%note,year,month,sday
120 FORMAT(1X,A9,1X,A1,1X,A1,1X,A1,1X,I4.4,1X,I2.2,1X,A)
  IF(obs%type == 'S') THEN
     tmp3=tmp1(1:33)
     tmp3(14:14)='s'
     SELECT CASE (obs%par_unit)
     CASE ('1')
        conv=aukm
     CASE ('2')
        conv=1.d0
     CASE DEFAULT
        STOP '**** write_rwo: internal error (04) ****'
     END SELECT
! Transformation of ecliptical coordinates into equatorial
!         CALL rotpn(rot,'ECLM','J2000',0.d0,'EQUM','J2000',0.d0)
     equpos=MATMUL(roteceq,obs%obspos)
     equpos=equpos*conv
     WRITE(unit,205) tmp3,obs%par_unit,equpos,obs%obscod_s(1:3)
205  FORMAT(A,1X,A1,3F12.4,1X,A)
  END IF
! Second record for roving observatory observations
  IF(obs%obscod_s == '247') THEN
     IF(obs%tech /= 'V') STOP '**** write_rwo: internal error (05) ****'
     IF(.NOT.ASSOCIATED(obs%geopos)) STOP '**** write_rwo: internal error (06) ****'
     tmp3=tmp1(1:33)
     tmp3(14:14)='v'
     WRITE(unit,206) tmp3,obs%geopos(1:2)*degrad,obs%geopos(3),obs%obscod_s(1:3)
206  FORMAT(A,1X,F10.6,1X,F10.6,1X,F8.1,1X,A)
  END IF
END SUBROUTINE write_rwo_pos

CHARACTER(LEN=11) FUNCTION radar_resid_string(value,defined)
  DOUBLE PRECISION, INTENT(IN) :: value
  LOGICAL,          INTENT(IN) :: defined
  LOGICAL :: fixed
  radar_resid_string=' '
  IF(.NOT.defined) RETURN
  fixed=.true.
  IF(value>99999.D0) fixed=.false.
  IF(value<-9999.d0) fixed=.false.
  IF(fixed) THEN
     WRITE(radar_resid_string,101) value
     101 FORMAT(F11.5)
  ELSE
     WRITE(radar_resid_string,102) value
     102 FORMAT(1P,E11.4)
  END IF
END FUNCTION radar_resid_string

SUBROUTINE write_rwo_header(unit, error_model, nrec, rms, rmsmag)
INTEGER,          INTENT(IN) :: unit !output unit
CHARACTER(LEN=*), INTENT(IN) :: error_model ! error model file name
INTEGER,       INTENT(INOUT) :: nrec ! record counter 
DOUBLE PRECISION, INTENT(IN), OPTIONAL :: rms         ! RMS of astrometric fit
DOUBLE PRECISION, INTENT(IN), OPTIONAL :: rmsmag      ! RMS of photometric fit
INTEGER lench, ll
INCLUDE 'parcmc.h90'     ! comment character
! write header

WRITE(unit,101) rwo_version
101 FORMAT('version = ',I3)
nrec=nrec+1
ll=lench(error_model)
IF(ll > 0) THEN
   WRITE(unit,102) error_model(1:ll)
102 FORMAT('errmod  = ''',A,'''')
ELSE
   WRITE(unit,102) ' '
END IF
nrec=nrec+1
IF(PRESENT(rms))THEN
   IF(rms >= 0.d0) THEN
      WRITE(unit,103) rms
103   FORMAT('RMSast  = ',1P,E13.5)
      nrec=nrec+1
   ENDIF
ENDIF
IF(PRESENT(rmsmag))THEN
   IF(rmsmag >= 0.d0) THEN
      WRITE(unit,104) rmsmag
104   FORMAT('RMSmag  = ',1P,E13.5)
      nrec=nrec+1
   ENDIF
ENDIF
WRITE(unit,105)
105 FORMAT('END_OF_HEADER')
WRITE(unit,106) comcha,comcha
nrec=nrec+3
106 FORMAT(A1,' Object   Obser ============= Date ============= ================== Right Ascension =================', &
         '  ================= Declination ===================== ==== Magnitude ==== Ast Obs  Residual SEL'/  &
      A1,' Design   K T N YYYY MM DD.dddddddddd   Accuracy HH MM SS.sss  Accuracy      RMS  F     Bias    Resid sDD MM SS.ss', &
         '  Accuracy      RMS  F     Bias    Resid Val  B   RMS  Resid Cat Cod       Chi A M')

END SUBROUTINE write_rwo_header

! input from rwo files - structured
SUBROUTINE read_rwo(file,obs,obsw,n,error_model,rms,rmsmag,unit0,eof0)
! INTERFACE
  CHARACTER(LEN=*),             INTENT(IN)            :: file        ! input file, to be opened
          ! unless the optional variable unit0 is present
  TYPE(ast_obs),  DIMENSION(:), INTENT(OUT)           :: obs         ! observations
  TYPE(ast_wbsr), DIMENSION(:), INTENT(OUT)           :: obsw        ! observations weights/residuals
  INTEGER,                      INTENT(OUT)           :: n           ! number of observations found
  CHARACTER(LEN=*),             INTENT(OUT)           :: error_model ! error model file name 
  DOUBLE PRECISION,             INTENT(OUT), OPTIONAL :: rms         ! RMS of astrometric fit
  DOUBLE PRECISION,             INTENT(OUT), OPTIONAL :: rmsmag      ! RMS of photometric fit
  INTEGER,                      INTENT(IN) , OPTIONAL :: unit0 ! input unit if file is already opened
          ! then the routine must read until the obj.designation changes, then leave open
  LOGICAL,                      INTENT(OUT), OPTIONAL :: eof0 ! end of file if file was opened
          ! note eof0 is required if unit0 is present
  INCLUDE 'parobx.h90' ! observation numbers: maximum
! END INTERFACE
  LOGICAL error, eof, stored
  INTEGER :: unit,nr,kr,lf,j
  TYPE(ast_obs)  :: obs1
  TYPE(ast_wbsr) :: obsw1
  INTEGER, EXTERNAL :: lench
  CHARACTER*(name_len) ::  name0, name 
  SAVE obs1,obsw1,stored
  DATA stored /.false./
! ==========================================================
! check operating mode
  lf=lench(file)
  IF(PRESENT(unit0))THEN
     IF(.not.PRESENT(eof0))THEN
        WRITE(*,302) file(1:lf)
302     FORMAT('read_rwo: missing eof0 to warn of end of file ',A)
        STOP
     ENDIF
     IF(PRESENT(rms))THEN
       rms=-99.d0
     ENDIF
     IF(PRESENT(rmsmag))THEN
       rmsmag=-99.d0
     ENDIF
     unit=unit0
  ELSE
! open file
     CALL filopn(unit,file,'OLD')
     n=0 ! observations counter
     CALL rdfnam(unit,file,nr)
!read header; should count records...
     CALL read_rwo_header(unit, error_model, error, nr, rms, rmsmag)
     IF(error)THEN
        WRITE(ierrou,301) file(1:lf)
301     FORMAT('read_rwo: input error in header of file ',A)
     ENDIF
   ENDIF
! main loop on observations
  DO j=1,nobx
     IF(.not.stored)THEN
        CALL read_rwo_rec(unit,nr,obs1,obsw1,error,eof)
! redifine rejection when tech is equal 'x'
        IF(obs1%tech.EQ.'x'.OR.obs1%tech.EQ.'X') THEN
           obsw1%sel_coord=0
        END IF
        IF(PRESENT(eof0))THEN
           eof0=eof
        ENDIF
        IF(eof)THEN
           stored=.false.
           EXIT
        ELSEIF(error)THEN 
           WRITE(*,202) nr,file(1:lf)
           IF(ierrou > 0) WRITE(ierrou,202) nr,file(1:lf)
202        FORMAT('read_rwo: input error at record',I6,' of file ',A)
           numerr=numerr+1
           stored=.false.
        ENDIF
     ELSE
        IF(PRESENT(eof0))THEN
           eof0=.false.
        ENDIF
     ENDIF
     IF(PRESENT(unit0))THEN
        IF(j.eq.1)THEN
           name0=obs1%objdes
           stored=.false.
        ELSE
! check if the designation has changed
           name=obs1%objdes
           IF(name.ne.name0)THEN
! store for next call
              stored=.true.                 
              EXIT
           ELSE
              stored=.false.
           ENDIF
        ENDIF
     ENDIF
! observation has been read
     n=n+1
     obs(n)=obs1
     obsw(n)=obsw1
  END DO
  IF(.not.PRESENT(unit0))CALL filclo(unit,' ')
END SUBROUTINE read_rwo

! ======================================
! read_rwo_rec
! ======================================
SUBROUTINE read_rwo_rec(unit,nr,obs1,obsw1,error, eof)
USE station_coordinates
USE reference_systems
  INTEGER,         INTENT(IN)   :: unit      ! input unit
  INTEGER,        INTENT(INOUT) :: nr        ! record counter
  TYPE(ast_obs),   INTENT(OUT)  :: obs1      ! observation
  TYPE(ast_wbsr),  INTENT(OUT)  :: obsw1     ! weights residuals etc.
  LOGICAL,         INTENT(OUT)  :: error,eof ! error flag, end of file
  INCLUDE 'parcmc.h90'                       ! comment character
  CHARACTER(LEN=200) :: record
  DOUBLE PRECISION                 :: day,sec,sect,ra_sec,dec_sec,coord1,acc_coord1,rms_coord1,conv
  DOUBLE PRECISION, DIMENSION(3)   :: equpos,geopos,bfpos0
  CHARACTER(LEN=1)   :: dec_sign,col1,col2
  CHARACTER(LEN=3)   :: scale,obscod_s1
  CHARACTER(LEN=5)   :: mag_field
  CHARACTER(LEN=8)   :: rmsa_str,rmsd_str
  CHARACTER(LEN=9)   :: resida_str,residd_str,chi_str
  CHARACTER(LEN=11)  :: resid_rad,bias_rad
  CHARACTER(LEN=19)  :: mag_str
  CHARACTER(LEN=33)  :: tmp3
  LOGICAL            :: force1
  INTEGER :: year,month,iday,mjd,mjdt,ra_h,ra_min,dec_deg,dec_min,len_mag_field,pos_point,n_dec_mag
  INTEGER :: hh,mm,iss,ipr0,ipr1,code_tx,code_rx
  DOUBLE PRECISION, EXTERNAL :: tjm1
  CHARACTER(LEN=16)  :: stname
! ======================================
! setup of flags
  error=.false.
  eof=.false.
! read one record
3 CONTINUE
  READ(unit,201,END=2) record
201 FORMAT(A)
  nr=nr+1
  IF(record(1:1) == comcha) GOTO 3
  obs1=undefined_ast_obs
  obsw1=undefined_ast_wbsr
  READ(record,110, ERR=1, END=2) obs1%objdes,obs1%type,obs1%tech,obs1%note
110 FORMAT(1X,A9,1X,A1,1X,A1,1X,A1)
  SELECT CASE (obs1%type)
  CASE ('O', 'S')                  ! OPTICAL/SATELLITE OBSERVATION
! time of observations
     READ(record(18:),120,ERR=1) year,month,day,obs1%acc_time
120  FORMAT(I4,1X,I2,1X,F13.10,1X,E10.3)
     iday=day
     sec=(day-iday)*86400.d0
     mjd=NINT(tjm1(iday,month,year,0.d0))
     IF(year < 1972) THEN
        scale='UT1'
     ELSE
        scale='UTC'
     ENDIF
     CALL cnvtim(mjd,sec,scale,mjdt,sect,'TDT')
     obs1%time_utc=mjd+sec/86400.d0
     obs1%time_tdt=mjdt+sect/86400.d0
! change_format='INPUT'
     IF(change_format.ne.'INPUT')THEN
! read alpha, delta with weights, bias, residuals... also magnitude, chi, selection flags
        READ(record(51:),521,ERR=1) ra_h,ra_min,ra_sec,obs1%acc_coord(1),rmsa_str,obsw1%force_w(1),&
&              obsw1%bias_coord(1),resida_str,                                              &
&              dec_sign,dec_deg,dec_min,dec_sec,obs1%acc_coord(2),rmsd_str,obsw1%force_w(2),&
&              obsw1%bias_coord(2),residd_str,                                              &
&              mag_str,obs1%catcodmpc,obs1%obscod_s,chi_str,obsw1%sel_coord,obsw1%sel_mag
521     FORMAT(I2,1X,I2,1X,F6.6,1X,E10.3,1X,A8,1X,L1,1X,F8.3,A9, &
&           1X,A1,I2,1X,I2,1X,F5.2,1X,E10.3,1X,A8,1X,L1,1X,F8.3,A9, &
&           1X,A19,2X,A2,1X,A3,1X,A9,1X,I1,1X,I1)
     ELSE
! read alpha, delta with weights, bias, residuals... also magnitude, chi, selection flags, no catcodmpc
        READ(record(51:),121,ERR=1) ra_h,ra_min,ra_sec,obs1%acc_coord(1),rmsa_str,obsw1%force_w(1),&
&              obsw1%bias_coord(1),resida_str,                                              &
&              dec_sign,dec_deg,dec_min,dec_sec,obs1%acc_coord(2),rmsd_str,obsw1%force_w(2),&
&              obsw1%bias_coord(2),residd_str,                                              &
&              mag_str,obs1%obscod_s,chi_str,obsw1%sel_coord,obsw1%sel_mag
121     FORMAT(I2,1X,I2,1X,F6.6,1X,E10.3,1X,A8,1X,L1,1X,F8.3,A9, &
&           1X,A1,I2,1X,I2,1X,F5.2,1X,E10.3,1X,A8,1X,L1,1X,F8.3,A9, &
&           1X,A19,1X,A3,1X,A9,1X,I1,1X,I1)
     ENDIF
     obs1%acc_coord=obs1%acc_coord/secrad
     IF(rmsa_str == ' ') THEN
        obsw1%rms_coord(1)=-1.d0
     ELSE
        READ(rmsa_str,*,ERR=1) obsw1%rms_coord(1)
        obsw1%rms_coord(1)=obsw1%rms_coord(1)/secrad
     END IF
     IF(rmsd_str == ' ') THEN
        obsw1%rms_coord(2)=-1.d0
     ELSE
        READ(rmsd_str,*,ERR=1) obsw1%rms_coord(2)
        obsw1%rms_coord(2)=obsw1%rms_coord(2)/secrad
     END IF
     obsw1%bias_coord=obsw1%bias_coord/secrad
     obs1%coord(1)=15.d0*(ra_h*3600.d0+ra_min*60.d0+ra_sec)/secrad
     obs1%coord(2)=(dec_deg*3600.d0+dec_min*60.d0+dec_sec)/secrad
! correct for cos(DEC)
     obs1%acc_coord(1)=obs1%acc_coord(1)/COS(obs1%coord(2))
     obsw1%bias_coord(1)=obsw1%bias_coord(1)/COS(obs1%coord(2))
     obsw1%rms_coord(1)=obsw1%rms_coord(1)/COS(obs1%coord(2))
     IF(dec_sign == '-') obs1%coord(2)=-obs1%coord(2)
     obsw1%resc_def=(resida_str /= ' ' .AND. residd_str /= ' ')
     IF(obsw1%resc_def) THEN
        READ(resida_str,*,ERR=1) obsw1%res_coord(1)
        READ(residd_str,*,ERR=1) obsw1%res_coord(2)
        obsw1%res_coord=obsw1%res_coord/secrad
        obsw1%res_coord(1)=obsw1%res_coord(1)/COS(obs1%coord(2))
        READ(chi_str,*,ERR=1) obsw1%chi
     ELSE
        obsw1%chi=-1.d0
     END IF
! station
     CALL statcode(obs1%obscod_s,obs1%obscod_i)
! magnitude
     obs1%mag_str=mag_str(1:6)
     mag_field=mag_str(1:5)
     CALL rmsp(mag_field,len_mag_field)
     mag_field=mag_str(1:5) !WARNING: to fix 644 mistake
     IF(len_mag_field>0) THEN
        READ(mag_field,*,ERR=1) obs1%mag
        obs1%mag_def=.true.
        pos_point=INDEX(mag_field,'.')
        IF(pos_point == 0) THEN
           obs1%acc_mag=1
        ELSE
           n_dec_mag=len_mag_field-pos_point
           obs1%acc_mag=10.0d0**(-n_dec_mag)
        END IF
     ELSE
        obs1%mag_def=.false.
        obs1%mag=0.d0
        obs1%acc_mag=99.d9
     END IF
     obs1%mag_band=mag_str(6:6)
     IF(mag_str(8:12) == ' ') THEN
        obsw1%rms_mag=-1.d0
     ELSE
        READ(mag_str(8:12),*,ERR=1) obsw1%rms_mag
     END IF
     IF(mag_str(14:19) == ' ') THEN
        obsw1%resm_def=.false.
        obsw1%res_mag=0.d0
     ELSE
        READ(mag_str(14:19),*,ERR=1) obsw1%res_mag
        obsw1%resm_def=.true.
     END IF
! Observer's position is computed for ground-based observations and read from the following record
! for roving observatory and satellite observations
     SELECT CASE (obs1%type)
     CASE ('O')
        SELECT CASE (obs1%obscod_s)
        CASE ('247')
           IF(obs1%tech /= 'V') GOTO 1
           READ(unit,206,ERR=1) tmp3,geopos,obscod_s1
206        FORMAT(A33,1X,F10.6,1X,F10.6,1X,F8.1,1X,A)
           nr=nr+1
           IF(tmp3(14:14) /= 'v') GOTO 1
           tmp3(14:14)='V'
           IF(tmp3 /= record(1:33)) GOTO 1
           IF(obscod_s1 /= obs1%obscod_s(1:3)) GOTO 1
           geopos(1:2)=geopos(1:2)*radeg
           IF(ASSOCIATED(obs1%geopos)) DEALLOCATE(obs1%geopos)
           ALLOCATE(obs1%geopos(3))
           obs1%geopos=geopos
           CALL geodetic_to_cartesian(geopos(1),geopos(2),geopos(3),bfpos0)
           bfpos0=bfpos0/(aukm*1d3)
           CALL observer_position(obs1%time_tdt,obs1%obspos,obs1%obsvel,BFPOS=bfpos0)
        CASE DEFAULT
           CALL observer_position(obs1%time_tdt,obs1%obspos,obs1%obsvel,OBSCODE=obs1%obscod_i)
        END SELECT
     CASE ('S')
        READ(unit,205,ERR=1) tmp3,obs1%par_unit,equpos,obscod_s1
205     FORMAT(A33,1X,A1,3F12.4,1X,A)
        nr=nr+1
        tmp3(14:14)='S'
        IF(tmp3 /= record(1:33)) GOTO 1
        IF(obscod_s1 /= obs1%obscod_s(1:3)) GOTO 1
        SELECT CASE (obs1%par_unit)
        CASE ('1')
           conv=1.d0/aukm
        CASE ('2')
           conv=1.d0
        CASE DEFAULT
           GOTO 1
        END SELECT
        equpos=equpos*conv
! Transformation of equatorial coordinates into ecliptical
        obs1%obspos=MATMUL(roteqec,equpos)
        obs1%obsvel=0.d0
     CASE DEFAULT
        STOP '**** read_rwo: internal error (01) ****'
     END SELECT

  CASE ('R', 'V')                  ! RADAR OBSERVATION (RANGE OR RANGE RATE)
     READ(record(18:),122,ERR=1) year,month,iday,hh,col1,mm,col2,iss
122  FORMAT(I4,1X,I2,1X,I2,1X,I2,A1,I2,A1,I2)
     IF(col1 /= ':') GOTO 1
     IF(col2 /= ':') GOTO 1
     mjd=NINT(tjm1(iday,month,year,0.d0))
     sec=iss+mm*60+hh*3600
     IF(year < 1972) THEN
        scale='UT1'
     ELSE
        scale='UTC'
     ENDIF
     CALL cnvtim(mjd,sec,scale,mjdt,sect,'TDT')
     obs1%time_utc=mjd+sec/86400.d0
     obs1%time_tdt=mjdt+sect/86400.d0
     READ(record(37:),123,ERR=1) coord1,acc_coord1,rms_coord1,force1, &
          &       bias_rad,resid_rad,obs1%obscod_s,chi_str,obsw1%sel_coord
123  FORMAT(F16.5,1X,F9.5,1X,F9.5,1X,L1,1X,A,1X,A,1X,A7,1X,A,1X,I1)
     CALL statcode(obs1%obscod_s(1:3),code_tx)
     CALL statcode(obs1%obscod_s(5:7),code_rx)
     obs1%obscod_i=code_tx*10000+code_rx
     CALL obscoo(code_tx,obs1%obspos,stname)
     CALL obscoo(code_rx,obs1%obsvel,stname)
     obs1%mag_def=.false.
     SELECT CASE (obs1%type)
     CASE ('R')
        ipr1=1
     CASE ('V')
        ipr1=2
     END SELECT
     ipr0=3-ipr1
     obs1%coord(ipr1)=coord1/aukm
     obs1%coord(ipr0)=0
     obs1%acc_coord(ipr1)=acc_coord1/aukm
     obsw1%rms_coord(ipr1)=rms_coord1/aukm
     obsw1%force_w(ipr1)=force1
     obsw1%force_w(ipr0)=.false.
     IF(bias_rad == ' ') THEN
        obsw1%bias_coord(ipr1)=0.d0
     ELSE
        READ(bias_rad,*,ERR=1) obsw1%bias_coord(ipr1)
        obsw1%bias_coord(ipr1)=obsw1%bias_coord(ipr1)/aukm
     END IF
     obsw1%bias_coord(ipr0)=0.d0
     obsw1%resc_def=(resid_rad /= ' ')
     IF(obsw1%resc_def) THEN
        READ(resid_rad,*,ERR=1) obsw1%res_coord(ipr1)
        obsw1%res_coord(ipr1)=obsw1%res_coord(ipr1)/aukm
     ELSE
        obsw1%res_coord(ipr1)=0.d0
     END IF
     obsw1%res_coord(ipr0)=0.d0
     IF(chi_str /= ' ') THEN
        READ(chi_str,*,ERR=1) obsw1%chi
     ELSE
        obsw1%chi=-9.d99
     END IF
     obs1%acc_time=10d-10
  CASE DEFAULT
     WRITE(*,203) obs1%type,nr
     WRITE(ierrou,203) obs1%type,nr
203  FORMAT('read_rwo: unknown observation type ',A,' at record',I6,' of file ',A)
     numerr=numerr+1
     error=.true.
  END SELECT
  RETURN
1 CONTINUE ! error cases
  error=.true.
  RETURN
2 CONTINUE ! regular end of file
  eof=.true.
END SUBROUTINE read_rwo_rec

! ========================================
!  read_rwo_header
! ========================================
SUBROUTINE  read_rwo_header(unit, error_model, error, nr, rms, rmsmag)
  INTEGER,          INTENT(IN)  :: unit               !output unit
  CHARACTER(LEN=*), INTENT(OUT) :: error_model        ! error model file name
  LOGICAL,          INTENT(OUT) :: error              ! error flag
  INTEGER,          INTENT(OUT) :: nr               ! record counter
  DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: rms      ! RMS of astrometric fit
  DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: rmsmag   ! RMS of photometric fit
  EXTERNAL :: rdfnam,rdfint,rdfcha,rdfrea,cnvtim,rmsp
  INTEGER :: version,kr
  LOGICAL :: found
! ========================================
  error=.false.
  nr=0
  CALL rdfint(unit,'version',.false.,version,found,kr)
  IF(found) THEN
     nr=nr+1
     IF(version == 1) THEN
        change_format='INPUT'
     ELSEIF(version == rwo_version) THEN
        change_format=' '
     ELSE
        WRITE(*,200) version
200     FORMAT('Unknown version ',I2)
        WRITE(ierrou,200) version
        numerr=numerr+1
        error=.true.
        RETURN
     END IF
  ELSE
     WRITE(*,*) 'no rwo version  in header'
     WRITE(*,*) 'may be old f77 format '
     WRITE(ierrou,*) 'no rwo version  in header'
     numerr=numerr+1
     error=.true.
     RETURN
  END IF
  error_model=' '
  CALL rdfcha(unit,'errmod',.false.,error_model,found,kr)
  IF(.not.found)THEN
     WRITE(*,222) version
222  FORMAT('No error_model, version ',I2)
     WRITE(ierrou,200) version
     numerr=numerr+1
     error=.true.
     RETURN
  ENDIF
  nr=nr+1
  IF(PRESENT(rms)) THEN
     rms=-99.d0
     CALL rdfrea(unit,'RMSast',.false.,rms,found,kr)
     If(found) nr=nr+1
  END IF
  IF(PRESENT(rmsmag)) THEN
     rmsmag=-99.d0
     CALL rdfrea(unit,'RMSmag',.false.,rmsmag,found,kr)
     If(found) nr=nr+1
  END IF
END SUBROUTINE read_rwo_header
! END SUBMODULE iorwo

! SUBMODULE iompc
! read file of MPC observation data

SUBROUTINE mpc_obs_input(obs_found,obs,nobs,unit,filnam,eof,ons_mode)

LOGICAL, INTENT(OUT) :: obs_found ! true if there is something
INTEGER, INTENT(OUT) :: nobs ! number of observations
TYPE(ast_obs), INTENT(OUT), DIMENSION(:) :: obs ! vector of observations
INTEGER, INTENT(IN), OPTIONAL :: unit ! input unit if file is already opened
          ! then the routine must read until the obj.designation changes, then leave open
CHARACTER(*), INTENT(IN), OPTIONAL :: filnam ! input file name if it is to be opened

          ! then the routine must open, read until end of file, then close
LOGICAL, INTENT(OUT), OPTIONAL :: eof ! end of file flag in case the file was open
LOGICAL, INTENT(IN), OPTIONAL :: ons_mode ! if present, reading an ONS
                ! requiring a new designations manufactured by the perl routine
!end interface

INTEGER iunit, nobsx, n, lf
INTEGER, EXTERNAL :: lench 

! record and logical controls for mpcrec_transform
CHARACTER*(80) mpcrec,mpcrec2
LOGICAL skip_rad, error, complete, rad_in_mpcfile

Character*9 name ! for message

! check that one or the other operations mode is selected (replacing either mpcin or mpcin3)

obs_found=.FALSE.
rad_in_mpcfile=.FALSE.

nobs=0
IF(PRESENT(unit).and.PRESENT(filnam))THEN
   WRITE(*,*) ' mpc_obs_input: both options present'
   STOP
ENDIF

nobsx=SIZE(obs,1)

! if necessary open file
IF(PRESENT(filnam))THEN
   CALL filopn(iunit,filnam,'old')
ELSEIF(PRESENT(unit))THEN
   iunit=unit
   IF(.not.PRESENT(eof))THEN
      WRITE(*,*)' mpc_obs_input: inconsistent optional arguments'
   ENDIF
ELSE
   WRITE(*,*) ' mpc_obs_input: no options present'
   STOP
ENDIF

nobs=0
DO n=1,nobsx
   nobs=nobs+1
   READ(iunit,'(a)',end=10, err=11) mpcrec
! blank line is like file termination
   IF(lench(mpcrec).eq.0) GOTO 10
! interpret one record
   IF(PRESENT(ons_mode))THEN 
      CALL mpcrec_transform(mpcrec,obs(nobs),complete,error,skip_rad,ons_mode)
   ELSE  
      CALL mpcrec_transform(mpcrec,obs(nobs),complete,error,skip_rad)
   ENDIF
   IF(n.eq.1)name=obs(1)%objdes
! if error skip one record and hope it can go ahead
   IF(error) THEN
      WRITE(ierrou,*)'object ',name,'error at MPC observation record no ',n
      numerr=numerr+1
      nobs=nobs-1
      CYCLE
   ENDIF
! skipping MPC radar data and issue a warning
   IF(skip_rad)THEN
      nobs=nobs-1
      rad_in_mpcfile=.TRUE.
      CYCLE
   ENDIF
! if asteroid is changed, go back of the only record read so far
   IF(PRESENT(unit))THEN
      IF(obs(nobs)%objdes.ne.obs(1)%objdes)THEN
         eof=.FALSE.
         nobs=nobs-1
         BACKSPACE(iunit)
         RETURN
      ENDIF
   ENDIF
! handle cases not including satellite and roving observer observations
   IF(complete)THEN
! find observer position
      CALL mpcrec_add_obspos(obs(nobs),error)
! if error skip the two records and hope it can go ahead
      IF(error) THEN
         WRITE(ierrou,*)'object ',name,'error in position at MPC record no ',n
         numerr=numerr+1
         nobs=nobs-1
         CYCLE
      ENDIF
   ELSE
! handle roving observer/satellite position
      IF(obs(nobs)%type.eq.'S'.or.obs(nobs)%obscod_s.eq.'247')THEN
! read second line with satellite/roving obs. position
         READ(iunit,'(a)',end=10) mpcrec2
         IF(mpcrec2(15:15).eq.'X')THEN
! ignore X record if inserted between S and s
            READ(iunit,'(a)',end=10) mpcrec2
            CALL mpcrec_add_obspos(obs(nobs),error,mpcrec,mpcrec2)
         ELSEIF(mpcrec2(15:15).eq.'s'.or.mpcrec2(15:15).eq.'v')THEN
!            WRITE(*,*)' processing satellite/roving observer data, t=',obs(nobs)%time_utc
            CALL mpcrec_add_obspos(obs(nobs),error,mpcrec,mpcrec2)
         ELSE
            error=.true.
         ENDIF
! if error skip the two records and hope it can go ahead
         IF(error) THEN
            IF(obs(nobs)%type.eq.'S')THEN
               WRITE(ierrou,*)'object ',name,'error in satellite position at MPC record no ',n
            ELSEIF(obs(nobs)%obscod_s.eq.'247')THEN
               WRITE(ierrou,*)'object ',name,'error in roving observer pos. at MPC record no ',n
            ENDIF
            numerr=numerr+1
            nobs=nobs-1
            CYCLE
         ENDIF
      ELSE
! logical error
         WRITE(ierrou,*)' object ', name,' not complete but not roving/sat, rec= ',n
      ENDIF
   ENDIF
   obs_found=.TRUE.
ENDDO
! too many observations
WRITE(ierrou, *)' too many MPC observations, record no ',n
numerr=numerr+1

! error
11 WRITE(*,*) 'mpc_obs_input: error in reading at nobs =',nobs
! end of file
10 nobs=nobs-1

! If necessary close file and issue warning for presence of radar observations
IF(PRESENT(filnam))THEN
   IF(rad_in_mpcfile) THEN
      lf=lench(filnam)
!      WRITE(ierrou,120) filnam(1:lf)
!      WRITE(ierrou,122)
!      numerr=numerr+1
   END IF
   CALL filclo(iunit,' ')
ELSE
   IF(rad_in_mpcfile) THEN
!      WRITE(ierrou,121)
!      WRITE(ierrou,122)
!      numerr=numerr+1
   END IF
   eof=.TRUE.
END IF
120 FORMAT('MPC file ',A,' contains radar observations')
121 FORMAT('Some input MPC file contains radar observations')
122 FORMAT(5X,'Radar observations cannot be handled correctly by MPC format:'/   &
           5X,'please use appropriate .rad file instead')

END SUBROUTINE mpc_obs_input
! Transformation of a MPC record into an observation (TYPE ast_obs)
! This routine works on a single record and therefore cannot handle correctly satellite and radar observations
SUBROUTINE mpcrec_transform(mpcrec,obs,complete,error,skip_rad,ons_mode)
USE station_coordinates
USE output_control

CHARACTER(LEN=*),                INTENT(IN)  :: mpcrec      ! MPC record
TYPE(ast_obs),                   INTENT(OUT) :: obs         ! astrometric observations
LOGICAL,                         INTENT(OUT) :: error       ! error conversion flag,
LOGICAL,                         INTENT(OUT) :: complete    ! input complete flag: if(.not.complete) the next record has to be read
LOGICAL,                         INTENT(OUT) :: skip_rad    ! radar skip flag: if(skip_rad) a radar data has been skipped
LOGICAL, INTENT(IN), OPTIONAL :: ons_mode ! if present, reading an ONS
                ! requiring a new designations manufactured by the perl routine
!end interface
CHARACTER(LEN=80) :: error_code          ! internal error code
LOGICAL           :: err_tmp             ! error codes for called subroutines
INTEGER           :: year,month,day
DOUBLE PRECISION  :: seconds_utc,seconds_tdt
CHARACTER(LEN=9)  :: chdate
INTEGER           :: len_chdate,pos_point,n_dec_day,mjd_utc,mjd_tdt,len_error_code,len_mag_field,n_dec_mag
CHARACTER(LEN=3)  :: timescale
CHARACTER(LEN=5)  :: mag_field

! Obsolete CALLs (to be substituted by modules)
DOUBLE PRECISION, EXTERNAL :: tjm1
INTEGER,          EXTERNAL :: lench
INTEGER le
error=.true.
error_code='undefined'

! Default values
obs=undefined_ast_obs

! Observation type and technology
obs%tech=mpcrec(15:15)
! error if it is a second record with observer position
IF(obs%tech.eq.'v'.or.obs%tech.eq.'s')THEN
   error_code='second record with obs. pos. given as first record'
   GOTO 10
ELSEIF(obs%tech.eq."'")THEN  ! anti-quote safety
   obs%tech=' '
ELSEIF(obs%tech.eq.'"')THEN
   obs%tech=' '
ELSEIF(obs%note.eq."\")THEN !" anti-backslash safety
   obs%note=' ' 
ENDIF
! Radar and satellite observations are not handled at once, another record has to be read
obs%type='O'
skip_rad=.FALSE.
IF(obs%tech.eq.'S') THEN
   complete=.FALSE.
   obs%type='S'
ELSE
   complete=.TRUE.
END IF

IF(PRESENT(ons_mode))THEN
   IF(ons_mode)THEN
! Object Name: read designation created by the perl script, given by the observer, etc.
      obs%objdes=mpcrec(4:12)
      CALL rmsp(obs%objdes,le)
   ELSE
! Object Name
      CALL iaucod(mpcrec(1:12),obs%objdes,err_tmp)
      IF (err_tmp) THEN
         error_code='conversion of object name'
         GOTO 10
      END IF
   ENDIF
ELSE
! Object Name
   CALL iaucod(mpcrec(1:12),obs%objdes,err_tmp)
   IF (err_tmp) THEN
      error_code='conversion of object name'
      GOTO 10
   END IF
ENDIF
! Radar observations are not read from MPC files
IF(obs%tech.eq.'R' .OR. obs%tech.eq.'r') THEN
   error=.FALSE.
   skip_rad=.TRUE.
   RETURN
END IF
! Default value for technology descriptor is 'P' (Photographic)
IF(obs%tech.eq.' ') obs%tech='P'

! Note on observation circumstances
obs%note=mpcrec(14:14)
IF(obs%note.eq."'")THEN ! anti-quote safety
   obs%note=' '
ELSEIF(obs%note.eq.'"')THEN
   obs%note=' '
ELSEIF(obs%note.eq.'\')THEN !' anti-backslash safety
   obs%note=' '
ENDIF
! Observatory code
error_code='input field: observatory code'

READ(mpcrec(78:80),*,ERR=10) obs%obscod_s
CALL statcode(obs%obscod_s,obs%obscod_i)
IF(obs%obscod_i.eq.247)THEN
  complete=.false.
ENDIF

! Read time fields
error_code='input field: year, month and date'
READ(mpcrec,100,ERR=10) year,month,chdate
100 FORMAT(15X,I4,1X,I2,1X,A9)
len_chdate=lench(chdate)
IF(len_chdate.LE.0) GOTO 10

! Time accuracy
pos_point=INDEX(chdate,'.')
IF(pos_point == 0) THEN
   error_code='input field: day'
   READ(chdate,*,ERR=10) day
   seconds_utc=0
   obs%acc_time=1
ELSE
   error_code='input field: day'
   READ(chdate(1:pos_point-1),*,ERR=10) day
   error_code='input field: seconds'
   READ(chdate(pos_point:),*,ERR=10) seconds_utc
   seconds_utc=seconds_utc*86400
   n_dec_day=len_chdate-pos_point
   obs%acc_time=10.0d0**(-n_dec_day)
END IF

! Conversion of time fields into MJD
IF(year.LT.1972) THEN
   timescale='UT1'
ELSE
   timescale='UTC'
END IF
mjd_utc=NINT(tjm1(day,month,year,0.d0))
CALL cnvtim(mjd_utc,seconds_utc,timescale,mjd_tdt,seconds_tdt,'TDT')
obs%time_utc=mjd_utc+seconds_utc/86400.d0
obs%time_tdt=mjd_tdt+seconds_tdt/86400.d0

! Conversion of RA and DEC
CALL rdanga(mpcrec(45:56),obs%coord(2),obs%acc_coord(2),err_tmp)
IF (err_tmp) THEN
   error_code='conversion of DEC'
   GOTO 10
ELSE
   obs%coord(2)=obs%coord(2)*radeg
   obs%acc_coord(2)=obs%acc_coord(2)*radeg
END IF
CALL rdanga(mpcrec(33:44),obs%coord(1),obs%acc_coord(1),err_tmp)
IF (err_tmp) THEN
   error_code='conversion of RA'
   GOTO 10
ELSE
   obs%coord(1)=obs%coord(1)*radh
   obs%acc_coord(1)=obs%acc_coord(1)*radh
   obs%acc_coord(1)=obs%acc_coord(1)/COS(obs%coord(2))
END IF

! Magnitude
obs%mag_str=mpcrec(66:71)
obs%mag_band=mpcrec(71:71)
mag_field=mpcrec(66:70) !WARNING: to fix 644 mistake
CALL rmsp(mag_field,len_mag_field)
mag_field=mpcrec(66:70)
IF(len_mag_field>0) THEN
   error_code='input field: magnitude'
   READ(mag_field,*,ERR=10) obs%mag
   obs%mag_def=.true.
   pos_point=INDEX(mag_field,'.')
   IF(pos_point == 0) THEN
      obs%acc_mag=1
   ELSE
      n_dec_mag=len_mag_field-pos_point
      obs%acc_mag=10.0d0**(-n_dec_mag)
   END IF
END IF
! Astrometric catalog code from MPC
obs%catcodmpc(2:2)=mpcrec(72:72)

! discovery records to be removed
!IF (obs%tech.eq.'X') THEN
!   error_code=' corrected discovery obs. removed'
!   GOTO 10
!ENDIF
! Normal termination
error=.false.
RETURN

! Error termination
10 CONTINUE
len_error_code=lench(error_code)
WRITE(ierrou,101)mpcrec(1:12),error_code(1:len_error_code)
101 FORMAT(A12,'  ERROR in mpcrec_transform: ',A)
WRITE(ierrou,'(A80)')mpcrec
END SUBROUTINE mpcrec_transform

! Add the geocentric position of the observer to an astrometric observation
! In case of a satellite or "roving observatory" observation, the two MPC records must be passed
SUBROUTINE mpcrec_add_obspos(obs,error,mpcrec1,mpcrec2)
USE station_coordinates
USE reference_systems

TYPE(ast_obs),      INTENT(INOUT)          :: obs        ! astrometric observation to be completed
LOGICAL,            INTENT(OUT)            :: error      ! error conversion flag
CHARACTER(LEN=*),   INTENT(IN),   OPTIONAL :: mpcrec1    ! first  MPC record
CHARACTER(LEN=*),   INTENT(IN),   OPTIONAL :: mpcrec2    ! second MPC record

CHARACTER(LEN=80)                :: error_code             ! internal error code
CHARACTER*3                      :: obscods
INTEGER                          :: len_error_code,sign,obscod,le
DOUBLE PRECISION                 :: longitude,latitude,altitude,conv
DOUBLE PRECISION, DIMENSION(3)   :: bfpos0,equpos
! DOUBLE PRECISION, DIMENSION(3,3) :: rot
CHARACTER*5 :: altstr

! astronomical unit au in km from fund_const

! Obsolete CALLs (to be substituted by modules)
INTEGER, EXTERNAL :: lench
Character*9 name ! for message

error=.false.
name=obs%objdes

SELECT CASE (obs%type)
   CASE ('O')
      IF(obs%obscod_s == '247') THEN
! Roving observatory
         IF(.NOT.PRESENT(mpcrec1)) STOP '**** mpcrec_add_obspos: internal error (01) ****'
         IF(.NOT.PRESENT(mpcrec2)) STOP '**** mpcrec_add_obspos: internal error (02) ****'
         IF(mpcrec1(1:14) /= mpcrec2(1:14)) THEN
            error_code='mismatch in col 1-14'
            GOTO 10
         END IF
         IF(mpcrec2(15:15) /= 'v') THEN
            error_code='missing code "v" in column 15'
            GOTO 10
         END IF
         IF(mpcrec1(16:32) /= mpcrec2(16:32)) THEN
            error_code='mismatch in col 16-32'
            GOTO 10
         END IF
         IF(mpcrec2(33:34) /= '1 ') THEN
            error_code='abnormal content of columns 33-34'
            GOTO 10
         END IF
         IF(mpcrec2(78:80) /= '247') THEN
            error_code='abnormal content of columns 78-80'
            GOTO 10
         END IF
! Input and transformation of coordinates
         error_code='input field: longitude'
         READ(mpcrec2(35:44),*,ERR=10) longitude
         SELECT CASE (mpcrec2(46:46))
            CASE (' ', '+')
               sign=+1
            CASE ('-')
               sign=-1
            CASE DEFAULT
               error_code='input field: sign of latitude'
               GOTO 10
         END SELECT
         error_code='input field: latitude'
         READ(mpcrec2(47:55),*,ERR=10) latitude
         latitude=sign*latitude
         error_code='input field: altitude'
         altstr=mpcrec2(57:61)
         CALL rmsp(altstr,le)
         IF(le.gt.0)THEN
            READ(mpcrec2(57:61),*,ERR=10) altitude
         ELSE
! altitude error case; does not matter
            WRITE(ierrou,102)mpcrec2(57:61),obs%objdes
102         FORMAT('mpcrec_add_obspos: altitude is ',a5,' for ',a9)
            numerr=numerr+1
            altitude=0.d0
         ENDIF         
         longitude=longitude*radeg
         latitude=latitude*radeg
         ! Store information on the geodetic position of the observatory (needed only for output)
         IF(ASSOCIATED(obs%geopos)) DEALLOCATE(obs%geopos)
         ALLOCATE(obs%geopos(3))
         obs%geopos(1)=longitude
         obs%geopos(2)=latitude
         obs%geopos(3)=altitude
         CALL geodetic_to_cartesian(longitude,latitude,altitude,bfpos0)
         bfpos0=bfpos0/(aukm*1d3)
         CALL observer_position(obs%time_tdt,obs%obspos,obs%obsvel,BFPOS=bfpos0)
      ELSE
! Fixed astronomical observatory
         CALL observer_position(obs%time_tdt,obs%obspos,obs%obsvel,OBSCODE=obs%obscod_i)
      END IF

   CASE ('S')
! Satellite observation
      IF(.NOT.PRESENT(mpcrec1)) STOP '**** mpcrec_add_obspos: internal error (03) ****'
      IF(.NOT.PRESENT(mpcrec2)) STOP '**** mpcrec_add_obspos: internal error (04) ****'
      IF(mpcrec1(1:12) /= mpcrec2(1:12)) THEN
         error_code='mismatch in col 1-12 '//name
         GOTO 10
      END IF
      IF(mpcrec2(15:15) /= 's') THEN
         error_code='missing code "s" in column 15 '//name
         GOTO 10
      END IF
      IF(mpcrec1(16:32) /= mpcrec2(16:32)) THEN
         error_code='mismatch in col 16-32 '//name
         GOTO 10
      END IF
! Conversion of unit of length
      SELECT CASE (mpcrec2(33:33))
      CASE ('1')
         conv=1.d0/aukm
      CASE ('2')
         conv=1.d0
      CASE DEFAULT
         error_code='input field: unit of parallax indicator '//name
         GOTO 10
      END SELECT
      obs%par_unit=mpcrec2(33:33)
      error_code='input field: parallax vector'
      READ(mpcrec2(36:45),*,ERR=10) equpos(1)
      IF(mpcrec2(35:35) == '-') equpos(1)=-equpos(1)
      READ(mpcrec2(48:57),*,ERR=10) equpos(2)
      IF(mpcrec2(47:47) == '-') equpos(2)=-equpos(2)
      READ(mpcrec2(60:69),*,ERR=10) equpos(3)
      IF(mpcrec2(59:59) == '-') equpos(3)=-equpos(3)
! Transformation of equatorial coordinates into ecliptical
!      CALL rotpn(rot,'EQUM','J2000',0.d0,'ECLM','J2000',0.d0)
      obs%obspos=MATMUL(roteqec,equpos)
      obs%obspos=conv*obs%obspos
      obs%obsvel=0.d0
      error_code='input field: observatory code '//name
      READ(mpcrec1(78:80),*,ERR=10) obscods
      IF(obscods /= obs%obscod_s) THEN
         error_code='abnormal content of columns 78-80 '//name
         GOTO 10
      END IF

   CASE ('R', 'V')
! Radar observations
      WRITE(*,*)' mpcrec_add_obspos: MPC radar not alllowed '
      STOP
      CALL observer_position(obs%time_tdt,obs%obspos,obs%obsvel,OBSCODE=obs%obscod_i,PRECISION=2)
   CASE DEFAULT
! Unknown type of observation
      error_code='observation type: ' // obs%type
      GOTO 10

END SELECT

RETURN

10 CONTINUE
len_error_code=lench(error_code)
WRITE(ierrou,101) error_code(1:len_error_code)
101 FORMAT('*** ERROR in mpcrec_add_obspos: ',A)
numerr=numerr+1
error=.true.
END SUBROUTINE mpcrec_add_obspos
! END SUBMODULE iompc

! SUBMODULE ioradar
SUBROUTINE  jpl_radarobs(obs_found,file,obs,nobs)

CHARACTER*(*), INTENT(IN) :: file
LOGICAL, INTENT(OUT) :: obs_found ! true if there is something
INTEGER, INTENT(OUT) :: nobs ! number of observations
TYPE(ast_obs), INTENT(OUT), DIMENSION(:) :: obs ! vector of observations
! end interface
LOGICAL error
CHARACTER*100 rec
INTEGER i,l,unit,nlef
external lench
integer lench
  nlef=SIZE(obs)
  call filopn(unit,file,'old')
  nobs=0
  DO  i=1,nlef
     read(unit,101,end=10) rec
101  format(a)
     l=lench(rec)
     IF(l.eq.0) CYCLE
     nobs=nobs+1
     radar_rec(i)=rec
! parse one jpl radar record
     CALL jplradar_transform(rec,obs(i),error)
     IF(error)THEN
        l=lench(file)
        write(*,102) file(1:l),nobs
102     format(' **** jpl_radarobs: input conversion error ****'/          &
             &        ' **** file = "',a,'", line',i5,' ****')
        nobs=0
        obs_found=.false.
        call filclo(unit,' ')
        RETURN
     ENDIF
 ENDDO
! too many records
   STOP ' **** jpl_radarobs: nobs > nlef ****'

! regular ending
10 continue
   call filclo(unit,' ')
   IF(nobs.eq.0)THEN
      obs_found=.false.
      write(*,*) 'Warning: Found ',nobs,' obs in ',file
   ELSE
      obs_found=.true.
      IF(verb_io.gt.9)write(*,*) 'Found ',nobs,' obs in ',file
   ENDIF

RETURN
END SUBROUTINE jpl_radarobs

! Revised by Genny on 5 November 2005 to deal with the new format of 
! radar_data_ast.txt to include numbered objects > 99999
SUBROUTINE jplradar_transform(rec,obs,error)
USE station_coordinates

CHARACTER*(*), INTENT(IN)  :: rec
TYPE(ast_obs), INTENT(OUT) :: obs
LOGICAL,       INTENT(OUT) :: error

      INTEGER errcod
! name
      character*6 number
      character*17 nametmp,stname
      integer lnum,lnam
! time
      integer year,month,day,hour,min,isec,mjd,mjdt
      double precision sec,sect
      character*3 scale
! funcs
      INTEGER ix
      DOUBLE PRECISION tjm1
      EXTERNAL tjm1
! measurements
      double precision observ,rms,freq
      DOUBLE PRECISION r,v, accr,accv
      character*2 unit
      character*3 surf
      logical range
! obscode
      character*(name_len) trxstr,recstr
      integer iotr,iore
      CHARACTER(LEN=3) :: code_trx,code_rcx

      error=.true.
      errcod=1
      obs=undefined_ast_obs
! ========== HANDLE OBJECT NAME =================
      READ(rec,101,ERR=10) number,nametmp
  101 FORMAT(a6,1x,a16)
      call rmsp(number,lnum)
      call rmsp(nametmp,lnam)
      if(lnam.gt.9)lnam=9
      if(lnum.eq.0)then
! remove identified object if present ('1991AQ=1994RD')
         ix=index(nametmp,'=')
         if(ix.gt.0)lnam=ix-1
         obs%objdes=nametmp(1:lnam)
      else
         obs%objdes=number(1:lnum)
      endif
! ========== HANDLE DATE AND TIME =================
      READ(rec,102,ERR=10) year,month,day,hour,min,isec
  102 FORMAT(24X,I4,5(1X,I2))
      IF(year.LT.1972) THEN
          scale='UT1'
      ELSE
          scale='UTC'
      END IF
      sec=(hour*60d0+min)*60d0+isec
      mjd=nint(tjm1(day,month,year,0.d0))
      CALL cnvtim(mjd,sec,scale,mjdt,sect,'TDT')
      obs%time_utc=mjd+sec/86400.d0
      obs%time_tdt=mjdt+sect/86400.d0
      obs%acc_time=10d-10
! ========== READ MEASUREMENT =================
      READ(rec,103,ERR=10) observ,rms,unit,surf,freq
  103 FORMAT(44x,f13.2,1x,f7.3,1x,a2,1x,a3,1x,f5.0)
      freq=freq*1d6
      if(unit.eq.'us')then
         obs%type='R'
         range=.true.
      elseif(unit.eq.'Hz')then
         obs%type='V'
         range=.false.
      else
         errcod=37
         WRITE(*,*)rec
         goto 10
      endif
      if(surf.eq.'COM')then
         obs%tech='c'
      elseif(surf.eq.'PP ')then
         obs%tech='s'
      else
         errcod=47
         goto 10
      endif
      if(range)then
         r=observ*1d-6/86400d0*vlight/2.d0
         accr=rms*1d-6/86400d0*vlight/2.d0
         obs%coord(1)=r
         obs%acc_coord(1)=accr
         obs%coord(2)=0.d0
         obs%acc_coord(2)=-1.d0
      else
         v=-observ*vlight/(freq*2.d0)
         accv=rms*vlight/(freq*2.d0)
         obs%coord(2)=v
         obs%acc_coord(2)=accv
         obs%coord(1)=0.d0
         obs%acc_coord(1)=-1.d0
      endif
! ======================= HANDLE OBSERVATORY CODES ====================
      READ(rec,104,ERR=10) trxstr,recstr
  104 FORMAT(79x,a9,1x,a9)
      errcod=5
      IF(trxstr.eq.'-73')THEN
         trxstr='EISCAT'
         recstr='EISCAT'
      END IF
      iotr=radar_station(trxstr)
      iore=radar_station (recstr)
      obs%obscod_i=iotr*10000+iore
      CALL codestat(iotr,code_trx)
      CALL codestat(iore,code_rcx)
      obs%obscod_s=code_trx//' '//code_rcx
! obs%obspos = body-fixed position of transmitting station
! obs%obsvel = body-fixed position of receiving station
      CALL obscoo(iotr,obs%obspos,stname)
      CALL obscoo(iore,obs%obsvel,stname)
! ======================= FILL DUMMY FIELDS ===========================
      obs%note=' '
      obs%mag=0.d0
      obs%mag_str='      '
      obs%mag_band=' '
      obs%acc_mag=0.d0
      obs%mag_def=.FALSE.
! ====================== CLEAN UP =====================================
      error=.false.
   10 CONTINUE
      IF(error) THEN
         WRITE(*,105) errcod,rec
  105    FORMAT(' jplrad: error code',I3,'rec:',A)
      ENDIF
      RETURN

END SUBROUTINE jplradar_transform

! ====================== station func =============================
! compute observatory code from string
INTEGER FUNCTION radar_station(stastr)
  character*(*) stastr
! We have the following stations:
!     'Arecibo  ' = 251
!     'DSS 13   ' = 252
!     'DSS 14   ' = 253
!     'DSS 25   ' = 257
!     'Haystack ' = 254
!     'Evpatoria' = 255 (not in MPC file...)
!     'Greenbank' = 256
!     'EISCAT'    = 273 (not in MPC file...)
  if(stastr.eq.'Arecibo  ')then
     radar_station=251
  elseif(stastr.eq.'DSS 13   ')then
     radar_station=252
  elseif(stastr.eq.'DSS 14   ')then
     radar_station=253
  elseif(stastr.eq.'DSS 25   ')then
     radar_station=257
  elseif(stastr.eq.'Haystack ')then
     radar_station=254
  elseif(stastr.eq.'Evpatoria')then
     radar_station=255
  elseif(stastr.eq.'Greenbank')then
     radar_station=256
  elseif(stastr.eq.'EISCAT   ')then
!     radar_station=273
     radar_station=259
  else
     radar_station=500
     write(*,*)'jplrad: unknown station: ',stastr
     STOP '***unknown radar station****'
  endif
END FUNCTION radar_station
! END SUBMODULE ioradar
! Input of observations: reads sequentially the .rwo, .obs, .rad file and combines
! the data according to the precedence rule specified by precob

SUBROUTINE input_obs(obsdir,astna0,precob,error_model,obs0,obs,obsw,m,iun20,change,rms,rmsmag)

CHARACTER(LEN=*),             INTENT(IN)  :: obsdir      ! input dir. where files astna0.rwo, astna0.obs and astna0.rad are
CHARACTER(LEN=*),             INTENT(IN)  :: astna0      ! asteroid name
INTEGER,                      INTENT(IN)  :: iun20       ! messages unit
LOGICAL,                      INTENT(IN)  :: precob      ! .TRUE. for overwriting .rwo, .FALSE. for updating .rwo
CHARACTER(LEN=*),             INTENT(IN)  :: error_model ! error model file name
LOGICAL,                      INTENT(OUT) :: obs0        ! successful input
LOGICAL,                      INTENT(OUT) :: change      ! existence of new data
INTEGER,                      INTENT(OUT) :: m           ! number of observations
TYPE(ast_obs),  DIMENSION(:), INTENT(OUT) :: obs         ! observations
TYPE(ast_wbsr), DIMENSION(:), INTENT(OUT) :: obsw        ! weights and possibly residuals
DOUBLE PRECISION,             INTENT(OUT), OPTIONAL :: rms    ! RMS of astrometric fit (-1 = undefined)
DOUBLE PRECISION,             INTENT(OUT), OPTIONAL :: rmsmag ! RMS of photometric fit (-1 = undefined)

! observation numbers: maximum, space left
INCLUDE 'parobx.h90'
INTEGER :: nlef
! file names
CHARACTER*77 file
INTEGER lfile
LOGICAL rwo,mpc,rad  ! logical flags for existence of input files
! new obs. number
INTEGER mnew,mr,nlefm
! ===== observational data: temporary copy===========
INTEGER mt ! observation number
TYPE(ast_obs),  DIMENSION(nobx) :: obst ! observations
TYPE(ast_wbsr), DIMENSION(nobx) :: obswt ! weights/residuals
CHARACTER*(20) :: error_modelt ! error model file name
LOGICAL init ! for observ_rms: if not initialized, obsw should be set to default 
! ===============================================
! sorting
INTEGER iperm(nobx)
! directory char
INCLUDE 'sysdep.h90'
! loop index
INTEGER i
INTEGER ld
INTEGER, EXTERNAL :: lench
LOGICAL changemod ! error model changed with respect to previous .rwo file
! =============EXECUTION BEGINS======================
IF(PRESENT(rms)) rms=-1.d0
IF(PRESENT(rmsmag)) rmsmag=-1.d0
changemod=.false.
!  compute file name
ld=lench(obsdir)
IF(ld.GT.0) THEN
   IF(obsdir(ld:ld).EQ.dircha) THEN
      file=obsdir(1:ld)//astna0
   ELSE
      file=obsdir(1:ld)//dircha//astna0
   END IF
ELSE
   file=astna0
END IF
CALL rmsp(file,lfile)
! existence of .rwo, .obs, .rad
INQUIRE(file=file(1:lfile)//'.rwo',exist=rwo)
INQUIRE(file=file(1:lfile)//'.obs',exist=mpc)
INQUIRE(file=file(1:lfile)//'.rad',exist=rad)
! For now we will not do radar only orbits
IF(.not.rwo .and. .not. mpc)THEN
   WRITE(*,*)'You must provide either a .obs or .rwo file in directory ',obsdir
   obs0=.false.
   RETURN
ENDIF
! check sizes of arrays
nlef=SIZE(obs)
IF(nlef.ne.SIZE(obsw))THEN
   WRITE(*,*)' input_obs: problems in dimensions of arrays obs, obsw ',nlef,SIZE(obsw)
   STOP
ENDIF
m=0
! select operations mode
IF(.not.rwo)THEN
   IF(verb_io.gt.9)WRITE(*,*) 'No .rwo file, reading .obs and/or .rad files.'
! there is no .rwo, so read .obs and/or .rad
   IF(mpc)THEN
! Input of astrometric observations from a file (MPC format)
      CALL mpc_obs_input(mpc,obs,m,FILNAM=file(1:lfile)//'.obs',ONS_MODE=ons_name)
      IF(verb_io.gt.9)WRITE(*,*)'input_obs: ',m,' obs in ',file(1:lfile)//'.obs'
      WRITE(iun20,*)'input_obs: ',m,' obs in ',file(1:lfile)//'.obs'
   ENDIF
! read radar jpl data
   IF(rad)THEN
      nlefm=nlef-m
      CALL jpl_radarobs(rad,file(1:lfile)//'.rad',obs(m+1:),mr)
      IF(verb_io.gt.9)WRITE(*,*)'input_obs:',mr,' obs from  ',file(1:lfile)//'.rad'
      WRITE(iun20,*)'input_obs:',mr,' obs  ',file(1:lfile)//'.rad'
      m=mr+m
   ENDIF
   obs0=mpc.or.rad
   IF(.not.obs0)RETURN
! find weights for these; a priori RMS of astrometric observations
   init=.true.
   CALL observ_rms(obs,error_model,init,obsw,m)
! output data for possible manual fixing: create weights file; no rms available
   CALL write_rwo(file(1:lfile)//'.rwo',obs,obsw,m,error_model)
   change=.true.
! select between update and overwrite of .rwo
ELSEIF(precob)THEN
   IF(verb_io.gt.9)WRITE(*,*) file(1:lfile),'.rwo found but ALL obs will come from .obs/.rad files.'
! give the precedence to the observation files .obs and .rad
! with respect to .rwo, which is overwritten
   IF(mpc)THEN
! Input of astrometric observations from a file (MPC format)
      CALL mpc_obs_input(mpc,obs,m,FILNAM=file(1:lfile)//'.obs',ONS_MODE=ons_name)
      IF(verb_io.gt.9)WRITE(*,*)'input_obs: ',m,' obs in ',file(1:lfile)//'.obs'
      WRITE(iun20,*)'input_obs: ',m,' obs in ',file(1:lfile)//'.obs'
   ELSE
      m=0
   ENDIF
! read radar jpl data
   IF(rad)THEN
      nlefm=nlef-m
      CALL jpl_radarobs(rad,file(1:lfile)//'.rad',obs(m+1:),mr)
      IF(verb_io.gt.9)WRITE(*,*)'input_obs:',mr,' obs from  ',file(1:lfile)//'.rad'
      WRITE(iun20,*)'input_obs:',mr,' obs  ',file(1:lfile)//'.rad'
      m=mr+m
   ENDIF
   obs0=rad.or.mpc
! If no obs then object does not "exist", so rwo should not be read:
   IF(.not.obs0)return
! find weights for these; a priori RMS of astrometric observations
   init=.true.
   CALL observ_rms(obs,error_model,init,obsw,m)
! give default selection flag of 1
   obsw(1:m)%sel_coord=1
! format transition:
!  CALL read_rwo_trans(file(1:lfile)//'.rwo',obst,obswt,mt,error_modelt,rms,rmsmag)
! read .rwo  anyway, but store in temporary array the data
   CALL read_rwo(file(1:lfile)//'.rwo',obst,obswt,mt,error_modelt,rms,rmsmag)
   IF(error_model.ne.error_modelt)THEN
      init=.false. ! structure is initialized already
      CALL observ_rms(obst,error_model,init,obswt,mt)
      changemod=.true.
   ELSE
      changemod=.false.
   ENDIF
! recover informations from .rwo (only selection flags)
   CALL addobs_rwo(obs,obsw,m,obst,obswt,mt,change)
! if error model is changed, the data have to be considered changed
   change=(change.or.changemod) 
   IF(change)THEN
      IF(verb_io.gt.9)WRITE(*,*)'There are new/changed obs. New numobs=',m
! output updated .rwo file, if there are new observations (erasing residuals)
      CALL write_rwo(file(1:lfile)//'.rwo',obs,obsw,m,error_model)
   ELSE
      IF(verb_io.gt.9)WRITE(*,*)'There are no updates in .obs or .rad files.'
   ENDIF
ELSE
! give the precedence to .rwo, the .obs and .rad files are intended
! as additional observations only; data in .rwo are not erased, can only be
! changed
! read .rwo, and store data in final array
   IF(verb_io.gt.9)WRITE(*,*)'Using .rwo file, but checking .obs,.rad for update.'
   CALL read_rwo(file(1:lfile)//'.rwo',obs,obsw,m,error_modelt,rms,rmsmag)
   IF(verb_io.gt.9)WRITE(*,*)'input_obs: ',m,' obs from  ',file(1:lfile)//'.rwo'
   WRITE(iun20,*)'input_obs: ',m,' obs from ',file(1:lfile)//'.rwo'
   IF(m.eq.0)THEN
      obs0=.false.
   ELSE
      obs0=.true.
   ENDIF
   IF(error_model.ne.error_modelt)THEN
      init=.false. ! structure is initialized already
      CALL observ_rms(obs,error_model,init,obsw,m)
      changemod=.true.
   ELSE
      changemod=.false.
   ENDIF
! if there are input data
   IF(mpc)THEN
! Input of astrometric observations into temporary from a file (MPC form
      CALL mpc_obs_input(mpc,obst,mt,FILNAM=file(1:lfile)//'.obs',ONS_MODE=ons_name)
!      IF(error_model.EQ.'cbm10') THEN
!         obs%catcodmpc=obst%catcodmpc
!      ENDIF
      IF(.not.mpc)THEN
         WRITE(*,*) file(1:lfile)//'.obs is possibly corrupt. Not using any data from this file.'
      ELSE
         IF(verb_io.gt.9)WRITE(*,*)'input_obs:',mt,' obs from  ',file(1:lfile)//'.obs'
         WRITE(iun20,*)'input_obs:',mt,' from ',file(1:lfile)//'.obs'
      ENDIF
   ENDIF
! read radar jpl data
   IF(rad)THEN
      nlefm=nobx-mt
      CALL jpl_radarobs(rad,file(1:lfile)//'.rad',obst(mt+1:),mr)
      IF(.not.rad)THEN
         WRITE(*,*) file(1:lfile)//'.rad is possibly corrupt. Not using any data from this file.'
      ELSE
         IF(verb_io.gt.9)WRITE(*,*)'input_obs:',mr,' obs from  ',file(1:lfile)//'.rad'
         WRITE(iun20,*)'input_obs:',mr,' obs  ',file(1:lfile)//'.rad'
         mt=mr+mt
      ENDIF
   ENDIF
! add information from .obs and .rad
   IF(mpc.or.rad)THEN
      obs0=.true.
! find weights for these; a priori RMS of astrometric observations
      init=.true.
      CALL observ_rms(obst,error_model,init,obswt,mt)
      CALL addobs_mpc(obs,obsw,m,obst,obswt,mt,mnew,change)
! if error model is changed, the data have to be considered changed
      change=(change.or.changemod)   
      IF(change)THEN
         IF(verb_io.gt.9)WRITE(*,*)'There are new/changed obs. New numobs=',mnew
! output updated .rwo file, if there are new observations (erasing residuals)
         CALL write_rwo(file(1:lfile)//'.rwo',obs,obsw,mnew,error_model)
         m=mnew
      ELSE
         IF(verb_io.gt.9)WRITE(*,*)'There are no updates in .obs or .rad files.'
      ENDIF
   else
      change=.true.
! check for new weights anyway???
   ENDIF
ENDIF

! get asteroid radius (if necessary) before returning
CALL aster_radius(obs%objdes,obs%type,m)
! =======  sort data before returning =========
CALL heapsort(obs%time_tdt,m,iperm)
! copy output into temp vectors
obst(1:m)=obs(1:m)
obswt(1:m)=obsw(1:m)
! copy back input output vectors in sorted order
Do i=1,m
   obs(i)=obst(iperm(i))
   obsw(i)=obswt(iperm(i))
ENDDO

END SUBROUTINE input_obs

! SUBMODULE addobs
! =========================================
!  A D D O B S _ M P C
!
! add information from .obs and .rad to the one available from a file .rwo
! The number of observations can increase, from m to a maximum
! which is m+mt; however, the vectors have a maximum dimension nlef
SUBROUTINE addobs_mpc(obs,obsw,m,obst,obswt,mt,mnew,change)
! ===============================================
! change flag, new obs. number
  LOGICAL, INTENT(OUT) :: change
  INTEGER, INTENT(OUT) ::  mnew
! ===== observational data ===========================
  INTEGER m ! observation number
  TYPE(ast_obs), INTENT(INOUT), DIMENSION(:) :: obs ! observations
  TYPE(ast_wbsr), INTENT(INOUT), DIMENSION(:) :: obsw ! weights and possibly residuals
! ===== observational data: temporary copy===========
  INTEGER mt ! observation number
  TYPE(ast_obs), INTENT(IN), DIMENSION(:) :: obst ! observations
  TYPE(ast_wbsr), INTENT(INOUT), DIMENSION(:) :: obswt ! weights and possibly residuals
! ===========================================
  INTEGER nlef ! length of obs,obsw vectors
  INTEGER j,mj,double
! ==================================================
  nlef=SIZE(obs)
  IF(nlef.ne.SIZE(obsw))THEN
     WRITE(*,*)' addobs_mpc: problems in dimensions of arrays obs, obsw ',nlef,SIZE(obsw)
     STOP
  ENDIF
! monitor changes
  mnew=m
  change=.false.
! scan supposedly new observations
  DO 1 j=1,mt
     mj=find_obs(obst(j),obs,m,double)
     IF(mj.ne.0.and.double.eq.0)THEN
!     the observation was already there
        IF(changed_obs(obst(j),obswt(j),obs(mj),obsw(mj)))THEN
! ... but it is changed
           change=.true.
           obs(mj)=obst(j)
           obsw(mj)=obswt(j)
        ELSEIF(change_format.EQ.'INPUT') THEN
           change=.true.
           obs(mj)%catcodmpc=obst(j)%catcodmpc
           obsw(mj)%bias_coord=obswt(j)%bias_coord
           IF(.not.obsw(mj)%force_w(1)) THEN
              obsw(mj)%rms_coord(1)=obswt(j)%rms_coord(1)
           ENDIF
           IF(.not.obsw(mj)%force_w(2)) THEN
              obsw(mj)%rms_coord(2)=obswt(j)%rms_coord(2)
           ENDIF
        ELSEIF(obs(mj)%catcodmpc.NE.obst(j)%catcodmpc) THEN
           change=.true.
           obs(mj)%catcodmpc=obst(j)%catcodmpc
           obsw(mj)%bias_coord=obswt(j)%bias_coord
        ENDIF
     ELSEIF(mj.ne.0.and.double.ne.0)THEN
! the observation was already there, in double copy!
        IF(changed_obs(obst(j),obswt(j),obs(mj),obsw(mj)))THEN
           IF(changed_obs(obst(j),obswt(j),obs(double),obsw(double)))THEN
              change=.true.
! double, and changed! human intervention required
              WRITE(*,*)'addobs_mpc: double and changed'
              WRITE(*,*)' records ',mj,' and ',double,' in .rwo'
              WRITE(*,*)' record ',j,' in .obs'
!             STOP
           ELSE
! OK, it is the double
           ENDIF
        ELSE
! OK, it is the first one
        ENDIF
     ELSEIF(mj.eq.0)THEN
! the observation is new: add it
        change=.true.
!       WRITE(*,*)'addobs: new observation at record ',j
        mnew=mnew+1
        obs(mnew)=obst(j)
        obsw(mnew)=obswt(j)
     ENDIF
1 ENDDO
END SUBROUTINE addobs_mpc
! =======================================================
! F I N D _ O B S
! find an observation from a list, matching the given one
INTEGER FUNCTION find_obs(obst,obs,m,double)
! ============= INPUT =====================
! number of obs. record to be scanned, time, time to be found
  INTEGER m ! number of obs
  TYPE(ast_obs), INTENT(IN) :: obst ! one new observation
  TYPE(ast_obs),  INTENT(IN), DIMENSION(m) :: obs ! observations
! ============OUTPUT (impure func!) ===========
  INTEGER, INTENT(OUT) :: double ! location of the duplicate observation
! =========END INTERFACE=====================
  INTEGER j
  DOUBLE PRECISION, PARAMETER ::  epst= 1.d-8 ! time control
  find_obs=0
  double=0
  DO 1 j=1,m
     IF(abs(obst%time_utc-obs(j)%time_utc).lt.epst.and.                                &
          &      obst%obscod_i.eq.obs(j)%obscod_i.and.obst%type.eq.obs(j)%type)THEN
        ! observation found; is it a double?
        IF(find_obs.ne.0)THEN
           IF(double.eq.0)THEN
              double=j
              IF(ierrou.gt.0)THEN
                 IF(obst%tech.eq.'X'.or.obs(find_obs)%tech.eq.'X'.or.obs(j)%tech.eq.'X')THEN
!                   X obs. is often a duplicate, no problem
                 ELSE
                    WRITE(ierrou,*)'findob: two same time',find_obs,j,obst%time_utc,obst%obscod_s,obst%tech 
                    numerr=numerr+1
                 ENDIF
              ELSE
                 IF(verb_io.gt.9)WRITE(*,*)'findob: two same time', find_obs,j,obst%time_utc, obst%obscod_s
              ENDIF
           ELSE
              IF(ierrou.gt.0)THEN ! disaster case: triple observation!!!
                 WRITE(ierrou,*)'findob: three same time', &
                      &                     find_obs,double, j,obst%time_utc, obst%obscod_s, obst%tech 
              ELSE
                 WRITE(*,*)'findob: three same time',                &
                      &                     find_obs,double, j,obst%time_utc, obst%obscod_s, obst%tech 
              ENDIF
              !             STOP
           ENDIF
        ELSE
           find_obs=j !normal case, no double
        ENDIF
     ENDIF
1 ENDDO
END FUNCTION find_obs
! =======================================================
! C H A N G E D _ O B S
! is an observation changed?
LOGICAL FUNCTION changed_obs(obst,obswt,obs,obsw)
  TYPE(ast_obs), INTENT(IN) :: obst ! one new observation
  TYPE(ast_wbsr), INTENT(IN) :: obswt ! one new observation weight
  TYPE(ast_obs),  INTENT(IN) :: obs ! one old observations
  TYPE(ast_wbsr), INTENT(IN) :: obsw ! one old observation weight
  DOUBLE PRECISION::  epsa, epsw, epsb ! control on angle, weight, bias
  changed_obs=.false.
  epsa=1d-8
  IF(abs(obs%coord(1)-obst%coord(1)).gt.epsa*abs(obs%coord(1)).or.             &
     &      abs(obs%coord(2)-obst%coord(2)).gt.epsa*abs(obs%coord(2)))changed_obs=.true.
! temporary fix
!      IF(obs%mag_str.ne.obst%mag_str.or.obs%tech.ne.obst%tech)changed_obs=.true.
! problem: if the weights are changed, should this be considered a changed observation??
! obsw, obswt are included to make this possible: but no
!       epsw=1.d-2
!       IF(abs(obsw%rms_coord(1)-obswt%rms_coord(1)).gt.epsw*abs(obsw%rms_coord(1)).or.      &
!      &   abs(obsw%rms_coord(2)-obswt%rms_coord(2)).gt.epsw*abs(obsw%rms_coord(2)))       &
!      & changed_obs=.true. 
!       epsb=1.d-8 ! control on bias is absolute, because bias can be small!!
!      IF(abs(obsw%bias_coord(1)-obswt%bias_coord(1)).gt.epsb.or.     &
!      &   abs(obsw%bias_coord(2)-obswt%bias_coord(2)).gt.epsb) changed_obs=.true. 
END FUNCTION changed_obs
! =========================================
!  A D D O B S _ R W O
!
! add information from .rwo file (only selection flags and manually fixed weights)
! to the observations from .obs and .rad
! The number of observations cannot increase, remains m
SUBROUTINE addobs_rwo(obs,obsw,m,obst,obswt,mt,change)
  USE output_control
! logical change flag
  LOGICAl change
! ===== observational data ===========================
  INTEGER m ! observation number
  TYPE(ast_obs), INTENT(INOUT), DIMENSION(:) :: obs ! observations
  TYPE(ast_wbsr), INTENT(INOUT), DIMENSION(:) :: obsw ! weights and possibly residuals
! ===== observational data: temporary copy===========
  INTEGER mt ! observation number
  TYPE(ast_obs), INTENT(IN), DIMENSION(:) :: obst ! observations
  TYPE(ast_wbsr), INTENT(IN), DIMENSION(:) :: obswt ! weights and possibly residuals
! ===========================================
!     INTEGER nlef ! length of obs,obsw vectors
! ===========================================
  INTEGER j,mj,double
! if all the same...
  change=.false.
! scan old weigts and selection flags, see if they match  observations
  DO 1 j=1,mt
     mj=find_obs(obst(j),obs,m,double)
     IF(mj.ne.0.and.double.eq.0)THEN
! this observation is still present in .obs, .rad files
        IF(changed_obs(obst(j),obswt(j),obs(mj),obsw(mj)))THEN
! ... but it is changed, thus leave selection flag=1 and default weight
           change=.true.
        ELSE
! if no change to observation then preserve the
! manually fixed weights and selection flags
           IF(obswt(j)%force_w(1))THEN
              obsw(mj)%rms_coord(1)=obswt(j)%rms_coord(1)
              obsw(mj)%force_w(1)=obswt(j)%force_w(1)
           ENDIF
           IF(obswt(j)%force_w(2))THEN
              obsw(mj)%rms_coord(2)=obswt(j)%rms_coord(2)
              obsw(mj)%force_w(2)=obswt(j)%force_w(2)
           ENDIF
           IF(obswt(j)%rms_mag.lt.0.d0)obsw(mj)%rms_mag=obswt(j)%rms_mag
! selection flags are preserved anyway
           obsw(mj)%sel_coord= obswt(j)%sel_coord
           obsw(mj)%sel_mag= obswt(j)%sel_mag
        ENDIF
     ELSEIF(mj.ne.0.and.double.ne.0)THEN
! the observation was already there, in double copy!
        IF(changed_obs(obst(j),obswt(j),obs(mj),obsw(mj)))THEN
           IF(changed_obs(obst(j),obswt(j),obs(double),obsw(double)))THEN
              change=.true.
! double, and changed! human intervention required
              WRITE(*,*)'addobs_rwo: double and changed'
              WRITE(ierrou,*)'addobs_rwo: double and changed'
              numerr=numerr+1
              WRITE(*,*)' records ',mj,' and ',double,' in .obs'
              WRITE(*,*)' record ',j,' in .rwo'
              WRITE(ierrou,*)' records ',mj,' and ',double,' in .obs'
              WRITE(ierrou,*)' record ',j,' in .rwo'
!             STOP
           ELSE
! OK, it is the double
! it is the same, so preserve the selection flags
! manually fixed weights and selection flags
              IF(obswt(j)%force_w(1))THEN
                 obsw(double)%rms_coord(1)=obswt(j)%rms_coord(1)
                 obsw(double)%force_w(1)=obswt(j)%force_w(1)
              ENDIF
              IF(obswt(j)%force_w(2))THEN
                 obsw(double)%rms_coord(2)=obswt(j)%rms_coord(2)
                 obsw(double)%force_w(2)=obswt(j)%force_w(2)
              ENDIF
              IF(obswt(j)%rms_mag.lt.0.d0)obsw(mj)%rms_mag=obswt(j)%rms_mag
              obsw(double)%sel_coord=obswt(j)%sel_coord
              obsw(double)%sel_mag=obswt(j)%sel_mag
           ENDIF
        ELSE
! OK, it is the first one
! it is the same, so preserve the selection flags
! manually fixed weights and selection flags
           IF(obswt(j)%force_w(1))THEN
              obsw(mj)%rms_coord(1)=obswt(j)%rms_coord(1)
              obsw(mj)%force_w(1)=obswt(j)%force_w(1)
           ENDIF
           IF(obswt(j)%force_w(2))THEN
              obsw(mj)%rms_coord(2)=obswt(j)%rms_coord(2)
              obsw(mj)%force_w(2)=obswt(j)%force_w(2)
           ENDIF
           IF(obswt(j)%rms_mag.lt.0.d0)obsw(mj)%rms_mag=obswt(j)%rms_mag
           obsw(mj)%sel_coord=obswt(j)%sel_coord
           obsw(mj)%sel_mag=obswt(j)%sel_mag
        ENDIF
     ELSEIF(mj.eq.0)THEN
! if it is not found in .rwo, leave the default weights and selection flags
        change=.true.
     ENDIF
1 ENDDO
! check if there are extra (added) observations in .obs file
! it might be better to loop on the .obs rather than the .rwo data
  IF (mt.lt.m) change=.true.
END SUBROUTINE addobs_rwo

! Copyright Orbfit Consortium 1999
! this routine determines an asteroid radius if needed for radar
! Revised by Genny on 5 November 2005 to deal with the new format of 
! astorb.dat to include numbered objects > 99999
SUBROUTINE aster_radius(objid,obstype,m)
  INTEGER,INTENT(IN) :: m
  CHARACTER*(name_len), INTENT(IN), DIMENSION(m) ::  objid  ! designation
  CHARACTER*1, INTENT(iN), DIMENSION(m) ::  obstype ! find if radar
! OUTPUT: radius of asteroid is public radius
! ====================================
  integer unit,ln,lnam,lnum,isav,i
  character*18 oid,number,name
  double precision hmag,diam,exponent
  logical needed

! First find out if we need the radius at all
  needed=.false.
  do i=1,m
     if(obstype(i).eq.'R'.or.obstype(i).eq.'V')then
        needed=.true.
        isav=i
        exit
     endif
  enddo
  if(.not.needed)then
! bail out, leaving as it is.
!         write (*,*) 'Asteroid radius not required.'
!         radius=-1d0
     return
  else
! get radius
     call filopl(unit,'astorb.rad')
     oid=objid(isav)
     call rmsp(oid,ln)
! read records in an infinite loop (it's F(UGLY)77)
1    continue
     read(unit,101,end=109) number,name,hmag,diam
     call rmsp(number,lnum)
     call rmsp(name,lnam)
     if(oid(1:ln).eq.number(1:lnum) .or.                       &
          &           oid(1:ln).eq.name(1:lnam)    )THEN
!              found object in file
        if(diam.gt.0d0)then
           radius=diam/2d0/1.4998d8
           IF(verb_io.gt.9)write(*,102)'RADAR: Using radius = ',                 &
 &                 radius*1.4998d8,' km from IRAS diameter.'
        else
           exponent=0.5d0*(6.3d0-log10(0.2d0)-0.4d0*hmag)
           radius=10**exponent/2d0/1.4998d8
           IF(verb_io.gt.9)write(*,102)'RADAR: Using radius = ',                 &
                &                 radius*1.4998d8,' km from absolute magnitude.'
        endif
! exit infinite loop
        goto 109
     endif
     goto 1
  endif
  write(*,*)'*** astrad warning: ',name,' not found in astorb.dat.'
  radius=1d0/1.4998d8
  write(*,*)'RADAR: Using radius = ',radius*1.4998d8,' km.'
109 call filclo(unit,' ')
! ==============================================================
!                  num   nam    comp   Hmag  Gmag    col   diam
101 format(a6,1X,A18,1X,15x,1X,f5.2,1X,5x,1X,4x,1X,f5.1)
102 format(a,f6.2,a)
! ==============================================================
END SUBROUTINE aster_radius
!
! Copyright (C) 1997 by Mario Carpino
! Version: February 24, 1997
!
!  *****************************************************************
!  *                                                               *
!  *                         R D A N G A                           *
!  *                                                               *
!  *    Read an angle with its accuracy from a character string    *
!  *                                                               *
!  *****************************************************************
!
! INPUT:    STRING    -  Character string
!
! OUTPUT:   ANGLE     -  Angle (no unit conversion is performed)
!           ACC       -  Angle accuracy
!           ERROR     -  Conversion error
!
SUBROUTINE rdanga(string,angle,acc,error)
  CHARACTER*(*) string
  DOUBLE PRECISION angle,acc
  LOGICAL error

! Max string length
  INTEGER lx
  PARAMETER (lx=200)

  CHARACTER*(lx) c1,c,field
  INTEGER l,isig,nf,i,pp,iv,ll
  LOGICAL nospli
  DOUBLE PRECISION fact,rv


  INTEGER lench,nitchs
  EXTERNAL lench,nitchs

  error=.true.
  IF(lench(string).GT.lx) STOP '**** rdanga: LEN(string) > lx ****'
  c1=string
  CALL norstr(c1,l)

! The sign may be separated from the value by blank space
  isig=1
  IF(c1(1:1).EQ.'+') THEN
     c=c1(2:)
     CALL norstr(c,l)
  ELSEIF(c1(1:1).EQ.'-') THEN
     isig=-1
     c=c1(2:)
     CALL norstr(c,l)
  ELSE
     c=c1
  END IF

  nf=nitchs(c)
  IF(nf.LT.1.OR.nf.GT.3) RETURN
  angle=0
  fact=1
  DO 1 i=1,nf
     CALL stspli(c,' ',field,nospli)
     IF(nospli) RETURN
     pp=INDEX(field,'.')
     IF(pp.GT.0) THEN
        IF(i.NE.nf) RETURN
        READ(field,*,err=10,end=10) rv
     ELSE
        READ(field,*,err=10,end=10) iv
        rv=iv
     END IF
     angle=angle+fact*rv
     IF(i.EQ.nf) THEN
        CALL norstr(field,ll)
        pp=INDEX(field,'.')
        IF(pp.EQ.0) THEN
           acc=1
        ELSE
           ll=lench(field)
           acc=10.0d0**(pp-ll)
        END IF
        acc=acc*fact
     ELSE
        fact=fact/60.d0
     END IF
1 END DO
  angle=angle*isig
  error=.false.

10 CONTINUE
END SUBROUTINE rdanga

!SUB-MODULE observ_bias_weight
! Version 4.2 Copyright OrbFit consortium 2010
! ====================================================================
! ASTROW_BIAS_CBM
!     Computation of bias from catalogs, and 
!          Computation of a-priori RMS of observations 
!          from Chesley, Baer and Monet (2010) error model
! ====================================================================
!
! INPUT:    OBS          - observations
!           ERROR_MODEL  - name of weighing scheme to be used
!           N            - number of observations
!           INIT         - reinitialize all the structure (when init=.false.,
!                          information different from weight/bias is preserved)
!
! OUTPUT:   OBSW         - observation weights
!
SUBROUTINE observ_rms(obs,error_model,init,obsw,n)
 INTEGER,                      INTENT(IN)    ::  n  ! number of obs
 LOGICAL,                      INTENT(IN)    :: init! if true reset obsw 
 TYPE(ast_obs), DIMENSION(n),  INTENT(INOUT) :: obs ! observations
 CHARACTER*(*),                INTENT(IN)    :: error_model ! code for errmod
 TYPE(ast_wbsr), DIMENSION(n), INTENT(INOUT) :: obsw ! when init=.false.,
! =================END INTERFACE========================================
 INTEGER i,j, k ,jj
 DOUBLE PRECISION rmsrmi,rmsvmi, rmstmp(2), t2
 INTEGER unit ! for input of RMS values
 INTEGER ncat,lench
! Get default magnitude weights from a function magrms_new in same module
 CHARACTER q*4, obscod*3, catcod*1, qq*3
 CHARACTER*6 catcodstr
 INTEGER ind(n),nob(n), nob1(n) !sorting by time, number obs in batch
 TYPE(ast_obs), DIMENSION(n) :: obs1
 DOUBLE PRECISION, PARAMETER :: gapmax=8.d0/24.d0 ! 8 hours max gap 
! ======================================================================
 LOGICAL, SAVE :: first
 DATA first/.true./
! ==============initialize arrays of RMS rules==========================
 IF(first)THEN
    IF(error_model.eq.'fcct14')THEN
       CALL filopl(unit,'weights/fcct14.rules')
       nrules=0
       DO j=1, nrulesx
          READ(unit, 101, END=97)obscod,catcodstr,rmstmp
101       FORMAT(A3,3X,A6,2X,F5.2,2X,F5.2)
          ncat=lench(catcodstr)-2
          DO i=1,ncat
             obscatcod(nrules+i)=obscod//catcodstr(2+i:2+i)
             rmsbin(1:2,nrules+i)=rmstmp
          ENDDO          
          nrules=nrules+ncat
       ENDDO
97     CALL filclo(unit,' ')
       CALL heapsortnachlo(obscatcod,nrules,4,q,indrules)
    ELSEIF(error_model.eq.'cbm10')THEN
! rules for observatory, whatever star catalog
       CALL filopl(unit,'weights/cbm10.sta')
       READ(unit,*) ! skip header
       nstarms=0
       DO j=1, nstarmsx
          READ(unit, *, END=98)obscod,rmsstabin(1:2,j)
! example: 244        0.20       0.20
          obscodrms(j)=obscod
          nstarms=nstarms+1
       ENDDO
98     CALL filclo(unit,' ')
       CALL heapsortnachlo(obscodrms,nstarms,3,qq,indstarms)
! rules for observatory and star catalog
       CALL filopl(unit,'weights/cbm10.rules')
       READ(unit,*) ! skip header
       nrules=0
       DO j=1, nrulesx
          READ(unit, 100, END=99)obscod,catcod,rmsbin(1:2,j)
          obscatcod(j)=obscod//catcod
! example: 704:cc @ 1.23, 1.19
100       FORMAT(A3,2X,A1,2X,F5.2,1X,F5.2)
          nrules=nrules+1
       ENDDO
99     CALL filclo(unit,' ')
       CALL heapsortnachlo(obscatcod,nrules,4,q,indrules)
    ENDIF
    first=.false.
 ENDIF
! for each observations build the corresponding weights and bias record
 DO  i=1,n
    IF(init) obsw(i)=undefined_ast_wbsr
    IF(obs(i)%type.eq.'O'.and.(error_model.eq.'cbm10'.OR.error_model.eq.'fcct14'))THEN
       CALL astrow_bias_cbm(error_model,obs(i),obsw(i))
! magnitude weigthing
       obsw(i)%rms_mag=magrms_new(obs(i)%mag_str,obs(i)%time_tdt,obs(i)%obscod_i,obs(i)%tech)
    ELSEIF(obs(i)%type.eq.'O'.or.obs(i)%type.eq.'S')THEN
! astrometric (telescope) observations; get weight
       CALL astrow_nobias(error_model,obs(i),obsw(i))
! magnitude weigthing
       obsw(i)%rms_mag=magrms_new(obs(i)%mag_str,obs(i)%time_tdt,obs(i)%obscod_i,obs(i)%tech)
    ELSEIF(obs(i)%type.eq.'R')THEN
! weighting of radar data is given with observations, but there is
! minimum credible (for a given acccuracy of the models)
       rmsrmi=1.d-3/aukm
       obsw(i)%rms_coord(1)=MAX(obs(i)%acc_coord(1),rmsrmi)
       obsw(i)%rms_coord(2)=9.d9
       obsw(i)%rms_mag=9.d9
       obsw(i)%bias_coord(2)=0.d0
    ELSEIF(obs(i)%type.eq.'V')THEN
       rmsvmi=1.d-3/aukm
       obsw(i)%rms_coord(2)=MAX(obs(i)%acc_coord(2),rmsvmi)
       obsw(i)%rms_coord(1)=9.d9
       obsw(i)%rms_mag=9.d9
       obsw(i)%bias_coord(1)=0.d0
    ELSE
       WRITE(*,*)'observ_rms: obs. type', obs(i)%type,' not known ', &
     &           'at time ',obs(i)%time_tdt,' number ',i
       STOP
    ENDIF
    IF(obs(i)%tech.EQ.'x'.OR.obs(i)%tech.EQ.'X') THEN
       obsw(i)%sel_coord=0
    END IF
 END DO
! apply sqrt(nobs) factor within one batch (ending with an 8 hour gap).
! sort by time first
  CALL heapsort(obs(1:n)%time_tdt,n,ind)
  DO i=1,n
     obs1(i)=obs(ind(i))
  ENDDO
! find number of observations in batch from same station
  nob1=0  
  DO  j=1,n
    IF(obs1(j)%type.eq.'O'.and.error_model.eq.'fcct14')THEN
       IF(nob1(j).gt.0)CYCLE ! already done
       t2=obs1(j)%time_tdt
       nob1(j)=1
       IF(obs1(j)%tech.eq.'X'.or.obs1(j)%tech.eq.'x')CYCLE
! find batch from same station
       DO k=j+1,n
          IF(obs1(k)%time_tdt-t2.gt.gapmax) EXIT ! batch finished
          IF(obs1(k)%obscod_s.eq.obs1(j)%obscod_s)THEN
             nob1(j)=nob1(j)+1
             t2=obs1(k)%time_tdt
          ENDIF
       ENDDO
! assign number of obs in batch to each obs in batch
       DO k=j+1,n
          IF(obs1(k)%time_tdt-t2.gt.gapmax) EXIT ! batch finished
          IF(obs1(k)%obscod_s.eq.obs1(j)%obscod_s)THEN
             nob1(k)=nob1(j)
          ENDIF
       ENDDO      
    ENDIF
  ENDDO
!  DO i=1,n
!     WRITE(*,*) nob1(i)
!  END DO
!  STOP
! go back to original array of obs., apply sqrt
  DO k=1,n
    jj=ind(k)
    nob(jj)=nob1(k)
    IF(nob(jj).ne.0)THEN
       IF(.NOT.obsw(jj)%force_w(1))THEN
          obsw(jj)%rms_coord(1)=obsw(jj)%rms_coord(1)*sqrt(1.d0*nob(jj))
       END IF
       IF(.NOT.obsw(jj)%force_w(2))THEN
          obsw(jj)%rms_coord(2)=obsw(jj)%rms_coord(2)*sqrt(1.d0*nob(jj))
       END IF
    ENDIF
  ENDDO

END SUBROUTINE observ_rms
! ====================================================================
! ASTROW_BIAS_CBM
!     Computation of bias from catalogs, and 
!          Computation of a-priori RMS of observations 
!          from Chesley, Baer and Monet (2010) error model
!          also from Farnocchia, Chesley, Chamberlin, Tholen (2014)
! ====================================================================
!
! INPUT:    OBS          - observations
!           ERROR_MODEL  - name of weighing scheme to be used
!           N            - number of observations
!           INIT         - reinitialize all the structure (when init=.false.,
!                          information different from weight/bias is preserved)
!
! OUTPUT:   OBSW         - observation weights
! ====================================================================
SUBROUTINE astrow_bias_cbm(error_model,obs,obsw)
  USE cat_debias
  IMPLICIT NONE
  TYPE(ast_obs),INTENT(INOUT)    :: obs
  CHARACTER*(*),INTENT(IN)    :: error_model
  TYPE(ast_wbsr), INTENT(INOUT) :: obsw ! when init=.false.,
                   ! information different from weight/bias is preserved
! ===============end interface===============================
  CHARACTER*2 catcode
  LOGICAL biasflag
  DOUBLE PRECISION bias_ra,bias_dec,pmRA,pmDEC,rmsmin,tdt,rmsa,rmsd,deltat
  CHARACTER observ*3, obscat*4
  LOGICAL ermrms
  INTEGER indobs
  DOUBLE PRECISION, PARAMETER   :: mjd2000=51544.d0
! ============================================================
  IF(error_model.NE.'cbm10'.AND.error_model.NE.'fcct14')THEN
     WRITE(*,*)  '**** WRONG ERROR MODEL ***** ',error_model
     STOP
  ENDIF
  tdt=obs%time_utc
  observ=obs%obscod_s
! years since 2000
  deltat=(tdt-mjd2000)/365.25d0
! find bias
  CALL cat_bias(obs%coord(1),obs%coord(2),obs%catcodmpc(2:2),bias_ra,bias_dec,pmRA,pmDEC,biasflag,error_model)
! record bias in ast_wbsr data type
  IF(biasflag)THEN
! pmRA and pmDEC [mas/y], bias_ra, bias_dec [arcsec]
    obsw%bias_coord(1)=(bias_ra+pmRA*deltat/1000.d0)*radsec/cos(obs%coord(2)) !!!Checked with Steve
    obsw%bias_coord(2)=(bias_dec+pmDEC*deltat/1000.d0)*radsec
  ELSE
    obsw%bias_coord(1:2)=0.d0
  ENDIF
! new weighting for fcct14
  IF(error_model.eq.'fcct14')THEN
! first set default based on observation technique
     IF(obs%tech.eq.'P'.or.obs%tech.eq.'A'.or.obs%tech.eq.'N')THEN
! Photografic: Default value of RMS (time dependent, in arcsec)
! Before 1890
        IF(tdt.LT.11368.d0) THEN
           rmsa=3.d0
           rmsd=3.d0
! From 1890 to 1950
        ELSE IF(tdt.LT.33282.d0) THEN
           rmsa=2.d0
           rmsd=2.d0
! After 1950
        ELSE
           rmsa=1.5d0
           rmsd=1.5d0
        ENDIF
!Occultations
     ELSEIF(obs%tech.eq.'E')THEN
        rmsa=0.2d0
        rmsd=0.2d0
!Hypparcos
     ELSEIF(obs%tech.eq.'H')THEN
        rmsa=0.4d0
        rmsd=0.4d0
!Transit circle
     ELSEIF(obs%tech.eq.'T')THEN
        rmsa=0.5d0
        rmsd=0.5d0
!Encoder
     ELSEIF(obs%tech.eq.'e')THEN
        rmsa=0.75d0
        rmsd=0.75d0
!Micrometer
     ELSEIF(obs%tech.eq.'M')THEN
        rmsa=3.0d0
        rmsd=3.0d0
!CCD - normal place - roving - satellites
     ELSEIF(obs%tech.eq.'c'.or.obs%tech.eq.'C'.or.obs%tech.eq.'n'.or.obs%tech.eq.'V'.or.obs%tech.eq.'S')THEN
        rmsa=1.0d0
        rmsd=1.0d0
!Discovery observations - normally duplicated or very poor
     ELSEIF(obs%tech.eq.'x'.or.obs%tech.eq.'X')THEN
        rmsa=999.d0
        rmsd=999.d0
     ELSE
        WRITE(*,*)' obs. type unknown', obs%type,' ', obs%tech,' ',obs%obscod_s
        STOP 
     ENDIF
! use rules to overwrite default, but only for CCD observations 
     IF(obs%tech.eq.'c'.or.obs%tech.eq.'C'.or.obs%tech.eq.'n'.or.obs%tech.eq.'V'.or.obs%tech.eq.'S')THEN
        obscat=observ//obs%catcodmpc(2:2) ! four character string for a rule
        CALL bin_search(obscat,nrules,obscatcod,indrules,indobs)
! use rms appropriate for observatory, star catalog (in arcsec)
        IF(indobs.gt.0.and.indobs.le.nrules)THEN
           rmsa=rmsbin(1,indobs)
           rmsd=rmsbin(2,indobs)
        ELSE
           obscat=observ//'*'
           CALL bin_search(obscat,nrules,obscatcod,indrules,indobs)
! use rms appropriate for observatory, any catalog
           IF(indobs.gt.0.and.indobs.le.nrules)THEN
              rmsa=rmsbin(1,indobs)
              rmsd=rmsbin(2,indobs)
           ELSE
              obscat='ALL'//obs%catcodmpc(2:2)
              CALL bin_search(obscat,nrules,obscatcod,indrules,indobs)
! use rms appropriate for catalog, any station
              IF(indobs.gt.0.and.indobs.le.nrules)THEN
                 rmsa=rmsbin(1,indobs)
                 rmsd=rmsbin(2,indobs)
              ENDIF
           ENDIF
        ENDIF
     ENDIF
! old weighting for cbm10
  ELSEIF(error_model.eq.'cbm10')THEN
     IF(biasflag)THEN
! use rms appropriate for observatory, whatever the star catalog
        CALL bin_search(observ,nstarms,obscodrms,indstarms,indobs)
        IF(indobs.gt.0.and.indobs.le.nstarms)THEN
           rmsa=rmsstabin(1,indobs)
           rmsd=rmsstabin(2,indobs)
           ermrms=.true.
        ELSE
! use rms appropriate for observatory, star catalog (in arcsec)
           obscat=observ//obs%catcodmpc(2:2) ! four character string for a rule
           CALL bin_search(obscat,nrules,obscatcod,indrules,indobs)
           IF(indobs.gt.0.and.indobs.le.nrules)THEN
              rmsa=rmsbin(1,indobs)
              rmsd=rmsbin(2,indobs)
              ermrms=.true.
           ELSE
              obscat='ALL'//obs%catcodmpc(2:2) ! rule for all other obs with given catalog
              CALL bin_search(obscat,nrules,obscatcod,indrules,indobs)
              IF(indobs.gt.0.and.indobs.le.nrules)THEN
                 rmsa=rmsbin(1,indobs)
                 rmsd=rmsbin(2,indobs)
                 ermrms=.true.
              ELSE
!              WRITE(iwarou,'(A,A,A,A,A,F11.5)')'astrow_bias_cbm: case without weights ',obscat, &
!    &     ' Observatory: ',observ, ' Time UTC: ',obs%time_utc
                 numwar=numwar+1
                 ermrms=.false.
! Default value of RMS (time dependent, in arcsec)
                 rmsa=1.5d0
                 rmsd=1.5d0
              ENDIF
           ENDIF
        ENDIF
     !!! for non-CCD default RMS
     !!IF(obs%tech.ne.'c'.or.obs%tech.ne.'C')THEN
     !!  ! Default value of RMS (time dependent, in arcsec)
     !!  rmsa=1.5d0
     !!  rmsd=1.5d0
     !!ENDIF
     ELSEIF(obs%tech.eq.'c'.or.obs%tech.eq.'C'.OR.&
!!! DF when the catalog is known, use proper weights
          obs%catcodmpc.EQ.' h'.OR.obs%catcodmpc.EQ.' i'&
          .OR.obs%catcodmpc.EQ.' j'.OR.obs%catcodmpc.EQ.' k'&
          .OR.obs%catcodmpc.EQ.' l'.OR.obs%catcodmpc.EQ.' m'&
          .OR.obs%catcodmpc.EQ.' n'.OR.obs%catcodmpc.EQ.' p'&
          .OR.obs%catcodmpc.EQ.' v'.OR.obs%catcodmpc.EQ.' w'&
          .OR.obs%catcodmpc.EQ.' x'.OR.obs%catcodmpc.EQ.' z')THEN
! only CCD observations have special model
! is there a rule for this observatory?
        CALL bin_search(observ,nstarms,obscodrms,indstarms,indobs)
        IF(indobs.gt.0.and.indobs.le.nrules)THEN
           rmsa=rmsstabin(1,indobs)
           rmsd=rmsstabin(2,indobs)
           ermrms=.true.
        ELSE
! use rms appropriate for observatory, star catalog (in arcsec)
           obscat=observ//obs%catcodmpc(2:2) ! four character string for a rule
           CALL bin_search(obscat,nrules,obscatcod,indrules,indobs)
           IF(indobs.gt.0.and.indobs.le.nrules)THEN
              rmsa=rmsbin(1,indobs)
              rmsd=rmsbin(2,indobs)
              ermrms=.true.
           ELSE
              obscat='ALL'//obs%catcodmpc(2:2) ! rule for all other obs with given catalog
              CALL bin_search(obscat,nrules,obscatcod,indrules,indobs)
              IF(indobs.gt.0.and.indobs.le.nrules)THEN
                 rmsa=rmsbin(1,indobs)
                 rmsd=rmsbin(2,indobs)
                 ermrms=.true.
              ELSE
!              WRITE(iwarou,'(A,A,A,A,A,F11.5)')'astrow_bias_cbm: case without weights ',obscat,  &
!&         ' Observatory: ',observ, ' Time UTC: ',obs%time_utc
                 numwar=numwar+1
                 ermrms=.false.
! Default value of RMS (time dependent, in arcsec)
                 rmsa=1.5d0
                 rmsd=1.5d0
              ENDIF
           ENDIF
        ENDIF
     ELSE
        ermrms=.false.
! non-CCD observations
!!! DF transit type
        IF(obs%tech.eq.'t'.or.obs%tech.eq.'T')THEN
           rmsa=0.5d0
           rmsd=0.5d0
! Default value of RMS (time dependent, in arcsec)
! Before 1890
        ELSEIF(tdt.LT.11368.d0) THEN
           rmsa=3.d0
           rmsd=3.d0
! From 1890 to 1950
        ELSE IF(tdt.LT.33282.d0) THEN
           rmsa=2.d0
           rmsd=2.d0
! After 1950
        ELSE
           rmsa=1.2d0
           rmsd=1.2d0
        ENDIF
     ENDIF
  ENDIF
! change rms value only if it is not forced
  IF(.not.obsw%force_w(1))THEN 
     obsw%rms_coord(1)=MAX(obs%acc_coord(1),rmsa)*radsec/cos(obs%coord(2))
  ENDIF
  IF(.not.obsw%force_w(2))THEN
     obsw%rms_coord(2)=MAX(obs%acc_coord(2),rmsd)*radsec
  ENDIF
END SUBROUTINE astrow_bias_cbm
! ====================================================================
! ASTROW_NOBIAS
! No bias from catalogs, assign RMS
! ====================================================================
! INPUT:    MPCTYP    -  Observation type (column 15 of MPC record)
!           TDT       -  Time of observation (MJD, TDT)
!           IDSTA     -  Observatory code 9numeric)
!           ACCA      -  Accuracy of right ascension (rad)
!           ACCD      -  Accuracy of declination (rad)
!
! OUTPUT:   RMSA      -  A-priori RMS of right ascension (rad)
!           RMSD      -  A-priori RMS of declination (rad)
!           BIASA     -  Bias in  right ascension (rad)
!           BIASD     -  Bias in  declination (rad)
!           DECSTEP   -  Decision steps of RMS assignment
!
SUBROUTINE astrow_nobias(error_model,obs,obsw)
  TYPE(ast_obs),  INTENT(IN)    :: obs         ! observation
  CHARACTER*(*),  INTENT(IN)    :: error_model ! code for errmod
  TYPE(ast_wbsr), INTENT(INOUT) :: obsw        ! assigned weight

! ==============END INTERFACE==================================
!mpctyp,tdt,idsta,acca,accd,        &
!     &        rmsa,rmsd,biasa,biasd,decstep)
  DOUBLE PRECISION tdt,acca,accd
  CHARACTER*(1) mpctyp
  DOUBLE PRECISION rmsa,rmsd,biasa,biasd, rmsmin
  LOGICAL error
! =============================================================
! no bias anyway
  obsw%bias_coord(1:2)=0.d0
! default, rough rule-of-thumb error model
! Default value of RMS (time dependent)
  tdt=obs%time_tdt
! Before 1890
  IF(tdt.LT.11368.d0) THEN
     rmsmin=3.d0
! From 1890 to 1950
  ELSE IF(tdt.LT.33282.d0) THEN
     rmsmin=2.d0
! After 1950
  ELSE
     rmsmin=1.d0
  END IF
  rmsmin=rmsmin*radsec
  acca=obs%acc_coord(1)
  accd=obs%acc_coord(2)
  rmsa=MAX(rmsmin,acca)/cos(obs%coord(2))
  rmsd=MAX(rmsmin,accd)
  obsw%rms_coord(1)=rmsa
  obsw%rms_coord(2)=rmsd
END SUBROUTINE astrow_nobias
!***********************************************************************
! 'magrms' returns the default magnitude rms based on the MPC obs string
!***********************************************************************
DOUBLE PRECISION FUNCTION magrms_new(magstr,tdt,idsta,typ)
  CHARACTER*6 magstr
! obs. type (from column 15 of .obs format), station code, time MJD
  CHARACTER*1 typ
  INTEGER idsta
  DOUBLE PRECISION tdt
  INTEGER ll,lench
! magnitude weighting for now is simple;
! should be based on digits given,color
  ll=lench(magstr)
  IF(ll.le.0)THEN
     magrms_new=-1.d0
     RETURN
  ENDIF
  IF(magstr(3:5).eq.'   ')THEN
     magrms_new=1.0
  ELSEIF(magstr(5:5).eq. ' ')THEN
     magrms_new=0.7
  ELSE
     magrms_new=0.5
  ENDIF
END  FUNCTION magrms_new

! END SUB_MODULE observ_bias_weight
END MODULE astrometric_observations
