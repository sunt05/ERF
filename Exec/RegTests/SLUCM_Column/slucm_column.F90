!===============================================================================
! Offline single-column driver for the SLUCM physics, urban().
!
! Purpose (dev/plan-slucm-coupling.md section 8, item 6): drive urban() with an
! explicit forcing file so that ERF's vendored module_sf_urban.F90 and the
! upstream WRF phys/module_sf_urban.F it came from can be run under *identical*
! inputs and their trajectories compared. Both are built from this same driver;
! only the module source and its two `use` dependencies differ.
!
! It deliberately does NOT go through NoahmpUrbanDriverMainMod: that would drag
! in the whole NoahmpIO type and with it the coupling this test is meant to hold
! constant. The initialisation below instead mirrors that driver's own call site
! (drivers/erf/NoahmpUrbanDriverMainMod.F90, around the `call urban(`), so the
! state urban() sees here is the state it sees inside ERF.
!
! Everything urban() reads as INTENT(IN) comes from the forcing file, so the two
! builds cannot diverge through the driver. The only scalars fixed here are the
! ones that describe the column rather than the weather.
!===============================================================================
program slucm_column

  use module_sf_urban, only : urban, urban_param_init, LU_DATA_TYPE, ICATE
#ifdef HAVE_KIND_URBAN_REAL
  ! ERF-local addition; the pristine upstream file does not export it, which is
  ! why this driver must not depend on it to compile.
  use module_sf_urban, only : KIND_URBAN_REAL
#endif

  implicit none

  ! ---- column definition (not weather; identical on both sides) --------------
  integer, parameter :: NL     = 4        ! roof/wall/road layers, WRF default
  integer, parameter :: UTYPE_ = 2        ! high-intensity residential; matches
                                          ! the 3-D smoke test's UTYPE = 2
  integer, parameter :: JMONTH_ = 7
  real,    parameter :: XLAT_  = 40.0     ! deg
  real,    parameter :: ZA_    = 20.0     ! reference height [m]

  ! ---- urban() dummy arguments ----------------------------------------------
  logical :: LSOLAR
  integer :: num_roof_layers, num_wall_layers, num_road_layers
  real, dimension(NL) :: DZR, DZB, DZG
  integer :: UTYPE
  real :: TA, QA, UA, U1, V1, SSG, SSGD, SSGQ, LLG, RAIN, RHOO
  real :: ZA, DECLIN, COSZ, OMG, XLAT, DELT, ZNT
  real :: CHS, CHS2
  real :: TR, TB, TG, TC, QC, UC
  real, dimension(NL) :: TRL, TBL, TGL
  real :: XXXR, XXXB, XXXG, XXXC
  real :: TS, QS, SH, LH, LH_KINEMATIC
  real :: SW, ALB, LW, G, RN, PSIM, PSIH
  real :: GZ1OZ0
  real :: CMR_URB, CHR_URB, CMC_URB, CHC_URB
  real :: U10, V10, TH2, Q2, UST
  real :: mh_urb, stdh_urb
  real, dimension(4) :: lf_urb
  real :: lp_urb, hgt_urb, frc_urb, lb_urb, zo_check
  real :: CMCR, TGR
  real, dimension(NL) :: TGRL, SMR
  real :: CMGR_URB, CHGR_URB
  integer :: jmonth
  real :: DRELR, DRELB, DRELG, FLXHUMR, FLXHUMB, FLXHUMG
  real :: TVG, TT, XXXVG
  real, dimension(NL) :: TVGL, SMG
  real :: CMCG, FLXHUMVG, FLXHUMT
  real :: lf_urb_s, z0_urb, vegfrac_in

  ! ---- driver bookkeeping ----------------------------------------------------
  integer, parameter :: FIN = 21, FOUT = 22
  character(len=512) :: line
  character(len=256) :: forcing_file, output_file
  integer :: ios, nstep, k, nargs
  real    :: t_hours

  nargs = command_argument_count()
  if (nargs >= 1) then
     call get_command_argument(1, forcing_file)
  else
     forcing_file = 'forcing_48h.csv'
  end if
  if (nargs >= 2) then
     call get_command_argument(2, output_file)
  else
     output_file = 'trajectory.csv'
  end if

  !-----------------------------------------------------------------------------
  ! The promotion assertion the vendored header advertises. urban() is written
  ! in bare REAL; ERF compiles it with -fdefault-real-8. If that flag is ever
  ! dropped, every argument this driver passes would be silently reinterpreted,
  ! so fail loudly instead.
  !-----------------------------------------------------------------------------
  ! Works on both sides: the driver's own default REAL is promoted by the same
  ! flag that promotes urban(), so if this is 4 the promotion was dropped.
  if (kind(TA) /= kind(1.0d0)) then
     write(*,'(a,i0,a)') 'FATAL: real promotion was not applied (kind = ', &
                         kind(TA), ', expected 8).'
     stop 1
  end if
#ifdef HAVE_KIND_URBAN_REAL
  ! Stronger, ERF only: confirm the *module* translation unit agrees, which is
  ! what KIND_URBAN_REAL exists to prove.
  if (KIND_URBAN_REAL /= kind(TA)) then
     write(*,'(a)') 'FATAL: driver and module_sf_urban disagree on REAL kind.'
     write(*,'(a,i0,a,i0)') '       KIND_URBAN_REAL = ', KIND_URBAN_REAL, &
                            ', driver = ', kind(TA)
     stop 1
  end if
#endif

  num_roof_layers = NL
  num_wall_layers = NL
  num_road_layers = NL

  ! LU_DATA_TYPE and ICATE are module state that WRF sets from the land-use
  ! dataset before urban_param_init reads the table; they are not arguments.
  LU_DATA_TYPE = 'USGS'
  ICATE        = 33

  ! Reads URBPARM.TBL from the working directory and fills the *_TBL arrays and
  ! the AHOPTION / ALHOPTION / TREEOPTION / GROPTION / IRI_SCHEME switches.
  ! sf_urban_physics = 1 (SLUCM), use_wudapt_lcz = 0, distributed drag off --
  ! the configuration ERF ships and asserts.
  DZR = 0.0 ; DZB = 0.0 ; DZG = 0.0
  call urban_param_init(DZR, DZB, DZG, NL, 1, 0, .false.)

  !-----------------------------------------------------------------------------
  ! Initial state. Mirrors NoahmpUrbanDriverMainMod's call site: an isothermal
  ! fabric at the initial air temperature, zero canopy water, and the in-canyon
  ! vegetation state seeded from the ground/canopy values so that urban()'s
  ! unconditional *P prologue copies defined memory (TREEOPTION == 0 here, so
  ! nothing downstream consumes it).
  !-----------------------------------------------------------------------------
  TR = 290.0 ; TB = 290.0 ; TG = 290.0 ; TC = 290.0
  QC = 0.010 ; UC = 1.0
  TRL = 290.0 ; TBL = 290.0 ; TGL = 290.0
  XXXR = 0.0 ; XXXB = 0.0 ; XXXG = 0.0 ; XXXC = 0.0

  TGR  = 290.0 ; TGRL = 290.0 ; SMR = 0.3 ; CMCR = 0.0
  DRELR = 0.0 ; DRELB = 0.0 ; DRELG = 0.0
  FLXHUMR = 0.0 ; FLXHUMB = 0.0 ; FLXHUMG = 0.0

  TVG = TG ; TT = TC ; XXXVG = XXXG
  TVGL = TGL ; SMG = SMR
  CMCG = 0.0 ; FLXHUMVG = 0.0 ; FLXHUMT = 0.0

  ! Exchange coefficients. urban() uses CHS/CHS2 only for the diagnostic
  ! TS = TA + FLXTH/CHS and QS = QA + FLXHUM/CHS; no prognostic state depends on
  ! them (see the note at the ERF call site).
  CHS  = 0.01 ; CHS2 = 0.01
  CMR_URB = 0.0 ; CHR_URB = 0.0 ; CMC_URB = 0.0 ; CHC_URB = 0.0
  CMGR_URB = 0.0 ; CHGR_URB = 0.0

  ! NUDAPT morphology absent, exactly as urban_var_init leaves it when no
  ! high-resolution urban dataset is supplied.
  mh_urb = 0.0 ; stdh_urb = 0.0 ; lp_urb = 0.0 ; hgt_urb = 0.0
  lf_urb = 0.0 ; lb_urb = 0.0 ; zo_check = 0.0
  frc_urb = 0.9                      ! matches the smoke test's FRC_URB2D
  lf_urb_s = 0.0 ; z0_urb = 0.0      ! distributed aerodynamics only
  vegfrac_in = 0.0

  ZNT    = 0.5
  LSOLAR = .false.                   ! urban() splits SSG into SSGD/SSGQ itself
  SSGD   = 0.0 ; SSGQ = 0.0
  UTYPE  = UTYPE_
  jmonth = JMONTH_
  XLAT   = XLAT_
  ZA     = ZA_

  !-----------------------------------------------------------------------------
  ! Run
  !-----------------------------------------------------------------------------
  open(unit=FIN, file=trim(forcing_file), status='old', action='read', iostat=ios)
  if (ios /= 0) then
     write(*,'(a)') 'FATAL: cannot open forcing file '//trim(forcing_file)
     stop 1
  end if

  open(unit=FOUT, file=trim(output_file), status='replace', action='write')
  ! 17 significant digits round-trips an IEEE double exactly, so two builds that
  ! agree bitwise produce byte-identical files.
  write(FOUT,'(a)') 'step,time_h,TS,TR,TB,TG,SH,LH,G,RN,TC,UST,ALB,LW'

  nstep = 0
  read(FIN,'(a)', iostat=ios) line     ! header
  do
     read(FIN,'(a)', iostat=ios) line
     if (ios /= 0) exit
     if (len_trim(line) == 0) cycle
     read(line,*) DELT, TA, QA, UA, U1, V1, SSG, LLG, RAIN, RHOO, COSZ, OMG, DECLIN

     call urban(LSOLAR,                                           &
                num_roof_layers, num_wall_layers, num_road_layers, &
                DZR, DZB, DZG,                                    &
                UTYPE, TA, QA, UA, U1, V1, SSG, SSGD, SSGQ, LLG, RAIN, RHOO, &
                ZA, DECLIN, COSZ, OMG, XLAT, DELT, ZNT,           &
                CHS, CHS2,                                        &
                TR, TB, TG, TC, QC, UC,                           &
                TRL, TBL, TGL,                                    &
                XXXR, XXXB, XXXG, XXXC,                           &
                TS, QS, SH, LH, LH_KINEMATIC,                     &
                SW, ALB, LW, G, RN, PSIM, PSIH,                   &
                GZ1OZ0,                                           &
                CMR_URB, CHR_URB, CMC_URB, CHC_URB,               &
                U10, V10, TH2, Q2, UST, mh_urb, stdh_urb, lf_urb, &
                lp_urb, hgt_urb, frc_urb, lb_urb, zo_check,       &
                CMCR, TGR, TGRL, SMR, CMGR_URB, CHGR_URB, jmonth, &
                DRELR, DRELB, DRELG, FLXHUMR, FLXHUMB, FLXHUMG,   &
                TVG, TT, XXXVG, TVGL, SMG, CMCG, FLXHUMVG, FLXHUMT, &
                lf_urb_s, z0_urb, vegfrac_in                      &
                )

     nstep   = nstep + 1
     t_hours = real(nstep) * DELT / 3600.0

     write(FOUT,'(i0,",",13(es26.17e3,:,","))') &
          nstep, t_hours, TS, TR, TB, TG, SH, LH, G, RN, TC, UST, ALB, LW
  end do

  close(FIN)
  close(FOUT)

  write(*,'(a,i0,a)') 'SLUCM_COLUMN_OK steps=', nstep, ''
  write(*,'(a,i0)')   'driver REAL kind=', kind(TA)
#ifdef HAVE_KIND_URBAN_REAL
  write(*,'(a,i0)')   'KIND_URBAN_REAL =', KIND_URBAN_REAL
#endif
  write(*,'(a)')      'final state:'
  write(*,'(a,es26.17e3)') '  TS = ', TS
  write(*,'(a,es26.17e3)') '  SH = ', SH
  write(*,'(a,es26.17e3)') '  LH = ', LH
  write(*,'(a,es26.17e3)') '  G  = ', G
  write(*,'(a,es26.17e3)') '  RN = ', RN
  do k = 1, NL
     write(*,'(a,i0,a,es26.17e3)') '  TGL(', k, ') = ', TGL(k)
  end do

end program slucm_column
