!
! (mdj) Replacement modules for common blocks and global variables
!  - should also subsume functions for better argument validation
!
!
module htvars
  implicit none
  !
  ! in case we want to change precisions at some point,
  !  for now sp seems to be enough
  !
  integer,parameter :: sp=KIND(1.0),  &
       & dp=SELECTED_REAL_KIND(2*PRECISION(1.0_sp))
  !
  ! constants for dimensioning - switch to allocatable arrays
  !
  integer,PARAMETER :: NVAR=24,NPAR=19
  real,PARAMETER :: PI=3.1415927
  logical :: flag_interp,flag_umbrella, flag_fieldtxt
end module htvars

module htpaths
  use htvars
  implicit none
  real(kind=sp) :: z
  real(kind=sp) :: grav,rv,ra,cv,ca,cp,m,rs,rhog,mu,sigma,fraction, &
       & rho0,ff,h,t0,en,v,sj,hj,alpha,beta,emall,u,b,rho,t,lmf,fo, &
       & rhol, bb, bmax, plumezwidth, plumehwidth, emallhb, bvf, nn ! hb,qhb in umb_mod
 ! bvf is a measured b-v freq from a local density gradient. nn is a given b-v freq**2
  real(kind=sp) :: vso(npar)
 ! replace common blocks with module
 ! COMMON /PATH/ Z
 ! COMMON /PATH1/ PI,GRAV,RV,RA,CV,CA,CP,M,RS,RHOG,MU,SIGMA, &
 !      &     RHO0,FF,H,T0,EN,V,SJ,HJ,ALPHA,BETA,EMALL,U,B,RHO,T,LMF,FO
 ! COMMON /PATH2/ VSO
end module htpaths

module htdata
  use htvars
  implicit none
  integer, parameter :: MAXNAMELEN=60,MAXLINELEN=200
  character(LEN=MAXNAMELEN) :: dataname,volc,eruptdate,plumeshape, &
       & phifrac,phidist,kindex,volclat,volclong
  integer :: index_height ! current height index, so we don't have to keep
                          ! searching on it
  integer :: len_dataname,col_h,col_t,col_v,col_rho,ndata
  real(kind=sp),allocatable :: dat_h(:),dat_t(:),dat_v(:),dat_rho(:)
end module htdata

!
! variables for umbrella contribution (strongwind.f)
!
module umb_mod
  use htvars
  implicit none
  !
  ! specific new additions to incorporate umbrella cloud,
  !  c.f. Marcus' strongwind.f code
  !
  integer,PARAMETER :: MAXSTEP=200
  integer,PARAMETER :: unit_mi=20,unit_umbr=21,unit_down=22
  integer :: index,fi,mind,iv,ido

  real(kind=sp),parameter :: g=9.82
  real(kind=sp),parameter :: n1=0.014,n1s=0.02,n2=0.001225,ps=0.226

  real(kind=sp) :: ekh,term,distto

  real(kind=sp) :: cd,cdold,dens,deltat,deltah0,diff,di,enu,hb, &
       ht,ratio1,ratio2,ratio3,lamda,re,Qhb,Rhb
  real(kind=sp) :: x_umb,xend,y_umb(nvar),yprime(nvar),step,incr,final
  real(kind=sp) :: Qo,Uo,wind
  real(kind=sp) :: p_umb(npar),scmm(npar),scfi(npar)

  real(kind=sp) :: xx(MAXSTEP),w_umb(MAXSTEP),v_umb(MAXSTEP),f(MAXSTEP),m_umb(NPAR,MAXSTEP),&
       & mprime(NPAR,MAXSTEP),a(NPAR,MAXSTEP),b_umb(MAXSTEP)

  external derivs_mod
end module umb_mod

!
! parse.f90
! $Author: jonesm $
! $Date: 2010-07-12 20:20:34 $
! $Id: parse.f90,v 1.7 2010-07-12 20:20:34 jonesm Exp $
! $Revision: 1.7 $
!
! <TB> support routine to parse free format files.
!
MODULE parse
  implicit none
  ! 
  ! Designed to be a self-contained module for the parse subroutine,
  !  which reads free format files, returning the contents of a line
  !  as an array of individual strings.  These strings can then be
  !  translated into integers or doubles as required.
  !
  ! Summer 2001 (mdj)
  ! Last Modified: 09/17/10 (mdj)
  !  - added modification for quoted strings that can contain spaces
  !
  ! MXDLEN - maximum line length
  ! MXDTERM - maximum number of terms per line
  ! MXDTERMLEN - maximum character length of each term
  ! cterms(iterm,*) - ith term on input line
  ! nterms - number of distinct terms
  integer,parameter :: idble=SELECTED_REAL_KIND(10,100)
  integer,parameter :: MXDLEN=100,MXDTERM=20,MXDTERMLEN=80
  character(LEN=MXDTERMLEN) :: cterms(MXDTERM)
  integer :: nterms

CONTAINS
  integer function lparse(iunit)
  implicit none

  ! Dummy arguments
  integer :: iunit
  
  ! Locals
  logical :: newterm,insidequotes
  integer :: istart,iterm,lenw,clen,indx,i,j,iline
  character(LEN=MXDLEN) :: cdata

  ! Read current line, parsing off any characters following a comment
  !  
  iline=0
  do
     indx=0
     clen=0
     istart=1
!     read (unit=iunit,fmt=*) cdata
     read (unit=iunit,fmt='(a100)',end=2,err=2) cdata
!     do ! if the entire line is blank, this fails, use ADJUSTL intrinsic
!        if (cdata(1:1).eq.' ') then
!           do i=2,MXDLEN
!              cdata(i-1:i-1)=cdata(i:i)
!           end do
!           cdata(MXDLEN:MXDLEN)=' '
!        else
!           exit
!        endif
!     end do
     cdata = ADJUSTL(cdata)
!     cdata = TRIM(cdata)
!     clen = LEN_TRIM(cdata)
     do i=MXDLEN,1,-1
        if (cdata(i:i).ne.' ') exit
     enddo
     clen=i
     iline=iline+1
     ! print*,'<lparse> line,clen,string = ',iline,clen,':',cdata(1:clen),':'
     insidequotes=.false.
     do i=1,clen
        if (cdata(i:i).eq.'"') insidequotes=.not.insidequotes
        if ((cdata(i:i).eq.'#'.or.cdata(i:i).eq.';'.or.cdata(i:i).eq.':').and.(.not.insidequotes)) exit
!           print*,' lparse: indx,cdata(1:1) :',indx,cdata(i:i)
!           exit
!        else
           indx=i
!        endif
     end do
     indx=indx+1
     cdata(indx:indx)=' '
     clen = indx
     ! if indx=1, then we have no non-comment terms, and need a new line
     if(indx.gt.1) exit
  end do
  ! write(*,*) '<lparse> indx,cdata = ',indx,cdata(1:indx)
 
  ! Real parsing - delimiters between terms are space or comma, or a combo
  i=0
  istart=1
  iterm=1
  newterm=.false.
  insidequotes=.false.
  do
     i = i+1
     if (i.gt.clen) exit
     !print*,'i,cdata(i:i)=',i,cdata(i:i) 
     if ((cdata(i:i).eq.'"')) insidequotes=.not.insidequotes
     if ((cdata(i:i).eq.' '.or.cdata(i:i).eq.',').and.(.not.insidequotes)) then

        if (.not.newterm) then
           lenw = i-istart
           if (lenw.gt.MXDTERMLEN) then
              write(*,*) '<lparse> MXDTERMLEN exceeded. ',lenw,MXDTERMLEN,cdata(istart:i)
              stop
           endif
           cterms(iterm) = cdata(istart:i-1)
!           do j=1,lenw
!              cterms(iterm,j:j) = cdata(istart+j-1:istart+j-1)
!           enddo
!           do j=lenw+1,MXDTERMLEN
!              cterms(iterm,j)=' '
!           enddo
           newterm=.true.
           ! write(*,*) '<lparse> term found, iterm, term = ',iterm,lenw,cterms(iterm),cdata(istart:i-1)
        endif
     else
        if (newterm) then ! new term
           istart=i
           newterm=.false.
           iterm=iterm+1
           if(iterm.gt.MXDTERM) then
              write(*,*) '<lparse> MXDTERM exceeded.',iterm,MXDTERM
              stop
           endif
        endif
     endif
  enddo
  nterms=iterm
  ! write(*,*) '<lparse>: line parsed.  ',iterm,' terms found.'
  lparse=0
  return
2 continue
!2 write(*,*) '<lparse>: EOF found.'
  lparse=1
  return
  end function lparse

  double precision function rlread(string)
  implicit none
  !
  ! Reads a double from the input string
  !
  character(*) string
  read(unit=string,fmt='(e30.20)',err=1) rlread
  return
1 write(*,'("<rlread>: error in conversion, string =",a30)') string
  STOP
  end function rlread

  integer function intread(string)
  implicit none
  !
  ! Reads an integer from the input string
  !
  character(*) string
  read(unit=string,fmt='(i30)',err=1) intread
  return
1 write(*,'("<intread>: error in conversion, string =",a30)') string
  STOP
  end function intread

  character(LEN=8) function itos(int)
    implicit none
    integer :: int
    write(unit=itos,fmt='(i8)',err=1) int
    itos = ADJUSTL(itos)
    itos = TRIM(itos)
    return
1   write(*,'("<itos>: error in conversion, int =",i8)') int
    STOP
  end function itos

END MODULE PARSE

PROGRAM PUFFIN
  !---------------------------------------------------------------C
  !     THIS IS Puffin (PLUME HEIGHT/Puff Input)
  !     Author: M. Bursik mib@buffalo.edu
  !     Modifications:
  !      - conversion to f90 (2010-05) (mdj)
  !      - added experimental radiosonde capability for temperature,
  !        density, and wind speed (2010-05) (mdj)
  !      - flexible input format by keyword (2010-08) (mdj)
  !
  ! Notes on radiosonde data:
  !      - verify temperature interpolation should be by atmosph. temp.
  !      - verify density as function substitution (rhoa)
  !      - completely ignore wind direction?  Looks like exp. data
  !        has shears/lots of noise
  !      - interpolation/extrapolation is still very crude
  !---------------------------------------------------------------C
  use htvars
  use htpaths
  use htdata
  use umb_mod
  use parse
  implicit none

  REAL(KIND=SP) :: YOUT(NVAR),Y(NVAR),DXF(NVAR),DY(NVAR),XFF(NVAR)
  REAL(KIND=SP) :: MV(NPAR,1),PV(NPAR),DYDX(NVAR),frac(NVAR)
  REAL(KIND=SP) :: P,K2

  integer,parameter :: unit_in=13,unit_height=14,unit_puff=15,unit_gui=16
  integer :: ieof,itype,i,j,k,kk,maxiter,np,b0,len_string1,count
  character(LEN=8) :: string1
  character(LEN=24) :: stype(4)
  character(LEN=400) :: puff_string,puff_numstring

  logical :: exist_in
  real(kind=sp) :: c,h0,ends,en0,phimin,phimax,r,s,theta1, &
       & starts,theta,y3,y2,y4,u0,uz,w,x,y1,y24,zmax,dy24, &
       & ashlogmean,ashlogsdev,mu_sigma,mu_rand,n_gsd, &
       & UZ_old,UZ_next,H_d0,zmax_0,count_hb
  real(kind=sp) :: rhoa,ta
  external rhoa,ta
  integer :: bsearch
  external bsearch
  EXTERNAL DERIVS,derivs_umb
  external init_random_seed
!      
!      BEGIN EXECUTABLES
!      DEFINE CONSTANTS		
!      UNIVERSAL:
!      
  stype(1) = 'PHREATOPLINIAN'
  stype(2) = 'VULCANIAN'
  stype(3) = 'PLINIAN'
  stype(4) = 'STROMBOLIAN'
!  PI=3.14159    
  GRAV=9.807
  nn=0.003333  ! b-v freq squared (N-squared) for standard atmosphere
! from "Handbook of Acoustics" by Crocker p. 334
! note that (26 July 2012) umbrella module has its own brunt-vaisala freq
!      HEAT CAPS, GAS CONSTS FOR VAPOR, AIR, PARTICLES (S97)
  CV=1860.
  CA=1000.
  CP=1100.  !print*,'derivs: on entry, y(1),y(2)=',y(1),y(2)
  RA=285.
  RV=460.
!      ENTRAINMENT CONSTANTS
  ALPHA=0.15
  BETA=1.0
  !
  !
  ! Open input file - fail gracefully if not found.
  !
  INQUIRE(file='puffin.inp',exist=exist_in)
  if (.not.exist_in) then
     write(*,*) "Unable to open main input file, puffin.inp"
     STOP
  end if
  OPEN(unit_in,FILE='puffin.inp',status='old')
  rewind(unit_in)
  OPEN(unit_height,FILE='height.out',status='unknown')
  OPEN(unit_puff,FILE='run.puff',status='unknown')
  OPEN(unit_gui,FILE='gui.puff',status='UNKNOWN')
  OPEN(3,FILE='FIELD.TXT')
!
! parse input file for keywords and associated data
!  new fields are the RAWDATA and UMBRELLA keywords
!  for radiosonde data and stratification, respectively
! 
  B0=100.0 ! default vent radius
  mu=0.0 ! test for changes
  do
     ieof = lparse(unit_in)
     !print*, 'Read line from input, return = ',ieof
     if (ieof.ne.0) exit
     !write(*,'("keyword: ",a20," nterms = ",i3)') cterms(1),nterms
     select case (cterms(1))
     case ("VOLCANONAME")
        volc=cterms(2)
        volc=TRIM(volc)
        WRITE(*,*)'VOLCANO NAME:', volc
     case ("VOLCANOLATLONG")
        volclat=cterms(2)
        volclong=cterms(3)
     case ("ERUPTIONDATE")
        eruptdate=cterms(2)
        WRITE(*,*)'ERUPTION DATE:', eruptdate
     case ("PLUMESHAPE")
        plumeshape=cterms(2)
        WRITE(*,*)'PLUME SHAPE:', plumeshape
     case ("VENTRADIUS") 
        B0 = rlread(cterms(2))
        WRITE(*,*)'VENT RADIUS: ',B0,' M'
     case ("AXIALVELOCITY","VENTVELOCITY","JETVELOCITY") 
        U0 = rlread(cterms(2))
        WRITE(*,*)'AXIAL VENT/JET VELOCITY:',U0,' M/S'
     case ("WATERFRACTION")
        EN0 = rlread(cterms(2))
        WRITE(*,*)'WATER CONTENT:',EN0*100.,' WT %'
     case ("TEMPERATURE")
        T0 = rlread(cterms(2))
        WRITE(*,*)'ERUPTION TEMPERATURE:',T0,' K'
     case ("PLUMEANGLE")
        THETA1 = rlread(cterms(2))
        WRITE(*,*)'ANGLE THE PLUME MAKES WITH THE HORIZON:',THETA1,' DEG'
     case ("VENTHEIGHT")
        H0 = rlread(cterms(2))
        WRITE (*,*)'HEIGHT OF VENT ABOVE SEA LEVEL:',H0,' M'
     case ("WINDSPEED")
        V = rlread(cterms(2))
        WRITE(*,*)'MAX HORIZONTAL WIND SPEED: ',V,' m/s'
     case ("JETHEIGHT")
        HJ = rlread(cterms(2))
        WRITE(*,*)'MAX JET HEIGHT: ',HJ,' m'
     case ("JETWIDTH")
        SJ = rlread(cterms(2))
        WRITE(*,*)'MAX JET WIDTH: ',SJ,' m'
     case ("ERUPTIONTYPE")
        itype = intread(cterms(2))
        select case (itype)
        case (1)  ! phreatoplinian
           mu_sigma=4.0   ! (mib) sign fix 2010-08
           sigma=2.0
           n_gsd=4  ! need to put right value from Woods, Bursik (1991)
        case (2)  ! vulcanian
           mu_sigma=1.5
           sigma=1.0
           n_gsd=3
        case (3)   ! plinian
           mu_sigma=0.0
           sigma=3.0
           n_gsd=6
        case (4)   ! strombolian
           mu_sigma=-4.0
           sigma=2.5
           n_gsd=3
        case default
           write(*,*) 'Invalid eruption type (1<=itype<=4):',itype
           STOP
        end select
        WRITE(*,*) 'ERUPTION TYPE: ',stype(itype)
        ! call init_random_seed
     case ("PARTICLEDENSITY")
        RS = rlread(cterms(2))
        WRITE(*,*)'PARTICLE DENSITY:',RS, ' KG/M^3'
     case ("PARTICLESHAPE")
        FF = rlread(cterms(2))
        WRITE(*,*)'PARTICLE SHAPE FACTOR:',FF
        !      MOST ELONGATED PARTICLES-FF=0.4.  EQUANT-FF=1
     case ("PARTICLEFRACTION")
        !! DON'T CHANGE FRACTION FROM 1!
        !!  THIS WOULD BE THE PARAMETER TO MODIFY IF 
        !!  IT WAS DESIRED TO TRACK A PARTICULAR PARTICLE
        !!  PHASE THAT MIGHT BE A FRACTION OF TOTAL PARTICLE POP.
        FRACTION = rlread(cterms(2))
        WRITE(*,*)'FRACTION OF TOTAL PARTICLE POPULATION '
        WRITE(*,*)'   SUCH PARTICLES REPRESENT?',FRACTION
     case ("GRAINRANDOM")
        call init_random_seed
        call random_number(mu_rand)
        ! unif or gauss(??) distr on [mu-sigma,mu+sigma]
        ! taken out random seed in favor of value from input (mdj) 2010-09-29
        mu=mu_sigma+2.0*(mu_rand-0.5)*sigma/sqrt(n_gsd)
        write(*,*)'Using random grain size around mean=2.0, mu,sigma = ',mu,sigma
     case ("GRAINMEAN")
        mu = rlread(cterms(2))
        write(*,*)'Grain size, mean (phi scale):',mu
     case ("GRAINSDEV")
        sigma = rlread(cterms(2))
        write(*,*)'Grain size, std dev (phi scale):',sigma
     case ("PRINTGRAINSIZES")
        phimin = rlread(cterms(2))
        phimax = rlread(cterms(3))
        WRITE(*,*)'GRAIN-SIZES TO PRINT (PHIMIN PHIMAX):',phimin,phimax
     case ("RAWDATAFILE")
        dataname = cterms(2)
        flag_interp=.true.
        WRITE(*,*) 'Using raw data file: ',TRIM(dataname)
     case ("RAWDATACOLUMNS")
        if (nterms.ne.5) then
           STOP 'Need 4 columns in RAWDATACOLUMNS'
        end if
        col_h = intread(cterms(2))
        col_t = intread(cterms(3))
        col_v = intread(cterms(4))
        col_rho = intread(cterms(5))
        write(*,'(" columns for height: ",i2," temperature: ",i2, &
             & " wind speed: ",i2," density: ",i2)') col_h,col_t, &
             & col_v,col_rho
     case ("UMBRELLA") ! two arguments, turb. diffusivity and runout time in seconds
        flag_umbrella=.true.
        if (nterms.ne.4) then
           STOP 'Usage: UMBRELLA (float)diffusivity (float)runout_time (float)max_distance'
        end if
        ekh = rlread(cterms(2))
        term = rlread(cterms(3))
        distto = rlread(cterms(4))
        write(*,*) 'Umbrella feature on, stratification influence.'
     case("FIELDTXT") ! TODO
        flag_fieldtxt = .true.
        WRITE(*,*) 'Outputting FIELD.txt'
     case default
        write(*,*) 'Error, unrecognized keyword in input: ',cterms(1)
        STOP
     end select
  end do
  if (flag_interp) then
     CALL getdata
  end if
  print*,'flags - interp,umbrella = ',flag_interp,flag_umbrella
  !
  ! more initialization and sanity checking
  !
  if (mu.eq.0.0) then
     mu = mu_sigma ! pre-loaded value according to eruption type
  end if
  count=0
  count_hb=0
  bb=0
  bmax=0
  plumehwidth=0
  plumezwidth=0
  hb=h0
  qhb=0
  rho=0
  rhol=0
 
  !STOP 'dbx'
!      
!      HERE BEGINS THE LOOPING, LOOPS FOR EACH VALUE OF B0
!      
!      
!  remove loop over vent radii - make input parameter?
!    
!  DO B0=100,105,10
  WRITE(*,*)'VENT RADIUS:',B0,' M'
!      
!      SET UP INTEGRATION DOMAIN, SPACE VARIABLES  
!      INTEGRATE BETWEEN STARTS AND ENDS, W/ STEP H
!      
!(mdj)
! modify for experimental range or extrapolate (with all the accompanying
!  hazard so implied - interpolation is relatively easy (T and V look to be
!  pretty noisy)?
  STARTS=0.0
  ENDS=100000.
!      ENDS=5000.
!      H IS THE STEP SIZE
  H=1.0
!(mdj)
!      FOR PRINT SKIPPING
  KK=0
!      INITIALIZE S
  S=STARTS
!      INITIAL VALUES  
!      INITIALIZE THE POSITION VARIABLES
  X=0.0
  Z=H0
  !(mdj) zmax was not being initialized?
  zmax=0
  zmax_0=0
  count=0
  count_hb=0
  bb=0
  bmax=0
  plumehwidth=0
  plumezwidth=0
  hb=h0
  qhb=0
  rhog=0
  rho=0
  rhol=0
  THETA=THETA1*PI/180.
  ! (mdj) set initial interpolation point for experimental data
  !       currently assumes fixed value for height < dat_h(1)
  if (flag_interp) index_height = bsearch(z,dat_h,ndata)
  RHO0=RHOA(Z)
  B=B0*1.00000
  EN=EN0
  W=EN0
  T=T0  ! (mdj) T0 is eruption temperature entered above
  U=U0*1.00000
  UZ_old=0
  UZ_next=0
  R=(W/EN)*RV+(1-W/EN)*RA
  C=W*CV+(EN-W)*CA+(1-EN)*CP
  RHOG=RHOA(Z)*(RA*TA(Z)/R/T)
  RHO=1/((1-EN)/RS+EN/RHOG)
  M=RHO*B*B*U
  K2=RHO*U*U*B*B
  P=M*C*T
  LMF=M**(3/4)*P**(-1/2)
  FO = P
  EMALL=M*(1-EN)
  !call random_number(mu_rand)
  ! unif or gauss(??) distr on [mu-sigma,mu+sigma]
  ! taken out random seed in favor of value from input (mdj) 2010-09-29
  ! mu=mu_sigma+2.0*(mu_rand-0.5)*sigma/sqrt(n_gsd)
  
  WRITE(*,*)'INIT ATM RHOA0, TA0:',RHOA(Z),TA(Z)
  WRITE(*,*)'PLUME INIT. DENSITY:',RHO,'  ',RHOG
  WRITE(*,*)'INIT MASS FLUX:',M*PI
  WRITE(*,*)'INIT MOMENTUM FLUX:',K2*PI
  WRITE(*,*)'INIT THERMAL FLUX:',P*PI
  WRITE(*,*)'INIT MASS FLUX OF PARTICLES:',EMALL*PI
  WRITE(*,'("    mu = ",f6.3," sigma = ",f6.3," mu_rand = ",f6.3)') mu,sigma,mu_rand
  !      INITIAL VALUES TO GET INTEGRATED BY R-K
  Y(1)=M
  Y(2)=K2
  Y(3)=THETA
  Y(4)=P
  UZ=U*sin(y(3))
  !      CALCULATE FALL SPEEDS OF THE GRAINS AT VENT
  CALL FALLSPEED(PV)
!      INITIAL MASSES IN GRAIN SIZE FRACTIONS DEPEND ON
!      BULK DENSITY, HENCE EXSOLVED WATER CONTENT
  DO I=1,19
     MV(I,1)=PV(I)*EMALL*FRACTION
  end DO
  EMALL=0.0
  DO I=1,19
     Y(I+4)=MV(I,1)
     EMALL=EMALL+Y(I+4)
     frac(i+4)=0.0
  end DO
  Y(24)=EMALL
!      
!      
!      THIS BEGINS THE NUMERICAL INTEGRATION 
!      
!      
!      INITIALIZE THE VALUES OF DY
  CALL DERIVS(S,Y,DYDX)
!      MAX ITERATION SAFETY CATCH
! TODO I REDUCED THIS FROM 100K BECAUSE PUFF COMPLAINS IF I GO TOO HIGH
  MAXITER=100000.
90 CONTINUE ! target of GOTO statement below
  if (mod(int(s),1000).eq.0) write(*,'(" current height,uz,rho = ",2f8.1,1f8.4)') s,uz,rho
  UZ_old=UZ
  UZ=Y(2)/Y(1)*SIN(Y(3))
  !IF (flag_interp) 
  UZ_next=UZ+(UZ-UZ_old)/H*50. ! Would UZ go to 0 in 50 m?
  !print*,'KK,UZ = ',KK,UZ
  IF (S .LT. ENDS .AND. S .LT. MAXITER .AND. UZ_next .GE. 0.0) THEN
     CALL FALLSPEED (PV)
     !print*, 'pre-rk4 y(1)-y(4) =',yout(1),yout(2),yout(3),yout(4)
     !print*,'pre-rk4 vso = ',vso
     CALL RK4(Y,DYDX,NVAR,S,H,YOUT,DERIVS)
     !print*, 'post-rk4 y(1)-y(4) =',yout(1),yout(2),yout(3),yout(4)
     DO J=1,24
        DY(J)=(Y(J)-YOUT(J))
        if (y(j) .gt. 0.0) then ! (mib) 2010-09-01
           Y(J)=YOUT(J)
        else
           y(j)=0.0
        end if
     end DO
     Y1=Y(1)
     Y2=Y(2)
     Y3=Y(3)
     Y4=Y(4)
     Y24=Y(24)
     ! (mib) DY24 is fallout of all particles per meter
     DY24=DY(24)/(sin(y3)*h)
     RHOL=RHO
     CALL UPDATE(Y1,Y2,Y4,Y24,U,RHO,EN,RS,RHOG,B, &
          &           W,M,R,RV,RA,C,CV,CA,CP,T)
     RHOL=(RHO-RHOL)/H
     !      BACK CALCULATE X AND Z
     X=X+(COS(Y3)*H)
     Z=Z+(SIN(Y3)*H)
! (mib) To calculate HB, need to find where density of plume changes sign
! Since density diff becomes irregular, we only trap the first time 
! that the density diff changes sign
        IF (RHOL .NE. 0 .AND. RHOL .GT. (RHOA(Z)-RHOA(Z-1))/H .AND. RHO .GT. RHOA(Z) .AND. COUNT_HB .EQ. 0) THEN
        HB=Z
        print*,'rhol: ',rhol,'drhoa/dz: ',(rhoa(Z)-rhoa(z-1))/H,'Hb: ',HB
        QHB=Y1*PI/RHO
        BB=B
        emallhb=0.0
        ashlogmean=0.0
        ashlogsdev=0.0
        DO i=5,23
           emallhb=emallhb+y(i)
        END DO
        DO j=5,23
           frac(j)=y(j)/emallhb*1.0
           ashlogmean=ashlogmean+(j-13)*frac(j)
           !              WRITE(*,*)'Z:  ',Z,'B:  ',B,'y, frac',j,':  ',y(j),frac(j)
        END DO
        DO k=5,23
           frac(k)=y(k)/emallhb*1.0
           ashlogsdev=ashlogsdev+((k-13)-ashlogmean)**2.0*frac(k)
        END DO
        COUNT_HB=1
!           print*,'crazy Nans again: ashlog mean,sdev = ',ashlogmean,ashlogsdev
        ashlogmean=log10(0.001*2.0**(-ashlogmean))
        ashlogsdev=abs(log10(sqrt(0.001*2.0**(ashlogsdev)))) ! + definite
     END IF
!(mdj)        IF(Z/Z.EQ.1.0) ZMAX=Z  ! note that z/z is not nec. unity in fp
     if (flag_interp) index_height = bsearch(z,dat_h,ndata)  ! unsure of position here
     !print*,'y3,h,z,index_height = ',y3,h,z,index_height
     zmax_0=MAX(zmax_0,z)
     bmax=MAX(bmax,b)
!(mdj)
!      FALLOUTCOR CORRECTS X TO ACCOUNT FOR THE PLUME RADIUS AND
!      HORIZONTAL DISTANCE TRAVELED AFTER LEAVING THE PLUME.
     !print*,'pre-falloutcor vso = ',vso
!     CALL FALLOUTCOR(X,Y,DXF,B,XFF,Z,H)
     S=S+H
!      CONVERT TO FALLOUT PER UNIT AREA INSTEAD OF PER SLICE
     DO J=5,23
92      DY(J)=DY(J)/(2*B)
     end DO
!      PRINT OUT CALCULATIONS FROM THIS STEP***
!      C     NP=10
     H_d0=-UZ*H/(UZ-UZ_old)  ! Distance to 0, if it's going there.
!     WRITE(unit_height,*)'Z: ',Z,'UZ: ',UZ,'RHO: ',RHO,'RHOA(Z): ',RHOA(Z),'H_d0: ',H_d0
     BVF = ((grav/rhoa(z-1)))**0.5*(-(rhoa(Z)-rhoa(z-1))/H)**0.5
     WRITE(unit_height,*)'Z: ',Z,'UZ: ',UZ,'RHO: ',RHO,'RHOA(Z): ',RHOA(Z),'H_d0: ',H_d0,'BVF: ',BVF
     NP=99
     KK=KK+1
     IF(KK.LE.NP) GO TO 95
!        WRITE(*,*)'Z:  ',Z
!        WRITE(*,*)'Z:  ',Z,'B:  ',B,'Y24:  ',Y24,'DY24:  ',DY24
     if (flag_fieldtxt) then
         CALL FIELD(X,B,Z,U,Y(3),Y(11))
     endif

     KK=0
95   CONTINUE
!      RECALCULATE DY'S (PREPARING TO CALL RK4 AGAIN)
     GO TO 90
  ELSE
     zmax=zmax_0+H_d0
     plumehwidth=(bb+bmax)/2.0
     plumezwidth=(zmax-hb)/2.0
     WRITE(*,*)'NEUTRAL BUOYANCY HEIGHT:  ',HB,' M'
     WRITE(*,*)'ERUPTION PLUME HEIGHT:  ',ZMAX,' M'
     WRITE(*,*)'RADIUS AT HB:  ',BB,' M'
     WRITE(*,*)'ATM DENSITY AT HB:  ',RHOA(Z),' KG/CU M'    
     WRITE(*,*)'PARTICLE FLUX AT HB:  ',EMALLHB,' KG/S'
     WRITE(*,*)'VOLUME FLUX AT HB:  ',QHB,' CU M/S'
     WRITE(*,*)grav,rhoa(z),rhoa(z-1),h
!     BVF = ((grav/rhoa(z-1)))**0.5*(-(rhoa(Z)-rhoa(z-1))/H)**0.5
!     WRITE(*,*)'BRUNT-VAISALA FREQ: ',BVF,' /S'
     WRITE(unit_height,*)U0,B0,M*PI,ZMAX
     WRITE(*,*)B0,M*PI,ZMAX
     ! Note that plumehwidth and plumezwidth are in km
     !        write(unit_height,*)' '
     !
     ! puff_numstring used to compose string with numerical values,
     !  then incorporated into puff_string for full command line
     !
     write(puff_numstring,fmt=666) "--plumeMax",zmax,"--plumeMin",HB, &
          "--plumeHwidth",plumehwidth/1000.,"--plumeZwidth",plumezwidth/1000., &
          "--ashLogMean",ashlogmean,"--ashLogSdev",ashlogsdev
! TODO MADE plumeMax 10.3 plumeHwidth f10.3, up from f8.3
666  format(a10,1x,f10.1,1x,a10,1x,f8.1,1x,a13,1x,f10.3,1x,a13,1x,f8.3,1x,a12,f9.5,1x,a12,f9.5)
     puff_string = "puff --ashOutput=false --model=reanalysis --rcfile puffrc --volc " // &
          & volc(1:len_trim(volc)) // " --volcLat " // volclat(1:len_trim(volclat)) // & 
          & " --volcLon " // volclong(1:len_trim(volclong)) // &
          & " --eruptDate " // eruptdate(1:len_trim(eruptdate)) // &
          & " --plumeShape " // plumeshape(1:len_trim(plumeshape)) // &
          & " --runHours 72 --nAsh 4000000 " // &
          & "--gridOutput=True --gridBox -14.8819:32.8058/37.3622:68.3775/0:16000 " // &
          & puff_numstring(1:len_trim(puff_string)) 
     write(*,'(a400)') puff_string
     write(*,'("Length of puff string: ",i4)') len_trim(puff_string)
     write(unit_puff,'(a400)') TRIM(ADJUSTL(puff_string(1:len_trim(puff_string))))
     write(unit_gui,*) TRIM(ADJUSTL(puff_numstring(1:len_trim(puff_numstring))))
     WRITE(*,*)'puff command line is written to: run.puff'
     WRITE(*,*)'  '
  END IF
!97 CONTINUE
  ! STOP 'dbx'
     
  !
  ! Umbrella section, some initialization overlap with htcalc
  !  above, rest reproduced here from strongwind.f
  !
  if (flag_umbrella) then
     deltat=900.0
     t=300.0
     wind=v ! (mdj) not sure about this setting ...
     if (wind.le.0.0) then
        write(*,'("WARNING: input wind speed = ",e12.4," <= 0!")') wind
        write(*,'("         resetting wind = 1.0 m/s")')
        wind=1.0
     end if
     !ps=0.226
     !
     ! wrong flux!  need flux at bottom of umbrella cloud
     Qo = PI*B0*B0*U0
     Uo = U0 ! (mdj)-maybe? was fixed at 20.0
     incr=500 ! (mdj) - was parameter
     final=60000.0 ! (mdj)-was parameter
     distto = 5.e5
     ! calculate density and viscosity of air at vent height
     !  sea level density from NASA 1962 Standard Atmosphere
     !  viscosity from Middleton & Wilcock
     rho0=1.225*exp(-H*Nn/g)
     !
     !enu=0.000018/rho0
     ! enu from Atham results
     enu=0.37/rhog
     
     ! Calculate values for p (total of 19) based on
     ! grain size and eruption type, memorize millimeter scale
     ! and fi-scale:
     
     index=1
     do fi = -9, 9
        di=2.0**(-fi)
        ! save millimeter scale:
        scmm(index)=di
        ! save fi-scale:
        scfi(index)=fi
        
        p_umb(index)=EXP(-((fi-mu)**2.0)/(2.0*sigma*sigma))/(s*sigma)
        ! get vent-height fall speeds
        cd=1.0
        do iv = 1, 10
           Vso(index)=sqrt(4.0*(0.001*di)*g*(Rs-rho0)/(3.0*rho0*cd))
           Re=Vso(index)*(0.001*di)/enu
           cdold=cd
           cd=(24.0/Re)*(Ff**(-0.32))+2*sqrt(1.07-Ff)
           diff=abs((cdold-cd)/cd)
           ! use this output to check that cd is converging, if wanted.
           ! write(*,*) 'diff of ', index, ' = ', diff
        end do
        ! use this output to see fall speeds, if wanted.
        write(*,*)index,di,Vso(index)
        if (diff.ge.0.01) then
           write(*,*)'error in fall speed (index=', index,') >1%',' = ',diff
        endif
        index=index+1
     end do
     
     ! Initial conditions for all 22 functions:
     w_umb(1)=Qo/3.14159
     v_umb(1)=sqrt(Qo*Uo/PI)
     ! f is derived from relats in Bursik et al. (1992)
     ! and Wilson et al. (1978).  (Using F*=(2/pi)F in derivs)
     f(1)=deltaT*5.0*g*Qo*2.0/(rho0*T*3.14159)
     do i=1,19
        ! initial masses in grain size fractions depend on
        ! water content, chosen such that initial bulk density=5kg/m3,
        ! as above for f(1). (This represents exs. vols. of 0.03-0.05).
        m_umb(i,1)=p_umb(i)*Qo*5.0*Fraction
     end do
     
     ! Open files for output. Values of m1 through m19 will be saved 
     ! to mi.dat
     ! Values of v,f,w and m1 will be displayed on the screen.
     ! Output from "umbrella" will be saved to umbr.dat
     string1 = itos(B0)
     len_string1 = LEN_TRIM(string1)
     open(unit=unit_mi,status='unknown',file='mi_'//string1(1:len_string1)//'.dat')
     open(unit=unit_umbr,status='unknown',file='umbr_'//string1(1:len_string1)//'.dat')
     open(unit=unit_down,status='unknown',file='downw_'//string1(1:len_string1)//'.dat')

     
     ! Initial conditions:
     x_umb=0.0
     y_umb(1)=w_umb(1)
     y_umb(2)=v_umb(1)
     y_umb(3)=f(1)
     do i=1,19
        y_umb(3+i)=m_umb(i,1)
     end do
     xx(1)=x_umb
     
     ! ratio1, ratio2 - to monitor the approaching singularity and
     ! stop execution before the program crashes. ratio2 (current one)
     ! and ratio1 (old one) are used to approximate what the ratio
     ! v*v/w at the next step would be. Linear extrapolation is used.
     ! If the result is close to zero (ratio3=15 is the current cutoff value
     ! but it may be increased if the program crashes prematurely)
     ! then execution is terminated.
     ratio1=0.0
     ratio2=0.0
     
     index=1
     ido=1
     step=0.0
     write (*,'(/4x,"x",10x,"w",11x,"v",11x,"f",15x,"speed")')
     
     !
     ! main integration loop - can we combine with the one above or
     !  are we doomed to do them separately?
     !
     do
        step=step+incr
        xend=step
        index=index+1
        ! call ivmrk (ido,n,deriv,x,xend,y,yprime)
        !  "ivmrk" is an IMSL subroutine for solution of a system
        !  of ordinary differential equations using Runge-Kutta method.
        !  "ido" - flag indicating the state of the computation
        !       "1" for initial entry, "3" - final call to release workspace
        !  "n" - number of equations
        !  "deriv" - name of subroutine to evaluate derivatives
        !  "x" - independent variable
        !  "xend" - value of x where the solution is required
        !  "y" - array of size n of dependent variables
        !  "yprime" - array of size n containing the values of 
        !        derivatives evaluated at (x,y)
        !  See IMSL documentation for more details.
        !print*,'pre-RK4, y_umb = ',y_umb
        CALL RK4(y_umb,yprime,22,x_umb,incr,y_umb,derivs_umb)
        !
        ! (mdj) N.B., RK4 does not advance x_umb by incr
        x_umb = xend
        !
        if ( (step.lt.final).and.(ido.ne.3) ) then
           xx(index)=x_umb
           w_umb(index)=y_umb(1)  ! (mdj) w is volume flux
           v_umb(index)=y_umb(2)
           f(index)=y_umb(3)
           ! b(index) is the radius of the visible edge of plume from MTT.
           b_umb(index)=2.0*y_umb(1)/y_umb(2)
           do i=1,19
              ! use y below to find particles still in plume
              ! yprime for particles that have left the plume
              m_umb(i,index)=y_umb(i+3)  ! m is mass flux for various particle sizes
              ! Use the first m' for kg/m, second for kg/m2 deposition.
              ! mprime(i,index)=-yprime(i+3)
              mprime(i,index)=-yprime(i+3)/2/pi/b_umb(index)
           end do
           
           ratio1=ratio2
           ratio2=y_umb(2)*y_umb(2)/y_umb(1)
           
           write(*,'(f7.1,3e13.5,f15.3)') x_umb,y_umb(1),y_umb(2),y_umb(3),ratio2
           
           if (ratio2 .lt. ratio1) then
              ratio3=(ratio1-ratio2)/incr/ratio2
           else
              ratio3=0.0000001
              ! (Arbitrary safe number)
           endif
           ! Final call to release workspace:
           if ((ratio3 .gt. 0.0007).or.(step .gt. final)) ido=3
           !
           ! was main loop variable
           ! go to 10
        else
           exit ! bail out on main loop
        end if
        ! STOP 'dbx'
     end do ! main integration loop
     
     write (*,*) "Remaining results saved to file mi_"//string1(1:len_string1)//".dat"
     ! Index of the last x for which functions were calculated:
     ! (useful in further calculations as upper limit for loop counters)
     count=index-1
     
     ! Find Hb (when "f" changes sign) and Ht
     ! Hb is the height of the bottom of the cloud (Height when cloud
     ! density = atmospheric density)
     ! Ht is for cloud top height
     ! Rhb - radius at Hb=const
     ! Qhb - flux at Hb
     do i=2,count
        if ((f(i-1) .gt. 0.0) .and. (f(i) .le. 0.0)) then
           Hb=xx(i-1)
           Qhb=w_umb(i-1)*PI
        end if
     end do
     
     ! Execution is stopped *before* the ratio v*v/w approaches zero
     ! to prevent the program crashing, therefore Ht is slightly
     ! underestimated.
     ! To correct that 1/ratio3 is added (linear extrapolation
     ! of final velocity trend to estimated top).
     ! deltah0 is initial cloud thickness for steady release
     Ht=xx(count)+1/ratio3
     deltah0=Ht-Hb
     Rhb=deltah0/2.0
     dens=1.225*exp(-Hb*N1s*N1s/g)
     lamda=0.8
     
     ! Determine which value of M-functions (already calculated) to 
     ! use - it will be one corresponding to Hb, or the first larger than Hb.
     ! "mind" is number of that value in the array
     mind=INT(Hb/incr)+1
     print*,'pre-umbrella, mind = ',mind
     call umbrella
     
     write (*,*) "Results for umbrella cloud saved to file umbr_"//string1(1:len_string1)//".dat"
     write(*,*)"Results for downwind plume saved to file downw_"//string1(1:len_string1)//".dat"
     
     ! To find fallout from column in wafers, use
     ! integration by Simpson's 1/3 rule; value of y' is 
     ! calculated at every 
     ! second x, that is, if output is printed for
     ! x=0, 100, 200, 300, 400, 500, etc then integrals will be
     ! calculated for x=100, 300, 500 etc. These values will correspond
     ! to middle points of intervals: 0-200, 200-400, etc. 
     
     do index=1,19
        do i=2,count-1,2
           a(index,i/2)=(mprime(index,i-1)+4.0*mprime(index,i)+  &
                & mprime(index,i+1) )*2.0*incr/3.0
        end do
     end do
     
     write(unit_mi,*)'The numbers in the m-functions are related to'
     write(unit_mi,*)'the grain size in the following way:'
     write(unit_mi,*)'   fi-size    size in mm     m-function '
     do i=1,19
        write(unit_mi,'(f9.1,f14.4,i12)') scfi(i),scmm(i),i
     end do
     write(unit_mi,'(/)')
     ! The user can select one of two output formats - uncomment 
     ! one of two blocks of instructions:
     !
     ! First output format - masses of particles left in the plume
     ! are printed separately for each size (in kg/s) together with 
     ! results of integration:
     ! do index=1,19
     !    write(unit_mi,*)'    M number',index
     !    write(unit_mi,*)'   x               m        integral'
     !    do i=1,count
     !       if (MOD(i,2).eq.1) then
     !          write(unit_mi,'(f7.1,f15.3)') xx(i),m(index,i)
     !       else
     !          write(unit_mi,'(f7.1,f15.3,f16.1)') xx(i),m(index,i), &
     !                & a(index,i/2)
     !       end if
     !    end do
     !    write(unit_mi,'(/)')
     ! end do
     !
     ! Alternative output format:
     ! To print without integral value in format convenient for
     ! further visualization (to be imported into Matlab or other
     ! graph-plotting program): fallouts for different sizes are
     ! printed together as a matrix, with the first column being
     ! the distance and the rest 19 columns being fallout values 
     ! in kg/m*s
     do i=1,count
        write(unit_mi,'(2f7.1,19(e17.9))') xx(i),b_umb(i), &
             & (mprime(index,i),index=4,8)
     end do
     
     close(unit=unit_mi,status='keep')
     close(unit=unit_umbr,status='keep')
     close(unit=unit_down,status='keep')
     
  end if ! flag_umbrella (strongwind.f)
  
  !  end DO ! b0=100.,1000.,100.
  !      CLOSE THE I/O FILES
  CLOSE(unit_height, STATUS='KEEP')
  CLOSE(unit_in, STATUS='KEEP')
  CLOSE(3, STATUS='KEEP')
!      THIS ENDS THE NUMERICAL INTEGRATION
!(mdj)     STOP  ! using stop for regular termination not so good
!(mdj)  END DO   ! superfluous?  should be end program
end program PUFFIN

! -----------------------------------------C
!      C
!      C
!      SUBROUTINES AND FUNCTIONS          C
!      C
!      C
! -----------------------------------------C
SUBROUTINE RK4(Y,DYDX,N,X,H,YOUT,DERIVS)
  use htvars
  implicit none
!  integer, parameter :: NMAX=26
  integer :: n
  real(kind=sp) :: x,h
!  real(kind=sp) :: Y(N),DYDX(N),YOUT(N),YT(NMAX),DYT(NMAX),DYM(NMAX)
  real(kind=sp) :: Y(N),DYDX(N),YOUT(N),YT(N),DYT(N),DYM(N)
  external derivs

  integer :: i
  real(kind=sp) :: hh,h6,xh
  HH=H*0.5
  H6=H/6.
  XH=X+HH
  DO I=1,N
     YT(I)=Y(I)+HH*DYDX(I)
  end DO
  CALL DERIVS(XH,YT,DYT)
  !print*,'rk4: 1st dyt(1-3)  ',dyt(1),dyt(2),dyt(3)
  !print*,'rk4: 1st yt(1-3)  ',yt(1),yt(2),yt(3)
  DO I=1,N       
     YT(I)=Y(I)+HH*DYT(I)
  end DO
  CALL DERIVS(XH,YT,DYM)
  !print*,'rk4: dym(1-3)  ',dym(1),dym(2),dym(3)
  !print*,'rk4: 2nd yt(1-3)  ',yt(1),yt(2),yt(3)
  DO I=1,N
     YT(I)=Y(I)+H*DYM(I)
     DYM(I)=DYT(I)+DYM(I)
  end do
  CALL DERIVS(X+H,YT,DYT)
  !print*,'rk4: dydx(1-3) ',dydx(1),dydx(2),dydx(3)
  !print*,'rk4: dyt(1-3)  ',dyt(1),dyt(2),dyt(3)
  !print*,'rk4: dym(1-3)  ',dym(1),dym(2),dym(3)
  !print*,'rk4: yt(1-3)  ',yt(1),yt(2),yt(3)
  DO I=1,N
     YOUT(I)=Y(I)+H6*(DYDX(I)+DYT(I)+2.*DYM(I))
  end do
  RETURN  

END subroutine rk4

!      
!      
!      DIFFERENTIAL EQUATIONS      
!      
!      
SUBROUTINE DERIVS(S,Y,DYDX)
  use htvars
  use htpaths
  implicit none

  real(kind=sp) :: s,Y(NVAR),DYDX(NVAR)

  integer :: j,jjj
  real(kind=sp) :: ps,eps,ydy,y1y2,y2y1,emallds,ue,phi2

  real(kind=sp) :: phi,rhoa,ta,vv
  external phi,rhoa,ta,vv

  ! PARAMETERS NOT DEPENDENT ON MI
  !  REENTRAINMENT AREA AND B-V FREQ
  !print*,'derivs: on entry, y(1),y(2),y(3)=',y(1),y(2),y(3)
  PS=0.226*(SIN(Y(3)))**2+0.5*(COS(Y(3)))**2
  EPS=(0.014)**2
  ! CALCULATE THE MASS OF PYROCLASTS
  !  AMOUNT OF FALLOUT
  !
  Y1Y2=SQRT(Y(1)**2/Y(2))
  ! (mdj) at the last step we seem to hit y(2) < 0,
  !  so test for that and set to zero
  !
  !if (y(2).le.0) then
  !   y1y2 = 1.0
  !else
  !   y1y2=SQRT(y(1)**2/y(2))
  !end if
  ! (mdj) end
  Y2Y1=Y(2)/Y(1)
  !print*,'derivs: on entry, y(1),y(2)=',y(1),y(2)
  EMALLDS=0.0
  DO JJJ=5,23
     J=JJJ-4
     PHI2=PHI(LMF,FO,VSO,J)
     DYDX(JJJ)=-2.0*PS*VSO(J)*Y(JJJ)/Y1Y2/Y2Y1  &
          &        +((6*PHI2*EPS*U)/(5*B**2))
     EMALLDS=DYDX(JJJ)+EMALLDS     
     !print*,'jjj,phi2,vso,y,dydx=',jjj,phi2,vso(j),y(jjj),dydx(jjj)
  end DO
  DYDX(24)=EMALLDS
  UE=ALPHA*ABS(U-VV(V,Z,SJ,HJ)*COS(Y(3))) &
       &     +BETA*ABS(VV(V,Z,SJ,HJ)*SIN(Y(3)))  
  DYDX(1)=2*RHOA(Z)*UE*B+EMALLDS
  YDY=DYDX(1) 
  DYDX(2)=-B**2*(RHO-RHOA(Z))*GRAV*SIN(Y(3)) &
       &     +VV(V,Z,SJ,HJ)*COS(Y(3))*YDY+U*EMALLDS
  DYDX(3)=-(B**2*(RHO-RHOA(Z))*GRAV*COS(Y(3)) &
       &     -VV(V,Z,SJ,HJ)*SIN(Y(3))*YDY)/Y(2)
  DYDX(4)=2*UE*RHOA(Z)*B*CA*TA(Z)-B**2*UE*RHOA(Z)*GRAV &
       &     +CP*T*EMALLDS
  !print*,'derivs: z,rhoa(z) = ',z,rhoa(z)
  !print*,'derivs: emallds,b,ue = ',emallds,b,ue
  !print*,'derivs: vso(j),phi2,u=',vso(j),phi2,u
  !print*,'derivs: on exit, y(1),y(2)=',y(1),y(2)
  RETURN
END subroutine derivs

!      '(e30.20)'
!      
!      UPDATE THE PARAMETERS BETWEEN STEPS
!      
!      
SUBROUTINE UPDATE(Y1,Y2,Y4,Y24,U,RHO,EN,RS,RHOG,B, &
     &     W,M,R,RV,RA,C,CV,CA,CP,T)
  use htvars
  use htpaths,only: z
  implicit none
!      UPDATE VALUES OF PARAMETERS

  real(kind=sp) :: y1,y2,y1y2,y4,y24,u,rho,en,rs,rhog,b,w,m, &
       & r,rv,ra,c,cv,ca,cp,t

  real(kind=sp) :: rhoa,ta
  external rhoa,ta

  Y1Y2=SQRT(Y1**2/Y2)
  U=Y2/Y1
  RHO=1/((1-EN)/RS+EN/RHOG)
  B=1/(SQRT(RHO))*Y1Y2
  ! EN is total gas fraction
  EN=(Y1-Y24)/Y1
! W = total - fraction air (ie, cur mass flux - orig)) - fraction particles 
! W = fraction water (bug fix 2010-10-06, 2nd term was -(Y1-M)/Y2
  W=1-(Y1-M)/Y1-Y24/Y1
  R=W/EN*RV+(1-W/EN)*RA
  C=W*CV+(EN-W)*CA+(1-EN)*CP
  T=Y4/C/Y1
  RHOG=RHOA(Z)*RA*TA(Z)/R/T
  RETURN
END SUBROUTINE UPDATE

!      
!      print*,'rk4: 1st dyt(1-3)  ',dyt(1),dyt(2),dyt(3)
!      PARTICLE SETTLING SPEEDS
!      
!      
SUBROUTINE FALLSPEED(PV)
  use htvars
  use htpaths
  implicit none

  real(kind=sp) :: PV(NPAR)
  REAL(KIND=SP) :: SCFI(NPAR),SCMM(NPAR)
  integer :: iv,index,i
  real(kind=sp) :: fi,di,diff,i2di,cd,cdold,enu,re,sv

17 FORMAT(5A10)
!      MEMORIZES THE PHI AND MM SCALE
  INDEX=1
!  DO 15 FI=-8,10,1.0
  DO I=-8,10,1
     fi=1.0_sp*i
     DI=2.0**(-FI)
     SCMM(INDEX)=DI
     SCFI(INDEX)=FI
     SV=SQRT(2*PI)
     I2DI=LOG10(DI)/LOG10(2.0)
     PV(INDEX)=EXP(-((FI-MU)**2)/(2.0*SIGMA**2))/(SV*SIGMA)
!      CALCULATE DENSITY AND VISCOSITY OF AIR AT VENT HEIGHT
!      SEA LEVEL DENSITY FROM NASA 1962 STANDARD ATMOSPHERE
!      VISCOSITY FROM MIDDLETON & WILCOCK
!     ENU=0.000018/RHOG
!      That doesn't work to well, use enu from Atham results:
     enu=0.37/rhog
!      GET FALL SPEEDS
     CD=1.0
     DO IV = 1, 20
        VSO(INDEX)=SQRT(4.0*(0.001*DI)*GRAV*(RS-RHO)/(3.0*RHO*CD))
        RE=VSO(INDEX)*(0.001*DI)/ENU
        CDOLD=CD
        CD=(24.0/RE)*(FF**(-0.32))+2*SQRT(1.07-FF)
        DIFF=ABS((CDOLD-CD)/CD)
!      USE THIS OUTPUT TO CHECK THAT CD IS CONVERGING
!      WRITE(*,*) 'DIFF OF ', INDEX, ' = ', DIFF
!14   CONTINUE
     end do
     IF (DIFF.GE.0.01) THEN
        WRITE(*,*)'ERROR IN FALL SPEED (INDEX=', INDEX,') >1%'
     ENDIF
     INDEX=INDEX+1
!15 CONTINUE
  end DO
!      WRITE(*,20)(RHOA(Z),I-9,VSO(I),I=1,19)
!      20      FORMAT(F10.3,I10,F20.5)
  RETURN
END subroutine fallspeed

!      
!      
!      CORRECT PARTICLE FALLOUT TO EDGE OF PLUME
!      
!      
SUBROUTINE FALLOUTCOR(X,Y,DXF,B,XFF,Z,H)
  use htvars
  use htpaths,only: nvar,npar,v,sj,hj,vso
  implicit none

  REAL(KIND=SP) :: DXF(NVAR),Y(NVAR),XFF(NVAR)
  real(kind=sp) :: x,b,z,h
  real(kind=sp) :: xf,zf
  real(kind=sp) :: vv
  external vv
  integer :: j

!  print*,'falloutcor - vso = ',vso
!      TRANSFORMS X AND Z POSTION TO THE PLUME EDGE
  XF=X+SIN(Y(3))*(B)
  ZF=Z-COS(Y(3))*(B)
!      INITIALIZES THE ARRAY DXF(J) (X POS FOR EACH GRAIN SIZE)
  DO J=1,19
50   DXF(J)=XF
  end DO
!      DXF(INDEX) IS THE CORRECTION WHICH
!      ACCOUNTS FOR HORIZONTAL MOVEMENT CAUSED BY THE WIND AFTER
!      THE PARTICLE LEAVES THE PLUME
!      DO WHILE THE PARTICLE HAS NOT HIT THE GROUND
  DO WHILE (ZF.GT.0)
!      CALCULATE THE DISTANCE EACH PARTICLE SIZE HAS MOVED
!      IN THIS STEP
     DO J=1,19
        !print*,'j,vso(j),vv = ',j,vso(j),VV(V,ZF,SJ,HJ)
        DXF(J)=DXF(J)+(H*VV(V,ZF,SJ,HJ))/(VSO(J))
        XFF(J)=DXF(J)
        IF (ZF.GT.10000) THEN
        !WRITE(*,*) J,ZF,VV(V,ZF,SJ,HJ),VSO(J),XFF(J)
        ENDIF
!100  CONTINUE
     end DO
!      CHANGE IN Z FOR THIS STEP
     ZF=ZF-H
!400 CONTINUE
  end DO
  RETURN
END subroutine falloutcor

!      
!      
!      TO CALCULATE A 2-D PARTICLE FIELD
!      
!      

SUBROUTINE FIELD(XC,BF,ZC,UF,THET,EMASS)
  use htvars
  implicit none
!      SUBROUTINE TO ESTIMATE INSTANTANEOUS FIELD OF PARTICLE MASS
!      WITHIN PLUME
  real(kind=sp) :: xc,bf,zc,uf,thet,emass
  integer :: i,ixf
  real(kind=sp) :: xf,ri,ca,tca,xi

  XF=BF
  IXF=INT(XF)
  DO I=-20000,20000,200
     RI=FLOAT(I)
     CA=(EMASS/UF)*EXP(-RI**2/(XF/2.)**2)/3.1416/XF**2
!     print*,'field: ca,emass,uf = ',ca,emass,uf
     TCA=LOG10(CA)
     XI=XC+RI
     WRITE(3,*)RI,ZC,CA,TCA
  end DO
  WRITE(3,*)

  RETURN
END subroutine field

!      
!      
!      ATMOSPHERIC TEMPERATURE
!      
!      
FUNCTION TA(Z)
  use htvars
  use htdata
  implicit none
  real(kind=sp) :: ta
  real(kind=sp) :: z

  if (flag_interp) then ! (mdj) we have experimental data to use
     ! index of height in height array should be set already from main routine
     ! index_height = bsearch(z,dat_h,ndata)
     if (z<dat_h(1)) then ! set to first value (crude)
        ta=dat_t(1)
     else if (z>dat_h(ndata)) then ! set to last value (also crude)
        ta=dat_t(ndata)
     else ! linearly interpolate
        ta = dat_t(index_height)+(z-dat_h(index_height))*(dat_t(index_height+1)-dat_t(index_height))/ &
             & (dat_h(index_height+1)-dat_h(index_height))
     end if
  else
!      FROM   KALEIDAGRAPH BEST FIT TO ATM FROM SPARKSETAL97
!      Y = M0 + M1*X + ... M8*X^8 + M9*X^9
!      TA=279.81-0.011335*Z+5.3298E-07*Z**2
!      &  -9.0085E-12*Z**3+5.4568E-17*Z**4
!      NEXT IS TEMP FOR ISOTHERMAL ATM WITH H=8500. M (WALLACE&HOBBS)
     TA=288.
  end if
  RETURN
END FUNCTION TA

!      
!      
!      ATMOSPHERIC DENSITY
!
FUNCTION RHOA(Z)
  use htvars
  use htpaths, only: nn,grav
  use htdata
  implicit none
  real(kind=sp) :: rhoa
  real(kind=sp) :: z,rhoa0

  if (flag_interp) then ! (mdj) we have experimental data to use
     ! index of height in height array should be set already from main routine
     ! index = bsearch(z,dat_h,ndata)
     
     ! fit to Apr24 El. eruption daa for rho, 1.5557 * exp(-0.00014536 * z)
     !  but overall fit is very poor, especially for low height
     !
     if (z<dat_h(1)) then ! set to first value (crude)
        rhoa=dat_rho(1)
     else if (z>dat_h(ndata)) then ! set to last value (also crude)
        rhoa=dat_rho(ndata)
!      Or really it should decrease exponentially from the last value
        rhoa=dat_rho(ndata)*exp(-Z*Nn/grav)
     else ! linearly interpolate
        rhoa = dat_rho(index_height)+(z-dat_h(index_height))*(dat_rho(index_height+1)-dat_rho(index_height))/ &
           & (dat_h(index_height+1)-dat_h(index_height))
     end if
  else
!      COMMON /PATH1/ PI,GRAV,RV,RA,CV,CA,CP,M,RS,RHOG,MU,SIGMA,
!      $ RHO0,FF,H,T0,EN,V,SJ,HJ,ALPHA,BETA,EMALL,U,B,RHO,T,LMF,FO
!      RHOA0 IS SEA-LEVEL AIR DENSITY
     RHOA0=1.225
!      RHOA=1.3
!      NEXT IS BASED ON TEMPERATURE STRUCTURE FROM FIT TO SPARKSETAL97
!      RHOA=RHOA0*TA(0.0)/TA(Z)*EXP(-GRAV/RA*54.3159*ATAN(0.0000143345
!      & *(-151262.+2*Z))+90.0861*ATAN(0.000019015*(-13826.1+2*Z))-
!      & 31.0212*LOG(6.93669E+09-151262.*Z+Z**2)-
!      & 31.0212*LOG(7.39218E+08-13826.1*Z+Z**2))
!     RHOA=RHOA0*EXP(-Z/8500.)
     rhoA=rhoa0*exp(-Z*Nn/grav)
  end if
  RETURN
END FUNCTION RHOA

!      
!      
!      WIND SPEED AS A FUNCTION OF HEIGHT
!      
!      
FUNCTION VV(V,Z,SJ,HJ)
  use htvars
  use htdata,only: index_height,ndata,dat_h,dat_v
  implicit none
  real(kind=sp) :: vv
  real(kind=sp) :: v,z,sj,hj
  if (flag_interp) then ! (mdj) we have experimental data to use
     ! index of height in height array should be set already from main routine
     ! index = binarysearch(dat_h,z,ndata)
     if (z<dat_h(1)) then ! set to first value (crude)
        vv=dat_v(1)
     else if (z>dat_h(ndata)) then ! set to last value (also crude)
        vv=dat_v(ndata)
     else ! linearly interpolate
        vv = dat_v(index_height)+(z-dat_h(index_height))*(dat_v(index_height+1)-dat_v(index_height))/ &
             & (dat_h(index_height+1)-dat_h(index_height))
     end if
  else
!      VV=V
     VV=V*EXP(-(Z-HJ)**2./SJ**2.)
  endif
  RETURN
END FUNCTION VV

!      
!      
!      RE-ENTRAINMENT PARAMETER
!      
!      
FUNCTION PHI(LMF,FO,VSO,I)
  use htvars
  implicit none

  REAL(KIND=SP) :: VSO(NPAR),PPHI,LMF
  real(kind=sp) :: phi

  integer :: i
  real(kind=sp) :: FO
  PPHI=(FO/LMF)**(1/3)/VSO(I)
  PHI=0.43/(1.+(0.78/PPHI)**6.)
  RETURN
END FUNCTION PHI

!
! new routines (mdj)
!
subroutine getdata
 use htvars
 use htdata
 implicit none

 integer,parameter :: datunit=33,MAXTERMS=20,MAXTERMLEN=30
 integer :: i,j,numlines,rawlines,len_cdata,istart,iterm,nterms,lenw
 character(LEN=MAXLINELEN) :: cdata
 character(LEN=MAXTERMLEN) :: cterms(MAXTERMS)
 logical flag_file,newterm
 !
 ! open data file, parse
 !
 dataname = TRIM(dataname)
 len_dataname = LEN_TRIM(dataname)
 if (len_dataname.le.0) then
    write(*,*) 'dataname has length <=0, ',len_dataname
    STOP 'bad length on data file name.'
 end if
 INQUIRE(file=dataname(1:len_dataname),exist=flag_file)
 if (.not.flag_file) then
    write(*,*) 'Unable to open data file: ',dataname(1:len_dataname)
    STOP 'Error in opening data file.'
 end if
 OPEN(file=dataname(1:len_dataname),status='old',unit=datunit)
 !
 ! parse file once, count significant lines
 !
 numlines=0
 do
    read(unit=datunit,fmt='(a200)',end=2) cdata
!    cdata = TRIM(cdata)
!    len_cdata = LEN(cdata)
    do j=1,MAXLINELEN ! shift leading blanks to end of string
       if (cdata(1:1).eq.' ') then
          do i=2,MAXLINELEN
             cdata(i-1:i-1)=cdata(i:i)
          end do
          cdata(MAXLINELEN:MAXLINELEN)=' '
       else
          exit
       endif
    end do
    do i=MAXLINELEN,1,-1
       if (cdata(i:i).ne.' ') exit
    enddo
    len_cdata=i
!    print*,'cdata = ',cdata
!    print*,'line,len(cdata) = ',numlines,len_cdata
    if ((len_cdata.gt.0).and.((cdata(1:1).ne.'#').and.(cdata(1:1).ne.'%'))) then
       numlines = numlines+1
    end if
 end do
2 continue
 write(*,*) 'Found ',numlines,' significant data points in ', &
      & dataname(1:len_dataname)
 ndata=numlines
 ALLOCATE(dat_h(ndata),dat_t(ndata),dat_v(ndata),dat_rho(ndata))
 !
 ! parse again, load data arrays
 !
 REWIND(datunit)
 numlines=0
 rawlines=0
 do
    cdata=""
    read(unit=datunit,fmt='(a200)',end=3) cdata
    rawlines=rawlines+1
!    cdata = TRIM(cdata)
!    len_cdata = LEN(cdata)  ! not sure why these are not doing the job ...
    do j=1,MAXLINELEN ! shift leading blanks to end of string
       if (cdata(1:1).eq.' ') then
          do i=2,MAXLINELEN
             cdata(i-1:i-1)=cdata(i:i)
          end do
          cdata(MAXLINELEN:MAXLINELEN)=' '
       else
          exit
       endif
    end do
    do i=MAXLINELEN,1,-1
       if (cdata(i:i).ne.' ') exit
    enddo
    len_cdata=i
    !print*,'cdata = ',cdata(1:len_cdata)
    !print*,'len(cdata) = ',len_cdata
    if ((len_cdata.gt.0).and.((cdata(1:1).ne.'#').and.(cdata(1:1).ne.'%'))) then
       numlines=numlines+1
       ! significant data
       !
       ! break into terms
       i=0
       istart=1
       iterm=1
       newterm=.false.
       do
          i = i+1
          if (i.gt.len_cdata+1) exit ! pad by 1 to avoid losing last character
          !print*,'i,cdata(i:i),char: ',i,cdata(i:i),IACHAR(cdata(i:i))
          if (cdata(i:i).eq.' '.or.cdata(i:i).eq.','.or.cdata(i:i).eq.char(9).or.i.eq.len_cdata+1) then 
             ! space , or tab delimited terms in the overall string
             if (.not.newterm) then
                lenw = i-istart
                if (lenw.gt.MAXTERMLEN) then
                   write(*,*) 'MAXTERMLEN exceeded. ',lenw,MAXTERMLEN,cdata(istart:i)
                   stop
                endif
                cterms(iterm) = cdata(istart:i-1)
                newterm=.true.
                !write(*,*) '<getdata> term found, iterm, term = ',iterm,lenw,cterms(iterm),cdata(istart:i-1)
             endif
          else
             if (newterm) then ! new term
                istart=i
                newterm=.false.
                iterm=iterm+1
                if(iterm.gt.MAXTERMS) then
                   write(*,*) 'MAXTERMS exceeded.',iterm,MAXTERMS
                   stop
                endif
             endif
          endif
       enddo
       nterms=iterm
       !write(*,*) '<getdata>: line parsed.  ',nterms,' terms found.'
       !do j=1,nterms
       !   print*,'j,cterms(j) = ',j,cterms(j)
       !end do
       if (nterms.lt.col_rho) then
          write(*,*) 'nterms = ',nterms,' < col_rho = ',col_rho,' on line ',rawlines
          STOP '<getdata> fewer columns than expected.'
       end if
       !
       ! now use input column pointers to load data
       !
       read(unit=cterms(col_h),fmt=*,err=4) dat_h(numlines)
       read(unit=cterms(col_t),fmt=*,err=4) dat_t(numlines)
       read(unit=cterms(col_v),fmt=*,err=4) dat_v(numlines)
       read(unit=cterms(col_rho),fmt=*,err=4) dat_rho(numlines)
       !print*,'line,h,t,v,rho: ',numlines,dat_h(numlines),dat_t(numlines),dat_v(numlines),dat_rho(numlines)
       ! any necessary conversions
       !  temp from C to K
       !  rho?  kg/m^3
       dat_t(numlines) = dat_t(numlines)+273.15
    end if ! significant line of data
 end do
3 continue
 write(*,*) 'Data file parsed, read ',numlines,' of total ',rawlines,' lines.'
 ndata=numlines
 write(*,'(" Height range in data [m]: ",2e12.4)') dat_h(1),dat_h(ndata)
 return
4 write(*,*) 'Error on line ',rawlines
 STOP 'Error in parsing data file.'

end subroutine getdata
  
integer FUNCTION bsearch(a,x,n)
  ! search for value a in x(:)
  USE htvars
  IMPLICIT NONE
  integer, intent(IN) :: n
  REAL(kind=SP), INTENT(IN) :: x(n)
  REAL(kind=SP), INTENT(IN) :: a

!  INTEGER :: bsearch
  INTEGER :: jl,jm,ju
  LOGICAL :: ascend
!  n=size(x)
!  print*,'bsearch: n,a,x(1),x(n): ',n,a,x(1),x(n)
  ascend = (x(n) >= x(1))
  jl=0
  ju=n+1
  do
     if (ju-jl <= 1) exit
     jm=(ju+jl)/2
     if (ascend .eqv. (a >= x(jm))) then
        jl=jm
     else
        ju=jm
     end if
  end do
  if (a == x(1)) then
     bsearch=1
  else if (a == x(n)) then
     bsearch=n-1
  else
     bsearch=jl
  end if
  !
  ! hack for height less than 1st data point - fixes values
  !  below that to first available data value rather than
  !  extrapolate
  !
  if (bsearch.le.0) then
     bsearch=1
  end if
END FUNCTION bsearch

!subroutine umbrella(Qhb,Rhb,term,deltah0,dens,lamda,mind,wind,eKh,distto)
subroutine umbrella
  use htvars
  use htpaths
  use umb_mod
  implicit none

  !real(kind=sp) :: Qhb,Rhb,term,deltah0,dens,lamda !,Ht,Hb
  !real n1,n1s,n2,ps,eKh,Vso(19),m(19,200)
  real time(100), R(100), deltaH(100),S(19,16),rcur(16)
  real(kind=sp) :: Vr(100),Sdw(19,16),Mdw(19,16)
  real(kind=sp) :: xcur(16),Mini(19),Mum(19,16),SdwN(19,16)
  real(kind=sp) :: timedw(16),width(16),thick(16),Xxx(16)
  
  integer :: i,k,l,kend
  real(kind=sp) :: temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8
  real(kind=sp) :: rstep,tstep,v2
  !common N1,N1s,N2,ps,Vso
  !common /um/ m,Ht,Hb

  ! Set up time stepping for umbrella cloud------------
  tstep=500.0
  kend=INT(term/tstep)+1

  do k=1,kend

     time(k)=tstep*(k-1)
     ! -------Calculation of umbrella spreading---------------
     !  Two possible ways to calculate R(k).  Second seems to work better
     !  wrt data.
     !       R(k)=(Rhb**3.0+3.0*lamda*N1s*Qhb*time(k)*time(k)/6.28)**0.3334

     R(k)=(Rhb**3.0+lamda*N1s*254.0*(Ht/1000.0)**5.26*time(k)*time(k))**0.3334

     if (k.gt.1) then
        Vr(k)=(R(k)-R(k-1))/tstep
        if (Vr(k).le.wind) then 
           !go to 204
           exit ! bail on do k=1,kend
        end if
     end if
     ! For model of decreasing cloud thickness with time.  Seems to yield
     ! unrealistic fallout patterns.  Constant cloud thickness is therefore
     ! assumed.  Uncomment this to see decreasing thicknesses, if necessary.
     ! deltaH(k)=(Qhb*time(k) + 3.14*Rhb*Rhb*deltah0*dens)/(3.14*R(k)*R(k)*dens)
     deltaH(k)=deltaH0

!200  Rs=R(k)
     Rs=R(k)
  end do
!204 kend=k
!  kend=k

  write(unit_umbr,*)'Stagnation point=', Rs
  write(unit_umbr,*)'Ht= ',Ht,'Hb= ',Hb,'Rhb= ',Rhb
  write(unit_umbr,*)'     Time          R             Vr'
  do k=1,kend
     ! write(9,*)R(k)
     write(unit_umbr,'(f11.1,f15.1,f12.2)') time(k),R(k),Vr(k)
  end do

  ! The following routine calculates sedimentation from
  ! umbrella cloud (fallout per unit area) from Rhb up to
  ! maximum R for a given time (calculated above) in 15 increments
  ! (which gives us 16 points).  
  ! Current R is called rcur(l), l=1...16
  ! Number of increments can be changed.
  
  ! Vsi(i)=Vso(i)*exp(x*N*N/(2*g)). Exponent is calculated here
  ! and is called v2
  v2=exp(Hb*N2/19.6)

  do k=2,kend
     write(unit_umbr,*) ' Time =',time(k)
     write(unit_umbr,*) ' Maximum radius for this time is ',R(k), &
          & '     deltaH=',deltaH(k)
     write(unit_umbr,*) ' Printing fallout values for points between ', &
          & 'Rhb and maximum R'
     rstep=(R(k)-Rhb)/15.0
     write(unit_umbr,*) 'rstep =',rstep
     ! rcur(l) are radial distances for which fallout from umbrella
     !  cloud should be calculated.
     !  rcur are different for every time step.
     do l=1,16
        rcur(l)=Rhb+rstep*(l-1)
     end do

     ! ---------Fallout for radially spreading umbrella
     !write(unit_umbr,'("(mdj) Qhb,mind = ",e13.4,i4)') Qhb,mind
     !write(unit_umbr,'("(mdj) m_umb(i,mind):",/4(5e13.4/))') (m_umb(i,mind),i=1,19)
     do l=1,16
        do i=1,19
           ! temp1, temp2, temp3 - temporary variables that hold different
           ! parts of formula for fallout, which can also be calculated 
           ! in one step 
           temp1= -3.14*Vso(i)*v2*(rcur(l)*rcur(l)-Rhb*Rhb)/Qhb
           ! pow=(deltaH(k)/deltah0)**temp1
           temp2=exp(temp1)
           temp3=m_umb(i,mind)*Vso(i)*v2/Qhb

           ! ***IMPORTANT*** Different versions of the formula exist.
           ! Select and uncomment one of them.
           Mum(i,l)=m_umb(i,mind)*temp2
           !        S(i,l)=temp3*pow*temp2
           S(i,l)=temp3*temp2
           Mini(i)=Mum(i,l) ! name is likely a problem here
        end do
     end do

     ! Primary format of output - prints all data
     !     do i=1,19
     !        write(unit_umbr,'(16f14.3)')(rcur(l),l=1,16)
     !        write(unit_umbr,'(16e14.6)') i,(S(i,l),l=1,16)
     !     end do
     !     write(unit_umbr,*)
     ! Alternative format of output - for 
     ! subsequent plotting of selected curves only:
     write(unit_umbr,'(16f14.3)')(rcur(l),l=1,16)
     write(unit_umbr,'(16f14.6)') (S(1,l),l=1,16)
     write(unit_umbr,'(16f14.6)') (S(3,l),l=1,16)
     write(unit_umbr,'(16f14.6)') (S(5,l),l=1,16)
     write(unit_umbr,'(16f14.6)') (S(7,l),l=1,16)
     write(unit_umbr,'(16f14.6)') (S(9,l),l=1,16)
     write(unit_umbr,'(16f14.6)') (S(11,l),l=1,16)
     write(unit_umbr,*)

  end do ! k=2,kend

  ! Check to see whether a downwind plume should be calculated.
  ! if (distto.le.Rs) then
  !     write(*,*)'Downwind plume cannot be calculated'
  !     Go to 304
  ! end if
  if (distto.gt.Rs) then
     !-----------This part of the subroutine handles the 
     !               downwind plume.
     !  
     !  xcur(l) are downwind distances for which fallout from downwind plume
     !   should be calculated.
     do l=1,16
        xcur(l)=Rs+(distto-Rs)*(l-1)/15.0
     end do
     
     ! ----------Fallout for downwind plume
     do l=1,16
        do i=1,19
           ! temp4, temp5, temp6, temp7 - temporary variables that hold different
           ! parts of formula for fallout from downwind plume, which can also be
           !  calculated in one step
           temp4=(2.0*lamda*N1s*Qhb+4.0*eKh)**0.5/wind
           temp5= -2.0*Vso(i)*v2*(xcur(l)**1.5-Rs**1.5)/3.0/Qhb
           temp6=Mini(i)*EXP(temp4*temp5)
           temp7=temp6*Vso(i)*v2/Qhb
           temp8=(xcur(l)/1000.0)**1.5-(Rs/1000.0)**1.5
           
           timedw(l)=time(kend)+(xcur(l)-Rs)/wind
           width(l)=temp4*(xcur(l)+2.0*Rs)**0.5
           thick(l)=Qhb/wind/width(l)
           Xxx(l)=temp8
           Mdw(i,l)=temp6
           Sdw(i,l)=temp7
           SdwN(i,l)=temp7/Sdw(i,1)
        end do
     end do
     
     ! Output, as examples above, should be adapted for 
     !  particular uses.
     
     write(unit_down,'("(mdj) rs,wind = ",2e14.4)') rs,wind
     write(unit_down,'("(mdj) timedw = ",/2(8e12.4/))') (timedw(l),l=1,16)
     write(unit_down,*)'   t    x            X           S'
     write(unit_down,'(3f12.1,4e12.3)')(timedw(l),xcur(l),width(l), &
          &  SdwN(14,l),SdwN(15,l),SdwN(16,l),SdwN(17,l),l=1,16)
     ! write(*,*)(Vso(i),Vso(i)*v2,i=1,19)
  else
     write(*,*)'Downwind plume cannot be calculated'
  end if
  !304 return
  return
end subroutine umbrella

subroutine derivs_umb(x,y,yprime)
  ! adapted directly from strongwind.f
  use htvars
  use htpaths,only: Vso
  use umb_mod,only: N2,ps
  implicit none
  
  real :: x,y(*),yprime(*)

  integer :: i
  real :: alfa,temp2

  alfa=0.093
  temp2=exp(x*N2/19.6)

  yprime(1)=alfa*y(2)*2.0
  yprime(2)=y(1)*y(3)/(y(2)**3)
  yprime(3)=N2*y(1)*(-2.0)
  do i=4,22
     yprime(i)=ps*(-2.0)*Vso(i-3)*temp2*y(i)/y(2)
  end do
  return
end subroutine derivs_umb

subroutine init_random_seed()
  ! (mdj) getting correlated seeds on near-simultaneous runs,
  !  moving to seeding from milliseconds
  integer :: clock,n,idate(8)
  integer, dimension(:), allocatable :: iseed

  call random_seed(size=n)
  allocate(iseed(n))

  !call system_clock(count=clock)
  !seed=clock
  CALL RANDOM_SEED(get=iseed)
  CALL DATE_AND_TIME(values=idate)
  iseed = iseed*(idate(5)+idate(6)+idate(7)+idate(8)-500)     ! 8 is ms, 7 is sec, 6 is m, 5 is h
  write(*,'(" Seeding random numbers with seed = ",4i12,/36x,4i12)') iseed
  write(*,'("  N.B., seed values determined from date_and_time,")')
  write(*,'("  can be correlated for multiple runs with similar start times.")')
  call random_seed(put=iseed)

  deallocate(iseed)
end subroutine init_random_seed

