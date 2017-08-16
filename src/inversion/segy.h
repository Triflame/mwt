c* segy.h - include file for SEGY traces
c*****************************************
c Byte
c 1-4	long tracl;	/* trace sequence number within line */

c 5-8	long tracr;	/* trace sequence number within reel */

c 9-12	long fldr;	/* field record number */

c13-16	long tracf;	/* trace number within field record */

c17-20	long ep;	/* energy source point number */

c21-24	long cdp;	/* CDP ensemble number */

c25-28	long cdpt;	/* trace number within CDP ensemble */

	integer*4 tracl,tracr,fldr,tracf,ep,cdp,cdpt

c29-30	short trid;	/* trace identification code:
c			1 = seismic data
c			2 = dead
c			3 = dummy
c			4 = time break
c			5 = uphole
c			6 = sweep
c			7 = timing
c			8 = water break
c			9---, N = optional use (N = 32,767)
c
c			Following are CWP id flags:
c
c			 9 = autocorrelation
c
c			10 = Fourier transformed - no packing
c			     xr[0],xi[0], ..., xr[N-1],xi[N-1]
c
c			11 = Fourier transformed - unpacked Nyquist
c			     xr[0],xi[0],...,xr[N/2],xi[N/2]
c
c			12 = Fourier transformed - packed Nyquist
c	 		     even N:
c			     xr[0],xr[N/2],xr[1],xi[1], ...,
c				xr[N/2 -1],xi[N/2 -1]
c				(note the exceptional second entry)
c			     odd N:
c			     xr[0],xr[(N-1)/2],xr[1],xi[1], ...,
c				xr[(N-1)/2 -1],xi[(N-1)/2 -1],xi[(N-1)/2]
c				(note the exceptional second & last entries)
c
c			13 = Complex signal in the time domain
c			     xr[0],xi[0], ..., xr[N-1],xi[N-1]
c
c			14 = Fourier transformed - amplitude/phase
c			     a[0],p[0], ..., a[N-1],p[N-1]
c
c			15 = Complex time signal - amplitude/phase
c			     a[0],p[0], ..., a[N-1],p[N-1]
c
c			16 = Real part of complex trace from 0 to Nyquist
c
c			17 = Imag part of complex trace from 0 to Nyquist
c
c			18 = Amplitude of complex trace from 0 to Nyquist
c
c			19 = Phase of complex trace from 0 to Nyquist
c
c			21 = Wavenumber time domain (k-t)
c
c			22 = Wavenumber frequency (k-omega)
c
c			30 = Depth-Range (z-x) traces
c
c			101 = Seismic data packed to bytes (by supack1)
c			
c			102 = Seismic data packed to 2 bytes (by supack2)
c			*/
c
c31-32	short nvs;	/* number of vertically summed traces (see vscode
c			   in bhed structure) */

c33-34	short nhs;	/* number of horizontally summed traces (see vscode
c			   in bhed structure) */

c35-36	short duse;	/* data use:
c				1 = production
c				2 = test */

	integer*2 trid,nvs,nhs,duse 

c37-40	long offset;	/* distance from source point to receiver
c			   group (negative if opposite to direction
c			   in which the line was shot) */

c41-44	long gelev;	/* receiver group elevation from sea level
c			   (above sea level is positive) */

c45-48	long selev;	/* source elevation from sea level
c			   (above sea level is positive) */

c49-52	long sdepth;	/* source depth (positive) */

c53-56	long gdel;	/* datum elevation at receiver group */

c57-60	long sdel;	/* datum elevation at source */

c61-64	long swdep;	/* water depth at source */

c65-68	long gwdep;	/* water depth at receiver group */

	integer*4 offset,gelev,selev,sdepth,gdel,sdel,swdep,gwdep

c69-70	short scalel;	/* scale factor for previous 7 entries
c			   with value plus or minus 10 to the
c			   power 0, 1, 2, 3, or 4 (if positive,
c			   multiply, if negative divide) */

c71-72	short scalco;	/* scale factor for next 4 entries
c			   with value plus or minus 10 to the
c			   power 0, 1, 2, 3, or 4 (if positive,
c			   multiply, if negative divide) */

	integer*2 scalel,scalco

c73-76	long  sx;	/* X source coordinate */

c77-80	long  sy;	/* Y source coordinate */

c81-84	long  gx;	/* X group coordinate */

c85-88	long  gy;	/* Y source coordinate */

	integer*4 sx,sy,gx,gy

c89-90	short counit;	/* coordinate units code:
c				for previoius four entries
c				1 = length (meters or feet)
c				2 = seconds of arc (in this case, the
c				X values are longitude and the Y values
c				are latitude, a positive value designates
c				the number of seconds east of Greenwich
c				or north of the equator */

c91-92	short wevel;	/* weathering velocity */

c93-94	short swevel;	/* subweathering velocity */

c95-96	short sut;	/* uphole time at source */

c97-98	short gut;	/* uphole time at receiver group */

c99-100	short sstat;	/* source static correction */


c101-102  short gstat;	/* group static correction */

c103-104  short tstat;	/* total static applied */

c105-106  short laga;	/* lag time A, time in ms between end of 240-
c			   byte trace identification header and time
c			   break, positive if time break occurs after
c			   end of header, time break is defined as
c			   the initiation pulse which maybe recorded
c			   on an auxiliary trace or as otherwise
c			   specified by the recording system */

c107-108  short lagb;	/* lag time B, time in ms between the time break
c		   and the initiation time of the energy source,
c		   may be positive or negative */

c109-110  short delrt;	/* delay recording time, time in ms between
c			   initiation time of energy source and time
c			   when recording of data samples begins
c			   (for deep water work if recording does not
c			   start at zero time) */

c111-112  short muts;	/* mute time--start */

c113-114  short mute;	/* mute time--end */

c115-116  unsigned short ns;	/* number of samples in this trace */

c117-118  unsigned short dt;	/* sample interval; in micro-seconds */

	integer*2 counit,wevel,swevel,sut,gut,sstat,gstat, 
     1            tstat,laga,lagb,delrt,muts,mute,ns,dt 

c119-120  short gain;	/* gain type of field instruments code:
c				1 = fixed
c				2 = binary
c				3 = floating point
c				4 ---- N = optional use */

c121-122  short igc;	/* instrument gain constant */

c123-124  short igi;	/* instrument early or initial gain */

c125-126  short corr;	/* correlated:
c				1 = no
c				2 = yes */

c127-128  short sfs;	/* sweep frequency at start */

c129-130  short sfe;	/* sweep frequency at end */

c131-132  short slen;	/* sweep length in ms */

c133-134  short styp;	/* sweep type code:
c				1 = linear
c				2 = cos-squared
c				3 = other */

c135-136  short stas;	/* sweep trace length at start in ms */

c137-138  short stae;	/* sweep trace length at end in ms */

c139-140  short tatyp;	/* taper type: 1=linear, 2=cos^2, 3=other */

c141-142  short afilf;	/* alias filter frequency if used */

c143-144  short afils;	/* alias filter slope */

c145-146  short nofilf;	/* notch filter frequency if used */

c147-148  short nofils;	/* notch filter slope */

c149-150  short lcf;	/* low cut frequency if used */

c151-152  short hcf;	/* high cut frequncy if used */

c153-154  short lcs;	/* low cut slope */

c155-156  short hcs;	/* high cut slope */

c157-158  short year;	/* year data recorded */

c159-160  short day;	/* day of year */

c161-162  short hour;	/* hour of day (24 hour clock) */

c163-164  short minute;	/* minute of hour */

c165-166  short sec;	/* second of minute */

c167-168  short timbas;	/* time basis code:
c				1 = local
c				2 = GMT
c				3 = other */

c169-170  short trwf;	/* trace weighting factor, defined as 1/2^N
c			   volts for the least sigificant bit */

c171-172  short grnors;	/* geophone group number of roll switch
c			   position one */

c173-174  short grnofr;	/* geophone group number of trace one within
c			   original field record */

c175-176  short grnlof;	/* geophone group number of last trace within
c			   original field record */

c177-178   short gaps;	/* gap size (total number of groups dropped) */

c179-180   short otrav;	/* overtravel taper code:
c				1 = down (or behind)
c				2 = up (or ahead) */

	integer*2 gain,igc,igi,corr,sfs,sfe,slen,styp,stas,
     1            stae,tatyp,afilf,afils,nofilf,nofils,lcf,hcf,
     1            lcs,hcs,year,day,hour,minute,sec,timbas,
     1            trwf,grnors,grnofr,grnlof,gaps,otrav

c	/* local assignments */
c181-184  float d1;	/* sample spacing for non-seismic data */

c185-188  float f1;	/* first sample location for non-seismic data */

c189-192  float d2;	/* sample spacing between traces */

c193-196  float f2;	/* first trace location */

c197-200  float ungpow;	/* negative of power used for dynamic
c	   range compression */

c201-204  float unscale;	/* reciprocal of scaling factor to normalize
c	   range */

	real*4 d1,f1,d2,f2,ungpow,unscale

c205-206  short mark;	/* mark selected traces */

c207-240  short unass[17];	/* unassigned--NOTE: last entry causes 
c			   a break in the word alignment, if we REALLY
c			   want to maintain 240 bytes, the following
c			   entry should be an odd number of short/UINT2
c			   OR do the insertion above the "mark" keyword
c			   entry */
	integer*2 mark,unass(17)	
