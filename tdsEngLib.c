/************************************************************/
/*															*/
/*	This DLL utility library 'tdsEngLib.dll' can be used	*/
/*  to execute various function from any language			*/
/*  that supports dynamic link  libraries.					*/
/*  Functions usually return an error number.				*/
/*  Error number of 0 means no error						*/
/*															*/
/*	Written by	: Nuri Cankurt								*/
/*	Version		: 1.0  Original version						*/
/*	Date		: March 13, 2004							*/
/*															*/
/*		2.00    Added new functions: polyValue, linearLookup*/
/*				Nuri Cankurt January 22, 2006				*/
/*															*/
/*		2.10	Added RTD capability to enter Callendar-	*/
/*				Van Dusen coeff obtained from calibration	*/
/*				Nuri Cankurt 9-28-2008						*/
/*															*/
/*      2.20    Added deg C versions for all RTD equations	*/
/*              Nuri Cankurt 9-2-2009						*/
/*															*/
/*      2.30    Added Type N, R, S, T thermocouple equations*/
/*              Nuri Cankurt 5-17-2010						*/
/*		3.00	Modified for Visual Studio 2010				*/
/*				added preprocessor definitions for			*/
/*				_CRT_SECURE_NO_WARNINGS						*/
/*				Nuri Cankurt  May 23, 2011					*/
/*		3.10	Modified for Visual Studio 2010				*/
/*				added preprocessor definitions for			*/
/*				_USE_32BIT_TIME_T							*/
/*				Nuri Cankurt  August 30, 2012				*/
/*		3.20	Preprocessor definition for					*/
/*				_USE_32BIT_TIME_T did not work				*/
/*				Had to define it explicitly in the source	*/
/*				tdsEng.c file								*/
/*				Nuri Cankurt  October 15, 2012				*/
/************************************************************/

#define _USE_32BIT_TIME_T	// Nuri Cankurt 10-15-2012

#include <stdio.h> 
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#ifdef _HPUX_SOURCE
#else
	#include <io.h>
	#include <windows.h>
	#include <sys/timeb.h>
#endif
#include <time.h>
#include <math.h>

#define MAX_STR 512

long error;
#define	DllExport extern __declspec( dllexport )

/************************************************************************/
/*																		*/
/* utcToLocalDateTimeStr: This function returns local date (mm/dd/yyyy) */
/*		and local time (hh:mm:ss) as string from UTC seconds given		*/
/*																		*/
/* If UTC seconds is in error, date and time strings are set to NULL	*/
/************************************************************************/
#ifdef _HPUX_SOURCE
	void utcToLocalDateTimeStr (long UTCsec, char *dateStr, char *timeStr)
#else
	DllExport void __cdecl utcToLocalDateTimeStr (long UTCsec, char *dateStr, char *timeStr)
#endif
{
	const	time_t *timer;
	const struct tm *timeptr;

	if (UTCsec <=0)
	{
		strcpy(dateStr,"");
		strcpy(timeStr,"");
		return;
	}

	timer = (const time_t *)&UTCsec;
	timeptr = localtime(timer);

	if (timeptr == NULL)
	{
		strcpy(dateStr,"");
		strcpy(timeStr,"");
		return;
	} 

	sprintf(dateStr,"%02d/%02d/%04d\0",(long) timeptr->tm_mon+1,(long) timeptr->tm_mday,(long) timeptr->tm_year+1900);
	sprintf(timeStr,"%02d:%02d:%02d\0",(long) timeptr->tm_hour,(long) timeptr->tm_min,(long) timeptr->tm_sec);

	return;
}

/************************************************************************/
/*																		*/
/* utcToLocalDateTimeStr2: This function returns local date (mm/dd/yyyy)*/
/*		and local time (hh:mm:ss) as string and millisec from fractional*/
/*		UTC seconds given												*/
/*																		*/
/* If UTC seconds is in error, date and time strings are set to NULL	*/
/************************************************************************/
#ifdef _HPUX_SOURCE
	void utcToLocalDateTimeStr2 (double fractionalUTCsec, char *dateStr, char *timeStr, long *millisec)
#else
	DllExport void __cdecl utcToLocalDateTimeStr2 (double fractionalUTCsec, char *dateStr, char *timeStr, long *millisec)
#endif
{
	long UTCsec;
	double xx;
	
	UTCsec=(long) fractionalUTCsec;

	if (UTCsec <=0)
	{
		strcpy(dateStr,"");
		strcpy(timeStr,"");
		*millisec=0;
		return;
	}

	utcToLocalDateTimeStr (UTCsec, dateStr, timeStr);

	xx=(fractionalUTCsec - (double) UTCsec) * 1000.0;
	*millisec = (long) xx;

	return;
}

/********************************************************************/
/*                                                                  */
/* linearLSCF: This function performs linear least square curve fit.*/
/* and returns slope, intercept and correlation coef				*/
/* Calling parameters:												*/
/*   xdata=vector containing independent variable data				*/
/*   ydata=vector containing dependent variable data				*/
/*   num=number of data points in xdata (and in ydata)				*/
/*																	*/
/* Returned parameters:												*/
/*   slope=slope of Y=aX + B (slope=a)								*/
/*   intercept=intercept of Y=aX + B (intercept=B)					*/
/*   corrCoef=Correlation coefficient R								*/
/*																	*/
/* Function returns zero if no error, -1 if error					*/
/********************************************************************/
#ifdef _HPUX_SOURCE
	long linearLSCF (double *xdata, double *ydata, long num, double *slope, double *intercept, double *corrCoef)
#else
	DllExport long __cdecl linearLSCF (double *xdata, double *ydata, long num, double *slope, double *intercept, double *corrCoef)
#endif
{
	double sumX;
	double sumX2;
	double sumY;
	double sumY2;
	double sumXY;
	double a, b, r;
	long ii;
	double xnum;

	*slope=0.0;
	*intercept=0.0;
	*corrCoef=0.0;

	if (num < 2) return -1;

	if (num == 2)
	{
		if (xdata[0] == xdata[1])return -1;
		a=(ydata[0] - ydata[1])/(xdata[0] - xdata[1]);

		b=ydata[0] - a*xdata[0];
		*slope=a;
		*intercept=b;
		*corrCoef=1.0;
		return 0;
	}

	sumX=0.0;
	sumX2=0.0;
	sumY=0.0;
	sumY2=0.0;
	sumXY=0.0;
	xnum=(double)num;

	for (ii=0;ii<num;ii++)
	{
		sumX=sumX + xdata[ii];
		sumX2=sumX2 + xdata[ii] * xdata[ii];
		sumY=sumY + ydata[ii];
		sumY2=sumY2 + ydata[ii] * ydata[ii];
		sumXY=sumXY + xdata[ii] * ydata[ii];
	}

	if (xnum * sumX2 == sumX * sumX)return -1;

	a=(xnum * sumXY - sumX * sumY)/(xnum * sumX2 - sumX * sumX);
	b=(sumY/xnum) - (a * (sumX/xnum));

	*slope=a;
	*intercept=b;

	r=(xnum * sumY2 - sumY * sumY) * (xnum * sumX2 - sumX * sumX);

	if (r <= 0.0)return -1;

	r=(xnum * sumXY - sumX * sumY)/sqrt(r);

	*corrCoef=r;

	return 0;
}

/* ============= sortDoublePtr ========================== */

/* Modified to run on PC platform. Nuri Cankurt July 2, 2003 */

/* This function sorts a given double array 'fdata' to create
   a pointer array 'index'.  It does not change 'fdata' in any way.

   The sorted values of 'fdata' would be:
   fdata[index[0]], fdata[index[1]], fdata[index[2]] ....

   assending or descending.

   Written by: N. Cankurt 10-30-97
   
   Input variables:
   fdata  = data
   num    = how many values are there in fdata.
   how   = how to sort 'a' for ascending and 'd' for descending

   Returned variables:
   index  = contains the index of fdata which would be in sorted order

   Function returns 0 if there is no error.  -1 if there is an error

   Bug fixed;  Nuri 4-24-2007 for descending flag

*/
#ifdef _HPUX_SOURCE
	long sortDoublePtr (double *fdata, long num, char how, long *index)
#else
	DllExport long __cdecl sortDoublePtr (double *fdata, long num, char how, long *index)
#endif
{
	long ii,jj,noswap;
	double temp;
	char chow;

	if (num < 1)return(-1);
//	if (how != 'A' && how != 'a' && how != 'B' && how != 'b')return -1;	  // this has a bug Nuri  4-24-2007
	if (how != 'A' && how != 'a' && how != 'D' && how != 'd')return -1;	  // fixed the bug  4-24-2007

	for (ii=0;ii<num;ii++)index[ii]=ii;
	if (num == 1)return(0);

	chow=how;
	if (chow == 'A')chow='a';
	if (chow == 'D')chow='d';

	switch (chow)
	{
		case 'a':   /* assending */
        for (ii=0;ii<num-1;ii++)
		{
			noswap=1;
			for (jj=num-1;jj>ii;jj--)
			{
				if (fdata[index[jj]]<fdata[index[jj-1]])
				{
					temp=(double)index[jj];
					index[jj]=index[jj-1];
					index[jj-1]=(int)temp;
					noswap=0;
				}
			}
			if (noswap) return(0);
        }
		return (0);

		case 'd':   /* descending */
		for (ii=0;ii<num-1;ii++)
		{
			noswap=1;
			for (jj=num-1;jj>ii;jj--)
			{
				if (fdata[index[jj]]>fdata[index[jj-1]])
				{
					temp=(double)index[jj];
					index[jj]=index[jj-1];
					index[jj-1]=(int)temp;
					noswap=0;
				}
			}
			if (noswap) return(0);
		}
		return (0);
	}
	return (-1);
}

/********************************************************************************************/
/*																							*/
/* ctimeToDbl: This function converts'timeString' to double value and returns it			*/
/*																			   				*/
/*  timeString must be a string containing DDDHHMMSSsss	(standard datum output)				*/
/*		where:																				*/
/*			DDD=julian day																	*/
/*			HH =hour																		*/
/*			MM =minute																		*/
/*			SS =second																		*/
/*			sss=millisecond																	*/
/*																							*/
/********************************************************************************************/
#ifdef _HPUX_SOURCE
	double ctimeToDbl (char *timeString)
#else
	DllExport double __cdecl ctimeToDbl (char *timeString)
#endif
{
	double val;
	sscanf(timeString,"%lg",&val);
	return val;
}

/********************************************************************************************/
/*																							*/
/* timeToDDDHHMMSS_MS: This function converts'timeVal' to double value as DDDHHMMSS.sss		*/
/*																			   				*/
/*  timeVal must be a double containing DDDHHMMSSsss	(standard datum output)				*/
/*		where:																				*/
/*			DDD=julian day																	*/
/*			HH =hour																		*/
/*			MM =minute																		*/
/*			SS =second																		*/
/*			sss=millisecond																	*/
/*																							*/
/********************************************************************************************/
#ifdef _HPUX_SOURCE
	double timeToDDDHHMMSS_MS (double timeVal)
#else
	DllExport double __cdecl timeToDDDHHMMSS_MS (double timeVal)
#endif
{
	return timeVal/1000.0;
}

/********************************************************************************************/
/*																							*/
/* timeToDDDHHMMSS: This function converts'timeVal' to double value as DDDHHMMSS			*/
/*																			   				*/
/*  timeVal must be a double containing DDDHHMMSSsss	(standard datum output)				*/
/*		where:																				*/
/*			DDD=julian day																	*/
/*			HH =hour																		*/
/*			MM =minute																		*/
/*			SS =second																		*/
/*			sss=millisecond																	*/
/*																							*/
/********************************************************************************************/
#ifdef _HPUX_SOURCE
	long timeToDDDHHMMSS (double timeVal)
#else
	DllExport long __cdecl timeToDDDHHMMSS (double timeVal)
#endif
{
	double val;
	long val2;

	val=timeVal/1000.0;
	val2=(long) val;

	return val2;
}

/********************************************************************************************/
/*																							*/
/* timeToDay: This function extracts day from 'timeVal'										*/
/*																			   				*/
/*  timeVal must be a double containing DDDHHMMSSsss	(standard datum output)				*/
/*		where:																				*/
/*			DDD=julian day																	*/
/*			HH =hour																		*/
/*			MM =minute																		*/
/*			SS =second																		*/
/*			sss=millisecond																	*/
/*																							*/
/********************************************************************************************/
#ifdef _HPUX_SOURCE
	long timeToDay (double timeVal)
#else
	DllExport long __cdecl timeToDay (double timeVal)
#endif
{
	double val;

	val=timeVal/1000000000.0;

	return (long) val;
}

/********************************************************************************************/
/*																							*/
/* timeToHour: This function extracts hour from 'timeVal'									*/
/*																			   				*/
/*  timeVal must be a double containing DDDHHMMSSsss	(standard datum output)				*/
/*		where:																				*/
/*			DDD=julian day																	*/
/*			HH =hour																		*/
/*			MM =minute																		*/
/*			SS =second																		*/
/*			sss=millisecond																	*/
/*																							*/
/********************************************************************************************/
#ifdef _HPUX_SOURCE
	long timeToHour (double timeVal)
#else
	DllExport long __cdecl timeToHour (double timeVal)
#endif
{
	double val;
	long val2;
	long val3;

	val=timeVal/1000000000.0;
	val2=(long)val;

	val=timeVal/10000000.0;
	val3=(long)val;

	val2=val3 - val2*100;

	return val2;
}

/********************************************************************************************/
/*																							*/
/* timeToMinute: This function extracts minute from 'timeVal'								*/
/*																			   				*/
/*  timeVal must be a double containing DDDHHMMSSsss	(standard datum output)				*/
/*		where:																				*/
/*			DDD=julian day																	*/
/*			HH =hour																		*/
/*			MM =minute																		*/
/*			SS =second																		*/
/*			sss=millisecond																	*/
/*																							*/
/********************************************************************************************/
#ifdef _HPUX_SOURCE
	long timeToMinute (double timeVal)
#else
	DllExport long __cdecl timeToMinute (double timeVal)
#endif
{
	double val;
	long val2;
	long val3;

	val=timeVal/10000000.0;
	val2=(long)val;

	val=timeVal/100000.0;
	val3=(long) val;

	val2=val3 - val2*100;

	return val2;
}

/********************************************************************************************/
/*																							*/
/* timeToSecond: This function extracts second from 'timeVal'								*/
/*																			   				*/
/*  timeVal must be a double containing DDDHHMMSSsss	(standard datum output)				*/
/*		where:																				*/
/*			DDD=julian day																	*/
/*			HH =hour																		*/
/*			MM =minute																		*/
/*			SS =second																		*/
/*			sss=millisecond																	*/
/*																							*/
/********************************************************************************************/
#ifdef _HPUX_SOURCE
	long timeToSecond (double timeVal)
#else
	DllExport long __cdecl timeToSecond (double timeVal)
#endif
{
	double val;
	long val2;
	long val3;

	val=timeVal/100000.0;
	val2=(long)val;

	val=timeVal/1000.0;
	val3=(long) val;

	val2=val3 - val2*100;

	return val2;
}

/********************************************************************************************/
/*																							*/
/* timeToMillisec: This function extracts millisecond from 'timeVal'						*/
/*																			   				*/
/*  timeVal must be a double containing DDDHHMMSSsss	(standard datum output)				*/
/*		where:																				*/
/*			DDD=julian day																	*/
/*			HH =hour																		*/
/*			MM =minute																		*/
/*			SS =second																		*/
/*			sss=millisecond																	*/
/*																							*/
/********************************************************************************************/
#ifdef _HPUX_SOURCE
	long timeToMillisec (double timeVal)
#else
	DllExport long __cdecl timeToMillisec (double timeVal)
#endif
{
	double val;
	long val2;

	val2=timeToDDDHHMMSS(timeVal);
	val=(double)val2;

	val=timeVal - val*1000;
	val2=(long) val;

	return val2;
}

/********************************************************************************************/
/*																							*/
/* timeToHHMMSS: This function extracts HHMMSS from 'timeVal'								*/
/*																			   				*/
/*  timeVal must be a double containing DDDHHMMSSsss	(standard datum output)				*/
/*		where:																				*/
/*			DDD=julian day																	*/
/*			HH =hour																		*/
/*			MM =minute																		*/
/*			SS =second																		*/
/*			sss=millisecond																	*/
/*																							*/
/********************************************************************************************/
#ifdef _HPUX_SOURCE
	long timeToHHMMSS (double timeVal)
#else
	DllExport long __cdecl timeToHHMMSS (double timeVal)
#endif
{
	double val;
	long val2;
	long val3;

	val=timeVal/1000.0;
	val2=(long) val;

	val=timeVal/1000000000.0;
	val3=(long) val;

	val2=val2 - val3*1000000;

	return val2;
}

/********************************************************************************************/
/*																							*/
/* timeToHHMMSS_MS: This function extracts HHMMSS.sss from 'timeVal'						*/
/*																			   				*/
/*  timeVal must be a double containing DDDHHMMSSsss	(standard datum output)				*/
/*		where:																				*/
/*			DDD=julian day																	*/
/*			HH =hour																		*/
/*			MM =minute																		*/
/*			SS =second																		*/
/*			sss=millisecond																	*/
/*																							*/
/********************************************************************************************/
#ifdef _HPUX_SOURCE
	double timeToHHMMSS_MS (double timeVal)
#else
	DllExport double __cdecl timeToHHMMSS_MS (double timeVal)
#endif
{
	double val;
	double val2;

	val=timeToDDDHHMMSS_MS(timeVal);
	val2=(double) timeToDay(timeVal);

	val2=val - val2*1000000.0;

	return val2;
}

/********************************************************************************************/
/*																							*/
/* timeToDecDay: This function converts 'timeVal' to decimal day and returns it				*/
/*																			   				*/
/*  time must be a double containing DDDHHMMSSsss	(standard datum output)					*/
/*		where:																				*/
/*			DDD=julian day																	*/
/*			HH =hour																		*/
/*			MM =minute																		*/
/*			SS =second																		*/
/*			sss=millisecond																	*/
/*																							*/
/********************************************************************************************/
#ifdef _HPUX_SOURCE
	double timeToDecDay (double timeVal)
#else
	DllExport double __cdecl timeToDecDay (double timeVal)
#endif
{
	double val;
	long day,hour,minute,second,millisec;
	double result;

	day=timeToDay(timeVal);
	hour=timeToHour(timeVal);
	minute=timeToMinute(timeVal);
	second=timeToSecond(timeVal);
	millisec=timeToMillisec(timeVal);

	result=(double)day;

	val=(double)hour;
	val=val/24.0;
	result=result+val;

	val=(double)minute;
	val=val/(24.0*60.0);
	result=result+val;

	val=(double)second;
	val=val/(24.0*60.0*60.0);
	result=result+val;

	val=(double)millisec;
	val=val/(24.0*60.0*60.0*1000.0);
	result=result+val;

	return result;
}

/********************************************************************************************/
/*																							*/
/* timeToYearSecond: This function converts 'timeVal' to seconds since the year began.		*/
/*  timeVal must be a double containing DDDHHMMSSsss	(standard datum output)				*/
/*		where:																				*/
/*			DDD=julian day																	*/
/*			HH =hour																		*/
/*			MM =minute																		*/
/*			SS =second																		*/
/*			sss=millisecond																	*/
/*																							*/
/********************************************************************************************/
#ifdef _HPUX_SOURCE
	double timeToYearSecond (double timeVal)
#else
	DllExport double __cdecl timeToYearSecond (double timeVal)
#endif
{
	double val;

	val=timeToDecDay(timeVal);
	if(val < 1.0)return 0.0;

	val=(val-1.0) * 24.0 * 60.0 * 60.0;

	return val;
}

/********************************************************************************************/
/*																							*/
/* timeToStdStr: This function converts 'timeVal' to 'strout' as 'DDD-HH:MM:SS.sss'			*/
/*  timeVal must be a double containing DDDHHMMSSsss	(standard datum output)				*/
/*		where:																				*/
/*			DDD=julian day																	*/
/*			HH =hour																		*/
/*			MM =minute																		*/
/*			SS =second																		*/
/*			sss=millisecond																	*/
/*																							*/
/*  Function always returns 0																*/
/********************************************************************************************/
#ifdef _HPUX_SOURCE
	long timeToStdStr (double timeVal, char *strout)
#else
	DllExport long __cdecl timeToStdStr (double timeVal, char *strout)
#endif
{
	sprintf(strout,"%03d-%02d:%02d:%02d.%03d\0",timeToDay(timeVal),
		timeToHour(timeVal),timeToMinute(timeVal),timeToSecond(timeVal),
		timeToMillisec(timeVal));

	return 0;
}

/********************************************************************************************/
/*																							*/
/* timeToCtime: This function converts 'timeVal' to 'strout' as 'DDDHHMMSSsss'			*/
/*  timeVal must be a double containing DDDHHMMSSsss	(standard datum output)				*/
/*		where:																				*/
/*			DDD=julian day																	*/
/*			HH =hour																		*/
/*			MM =minute																		*/
/*			SS =second																		*/
/*			sss=millisecond																	*/
/*																							*/
/*  Function always returns 0																*/
/********************************************************************************************/
#ifdef _HPUX_SOURCE
	long timeToCtime (double timeVal, char *strout)
#else
	DllExport long __cdecl timeToCtime (double timeVal, char *strout)
#endif
{
	sprintf(strout,"%03d%02d%02d%02d%03d\0",timeToDay(timeVal),
		timeToHour(timeVal),timeToMinute(timeVal),timeToSecond(timeVal),
		timeToMillisec(timeVal));

	return 0;
}

/********************************************************************************************/
/*																							*/
/* tck_V32toF: This function calculates temperature in deg F from volts when cold junction	*/
/*	reference temperature is 32 F (standard temprature) for type K TC. The function returns */
/*	calculated temperature.																	*/
/*																							*/
/*							This function is valid for TYPE K thermocouple					*/
/* 
   In this function, NIST equations refer to ITS-90 Table.

   This function replaces previous functions 'typek_degfC.c and typek_degf.f".

   This function calculates temperature in deg F referenced to 32 deg F from type K thermocouple. 
   Th input is variable 'volts' in volts.

      The quations were obtained from N.I.S.T web site on July 16, 2002.  The NIST equations are valid
      for temperatures from:
       -200 deg C to 1372 deg C (-328 deg F to 2501.6 deg F) or -5.891 millivolts to 54.886 millivolts.

      Accuracy of these NIST equations are as follows:
      -200 to    0 deg C (-328 to   32 deg F)   Accuracy= -0.02 to 0.04 deg C (-0.036 to  0.072 deg F)  
	 0 to  500 deg C (  32 to  932 deg F)   Accuracy= -0.05 to 0.04 deg C (-0.09  to  0.072 deg F)
       500 to 1372 deg C ( 932 to 2501.6 deg F) Accuracy= -0.05 to 0.06 deg C (-0.09  to  0.108 deg F)

   If 'volts' is outside of temperature range, function returns 0 deg F.

   'volts' is in volts
   'typek_degf2' is in deg F

   'typek_degf2' Written by: Nuri Cankurt July 23, 2002  SRC-330
   dll version 'tck_V32toF': written by Nuri Cankurt Jan 12, 2004
*/
/********************************************************************************************/
#ifdef _HPUX_SOURCE
	double tck_V32toF (double volts)
#else
	DllExport double __cdecl tck_V32toF (double volts)
#endif
{
   double tf;  /* deg F */
   double tc;  /* deg C */
   double milli_volt;
   double emf;

   tf=0.0;
   tc=0.0;	// Nuri Cankurt 5-17-2010
   if (volts < -5.891e-3) return tf;
   if (volts > 54.886e-3) return tf;
   if (volts == 0.0) return tf;

   milli_volt= volts * 1000.0;
   emf=milli_volt;

   if (emf >= -5.891 &&  emf <= 0.0)  /* -200 deg C to 0 deg C (-328 to   32 deg F) */
   {
      tc=     0.0000000E+00 
             +2.5173462E+01*emf  
             -1.1662878E+00*emf*emf
             -1.0833638E+00*emf*emf*emf
             -8.9773540E-01*emf*emf*emf*emf
             -3.7342377E-01*emf*emf*emf*emf*emf
             -8.6632643E-02*emf*emf*emf*emf*emf*emf
             -1.0450598E-02*emf*emf*emf*emf*emf*emf*emf
             -5.1920577E-04*emf*emf*emf*emf*emf*emf*emf*emf; 
   }
   else if (emf > 0.0 && emf <= 20.644)  /* 0 deg C to 500 deg C (  32 to 932 deg F) */
   {
      tc=   0.000000E+00 
	   +2.508355E+01*emf 
           +7.860106E-02*emf*emf
           -2.503131E-01*emf*emf*emf 
           +8.315270E-02*emf*emf*emf*emf
           -1.228034E-02*emf*emf*emf*emf*emf
           +9.804036E-04*emf*emf*emf*emf*emf*emf 
           -4.413030E-05*emf*emf*emf*emf*emf*emf*emf 
           +1.057734E-06*emf*emf*emf*emf*emf*emf*emf*emf 
           -1.052755E-08*emf*emf*emf*emf*emf*emf*emf*emf*emf; 
   }
   else if (emf > 20.644 && emf <= 54.886)  /* 500 deg C to 1372 deg C ( 932 to 2501.6 degF) */
   {
      tc= -1.318058E+02
          +4.830222E+01*emf
          -1.646031E+00*emf*emf
          +5.464731E-02*emf*emf*emf
          -9.650715E-04*emf*emf*emf*emf
          +8.802193E-06*emf*emf*emf*emf*emf
          -3.110810E-08*emf*emf*emf*emf*emf*emf;
   }
   tf=tc*1.8 + 32.0;  /* deg F */

   return tf;
}

/********************************************************************************************/
/*																							*/
/* tck_FtoV32: This function calculates the voltage (in volts) for a Type K thermocouple	*/
/*	when the cold juction reference temperature is at 32 F (standard temperature).			*/
/*	The function returns the calculated voltage.											*/
/*																							*/
/*							This function is valid for TYPE K thermocouple					*/
/* 
/* 
   In this function, NIST equations refer to ITS-90 Table.

   This function replaces previous functions 'typek_voltC.c and typek_volt.f".
	 
   This function calculates volts referenced to 32 deg F from type K thermocouple at 'temp' deg F.
   the quation was obtained from N.I.S.T web site on July 16, 2002.

   This function is valid for temperatures from -270 deg C  to 1372 deg C (-454 deg F to 2501.6 deg F)
   If 'temp' is outside of this range, function returns 0 volts.

   'temp' is in deg F
   'typek_volts2' is in volts

   'typek_volts2' Written by: Nuri Cankurt July 23, 2002  SRC-330
    dll version 'tck_FtoV32': written by Nuri Cankurt Jan 12, 2004
*/
#ifdef _HPUX_SOURCE
	double tck_FtoV32 (double temp)
#else
	DllExport double __cdecl tck_FtoV32 (double temp)
#endif
{
   double volt;
   double tf;  /* deg F */
   double tc;  /* deg C */
   double milli_volt;

   volt=0.0;
   milli_volt=0.0;	// Nuri Cankurt 5-17-2010

   if (temp < -454.0 || temp > 2501.6) return volt;
   if (temp == 0.0) return volt;
   if (temp == 32.0)return volt;

   tf= temp;
   tc=(tf-32.0)/1.8;

   if (tc >= -270.0 && tc <= 0.0)  /* For temps from -270 deg C to 0 deg C (-454 deg F to 32 deg F) */
   {
      milli_volt     =   0.000000000000E+00
                        +0.394501280250E-01*tc
                        +0.236223735980E-04*tc*tc
                        -0.328589067840E-06*tc*tc*tc
                        -0.499048287770E-08*tc*tc*tc*tc
                        -0.675090591730E-10*tc*tc*tc*tc*tc
                        -0.574103274280E-12*tc*tc*tc*tc*tc*tc
                        -0.310888728940E-14*tc*tc*tc*tc*tc*tc*tc
                        -0.104516093650E-16*tc*tc*tc*tc*tc*tc*tc*tc
                        -0.198892668780E-19*tc*tc*tc*tc*tc*tc*tc*tc*tc
                        -0.163226974860E-22*tc*tc*tc*tc*tc*tc*tc*tc*tc*tc;
    }
    else if (tc >0.0 && tc <= 1372.0) /* For temps from 0 degC to 1372 deg C ( 32 deg F to 2501.6 deg F) */
    {
       milli_volt     = -0.176004136860E-01
                        +0.389212049750E-01*tc
                        +0.185587700320E-04*tc*tc
                        -0.994575928740E-07*tc*tc*tc
                        +0.318409457190E-09*tc*tc*tc*tc
                        -0.560728448890E-12*tc*tc*tc*tc*tc
                        +0.560750590590E-15*tc*tc*tc*tc*tc*tc
                        -0.320207200030E-18*tc*tc*tc*tc*tc*tc*tc
                        +0.971511471520E-22*tc*tc*tc*tc*tc*tc*tc*tc
                        -0.121047212750E-25*tc*tc*tc*tc*tc*tc*tc*tc*tc
                        +0.118597600000E+00*exp(-0.118343200000E-03*(tc-126.9686)*(tc-126.9686));
   }
   volt= (double) milli_volt / 1000.0;
   return volt;
}

/********************************************************************************************/
/*																							*/
/* tck_VtoF: This function calculates temperature in deg F from volts at cold junction		*/
/*  reference temperature 'refTemp' (F).  The function returns the calculated temperature	*/
/*																							*/
/*							This function is valid for TYPE K thermocouple					*/
/*																							*/
/********************************************************************************************/
#ifdef _HPUX_SOURCE
	double tck_VtoF (double volts, double refTemp)
#else
	DllExport double __cdecl tck_VtoF (double volts, double refTemp)
#endif
{
	double refVolts;

	if (refTemp == 0.0) return 0.0;
	if (volts == 0.0) return 0.0;
	if (refTemp == 32.0)return tck_V32toF (volts);

	refVolts=tck_FtoV32(refTemp);
	if (refVolts == 0.0)return 0.0;

	return tck_V32toF((refVolts+volts));
}

/********************************************************************************************/
/*																							*/
/* tcb_V32toF: This function calculates temperature in deg F from volts when cold junction	*/
/*	reference temperature is 32 F (standard temprature) for type B TC. The function returns */
/*	calculated temperature.																	*/
/*																							*/
/*							This function is valid for TYPE B thermocouple					*/
/* 
   In this function, NIST equations refer to ITS-90 Table.

   This function replaces previous functions 'typeb_degfC.c and typeb_degf.f".

   This function calculates temperature in deg F referenced to 32 deg F from type B thermocouple. 
   Th input is variable 'volts' in volts.

   This function calculates temperatures using various segments:

   1) For temperatures from 250 deg C  to 1820 deg C (482 deg F to 3308 deg F), the quation was obtained 
      from N.I.S.T web site on July 16, 2002.  This NIST equation is valid for temperatures from
      250 deg C to 1820 deg C (482 deg F to 3308 deg F) or 0.291 millivolts to 13.820 millivolts.

      Accuracy of these NIST equations are as follows:
      250 to  700 deg C ( 482 to 1292 deg F)  Accuracy= -0.02 to 0.03 deg C (-0.036 to  0.54 deg F)
      700 to 1820 deg C (1292 to 3308 deg F)  Accuracy= -0.01 to 0.02 deg C (-0.018 to 0.036 deg F)

   2) For 22.2 deg C to 250 deg C (72 to 482 deg F), GTTL curve fitted NIST voltage data
      
      Accuracy of GTTL curves are as follows:
       72 to 120 deg F     Accuracy= -0.52 to 1.00 deg F
      120 to 200 deg F     Accuracy= -0.22 to 0.01 deg F
      200 to 482 deg F     Accuracy= -0.19 to 0.10 deg F

   3) Below 22.2 deg C (72 deg F), there is no equation.  The function returns extrapolation of 72 to 120 deg F equation
      which is totally in error.


   If 'volts' is outside of temperature range, function returns 0 deg F.

   'volts' is in volts
   'typeb_degf2' is in deg F

   'typeb_degf2' Written by: Nuri Cankurt July 17, 2002  SRC-327
   dll version 'tcb_V32toF': written by Nuri Cankurt Jan 12, 2004
*/
/********************************************************************************************/
#ifdef _HPUX_SOURCE
	double tcb_V32toF (double volts)
#else
	DllExport double __cdecl tcb_V32toF (double volts)
#endif
{
   double tf;  /* deg F */
   double tc;  /* deg C */
   double milli_volt;
   double emf;

   tf=0.0;
   tc=0.0;	// Nuri Cankurt 5-17-2010

   if (volts == 0.0)return tf;
   if (volts < -0.002585e-3) return tf;
   if (volts > 13.820e-3) return tf;

   milli_volt= volts * 1000.0;
   emf=milli_volt;

   if (emf <= 0.00191344 )              /* GTTL curve fit valid for 72 to 120 F.  However due to double valued */
   {                                    /* nature of type B TC at low temps, function will return values even  */
				        /* below 72 F with large errors */
      tf = 1.081395E+02                 /*   72 to 120 deg F     Accuracy= -0.52 to 1.00 deg F                 */
          + 6.970012E+03*emf
          - 2.065585E+06*emf*emf
          + 6.670762E+08*emf*emf*emf
          + 8.525493E+11*emf*emf*emf*emf
          - 1.135224E+14*emf*emf*emf*emf*emf
          - 1.477129E+17*emf*emf*emf*emf*emf*emf;
      tc=(tf-32.0)/1.8;
   }

   else if (emf <= 0.0274527 )         /* GTTL curve fit valid for 120 to 200 F.   */
   {                                   /* 120 to 200 deg F     Accuracy= -0.22 to 0.01 deg F  */
      tf = 1.083872E+02
          + 6.703680E+03*emf
          - 3.659377E+05*emf*emf
          + 2.219400E+07*emf*emf*emf
          - 8.965692E+08*emf*emf*emf*emf
          + 2.016542E+10*emf*emf*emf*emf*emf
          - 1.892776E+11*emf*emf*emf*emf*emf*emf;
      tc=(tf-32.0)/1.8;
   }

   else if (emf < 0.291  )         /* GTTL curve fit valid for 200 to 482 F.   */
   {                               /* 200 to 482 deg F     Accuracy= -0.19 to 0.10 deg F  */ 
      tf =  1.306166E+02
          + 3.043443E+03*emf
          - 2.203554E+04*emf*emf
          + 1.432234E+05*emf*emf*emf
          - 5.762741E+05*emf*emf*emf*emf
          + 1.254962E+06*emf*emf*emf*emf*emf
          - 1.126887E+06*emf*emf*emf*emf*emf*emf;
      tc=(tf-32.0)/1.8;
   }

   else if (emf >= 0.291 &&  emf <= 2.431)  /* NIST eqn. Temp is between 250 deg C to 700 deg C (482 deg F to 1292 deg F) */
   {                                        /* Accuracy= -0.02 to 0.03 deg C (-0.036 to  0.54 deg F)            */
      tc=  9.8423321e+1
	  +6.9971500e+2*emf
	  -8.4765304e+2*emf*emf
	  +1.0052644e+3*emf*emf*emf
	  -8.3345952e+2*emf*emf*emf*emf
	  +4.5508542e+2*emf*emf*emf*emf*emf
	  -1.5523037e+2*emf*emf*emf*emf*emf*emf
	  +2.9886750e+1*emf*emf*emf*emf*emf*emf*emf
	  -2.4742860   *emf*emf*emf*emf*emf*emf*emf*emf;
   }
   else if (emf > 2.431 && emf <= 13.820 )  /* NIST eqn. 700 deg C to 1820 deg C (1292 to 3308 deg F) */
   {                                        /* Accuracy= -0.01 to 0.02 deg C (-0.018 to 0.036 deg F)  */
      tc=  2.1315071e+2
          +2.8510504e+2*emf
	  -5.2742887e+1*emf*emf
	  +9.9160804e+0*emf*emf*emf
	  -1.2965303e+0*emf*emf*emf*emf
	  +1.1195870e-1*emf*emf*emf*emf*emf
	  -6.0625199e-3*emf*emf*emf*emf*emf*emf
	  +1.8661696e-4*emf*emf*emf*emf*emf*emf*emf
	  -2.4878585e-6*emf*emf*emf*emf*emf*emf*emf*emf;
   }
   tf=tc*1.8 + 32.0;  /* deg F */

   return tf;
}

/********************************************************************************************/
/*																							*/
/* tcb_FtoV32: This function calculates the voltage (in volts) for a Type B thermocouple	*/
/*	when the cold juction reference temperature is at 32 F (standard temperature).			*/
/*	The function returns the calculated voltage.											*/
/*																							*/
/*							This function is valid for TYPE B thermocouple					*/
/* 
/* 
   In this function, NIST equations refer to ITS-90 Table.

   This function replaces previous functions 'typeb_voltC.c and typeb_volt.f".
	 
   This function calculates volts referenced to 32 deg F from type B thermocouple at 'temp' deg F.
   the quation was obtained from N.I.S.T web site on July 16, 2002.

   This function is valid for temperatures from 0 deg C  to 1820 deg C (32 deg F to 3308 deg F)
   If 'temp' is outside of this range, function returns 0 volts.

   'temp' is in deg F
   'typeb_volts2' is in volts

   'typeb_volts2' Written by: Nuri Cankurt July 17, 2002  SRC-327
    dll version 'tcb_FtoV32': written by Nuri Cankurt Jan 12, 2004
*/
#ifdef _HPUX_SOURCE
	double tcb_FtoV32 (double temp)
#else
	DllExport double __cdecl tcb_FtoV32 (double temp)
#endif
{
   double volt;
   double tf;  /* deg F */
   double tc;  /* deg C */
   double milli_volt;

   volt=0.0;
   milli_volt=0.0;	// Nuri Cankurt 5-17-2010

   if (temp < 32.0 || temp > 3308.0) return volt;
   if (temp == 32.0)return volt;

   tf= temp;
   tc=(tf-32.0)/1.8;

   if (tc >= 0.0 && tc <= 630.615)  /* Temp is between 0 deg C to 630.615 deg C (32 deg F and 1167.107 deg F) */
   {
      milli_volt= 0.0
		       -0.246508183460e-3*tc 
		       +0.590404211710e-5*tc*tc
		       -0.132579316360e-8*tc*tc*tc
		       +0.156682919010e-11*tc*tc*tc*tc
		       -0.169445292400e-14*tc*tc*tc*tc*tc
		       +0.629903470940e-18*tc*tc*tc*tc*tc*tc;
   }
   else if (tc >630.615 && tc <= 1820.0)  /* Temp is between 630.615 deg C to 1820 deg C (1167.107 deg F and 3308 deg F */  
   {
      milli_volt= -0.389381686210e+1
            +0.285717474700e-1*tc
            -0.848851047850e-4*tc*tc
			+0.157852801640e-6*tc*tc*tc
			-0.168353448640e-9*tc*tc*tc*tc
			+0.111097940130e-12*tc*tc*tc*tc*tc
			-0.445154310330e-16*tc*tc*tc*tc*tc*tc
			+0.989756408210e-20*tc*tc*tc*tc*tc*tc*tc
			-0.937913302890e-24*tc*tc*tc*tc*tc*tc*tc*tc;
   }
   volt=  milli_volt / 1000.0;
   return volt;
}

/********************************************************************************************/
/*																							*/
/* tcb_VtoF: This function calculates temperature in deg F from volts at cold junction		*/
/*  reference temperature 'refTemp' (F).  The function returns the calculated temperature	*/
/*																							*/
/*							This function is valid for TYPE B thermocouple					*/
/*																							*/
/********************************************************************************************/
#ifdef _HPUX_SOURCE
	double tcb_VtoF (double volts, double refTemp)
#else
	DllExport double __cdecl tcb_VtoF (double volts, double refTemp)
#endif
{
	double refVolts;

	if (refTemp == 0.0)return 0.0;
	if (volts == 0.0)return 0.0;
	if (refTemp == 32.0)return tcb_V32toF (volts);

	refVolts=tcb_FtoV32(refTemp);
	if (refVolts == 0.0)return 0.0;

	return tcb_V32toF((refVolts+volts));
}

/********************************************************************************************/
/*																							*/
/* tce_V32toF: This function calculates temperature in deg F from volts when cold junction	*/
/*	reference temperature is 32 F (standard temprature) for type E TC. The function returns */
/*	calculated temperature.																	*/
/*																							*/
/*							This function is valid for TYPE E thermocouple					*/
/* 
   In this function, NIST equations refer to ITS-90 Table.

   This function calculates temperature in deg F referenced to 32 deg F from type E thermocouple. 
   Th input is variable 'volts' in volts.

      The quations were obtained from N.I.S.T web site on Jan 27, 2004.  The NIST equations are valid
      for temperatures from:
       -200 deg C to 1000 deg C (-328 deg F to 1832 deg F) or -8.825 millivolts to 76.373 millivolts.

      Accuracy of these NIST equations are as follows:
      -200 to    0 deg C (-328 to   32 deg F)   Accuracy= -0.01 to 0.03 deg C (-0.018 to  0.054 deg F)  
	     0 to 1000 deg C (  32 to 1832 deg F)   Accuracy= -0.02 to 0.02 deg C (-0.036  to 0.036 deg F)

   If 'volts' is outside of temperature range, function returns 0 deg F.

   'volts' is in volts

   dll version 'tce_V32toF': written by Nuri Cankurt March 12, 2004
*/
/********************************************************************************************/
#ifdef _HPUX_SOURCE
	double tce_V32toF (double volts)
#else
	DllExport double __cdecl tce_V32toF (double volts)
#endif
{
   double tf;  /* deg F */
   double tc;  /* deg C */
   double milli_volt;
   double emf;

   tf=0.0;
   tc=0.0;	// Nuri Cankurt 5-17-2010
   if (volts < -8.825e-3) return tf;
   if (volts > 76.373e-3) return tf;
   if (volts == 0.0) return tf;

   milli_volt= volts * 1000.0;
   emf=milli_volt;

   if (emf >= -8.825 &&  emf <= 0.0)  /* -200 deg C to 0 deg C (-328 to   32 deg F) */
   {
      tc=    0.0000000E+00 
            +1.6977288E+01*emf  
            -4.3514970E-01*emf*emf
            -1.5859697E-01*emf*emf*emf
            -9.2502871E-02*emf*emf*emf*emf
            -2.6084314E-02*emf*emf*emf*emf*emf
            -4.1360199E-03*emf*emf*emf*emf*emf*emf
            -3.4034030E-04*emf*emf*emf*emf*emf*emf*emf
            -1.1564890E-05*emf*emf*emf*emf*emf*emf*emf*emf; 
   }
   else if (emf > 0.0 && emf <=  76.373)  /* 0 deg C to 1000 deg C (  32 to 1832 deg F) */
   {
      tc=    0.000000E+00 
			+1.7057035E+01*emf 
            -2.3301759E-01*emf*emf
            +6.5435585E-03*emf*emf*emf 
            -7.3562749E-05*emf*emf*emf*emf
            -1.7896001E-06*emf*emf*emf*emf*emf
            +8.4036165E-08*emf*emf*emf*emf*emf*emf 
            -1.3735879E-09*emf*emf*emf*emf*emf*emf*emf 
            +1.0629823E-11*emf*emf*emf*emf*emf*emf*emf*emf 
            -3.2447087E-14*emf*emf*emf*emf*emf*emf*emf*emf*emf; 
   }

   tf=tc*1.8 + 32.0;  /* deg F */

   return tf;
}

/********************************************************************************************/
/*																							*/
/* tce_FtoV32: This function calculates the voltage (in volts) for a Type E thermocouple	*/
/*	when the cold juction reference temperature is at 32 F (standard temperature).			*/
/*	The function returns the calculated voltage.											*/
/*																							*/
/*							This function is valid for TYPE E thermocouple					*/
/* 
/* 
   In this function, NIST equations refer to ITS-90 Table.
 
   This function calculates volts referenced to 32 deg F from type K thermocouple at 'temp' deg F.
   the quation was obtained from N.I.S.T web site on Jan 27, 2004.

   This function is valid for temperatures from -270 deg C  to 1000 deg C (-454 deg F to 1832 deg F)
   If 'temp' is outside of this range, function returns 0 volts.

   'temp' is in deg F

    dll version 'tce_FtoV32': written by Nuri Cankurt March 12, 2004
*/
#ifdef _HPUX_SOURCE
	double tce_FtoV32 (double temp)
#else
	DllExport double __cdecl tce_FtoV32 (double temp)
#endif
{
   double volt;
   double tf;  /* deg F */
   double tc;  /* deg C */
   double milli_volt;

   volt=0.0;
   milli_volt=0.0;	// Nuri Cankurt 5-17-2010

   if (temp < -454.0 || temp > 1832.0) return volt;
   if (temp == 0.0) return volt;
   if (temp == 32.0)return volt;

   tf= temp;
   tc=(tf-32.0)/1.8;

   if (tc >= -270.0 && tc <= 0.0)  /* For temps from -270 deg C to 0 deg C (-454 deg F to 32 deg F) */
   {
      milli_volt     =   0.000000000000E+00
                        +0.586655087080E-01*tc
                        +0.454109771240E-04*tc*tc
                        -0.779980486860E-06*tc*tc*tc
                        -0.258001608430E-07*tc*tc*tc*tc
                        -0.594525830570E-09*tc*tc*tc*tc*tc
                        -0.932140586670E-11*tc*tc*tc*tc*tc*tc
                        -0.102876055340E-12*tc*tc*tc*tc*tc*tc*tc
                        -0.803701236210E-15*tc*tc*tc*tc*tc*tc*tc*tc
                        -0.439794973910E-17*tc*tc*tc*tc*tc*tc*tc*tc*tc
                        -0.164147763550E-19*tc*tc*tc*tc*tc*tc*tc*tc*tc*tc
						-0.396736195160E-22*tc*tc*tc*tc*tc*tc*tc*tc*tc*tc*tc
						-0.558273287210E-25*tc*tc*tc*tc*tc*tc*tc*tc*tc*tc*tc*tc
						-0.346578420130E-28*tc*tc*tc*tc*tc*tc*tc*tc*tc*tc*tc*tc*tc;
    }
    else if (tc >0.0 && tc <= 1000.0) /* For temps from 0 degC to 1000 deg C ( 32 deg F to 1832.0 deg F) */
    {
       milli_volt     =  0.000000000000E+00
                        +0.586655087100E-01*tc
                        +0.450322755820E-04*tc*tc
                        +0.289084072120E-07*tc*tc*tc
                        -0.330568966520E-09*tc*tc*tc*tc
                        +0.650244032700E-12*tc*tc*tc*tc*tc
                        -0.191974955040E-15*tc*tc*tc*tc*tc*tc
                        -0.125366004970E-17*tc*tc*tc*tc*tc*tc*tc
                        +0.214892175690E-20*tc*tc*tc*tc*tc*tc*tc*tc
                        -0.143880417820E-23*tc*tc*tc*tc*tc*tc*tc*tc*tc
                        +0.359608994810E-27*tc*tc*tc*tc*tc*tc*tc*tc*tc*tc;
   }
   volt= (double) milli_volt / 1000.0;
   return volt;
}

/********************************************************************************************/
/*																							*/
/* tce_VtoF: This function calculates temperature in deg F from volts at cold junction		*/
/*  reference temperature 'refTemp' (F).  The function returns the calculated temperature	*/
/*																							*/
/*							This function is valid for TYPE E thermocouple					*/
/*																							*/
/********************************************************************************************/
#ifdef _HPUX_SOURCE
	double tce_VtoF (double volts, double refTemp)
#else
	DllExport double __cdecl tce_VtoF (double volts, double refTemp)
#endif
{
	double refVolts;

	if (refTemp == 0.0) return 0.0;
	if (volts == 0.0) return 0.0;
	if (refTemp == 32.0)return tce_V32toF (volts);

	refVolts=tce_FtoV32(refTemp);
	if (refVolts == 0.0)return 0.0;

	return tce_V32toF((refVolts+volts));
}

/********************************************************************************************/
/*																							*/
/* tcj_V32toF: This function calculates temperature in deg F from volts when cold junction	*/
/*	reference temperature is 32 F (standard temprature) for type J TC. The function returns */
/*	calculated temperature.																	*/
/*																							*/
/*							This function is valid for TYPE J thermocouple					*/
/* 
   In this function, NIST equations refer to ITS-90 Table.

   This function calculates temperature in deg F referenced to 32 deg F from type J thermocouple. 
   Th input is variable 'volts' in volts.

      The quations were obtained from N.I.S.T web site on Jan 27, 2004.  The NIST equations are valid
      for temperatures from:
       -210 deg C to 1200 deg C (-346 deg F to 2192 deg F) or -8.095 millivolts to 69.553 millivolts.

      Accuracy of these NIST equations are as follows:
      -210 to    0 deg C (-346 to   32 deg F)   Accuracy= -0.05 to 0.03 deg C (-0.09  to  0.054 deg F)  
	     0 to  760 deg C (  32 to 1400 deg F)   Accuracy= -0.04 to 0.04 deg C (-0.072 to  0.072 deg F)
	   760 to 1200 deg C (1400 to 2192 deg F)   Accuracy= -0.04 to 0.03 deg C (-0.072 to  0.054 deg F)
 
   If 'volts' is outside of temperature range, function returns 0 deg F.

   'volts' is in volts

   dll version 'tcj_V32toF': written by Nuri Cankurt March 12, 2004
*/
/********************************************************************************************/
#ifdef _HPUX_SOURCE
	double tcj_V32toF (double volts)
#else
	DllExport double __cdecl tcj_V32toF (double volts)
#endif
{
   double tf;  /* deg F */
   double tc;  /* deg C */
   double milli_volt;
   double emf;

   tf=0.0;
   tc=0.0;	// Nuri Cankurt 5-17-2010
   if (volts < -8.09538e-3) return tf;
   if (volts > 69.55318e-3) return tf;
   if (volts == 0.0) return tf;

   milli_volt= volts * 1000.0;
   emf=milli_volt;

   if (emf >= -8.09538 &&  emf <= 0.0)  /* -210 deg C to 0 deg C (-346 to   32 deg F) */
   {
      tc=    0.0000000E+00 
            +1.9528268E+01*emf  
            -1.2286185E+00*emf*emf
            -1.0752178E+00*emf*emf*emf
            -5.9086933E-01*emf*emf*emf*emf
            -1.7256713E-01*emf*emf*emf*emf*emf
            -2.8131513E-02*emf*emf*emf*emf*emf*emf
            -2.3963370E-03*emf*emf*emf*emf*emf*emf*emf
            -8.3823321E-05*emf*emf*emf*emf*emf*emf*emf*emf; 
   }
   else if (emf > 0.0 && emf <=  42.919)  /* 0 deg C to 760 deg C (  32 to 1400 deg F) */
   {
      tc=    0.000000E+00 
			+1.978425E+01*emf 
            -2.001204E-01*emf*emf
            +1.036969E-02*emf*emf*emf 
            -2.549687E-04*emf*emf*emf*emf
            +3.585153E-06*emf*emf*emf*emf*emf
            -5.344285E-08*emf*emf*emf*emf*emf*emf 
            +5.099890E-10*emf*emf*emf*emf*emf*emf*emf; 
   }
   else if (emf > 42.919 && emf <=  69.55318)  /* 760 deg C to 1200 deg C ( 1400 to 2192 deg F) */
   {
      tc=   -3.11358187E+03 
			+3.00543684E+02*emf 
            -9.94773230E+00*emf*emf
            +1.70276630E-01*emf*emf*emf 
            -1.43033468E-03*emf*emf*emf*emf
            +4.73886084E-06*emf*emf*emf*emf*emf;
   }

   tf=tc*1.8 + 32.0;  /* deg F */

   return tf;
}

/********************************************************************************************/
/*																							*/
/* tcj_FtoV32: This function calculates the voltage (in volts) for a Type J thermocouple	*/
/*	when the cold juction reference temperature is at 32 F (standard temperature).			*/
/*	The function returns the calculated voltage.											*/
/*																							*/
/*							This function is valid for TYPE J thermocouple					*/
/* 
/* 
   In this function, NIST equations refer to ITS-90 Table.
 
   This function calculates volts referenced to 32 deg F from type K thermocouple at 'temp' deg F.
   the quation was obtained from N.I.S.T web site on Jan 27, 2004.

   This function is valid for temperatures from -210 deg C  to 1200 deg C (-346 deg F to 2192 deg F)
   If 'temp' is outside of this range, function returns 0 volts.

   'temp' is in deg F

    dll version 'tcj_FtoV32': written by Nuri Cankurt March 12, 2004
*/
#ifdef _HPUX_SOURCE
	double tcj_FtoV32 (double temp)
#else
	DllExport double __cdecl tcj_FtoV32 (double temp)
#endif
{
   double volt;
   double tf;  /* deg F */
   double tc;  /* deg C */
   double milli_volt;

   volt=0.0;
   milli_volt=0.0;	// Nuri Cankurt 5-17-2010

   if (temp < -346.0 || temp > 2192.0) return volt;
   if (temp == 0.0) return volt;
   if (temp == 32.0)return volt;

   tf= temp;
   tc=(tf-32.0)/1.8;

   if (tc >= -210.0 && tc <= 760.0)  /* For temps from -210 deg C to 760 deg C (-346 deg F to 1400 deg F) */
   {
      milli_volt     =   0.000000000000E+00
                        +0.503811878150E-01*tc
                        +0.304758369300E-04*tc*tc
                        -0.856810657200E-07*tc*tc*tc
                        +0.132281952950E-09*tc*tc*tc*tc
                        -0.170529583370E-12*tc*tc*tc*tc*tc
                        +0.209480906970E-15*tc*tc*tc*tc*tc*tc
                        -0.125383953360E-18*tc*tc*tc*tc*tc*tc*tc
                        +0.156317256970E-22*tc*tc*tc*tc*tc*tc*tc*tc;
    }
    else if (tc >760.0 && tc <= 1200.0) /* For temps from 760 degC to 1200 deg C ( 1400 deg F to 2192.0 deg F) */
    {
       milli_volt     =  0.296456256810E+03
                        -0.149761277860E+01*tc
                        +0.317871039240E-02*tc*tc
                        -0.318476867010E-05*tc*tc*tc
                        +0.157208190040E-08*tc*tc*tc*tc
                        -0.306913690560E-12*tc*tc*tc*tc*tc;
   }
   volt= (double) milli_volt / 1000.0;
   return volt;
}

/********************************************************************************************/
/*																							*/
/* tcj_VtoF: This function calculates temperature in deg F from volts at cold junction		*/
/*  reference temperature 'refTemp' (F).  The function returns the calculated temperature	*/
/*																							*/
/*							This function is valid for TYPE J thermocouple					*/
/*																							*/
/********************************************************************************************/
#ifdef _HPUX_SOURCE
	double tcj_VtoF (double volts, double refTemp)
#else
	DllExport double __cdecl tcj_VtoF (double volts, double refTemp)
#endif
{
	double refVolts;

	if (refTemp == 0.0) return 0.0;
	if (volts == 0.0) return 0.0;
	if (refTemp == 32.0)return tcj_V32toF (volts);

	refVolts=tcj_FtoV32(refTemp);
	if (refVolts == 0.0)return 0.0;

	return tcj_V32toF((refVolts+volts));
}

/********************************************************************************************/
/*																							*/
/* tcn_V32toF: This function calculates temperature in deg F from volts when cold junction	*/
/*	reference temperature is 32 F (standard temprature) for type N TC. The function returns */
/*	calculated temperature.																	*/
/*																							*/
/*							This function is valid for TYPE N thermocouple					*/
/* 
   In this function, NIST equations refer to ITS-90 Table.

   This function calculates temperature in deg F referenced to 32 deg F from type N thermocouple. 
   Th input is variable 'volts' in volts.

      The quations were obtained from N.I.S.T web site on Jan 27, 2004.  The NIST equations are valid
      for temperatures from:
       -200 deg C to 1300 deg C (-328 deg F to 2372 deg F) or -3.990 millivolts to 47.513 millivolts.

      Accuracy of these NIST equations are as follows:
      -200 to    0 deg C (-328 to   32 deg F)   Accuracy= -0.02 to 0.03 deg C (-0.036 to  0.054 deg F) 
	     0 to  600 deg C (  32 to 1112 deg F)   Accuracy= -0.02 to 0.03 deg C (-0.036 to  0.054 deg F) 
	   600 to 1300 deg C (1112 to 2372 deg F)   Accuracy= -0.04 to 0.02 deg C (-0.072  to 0.036 deg F)

   If 'volts' is outside of temperature range, function returns 0 deg F.

   'volts' is in volts

   dll version 'tcn_V32toF': written by Nuri Cankurt May 17, 2010
*/
/********************************************************************************************/
#ifdef _HPUX_SOURCE
	double tcn_V32toF (double volts)
#else
	DllExport double __cdecl tcn_V32toF (double volts)
#endif
{
   double tf;  /* deg F */
   double tc;  /* deg C */
   double milli_volt;
   double emf;

   tf=0.0;
   tc=0.0;
   if (volts < -3.990e-3) return tf;
   if (volts > 47.513e-3) return tf;
   if (volts == 0.0) return tf;

   milli_volt= volts * 1000.0;
   emf=milli_volt;

   if (emf >= -3.990 &&  emf <= 0.0)  /* -200 deg C to 0 deg C (-328 to   32 deg F) */
   {
      tc=    0.0000000E+00 
            +3.8436847E+01*emf  
            +1.1010485E+00*emf*emf
            +5.2229312E+00*emf*emf*emf
            +7.2060525E+00*emf*emf*emf*emf
            +5.8488586E+00*emf*emf*emf*emf*emf
            +2.7754916E+00*emf*emf*emf*emf*emf*emf
            +7.7075166E-01*emf*emf*emf*emf*emf*emf*emf
            +1.1582665E-01*emf*emf*emf*emf*emf*emf*emf*emf
			+7.3138868E-03*emf*emf*emf*emf*emf*emf*emf*emf*emf;
   }
   else if (emf > 0.0 && emf <=  20.613)  /* 0 deg C to 600 deg C (  32 to 1112 deg F) */
   {
      tc=    0.000000E+00 
			+3.86896E+01*emf 
            -1.08267E+00*emf*emf
            +4.70205E-02*emf*emf*emf
			-2.12169E-06*emf*emf*emf*emf
            -1.17272E-04*emf*emf*emf*emf*emf
            +5.39280E-06*emf*emf*emf*emf*emf*emf
            -7.98156E-08*emf*emf*emf*emf*emf*emf*emf; 
   }
   else if (emf > 20.613 && emf <=  47.513)  /* 600 deg C to 1300 deg C ( 1112 to 2372 deg F) */
   {
      tc=    1.972485E+01 
			+3.300943E+01*emf 
            -3.915159E-01*emf*emf
            +9.855391E-03*emf*emf*emf
			-1.274371E-04*emf*emf*emf*emf
            +7.767022E-07*emf*emf*emf*emf*emf;
   }

   tf=tc*1.8 + 32.0;  /* deg F */

   return tf;
}

/********************************************************************************************/
/*																							*/
/* tcn_FtoV32: This function calculates the voltage (in volts) for a Type N thermocouple	*/
/*	when the cold juction reference temperature is at 32 F (standard temperature).			*/
/*	The function returns the calculated voltage.											*/
/*																							*/
/*							This function is valid for TYPE N thermocouple					*/
/* 
/* 
   In this function, NIST equations refer to ITS-90 Table.
 
   This function calculates volts referenced to 32 deg F from type N thermocouple at 'temp' deg F.
   the quation was obtained from N.I.S.T web site on Jan 27, 2004.

   This function is valid for temperatures from -270 deg C  to 1300 deg C (-454 deg F to 2372 deg F)
   If 'temp' is outside of this range, function returns 0 volts.

   'temp' is in deg F

    dll version 'tcn_FtoV32': written by Nuri Cankurt May 17, 2010
*/
#ifdef _HPUX_SOURCE
	double tcn_FtoV32 (double temp)
#else
	DllExport double __cdecl tcn_FtoV32 (double temp)
#endif
{
   double volt;
   double tf;  /* deg F */
   double tc;  /* deg C */
   double milli_volt;

   volt=0.0;
   milli_volt=0.0;

   if (temp < -454.0 || temp > 2372.0) return volt;
   if (temp == 0.0) return volt;
   if (temp == 32.0)return volt;

   tf= temp;
   tc=(tf-32.0)/1.8;

   if (tc >= -270.0 && tc <= 0.0)  /* For temps from -270 deg C to 0 deg C (-454 deg F to 32 deg F) */
   {
      milli_volt     =   0.000000000000E+00
                        +0.261591059620E-01*tc
                        +0.109574842280E-04*tc*tc
                        -0.938411115540E-07*tc*tc*tc
                        -0.464120397590E-10*tc*tc*tc*tc
                        -0.263033577160E-11*tc*tc*tc*tc*tc
                        -0.226534380030E-13*tc*tc*tc*tc*tc*tc
                        -0.760893007910E-16*tc*tc*tc*tc*tc*tc*tc
                        -0.934196678350E-19*tc*tc*tc*tc*tc*tc*tc*tc;
    }
    else if (tc >0.0 && tc <= 1300.0) /* For temps from 0 degC to 1300 deg C ( 32 deg F to 2372.0 deg F) */
    {
       milli_volt     =  0.000000000000E+00
                        +0.259293946010E-01*tc
                        +0.157101418800E-04*tc*tc
                        +0.438256272370E-07*tc*tc*tc
                        -0.252611697940E-09*tc*tc*tc*tc
                        +0.643118193390E-12*tc*tc*tc*tc*tc
                        -0.100634715190E-14*tc*tc*tc*tc*tc*tc
                        +0.997453389920E-18*tc*tc*tc*tc*tc*tc*tc
                        -0.608632456070E-21*tc*tc*tc*tc*tc*tc*tc*tc
                        +0.208492293390E-24*tc*tc*tc*tc*tc*tc*tc*tc*tc
                        -0.306821961510E-28*tc*tc*tc*tc*tc*tc*tc*tc*tc*tc;
   }
   volt= (double) milli_volt / 1000.0;
   return volt;
}

/********************************************************************************************/
/*																							*/
/* tcn_VtoF: This function calculates temperature in deg F from volts at cold junction		*/
/*  reference temperature 'refTemp' (F).  The function returns the calculated temperature	*/
/*																							*/
/*							This function is valid for TYPE N thermocouple					*/
/*																							*/
/********************************************************************************************/
#ifdef _HPUX_SOURCE
	double tcn_VtoF (double volts, double refTemp)
#else
	DllExport double __cdecl tcn_VtoF (double volts, double refTemp)
#endif
{
	double refVolts;

	if (refTemp == 0.0) return 0.0;
	if (volts == 0.0) return 0.0;
	if (refTemp == 32.0)return tcn_V32toF (volts);

	refVolts=tcn_FtoV32(refTemp);
	if (refVolts == 0.0)return 0.0;

	return tcn_V32toF((refVolts+volts));
}

/********************************************************************************************/
/*																							*/
/* tcr_V32toF: This function calculates temperature in deg F from volts when cold junction	*/
/*	reference temperature is 32 F (standard temprature) for type R TC. The function returns */
/*	calculated temperature.																	*/
/*																							*/
/*							This function is valid for TYPE R thermocouple					*/
/* 
   In this function, NIST equations refer to ITS-90 Table.

   This function calculates temperature in deg F referenced to 32 deg F from type R thermocouple. 
   Th input is variable 'volts' in volts.

      The quations were obtained from N.I.S.T web site on Jan 27, 2004.  The NIST equations are valid
      for temperatures from:
       -50 deg C to 1768.1 deg C (-58 deg F to 3214.58 deg F) or -0.226 millivolts to 21.103 millivolts.

      Accuracy of these NIST equations are as follows:
       -50   to   250   deg C ( -58   to  482    deg F)   Accuracy= -0.02   to 0.02  deg C (-0.036  to 0.036  deg F) 
	   250   to  1200   deg C ( 482   to 2192    deg F)   Accuracy= -0.005  to 0.005 deg C (-0.009  to 0.009  deg F) 
	  1064   to  1664.5 deg C (1947.2 to 3028.1  deg F)   Accuracy= -0.0005 to 0.001 deg C (-0.0009 to 0.0018 deg F)
	  1664.5 to  1768.1 deg C (3028.1 to 3214.58 deg F)   Accuracy= -0.001  to 0.002 deg C (-0.0018 to 0.0036 deg F)

   If 'volts' is outside of temperature range, function returns 0 deg F.

   'volts' is in volts

   dll version 'tcr_V32toF': written by Nuri Cankurt May 17, 2010
*/
/********************************************************************************************/
#ifdef _HPUX_SOURCE
	double tcr_V32toF (double volts)
#else
	DllExport double __cdecl tcr_V32toF (double volts)
#endif
{
   double tf;  /* deg F */
   double tc;  /* deg C */
   double milli_volt;
   double emf;

   tf=0.0;
   tc=0.0;
   if (volts < -0.226e-3) return tf;
   if (volts > 21.103e-3) return tf;
   if (volts == 0.0) return tf;

   milli_volt= volts * 1000.0;
   emf=milli_volt;

   if (emf >= -0.226 &&  emf <= 1.923)  /* -50 deg C to 250 deg C (-58 to   482 deg F) */
   {
      tc=    0.0000000E+00 
            +1.8891380E+02*emf  
            -9.3835290E+01*emf*emf
            +1.3068619E+02*emf*emf*emf
            -2.2703580E+02*emf*emf*emf*emf
            +3.5145659E+02*emf*emf*emf*emf*emf
            -3.8953900E+02*emf*emf*emf*emf*emf*emf
            +2.8239471E+02*emf*emf*emf*emf*emf*emf*emf
            -1.2607281E+02*emf*emf*emf*emf*emf*emf*emf*emf
			+3.1353611E+01*emf*emf*emf*emf*emf*emf*emf*emf*emf
	        -3.3187769E+00*emf*emf*emf*emf*emf*emf*emf*emf*emf*emf;
   }
   else if (emf > 1.923 && emf <=  11.361)  /* 250 deg C to 1064 deg C ( 482 to 1947.2 deg F) */
   {
      tc=    1.334584505E+01
			+1.472644573E+02*emf 
            -1.844024844E+01*emf*emf
            +4.031129726E+00*emf*emf*emf
			-6.249428360E-01*emf*emf*emf*emf
            +6.468412046E-02*emf*emf*emf*emf*emf
            -4.458750426E-03*emf*emf*emf*emf*emf*emf
            +1.994710149E-04*emf*emf*emf*emf*emf*emf*emf 
			-5.313401790E-06*emf*emf*emf*emf*emf*emf*emf*emf
			+6.481976217E-08*emf*emf*emf*emf*emf*emf*emf*emf*emf; 
   }
   else if (emf > 11.361 && emf <=  19.739)  /* 1064 deg C to 1664.5 deg C ( 1947.2 to 3028.1 deg F) */
   {
      tc=   -8.199599416E+01
			+1.553962042E+02*emf 
            -8.342197663E+00*emf*emf
            +4.279433549E-01*emf*emf*emf
			-1.191577910E-02*emf*emf*emf*emf
            +1.492290091E-04*emf*emf*emf*emf*emf;
   }
   else if (emf > 19.739 && emf <=  21.103)  /* 1664.5 deg C to 1768.1 deg C ( 3028.1 to 3214.58 deg F) */
   {
      tc=    3.406177836E+04
			-7.023729171E+03*emf 
            +5.582903813E+02*emf*emf
            -1.952394635E+01*emf*emf*emf
			+2.560740231E-01*emf*emf*emf*emf;
   }

   tf=tc*1.8 + 32.0;  /* deg F */

   return tf;
}

/********************************************************************************************/
/*																							*/
/* tcr_FtoV32: This function calculates the voltage (in volts) for a Type R thermocouple	*/
/*	when the cold juction reference temperature is at 32 F (standard temperature).			*/
/*	The function returns the calculated voltage.											*/
/*																							*/
/*							This function is valid for TYPE R thermocouple					*/
/* 
/* 
   In this function, NIST equations refer to ITS-90 Table.
 
   This function calculates volts referenced to 32 deg F from type R thermocouple at 'temp' deg F.
   the quation was obtained from N.I.S.T web site on Jan 27, 2004.

   This function is valid for temperatures from -50 deg C  to 1768.1 deg C (-58 deg F to 3214.58 deg F)
   If 'temp' is outside of this range, function returns 0 volts.

   'temp' is in deg F

    dll version 'tcr_FtoV32': written by Nuri Cankurt May 17, 2010
*/
#ifdef _HPUX_SOURCE
	double tcr_FtoV32 (double temp)
#else
	DllExport double __cdecl tcr_FtoV32 (double temp)
#endif
{
   double volt;
   double tf;  /* deg F */
   double tc;  /* deg C */
   double milli_volt;

   volt=0.0;
   milli_volt=0.0;

   if (temp < -58.0 || temp > 3214.58) return volt;
   if (temp == 0.0) return volt;
   if (temp == 32.0)return volt;

   tf= temp;
   tc=(tf-32.0)/1.8;

   if (tc >= -50.0 && tc <= 1064.18)  /* For temps from -50 deg C to 1064.18 deg C (-58 deg F to 1947.52 deg F) */
   {
      milli_volt     =   0.000000000000E+00
                        +0.528961729765E-02*tc
                        +0.139166589782E-04*tc*tc
                        -0.238855693017E-07*tc*tc*tc
                        +0.356916001063E-10*tc*tc*tc*tc
                        -0.462347666298E-13*tc*tc*tc*tc*tc
                        +0.500777441034E-16*tc*tc*tc*tc*tc*tc
                        -0.373105886191E-19*tc*tc*tc*tc*tc*tc*tc
                        +0.157716482367E-22*tc*tc*tc*tc*tc*tc*tc*tc
						-0.281038625251E-26*tc*tc*tc*tc*tc*tc*tc*tc*tc;
    }
    else if (tc >1064.18 && tc <= 1664.5) /* For temps from 1064.18 degC to 1664.5 deg C ( 1947.52 deg F to 3028.1 deg F) */
    {
       milli_volt     =  0.295157925316E+01
                        -0.252061251332E-02*tc
                        +0.159564501865E-04*tc*tc
                        -0.764085947576E-08*tc*tc*tc
                        +0.205305291024E-11*tc*tc*tc*tc
                        -0.293359668173E-15*tc*tc*tc*tc*tc;
   }
    else if (tc >1664.5 && tc <= 1768.1) /* For temps from 1664.5 degC to 1768.1 deg C ( 3028.1 deg F to 3214.58 deg F) */
    {
       milli_volt     =  0.152232118209E+03
                        -0.268819888545E+00*tc
                        +0.171280280471E-03*tc*tc
                        -0.345895706453E-07*tc*tc*tc
                        -0.934633971046E-14*tc*tc*tc*tc;
   }

   volt= (double) milli_volt / 1000.0;
   return volt;
}

/********************************************************************************************/
/*																							*/
/* tcr_VtoF: This function calculates temperature in deg F from volts at cold junction		*/
/*  reference temperature 'refTemp' (F).  The function returns the calculated temperature	*/
/*																							*/
/*							This function is valid for TYPE R thermocouple					*/
/*																							*/
/********************************************************************************************/
#ifdef _HPUX_SOURCE
	double tcr_VtoF (double volts, double refTemp)
#else
	DllExport double __cdecl tcr_VtoF (double volts, double refTemp)
#endif
{
	double refVolts;

	if (refTemp == 0.0) return 0.0;
	if (volts == 0.0) return 0.0;
	if (refTemp == 32.0)return tcr_V32toF (volts);

	refVolts=tcr_FtoV32(refTemp);
	if (refVolts == 0.0)return 0.0;

	return tcr_V32toF((refVolts+volts));
}

/********************************************************************************************/
/*																							*/
/* tcs_V32toF: This function calculates temperature in deg F from volts when cold junction	*/
/*	reference temperature is 32 F (standard temprature) for type S TC. The function returns */
/*	calculated temperature.																	*/
/*																							*/
/*							This function is valid for TYPE S thermocouple					*/
/* 
   In this function, NIST equations refer to ITS-90 Table.

   This function calculates temperature in deg F referenced to 32 deg F from type S thermocouple. 
   Th input is variable 'volts' in volts.

      The quations were obtained from N.I.S.T web site on Jan 27, 2004.  The NIST equations are valid
      for temperatures from:
       -50 deg C to 1768.1 deg C (-58 deg F to 3214.58 deg F) or -0.226 millivolts to 21.103 millivolts.

      Accuracy of these NIST equations are as follows:
       -50   to   250   deg C ( -58   to  482    deg F)   Accuracy= -0.02   to 0.02   deg C (-0.036  to 0.036  deg F) 
	   250   to  1200   deg C ( 482   to 2192    deg F)   Accuracy= -0.01   to 0.01   deg C (-0.018  to 0.018  deg F) 
	  1064   to  1664.5 deg C (1947.2 to 3028.1  deg F)   Accuracy= -0.0002 to 0.0002 deg C (-0.0004 to 0.0004 deg F)
	  1664.5 to  1768.1 deg C (3028.1 to 3214.58 deg F)   Accuracy= -0.002  to 0.002  deg C (-0.0036 to 0.0036 deg F)

   If 'volts' is outside of temperature range, function returns 0 deg F.

   'volts' is in volts

   dll version 'tcs_V32toF': written by Nuri Cankurt May 17, 2010
*/
/********************************************************************************************/
#ifdef _HPUX_SOURCE
	double tcs_V32toF (double volts)
#else
	DllExport double __cdecl tcs_V32toF (double volts)
#endif
{
   double tf;  /* deg F */
   double tc;  /* deg C */
   double milli_volt;
   double emf;

   tf=0.0;
   tc=0.0;
   if (volts < -0.235e-3) return tf;
   if (volts > 18.693e-3) return tf;
   if (volts == 0.0) return tf;

   milli_volt= volts * 1000.0;
   emf=milli_volt;

   if (emf >= -0.235 &&  emf <= 1.874)  /* -50 deg C to 250 deg C (-58 to   482 deg F) */
   {
      tc=    0.00000000E+00 
            +1.84949460E+02*emf  
            -8.00504062E+01*emf*emf
            +1.02237430E+02*emf*emf*emf
            -1.52248592E+02*emf*emf*emf*emf
            +1.88821343E+02*emf*emf*emf*emf*emf
            -1.59085941E+02*emf*emf*emf*emf*emf*emf
            +8.23027880E+01*emf*emf*emf*emf*emf*emf*emf
			-2.34181944E+01*emf*emf*emf*emf*emf*emf*emf*emf
	        +2.79786260E+00*emf*emf*emf*emf*emf*emf*emf*emf*emf;
   }
   else if (emf > 1.874 && emf <=  10.332)  /* 250 deg C to 1064 deg C ( 482 to 1947.2 deg F) */
   {
      tc=    1.291507177E+01
			+1.466298863E+02*emf 
            -1.534713402E+01*emf*emf
            +3.145945973E+00*emf*emf*emf
			-4.163257839E-01*emf*emf*emf*emf
            +3.187963771E-02*emf*emf*emf*emf*emf
            -1.291637500E-03*emf*emf*emf*emf*emf*emf
            +2.183475087E-05*emf*emf*emf*emf*emf*emf*emf 
			-1.447379511E-07*emf*emf*emf*emf*emf*emf*emf*emf
			+8.211272125E-09*emf*emf*emf*emf*emf*emf*emf*emf*emf; 
   }
   else if (emf > 10.332 && emf <=  17.536)  /* 1064 deg C to 1664.5 deg C ( 1947.2 to 3028.1 deg F) */
   {
      tc=   -8.087801117E+01
			+1.621573104E+02*emf 
            -8.536869453E+00*emf*emf
            +4.719686976E-01*emf*emf*emf
			-1.441693666E-02*emf*emf*emf*emf
            +2.081618890E-04 *emf*emf*emf*emf*emf;
   }
   else if (emf > 17.536 && emf <=  18.693)  /* 1664.5 deg C to 1768.1 deg C ( 3028.1 to 3214.58 deg F) */
   {
      tc=    5.333875126E+04
			-1.235892298E+04*emf 
            +1.092657613E+03*emf*emf
            -4.265693686E+01*emf*emf*emf
			+6.247205420E-01*emf*emf*emf*emf;
   }

   tf=tc*1.8 + 32.0;  /* deg F */

   return tf;
}

/********************************************************************************************/
/*																							*/
/* tcs_FtoV32: This function calculates the voltage (in volts) for a Type S thermocouple	*/
/*	when the cold juction reference temperature is at 32 F (standard temperature).			*/
/*	The function returns the calculated voltage.											*/
/*																							*/
/*							This function is valid for TYPE S thermocouple					*/
/* 
/* 
   In this function, NIST equations refer to ITS-90 Table.
 
   This function calculates volts referenced to 32 deg F from type S thermocouple at 'temp' deg F.
   the quation was obtained from N.I.S.T web site on Jan 27, 2004.

   This function is valid for temperatures from -50 deg C  to 1768.1 deg C (-58 deg F to 3214.58 deg F)
   If 'temp' is outside of this range, function returns 0 volts.

   'temp' is in deg F

    dll version 'tcs_FtoV32': written by Nuri Cankurt May 17, 2010
*/
#ifdef _HPUX_SOURCE
	double tcs_FtoV32 (double temp)
#else
	DllExport double __cdecl tcs_FtoV32 (double temp)
#endif
{
   double volt;
   double tf;  /* deg F */
   double tc;  /* deg C */
   double milli_volt;

   volt=0.0;
   milli_volt=0.0;

   if (temp < -58.0 || temp > 3214.58) return volt;
   if (temp == 0.0) return volt;
   if (temp == 32.0)return volt;

   tf= temp;
   tc=(tf-32.0)/1.8;

   if (tc >= -50.0 && tc <= 1064.18)  /* For temps from -50 deg C to 1064.18 deg C (-58 deg F to 1947.52 deg F) */
   {
      milli_volt     =   0.000000000000E+00
                        +0.540313308631E-02*tc
                        +0.125934289740E-04*tc*tc
                        -0.232477968689E-07*tc*tc*tc
                        +0.322028823036E-10*tc*tc*tc*tc
                        -0.331465196389E-13*tc*tc*tc*tc*tc
                        +0.255744251786E-16*tc*tc*tc*tc*tc*tc
                        -0.125068871393E-19*tc*tc*tc*tc*tc*tc*tc
                        +0.271443176145E-23*tc*tc*tc*tc*tc*tc*tc*tc;
    }
    else if (tc >1064.18 && tc <= 1664.5) /* For temps from 1064.18 degC to 1664.5 deg C ( 1947.52 deg F to 3028.1 deg F) */
    {
       milli_volt     =  0.132900444085E+01
                        +0.334509311344E-02*tc
                        +0.654805192818E-05*tc*tc
                        -0.164856259209E-08*tc*tc*tc
                        +0.129989605174E-13*tc*tc*tc*tc;
   }
    else if (tc >1664.5 && tc <= 1768.1) /* For temps from 1664.5 degC to 1768.1 deg C ( 3028.1 deg F to 3214.58 deg F) */
    {
       milli_volt     =  0.146628232636E+03
                        -0.258430516752E+00*tc
                        +0.163693574641E-03*tc*tc
                        -0.330439046987E-07*tc*tc*tc
                        -0.943223690612E-14*tc*tc*tc*tc;
   }

   volt= (double) milli_volt / 1000.0;
   return volt;
}

/********************************************************************************************/
/*																							*/
/* tcs_VtoF: This function calculates temperature in deg F from volts at cold junction		*/
/*  reference temperature 'refTemp' (F).  The function returns the calculated temperature	*/
/*																							*/
/*							This function is valid for TYPE S thermocouple					*/
/*																							*/
/********************************************************************************************/
#ifdef _HPUX_SOURCE
	double tcs_VtoF (double volts, double refTemp)
#else
	DllExport double __cdecl tcs_VtoF (double volts, double refTemp)
#endif
{
	double refVolts;

	if (refTemp == 0.0) return 0.0;
	if (volts == 0.0) return 0.0;
	if (refTemp == 32.0)return tcs_V32toF (volts);

	refVolts=tcs_FtoV32(refTemp);
	if (refVolts == 0.0)return 0.0;

	return tcs_V32toF((refVolts+volts));
}

/********************************************************************************************/
/*																							*/
/* tct_V32toF: This function calculates temperature in deg F from volts when cold junction	*/
/*	reference temperature is 32 F (standard temprature) for type T TC. The function returns */
/*	calculated temperature.																	*/
/*																							*/
/*							This function is valid for TYPE T thermocouple					*/
/* 
   In this function, NIST equations refer to ITS-90 Table.

   This function calculates temperature in deg F referenced to 32 deg F from type T thermocouple. 
   Th input is variable 'volts' in volts.

      The quations were obtained from N.I.T.T web site on Jan 27, 2004.  The NIST equations are valid
      for temperatures from:
       -200 deg C to 400 deg C (-328 deg F to 752 deg F) or -5.603 millivolts to 20.872 millivolts.

      Accuracy of these NIST equations are as follows:
       -200 to   0 deg C ( -328 to  32 deg F)   Accuracy= -0.02 to 0.04 deg C (-0.036  to 0.072  deg F) 
	      0 to 400 deg C (   32 to 752 deg F)   Accuracy= -0.03 to 0.03 deg C (-0.054  to 0.054  deg F)

   If 'volts' is outside of temperature range, function returns 0 deg F.

   'volts' is in volts

   dll version 'tct_V32toF': written by Nuri Cankurt May 17, 2010
*/
/********************************************************************************************/
#ifdef _HPUX_SOURCE
	double tct_V32toF (double volts)
#else
	DllExport double __cdecl tct_V32toF (double volts)
#endif
{
   double tf;  /* deg F */
   double tc;  /* deg C */
   double milli_volt;
   double emf;

   tf=0.0;
   tc=0.0;
   if (volts < -5.603e-3) return tf;
   if (volts > 20.872e-3) return tf;
   if (volts == 0.0) return tf;

   milli_volt= volts * 1000.0;
   emf=milli_volt;

   if (emf >= -5.603 &&  emf <= 0.0)  /* -200 deg C to 0 deg C (-328 to  32 deg F) */
   {
      tc=    0.0000000E+00 
            +2.5949192E+01*emf  
            -2.1316967E-01*emf*emf
            +7.9018692E-01*emf*emf*emf
            +4.2527777E-01*emf*emf*emf*emf
            +1.3304473E-01*emf*emf*emf*emf*emf
            +2.0241446E-02*emf*emf*emf*emf*emf*emf
            +1.2668171E-03*emf*emf*emf*emf*emf*emf*emf;
   }
   else if (emf > 0.0 && emf <=  20.872)  /* 0 deg C to 400 deg C ( 32 to 752 deg F) */
   {
      tc=    0.000000E+00
			+2.592800E+01*emf 
            -7.602961E-01*emf*emf
            +4.637791E-02*emf*emf*emf
			-2.165394E-03*emf*emf*emf*emf
            +6.048144E-05*emf*emf*emf*emf*emf
            -7.293422E-07*emf*emf*emf*emf*emf*emf;
   }

   tf=tc*1.8 + 32.0;  /* deg F */

   return tf;
}

/********************************************************************************************/
/*																							*/
/* tct_FtoV32: This function calculates the voltage (in volts) for a Type T thermocouple	*/
/*	when the cold juction reference temperature is at 32 F (standard temperature).			*/
/*	The function returns the calculated voltage.											*/
/*																							*/
/*							This function is valid for TYPE T thermocouple					*/
/* 
/* 
   In this function, NIST equations refer to ITS-90 Table.
 
   This function calculates volts referenced to 32 deg F from type T thermocouple at 'temp' deg F.
   the quation was obtained from N.I.T.T web site on Jan 27, 2004.

   This function is valid for temperatures from -270 deg C  to 400 deg C (-454 deg F to 752 deg F)
   If 'temp' is outside of this range, function returns 0 volts.

   'temp' is in deg F

    dll version 'tct_FtoV32': written by Nuri Cankurt May 17, 2010
*/
#ifdef _HPUX_SOURCE
	double tct_FtoV32 (double temp)
#else
	DllExport double __cdecl tct_FtoV32 (double temp)
#endif
{
   double volt;
   double tf;  /* deg F */
   double tc;  /* deg C */
   double milli_volt;

   volt=0.0;
   milli_volt=0.0;

   if (temp < -454.0 || temp > 752.0) return volt;
   if (temp == 0.0) return volt;
   if (temp == 32.0)return volt;

   tf= temp;
   tc=(tf-32.0)/1.8;

   if (tc >= -270.0 && tc <= 0.0)  /* For temps from -270 deg C to 0 deg C (-454 deg F to 32 deg F) */
   {
      milli_volt     =   0.000000000000E+00
                        +0.387481063640E-01*tc
                        +0.441944343470E-04*tc*tc
                        +0.118443231050E-06*tc*tc*tc
                        +0.200329735540E-07*tc*tc*tc*tc
                        +0.901380195590E-09*tc*tc*tc*tc*tc
                        +0.226511565930E-10*tc*tc*tc*tc*tc*tc
                        +0.360711542050E-12*tc*tc*tc*tc*tc*tc*tc
                        +0.384939398830E-14*tc*tc*tc*tc*tc*tc*tc*tc
						+0.282135219250E-16*tc*tc*tc*tc*tc*tc*tc*tc*tc
						+0.142515947790E-18*tc*tc*tc*tc*tc*tc*tc*tc*tc*tc
						+0.487686622860E-21*tc*tc*tc*tc*tc*tc*tc*tc*tc*tc*tc
						+0.107955392700E-23*tc*tc*tc*tc*tc*tc*tc*tc*tc*tc*tc*tc
						+0.139450270620E-26*tc*tc*tc*tc*tc*tc*tc*tc*tc*tc*tc*tc*tc
						+0.797951539270E-30*tc*tc*tc*tc*tc*tc*tc*tc*tc*tc*tc*tc*tc*tc;
    }
    else if (tc >0.0 && tc <= 400.0) /* For temps from 0 degC to 400 deg C ( 32 deg F to 752 deg F) */
    {
       milli_volt     =  0.000000000000E+00
                        +0.387481063640E-01*tc
                        +0.332922278800E-04*tc*tc
                        +0.206182434040E-06*tc*tc*tc
                        -0.218822568460E-08*tc*tc*tc*tc
                        +0.109968809280E-10*tc*tc*tc*tc*tc
						-0.308157587720E-13*tc*tc*tc*tc*tc*tc
						+0.454791352900E-16*tc*tc*tc*tc*tc*tc*tc
						-0.275129016730E-19*tc*tc*tc*tc*tc*tc*tc*tc;
   }

   volt= (double) milli_volt / 1000.0;
   return volt;
}

/********************************************************************************************/
/*																							*/
/* tct_VtoF: This function calculates temperature in deg F from volts at cold junction		*/
/*  reference temperature 'refTemp' (F).  The function returns the calculated temperature	*/
/*																							*/
/*							This function is valid for TYPE T thermocouple					*/
/*																							*/
/********************************************************************************************/
#ifdef _HPUX_SOURCE
	double tct_VtoF (double volts, double refTemp)
#else
	DllExport double __cdecl tct_VtoF (double volts, double refTemp)
#endif
{
	double refVolts;

	if (refTemp == 0.0) return 0.0;
	if (volts == 0.0) return 0.0;
	if (refTemp == 32.0)return tct_V32toF (volts);

	refVolts=tct_FtoV32(refTemp);
	if (refVolts == 0.0)return 0.0;

	return tct_V32toF((refVolts+volts));
}

/********************************************************************************************/
/*
	Written by: Nuri Cankurt September 25, 2008

	rtd_PTgen_CtoOhm: This function calculates resistance in Ohms for Platinum RTD
	using The Callendar - van Dusen coefficients which are obtained by calibration.

    The equations are valid for -200 degC (-328 degF) to +850 degC (1562 degF)

    The equation were obtained as follows:

	The resistance values are calculated using DIN EN 60751 from Omega handbook
	RTD Tables  Omega Book Volume MM Page Z-251

			T >= 0 deg C					T< 0 deg C		
							
		R=R0 x (1+A x T + B x T^2)			R=R0 x [1 + A x T + B x T^2 + C x (T-100) x T^3]		
	Valid from -200 degC to +850 degC							
	Class A accuracy=	+/- (0.15 + 0.002 x |T|) deg C						
	Class B accuracy=	+/- (0.3 + 0.005 x |T|) deg C						
	Note=Omega table has error in coef C. It should be -4.183E-12 instead of -4.183E-13 							
	Where							
	T is temperature in deg C							
	R is resistance in OHMs							
	R0 is resistance in OHMs at 0 deg C	obtained from calibration
	A, B and C are constants and obtained from calibration

	NOTE: if T>= 0 deg C, Coeficient C is not used

	Calling parameters:

	tempInDegC:	Temperature in deg C
	R0 is resistance in OHMs at 0 deg C	obtained from calibration
	A, B and C are constants and obtained from calibration	

	if A, B and C are all set to zero, then standard equation is used for Pt RTD's with Alpha coefficient of 0.00385
	with the following parameters:
						
		T>=0					T<0
	A=	3.908300E-03			3.908300E-03	
	B=	-5.775000E-07			-5.775000E-07	
	C=	N/A						-4.183000E-12	

    R0 must still be provided. Normally R0 is 100, 500 or 1000 ohms

	Function returns resistance in Ohms.  If the values are outside the valid range, function returns 0 ohms.

*/
/********************************************************************************************/
#ifdef _HPUX_SOURCE
	double rtd_PTgen_CtoOhm (double tempInDegC, double R0, double A, double B, double C)
#else
	DllExport double __cdecl rtd_PTgen_CtoOhm (double tempInDegC, double R0, double A, double B, double C)
#endif
{
   double tc;  /* deg C */
   double ohm;
   double AA, BB, CC;

   ohm=0.0;
   tc=tempInDegC;

   if (tc == 0.0)return R0;
   if (tc < -200.0 || tc > 850.0) return 0.0;
   if (R0 <= 0.0) return 0.0;

    AA=A;
	BB=B;
	CC=C;

   if (A == 0.0 && B == 0.0 && C == 0.0)	// Use standard equation
   {
		AA=  3.908300E-03;
		BB= -5.775000E-07;
		CC= -4.183000E-12;
   }

   if (CC == 0.0)	// temps greater than 0 degC only
   {
		if (tc < 0.0)  /* -200 deg C to 0 deg C (-328 to   32 deg F) */
		{
			return 0.0;
		}
		else  /* 0 deg C to 850 deg C (  32 to 1562 deg F) */
		{
			ohm= 1.0 + AA*tc + BB*tc*tc;
			ohm=ohm*R0;
		}
   }
   else				// all temperatures
   {
		if (tc < 0.0)  /* -200 deg C to 0 deg C (-328 to   32 deg F) */
		{
			ohm= 1.0 + AA*tc + BB*tc*tc + CC*(tc-100.0)*tc*tc*tc;
			ohm=ohm*R0; 
		}
		else  /* 0 deg C to 850 deg C (  32 to 1562 deg F) */
		{
			ohm= 1.0 + AA*tc + BB*tc*tc;
			ohm=ohm*R0;
		}
   }

   return ohm;
}

/********************************************************************************************/
/*
	Written by: Nuri Cankurt March 10, 2004

	rtd_PTgenA385_CtoOhm2: This function calculates resistance in Ohms for Platinum RTD
	with Alpha coefficient of 0.00385

    The equations are valid for -200 degC (-328 degF) to +850 degC (1562 degF)

    The equation were obtained as follows:

	The resistance values are calculated using DIN EN 60751 from Omega handbook
	RTD Tables  Omega Book Volume MM Page Z-251

			T >= 0 deg C					T< 0 deg C		
							
		R=R0 x (1+A x T + B x T^2)			R=R0 x [1 + A x T + B x T^2 + C x (T-100) x T^3]		
							
	Class A accuracy=	+/- (0.15 + 0.002 x |T|) deg C						
	Class B accuracy=	+/- (0.3 + 0.005 x |T|) deg C						
	Note=Omega table has error in coef C. It should be -4.183E-12 instead of -4.183E-13 							
	Where							
	T is temperature in deg C							
	R is resistance in OHMs							
	R0 is resistance in OHMs at 0 deg C							
	Valid from -200 degC to +850 degC							
		T>=0					T<0
	A=	3.908300E-03			3.908300E-03	
	B=	-5.775000E-07			-5.775000E-07	
	C=	N/A						-4.183000E-12	
	R0=	100						100	

	Calling parameters:

	tempInDegC:	Temperature in deg C
	refResistanceInOhms:RTD reference resistance value at 0 degC.  Normally 100, 500 or 1000 ohms

	Function returns resiatance in Ohms.  If the values are outside the valid range, function returns 0 ohms.

*/
/********************************************************************************************/
#ifdef _HPUX_SOURCE
	double rtd_PTgenA385_CtoOhm2 (double tempInDegC, double refResistanceInOhms)
#else
	DllExport double __cdecl rtd_PTgenA385_CtoOhm2 (double tempInDegC, double refResistanceInOhms)
#endif
{
   double tc;  /* deg C */
   double ohm;
   double AA,BB,CC;

   ohm=0.0;
   tc=tempInDegC;

   if (tc == 0.0)return refResistanceInOhms;
   if (tc < -200.0 || tc > 850.0) return 0.0;
   if (refResistanceInOhms <= 0.0) return 0.0;

	AA=  3.908300E-03;
	BB= -5.775000E-07;
	CC= -4.183000E-12;

   if (tc < 0.0)  /* -200 deg C to 0 deg C (-328 to   32 deg F) */
   {
		ohm= 1.0 + AA*tc + BB*tc*tc + CC*(tc-100.0)*tc*tc*tc;
		ohm=ohm*refResistanceInOhms; 
   }
   else  /* 0 deg C to 850 deg C (  32 to 1562 deg F) */
   {
		ohm= 1.0 + AA*tc + BB*tc*tc;
		ohm=ohm*refResistanceInOhms;
   }

   return ohm;
}

/********************************************************************************************/
/*
	Written by: Nuri Cankurt March 10, 2004

	rtd_PTgenA385_OhmToC2: This function calculates temperature in deg C for Platinum RTD
	with Alpha coefficient of 0.00385

    The equations are valid for -200 degC (-328 degF) to +850 degC (1562 degF)

    The equation were obtained as follows:

	The resistance values were first calculated using DIN EN 60751 from Omega handbook
	RTD Tables  Omega Book Volume MM Page Z-251

			T >= 0 deg C					T< 0 deg C		
							
		R=R0 x (1+A x T + B x T^2)			R=R0 x [1 + A x T + B x T^2 + C x (T-100) x T^3]		
							
	Class A accuracy=	+/- (0.15 + 0.002 x |T|) deg C						
	Class B accuracy=	+/- (0.3 + 0.005 x |T|) deg C						
	Note=Omega table has error in coef C. It should be -4.183E-12 instead of -4.183E-13 							
	Where							
	T is temperature in deg C							
	R is resistance in OHMs							
	R0 is resistance in OHMs at 0 deg C							
	Valid from -200 degC to +850 degC							
		T>=0					T<0
	A=	3.908300E-03			3.908300E-03	
	B=	-5.775000E-07			-5.775000E-07	
	C=	N/A						-4.183000E-12	
	R0=	100						100	

    Then from the values above R/R0 were calculated then data was curve fitted to 5th order polynomial
	as:

	Y=Temperature (deg C)
	X=R/R0

    Y=C0 + C1*X + C2*X^2 + C3*X^3 + C4*X^4 + C5*X^5

    The maximum absolute error is 0.002 degC (0.004 degF) for temperatures from -200 degC (-328 degF) to +850 degC (1562 degF)


	Calling parameters:

	resistanceInOhms:		Resistance value for which temperature is desired
	refResistanceInOhms:	RTD reference resistance value at 0 degC.  Normally 100, 500 or 1000 ohms

	Function returns temperature in deg C.  If the values are outside the valid range, function returns 0 deg C

*/
/********************************************************************************************/
#ifdef _HPUX_SOURCE
	double rtd_PTgenA385_OhmToC2 (double resistanceInOhms, double refResistanceInOhms)
#else
	DllExport double __cdecl rtd_PTgenA385_OhmToC2 (double resistanceInOhms, double refResistanceInOhms)
#endif
{
   double tc;  /* deg C */
   double rr;

   tc=0.0;
   if (resistanceInOhms <= 0.0) return tc;
   if (refResistanceInOhms <= 0.0) return tc;

   rr=resistanceInOhms/refResistanceInOhms;

   if(rr < 0.1852007 || rr > 3.90481126) return tc;

   if(rr == 1.0)return 0.0;

   if (rr < 1.0)  /* -200 deg C to 0 deg C (-328 to   32 deg F) */
   {
      tc=    -242.01992875749800 
             +222.28124921393300*rr  
             +2.58588550459902E+01*rr*rr
             -4.82604183294947E+00*rr*rr*rr
             -2.81833863095458E+00*rr*rr*rr*rr
             +1.52425906434655E+00*rr*rr*rr*rr*rr;
 
   }
   else  /* 0 deg C to 850 deg C (  32 to 1562 deg F) */
   {
      tc=    -247.20591701840700 
             +2.39445086335996E+02*rr  
             +6.74989295331761E+00*rr*rr
             +1.10818800562992E+00*rr*rr*rr
             -1.23766399221495E-01*rr*rr*rr*rr
             +2.43332658428699E-02*rr*rr*rr*rr*rr;
   }

   return tc;
}

/********************************************************************************************/
/*
	Written by: Nuri Cankurt September 25, 2008

	rtd_PTgen_OhmToC: This function calculates temperature in deg C for Platinum RTD
	using The Callendar - van Dusen coefficients which are obtained by calibration.

    The equations are valid for -200 degC (-328 degF) to +850 degC (1562 degF)

    The equation were obtained as follows:

	The resistance values are calculated using DIN EN 60751 from Omega handbook
	RTD Tables  Omega Book Volume MM Page Z-251

			T >= 0 deg C					T< 0 deg C		
							
		R=R0 x (1+A x T + B x T^2)			R=R0 x [1 + A x T + B x T^2 + C x (T-100) x T^3]		
	Valid from -200 degC to +850 degC							
	Class A accuracy=	+/- (0.15 + 0.002 x |T|) deg C						
	Class B accuracy=	+/- (0.3 + 0.005 x |T|) deg C						
	Note=Omega table has error in coef C. It should be -4.183E-12 instead of -4.183E-13 							
	Where							
	T is temperature in deg C							
	R is resistance in OHMs							
	R0 is resistance in OHMs at 0 deg C	obtained from calibration
	A, B and C are constants and obtained from calibration

	NOTE: if T>= 0 deg C, Coeficient C is not used

	Calling parameters:

	resistanceInOhms:		Resistance value for which temperature is desired
	R0 is resistance in OHMs at 0 deg C	obtained from calibration
	A, B and C are constants and obtained from calibration	

	if A, B and C are all set to zero, then standard equation is used for Pt RTD's with Alpha coefficient of 0.00385
	with the following parameters:
						
		T>=0					T<0
	A=	3.908300E-03			3.908300E-03	
	B=	-5.775000E-07			-5.775000E-07	
	C=	N/A						-4.183000E-12

	If C is zero, temperatures greater than 0 deg C are calculated (R >= R0). In this case the temperatures are exact solution to
	quadratic equation.  When C is zero and temperatures are less than 0 deg C ( R < R0), there is no solution and zero deg F is returned.

	When C is non-zero, the temperatures are calculated by trial error.  In this case the maximum error is +/- 0.00032 deg F (+/- 0.00018 deg C)

    R0 must still be provided. Normally R0 is 100, 500 or 1000 ohms

	Function returns temperature in deg C.  If the values are outside the valid range, function returns 0 deg C

*/
/********************************************************************************************/
#ifdef _HPUX_SOURCE
	double rtd_PTgen_OhmToC (double resistanceInOhms, double R0, double A, double B, double C)
#else
	DllExport double __cdecl rtd_PTgen_OhmToC (double resistanceInOhms, double R0, double A, double B, double C)
#endif
{
   double tc;  /* deg C */
   double rr;
   double AA, BB, CC;
   double xx;
   double xhi, xlo, xguess;
   double ycalc, diff;
   long	  icount;

   tc=0.0;
   if (resistanceInOhms <= 0.0) return 0.0;
   if (R0 <= 0.0) return 0.0;

   rr=resistanceInOhms/R0;

   if(rr < 0.1852007 || rr > 3.90481126) return 0.0;

   if(rr == 1.0)return 0.0;

	AA=A;
	BB=B;
	CC=C;

   if (A == 0.0 && B == 0.0 && C == 0.0)
   {
		AA=  3.908300E-03;
		BB= -5.775000E-07;
		CC= -4.183000E-12;
   }

   if (CC == 0.0)		// only temps greater than 0.0 degC
   {
		if (rr < 1.0)  /* -200 deg C to 0 deg C (-328 to   32 deg F) */
		{
			return 0.0;	// not valid for temps less than 0 deg C
			
		}
		else  /* 0 deg C to 850 deg C (  32 to 1562 deg F) */
		{
			if(BB == 0.0)return 0.0;
			xx=AA*AA - 4.0*BB + 4.0*BB*rr;
			if(xx < 0.0)return 0.0;

			tc= ( -AA + sqrt(xx))/(2.0*BB);
			if(tc < 0.0)return 0.0;
			return tc;
		}
   }
   else	// all temperatures
   {
		if (rr < 1.0)  // -200 deg C to 0 deg C (-328 to   32 deg F).... Itereative solution
		{
			xhi=850.0;	// deg C
			xlo=-200.0;	// deg C
			icount=0;
			while (icount < 30)
			{
				icount=icount+1;
				xguess=(xhi+xlo)/2.0;	// deg C 
				ycalc=rtd_PTgen_CtoOhm (xguess, R0, AA, BB, CC);	// ohms
				if (ycalc == 0.0)return 0.0;
				diff=ycalc-resistanceInOhms;
				if (fabs(diff/resistanceInOhms) <= 1.0e-6)	// error is less than 0.0003 deg C (< 0.0005 deg F)
				{
						tc=xguess;		// deg C 
						return tc;
				}
				if (diff < 0.0)
				{
					xlo=xguess;	// deg C
				}
				else
				{
					xhi=xguess;	// deg C
				}
			}
		 
		}
		else  /* 0 deg C to 850 deg C (  32 to 1562 deg F) */
		{
			if(BB == 0.0)return 0.0;
			xx=AA*AA - 4.0*BB + 4.0*BB*rr;
			if(xx < 0.0)return 0.0;

			tc= ( -AA + sqrt(xx))/(2.0*BB);
			if(tc < 0.0)return 0.0;
			return tc;
		}
   }
 
   return 0.0;
}

/********************************************************************************************/
/*
	Written by: Nuri Cankurt March 10, 2004

	rtd_PT100A385_OhmToC: This function calculates temperature in deg C for 100 ohm Platinum RTD
	with Alpha coefficient of 0.00385

    The equations are valid for -200 degC (-328 degF) to +850 degC (1562 degF)


	Calling parameters:

	resistanceInOhms:		Resistance value for which temperature is desired

	Function returns temperature in deg C.  If the values are outside the valid range, function returns 0 deg C

*/
/********************************************************************************************/
#ifdef _HPUX_SOURCE
	double rtd_PT100A385_OhmToC (double resistanceInOhms)
#else
	DllExport double __cdecl rtd_PT100A385_OhmToC (double resistanceInOhms)
#endif
{
   double tc;  /* deg C */

   if (resistanceInOhms == 100.0) return 0.0;

//   tc=rtd_PTgenA385_OhmToC2 (resistanceInOhms, 100.0);
	 tc=rtd_PTgen_OhmToC (resistanceInOhms, 100.0, 0.0, 0.0, 0.0);	// Nuri Cankurt 9-29-2008

   return tc;
}

/********************************************************************************************/
/*
	Written by: Nuri Cankurt March 10, 2004

	rtd_PT500A385_OhmToC: This function calculates temperature in deg C for 500 ohm Platinum RTD
	with Alpha coefficient of 0.00385

    The equations are valid for -200 degC (-328 degF) to +850 degC (1562 degF)


	Calling parameters:

	resistanceInOhms:		Resistance value for which temperature is desired

	Function returns temperature in deg C.  If the values are outside the valid range, function returns 0 deg C

*/
/********************************************************************************************/
#ifdef _HPUX_SOURCE
	double rtd_PT500A385_OhmToC (double resistanceInOhms)
#else
	DllExport double __cdecl rtd_PT500A385_OhmToC (double resistanceInOhms)
#endif
{
   double tc;  /* deg C */

   if (resistanceInOhms == 500.0) return 0.0;

//   tc=rtd_PTgenA385_OhmToC2 (resistanceInOhms, 500.0);
   tc=rtd_PTgen_OhmToC (resistanceInOhms, 500.0, 0.0, 0.0, 0.0);	// Nuri Cankurt 9-29-2008

   return tc;
}

/********************************************************************************************/
/*
	Written by: Nuri Cankurt March 10, 2004

	rtd_PT1000A385_OhmToC: This function calculates temperature in deg C for 1000 ohm Platinum RTD
	with Alpha coefficient of 0.00385

    The equations are valid for -200 degC (-328 degF) to +850 degC (1562 degF)


	Calling parameters:

	resistanceInOhms:		Resistance value for which temperature is desired

	Function returns temperature in deg C.  If the values are outside the valid range, function returns 0 deg C

*/
/********************************************************************************************/
#ifdef _HPUX_SOURCE
	double rtd_PT1000A385_OhmToC (double resistanceInOhms)
#else
	DllExport double __cdecl rtd_PT1000A385_OhmToC (double resistanceInOhms)
#endif
{
   double tc;  /* deg C */

   if (resistanceInOhms == 1000.0) return 0.0;

//   tc=rtd_PTgenA385_OhmToC2 (resistanceInOhms, 1000.0);
   tc=rtd_PTgen_OhmToC (resistanceInOhms, 1000.0, 0.0, 0.0, 0.0);	// Nuri Cankurt 9-29-2008

   return tc;
}

/********************************************************************************************/
/*
	Written by: Nuri Cankurt March 10, 2004

	rtd_PT100A385_CtoOhm: This function calculates resistance in ohm for 100 ohm Platinum RTD
	with Alpha coefficient of 0.00385

    The equations are valid for -200 degC (-328 degF) to +850 degC (1562 degF)


	Calling parameters:

	tempInDegC:		temperature in deg C

	Function return resistance in ohm.  If the values are outside the valid range, function returns 0 ohm

*/
/********************************************************************************************/
#ifdef _HPUX_SOURCE
	double rtd_PT100A385_CtoOhm (double tempInDegC)
#else
	DllExport double __cdecl rtd_PT100A385_CtoOhm (double tempInDegC)
#endif
{
   double ohm;  

   if (tempInDegC == 0.0) return 100.0;

//   ohm=rtd_PTgenA385_CtoOhm2 (tempInDegC, 100.0);
   ohm=rtd_PTgen_CtoOhm (tempInDegC, 100.0, 0.0, 0.0, 0.0);	// Nuri Cankurt 9-29-2008

   return ohm;
}

/********************************************************************************************/
/*
	Written by: Nuri Cankurt March 10, 2004

	rtd_PT500A385_CtoOhm: This function calculates resistance in ohm for 500 ohm Platinum RTD
	with Alpha coefficient of 0.00385

    The equations are valid for -200 degC (-328 degF) to +850 degC (1562 degF)


	Calling parameters:

	tempInDegC:		temperature in deg C

	Function return resistance in ohm.  If the values are outside the valid range, function returns 0 ohm

*/
/********************************************************************************************/
#ifdef _HPUX_SOURCE
	double rtd_PT500A385_CtoOhm (double tempInDegC)
#else
	DllExport double __cdecl rtd_PT500A385_CtoOhm (double tempInDegC)
#endif
{
   double ohm;  

   if (tempInDegC == 0.0) return 500.0;

//   ohm=rtd_PTgenA385_CtoOhm2 (tempInDegC, 500.0);
   ohm=rtd_PTgen_CtoOhm (tempInDegC, 500.0, 0.0, 0.0, 0.0);	// Nuri Cankurt 9-29-2008

   return ohm;
}

/********************************************************************************************/
/*
	Written by: Nuri Cankurt March 10, 2004

	rtd_PT1000A385_CtoOhm: This function calculates resistance in ohm for 1000 ohm Platinum RTD
	with Alpha coefficient of 0.00385

    The equations are valid for -200 degC (-328 degF) to +850 degC (1562 degF)


	Calling parameters:

	tempInDegC:		temperature in deg C

	Function return resistance in ohm.  If the values are outside the valid range, function returns 0 ohm

*/
/********************************************************************************************/
#ifdef _HPUX_SOURCE
	double rtd_PT1000A385_CtoOhm (double tempInDegC)
#else
	DllExport double __cdecl rtd_PT1000A385_CtoOhm (double tempInDegC)
#endif
{
   double ohm;  

   if (tempInDegC == 0.0) return 1000.0;

//   ohm=rtd_PTgenA385_CtoOhm2 (tempInDegC, 1000.0);
   ohm=rtd_PTgen_CtoOhm (tempInDegC, 1000.0, 0.0, 0.0, 0.0);	// Nuri Cankurt 9-29-2008

   return ohm;
}

/********************************************************************************************/
/*
	Written by: Nuri Cankurt September 25, 2008

	rtd_PTgen_FtoOhm: This function calculates resistance in Ohms for Platinum RTD
	using The Callendar - van Dusen coefficients which are obtained by calibration.

    The equations are valid for -200 degC (-328 degF) to +850 degC (1562 degF)

    The equation were obtained as follows:

	The resistance values are calculated using DIN EN 60751 from Omega handbook
	RTD Tables  Omega Book Volume MM Page Z-251

			T >= 0 deg C					T< 0 deg C		
							
		R=R0 x (1+A x T + B x T^2)			R=R0 x [1 + A x T + B x T^2 + C x (T-100) x T^3]		
	Valid from -200 degC to +850 degC							
	Class A accuracy=	+/- (0.15 + 0.002 x |T|) deg C						
	Class B accuracy=	+/- (0.3 + 0.005 x |T|) deg C						
	Note=Omega table has error in coef C. It should be -4.183E-12 instead of -4.183E-13 							
	Where							
	T is temperature in deg C							
	R is resistance in OHMs							
	R0 is resistance in OHMs at 0 deg C	obtained from calibration
	A, B and C are constants and obtained from calibration

	NOTE: if T>= 0 deg C, Coeficient C is not used

	Calling parameters:

	tempInDegF:	Temperature in deg 
	R0 is resistance in OHMs at 0 deg C	obtained from calibration
	A, B and C are constants and obtained from calibration	

	if A, B and C are all set to zero, then standard equation is used for Pt RTD's with Alpha coefficient of 0.00385
	with the following parameters:
						
		T>=0					T<0
	A=	3.908300E-03			3.908300E-03	
	B=	-5.775000E-07			-5.775000E-07	
	C=	N/A						-4.183000E-12	

    R0 must still be provided. Normally R0 is 100, 500 or 1000 ohms

	Function returns resistance in Ohms.  If the values are outside the valid range, function returns 0 ohms.

*/
/********************************************************************************************/
#ifdef _HPUX_SOURCE
	double rtd_PTgen_FtoOhm (double tempInDegF, double R0, double A, double B, double C)
#else
	DllExport double __cdecl rtd_PTgen_FtoOhm (double tempInDegF, double R0, double A, double B, double C)
#endif
{
   double tf;  /* deg F */
   double tc;  /* deg C */
   double ohm;
   double AA, BB, CC;

   ohm=0.0;
   tf=tempInDegF;
   tc=(tf-32.0)/1.8;

   if (tf == 0.0)return 0.0;
   if (tc < -200.0 || tc > 850.0) return 0.0;
   if (R0 <= 0.0) return 0.0;
   if (tf == 32.0)return R0;

    AA=A;
	BB=B;
	CC=C;

   if (A == 0.0 && B == 0.0 && C == 0.0)	// Use standard equation
   {
		AA=  3.908300E-03;
		BB= -5.775000E-07;
		CC= -4.183000E-12;
   }

   if (CC == 0.0)	// temps greater than 0 degC only
   {
		if (tc < 0.0)  /* -200 deg C to 0 deg C (-328 to   32 deg F) */
		{
			return 0.0;
		}
		else  /* 0 deg C to 850 deg C (  32 to 1562 deg F) */
		{
			ohm= 1.0 + AA*tc + BB*tc*tc;
			ohm=ohm*R0;
		}
   }
   else				// all temperatures
   {
		if (tc < 0.0)  /* -200 deg C to 0 deg C (-328 to   32 deg F) */
		{
			ohm= 1.0 + AA*tc + BB*tc*tc + CC*(tc-100.0)*tc*tc*tc;
			ohm=ohm*R0; 
		}
		else  /* 0 deg C to 850 deg C (  32 to 1562 deg F) */
		{
			ohm= 1.0 + AA*tc + BB*tc*tc;
			ohm=ohm*R0;
		}
   }

   return ohm;
}

/********************************************************************************************/
/*
	Written by: Nuri Cankurt March 10, 2004

	rtd_PTgenA385_FtoOhm2: This function calculates resistance in Ohms for Platinum RTD
	with Alpha coefficient of 0.00385

    The equations are valid for -200 degC (-328 degF) to +850 degC (1562 degF)

    The equation were obtained as follows:

	The resistance values are calculated using DIN EN 60751 from Omega handbook
	RTD Tables  Omega Book Volume MM Page Z-251

			T >= 0 deg C					T< 0 deg C		
							
		R=R0 x (1+A x T + B x T^2)			R=R0 x [1 + A x T + B x T^2 + C x (T-100) x T^3]		
							
	Class A accuracy=	+/- (0.15 + 0.002 x |T|) deg C						
	Class B accuracy=	+/- (0.3 + 0.005 x |T|) deg C						
	Note=Omega table has error in coef C. It should be -4.183E-12 instead of -4.183E-13 							
	Where							
	T is temperature in deg C							
	R is resistance in OHMs							
	R0 is resistance in OHMs at 0 deg C							
	Valid from -200 degC to +850 degC							
		T>=0					T<0
	A=	3.908300E-03			3.908300E-03	
	B=	-5.775000E-07			-5.775000E-07	
	C=	N/A						-4.183000E-12	
	R0=	100						100	

	Calling parameters:

	tempInDegF:	Temperature in deg 
	refResistanceInOhms:RTD reference resistance value at 0 degC.  Normally 100, 500 or 1000 ohms

	Function returns resiatance in Ohms.  If the values are outside the valid range, function returns 0 ohms.

*/
/********************************************************************************************/
#ifdef _HPUX_SOURCE
	double rtd_PTgenA385_FtoOhm2 (double tempInDegF, double refResistanceInOhms)
#else
	DllExport double __cdecl rtd_PTgenA385_FtoOhm2 (double tempInDegF, double refResistanceInOhms)
#endif
{
   double tf;  /* deg F */
   double tc;  /* deg C */
   double ohm;
   double AA,BB,CC;

   ohm=0.0;
   tf=tempInDegF;
   tc=(tf-32.0)/1.8;

   if (tf == 0.0)return 0.0;
   if (tc < -200.0 || tc > 850.0) return 0.0;
   if (refResistanceInOhms <= 0.0) return 0.0;
   if (tf == 32.0)return refResistanceInOhms;

	AA=  3.908300E-03;
	BB= -5.775000E-07;
	CC= -4.183000E-12;

   if (tc < 0.0)  /* -200 deg C to 0 deg C (-328 to   32 deg F) */
   {
		ohm= 1.0 + AA*tc + BB*tc*tc + CC*(tc-100.0)*tc*tc*tc;
		ohm=ohm*refResistanceInOhms; 
   }
   else  /* 0 deg C to 850 deg C (  32 to 1562 deg F) */
   {
		ohm= 1.0 + AA*tc + BB*tc*tc;
		ohm=ohm*refResistanceInOhms;
   }

   return ohm;
}

/********************************************************************************************/
/*
	Written by: Nuri Cankurt March 10, 2004

	rtd_PTgenA385_OhmToF2: This function calculates temperature in deg F for Platinum RTD
	with Alpha coefficient of 0.00385

    The equations are valid for -200 degC (-328 degF) to +850 degC (1562 degF)

    The equation were obtained as follows:

	The resistance values were first calculated using DIN EN 60751 from Omega handbook
	RTD Tables  Omega Book Volume MM Page Z-251

			T >= 0 deg C					T< 0 deg C		
							
		R=R0 x (1+A x T + B x T^2)			R=R0 x [1 + A x T + B x T^2 + C x (T-100) x T^3]		
							
	Class A accuracy=	+/- (0.15 + 0.002 x |T|) deg C						
	Class B accuracy=	+/- (0.3 + 0.005 x |T|) deg C						
	Note=Omega table has error in coef C. It should be -4.183E-12 instead of -4.183E-13 							
	Where							
	T is temperature in deg C							
	R is resistance in OHMs							
	R0 is resistance in OHMs at 0 deg C							
	Valid from -200 degC to +850 degC							
		T>=0					T<0
	A=	3.908300E-03			3.908300E-03	
	B=	-5.775000E-07			-5.775000E-07	
	C=	N/A						-4.183000E-12	
	R0=	100						100	

    Then from the values above R/R0 were calculated then data was curve fitted to 5th order polynomial
	as:

	Y=Temperature (deg C)
	X=R/R0

    Y=C0 + C1*X + C2*X^2 + C3*X^3 + C4*X^4 + C5*X^5

    The maximum absolute error is 0.002 degC (0.004 degF) for temperatures from -200 degC (-328 degF) to +850 degC (1562 degF)


	Calling parameters:

	resistanceInOhms:		Resistance value for which temperature is desired
	refResistanceInOhms:	RTD reference resistance value at 0 degC.  Normally 100, 500 or 1000 ohms

	Function returns temperature in deg F.  If the values are outside the valid range, function returns 0 deg F

*/
/********************************************************************************************/
#ifdef _HPUX_SOURCE
	double rtd_PTgenA385_OhmToF2 (double resistanceInOhms, double refResistanceInOhms)
#else
	DllExport double __cdecl rtd_PTgenA385_OhmToF2 (double resistanceInOhms, double refResistanceInOhms)
#endif
{
   double tf;  /* deg F */
   double tc;  /* deg C */
   double rr;

   tf=0.0;
   if (resistanceInOhms <= 0.0) return tf;
   if (refResistanceInOhms <= 0.0) return tf;

   rr=resistanceInOhms/refResistanceInOhms;

   if(rr < 0.1852007 || rr > 3.90481126) return tf;

   if(rr == 1.0)return 32.0;

   if (rr < 1.0)  /* -200 deg C to 0 deg C (-328 to   32 deg F) */
   {
      tc=    -242.01992875749800 
             +222.28124921393300*rr  
             +2.58588550459902E+01*rr*rr
             -4.82604183294947E+00*rr*rr*rr
             -2.81833863095458E+00*rr*rr*rr*rr
             +1.52425906434655E+00*rr*rr*rr*rr*rr;
 
   }
   else  /* 0 deg C to 850 deg C (  32 to 1562 deg F) */
   {
      tc=    -247.20591701840700 
             +2.39445086335996E+02*rr  
             +6.74989295331761E+00*rr*rr
             +1.10818800562992E+00*rr*rr*rr
             -1.23766399221495E-01*rr*rr*rr*rr
             +2.43332658428699E-02*rr*rr*rr*rr*rr;
   }

   tf=tc*1.8 + 32.0;  /* deg F */

   return tf;
}

/********************************************************************************************/
/*
	Written by: Nuri Cankurt September 25, 2008

	rtd_PTgen_OhmToF: This function calculates temperature in deg F for Platinum RTD
	using The Callendar - van Dusen coefficients which are obtained by calibration.

    The equations are valid for -200 degC (-328 degF) to +850 degC (1562 degF)

    The equation were obtained as follows:

	The resistance values are calculated using DIN EN 60751 from Omega handbook
	RTD Tables  Omega Book Volume MM Page Z-251

			T >= 0 deg C					T< 0 deg C		
							
		R=R0 x (1+A x T + B x T^2)			R=R0 x [1 + A x T + B x T^2 + C x (T-100) x T^3]		
	Valid from -200 degC to +850 degC							
	Class A accuracy=	+/- (0.15 + 0.002 x |T|) deg C						
	Class B accuracy=	+/- (0.3 + 0.005 x |T|) deg C						
	Note=Omega table has error in coef C. It should be -4.183E-12 instead of -4.183E-13 							
	Where							
	T is temperature in deg C							
	R is resistance in OHMs							
	R0 is resistance in OHMs at 0 deg C	obtained from calibration
	A, B and C are constants and obtained from calibration

	NOTE: if T>= 0 deg C, Coeficient C is not used

	Calling parameters:

	resistanceInOhms:		Resistance value for which temperature is desired
	R0 is resistance in OHMs at 0 deg C	obtained from calibration
	A, B and C are constants and obtained from calibration	

	if A, B and C are all set to zero, then standard equation is used for Pt RTD's with Alpha coefficient of 0.00385
	with the following parameters:
						
		T>=0					T<0
	A=	3.908300E-03			3.908300E-03	
	B=	-5.775000E-07			-5.775000E-07	
	C=	N/A						-4.183000E-12

	If C is zero, temperatures greater than 0 deg C are calculated (R >= R0). In this case the temperatures are exact solution to
	quadratic equation.  When C is zero and temperatures are less than 0 deg C ( R < R0), there is no solution and zero deg F is returned.

	When C is non-zero, the temperatures are calculated by trial error.  In this case the maximum error is +/- 0.00032 deg F (+/- 0.00018 deg C)

    R0 must still be provided. Normally R0 is 100, 500 or 1000 ohms

	Function returns temperature in deg F.  If the values are outside the valid range, function returns 0 deg F

*/
/********************************************************************************************/
#ifdef _HPUX_SOURCE
	double rtd_PTgen_OhmToF (double resistanceInOhms, double R0, double A, double B, double C)
#else
	DllExport double __cdecl rtd_PTgen_OhmToF (double resistanceInOhms, double R0, double A, double B, double C)
#endif
{
   double tf;  /* deg F */
   double tc;  /* deg C */
   double rr;
   double AA, BB, CC;
   double xx;
   double xhi, xlo, xguess;
   double ycalc, diff;
   long	  icount;

   tf=0.0;
   if (resistanceInOhms <= 0.0) return 0.0;
   if (R0 <= 0.0) return 0.0;

   rr=resistanceInOhms/R0;

   if(rr < 0.1852007 || rr > 3.90481126) return 0.0;

   if(rr == 1.0)return 32.0;

	AA=A;
	BB=B;
	CC=C;

   if (A == 0.0 && B == 0.0 && C == 0.0)
   {
		AA=  3.908300E-03;
		BB= -5.775000E-07;
		CC= -4.183000E-12;
   }

   if (CC == 0.0)		// only temps greater than 0.0 degC
   {
		if (rr < 1.0)  /* -200 deg C to 0 deg C (-328 to   32 deg F) */
		{
			return 0.0;	// not valid for temps less than 0 deg C
			
		}
		else  /* 0 deg C to 850 deg C (  32 to 1562 deg F) */
		{
			if(BB == 0.0)return 0.0;
			xx=AA*AA - 4.0*BB + 4.0*BB*rr;
			if(xx < 0.0)return 0.0;

			tc= ( -AA + sqrt(xx))/(2.0*BB);
			if(tc < 0.0)return 0.0;
			tf=tc*1.8 + 32.0;  /* deg F */
			return tf;
		}
   }
   else	// all temperatures
   {
		if (rr < 1.0)  // -200 deg C to 0 deg C (-328 to   32 deg F).... Itereative solution
		{
			xhi=850.0;	// deg C
			xlo=-200.0;	// deg C
			icount=0;
			while (icount < 30)
			{
				icount=icount+1;
				xguess=((xhi+xlo)/2.0)*1.8+32.0;	// deg F 
				ycalc=rtd_PTgen_FtoOhm (xguess, R0, AA, BB, CC);	// ohms
				if (ycalc == 0.0)return 0.0;
				diff=ycalc-resistanceInOhms;
				if (fabs(diff/resistanceInOhms) <= 1.0e-6)	// error is less than 0.0003 deg C (< 0.0005 deg F)
				{
						tf=xguess;		// deg F 
						return tf;
				}
				if (diff < 0.0)
				{
					xlo=(xguess-32.0)/1.8;	// deg C
				}
				else
				{
					xhi=(xguess-32.0)/1.8;	// deg C
				}
			}
		 
		}
		else  /* 0 deg C to 850 deg C (  32 to 1562 deg F) */
		{
			if(BB == 0.0)return 0.0;
			xx=AA*AA - 4.0*BB + 4.0*BB*rr;
			if(xx < 0.0)return 0.0;

			tc= ( -AA + sqrt(xx))/(2.0*BB);
			if(tc < 0.0)return 0.0;
			tf=tc*1.8 + 32.0;  /* deg F */
			return tf;
		}
   }
 
   return 0.0;
}

/********************************************************************************************/
/*
	Written by: Nuri Cankurt March 10, 2004

	rtd_PT100A385_OhmToF: This function calculates temperature in deg F for 100 ohm Platinum RTD
	with Alpha coefficient of 0.00385

    The equations are valid for -200 degC (-328 degF) to +850 degC (1562 degF)


	Calling parameters:

	resistanceInOhms:		Resistance value for which temperature is desired

	Function returns temperature in deg F.  If the values are outside the valid range, function returns 0 deg F

*/
/********************************************************************************************/
#ifdef _HPUX_SOURCE
	double rtd_PT100A385_OhmToF (double resistanceInOhms)
#else
	DllExport double __cdecl rtd_PT100A385_OhmToF (double resistanceInOhms)
#endif
{
   double tf;  /* deg F */

   if (resistanceInOhms == 100.0) return 32.0;

//   tf=rtd_PTgenA385_OhmToF2 (resistanceInOhms, 100.0);
	 tf=rtd_PTgen_OhmToF (resistanceInOhms, 100.0, 0.0, 0.0, 0.0);	// Nuri Cankurt 9-29-2008

   return tf;
}

/********************************************************************************************/
/*
	Written by: Nuri Cankurt March 10, 2004

	rtd_PT500A385_OhmToF: This function calculates temperature in deg F for 500 ohm Platinum RTD
	with Alpha coefficient of 0.00385

    The equations are valid for -200 degC (-328 degF) to +850 degC (1562 degF)


	Calling parameters:

	resistanceInOhms:		Resistance value for which temperature is desired

	Function returns temperature in deg F.  If the values are outside the valid range, function returns 0 deg F

*/
/********************************************************************************************/
#ifdef _HPUX_SOURCE
	double rtd_PT500A385_OhmToF (double resistanceInOhms)
#else
	DllExport double __cdecl rtd_PT500A385_OhmToF (double resistanceInOhms)
#endif
{
   double tf;  /* deg F */

   if (resistanceInOhms == 500.0) return 32.0;

//   tf=rtd_PTgenA385_OhmToF2 (resistanceInOhms, 500.0);
   tf=rtd_PTgen_OhmToF (resistanceInOhms, 500.0, 0.0, 0.0, 0.0);	// Nuri Cankurt 9-29-2008

   return tf;
}

/********************************************************************************************/
/*
	Written by: Nuri Cankurt March 10, 2004

	rtd_PT1000A385_OhmToF: This function calculates temperature in deg F for 1000 ohm Platinum RTD
	with Alpha coefficient of 0.00385

    The equations are valid for -200 degC (-328 degF) to +850 degC (1562 degF)


	Calling parameters:

	resistanceInOhms:		Resistance value for which temperature is desired

	Function returns temperature in deg F.  If the values are outside the valid range, function returns 0 deg F

*/
/********************************************************************************************/
#ifdef _HPUX_SOURCE
	double rtd_PT1000A385_OhmToF (double resistanceInOhms)
#else
	DllExport double __cdecl rtd_PT1000A385_OhmToF (double resistanceInOhms)
#endif
{
   double tf;  /* deg F */

   if (resistanceInOhms == 1000.0) return 32.0;

//   tf=rtd_PTgenA385_OhmToF2 (resistanceInOhms, 1000.0);
   tf=rtd_PTgen_OhmToF (resistanceInOhms, 1000.0, 0.0, 0.0, 0.0);	// Nuri Cankurt 9-29-2008

   return tf;
}

/********************************************************************************************/
/*
	Written by: Nuri Cankurt March 10, 2004

	rtd_PT100A385_FtoOhm: This function calculates resistance in ohm for 100 ohm Platinum RTD
	with Alpha coefficient of 0.00385

    The equations are valid for -200 degC (-328 degF) to +850 degC (1562 degF)


	Calling parameters:

	tempInDegF:		temperature in deg F

	Function return resistance in ohm.  If the values are outside the valid range, function returns 0 ohm

*/
/********************************************************************************************/
#ifdef _HPUX_SOURCE
	double rtd_PT100A385_FtoOhm (double tempInDegF)
#else
	DllExport double __cdecl rtd_PT100A385_FtoOhm (double tempInDegF)
#endif
{
   double ohm;  

   if (tempInDegF == 32.0) return 100.0;

//   ohm=rtd_PTgenA385_FtoOhm2 (tempInDegF, 100.0);
   ohm=rtd_PTgen_FtoOhm (tempInDegF, 100.0, 0.0, 0.0, 0.0);	// Nuri Cankurt 9-29-2008

   return ohm;
}

/********************************************************************************************/
/*
	Written by: Nuri Cankurt March 10, 2004

	rtd_PT500A385_FtoOhm: This function calculates resistance in ohm for 500 ohm Platinum RTD
	with Alpha coefficient of 0.00385

    The equations are valid for -200 degC (-328 degF) to +850 degC (1562 degF)


	Calling parameters:

	tempInDegF:		temperature in deg F

	Function return resistance in ohm.  If the values are outside the valid range, function returns 0 ohm

*/
/********************************************************************************************/
#ifdef _HPUX_SOURCE
	double rtd_PT500A385_FtoOhm (double tempInDegF)
#else
	DllExport double __cdecl rtd_PT500A385_FtoOhm (double tempInDegF)
#endif
{
   double ohm;  

   if (tempInDegF == 32.0) return 500.0;

//   ohm=rtd_PTgenA385_FtoOhm2 (tempInDegF, 500.0);
   ohm=rtd_PTgen_FtoOhm (tempInDegF, 500.0, 0.0, 0.0, 0.0);	// Nuri Cankurt 9-29-2008

   return ohm;
}

/********************************************************************************************/
/*
	Written by: Nuri Cankurt March 10, 2004

	rtd_PT1000A385_FtoOhm: This function calculates resistance in ohm for 1000 ohm Platinum RTD
	with Alpha coefficient of 0.00385

    The equations are valid for -200 degC (-328 degF) to +850 degC (1562 degF)


	Calling parameters:

	tempInDegF:		temperature in deg F

	Function return resistance in ohm.  If the values are outside the valid range, function returns 0 ohm

*/
/********************************************************************************************/
#ifdef _HPUX_SOURCE
	double rtd_PT1000A385_FtoOhm (double tempInDegF)
#else
	DllExport double __cdecl rtd_PT1000A385_FtoOhm (double tempInDegF)
#endif
{
   double ohm;  /* deg F */

   if (tempInDegF == 32.0) return 1000.0;

//   ohm=rtd_PTgenA385_FtoOhm2 (tempInDegF, 1000.0);
   ohm=rtd_PTgen_FtoOhm (tempInDegF, 1000.0, 0.0, 0.0, 0.0);	// Nuri Cankurt 9-29-2008

   return ohm;
}

/********************************************************************************/
/*																				*/
/* polyValue: This function returns result of a polynomial equation				*/
/*																				*/
/*	xval: X value for which polynomial function is calculated					*/
/*  coef: Coeefients for polynomial function 0th, 1st, 2nd,.....Nth	order		*/
/*  order:  Polynomial function order (i.e. 9 for 9th order polynomial			*/
/*	if (order=0)polyValue=coef[0]												*/
/*  if (order>0)polyValue=coef[0]+coef[1]*xval+coef[2]*xval^2 +..+coef[N]*xval^N*/
/*  if (order<0)polyValue=coef[0]+coef[1]/xval+coef[2]/xval^2 +..+coef[N]/xval^N*/
/********************************************************************************/
#ifdef _HPUX_SOURCE
	double polyValue (double *coef, double xval, long order)
#else
	DllExport double __cdecl polyValue (double *coef, double xval, long order)
#endif
{
	long	orderUsed;
	long	ii;
	double	result;	
	double	lastFactor;
	double	factor;

	orderUsed=order;

	if (order == 0)
	{
		return coef[0];
	}
	else
	if (order > 0)
	{
		if (xval == 0.0)return coef[0];
	}
	else
	if (order < 0)
	{
		if (xval == 0.0)return 0.0; /* error */
		orderUsed=(-1)*order;
	}

	result=coef[0];
	lastFactor=1.0;

	if (order > 0)
	{
		for (ii=1;ii<=orderUsed;ii++)
		{
			factor=lastFactor*xval;
			result=result + coef[ii]*factor;
			lastFactor=factor;
		}
	}
	else
	if (order < 0)
	{
		for (ii=1;ii<=orderUsed;ii++)
		{
			factor=lastFactor*xval;
			result=result + coef[ii]/factor;
			lastFactor=factor;
		}
	}

	return result;
}

/********************************************************************************/
/*																				*/
/* linearLookup: This function interpolates set of xdata and ydata values		*/
/*	linearly and returns yval based on xval.  Data need not be in sorted order	*/
/*																				*/
/*	xdata and ydata: Array of data sets											*/
/*  xval: X value at which y value will be calculated							*/
/*	yval: Y value calculated from xval											*/
/*	numData: Number data points in xdata and ydata								*/
/*																				*/
/*  Function returns 0 if no error, or the following error codes:				*/
/*	-1 if numData < 2  (not enough data points)									*/
/*	-2 if id xval is outside of the scope on the low end						*/
/*	-3 if id xval is outside of the scope on the high end						*/
/* `-4 if error in sort function												*/
/*	-5 Other error																*/
/********************************************************************************/
#ifdef _HPUX_SOURCE
	long linearLookup (double *xdata, double *ydata, double xval, double *yval, long numData)
#else
	DllExport long __cdecl linearLookup (double *xdata, double *ydata, double xval, double *yval, long numData)
#endif
{
	long	ii;
	long	error;
	long	*sortedIndex;
	double	x1,x2,y1,y2;

	*yval=0.0;

	if (numData < 2)return -1;

	sortedIndex= (long *)malloc (numData*sizeof(long));

	error= sortDoublePtr (xdata, numData, 'a', sortedIndex);
	if (error !=0)
	{
		free(sortedIndex);
		return -4;
	}

	if (xval < xdata[sortedIndex[0]])
	{
		free(sortedIndex);
		return -2;
	}
	if (xval > xdata[sortedIndex[numData-1]])
	{
		free(sortedIndex);
		return -3;
	}

	for (ii=0;ii<numData;ii++)
	{
		if (xval <= xdata[sortedIndex[ii]])
		{
			if (ii == 0)
			{
				free(sortedIndex);
				return -5;
			}

			x1=xdata[sortedIndex[ii-1]];
			y1=ydata[sortedIndex[ii-1]];

			x2=xdata[sortedIndex[ii]];
			y2=ydata[sortedIndex[ii]];

			if (fabs(x2-x1) < 1.0e-40)
			{
				free(sortedIndex);
				return -5;
			}

			*yval =((xval-x1)/(x2-x1))*(y2-y1) + y1;
			free(sortedIndex);
			return 0;
		}
	}

	free(sortedIndex);

	return -5;
}