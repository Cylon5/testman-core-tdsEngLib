/* 
	This is tdsEngLib.dll external references file.  This DLL must exist in the same path
	as the calling program executable or in a directory within windows search path.
	(ie. c:\winnt\system32)

	File name	: 'tdsEngLib.h'
	Version		: 1.0 original version
	Written by	: Nuri Cankurt
	Date		: March 13, 2004
*/

#ifdef _HPUX_SOURCE
	extern long  linearLSCF (double *xdata, double *ydata, long num, double *slope, double *intercept, double *corrCoef);
	extern long  sortDoublePtr (double *fdata, long num, char how, long *index);

	extern void  utcToLocalDateTimeStr (long UTCsec, char *dateStr, char *timeStr);
	extern void  utcToLocalDateTimeStr2 (double fractionalUTCsec, char *dateStr, char *timeStr, long *millisec);

	extern double  ctimeToDbl (char *timeString);
	extern long  timeToDDDHHMMSS (double timeVal);
	extern double  timeToDDDHHMMSS_MS (double timeVal);
	extern long  timeToDay (double timeVal);
	extern long  timeToHour (double timeVal);
	extern long  timeToMinute (double timeVal);
	extern long  timeToSecond (double timeVal);
	extern long  timeToMillisec (double timeVal);
	extern long  timeToHHMMSS (double timeVal);
	extern double  timeToHHMMSS_MS (double timeVal);
	extern double  timeToDecDay (double timeVal);
	extern double  timeToYearSecond (double timeVal);
	extern long  timeToStdStr (double timeVal, char *strout);
	extern long  timeToCtime (double timeVal, char *strout);

	extern double  tck_V32toF (double volts);
	extern double  tck_FtoV32 (double temp);
	extern double  tck_VtoF (double volts, double refTemp);

	extern double  tcb_V32toF (double volts);
	extern double  tcb_FtoV32 (double temp);
	extern double  tcb_VtoF (double volts, double refTemp);

	extern double  tce_V32toF (double volts);
	extern double  tce_FtoV32 (double temp);
	extern double  tce_VtoF (double volts, double refTemp);

	extern double  tcj_V32toF (double volts);
	extern double  tcj_FtoV32 (double temp);
	extern double  tcj_VtoF (double volts, double refTemp);

	extern double  tcn_V32toF (double volts);
	extern double  tcn_FtoV32 (double temp);
	extern double  tcn_VtoF (double volts, double refTemp);

	extern double  tcr_V32toF (double volts);
	extern double  tcr_FtoV32 (double temp);
	extern double  tcr_VtoF (double volts, double refTemp);

	extern double  tcs_V32toF (double volts);
	extern double  tcs_FtoV32 (double temp);
	extern double  tcs_VtoF (double volts, double refTemp);

	extern double  tct_V32toF (double volts);
	extern double  tct_FtoV32 (double temp);
	extern double  tct_VtoF (double volts, double refTemp);

	extern double  rtd_PTgenA385_OhmToC2 (double resistanceInOhms, double refResistanceInOhms);
	extern double  rtd_PTgen_OhmToC (double resistanceInOhms, double R0, double A, double B, double C);
	extern double  rtd_PT100A385_OhmToC (double resistanceInOhms);
	extern double  rtd_PT500A385_OhmToC (double resistanceInOhms);
	extern double  rtd_PT1000A385_OhmToC (double resistanceInOhms);

	extern double  rtd_PTgenA385_CtoOhm2 (double tempInDegC, double refResistanceInOhms);
	extern double  rtd_PTgen_CtoOhm (double tempInDegC, double R0, double A, double B, double C);
	extern double  rtd_PT100A385_CtoOhm (double tempInDegC);
	extern double  rtd_PT500A385_CtoOhm (double tempInDegC);
	extern double  rtd_PT1000A385_CtoOhm (double tempInDegC);

	extern double  rtd_PTgenA385_OhmToF2 (double resistanceInOhms, double refResistanceInOhms);
	extern double  rtd_PTgen_OhmToF (double resistanceInOhms, double R0, double A, double B, double C);
	extern double  rtd_PT100A385_OhmToF (double resistanceInOhms);
	extern double  rtd_PT500A385_OhmToF (double resistanceInOhms);
	extern double  rtd_PT1000A385_OhmToF (double resistanceInOhms);

	extern double  rtd_PTgenA385_FtoOhm2 (double tempInDegF, double refResistanceInOhms);
	extern double  rtd_PTgen_FtoOhm (double tempInDegF, double R0, double A, double B, double C);
	extern double  rtd_PT100A385_FtoOhm (double tempInDegF);
	extern double  rtd_PT500A385_FtoOhm (double tempInDegF);
	extern double  rtd_PT1000A385_FtoOhm (double tempInDegF);

	double polyValue (double *coef, double xval, long order);
	long linearLookup (double *xdata, double *ydata, double xval, double *yval, long numData);
#else
	extern long __cdecl linearLSCF (double *xdata, double *ydata, long num, double *slope, double *intercept, double *corrCoef);
	extern long __cdecl sortDoublePtr (double *fdata, long num, char how, long *index);

	extern void __cdecl utcToLocalDateTimeStr (long UTCsec, char *dateStr, char *timeStr);
	extern void __cdecl utcToLocalDateTimeStr2 (double fractionalUTCsec, char *dateStr, char *timeStr, long *millisec);

	extern double __cdecl ctimeToDbl (char *timeString);
	extern long __cdecl timeToDDDHHMMSS (double timeVal);
	extern double __cdecl timeToDDDHHMMSS_MS (double timeVal);
	extern long __cdecl timeToDay (double timeVal);
	extern long __cdecl timeToHour (double timeVal);
	extern long __cdecl timeToMinute (double timeVal);
	extern long __cdecl timeToSecond (double timeVal);
	extern long __cdecl timeToMillisec (double timeVal);
	extern long __cdecl timeToHHMMSS (double timeVal);
	extern double __cdecl timeToHHMMSS_MS (double timeVal);
	extern double __cdecl timeToDecDay (double timeVal);
	extern double __cdecl timeToYearSecond (double timeVal);
	extern long __cdecl timeToStdStr (double timeVal, char *strout);
	extern long __cdecl timeToCtime (double timeVal, char *strout);

	extern double __cdecl tck_V32toF (double volts);
	extern double __cdecl tck_FtoV32 (double temp);
	extern double __cdecl tck_VtoF (double volts, double refTemp);

	extern double __cdecl tcb_V32toF (double volts);
	extern double __cdecl tcb_FtoV32 (double temp);
	extern double __cdecl tcb_VtoF (double volts, double refTemp);

	extern double __cdecl tce_V32toF (double volts);
	extern double __cdecl tce_FtoV32 (double temp);
	extern double __cdecl tce_VtoF (double volts, double refTemp);

	extern double __cdecl tcj_V32toF (double volts);
	extern double __cdecl tcj_FtoV32 (double temp);
	extern double __cdecl tcj_VtoF (double volts, double refTemp);

	extern double __cdecl tcn_V32toF (double volts);
	extern double __cdecl tcn_FtoV32 (double temp);
	extern double __cdecl tcn_VtoF (double volts, double refTemp);

	extern double __cdecl tcr_V32toF (double volts);
	extern double __cdecl tcr_FtoV32 (double temp);
	extern double __cdecl tcr_VtoF (double volts, double refTemp);

	extern double __cdecl tcs_V32toF (double volts);
	extern double __cdecl tcs_FtoV32 (double temp);
	extern double __cdecl tcs_VtoF (double volts, double refTemp);

	extern double __cdecl tct_V32toF (double volts);
	extern double __cdecl tct_FtoV32 (double temp);
	extern double __cdecl tct_VtoF (double volts, double refTemp);

	extern double __cdecl rtd_PTgenA385_OhmToC2 (double resistanceInOhms, double refResistanceInOhms);
	extern double __cdecl rtd_PTgen_OhmToC (double resistanceInOhms, double R0, double A, double B, double C);
	extern double __cdecl rtd_PT100A385_OhmToC (double resistanceInOhms);
	extern double __cdecl rtd_PT500A385_OhmToC (double resistanceInOhms);
	extern double __cdecl rtd_PT1000A385_OhmToC (double resistanceInOhms);

	extern double __cdecl rtd_PTgenA385_CtoOhm2 (double tempInDegC, double refResistanceInOhms);
	extern double __cdecl rtd_PTgen_CtoOhm (double tempInDegC, double R0, double A, double B, double C);
	extern double __cdecl rtd_PT100A385_CtoOhm (double tempInDegC);
	extern double __cdecl rtd_PT500A385_CtoOhm (double tempInDegC);
	extern double __cdecl rtd_PT1000A385_CtoOhm (double tempInDegC);

	extern double __cdecl rtd_PTgenA385_OhmToF2 (double resistanceInOhms, double refResistanceInOhms);
	extern double __cdecl rtd_PTgen_OhmToF (double resistanceInOhms, double R0, double A, double B, double C);
	extern double __cdecl rtd_PT100A385_OhmToF (double resistanceInOhms);
	extern double __cdecl rtd_PT500A385_OhmToF (double resistanceInOhms);
	extern double __cdecl rtd_PT1000A385_OhmToF (double resistanceInOhms);

	extern double __cdecl rtd_PTgenA385_FtoOhm2 (double tempInDegF, double refResistanceInOhms);
	extern double __cdecl rtd_PTgen_FtoOhm (double tempInDegF, double R0, double A, double B, double C);
	extern double __cdecl rtd_PT100A385_FtoOhm (double tempInDegF);
	extern double __cdecl rtd_PT500A385_FtoOhm (double tempInDegF);
	extern double __cdecl rtd_PT1000A385_FtoOhm (double tempInDegF);

	extern double __cdecl polyValue (double *coef, double xval, long order);
	extern long __cdecl linearLookup (double *xdata, double *ydata, double xval, double *yval, long numData);
#endif
