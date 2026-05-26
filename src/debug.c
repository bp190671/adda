/* Functions for printing debugging information when compiling with option -DDEBUGFULL.
 * Should only be used with this macro, but IDE syntax checking works fine even without it.
 *
 * Copyright (C) ADDA contributors
 * This file is part of ADDA.
 *
 * ADDA is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * ADDA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with ADDA. If not, see
 * <http://www.gnu.org/licenses/>.
 */
#include "const.h" // keep this first
#include "debug.h" // corresponding header
// project headers
#include "cmplx.h"
#include "comm.h"
#include "io.h"
#include "types.h"
#include "vars.h"
// system headers
#include <stdarg.h>
#include <stdio.h>
#include <string.h>

//======================================================================================================================

void DebugPrintf(ERR_LOC_DECL,const char * restrict fmt, ... )
/* Prints anything to stdout with additional debug information
 *
 * Not thread-safe! Should not be called in parallel from multiple threads (e.g. OpenMP)
 */
{
	va_list args;
	static char msg[MAX_PARAGRAPH]; // not to allocate it at every call

	if (who==ALL || IFROOT) { // controls whether output should be produced
		va_start(args,fmt);
		VsnprintfErr(ERR_LOC_CALL,msg,MAX_PARAGRAPH,fmt,args);
#ifdef PARALLEL
		if (who==ALL) printf("(ringID=%i) DEBUG: %s:%d: %s \n",ringid,srcfile,srcline,msg);
		else
#endif
		printf("DEBUG: %s:%d: %s \n",srcfile,srcline,msg);
		fflush(stdout);
		va_end(args);
	}
}

//======================================================================================================================

void FieldPrint(doublecomplex * restrict x)
/* print current field at certain dipole -- not used; left for deep debug; NOT ROBUST, since DipoleCoord is not always
 * available
 */
{
	int i=9810;

	i*=3;
	fprintf(logfile,"Dipole coordinates = "GFORM3V"\n",COMP3V(DipoleCoord+i));
	fprintf(logfile,"E = "CFORM3V,REIM3V(x+i));
}

//======================================================================================================================

void PrintScalar(const enum incpol which,doublecomplex * restrict cmplxF,
	double * restrict realF,const char * restrict fname_preffix,const char * restrict tmpl UOIP,
	const char * restrict field_name,const char * restrict fullname)
{
	FILE * restrict file;
	size_t j;
	char fname[MAX_FNAME],fname_sh[MAX_FNAME_SH];
	bool cmplx_mode; 
	if ((cmplxF==NULL) ^ (realF==NULL)) cmplx_mode=(realF==NULL);
	else LogError(ONE_POS,"One field (either real or complex) must be given to StoreFields");
	strcpy(fname_sh,fname_preffix);
	if (which==INCPOL_Y) strcat(fname_sh,F_YSUF);
	else strcat(fname_sh,F_XSUF);
#ifdef PARALLEL
	size_t shift=SnprintfErr(ALL_POS,fname,MAX_FNAME,"%s/",directory);
	SnprintfShiftErr(ALL_POS,shift,fname,MAX_FNAME,tmpl,ringid);
#else
	SnprintfErr(ALL_POS,fname,MAX_FNAME,"%s/%s",directory,fname_sh);
#endif
	file=FOpenErr(fname,"w",ALL_POS);
#ifdef PARALLEL
	if (ringid==0) {
#endif
		if (cmplx_mode) fprintf(file,"x y z Re(%s) Im(%s)\n",
			field_name,field_name);
		else fprintf(file,"x y z %s\n",field_name);
#ifdef PARALLEL
	}
#endif
	if (cmplx_mode) for (j=0;j<local_nvoid_Ndip;j++) fprintf(file,GFORM5L"\n",COMP3V(DipoleCoord+3*j),
		REIM(cmplxF[j]));
	else for (j=0;j<local_nvoid_Ndip;j++) fprintf(file,GFORM4L"\n",COMP3V(DipoleCoord+3*j),realF[j]);
	FCloseErr(file,fname,ALL_POS);
#ifdef PARALLEL
	Synchronize();
	if (IFROOT) CatNFiles(directory,tmpl,fname_sh);
#endif
	if (IFROOT) PRINTFB("%s saved to file\n",fullname);
}

//======================================================================================================================

void PrintVector(const enum incpol which,doublecomplex * restrict cmplxF,
	double * restrict realF,const char * restrict fname_preffix,const char * restrict tmpl UOIP,
	const char * restrict field_name,const char * restrict fullname)
{
	FILE * restrict file;
	size_t j;
	char fname[MAX_FNAME],fname_sh[MAX_FNAME_SH];
	bool cmplx_mode; 
	if ((cmplxF==NULL) ^ (realF==NULL)) cmplx_mode=(realF==NULL);
	else LogError(ONE_POS,"One field (either real or complex) must be given to StoreFields");
	strcpy(fname_sh,fname_preffix);
	if (which==INCPOL_Y) strcat(fname_sh,F_YSUF);
	else strcat(fname_sh,F_XSUF);
#ifdef PARALLEL
	size_t shift=SnprintfErr(ALL_POS,fname,MAX_FNAME,"%s/",directory);
	SnprintfShiftErr(ALL_POS,shift,fname,MAX_FNAME,tmpl,ringid);
#else
	SnprintfErr(ALL_POS,fname,MAX_FNAME,"%s/%s",directory,fname_sh);
#endif
	file=FOpenErr(fname,"w",ALL_POS);
#ifdef PARALLEL
	if (ringid==0) {
#endif
		if (cmplx_mode) fprintf(file,"x y z Re(%s.x) Im(%s.x) Re(%s.y) Im(%s.y) Re(%s.z) Im(%s.z)\n",
			field_name,field_name,field_name,field_name,field_name,field_name);
		else fprintf(file,"x y z %s.x %s.y %s.z\n",field_name,field_name,field_name);
#ifdef PARALLEL
	}
#endif
	if (cmplx_mode) for (j=0;j<local_nvoid_Ndip;j++) fprintf(file,GFORM9L"\n",COMP3V(DipoleCoord+3*j),
		REIM3V(cmplxF+3*j));
	else for (j=0;j<local_nvoid_Ndip;j++) fprintf(file,GFORM6L"\n",COMP3V(DipoleCoord+3*j),
		COMP3V(realF+3*j));
	FCloseErr(file,fname,ALL_POS);
#ifdef PARALLEL
	Synchronize();
	if (IFROOT) CatNFiles(directory,tmpl,fname_sh);
#endif
	if (IFROOT) PRINTFB("%s saved to file\n",fullname);
}

//======================================================================================================================

void PrintTensor(const enum incpol which,doublecomplex * restrict cmplxF,
	double * restrict realF,const char * restrict fname_preffix,const char * restrict tmpl UOIP,
	const char * restrict field_name,const char * restrict fullname)
{
	FILE * restrict file;
	size_t j;
	char fname[MAX_FNAME],fname_sh[MAX_FNAME_SH];
	bool cmplx_mode; 
	if ((cmplxF==NULL) ^ (realF==NULL)) cmplx_mode=(realF==NULL);
	else LogError(ONE_POS,"One field (either real or complex) must be given to StoreFields");
	strcpy(fname_sh,fname_preffix);
	if (which==INCPOL_Y) strcat(fname_sh,F_YSUF);
	else strcat(fname_sh,F_XSUF);
#ifdef PARALLEL
	size_t shift=SnprintfErr(ALL_POS,fname,MAX_FNAME,"%s/",directory);
	SnprintfShiftErr(ALL_POS,shift,fname,MAX_FNAME,tmpl,ringid);
#else
	SnprintfErr(ALL_POS,fname,MAX_FNAME,"%s/%s",directory,fname_sh);
#endif
	file=FOpenErr(fname,"w",ALL_POS);
#ifdef PARALLEL
	if (ringid==0) {
#endif
		if (cmplx_mode) fprintf(file,"x y z Re(%s.xx) Im(%s.xx) Re(%s.xy) Im(%s.xy) Re(%s.xz) Im(%s.xz) Re(%s.yx) Im(%s.yx) Re(%s.yy) Im(%s.yy) Re(%s.yz) Im(%s.yz) Re(%s.zx) Im(%s.zx) Re(%s.zy) Im(%s.zy) Re(%s.zz) Im(%s.zz)\n",
			field_name,field_name,field_name,field_name,field_name,field_name,field_name,field_name,field_name,field_name,field_name,field_name,field_name,field_name,field_name,field_name,field_name,field_name);
		else fprintf(file,"x y z %s.xx %s.xy %s.xz %s.yx %s.yy %s.yz %s.zx %s.zy %s.zz\n",field_name,field_name,field_name,field_name,field_name,field_name,field_name,field_name,field_name);
#ifdef PARALLEL
	}
#endif
	if (cmplx_mode) for (j=0;j<local_nvoid_Ndip;j++) fprintf(file,GFORM21L"\n",COMP3V(DipoleCoord+3*j),
		REIM9V(cmplxF+9*j));
	else for (j=0;j<local_nvoid_Ndip;j++) fprintf(file,GFORM12L"\n",COMP3V(DipoleCoord+3*j),
		COMP9V(realF+9*j));
	FCloseErr(file,fname,ALL_POS);
#ifdef PARALLEL
	Synchronize();
	if (IFROOT) CatNFiles(directory,tmpl,fname_sh);
#endif
	if (IFROOT) PRINTFB("%s saved to file\n",fullname);
}

