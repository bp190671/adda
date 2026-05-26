/* All the initialization is done here before actually calculating internal fields,
 * includes calculation of couple constants
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
#define _USE_MATH_DEFINES
#include "const.h" // keep this first
// project headers
#include "cmplx.h"
#include "comm.h"
#include "crosssec.h"
#include "debug.h"
#include "fft.h"
#include "interaction.h"
#include "io.h"
#include "memory.h"
#include "oclcore.h"
#include "Romberg.h"
#include "timing.h"
#include "vars.h"
#include "WD.h"

// system headers
#include <math.h>
#include <stdlib.h>
#include <string.h>

//======================================================================================================================
/*Each voxel is defined as the unit cube [-0.5,0.5]x[-0.5,0.5]x[-0.5,0.5] with the origin in the center*/
const double cubeOrderedPoints[8][3] = {{-0.5,-0.5,-0.5},{0.5,-0.5,-0.5},{-0.5,0.5,-0.5},{-0.5,-0.5,0.5},{0.5,-0.5,0.5},{0.5,0.5,-0.5},{-0.5,0.5,0.5},{0.5,0.5,0.5}};
const int cubeOrderedEdges[6][4] = {{0,2,5,1}, {0,3,6,2}, {3,4,7,6}, {4,1,5,7}, {0,1,4,3}, {2,6,7,5}};
const double CubeCenter[3] = {0,0,0};
const double CubeEdgeNorm[6][3] = {{0,0,-1},{-1,0,0},{0,0,1},{1,0,0},{0,-1,0},{0,1,0}};

//======================================================================================================================

bool vEq(double v[3], double u[3]) {
	/*Used to compare two vectors which can be considered as equivalent*/
	double diff[3];
	vSubtr(v, u, diff);
	return DotProd(diff, diff) < __DBL_EPSILON__;
}

//==========================================================

void BubbleSort(int* indxs, double* values, int n) {
	/*Bubble sorting of the values of an array*/
	int boof;
	for (int i = 0; i<n; i++){
		for (int j=i+1; j<n; j++){
			double vi = values[indxs[i]], vj = values[indxs[j]];
			if (vi<vj) {
				boof=indxs[i];
				indxs[i]=indxs[j];
				indxs[j]=boof;
			}
		}
	}
}

//======================================================================================================================

bool PointUpperPlane(const double p[static 3], const double n[static 3]) {
	/*Checks whether the point belongs to the secondary domain or not : this condition is never verified when the point
	corresponds to the voxel center.*/
	return DotProd(p, n) > 1.0;
}

//======================================================================================================================

void IntersectionParams(double nv[8], double n[3]) {
	/*Helper for storing the intersections of vertices with the plane in an array*/
	for (int i = 0; i < 8; i++) nv[i] = DotProd(cubeOrderedPoints[i], n);
}

//======================================================================================================================

int CubePlaneIntersect(const double n[3], double intersections[24][3], int edgeCounts[6]) {
	/*Determines the number of intersections per edge within the cube*/
	int i_idx = 0;
	double nv[8];

	IntersectionParams(nv, n);
	for (int nEdge = 0; nEdge < 6; nEdge++) {
		int edgeIntCount = 0;
		for (int j = 0; j < 4; j++) {
			int k = (j+1) % 4;
			int v1 = cubeOrderedEdges[nEdge][j];
			int v2 = cubeOrderedEdges[nEdge][k];
			if (EdgeIn(v1,v2,nv,intersections[i_idx])) {
				i_idx++;
				edgeIntCount++;
			}

		}

		edgeCounts[nEdge] = edgeIntCount;
	}

	return i_idx;
}

//======================================================================================================================

void PointsCenter(const double p[24][3], int n, double c[3]) {
	/*Computes the centroid from a set of vertices*/
	vInit(c);
	for (int i = 0; i < n; i++)
		vAdd(c, p[i], c);

	vMultScal(1.0 / (double) n, c, c);
}

//======================================================================================================================

double PointWalkOrd(const double u0[3], const double u[3], const double n[3]) {
	/*Ordering function Ord() for points around a normal vector : it is based on some reference point*/
	double vecProd[3];

	if (vEq(u0, u))
		return 2;

	CrossProd(u0, u, vecProd);
	double s = DotProd(vecProd, n);
	if (s > 0) s = 1;
	if (s < 0) s = -1;

	return s * (1.0 + AngleCos(u0, u));
}

//======================================================================================================================

int ReorderPoints(const double p[24][3], int k, const double n[3], double pReordered[24][3]) {
	/*Reorders points around their centroid, sorts them based on their values of Ord() and removes duplicate points*/
	double c[3];
	for(int l=0;l<3;l++) c[l]=CubeCenter[l];
	double u[24][3];
	double ords[24];
	int walkIndexes[24];
	if (k == 0) return 0;

	// switch to figure center coords
	PointsCenter(p, k, c);
	for (int i = 0; i < k; i++)
		vSubtr(p[i], c, u[i]);

	// calc ords for all vectors
	ords[0] = 2.0; walkIndexes[0] = 0;
	for (int i = 1; i < k; i++) {
		ords[i] = PointWalkOrd(u[0], u[i], n);
		walkIndexes[i] = i;
	}

	// reorder based on ords[], descending
	BubbleSort(walkIndexes, ords, k);

	// remove duplicate points and place to pReordered
	int prevApproved = 0, countPlaced = 1;
	vCopy(u[walkIndexes[0]], pReordered[0]);

	for (int i = 1; i < k; i++) {
		if (!vEq(u[walkIndexes[i]], u[prevApproved])){
			vCopy(u[walkIndexes[i]], pReordered[countPlaced]);
			prevApproved = walkIndexes[i];
			countPlaced++;
		}
	}

	// switch back to 0,0,0 coords
	for (int i = 0; i < countPlaced; i++)
		vAdd(pReordered[i], c, pReordered[i]);

	return countPlaced;
}

//======================================================================================================================
//added on 5.01.22.
bool EdgeIn(int i,int j,double nv[8],double res[3])
/* Finds intersection of plane n.r=1 with edges of the unit cube [-0.5,0.5]x[-0.5,0.5]x[-0.5,0.5]. Coordinates of 
 * cube vertices are given by static array v below. n is provided by its scalar products with v - by vector nv.
 * i and j are indices of adjacent vertices, defining the edge. Returns true if intersection is
 * inside the edge. Then (only if true) coordinates of the intersection are stored in res.
 * This code is based on https://cococubed.com/code_pages/raybox.shtml by F.X.Timmes
 */
{
	static const double v[8][3]={{-0.5,-0.5,-0.5},{0.5,-0.5,-0.5},{-0.5,0.5,-0.5},{-0.5,-0.5,0.5},{0.5,-0.5,0.5},
	{0.5,0.5,-0.5},{-0.5,0.5,0.5},{0.5,0.5,0.5}};
	bool cond;
	double t;

	if (nv[i]==nv[j]) cond=false;
	else {
		t=(1-nv[i])/(nv[j]-nv[i]);
		cond=(t>=0 && t<=1);
		if (cond) LinComb(v[i],v[j],1-t,t,res);
	}
	return cond;
}

//==========================================================

double VolumeFraction(double a,double b,double c)
/*Computes the volume fraction of the principal part (p) for the [-0.5,0.5]x[-0.5,0.5]x[-0.5,0.5] cube. A ready-to-use 
formula already exists for the evaluation of volume fraction based on the [0,1]x[0,1]x[0,1] cube and its intersection
vertices given by double nv[8]={0,a,b,c,a+c,a+b,b+c,a+b+c}. It requires some transformation of plane coefficients a,b,c 
to ensure that volume fraction is consistent with the change of origin.
The general formula (works for a,b,c non-negative) is f0=[1-h(a)-h(b)-h(c)+h(a+b)+h(a+c)+h(b+c)-h(a+b+c)]/(6abc), where
 h(r)=(1-r)^3, 0<r<1 and 0 otherwise.
The function optimizes the formula for speed and takes care of coefficients close to zero through some sorting algorithm
for plane coefficients and explicit expressions for the pyramidal decomposition.*/
{
	double A = fabs(a); double B = fabs(b); double C = fabs(c);
	// Transformation coefficient :
	double nAbs[3]={A,B,C}; // Modified normal vector (|a|,|b|,|c|)
	double shift[3]={-0.5,-0.5,-0.5};
	double gamma = 1.0-DotProd(nAbs,shift);
	if (A+B+C <= gamma) return 1.0;
	if (gamma>1.0/__DBL_EPSILON__) gamma=1.0;
	A *= 1.0/gamma;
	B *= 1.0/gamma;
	C *= 1.0/gamma;
	
	double tmp;
	// sort A,B,C in ascending order (A<=B<=C)
	if (A>B) {
		if (A>C) {
			tmp=C;
			C=A;
			if (B>tmp) A=tmp;
			else {
				A=B;
				B=tmp;
			}
		}
		else {
			tmp=A;
			A=B;
			B=tmp;
		}
	}
	else if (B>C) {
		tmp=C;
		C=B;
		if (A>tmp) {
			B=A;
			A=tmp;
		}
		else B=tmp;
	}
	// Pyramidal decomposition takes care of A,B,C->0 for f0
	// Tetrahedron of volume V0=1/6[((1/A,0,0)x(0,1/B,0)).(0,0,1/C)] formed between origin and plane
	if (A>=1) tmp=1/(2*A*B); // A>=1,B>=1,C>=1 : f0=1/(6abc)
	// below A<1
	else if (B>=1) tmp=(3-3*A+A*A)/(2*B); // A<1,B>=1,C>=1 : f0=(1-h(A))/(6abc)
	else { // below B<1
		if (A+B>=1) { // A<1,B<1,C>=1,A+B>=1 : f0=(1-h(A)-h(B))/(6abc)
			tmp=A*(3-3*A+A*A)-(1-B)*(1-B)*(1-B);
			if (C<1) tmp-=(1-C)*(1-C)*(1-C); // A<1,B<1,C<1,A+B>=1 : f0=(1-h(A)-h(B)-h(C))/(6abc)
			tmp/=2*A*B;
		}
		// below A+B<1
		else if (C>=1) tmp=3-1.5*(A+B); // A<1,B<1,C>=1,A+B<1 : f0=(1-h(A)-h(B)+h(A+B))/(6abc)
		// below C<1
		else if (A+C>=1) tmp=3-1.5*(A+B)-(1-C)*(1-C)*(1-C)/(2*A*B); // A<1,B<1,C<1,A+B<1,A+C>=1 : f0=(1-h(A)-h(B)-h(C)+h(A+B))/(6abc)
		// below A+C<1
		else if (B+C>=1) { // A<1,B<1,C<1,A+B<1,A+C<1,B+C>=1 : f0=(1-h(A)-h(B)-h(C)+h(A+B)+h(A+C))/(6abc)
			tmp=3-1.5*(A+B)+C/B*(3-1.5*(A+C))-(3-3*A+A*A)/(2*B);
		}
		else { // A<1,B<1,C<1,A+B<1,A+C<1,B+C<1,A+B+C>=1 : f0=(1-h(A)-h(B)-h(C)+h(A+B)+h(A+C)+h(B+C))/(6abc)
			tmp=3-1.5*(A+B)+C/B*(3-1.5*(A+C))-(3-3*A+A*A)/(2*B)+(1-B-C)*(1-B-C)*(1-B-C)/(2*A*B);
		}
	} // A<1,B<1,C<1,A+B<1,A+C<1,B+C<1,A+B+C<1 : f0=(1-h(A)-h(B)-h(C)+h(A+B)+h(A+C)+h(B+C)-h(A+B+C))/(6abc)
	if (C==0) return 0.0;
	double vf=tmp/(3*C); //Here, C=0 is always avoided because of sorting A<=B<=C, leading to A=B=C=0 (no plane)
	extern const bool curvcor;
	if(curvcor){
	/*The following term is a correction to the plane approximation which accounts the curvature effects. These effects
	cannot be taken into account into scattering quantities since a plane interface is needed to warrant constant fields 
	in each subvoxel: in this respect, this correction is used only to decrease errors related to the scatterer volume*/
		extern enum sh shape;
		double n[3]={a,b,c};
		double nv[8];
		double p[6][3];
		double Rc; //Radius of curvature
		IntersectionParams(nv, n);
		int vN=0;
		if (EdgeIn(0,1,nv,p[vN]) || EdgeIn(1,4,nv,p[vN]) || EdgeIn(4,7,nv,p[vN])) vN++; // 0->1->4->7
		if (EdgeIn(1,5,nv,p[vN])) vN++; // 1->5
		if (EdgeIn(0,2,nv,p[vN]) || EdgeIn(2,5,nv,p[vN]) || EdgeIn(5,7,nv,p[vN])) vN++; // 0->2->5->7
		if (EdgeIn(2,6,nv,p[vN])) vN++; // 2->6
		if (EdgeIn(0,3,nv,p[vN]) || EdgeIn(3,6,nv,p[vN]) || EdgeIn(6,7,nv,p[vN])) vN++; // 0->3->6->7
		if (EdgeIn(3,4,nv,p[vN])) vN++; // 3->4
		// now vN is the number of vertices
		/* compute moment of inertia of the polygon (<r^2>*S) respective to the axis defined by n
		* formula is based on http://en.wikipedia.org/wiki/List_of_moments_of_inertia with correction,
		* so that the sum includes i=N as well
		*/
		double v[3]; //Intersection point as n.v=1
		vMultScal(1/DotProd(n,n),n,v); // |v|=1/|n| -> v=n/|n|²
		int i,j;
		for (i=0;i<vN;i++) for(j=0;j<3;j++) p[i][j]-=v[j]; // shift vertices vectors into the plane
		double iner=0;
		for (i=0,j=1;i<vN;i++,j++) { // main sum for inertia moment
			if (j==vN) j=0; // cycle for the last vertex
			iner+=AbsOutProd(p[i],p[j])*(DotProd(p[i],p[i])+DotProd(p[i],p[j])+DotProd(p[j],p[j]));
		} // inertia moment is iner/12
		// this formula is based on expression for height between tangent plane and sphere, as r^2/2*R
		if(shape==SH_SPHERE){
			Rc=boxX/2; //radius of curvature is sphere radius
			vf-=iner/(24*Rc); //Information on the local curvature is needed to extend it to other shapes
		}
	}
	return vf;
}

//======================================================================================================================

void CalculationOfLsTensor(const double p[24][3], int k, const double n[3], double L[9]) {
	/*Determines the self-term dyadic of a single facet through the vector function h associated to some contour integral 
	on a polygon and the solid angle Ω as viewed from the origin.*/
	if (k<3) return;
	//Contour integral h over the polygon edges
	double h[3];
	double q[3]= {0,0,0}, du[3]={0,0,0};//Edge vector (and normalized)
	double lambda; //Weighting coefficients
	double u[24][3]; //Vertices of the cube for all the facets

	for (int j=0; j<k; j++)
			vSubtr(p[j],CubeCenter,u[j]);
	vSubtr(p[0],CubeCenter,u[k]); //Ensure that the formed polyhedron is closed

	for(int i=0; i<3; i++) h[i]=0.0;
	for (int j=0; j<k; j++){
		vSubtr(u[j+1], u[j], du);
		vMultScal(1.0 / vNorm(du), du, q);
		bool vertex=false; //Whether coefficients represent vertices contribution or not
		if (vertex==true) lambda=atanh(DotProd(u[j+1],q)/vNorm(u[j+1]))-atanh(DotProd(u[j],q)/vNorm(u[j]));
		else lambda=log((DotProd(u[j+1],q)+vNorm(u[j+1]))/(DotProd(u[j],q)+vNorm(u[j])));
		vMultScal(lambda,q,q);
		vAdd(h,q,h);
	}

	//calculation of omega
	double Omega;
	Omega=0.0;
	for (int i=0; i<k; i++) vMultScalSelf(1.0 / vNorm(u[i]), u[i]);
    for (int j=1; j<k-1; j++){
    	double prod[3];
    	CrossProd(u[j], u[j+1], prod);
    	double f = DotProd(u[0], prod);
    	double g = 1+DotProd(u[0],u[j])+DotProd(u[j],u[j+1])+DotProd(u[j+1],u[0]);
			Omega+=2*atan2(f,g);
    }
		//Calculation of Ls=Ωn⊗n+(n×h)⊗n
    double prod[3],cross[3];
    CrossProd(n,h,cross);
    vMultScal(Omega, n, prod);
    vAdd(prod, cross, prod);

	for (int i=0;i<3;i++){
		for (int j=0;j<3;j++){
			L[3*i+j]+=prod[i]*n[j];
		}
	}
}
//======================================================================================================================
void FiniteSizeCorrection(doublecomplex m,double vf,doublecomplex M0[9])
{
		doublecomplex M[9];
		double ka,kd2,S;
		int i;
		bool asym;
		const double *incPol;
		bool pol_avg=true;
		const enum incpol which=INCPOL_Y;
		doublecomplex RR=I*2*kd*kd*kd/3;
		extern const bool avg_inc_pol;
		extern const double polNlocRp;

		for (int i=0;i<9;i++) M[i]=0.0+0.0*I;

		asym = (PolRelation==POL_CLDR);
		if (asym && anisotropy) LogError(ONE_POS,"Incompatibility error in CoupleConstant");

		kd2=kd*kd;
		// Diagonal M
		if (asym) for (i=0;i<3;i++){ // loop over components of polarizability (for scalar input m)
			switch (PolRelation) {
				case POL_CLDR: M[4*i]=(LDR_B1+(LDR_B2+LDR_B3*prop[i]*prop[i])*m*m)*kd*kd; break;
				default: LogError(ONE_POS,"Incompatibility error in CoupleConstant");
			}
		}
		// Scalar M
		else for (i=0;i<3;i++) {
			switch (PolRelation) {
				case POL_CM: M[4*i]=0.0+0.0*I; break;
				case POL_DGF: M[4*i]=DGF_B1*kd2+RR; break;
				case POL_FCD: M[4*i]=2*ONE_THIRD*kd2*(2+kd*INV_PI*log((PI-kd)/(PI+kd)))+RR;
					break;
				case POL_IGT_SO: M[4*i]=SO_B1*kd2+RR; break;
				case POL_LAK: M[4*i]=2*FOUR_PI_OVER_THREE*((1-I*LAK_C*kd)*imExp(LAK_C*kd)-1); break;
				case POL_LDR:
					if (avg_inc_pol) S=0.5*(1-DotProdSquare(prop,prop));
					else {
						if (which==INCPOL_Y) incPol=incPolY;
						else incPol=incPolX;
						S = DotProdSquare(prop,incPol);
					}
					M[4*i]=(LDR_B1+(LDR_B2+LDR_B3*S)*m*m)*kd*kd+RR;
					break;
				case POL_NLOC:
					if (polNlocRp==0) M[4*i]=0.0+0.0*I;
					else M[4*i]=FOUR_PI_OVER_THREE*theta3(SQRT1_2PI*gridspace/polNlocRp);
					break;
				case POL_NLOC_AV:
					if (polNlocRp==0) M[4*i]=0.0+0.0*I;
					else {
						double x=gridspace/(2*SQRT2*polNlocRp);
						double g0,t;
						if (x<1) {
							t=erf(x);
							g0=1-t*t*t;
						}
						else {
							t=erfc(x);
							g0=t*(3-3*t+t*t);
						}
						// !!! dynamic part should be added here
						M[4*i]=FOUR_PI_OVER_THREE*g0;
					}
					break;
				case POL_RRC: M[4*i]=RR; break;
				default: LogError(ONE_POS,"Incompatibility error in CoupleConstant");
					// no break
			}
		}
		for (int i=0;i<9;i++) M0[i]=vf*M[i]; //First-order approximation in volume
}
//======================================================================================================================
void PolarizabilityCalc(doublecomplex mp, doublecomplex ms, double vf, double n[3],doublecomplex chi[3][3],doublecomplex alpha[3][3])
	/*Computes the effective polarizability based on both refractive indices for principal and secondary domains, volume
	fraction of the principal domain and the vector normal to the plane. Conditions on domain assignment are avoided in 
	this function. Total self-term dyadic (cube facets + plane) and other WD quantities are also computed here.*/
{
	double Eye[9] = {1,0,0,0,1,0,0,0,1};
	doublecomplex chi_p=(mp*mp-1.0)/FOUR_PI; // Susceptibility in principal domain
	doublecomplex chi_s=(ms*ms-1.0)/FOUR_PI; // Susceptibility in secondary domain
	doublecomplex T[3][3]; //Boundary condition tensor
	double Lp[9], LpMatr[3][3]; // Self-term dyadic in principal domain
	double Ls[9], LsMatr[3][3]; // Self-term dyadic in secondary domain
	doublecomplex Mp[9], MpMatr[3][3]; // Finite-size correction in principal domain
	doublecomplex Ms[9], MsMatr[3][3]; // Finite-size correction in secondary domain
	doublecomplex alphaT[3][3]; // Transpose of the polarizability tensor
	doublecomplex tmp[3][3], tempT[3][3], temp[3][3], tempP[3][3], tempS[3][3]; //Some 3x3 temporary matrices
	double trLp,trLs; //Normalized tr(Ls) and tr(Lp)
	double UnsortedEdgePoints[6][24][3], SortedEdgePoints[6][24][3];
	double intersections[24][3], pOrdered[24][3], pOrdered_inv[24][3];
	int counting = 0, accumCounting = 0, NumberOfPoints[6], edgeIntCounts[6]; // separate counts of intersections (by edges)

	/*The total self-term dyadic is computed as the sum of facet contributions (including the plane). For each facet, the 
	number of intersections is determined and both cube and intersection vertices are reordered. The self-term dyadic is
	then computed using these sorted points and the normal vector of the plane.*/
	
	if(vf==1.0){
		for (int l=0; l<9; l++){
			Ls[l]=0.0;
			Lp[l]=FOUR_PI_OVER_THREE*Eye[l];
		}
	}else{
		int k = CubePlaneIntersect(n, intersections, edgeIntCounts);
		k = ReorderPoints(intersections, k, n, pOrdered);
		for (int j=0; j<6; j++){
			counting=0;
			for (int i=0; i<4; i++){
				double* p = cubeOrderedPoints[cubeOrderedEdges[j][i]];
				if (DotProd(p,n)>1.0) vCopy(p, UnsortedEdgePoints[j][counting++]);
			}
			for (int l=0; l<edgeIntCounts[j]; l++) vCopy(intersections[accumCounting + l], UnsortedEdgePoints[j][counting + l]);

			counting += edgeIntCounts[j];
			accumCounting += edgeIntCounts[j];
			NumberOfPoints[j]=counting;
		}
		for (int i=0; i<6; i++)
			NumberOfPoints[i] = ReorderPoints(UnsortedEdgePoints[i], NumberOfPoints[i], CubeEdgeNorm[i], SortedEdgePoints[i]);
		/* Self-term dyadic Ls for each face of the cube...*/
		for (int l=0; l<9; l++) Ls[l]=0.0;
		for (int i=0; i<6; i++) 
			CalculationOfLsTensor(SortedEdgePoints[i], NumberOfPoints[i], CubeEdgeNorm[i], Ls);
		/*...and for the slice generated by the intersecting plane*/
		vNormalize(n);
		/*Normal vector of the facet associated to the intersecting plane is always outward to (s) while the computed normal
		vector of the plane always points toward (s) : the vertices are always CCW ordered here */
		for (int i=0; i<k; i++) vCopy(pOrdered[i],pOrdered_inv[(k-1)-i]);
		CalculationOfLsTensor(pOrdered_inv, k, n, Ls);
		for (int l=0; l<9; l++) Lp[l]=FOUR_PI_OVER_THREE*Eye[l]-Ls[l];
	}
	FiniteSizeCorrection(mp,vf,Mp);
	FiniteSizeCorrection(ms,1-vf,Ms);
	/*Here follows the routine for the calculation of effective polarizability tensor, based on principal domain p. The 
	domain assignment that is performed in CoupleConstant holds for volume fraction fp and refractive indices mp & ms.*/
	for (int i=0; i<3; i++){
		for (int j=0; j<3; j++){
			LsMatr[i][j]=Ls[3*i+j];
			LpMatr[i][j]=Lp[3*i+j];
			MsMatr[i][j]=Ms[3*i+j];
			MpMatr[i][j]=Mp[3*i+j];
		}
	}
	/*Computes the boundary condition tensor : T = I+(mp^2/ms^2-1)*nn'*/
	DyadProd(n,T);
	MatrMul(T, (mp*mp)/(ms*ms) - 1.0);
	MatrSum(T, Eye3);
	/*Computes the effective susceptibility tensor : χe=vf*χp*I+(1-vf)*χs*T*/
	MatrCopy(3, chi, Eye3);
	MatrCopy(3,tempT,T);
	MatrMul(chi, vf*chi_p);
	MatrMul(tempT,(1-vf)*chi_s);
	MatrSum(chi,tempT);
	/*Computes the effective polarizability tensor : αe = V*χe/[I+(Lp-Mp)*χp+(Ls-Ms)*χs*T]*/
	MatrCopy(3,tmp, Eye3);
	AlterMatrCopy(3, tempP, LpMatr);
	MatrDiff(tempP, MpMatr);
	MatrMul(tempP,chi_p);
	AlterMatrCopy(3,tempS,LsMatr);
	MatrDiff(tempS, MsMatr);
	MatrMul(tempS,chi_s);
	MatrProd(3,tempS,T,tempT);
	MatrSum(tmp,tempP);
	MatrSum(tmp,tempT);
	MatrInv(tmp,temp);
	MatrProd(3, chi, temp, alpha);
	MatrMul(alpha,dipvol);

	trLs=(LsMatr[0][0]+LsMatr[1][1]+LsMatr[2][2])/(4.0*PI);
	trLp=(LpMatr[0][0]+LpMatr[1][1]+LpMatr[2][2])/(4.0*PI);
	
}
//======================================================================================================================
void CoupleConstantWD(doublecomplex *mrel,const enum incpol which,doublecomplex chi[3][3],doublecomplex res[static 6], int index)
{
	doublecomplex alpha[3][3];
	doublecomplex mp,ms;
	double vf=volfrac[index];
	double r[3]={DipoleCoord[3*index],DipoleCoord[3*index+1],DipoleCoord[3*index+2]};
	double n[3]={plSec[3*index],plSec[3*index+1],plSec[3*index+2]};
	/*Domain assignment for scatterer : volume fraction, refractive indices*/
	bool cond=(volfrac[index]<1.0 && DotProd(n,r)<=gridspace); // Condition to fulfill for s=scatterer
	if(cond){
		mp=1.0+0.0*I;
		ms=ref_index[0];
	}else{
		mp=ref_index[0];
		ms=1.0+0.0*I;
	}
	PolarizabilityCalc(mp,ms,vf,n,chi,alpha);
	/*Symmetrization through upper/lower triangular matrix components*/
	res[0]=alpha[0][0]; res[1]=(alpha[0][1]+alpha[1][0])/2.0; res[2]=(alpha[2][0]+alpha[0][2])/2.0;
	res[3]=alpha[1][1]; res[4]=(alpha[1][2]+alpha[2][1])/2.0; res[5]=alpha[2][2];
}