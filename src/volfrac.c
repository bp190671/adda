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
#include "volfrac.h"
#include "QR.h"
#include "tagaki_factor.h"

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
	return DotProd(diff, diff) < DBL_EPS;
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

void FillCubePlaneIntHelperArray(double nv[8], double n[3]) {
	/*Helper for storing the intersections of vertices with the plane in an array*/
	for (int i = 0; i < 8; i++) nv[i] = DotProd(cubeOrderedPoints[i], n);
}

//======================================================================================================================

int CubePlaneRawIntersections(const double n[3], double intersections[24][3], int edgeCounts[6]) {
	/*Determines the number of intersections per edge within the cube*/
	int i_idx = 0;
	double nv[8];

	FillCubePlaneIntHelperArray(nv, n);
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

void CalculationOfLsTensor(const double p[24][3], int k, const double n[3], double L[9], double h[3], double *Omega) {
	/*Determines the self-term dyadic of a single facet through the vector function h associated to some contour integral 
	on a polygon and the solid angle Ω as viewed from the origin.*/
	if (k<3) return;
	//Contour integral h over the polygon edges

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
	*Omega=0.0;
	for (int i=0; i<k; i++) vMultScalSelf(1.0 / vNorm(u[i]), u[i]);
    for (int j=1; j<k-1; j++){
    	double prod[3];
    	CrossProd(u[j], u[j+1], prod);
    	double f = DotProd(u[0], prod);
    	double g = 1+DotProd(u[0],u[j])+DotProd(u[j],u[j+1])+DotProd(u[j+1],u[0]);
			*Omega+=2*atan2(f,g);
    }
		//Calculation of Ls=Ωn⊗n+(n×h)⊗n
    double prod[3],cross[3];
    CrossProd(n,h,cross);
    vMultScal(*Omega, n, prod);
    vAdd(prod, cross, prod);

	for (int i=0;i<3;i++){
		for (int j=0;j<3;j++){
			L[3*i+j]+=prod[i]*n[j];
		}
	}
}
//======================================================================================================================
void PolarizabilityCalc(doublecomplex mp, doublecomplex ms, double vf, double n[3],doublecomplex alpha[3][3])
	/*Computes the effective polarizability based on both refractive indices for principal and secondary domains, volume
	fraction of the principal domain and the vector normal to the plane. Conditions on domain assignment are avoided in 
	this function. Total self-term dyadic (cube facets + plane) and other WD quantities are also computed here.*/
{
	double Eye[9] = {1,0,0,0,1,0,0,0,1};
	doublecomplex chi_p=(mp*mp-1.0)/FOUR_PI; // Susceptibility in principal domain
	doublecomplex chi_s=(ms*ms-1.0)/FOUR_PI; // Susceptibility in secondary domain
	doublecomplex T[3][3]; //Boundary condition tensor
	doublecomplex chi[3][3]; //Effective susceptibility tensor
	double Lp[9], LpMatr[3][3]; // Self-term dyadic in principal domain
	double Ls[9], LsMatr[3][3]; // Self-term dyadic in secondary domain
	doublecomplex M=0.0; //Polarizability prescription (fixed to CM for now)
	doublecomplex Mp[9], MpMatr[3][3]; // Finite-size correction in principal domain
	doublecomplex Ms[9], MsMatr[3][3]; // Finite-size correction in secondary domain
	doublecomplex alphaT[3][3]; // Transpose of the polarizability tensor
	doublecomplex tmp[3][3], tempT[3][3], temp[3][3], tempP[3][3], tempS[3][3]; //Some 3x3 temporary matrices
	double h[3]; //Contour integral on the polygon
	double Omega; //Solid angle of the facet
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
		int k = CubePlaneRawIntersections(n, intersections, edgeIntCounts);
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
			CalculationOfLsTensor(SortedEdgePoints[i], NumberOfPoints[i], CubeEdgeNorm[i], Ls, h, &Omega);
		/*...and for the slice generated by the intersecting plane*/
		vNormalize(n);
		/*Normal vector of the facet associated to the intersecting plane is always outward to (s) while the computed normal
		vector of the plane always points toward (s) : the vertices are always CCW ordered here */
		for (int i=0; i<k; i++) vCopy(pOrdered[i],pOrdered_inv[(k-1)-i]);
		CalculationOfLsTensor(pOrdered_inv, k, n, Ls, h, &Omega);
		for (int l=0; l<9; l++) Lp[l]=FOUR_PI_OVER_THREE*Eye[l]-Ls[l];
	}
	for (int l=0; l<9; l++){
		Mp[l]=M*Eye[l];
		Ms[l]=M*Eye[l];
	}
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

	/*Print weighted discretization quantities to file*/
	if(print_wd){
		fprintf(voxel_wd,"Computational geometry : \n");
		fprintf(voxel_wd,"Normal vector n = (%.6e, %.6e, %.6e)\n", n[0], n[1], n[2]);
		fprintf(voxel_wd,"Volume fraction fp = %.6e\n\n", vf);
		fprintf(voxel_wd,"Principal domain (p) : \n");
		fprintf(voxel_wd,"Refractive index mp = %.6e+%.6ei\n", creal(mp), cimag(mp));
		fprintf(voxel_wd,"Electric susceptibility χp = %.6e+%.6ei\n\n", creal(chi_p), cimag(chi_p));
		fprintf(voxel_wd,"Self-term dyadic Lp = \n");
		DblDebugMatr(3,voxel_wd,LpMatr);
		fprintf(voxel_wd,"Normalized trace : tr(Lp)/(4π)=%.10f\n\n",trLp);
		fprintf(voxel_wd,"Finite-size correction Mp = \n");
		DebugMatr(3,voxel_wd,MpMatr);	
		fprintf(voxel_wd,"Secondary domain (s) : \n");
		fprintf(voxel_wd,"Refractive index ms = %.6e+%.6ei\n", creal(ms), cimag(ms));
		fprintf(voxel_wd,"Electric susceptibility χs = %.6e+%.6ei\n\n", creal(chi_s), cimag(chi_s));
		fprintf(voxel_wd,"Self-term dyadic Ls = \n");
		DblDebugMatr(3,voxel_wd,LsMatr);
		fprintf(voxel_wd,"Normalized trace : tr(Ls)/(4π)=%.10f\n\n",trLs);
		fprintf(voxel_wd,"Finite-size correction Ms = \n");
		DebugMatr(3,voxel_wd,MsMatr);
		fprintf(voxel_wd,"Weighted discretization : \n");
		fprintf(voxel_wd,"Boundary condition tensor T = \n");
		DebugMatr(3,voxel_wd,T);
		fprintf(voxel_wd,"Effective susceptibility χe = \n");
		DebugMatr(3,voxel_wd,chi);
		fprintf(voxel_wd,"Effective polarizability αe = \n");
		DebugMatr(3,voxel_wd,alpha);
		fprintf(voxel_wd,"-----------------------------------------------------------------------------------------\n\n");
	}
}