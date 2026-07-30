/*
 * WD.h
 *
 *  Created on: 20 ���. 2022 �.
 *      Author: Workstation
 */

#ifndef SRC_WD_H_
#define SRC_WD_H_

bool PointUpperPlane(const double p[static 3], const double plane_n[static 3]);
void IntersectionParams(double nv[8], double n[3]);
int CubePlaneIntersect(double n[3], double intersections[24][3], int edgeCounts[6]);
void PointsCenter(double p[24][3], int n, double c[3]);
double PointWalkOrd(double u0[3], double u[3], const double n[3]);
int ReorderPoints(double p[24][3], int k, const double n[3], double pReordered[24][3]);
bool EdgeIn(int i,int j,double nv[8],double res[3]);
double VolumeFraction(double a,double b,double c);
void SelfTermDyadic(double p[24][3], int k, const double n[3], double L[9]);
void FiniteSizeCorrection(doublecomplex m,double vf,doublecomplex Mp[9],doublecomplex Ms[9]);
void PolarizabilityCalc(doublecomplex mp, doublecomplex ms, double vf, double n[3],doublecomplex chi[3][3],doublecomplex alpha[3][3]);
void CoupleConstantWD(doublecomplex *mrel,doublecomplex chi[3][3],doublecomplex res[static 6], int index);

#endif /* SRC_WD_H_ */
