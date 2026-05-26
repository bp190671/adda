/*
 * WD.h
 *
 *  Created on: 20 ���. 2022 �.
 *      Author: Workstation
 */

#ifndef SRC_WD_H_
#define SRC_WD_H_

bool PointUpperPlane(const double p[static 3], const double plane_n[static 3]);
int CountEdgePlaneIntersections(int edgeIdx, const double n[3], const double intersections[2][3]);
void CalculationOfLsTensor(const double p[24][3], int k, const double n[3], double L[9]);
void PolarizabilityCalc(doublecomplex mp, doublecomplex ms, double vf, double n[3],doublecomplex chi[3][3],doublecomplex alpha[3][3]);
void CoupleConstantWD(doublecomplex *mrel,const enum incpol which,doublecomplex chi[3][3],doublecomplex res[static 6], int index);
bool EdgeIn(int i,int j,double nv[8],double res[3]);
double VolumeFraction(double a,double b,double c);

#endif /* SRC_WD_H_ */
