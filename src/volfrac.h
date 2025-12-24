/*
 * volfrac.h
 *
 *  Created on: 20 ���. 2022 �.
 *      Author: Workstation
 */

#ifndef SRC_VOLFRAC_H_
#define SRC_VOLFRAC_H_

bool PointUpperPlane(const double p[static 3], const double plane_n[static 3]);
int CountEdgePlaneIntersections(int edgeIdx, const double n[3], const double intersections[2][3]);
void CalculationOfLsTensor(const double p[24][3], int k, const double n[3], double L[9], double h[3], double *Omega);
void PolarizabilityCalc(doublecomplex mp, doublecomplex ms, double vf, double n[3],doublecomplex alpha[3][3]);
bool EdgeIn(int i,int j,double nv[8],double res[3]);

#endif /* SRC_VOLFRAC_H_ */
