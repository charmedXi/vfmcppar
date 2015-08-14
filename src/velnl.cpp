/* 	Calculates the non-local contributions to the velocities of every point.
	Adapted from CalcNonLocalVel_OtherFilament.m by Paul Walmsley. */

#include "tangle.h"

using namespace std;

const double kappa = 9.98e-8;

void Tangle::CalcVelocityNL(Point* pField){
	/* iterate over "source points" */
	vector <double> cumuvel(3);
	for(int Q(0); Q!=mTangle.size(); Q++){
		for(int l(0); l!=mTangle[Q]->mN; l++){
			Point* pSource = mTangle[Q]->mPoints[l];
			if(pSource==pField||pSource==pField->mPrev){continue;}
			else{
				// p = s_l - s_k, q = s_l+1 - s_l
				vector <double> q(3), p(3), pxq(3);	 
				// p.p, q.q, p.q
				double pp(0), qq(0), pq(0);	
				/* calculate p and q */
				for(int m(0);m<3;m++){
					p[m] = pSource->mPos[m] - pField->mPos[m];
					q[m] = pSource->mNext->mPos[m] - pField->mPos[m];
					pp += p[m]*p[m];
					qq += q[m]*q[m];
					pq += p[m]*q[m];
				}
				/* calculate pxq and assign temp variables */
				pxq[0] = p[1]*q[2] - p[2]*q[1];
				pxq[1] = p[2]*q[0] - p[0]*q[2];
				pxq[2] = p[0]*q[1] - p[1]*q[0];
				double D = (sqrt(pp)+sqrt(qq))/(sqrt(pp)*sqrt(qq)*(sqrt(pp)*sqrt(qq)+pq));
				/* assign values to mVelNL */
				for(int j(0);j<3;j++){
					cumuvel[j] += (kappa/(4*M_PI)) * D * pxq[j];
				}
			}
		}
	}
	pField->mVelNL = cumuvel;
}
