/* 	Calculates the non-local contributions to the velocities of every point.
	Adapted from CalcNonLocalVel_OtherFilament.m by Paul Walmsley. */

#include "tangle.h"

using namespace std;

const double kappa = 9.98e-8;

void Tangle::CalcVelocityNL(Point* pField){
	/* iterate over "source points" */
	//vector <double> cumuvel(3);
        vec3d cumuvel={0,0,0};
	for(int Q(0); Q!=mTangle.size(); Q++){
		for(int l(0); l!=mTangle[Q]->mN; l++){
			Point* pSource = mTangle[Q]->mPoints[l];
			if(pSource==pField||pSource==pField->mPrev){continue;}
			else{
				// p = s_l - s_k, q = s_l+1 - s_l
				vec3d q, p, pxq;	 
				/* calculate p and q */
                                vec3d fpos = pField->mPos;
                                vec3d ppos = pSource->mPos;
                                vec3d npos = pSource->mNext->mPos;

				// p.p, q.q, p.q
				double pp(0), qq(0), pq(0);	
				for(int m(0);m<3;m++){
					p[m] = ppos[m] - fpos[m];
					q[m] = npos[m] - fpos[m];
					pp += p[m]*p[m];
					qq += q[m]*q[m];
					pq += p[m]*q[m];
				}

				/* calculate pxq and assign temp variables */
				pxq[0] = p[1]*q[2] - p[2]*q[1];
				pxq[1] = p[2]*q[0] - p[0]*q[2];
				pxq[2] = p[0]*q[1] - p[1]*q[0];

                                double spp = sqrt(pp);
                                double sqq = sqrt(qq);
                                double sppqq = spp*sqq;
				double D = (spp+sqq)/(sppqq*(sppqq+pq));
				/* assign values to mVelNL */
                                D *= kappa/(4*M_PI);
				for(int j(0);j<3;j++){
					cumuvel[j] += D * pxq[j];
				}
			}
		}
	}
	pField->mVelNL = cumuvel;
}
