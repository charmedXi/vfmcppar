/* Calculate velocities of points in mPoints of a filament.
   Adapted from CalcVelMaster.m by Paul Walmsley */

#include "filament.h"
#include "tangle.h"

using namespace std;
/* circulation quantum, core radius */
const double	kappa = 9.98e-8, a0=1.3e-10, a1=exp(0.5)*a0;

/* calculate velocity at each point from s' and s'' eq(2) from Hanninen and Baggaley PNAS 111 p4667 (2014) */
void Tangle::CalcVelocity(Point* pField){
	// set velocity at previous time steps
	for(int j(0);j!=3;j++){
		pField->mVel3[j] = pField->mVel2[j];
		pField->mVel2[j] = pField->mVel1[j];
		pField->mVel1[j] = pField->mVel[j];
	}
	// calculate SPrime for just the point in question
	double A, B, C, D, E;
	double l, l1, l2, lm1;
	l = pField->mSegLength; l1 = pField->mNext->mSegLength;
	l2 = pField->mNext->mNext->mSegLength; lm1 = pField->mPrev->mSegLength;
	A = l * l1 * l1 + l * l1 * l2;
	A /= (lm1 * (lm1 + l) * (lm1 + l + l1) * (lm1 + l + l1 +l2));
	B = -lm1 * l1 * l1 - l * l1 * l1 - lm1 * l1 * l2 - l * l1 * l2;
	B /= (lm1 * l * (l + l1) * (l + l1 + l2));
	D = lm1 * l * l1 + l * l * l1 + lm1 * l * l2 + l * l * l2;
	D /= (l1 * l2 * (l + l1) * (lm1 + l + l1));
	E = -l1 * l * l - lm1 * l * l1;
	E /= (l2 * (l1 + l2) * (l + l1 + l2) * (lm1 + l + l1 + l2));
	C = -(A + B + D + E);

	for(int q=0;q<3;q++){
		pField->mSPrime[q]  = A*pField->mPrev->mPrev->mPos[q];
		pField->mSPrime[q] += B*pField->mPrev->mPos[q];
		pField->mSPrime[q] += C*pField->mPos[q];
		pField->mSPrime[q] += D*pField->mNext->mPos[q];
		pField->mSPrime[q] += E*pField->mNext->mNext->mPos[q];
	}
	// calculate S2Prime for just the point in question
	A = 2 * (-2 * l * l1  +  l1 * l1  - l * l2  +  l1 * l2 );
	A = A / (lm1 * (lm1 +l)*(lm1 + l + l1)*(lm1 + l + l1 + l2));
	B = 2 * (2 * lm1 * l1  + 2 * l * l1  -  l1 * l1  +  lm1 * l2  + l * l2  -  l1 * l2);
	B = B / (lm1 * l *(l + l1) * (l + l1 + l2));
	D = 2*(-lm1 * l - l * l +  lm1 * l1  + 2 * l * l1  +  lm1 * l2  + 2 * l * l2);
	D = D / (l1 * l2 *(l + l1) * (lm1 + l + l1));
	E = 2*(lm1 * l + l * l -  lm1 * l1  - 2 * l * l1 );
	E = E / ( l2  * (l1 + l2) * (l+ l1 + l2) * (lm1 + l + l1 + l2));
	C = -(A + B + D + E);

	for(int q=0;q<3;q++){
		pField->mS2Prime[q]  = A*pField->mPrev->mPrev->mPos[q];
		pField->mS2Prime[q] += B*pField->mPrev->mPos[q];
		pField->mS2Prime[q] += C*pField->mPos[q];
		pField->mS2Prime[q] += D*pField->mNext->mPos[q];
		pField->mS2Prime[q] += E*pField->mNext->mNext->mPos[q];
	}

	if(pField->mFlagFilled!=5){
		if(pField->mFlagFilled==3){pField->mFlagFilled++;}
		if(pField->mFlagFilled==2){pField->mFlagFilled++;}
		if(pField->mFlagFilled==1){pField->mFlagFilled++;}
		if(pField->mFlagFilled==0){pField->mFlagFilled++;}
		pField->mVel[0] = ((pField->mSPrime[1])*(pField->mS2Prime[2]) - (pField->mSPrime[2])*(pField->mS2Prime[1]));
		pField->mVel[1] = ((pField->mSPrime[2])*(pField->mS2Prime[0]) - (pField->mSPrime[0])*(pField->mS2Prime[2]));
		pField->mVel[2] = ((pField->mSPrime[0])*(pField->mS2Prime[1]) - (pField->mSPrime[1])*(pField->mS2Prime[0]));
		/* apply prefactor and add non-local contributions, resetting mVelNL after */
		for(int q=0;q<3;q++){
			pField->mVel[q] *= kappa*log(2*sqrt(pField->mSegLength * pField->mNext->mSegLength)/a1)/(4*PI);
			pField->mVel[q] += pField->mVelNL[q];   
			pField->mVelNL[q] = 0;		
		}
	}
}

// calculate s' using coefficients from Baggaley & Barenghi JLT 166:3-20 (2012)
void Filament::CalcSPrime(){
	vector <double> A(mN), B(mN), C(mN), D(mN), E(mN);
	double l, l1, l2, lm1;
	for(int i=0;i<mN;i++){
		l = mPoints[i]->mSegLength; l1 = mPoints[i]->mNext->mSegLength;
		l2 = mPoints[i]->mNext->mNext->mSegLength; lm1 = mPoints[i]->mPrev->mSegLength;
		A[i] = l * l1 * l1 + l * l1 * l2;
		A[i] /= (lm1 * (lm1 + l) * (lm1 + l + l1) * (lm1 + l + l1 +l2));

		B[i] = -lm1 * l1 * l1 - l * l1 * l1 - lm1 * l1 * l2 - l * l1 * l2;
		B[i] /= (lm1 * l * (l + l1) * (l + l1 + l2));

		D[i] = lm1 * l * l1 + l * l * l1 + lm1 * l * l2 + l * l * l2;
		D[i] /= (l1 * l2 * (l + l1) * (lm1 + l + l1));

		E[i] = -l1 * l * l - lm1 * l * l1;
		E[i] /= (l2 * (l1 + l2) * (l + l1 + l2) * (lm1 + l + l1 + l2));

		C[i] = -(A[i] + B[i] + D[i] + E[i]);
	}
	for(int p=0;p<mN;p++){
		for(int q=0;q<3;q++){
			mPoints[p]->mSPrime[q]  = A[p]*mPoints[p]->mPrev->mPrev->mPos[q];
			mPoints[p]->mSPrime[q] += B[p]*mPoints[p]->mPrev->mPos[q];
			mPoints[p]->mSPrime[q] += C[p]*mPoints[p]->mPos[q];
			mPoints[p]->mSPrime[q] += D[p]*mPoints[p]->mNext->mPos[q];
			mPoints[p]->mSPrime[q] += E[p]*mPoints[p]->mNext->mNext->mPos[q];
		}
	}
}

// calculate s'' using coefficients from Baggaley & Barenghi JLT 166:3-20 (2012)
void Filament::CalcS2Prime(){
	vector <double> A2(mN), B2(mN), C2(mN), D2(mN), E2(mN);
	double l, l1, l2, lm1;
	for(int i=0;i<mN;i++){

		l = mPoints[i]->mSegLength; l1 = mPoints[i]->mNext->mSegLength;
		l2 = mPoints[i]->mNext->mNext->mSegLength; lm1 = mPoints[i]->mPrev->mSegLength;

		A2[i] = 2 * (-2 * l * l1  +  l1 * l1  - l * l2  +  l1 * l2 );
		A2[i] = A2[i] / (lm1 * (lm1 +l)*(lm1 + l + l1)*(lm1 + l + l1 + l2));

		B2[i] = 2 * (2 * lm1 * l1  + 2 * l * l1  -  l1 * l1  +  lm1 * l2  + l * l2  -  l1 * l2);
		B2[i] = B2[i] / (lm1 * l *(l + l1) * (l + l1 + l2));

		D2[i] = 2*(-lm1 * l - l * l +  lm1 * l1  + 2 * l * l1  +  lm1 * l2  + 2 * l * l2);
		D2[i] = D2[i] / (l1 * l2 *(l + l1) * (lm1 + l + l1));

		E2[i] = 2*(lm1 * l + l * l -  lm1 * l1  - 2 * l * l1 );
		E2[i] = E2[i] / ( l2  * (l1 + l2) * (l+ l1 + l2) * (lm1 + l + l1 + l2));

		C2[i] = -(A2[i] + B2[i] + D2[i] + E2[i]);
	}
	for(int p=0;p<mN;p++){
		for(int q=0;q<3;q++){
			mPoints[p]->mS2Prime[q]  = A2[p]*mPoints[p]->mPrev->mPrev->mPos[q];
			mPoints[p]->mS2Prime[q] += B2[p]*mPoints[p]->mPrev->mPos[q];
			mPoints[p]->mS2Prime[q] += C2[p]*mPoints[p]->mPos[q];
			mPoints[p]->mS2Prime[q] += D2[p]*mPoints[p]->mNext->mPos[q];
			mPoints[p]->mS2Prime[q] += E2[p]*mPoints[p]->mNext->mNext->mPos[q];
		}
	}
}
