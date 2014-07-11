
#include "filament.h"

void Filament::CalcImpulse(){
//	cout << "Begin impulse" << endl;
	double r(0);
	vector < vector <double> > eps(mN), rxeps(mN);
	vector < double > r_eff(3);
	for (int k = 0; k < mN; k++){
		eps[k].resize(3); rxeps[k].resize(3);
		for (int j = 0; j < 3; j++){
			eps[k][j] = (mPoints[k]->mSegLength) * (mPoints[k]->mSPrime[j]) / (sqrt(mPoints[k]->mSPrime[0] * mPoints[k]->mSPrime[0] + mPoints[k]->mSPrime[1] * mPoints[k]->mSPrime[1] + mPoints[k]->mSPrime[2] * mPoints[k]->mSPrime[2]));
		}
		rxeps[k][0] = mPoints[k]->mPos[1] * eps[k][2] - mPoints[k]->mPos[2] * eps[k][1];
		rxeps[k][1] = mPoints[k]->mPos[2] * eps[k][0] - mPoints[k]->mPos[0] * eps[k][2];
		rxeps[k][2] = mPoints[k]->mPos[0] * eps[k][1] - mPoints[k]->mPos[1] * eps[k][0];
		for (int j = 0; j < 3; j++){
			r_eff[j] += rxeps[k][j];
		}
	}
	for (int j = 0; j < 3; j++){
		r_eff[j] = sqrt(pow(((1 / (2 * PI))*r_eff[j]), 2));
	}
	r = sqrt(r_eff[0]+r_eff[1]+r_eff[2]);
//	cout << "Effective radius = " << r << endl;
}