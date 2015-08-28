// calculate new positions after time step
// adapted from CalcVelMaster.m by Paul Walmsley

#include "tangle.h"
using namespace std;

// propagate positions using appropriate order method
void Tangle::PropagatePos(Point* pField){
	/* if velocities have already been filled 3 steps back, do AB4 */
	if(pField->mFlagFilled==4){
		for(int j(0);j<3; j++){
			pField->mPos[j] += (mDt/24)*(55*pField->mVel[j] - 59*pField->mVel1[j] + 37*pField->mVel2[j] - 9*pField->mVel3[j]);
		}
	}
	/* if velocities have already been filled 2 steps back, do AB3 */
	if(pField->mFlagFilled==3){
		for(int j(0);j<3; j++){
			pField->mPos[j] += (mDt/12)*(23*pField->mVel[j] - 16*pField->mVel1[j] + 5*pField->mVel2[j]);
		}
	}
	/* if velocities have only been filled 1 step back, do AB2 */
	else if(pField->mFlagFilled==2){
		for(int j(0);j<3; j++){
			pField->mPos[j] += (mDt/2)*(3*pField->mVel[j] - pField->mVel1[j]);
		}
	}
	/* if no previous velocities, use Euler */
	else if(pField->mFlagFilled==1){
		for(int j(0);j<3; j++){
			pField->mPos[j] += mDt * pField->mVel[j];
		}
	}

}



