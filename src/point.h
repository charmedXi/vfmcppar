#ifndef GUARD_POINTS_H
#define GUARD_POINTS_H

#include <array>
#include <cmath>
#include <iostream>

typedef std::array<double,3>  vec3d;
inline void zero(vec3d &v) {
    v[0]=0.0; v[1]=0.0; v[2]=0.0;
}

using namespace std;

class Point{
public:
	/* member data */
	Point 	*mPrev;				// pointer to previous point in filament
	Point 	*mNext;				// pointer to next point in filament
	vec3d 	mPos;				// position
	vec3d 	mVel; 				// current velocity
	vec3d 	mVel1;				// velocity last time step
	vec3d 	mVel2;  			// velocity 2 time steps ago
	vec3d 	mVel3;  			// velocity 3 time steps ago
	vec3d 	mVelNL; 			// non-local contributions to velocity
	vec3d 	mSPrime;			// tangent at point
	vec3d 	mS2Prime;			// binormal at point
	vec3d 	mSegLast; 			// vector to last point
	double	mSegLength; 		// distance to last point
	double 	mCharge; 			// charge at point
	int 	mFlagFilled;		// flag showing how many velocity steps back are present
	bool	mMarkedForRecon;	// flag used in reconnection
	bool	mFlagDummy;			// flag indicating dummy points used for strings
	bool	mMarkedForDeletion; // flag telling reconnect whether this point needs to be removed
	
	/* member functions */
	Point(){
		/* default constructor just reserves memory */
		mCharge = 0; mSegLength = 0; mFlagFilled = 0; mMarkedForDeletion = false; mMarkedForRecon = false;
		mFlagDummy = false;
	}
	explicit Point(Point* occ){
		/* parameterised constructor copies a point, used in reconnection */
		mVel = occ->mVel; mVel1 = occ->mVel1; mVel2 = occ->mVel2; mVel3 = occ->mVel3;
		zero(mVelNL);
		mPos = occ->mPos;
		zero(mSPrime); zero(mS2Prime); zero(mSegLast);
		mCharge = 0; mSegLength = 0; mFlagFilled = occ->mFlagFilled; mMarkedForDeletion = false; mMarkedForRecon = occ->mMarkedForRecon;
		mFlagDummy = false;
	}
	Point(const vec3d &CurrentPos) {
		/* parameterised constructor copies a point, used in reconnection */
		mPos=CurrentPos;
		zero(mVel); zero(mVel1); zero(mVel2); zero(mVel3); zero(mVelNL);
		zero(mSPrime); zero(mS2Prime); zero(mSegLast);
		mCharge = 0; mSegLength = 0; mFlagFilled = 0; mMarkedForDeletion = false; mMarkedForRecon = false;
		mFlagDummy = false;
	}
	~Point(){};

	double Disp2(Point* p2){
		double disp2(0);
		for(int i(0); i!=3; i++){
			disp2 += pow(mPos[i] - p2->mPos[i],2);
		}
		return disp2;
	}

};
#endif
