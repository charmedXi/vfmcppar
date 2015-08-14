#include "filament.h"
#include "tangle.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <iomanip>
#include <ctime>
#include <omp.h>

using namespace std;

	/* circulation quantum, core radii, mutual friction */
	double		kappa = 9.98e-8, a0=1.3e-10, a1=exp(0.5)*a0, alpha=0, q_e=1.6e-19;
	
int main(int argc, char* argv[]){

	/* initialise tangle  */
	Tangle Tangle;
	
	string runfile;
	if(argc!=1){runfile = argv[1];}
	else runfile = "init/ring_line.in";
	
	string filename = Tangle.Initialise(runfile);

	/* set number of timesteps and number of steps per save */
	int N_t(Tangle.mTotalTime/Tangle.mDt); // number of time steps
	Tangle.mN_f = 10000; // number of time steps per save
	Tangle.mN_slow = 0; // counts how many steps have occurred at slow-mo

	vector <Filament*>::iterator begin, current, end;

	/* prepare to time calculations */
	double percent;
	double us_Dt(Tangle.mDt * 1e6);
	clock_t t, t_temp;
	t=clock();
	int file_no(0);
	/* begin time-stepping */
	int i(0);
	cout << "\n\t - - - - - - - -    BEGINNING SIMULATION    - - - - - - - -\n\n";
	Tangle.mLog << Tangle.StringTime() << "\t\t\t\t\tsimulation begins" << endl;

	#pragma omp parallel
	{
		while(i < N_t){
			#pragma omp master
			{
				Tangle.mStep = i;

				/* save positions to file every mN_f steps */
				
				if(i%Tangle.mN_f==0){
					Tangle.Output(filename, i, file_no);
					t_temp = clock() -t;
					printf("\t\t wrote step %6u", i);
					Tangle.mLog << Tangle.StringTime() << "\t" << setw(10) << Tangle.mStep;
					Tangle.mLog << "\telapsed: " << ((float)t_temp)/CLOCKS_PER_SEC << " s:\t\twrote to file " << file_no << " for time " << i*us_Dt << " us" << endl;
					file_no++; 
				}
				if(i%100==0 || i%Tangle.mN_f==0){
					percent = (100*i/N_t);
					printf("\t\t\r %6.2f %% \t",percent); // output percentage completion
				}

				/* check for and perform reconnections if required */
				Tangle.Reconnection();
				/* adjust mesh until finished */
				bool MeshFinished(false);	
				int MeshCount(0);
				while(MeshFinished==false){
					MeshFinished = Tangle.MeshAdjust(); MeshCount++; 
					if(MeshCount == 10){
						cout << "Mesh too unruly! Exiting..." << endl;
						Tangle.mLog	<< Tangle.StringTime() << "Mesh too unruly! Program terminated." << endl;
						exit;
					}
				}
				/* remove rings smaller than 6 points and count them */
				bool LoopKilled = Tangle.LoopKill(); 
				if(LoopKilled == true) Tangle.mN_loopkills++;
			}
			/* calculate velocities and propagate positions */
			for(int P=0; P<Tangle.mTangle.size(); P++){
				#pragma omp for schedule(static) nowait
				for(int k=0; k<Tangle.mTangle[P]->mN; k++){
					Point* pField = Tangle.mTangle[P]->mPoints[k];
					Tangle.CalcVelocityNL(pField);	// calculates non-local contributions to velocity
				}
			}
			#pragma omp master
			{
				Tangle.CalcVelocity(); 		// calculates local contributions to velocity
				Tangle.PropagatePos(Tangle.mDt);	// propagate positions
				i++;	// step forward
			}
		}
	}

	cout << "\n\t - - - - - - -    SIMULATION FINISHED    - - - - - - - -" << endl;
	cout << "Leaf has size = " << sizeof(Tangle) << endl;
	Tangle.mLog << Tangle.StringTime() << "\t\t\t\tsimulation finished" << endl;
	ofstream timefile(filename+"/time.dat");
	t = clock()-t;
	timefile << "time elapsed = " << ((float)t)/CLOCKS_PER_SEC << " s " << endl;
	timefile << "number of recons = " << Tangle.mN_recon << endl;
	timefile << "number of loop kills = " << Tangle.mN_loopkills << endl;
	Tangle.mLog	<< Tangle.StringTime() << "\t\t\t\ttime elapsed = " << ((float)t)/CLOCKS_PER_SEC << " s " << endl;
	Tangle.mLog << Tangle.StringTime() << "\t\t\t\tnumber of recons = " << Tangle.mN_recon << endl;
	Tangle.mLog << Tangle.StringTime() << "\t\t\t\tnumber of loop kills = " << Tangle.mN_loopkills << endl;
	timefile.close(); 
	Tangle.mLog.close();
	return 0;
}
