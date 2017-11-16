#include "WaveField.h"

#include <fstream>
using namespace std;
int main(int argc, char const *argv[])
{

	// *********************************************
	// Input
	// *********************************************
	string acc_fileName = "ricker_acc.dat";
	string dis_fileName = "ricker_dis.dat";
	string soil_profile_fileName = "soil_profile.dat";
	double motion_depth = 250. ;
	WaveField theField ; 
	theField.set_motions(acc_fileName, dis_fileName);
	theField.set_motion_depth(motion_depth);
	theField.set_soil_profile(soil_profile_fileName);

	// For equivalent rock outcropping.
	// theField.deconvolution2bedrock();
	
	// *********************************************
	// Compute
	// *********************************************
	for (int i = 0; i <= 250 ; i+=25 ){
		theField.add_compute_depth(i);
	}


	theField.compute_upward();
	// theField.compute();




	// *********************************************
	// Ouput
	// *********************************************
	for (int i = 0; i <= 250 ; i+=25 ){
		theField.write_wave_at_depth(i);
	}

	return 0;
}



// {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240}
// 260
// 270
// 280
// 290
// 300
// 310
// 320
// 330
// 340
// 350
// 360
// 370
// 380
// 390
// 400
// 410
// 420
// 430
// 440
// 450
// 460
