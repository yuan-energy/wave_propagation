#include "WaveField.h"

#include <fstream>
using namespace std;
int main(int argc, char const *argv[])
{

	// *********************************************
	// Input
	// *********************************************
	string acc_fileName = "scaled_northridge_acc.dat";
	string dis_fileName = "scaled_northridge_dis.dat";
	string soil_profile_fileName = "soil_profile.txt";
	double motion_depth = 0. ;
	WaveField theField ; 
	theField.set_motions(acc_fileName, dis_fileName);
	theField.set_motion_depth(motion_depth);
	theField.set_soil_profile(soil_profile_fileName);

	// For equivalent rock outcropping.
	// theField.deconvolution2bedrock();
	
	// *********************************************
	// Compute
	// *********************************************
	vector<double> request_depths = {10, 30, 50, 70};
	for (auto& d: request_depths){
		theField.add_compute_depth(d);
	}


	theField.compute();




	// *********************************************
	// Ouput
	// *********************************************
	double request_depth = 70 ;
	theField.write_wave_at_depth(request_depth);

	return 0;
}

