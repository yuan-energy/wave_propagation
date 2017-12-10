#pragma once
#include <iostream>
#include <fstream>
#include <complex>
#include <vector>
#include <numeric>
#include <limits>
#include <limits.h>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <unordered_map>
#include <exception>
#include <sstream>
#include <cstring>
#include <streambuf>
#include "__assert_msg.h"

// #include <Channel.h>
// #include <FEM_ObjectBroker.h>
// #include <Vector.h>
// #include <ID.h>

using std::vector;
using std::complex;
using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::real;
using std::size_t;
using std::pair;
using std::unordered_map;
using std::stringstream;
#define SOILD_THICKNESS_TOLERANCE 1E-5
class WaveField
{
public:
	const double pi = 3.14159265359 ;
	const double machina_epsilon = std::numeric_limits<double>::min() ; 
	typedef	std::vector<std::vector<std::complex<double>>> mat_complex;
	typedef	std::vector<std::complex<double>> vec_complex;

	WaveField();
	WaveField(int tag);
	int set_motions(
		string const& filename_acc 
		, string const& filename_dis
		, double acc_unit_scale = 1.0
		, double dis_unit_scale = 1.0
		);
	int set_motions(vector<double> const& times, vector<double> const& acc, vector<double> const& dis );
	int set_motion_depth(double motion_depth, double depth_unit_scale = 1.0 );

	int set_motion_compensation_time(double add_compensation_time);

	int set_soil_profile(string const& filename
		, double Vs_unit_scale = 1.0
		, double rho_unit_scale = 1.0
		, double damp_unit_scale = 1.0
		, double thickness_unit_scale = 1.0
		);
	int set_soil_profile(vector<double> const& Vs
		, vector<double> const& rho
		, vector<double> const& damp
		, vector<double> const& thicks
		);

	int add_compute_depth(double new_depth);


	int compute();               // conduct the wave propagation. Take the input motion as the total wave (upward + downward)
	int compute_upward();        // conduct the upward wave propagation. Take the input motion as the upward wave.
	int	deconvolution2bedrock(); // for equivalent rock outcropping: deconvolution to bedrock and reset motion.


	vector<double> const& get_times() ;
	vector<double> const& get_acc_by_depth(double depth) ;
	vector<double> const& get_vel_by_depth(double depth) ;
	vector<double> const& get_dis_by_depth(double depth) ;
	vector<double> get_acc_by_depth(double depth, size_t pos, size_t len ) ;
	vector<double> get_vel_by_depth(double depth, size_t pos, size_t len ) ;
	vector<double> get_dis_by_depth(double depth, size_t pos, size_t len ) ;
	vector<double> const& get_freqs() ;
	vector<double> const& get_acc_freq_by_depth(double depth) ;
	vector<double> const& get_vel_freq_by_depth(double depth) ;
	vector<double> const& get_dis_freq_by_depth(double depth) ;

	int write_wave_at_depth(double depth, string filename_prefix="wave") ;

	pair<vector<double>, vector<double>> time2freq(double time_step, vector<double> const& time_series);
	
	int get_Ntime() const;
	int getTag() const;
	void setTag(int newTag) ;

	// int sendSelf( int commitTag, Channel &theChannel);
	// int receiveSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker
                      // &theBroker  );
	
private:
	int _theTag;
	// **************************************
	// Input 
	// **************************************
	string _soil_profile_filename; 
	string _motion_filename; 

	// Soil Profile
	vector<double> _soil_Vs;
	vector<double> _soil_rho;
	vector<double> _soil_damp;
	vector<double> _soil_thick;
	vector<double> _soil_depth;

	// Input Motions (Acceleration)
	vector<double> _times;
	vector<double> _input_acc;
	vector<double> _input_dis;
	double _time_step ;
	double _motion_depth ;
	double _max_time ;
	double _max_freq ;
	int _Ntime;

	// **************************************
	// Output
	// **************************************
	vector<vector<double>> _soil_acc ; //[depthID][timeID]
	vector<vector<double>> _soil_vel ; //[depthID][timeID]
	vector<vector<double>> _soil_dis ; //[depthID][timeID]

	vector<double> _freqs;
	vector<vector<double>> _soil_acc_freq ; //[depthID][freqID] soil layers, half freqs. magnitude.
	vector<vector<double>> _soil_vel_freq ; //[depthID][freqID] soil layers, half freqs. magnitude.
	vector<vector<double>> _soil_dis_freq ; //[depthID][freqID] soil layers, half freqs. magnitude.


	// **************************************
	// Subroutines
	// **************************************
	template<typename val_type>
	int print_vec(std::string name, std::vector<val_type> const& v){
		std::cout<< name << "\n";
		for (int i = 0; i < (int)v.size() ; ++i){
			std::cout<<v[i] << " ";
		}
		std::cout<<std::endl;
		return 1;
	}

	template<typename val_type>
	string to_short_string(val_type val){
		std::stringstream out ; 
		out.precision(5);
		out << val;
		return out.str();
	}
	int write_series_to_file(
		vector<double> const& c1, 
		vector<double> const& c2, 
		string filename) const;
	vector<double> integrate_vel(double time_step, vector<double> const& acc);
	vector<double> integrate_dis(double time_step, vector<double> const& acc, vector<double> const& vel);
	
	int depth2layer(double depth);
	int wave_propagation(double max_freq, int N_freq
		, std::vector<double> const& Vs
		, std::vector<double> const& rho
		, std::vector<double> const& damp
		, std::vector<double> const& layer_thick
		, mat_complex* disp
	);

	int wave_propagation(double max_freq, int N_freq
		, std::vector<double> const& Vs
		, std::vector<double> const& rho
		, std::vector<double> const& damp
		, std::vector<double> const& layer_thick
		, mat_complex* disp
		, mat_complex* up
	);

	/* 
	 * Computes the discrete Fourier transform (DFT) of the given complex vector, storing the result back into the vector.
	 * The vector can have any length. This is a wrapper function.
	 */
	void fft(std::vector<std::complex<double> > &vec);
	
	
	/* 
	 * Computes the inverse discrete Fourier transform (IDFT) of the given complex vector, storing the result back into the vector.
	 * The vector can have any length. This is a wrapper function. This transform does not perform scaling, so the inverse is not a true inverse.
	 */
	void ifft_inner(std::vector<std::complex<double> > &vec);
	void ifft(std::vector<std::complex<double> > &vec);
	
	
	/* 
	 * Computes the discrete Fourier transform (DFT) of the given complex vector, storing the result back into the vector.
	 * The vector's length must be a power of 2. Uses the Cooley-Tukey decimation-in-time radix-2 algorithm.
	 */
	void transformRadix2(std::vector<std::complex<double> > &vec);
	
	
	/* 
	 * Computes the discrete Fourier transform (DFT) of the given complex vector, storing the result back into the vector.
	 * The vector can have any length. This requires the convolution function, which in turn requires the radix-2 FFT function.
	 * Uses Bluestein's chirp z-transform algorithm.
	 */
	void transformBluestein(std::vector<std::complex<double> > &vec);
	
	
	/* 
	 * Computes the circular convolution of the given complex vectors. Each vector's length must be the same.
	 */
	void fft_convolve(
		const std::vector<std::complex<double> > &vecx,
		const std::vector<std::complex<double> > &vecy,
		std::vector<std::complex<double> > &vecout);
	size_t reverseBits(size_t x, int n);
	string removeComments(string prgm);
	// int send_vector(int commitTag, Channel & theChannel, std::vector<double> const& data, 
	// 	std::string const& data_name );
	// int receive_vector(int commitTag, Channel & theChannel, std::vector<double> & data, 
	// 	std::string const& data_name );
};