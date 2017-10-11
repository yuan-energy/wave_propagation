#include "WaveField.h"

using std::vector;
using std::complex;
using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::real;
using std::size_t;
using std::pair;
using std::make_pair;
using std::ifstream;
using std::ofstream;
using std::stringstream;
using std::to_string;


WaveField::WaveField(){}
WaveField::WaveField(int tag)
: _theTag(tag)
{
}

int WaveField::set_motion_depth(double motion_depth
	, double depth_unit_scale)
{
	if( motion_depth < - machina_epsilon ){
		cerr << " ERROR!!! in WaveField::add_compute_depth.  \n";
		cerr << " The motion_depth cannot be smaller than zero! \n";
		cerr << " motion_depth = " << motion_depth << endl;
		return -1;
	}
	_motion_depth = motion_depth * depth_unit_scale;
	return 1;
}

int WaveField::set_motions(string const& filename
	, double acc_unit_scale
	){
	ifstream the_motionFile;
	the_motionFile.open(filename.c_str());
	if (!the_motionFile.is_open()){
	    cerr << "ERROR!!! - WaveField::set_motions()";
	    cerr << " - could not open motion file: " << filename << endl;
	    return -1;
	}
	string word; 
	int count = 0;
	vector<double> times, acceleration ;
	while(the_motionFile >> word){
		if(count&1){
			acceleration.push_back(stof(word));
		}else{
			times.push_back(stof(word));
		}
		count++;
	}
	if(acceleration.size() != times.size()){
		cerr<< "ERROR!!! in User input file: "<< filename <<" \n\t time - acc length mismatch! " << endl;
	}
	for (auto& item: acceleration){
		item *= acc_unit_scale;
	}
	return this->set_motions(times, acceleration);
	
}

int WaveField::set_soil_profile(string const& filename
	, double Vs_unit_scale
	, double rho_unit_scale
	, double damp_unit_scale
	, double thickness_unit_scale
	){
	std::ifstream theFile(filename);
	if (!theFile.is_open()){
	    cerr << "ERROR!!! - WaveField::set_soil_profile()";
	    cerr << " - could not open soil profile file " << filename << endl;
	    return -1;
	}
	// convert the whole file into a big string
	std::stringstream buffer;
	buffer << theFile.rdbuf();
	string input = buffer.str();

	// remove comments
	string cleaned_input = this->removeComments(input);

	// save results to vector
	stringstream ss(cleaned_input);
	string word, line; 
	vector<double> Vs, rho, damp, thickness;
	int count = 0;
	while(ss >> word){
		switch ( (count++)%4 ){
			case 0 :{
				Vs.push_back(stof(word));
				break;
			}
			case 1:{
				rho.push_back(stof(word));
				break;
			}
			case 2 :{
				damp.push_back(stof(word));
				break;
			}
			case 3:{
				thickness.push_back(stof(word));
				break;
			}
		}
	}

	for (auto& item: Vs       ){ item *= Vs_unit_scale;	       }
	for (auto& item: rho      ){ item *= rho_unit_scale;	   }
	for (auto& item: damp     ){ item *= damp_unit_scale;	   }
	for (auto& item: thickness){ item *= thickness_unit_scale; }

	// print_vec("Vs = ", Vs);
	// print_vec("rho = ", rho);
	// print_vec("damp = ", damp);
	// print_vec("thickness = ", thickness);

	return this->set_soil_profile(Vs, rho, damp, thickness) ;
}


 
int WaveField::set_soil_profile(vector<double> const& Vs
		, vector<double> const& rho
		, vector<double> const& damp
		, vector<double> const& thick
	)
{
	if (Vs.empty() || rho.empty() || damp.empty() || thick.empty()){
		cerr<<" ERROR!! in WaveField::set_soil_profile: at least one soil properties is empty" << endl;
		return -1;
	}
	if(!(
		 Vs.size() == rho.size() && rho.size() == damp.size()
		 && damp.size() == (thick.size()+1)
		)
	){
		cerr << " ERROR!!! in WaveField::set_soil_profile. Soil properties Lengths mismatch!"<<"\n" ;
		cerr << " Vs.size()     = " << Vs.size() << "\n" ;
		cerr << " rho.size()    = " << rho.size() << "\n" ;
		cerr << " damp.size()   = " << damp.size() << "\n" ;
		cerr << " thick.size() = " << thick.size() << "\n" ;
		return -1;
	}
	if (fabs(thick[0]) < machina_epsilon){
		cerr << " ERROR!!! in WaveField::set_soil_profile. The first thickness cannot be zero! \n";
		cerr << " thick[0] = " << thick[0] <<endl;
		return -1;
	}
	for (int i = 0; i < (int)thick.size() ; ++i){
		if( thick [i] < 0 ){
			cerr << " ERROR!!! in WaveField::set_soil_profile. The thickness cannot be zero! \n";
			cerr << " thick["<<i<<"] = " << thick[i] << endl;
			return -1;
		}
	}
	_soil_Vs = Vs;
	_soil_rho = rho;
	_soil_damp = damp;
	_soil_thick = thick;

	int layernum = _soil_thick.size();
	_soil_depth.resize(layernum) ;
	_soil_depth[0] = _soil_thick[0];
	for(int i=1; i<layernum; ++i){
		_soil_depth[i] = _soil_thick[i] + _soil_thick[i-1];
	}
	_soil_depth.insert(_soil_depth.begin(), 0);
	if(fabs(_motion_depth) > machina_epsilon ){
		this->add_compute_depth(_motion_depth);
	}
	return 1;
}

int WaveField::set_motions(vector<double> const& times, vector<double> const& acceleration)
{
	_times = times;
	_input_acc = acceleration;
	_Ntime = times.size() ;

	if( _Ntime < 2 ){
		cerr<<" ERROR in WaveField::set_motions. Too short time series _Ntime = "<< _Ntime <<"\n";
		cerr<<" input times is " ; 
		for (int i = 0; i < (int)times.size() ; ++i){ cout<< times[i] << " ";}
		cout<<endl;
		return -1;
	}
	_time_step = times[1] - times[0] ; 
	if (_time_step < machina_epsilon ){
		cerr<<" ERROR in WaveField::set_motions. Too short _time_step = "<< _time_step <<"\n";
		return -1;
	}
	_max_freq = 1./_time_step ; 
	_max_time = times.back() ; 

	if( fabs(_motion_depth) > machina_epsilon  && !_soil_Vs.empty()){
		this->add_compute_depth(_motion_depth);
	}

	return 1;
}

int WaveField::add_compute_depth(double new_depth){
	for (auto d: _soil_depth){
		if ( fabs(new_depth - d) <= SOILD_THICKNESS_TOLERANCE )	{
			// cout<<"hi" <<endl;
			// print_vec<double>("_soil_Vs ", _soil_Vs);
			// print_vec<double>("_soil_rho ", _soil_rho);
			// print_vec<double>("_soil_damp ", _soil_damp);
			// print_vec<double>("_soil_thick ", _soil_thick);
			// print_vec<double>("_soil_depth ", _soil_depth);
			return 0;
		}
	}
	if( new_depth < 0 ){
		cerr << " ERROR!!! in WaveField::add_compute_depth.  \n";
		cerr << " The depth cannot be smaller than zero! \n";
		cerr << " new_depth = " << new_depth << endl;
		return -1;
	}
	

	_soil_depth.push_back(INT_MAX); 
	
	auto it = lower_bound(_soil_depth.begin(), _soil_depth.end(), new_depth) ;
	int layer = it - _soil_depth.begin() ;
	double pre_depth = _soil_depth[layer-1] ; 
	
	_soil_depth.insert(it, new_depth);
	_soil_Vs.insert(_soil_Vs.begin()+layer, _soil_Vs[layer-1]);
	_soil_rho.insert(_soil_rho.begin()+layer, _soil_rho[layer-1]);
	_soil_damp.insert(_soil_damp.begin()+layer, _soil_damp[layer-1]);

	double whole_thick = _soil_thick[layer-1] ; 
	double thick1 = new_depth - pre_depth ; 
	double thick2 = pre_depth + whole_thick - new_depth; 
	
	_soil_thick[layer-1] = thick2 ; 
	_soil_thick.insert(_soil_thick.begin()+layer-1 , thick1);

	// _soil_thick[0] = _soil_depth[0] ;
	// for (int i = 1; i < (int)_soil_thick.size() ; ++i){
		// _soil_thick[i] = _soil_depth[i] - _soil_depth[i-1];
	// }

	_soil_depth.pop_back();

	// cout<<"--------------------------------"<<endl;
	// cout<<" new_depth = " << new_depth <<endl;
	// cout<< " pre_depth = " << pre_depth <<endl;
	// cout<< " whole_thick = " << whole_thick <<endl;
	// cout<< " thick1 = " << thick1 <<endl;
	// cout<< " thick2 = " << thick2 <<endl;
	// print_vec<double>("_soil_Vs ", _soil_Vs);
	// print_vec<double>("_soil_rho ", _soil_rho);
	// print_vec<double>("_soil_damp ", _soil_damp);
	// print_vec<double>("_soil_thick ", _soil_thick);
	// print_vec<double>("_soil_depth ", _soil_depth);

	return 1;
}


int WaveField::compute(){
	if (_soil_Vs.empty() || _soil_rho.empty() || _soil_damp.empty() || _soil_thick.empty()){
		cerr<<" ERROR!! in WaveField::set_soil_profile: at least one soil properties is empty" << endl;
		return -1;
	}
	if(!(
		 _soil_Vs.size() == _soil_rho.size() 
		 && _soil_rho.size() == _soil_damp.size()
		 && _soil_damp.size() == (_soil_thick.size()+1)
		)
	){
		cerr << " ERROR!!! in WaveField::compute. Soil properties Lengths mismatch!"<<"\n" ;
		cerr << " Vs.size()     = " << _soil_Vs.size() << "\n" ;
		cerr << " rho.size()    = " << _soil_rho.size() << "\n" ;
		cerr << " damp.size()   = " << _soil_damp.size() << "\n" ;
		cerr << " thicks.size() = " << _soil_thick.size() << "\n" ;
		cerr << " ERROR may happened in WaveField::update_compute_layers \n" ;
		return -1;
	}
	for (auto thick : _soil_thick){
		if(thick <= SOILD_THICKNESS_TOLERANCE){
			cerr<<"WaveField::compute() soil layer thickness is too small! \n " ; 
			cerr<<" thick = " << thick  << endl;
			cerr<<" STOP before cause numeric ERROR!. \n" ; 
			return -1;
		}
	}
	// cout<< "_max_freq = "<< _max_freq << endl;
	// cout<< "_Ntime = "<< _Ntime << endl;
	// print_vec<double>("_soil_Vs ", _soil_Vs);
	// print_vec<double>("_soil_rho ", _soil_rho);
	// print_vec<double>("_soil_damp ", _soil_damp);
	// print_vec<double>("_soil_thick ", _soil_thick);
	// print_vec<double>("_soil_depth ", _soil_depth);

	int layernum = _soil_Vs.size();
	auto field = new mat_complex(layernum, vec_complex(_Ntime)) ;

	this->wave_propagation(_max_freq, _Ntime
		, _soil_Vs
		, _soil_rho
		, _soil_damp
		, _soil_thick
		, field
	);

	auto DepthFreqAccAmp = mat_complex(layernum, vec_complex(_Ntime));
	_soil_acc = vector<vector<double>>(layernum, vector<double>(_Ntime));
	vec_complex input_acc_freq(_Ntime);
	for (int i = 0; i < _Ntime; ++i){
		input_acc_freq[i] = _input_acc[i] ;
	}
	this->fft(input_acc_freq);

	auto it = lower_bound(_soil_depth.begin(), _soil_depth.end(), _motion_depth);
	int motion_layer = it - _soil_depth.begin() ; 

	// convolution for shallow soils
	for (int i = 0; i < motion_layer ; ++i){
		for (int j = 0; j < _Ntime ; ++j){
			DepthFreqAccAmp[i][j] = (*field)[motion_layer][j] / (*field)[i][j] * input_acc_freq[j] ; 
		}
	}
	// deconvolution for deeper soils (include itself, do nothing)
	for (int i = motion_layer; i < layernum ; ++i){
		for (int j = 0; j < _Ntime ; ++j){
			DepthFreqAccAmp[i][j] = (*field)[i][j] / (*field)[motion_layer][j] * input_acc_freq[j] ; 
		}
	}
	delete field;
	for (int i = 0; i < layernum ; ++i){
		this->ifft( DepthFreqAccAmp[i] ) ; 
	}

	// change back to time domain
	for (int i = 0; i < layernum ; ++i){
		for (int j = 0; j < _Ntime; ++j){
			_soil_acc[i][j] = real(DepthFreqAccAmp[i][j]) ; 
		}
	}

	// calculate the velocity
	_soil_vel = vector<vector<double>>(layernum, vector<double>(_Ntime));
	_soil_dis = vector<vector<double>>(layernum, vector<double>(_Ntime));
	_soil_acc_freq = vector<vector<double>>(layernum);
	_soil_vel_freq = vector<vector<double>>(layernum);
	_soil_dis_freq = vector<vector<double>>(layernum);
	pair<vector<double>, vector<double>> p ;
	for (int i = 0; i < layernum; ++i){
		_soil_vel[i] = integrate_vel(_time_step, _soil_acc[i]);
		_soil_dis[i] = integrate_dis(_time_step, _soil_acc[i], _soil_vel[i]);
		p = time2freq(_time_step, _soil_acc[i]); _soil_acc_freq[i] = p.second; 
		p = time2freq(_time_step, _soil_vel[i]); _soil_vel_freq[i] = p.second; 
		p = time2freq(_time_step, _soil_dis[i]); _soil_dis_freq[i] = p.second; 
	}

	return 1;
}

/********************************************************************************
*	Purpose of deconvolution2bedrock() is for equivalent rock outcropping: 
*	1. which deconvolution to bedrock first.
*	2. and change _motion_depth to the bedrock.
*	3. and change _input_acc to the motion at bedrock.
*	4. Later, this->Compute() will do the convolution to soil layers.
*********************************************************************************/
int
WaveField::deconvolution2bedrock(){
	// **************************************************************************
	// 1. deconvolution to bedrock first.
	// **************************************************************************
	vector<double> rock_Vs = {_soil_Vs.back(), _soil_Vs.back(), _soil_Vs.back()}  ;
	vector<double> rock_rho = {_soil_rho.back(), _soil_rho.back(), _soil_rho.back()}  ;
	vector<double> rock_damp = {_soil_damp.back(), _soil_damp.back(), _soil_damp.back()}  ;
	vector<double> rock_thick(2,0.) ;

	if( fabs(_motion_depth) > machina_epsilon ){
		rock_thick = {_motion_depth , _soil_depth.back() - _motion_depth } ;
	}else{
		rock_thick = {_soil_depth.back() / 2. , _soil_depth.back() - _soil_depth.back() / 2. }; 
	}

	int layernum = rock_Vs.size() ;
	auto field = new mat_complex(layernum, vec_complex(_Ntime)) ;
	this->wave_propagation(_max_freq, _Ntime
		, rock_Vs
		, rock_rho
		, rock_damp
		, rock_thick
		, field
	);

	auto DepthFreqAccAmp = mat_complex(layernum, vec_complex(_Ntime));
	auto rock_acc = vector<vector<double>>(layernum, vector<double>(_Ntime));
	vec_complex input_acc_freq(_Ntime);
	for (int i = 0; i < _Ntime; ++i){
		input_acc_freq[i] = _input_acc[i] ;
	}
	this->fft(input_acc_freq);

	int motion_layer = 0 ;
	if( fabs(_motion_depth) > machina_epsilon ){
		motion_layer = 1;
	}else{
		motion_layer = 0;
	}
	int bedrock_layer = 2;
	// deconvolution to bedrock from the motion_layer.
	for (int j = 0; j < _Ntime ; ++j){
		DepthFreqAccAmp[bedrock_layer][j] = (*field)[bedrock_layer][j] / (*field)[motion_layer][j] * input_acc_freq[j] ; 
	}

	delete field;
	this->ifft( DepthFreqAccAmp[bedrock_layer] ) ; 

	// ************************************************************************** 
	// 2.  change _motion_depth to the bedrock.
	// **************************************************************************
	_motion_depth = _soil_depth.back();

	// **************************************************************************
	// 3.  change _input_acc to the motion at bedrock.
	// **************************************************************************
	for (int j = 0; j < _Ntime; ++j){
		_input_acc[j] = real(DepthFreqAccAmp[bedrock_layer][j]) ; 
	}

	return 1;
}




vector<double>  
WaveField::integrate_vel(double dt, vector<double> const& acc){
	if (acc.empty() || !dt){
		cerr<<"ERROR !! in WaveField::integrate_vel Input acc is empty() or dt is zero" << endl;
	}
	int len = acc.size();
	vector<double> vel(len, 0.) ; 
	for(int i = 1; i<len; ++i){
		vel[i] = vel[i-1] + dt * acc[i-1]; 
	}
	return vel;
}
		
vector<double>  
WaveField::integrate_dis(double dt, vector<double> const& acc, vector<double> const& vel){
	if (vel.empty() || !dt){
		cerr<<"ERROR !! in WaveField::integrate_dis Input vel is empty() or dt is zero" << endl;
	}
	int len = vel.size();
	vector<double> dis(len, 0.) ; 
	for(int i = 1; i<len; ++i){
		dis[i] = dis[i-1] + dt * vel[i-1] + 0.5 * acc[i-1] * dt * dt;
	}
	return dis;
}

pair<vector<double>, vector<double>>  
WaveField::time2freq(double dt, vector<double> const& series){
	if (series.size()<2 || !dt){
		cerr<<"ERROR !! in WaveField::time2freq Input series is too short or dt is zero" << endl;
	}
	int len = series.size();
	vec_complex series_complex(len);
	for (int i = 0; i < len; ++i){
		series_complex[i] = series[i] ;
	}
	this->fft(series_complex);
	vector<double> full_freqs(len);
	for (int i = 0; i < len; ++i){
		full_freqs[i] = fabs(series_complex[i])/len ;
	}
	int Nfreq = len/2 + 1 ; 
	vector<double> mag_of_half_freqs(full_freqs.begin(), full_freqs.begin()+Nfreq);
	for (int i = 1; i < Nfreq-1 ; ++i){
		mag_of_half_freqs[i] *= 2. ; 
	}

	_freqs.resize(Nfreq);
	double max_freq = 1. / dt ;
	std::iota(_freqs.begin(), _freqs.end(), 0);
	for(auto& item: _freqs){
		item = item / (Nfreq-1) * max_freq ; 
	}

	return make_pair(_freqs, mag_of_half_freqs) ;
}





vector<double> const& WaveField::get_times() {
	if (_times.empty()){
		cerr<<"ERROR!!! in WaveField::get_times Empty Results\n";
		cerr<<"Please set the motions, soil profile and compute first!"<<endl;
	}
	return _times;
}
vector<double> const& WaveField::get_acc_by_depth(double depth) {
	if (_soil_acc.empty()){
		cerr<<"ERROR!!! in WaveField::get_acc_by_depth Empty Results\n";
		cerr<<"Please set the motions, soil profile and compute first!"<<endl;
	}
	return _soil_acc[depth2layer(depth)] ;
}
vector<double> const& WaveField::get_vel_by_depth(double depth) {
	if (_soil_vel.empty()){
		cerr<<"ERROR!!! in WaveField::get_vel_by_depth Empty Results\n";
		cerr<<"Please set the motions, soil profile and compute first!"<<endl;
	}
	return _soil_vel[depth2layer(depth)] ;
}
vector<double> const& WaveField::get_dis_by_depth(double depth) {
	if (_soil_dis.empty()){
		cerr<<"ERROR!!! in WaveField::get_dis_by_depth Empty Results\n";
		cerr<<"Please set the motions, soil profile and compute first!"<<endl;
	}
	return _soil_dis[depth2layer(depth)] ;
}

vector<double> WaveField::get_acc_by_depth(double depth, size_t pos, size_t len) {
	if (_soil_acc.empty()){
		cerr<<"ERROR!!! in WaveField::get_acc_by_depth Empty Results\n";
		cerr<<"Please set the motions, soil profile and compute first!"<<endl;
	}
	if (pos + len > _soil_acc[0].size())	{
		cerr<<"ERROR!!! in WaveField::get_acc_by_depth exceed vector boundary \n";
		cerr<<" pos + len = " << pos + len << "\n";
		cerr<<"  _soil_acc[0].size() = " <<  _soil_acc[0].size() <<endl;
	}
	auto result =  _soil_acc[depth2layer(depth)] ; 
	return vector<double>(result.begin()+pos, result.begin()+pos+len) ;
}
vector<double> WaveField::get_vel_by_depth(double depth, size_t pos, size_t len) {
	if (_soil_vel.empty()){
		cerr<<"ERROR!!! in WaveField::get_vel_by_depth Empty Results\n";
		cerr<<"Please set the motions, soil profile and compute first!"<<endl;
	}
	if (pos + len > _soil_vel[0].size())	{
		cerr<<"ERROR!!! in WaveField::get_vel_by_depth exceed vector boundary \n";
		cerr<<" pos + len = " << pos + len << "\n";
		cerr<<"  _soil_vel[0].size() = " <<  _soil_vel[0].size() <<endl;
	}
	auto result =  _soil_vel[depth2layer(depth)] ; 
	return vector<double>(result.begin()+pos, result.begin()+pos+len) ;
}
vector<double> WaveField::get_dis_by_depth(double depth, size_t pos, size_t len) {
	if (_soil_dis.empty()){
		cerr<<"ERROR!!! in WaveField::get_dis_by_depth Empty Results\n";
		cerr<<"Please set the motions, soil profile and compute first!"<<endl;
	}
	if (pos + len > _soil_dis[0].size())	{
		cerr<<"ERROR!!! in WaveField::get_dis_by_depth exceed vector boundary \n";
		cerr<<" pos + len = " << pos + len << "\n";
		cerr<<"  _soil_dis[0].size() = " <<  _soil_dis[0].size() <<endl;
	}
	auto result =  _soil_dis[depth2layer(depth)] ; 
	return vector<double>(result.begin()+pos, result.begin()+pos+len) ;
}

vector<double> const& WaveField::get_freqs() {
	if (_freqs.empty()){
		cerr<<"ERROR!!! in WaveField::get_freqs Empty Results\n";
		cerr<<"Please set the motions, soil profile and compute first!"<<endl;
	}
	return _freqs;
}
vector<double> const& WaveField::get_acc_freq_by_depth(double depth) {
	if (_soil_acc_freq.empty()){
		cerr<<"ERROR!!! in WaveField::get_acc_freq_by_depth Empty Results\n";
		cerr<<"Please set the motions, soil profile and compute first!"<<endl;
	}
	return _soil_acc_freq[depth2layer(depth)] ;
}
vector<double> const& WaveField::get_vel_freq_by_depth(double depth) {
	if (_soil_vel_freq.empty()){
		cerr<<"ERROR!!! in WaveField::get_vel_freq_by_depth Empty Results\n";
		cerr<<"Please set the motions, soil profile and compute first!"<<endl;
	}
	return _soil_vel_freq[depth2layer(depth)] ;
}
vector<double> const& WaveField::get_dis_freq_by_depth(double depth) {
	if (_soil_dis_freq.empty()){
		cerr<<"ERROR!!! in WaveField::get_dis_freq_by_depth Empty Results\n";
		cerr<<"Please set the motions, soil profile and compute first!"<<endl;
	}
	return _soil_dis_freq[depth2layer(depth)] ;
}


int WaveField::get_Ntime() const{
	return _Ntime;
}
int WaveField::getTag() const{
	return _theTag;
}

void WaveField::setTag(int newTag){
	_theTag = newTag;
}

int WaveField::depth2layer(double depth){
	auto it = lower_bound(_soil_depth.begin(), _soil_depth.end(), depth);
	if ( *it - depth > SOILD_THICKNESS_TOLERANCE ){
		if(depth<0){
			cerr << " ERROR!!! in WaveField::depth2layer.  \n";
			cerr << " The request depth cannot be smaller than zero! \n";
			cerr << " Trying to get motion from depth = " << depth << endl;
			return 0;
		}

		auto pre = it - 1;
		if( depth - *pre <= SOILD_THICKNESS_TOLERANCE ){
			return pre - _soil_depth.begin() ; 
		}

		this->add_compute_depth(depth);
		this->compute();
		return depth2layer(depth);
	}
	return it - _soil_depth.begin() ; 
}


int WaveField::write_wave_at_depth(double depth, string prefix){
	auto acc = get_acc_by_depth(depth);
	auto vel = get_vel_by_depth(depth);
	auto dis = get_dis_by_depth(depth);

	this->write_series_to_file(_times, acc, prefix+to_short_string(depth)+"_acc.txt");
	this->write_series_to_file(_times, vel, prefix+to_short_string(depth)+"_vel.txt");
	this->write_series_to_file(_times, dis, prefix+to_short_string(depth)+"_dis.txt");

	return 1;
}

int WaveField::write_series_to_file(
	vector<double> const& c1, 
	vector<double> const& c2, 
	string filename) const
{
	ofstream ofs(filename);
	for (int i = 0; i < (int)c1.size() ; ++i){
		ofs << c1[i] << " \t " << c2[i] <<endl;
	}
	return 1;
}

// **********************************************************************************
// Vs          : vector of shear wave velocities of the soil layers, 
//                  top one first - bedrock last 
//                  [Vstop ... Vsbottom Vsbedrock]
// rho         : vector of the unit weight of the soil layers, 
//                  top one first - bedrock last
//                  [rhotop ... rhobottom rhobedrock]
// damp        : vector of the material damping of the soil layers, 
//                  top one first - bedrock last 
//                  [ksitop ... ksibottom ksibedrock]
// freq        : vector of the frequencies of interest
//                  [0:df:(N-1)*df] 
// layer_thick : vector with the thickness of the layers, top one first    
//                  [Htop ... Hbottom]
// % EXAMPLE:
// f0   = 0;
// Fs   = 100;  % (in Hz) Frequency sample
// NFFT = 2^12; % number of frequencies for discretization
// 
// input.Vs          = [200 300 2000];      % (m/s)
// input.rho         = [2000 2100 2400];    % (kgr/m3)
// input.damp        = [0.04 0.03 0.01];    % 
// input.layer_thick = [10 10];             % (m) ! no thickness for bedrock!
// **********************************************************************************
int WaveField::wave_propagation(double max_freq, int N_freq
	, std::vector<double> const& Vs
	, std::vector<double> const& rho
	, std::vector<double> const& damp
	, std::vector<double> const& layer_thick
	, mat_complex* disp
)
{

	int layernum = Vs.size();

	std::vector<double> freq(N_freq);
	std::iota(freq.begin(), freq.end(), 0);
	for(auto& item: freq){
		item = item / (N_freq - 1) * max_freq ; 
	}

	std::vector<double> omega(N_freq)  ;
	for (int i = 0; i < N_freq; ++i){
		omega[i] = 2 * WaveField::pi * freq[i] ; 
	}
	std::complex<double> i1(0, 1);

	vec_complex Vsstar(layernum);
	for (int i = 0; i < layernum; ++i){
		Vsstar[i] = Vs[i] * ( 1. + i1 * damp[i] );
	}

	vec_complex az(layernum-1);
	for (int i = 0; i < layernum-1 ; ++i){
		az[i] = (rho[i]/rho[i+1]) * (Vsstar[i]/Vsstar[i+1]) ; 
	}

	auto kstar = new mat_complex(layernum, vec_complex(N_freq)) ;
	auto up = new mat_complex(layernum, vec_complex(N_freq)) ;
	auto down = new mat_complex(layernum, vec_complex(N_freq)) ;
	

	for (int i = 0; i < layernum; ++i){
		for (int j = 0; j < N_freq; ++j){
			(*kstar)[i][j] = omega[j] / Vsstar[i] ; 

			if( i == 0 ){
				(*up)[i][j] = 0.5 * exp( i1 * (*kstar)[i][j] * layer_thick[i])
					+ 0.5 * exp( -i1 * (*kstar)[i][j] * layer_thick[i]) ;
				(*down)[i][j] = (*up)[i][j] ; 
			}else{
				(*up)[i][j] = 0.5 * (*up)[i-1][j] * (1. + az[i-1]) * exp( i1 * (*kstar)[i-1][j] * layer_thick[i-1])
					+ 0.5 * (*down)[i-1][j] * (1. - az[i-1]) * exp( -i1 * (*kstar)[i-1][j] * layer_thick[i-1]);
				(*down)[i][j] = 0.5 * (*up)[i-1][j] * (1. - az[i-1]) * exp( i1 * (*kstar)[i-1][j] * layer_thick[i-1])
					+ 0.5 * (*down)[i-1][j] * (1. + az[i-1]) * exp( -i1 * (*kstar)[i-1][j] * layer_thick[i-1]);
			}
			(*disp)[i][j] = (*up)[i][j] + (*down)[i][j] ; 
		}
	}

	int len ;
	if( N_freq & 1 ){
		len = ( N_freq -1 ) / 2 ;
	}else{
		len = ( N_freq -2 ) / 2 ;
	}
	std::vector<int> ia(len) ;
	std::vector<int> ib(len) ;
	for (int i = 0; i < len; ++i){
		ia[i] = 1 + i ;
		ib[i] = N_freq - i - 1 ;
	}

	for (int i = 0; i < layernum; ++i){
		for (int j = 0; j < len; ++j){
			// (*up)[i][ib[j]] = std::conj((*up)[i][ia[j]]) ; 
			// (*down)[i][ib[j]] = std::conj((*down)[i][ia[j]]) ; 
			(*disp)[i][ib[j]] = std::conj((*disp)[i][ia[j]]) ; 
		}
	}

	delete kstar ; 
	delete up;
	delete down;

	return 1;
}





// **************************************************************
// fft subroutines
// **************************************************************
void WaveField::fft(vector<complex<double> > &vec) {
	size_t n = vec.size();
	if (n == 0)
		return;
	else if ((n & (n - 1)) == 0)  // Is power of 2
		transformRadix2(vec);
	else  // More complicated algorithm for arbitrary sizes
		transformBluestein(vec);
}


void WaveField::ifft_inner(vector<complex<double> > &vec) {
	std::transform(vec.cbegin(), vec.cend(), vec.begin(),
		static_cast<complex<double> (*)(const complex<double> &)>(std::conj));
	fft(vec);
	std::transform(vec.cbegin(), vec.cend(), vec.begin(),
		static_cast<complex<double> (*)(const complex<double> &)>(std::conj));
}

void WaveField::ifft(vector<complex<double> > &vec){
	ifft_inner(vec);
	for_each(vec.begin(), vec.end(), [&](complex<double>& a){a /= (int)vec.size() ;} ) ; 
}

void WaveField::transformRadix2(vector<complex<double> > &vec) {
	// Length variables
	size_t n = vec.size();
	int levels = 0;  // Compute levels = floor(log2(n))
	for (size_t temp = n; temp > 1U; temp >>= 1)
		levels++;
	if (static_cast<size_t>(1U) << levels != n)
		throw "Length is not a power of 2";
	
	// Trignometric table
	vector<complex<double> > expTable(n / 2);
	for (size_t i = 0; i < n / 2; i++)
		expTable[i] = std::exp(complex<double>(0, -2 * M_PI * i / n));
	
	// Bit-reversed addressing permutation
	for (size_t i = 0; i < n; i++) {
		size_t j = reverseBits(i, levels);
		if (j > i)
			std::swap(vec[i], vec[j]);
	}
	
	// Cooley-Tukey decimation-in-time radix-2 FFT
	for (size_t size = 2; size <= n; size *= 2) {
		size_t halfsize = size / 2;
		size_t tablestep = n / size;
		for (size_t i = 0; i < n; i += size) {
			for (size_t j = i, k = 0; j < i + halfsize; j++, k += tablestep) {
				complex<double> temp = vec[j + halfsize] * expTable[k];
				vec[j + halfsize] = vec[j] - temp;
				vec[j] += temp;
			}
		}
		if (size == n)  // Prevent overflow in 'size *= 2'
			break;
	}
}


void WaveField::transformBluestein(vector<complex<double> > &vec) {
	// Find a power-of-2 convolution length m such that m >= n * 2 + 1
	size_t n = vec.size();
	size_t m = 1;
	while (m / 2 <= n) {
		if (m > SIZE_MAX / 2)
			throw "Vector too large";
		m *= 2;
	}
	
	// Trignometric table
	vector<complex<double> > expTable(n);
	for (size_t i = 0; i < n; i++) {
		unsigned long long temp = static_cast<unsigned long long>(i) * i;
		temp %= static_cast<unsigned long long>(n) * 2;
		double angle = M_PI * temp / n;
		// Less accurate alternative if long long is unavailable: double angle = M_PI * i * i / n;
		expTable[i] = std::exp(complex<double>(0, -angle));
	}
	
	// Temporary vectors and preprocessing
	vector<complex<double> > av(m);
	for (size_t i = 0; i < n; i++)
		av[i] = vec[i] * expTable[i];
	vector<complex<double> > bv(m);
	bv[0] = expTable[0];
	for (size_t i = 1; i < n; i++)
		bv[i] = bv[m - i] = std::conj(expTable[i]);
	
	// Convolution
	vector<complex<double> > cv(m);
	fft_convolve(av, bv, cv);
	
	// Postprocessing
	for (size_t i = 0; i < n; i++)
		vec[i] = cv[i] * expTable[i];
}


void WaveField::fft_convolve(
		const vector<complex<double> > &xvec,
		const vector<complex<double> > &yvec,
		vector<complex<double> > &outvec) {
	
	size_t n = xvec.size();
	if (n != yvec.size() || n != outvec.size())
		throw "Mismatched lengths";
	vector<complex<double> > xv(xvec);
	vector<complex<double> > yv(yvec);
	fft(xv);
	fft(yv);
	for (size_t i = 0; i < n; i++)
		xv[i] *= yv[i];
	ifft_inner(xv);
	for (size_t i = 0; i < n; i++)  // Scaling (because this FFT implementation omits it)
		outvec[i] = xv[i] / static_cast<double>(n);
}


size_t WaveField::reverseBits(size_t x, int n) {
	size_t result = 0;
	for (int i = 0; i < n; i++, x >>= 1)
		result = (result << 1) | (x & 1U);
	return result;
}



string WaveField::removeComments(string prgm)
{
    int n = prgm.length();
    string res;
 
    // Flags to indicate that single line and multpile line comments
    // have started or not.
    bool s_cmt = false;
    bool m_cmt = false;
 
 
    // Traverse the given program
    for (int i=0; i<n; i++)
    {
        // If single line comment flag is on, then check for end of it
        if (s_cmt == true && prgm[i] == '\n')
            s_cmt = false;
 
        // If multiple line comment is on, then check for end of it
        else if  (m_cmt == true && prgm[i] == '*' && prgm[i+1] == '/')
            m_cmt = false,  i++;
 
        // If this character is in a comment, ignore it
        else if (s_cmt || m_cmt)
            continue;
 
        // Check for beginning of comments and set the approproate flags
        else if (prgm[i] == '/' && prgm[i+1] == '/')
            s_cmt = true, i++;
        else if (prgm[i] == '/' && prgm[i+1] == '*')
            m_cmt = true,  i++;
 
        // If current character is a non-comment character, append it to res
        else  res += prgm[i];
    }
    return res;
}

