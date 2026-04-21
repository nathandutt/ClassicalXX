#include <fstream>
#include <sstream>
#include <iostream>

#include <cassert>
#include <stdexcept>
#include <utility> 
#include <stdlib.h>
#include <time.h>

#include <complex>
#include <vector>
#include <string>

#include <cmath>
#include <algorithm>

using namespace std;

// Class to store simulation and initial data parameters
class Config {
    // Simulation Parameters
    public:
	double x_range;
	double dx;
	double dt;
	double smoothing;
	int domain_size;
	double t_max;
	double save_every;

	// Initial Data Parameters
	//This is for evolving a step
	double theta_l;
	double theta_r;
	double v_l;
	double v_r;

	Config(){
	}

	//Parser for simu.params
	Config(std::string filename) {

	    //Check if config file is valid
	    std::ifstream in(filename);
	    if (!in) {
		throw std::runtime_error("Could not open config file: " + filename);
	    }

	    std::string line;
	    while (std::getline(in, line)) {

		std::stringstream ss(line);
		std::string key;
		double value;

		// Skip comments or empty lines
		if (line.empty() || line[0] == '#') {
		    continue;
		}

		// Read the key-value pair
		if (ss >> key >> value) {
		    if (key == "x_range") x_range = value;
		    else if (key == "dx") dx = value;
		    else if (key == "dt") dt = value;
		    else if (key == "domain_size") domain_size = static_cast<int>(value);
		    else if (key == "t_max") t_max = value;
		    else if (key == "smoothing") smoothing = value;
		    else if (key == "save_every") save_every = value;
		    else if (key == "theta_l") theta_l = M_PI*value;
		    else if (key == "theta_r") theta_r = M_PI*value;
		    else if (key == "v_l") v_l = value;
		    else if (key == "v_r") v_r = value;

		    else {
			std::cerr << "Warning: Unknown config key '" << key << "'\n";
		    }
		}
	    }
	    //Get domain size from config.
	    domain_size = static_cast<int>(2 * x_range / dx);
	}
};


//Most basic object : single spin 
struct Spin{

    //Cartesian components
    double x, y, z;

    Spin(){};

    //Constructor given two angles
    Spin(double theta, double phi){
	x = std::sin(theta)*std::cos(phi);
	y = std::sin(theta)*std::sin(phi);
	z = std::cos(theta);
    }

    //Basic constructor
    Spin(double x, double y, double z) :
	x(x), y(y), z(z) {}

    //Arithmetic operations
    Spin operator-(const Spin& other) const{
	return Spin(x - other.x, y-other.y, z - other.z);
    }
    Spin operator+(const Spin& other) const{
	return Spin(x + other.x, y+other.y, z + other.z);
    }

    //Vector product
    Spin operator*(const Spin& other) const{
	double n_x = y*other.z - z*other.y;
	double n_y = z*other.x - x*other.z;
	double n_z = x*other.y - y*other.x;
	return Spin(n_x, n_y, n_z);
    }
};

//Multiplication by scalar
Spin operator*(const double l, const Spin& s){
    return Spin(l*s.x, l*s.y, l*s.z);
}

//Whole object : spin chain
//Wrapper for a vector of spins
struct SpinChain{

    std::vector<Spin> spins;

    //Constructors
    SpinChain(){};

    SpinChain(std::vector<Spin> spins) :
	spins(spins){} 

    //Construct from two arrays of angles.
    SpinChain(const std::vector<double>& thetas, const std::vector<double>& phis){
	assert((thetas.size() == phis.size()));
	int N = thetas.size();
	spins = std::vector<Spin>{};
	spins.reserve(N);

	for(int i = 0; i < N; i++){
	    Spin n_Spin = Spin(thetas[i], phis[i]);
	    spins.emplace_back(n_Spin);
	}
    }

    //Arithmetic operations
    //Note: for some speed gains it could be preferable to work
    //in place and not copy. But this is clearer.
    SpinChain operator-(const SpinChain& other) const{
	std::vector<Spin> chain;
	int N = spins.size();
	chain.reserve(N);

	for(int i =0; i < N; i++){
	    Spin n_s = spins[i] - other.spins[i];
	    chain.emplace_back(n_s);
	}

	return SpinChain(chain);
    }

    SpinChain operator+(const SpinChain& other) const{
	std::vector<Spin> chain;
	int N = spins.size();
	chain.reserve(N);

	for(int i =0; i < N; i++){
	    Spin n_s = spins[i] + other.spins[i];
	    chain.emplace_back(n_s);
	}

	return SpinChain(chain);
    }
};

SpinChain operator*(const double l, const SpinChain& sc){
    std::vector<Spin> chain;
    int N = sc.spins.size();
    chain.reserve(N);

    for(int i =0; i < N; i++){
	Spin n_s = l*sc.spins[i];
	chain.emplace_back(n_s);
    }
    return SpinChain(chain);
}

double IdxToCoord(const Config& config, int idx) {
    return -config.x_range + idx * config.dx;
}

// ---------------------------DIFFERENT INITAL CONDITION SHAPES-------------------------------------------

double StepTheta(const Config& config, double x){
    double width = config.smoothing*config.dx;
    return 0.5*(config.theta_l + config.theta_r) + 0.5*(config.theta_r - config.theta_l)*std::tanh(x/width);
}
double StepPhi(const Config& config, double x){
    double width = config.smoothing*config.dx;
    double v_l = config.v_l;
    double v_r = config.v_r;
    return 0.5 * (v_l + v_r)*x + 0.5*(v_r-v_l)*width * std::log(2.*std::cosh(x/width));
}

double Slab(const double x, const double width, const double center, const double smoothing){
    //Outputs a slab with fixed center and width/2. and smoothing, of amplitude 1 
    return 0.5*tanh((x-center + width/2.)/smoothing) - 0.5*tanh((x-center - width/2.)/smoothing);
}

double IntSlab(const double x, const double width, const double center, const double smoothing){
    //Outputs primitive of a slab:
    return 0.5*smoothing*(log(cosh((x-center + width/2.)/smoothing)) - log(cosh((x-center-width/2.)/smoothing)));
}

double width1 = 30.;
double width2 = 50.;
double separation = 100.;
double DoubleSlabTheta(const Config& config, const double x){
    double background_theta = 0.;
    double res = background_theta + (config.theta_l - background_theta)*Slab(x, width1, -1.*separation/2., config.smoothing*config.dx)
	+(config.theta_r - background_theta)*Slab(x, width2, separation/2., config.smoothing*config.dx);
    return res;
}
double DoubleSlabPhi(const Config& config, const double x){
    double background_v = 0.3;
    double res = background_v*x + (config.v_l - background_v)*IntSlab(x, width1, -1.*separation/2., config.smoothing*config.dx)
	+(config.v_r - background_v)*IntSlab(x, width2, separation/2., config.smoothing*config.dx);
    return res;
}
//---------------------------------------------------------------------------------------------------------

//To actually initialize thetas, phi, with given shape
std::vector<double> InitialTheta(const Config& config){
    std::vector<double> density;
    density.reserve(config.domain_size);

    for(int i =0; i<config.domain_size; i++){

	/* UNCOMMENT FOR SLAB INITIAL CONFIG
	 * 
	 * double x = IdxToCoord(config, i);
	 * double dens = DoubleSlabTheta(config, x);
	 *
	 */

	//--CONSTANT INITIAL DATA WITH NOISE---
	double dens = M_PI * (1./2.);
	dens += 0.1*drand48();
	//-------------------------------------
	
	density.emplace_back(dens);
    } 

    return density;
}

std::vector<double> InitialPhi(const Config& config){
    std::vector<double> u;
    u.reserve(config.domain_size);

    for(int i =0; i<config.domain_size; i++){

	/* UNCOMMENT FOR SLAB INITIAL CONFIG
	 *
	 * double x = IdxToCoord(config, i);
	 * double u_i = DoubleSlabPhi(config, x);
	 *
	 */

	double u_i = 0.;
	u.emplace_back(u_i);
    } 
    return u;
}

//This takes laplacian enforcing zero at points within int limit of boundaries
int limit = 5;
SpinChain DoubleDeriv(const Config& config, const SpinChain& chain) {
    std::vector<Spin> gradient(config.domain_size);

    const std::vector<Spin>& spins = chain.spins;

    for (int i = 0; i < config.domain_size; i++) {
	if(i<limit){
	    gradient[i] = Spin(0., 0., 0.);
	    continue;
	}
	if(i+limit > config.domain_size-1){
	    gradient[i] = Spin(0., 0., 0.);
	    continue;
	}

	gradient[i] = (1./ (config.dx * config.dx))*(spins[i-1] - 2.*spins[i] + spins[i+1]);
    }

    return SpinChain(gradient);
}

//This gives laplacian enforcing periodic boundary conditions
SpinChain DoubleDerivPeriodic(const Config& config, const SpinChain& chain) {
    std::vector<Spin> gradient(config.domain_size);
    const std::vector<Spin>& spins = chain.spins;

    int N = config.domain_size; 
    for (int i = 0; i < N; i++) {
	if((i>0) && (i<(N-1)))
	    gradient[i] = (1./ (config.dx * config.dx))*(spins[i-1] - 2.*spins[i] + spins[i+1]);

	else if(i==0)
	    gradient[0] = 1./(config.dx * config.dx)*(spins[N-1] -2.*spins[i] + spins[i+1]);

	else
	    gradient[N-1] = 1./(config.dx*config.dx)*(spins[N-2] - 2.*spins[N-1]+spins[0]);
    }
    return SpinChain(gradient);
}

//Calculate force term at each site
SpinChain F(const Config& config, const SpinChain& chain){
    //Implementing PDE: \partial_t S = S \times (\partial_x^2 S_{perp} - S_z e_z)
    vector<Spin> result;
    result.reserve(config.domain_size);

    //Taking periodic BC to observe coarsening
    SpinChain deriv = DoubleDerivPeriodic(config, chain);

    const std::vector<Spin>& S_xx = deriv.spins;
    const std::vector<Spin>& spins = chain.spins;

    for (int i = 0; i < config.domain_size; i++) {
	Spin z_proj = Spin{0., 0., spins[i].z};

	/* Uncomment for classical XX instead of easy axis LL
	 *
	 * Spin S_perp = Spin(S_xx[i].x, S_xx[i].y, 0.);
	 * Spin B_eff = S_perp - z_proj;
         */

	//--LL easy axis B_eff ---
	Spin B_eff = S_xx[i]  + z_proj;
	//------------------------

	Spin f_i = (spins[i]*B_eff);
	result.emplace_back(f_i);
    }
    SpinChain spin_chain_res(result);
    return spin_chain_res;
}

const complex I(0., 1.);
// Psi class for evolving the state using Runge-Kutta
class Psi {
    double time;
    SpinChain chain;
    Config config;

    public:
    //Construct from config
    Psi(Config config_) {
	time = 0.;
	config = config_;

	//Get intial thetas and phis
	std::vector<double> theta = InitialTheta(config_);
	std::vector<double> phi = InitialPhi(config_);

	//Construct inital spin chain from angles
	chain = SpinChain(theta, phi);
    };

    //Evolve chain with RK4
    //Note: maybe try a symplectic solver, and avoid so many copies
    void EvolveRK(){
	SpinChain k_1 = F(config, chain);
	SpinChain k_2 = F(config, chain + (config.dt/2.)*k_1);
	SpinChain k_3 = F(config, chain + (config.dt/2.)*k_2);
	SpinChain k_4 = F(config, chain + (config.dt)*k_3);
	SpinChain evolved = chain + (config.dt/6.)*(k_1 + 2.*k_2 + 2.*k_3 + k_4);
	chain = evolved;
	time = time + config.dt;
    }

    //Write spins to file
    //Format is 3 lines per timestep
    // t spin1_x spin2_x ..
    // t spin1_y spin2_y ..
    // t spin1_z spin2_z ...
    void Write(std::fstream& file) const {

	file << time;
	for (const Spin& c : chain.spins) {
	    file << "," << c.x;
	}
	file << "\n";
	file << time;
	for (const Spin& c : chain.spins){ 
	    file << "," << c.y;
	}
	file << "\n";
	file << time;
	for (const Spin& c : chain.spins){ 
	    file << "," << c.z;
	}
	file << "\n";

    }
};

// Main function to run the simulation
int main() {
    //Initialize seed for random initial profile (For studying coarsening)
    srand48(time(NULL));

    //Read config
    Config config("simu.params");

    //Calculate evolution timesteps
    int write_step = static_cast<int>(config.save_every/config.dt);
    int timesteps = static_cast<int>(config.t_max / config.dt);

    //Intialize RK4 evolver
    Psi state(config);

    //Open output file
    fstream file;
    string filename = "evolution.csv";
    file.open(filename, std::ios::out);

    // Evolve system and save data to file
    for (int t = 0; t < timesteps; t++) {

	//Only write if step is good
	if(t%write_step == 0){
	    cout << "t = " << config.dt*t << "\n";
	    state.Write(file);
	}
	//Evolve
	state.EvolveRK();
    }

    return 0;
}

