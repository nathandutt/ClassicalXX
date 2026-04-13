#include <cassert>
#include <fstream>
#include <complex>
#include <vector>
#include <string>
#include <sstream>
#include <stdexcept>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <utility> 
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
    double theta_l;
    double theta_r;
    double v_l;
    double v_r;

    Config(){
    }

    Config(std::string filename) {

	//Check if config file is valid
        std::ifstream in(filename);
        if (!in) {
            throw std::runtime_error("Could not open config file: " + filename);
        }

        std::string line;
        while (std::getline(in, line)) {
            // Remove leading and trailing whitespaces
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

        // Derived values
        domain_size = static_cast<int>(2 * x_range / dx);
    }
};

struct Spin{

    double x, y, z;

    Spin(){};
    Spin(double theta, double phi){
        x = std::sin(theta)*std::cos(phi);
        y = std::sin(theta)*std::sin(phi);
        z = std::cos(theta);
    }
    Spin(double x_, double y_, double z_) :
	x(x_), y(y_), z(z_) {}
    Spin operator-(const Spin& other) const{
        return Spin(x - other.x, y-other.y, z - other.z);
    }
    Spin operator+(const Spin& other) const{
        return Spin(x + other.x, y+other.y, z + other.z);
    }
    Spin operator*(const Spin& other) const{
        auto n_x = y*other.z - z*other.y;
        auto n_y = z*other.x - x*other.z;
        auto n_z = x*other.y - y*other.x;
        return Spin(n_x, n_y, n_z);
    }
};

Spin operator*(const double l, const Spin& s){
    return Spin(l*s.x, l*s.y, l*s.z);
}

struct SpinChain{
    std::vector<Spin> spins;

    //Constructors
    SpinChain(){};

    SpinChain(std::vector<Spin> spins_) :
	spins(spins_){} 

    SpinChain(const std::vector<double>& thetas, const std::vector<double>& phis){
	assert((thetas.size() == phis.size()));
        spins = std::vector<Spin>{};
	spins.reserve(thetas.size());
        for(int i = 0; i < thetas.size(); i++){
            auto n_Spin = Spin(thetas[i], phis[i]);
            spins.emplace_back(n_Spin);
        }
    }

    //Arithmetic operation overloads
    SpinChain operator-(const SpinChain& other) const{
        auto chain = std::vector<Spin>{};
        for(int i =0; i < spins.size(); i++){
            auto n_s = spins[i] - other.spins[i];
            chain.emplace_back(n_s);
        }
        return SpinChain(chain);
    }
    SpinChain operator+(const SpinChain& other) const{
        auto chain = std::vector<Spin>{};
        for(int i =0; i < spins.size(); i++){
            auto n_s = spins[i] + other.spins[i];
            chain.emplace_back(n_s);
        }
        return SpinChain(chain);
    }


};
SpinChain operator*(const double l, const SpinChain& sc){
    auto chain = std::vector<Spin>{};
    chain.reserve(sc.spins.size());
    for(int i =0; i < sc.spins.size(); i++){
	auto n_s = l*sc.spins[i];
	chain.emplace_back(n_s);
    }
    return SpinChain(chain);
}

auto IdxToCoord(const Config& config, int idx) {
    return -config.x_range + idx * config.dx;
}

// Initial conditions shapes
auto StepTheta(const Config& config, double x){
    auto width = config.smoothing*config.dx;
    return 0.5*(config.theta_l + config.theta_r) + 0.5*(config.theta_r - config.theta_l)*std::tanh(x/width);
}
auto StepPhi(const Config& config, double x){
    auto width = config.smoothing*config.dx;
    auto v_l = config.v_l;
    auto v_r = config.v_r;
    return 0.5 * (v_l + v_r)*x + 0.5*(v_r-v_l)*width * std::log(2.*std::cosh(x/width));
}

auto GaussianTheta(const Config& config, double x){

    return std::acos(0.45*std::exp(-x*x/100.));

}


//To actually initialize thetas, phi, with given shape
auto InitialTheta(const Config& config){
    auto density = vector<double>{};
    for(int i =0; i<config.domain_size; i++){
        auto x = IdxToCoord(config, i);
        auto dens = StepTheta(config, x);
        density.emplace_back(dens);
    } 
    return density;
}

auto InitialPhi(const Config& config){
    auto u = vector<double>{};
    for(int i =0; i<config.domain_size; i++){
        auto x = IdxToCoord(config, i);
        auto u_i = StepPhi(config, x);
        u.emplace_back(u_i);
    } 
    return u;
}

auto limit = 5;
auto DoubleDeriv(const Config& config, const SpinChain& chain) {
    auto gradient = std::vector<Spin>(config.domain_size);
    const auto& spins = chain.spins;
    
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

auto F(const Config& config, const SpinChain& chain){
    //Implementing PDE: \partial_t S = S \times (\partial_x^2 S_{perp} - S_z e_z)
    auto result = vector<Spin>{};
    result.reserve(config.domain_size);

    auto deriv = DoubleDeriv(config, chain);
    const auto& S_xx = deriv.spins;
    const auto& spins = chain.spins;
    for (int i = 0; i < config.domain_size; i++) {
        auto z_proj = Spin{0., 0., spins[i].z};
        auto S_perp = Spin(S_xx[i].x, S_xx[i].y, 0.);
        auto B_eff = S_perp - z_proj;
        auto f_i = (spins[i]*B_eff);
        result.emplace_back(f_i);
    }
    auto spin_chain_res = SpinChain(result);
    return spin_chain_res;
}

const auto I = complex(0., 1.);
// Psi class for evolving the state using Runge-Kutta
class Psi {
    double time;
    SpinChain chain;
    Config config;

public:
    Psi(Config config_) {
        time = 0.;
        config = config_;
        auto theta = InitialTheta(config_);
        auto phi = InitialPhi(config_);
        chain = SpinChain(theta, phi);
    };
    void EvolveRK(){
        auto k_1 = F(config, chain);
        auto k_2 = F(config, chain + (config.dt/2.)*k_1);
        auto k_3 = F(config, chain + (config.dt/2.)*k_2);
        auto k_4 = F(config, chain + (config.dt)*k_3);
        auto evolved = chain + (config.dt/6.)*(k_1 + 2.*k_2 + 2.*k_3 + k_4);
        chain = evolved;
        time = time + config.dt;
    }

      void Write(std::fstream& file) const {

        file << time;
        for (const auto& c : chain.spins) {
            file << "," << c.x;
        }
        file << "\n";
        file << time;
        for (const auto& c : chain.spins){ 
            file << "," << c.y;
        }
        file << "\n";
        file << time;
        for (const auto& c : chain.spins){ 
            file << "," << c.z;
        }
        file << "\n";

    }
};

// Main function to run the simulation
int main() {
    string filename = "evolution.csv";
    std::fstream file;
    // Load configuration from the file
    Config config("config.txt");
    std::cout << config.v_r << " "<< config.v_l << " Yepa \n";
    auto step = static_cast<int>(config.save_every/config.dt);
    // Calculate timesteps
    auto timesteps = static_cast<int>(config.t_max / config.dt);
    auto state = Psi(config);
    // Evolve system and save data to file
    file.open(filename, std::ios::out);
    for (int t = 0; t < timesteps; t++) {
        if(t%step == 0){
	cout << "t = " << config.dt*t << "\n";
	state.Write(file);
	}
	state.EvolveRK();
    }

    return 0;
}

