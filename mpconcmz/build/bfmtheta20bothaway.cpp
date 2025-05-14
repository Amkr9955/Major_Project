#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <random>
#include <numeric> 
#include <set>     
#include <utility> 
#include <algorithm> 


using namespace std;

// Constants
const int LATTICE_SIZE = 500;
const int CHAIN_LENGTH = 20;  
const double MIN_BOND_LENGTH_SQUARED = 4.00;   
const double MAX_BOND_LENGTH_SQUARED = 13.00;  
const double WALL_DISTANCE = 2.00;
const double WALL_DISTANCE_SQUARED = WALL_DISTANCE * WALL_DISTANCE; 
const double NB_CONTACT_DIST_L1 = 2.0;

// Temperature range
const double T_START = 0.01;
const double T_END = 2.00;
const double T_STEP = 0.01;
const double k = 1.00; // Boltzmann constant (set to 1 for reduced units)

// Allowed bond lengths for 4 site lattice model
const double allowed_bond_lengths_squared[] = {4.00, 5.00, 8.00, 9.00, 10.00,13.00};

const int Z_MIN = 2;      // Z position must be >=2

struct Monomer {
    int x, y, z;
    // int nucleotide; // 0-A, 1-T, 2-G, 3-C
};

ranlux48 rng(random_device{}());
uniform_real_distribution<double> dist(0.0, 1.0);

// Function prototypes
double count_contacts(const vector<Monomer>& chain1, const vector<Monomer>& chain2);

// === Modified file input with nucleotide info ===
void initialize_chains_from_file(vector<Monomer>& chain1, vector<Monomer>& chain2, const string& filename) {
    ifstream infile(filename);
    if (!infile) {
        cerr << "Error: Cannot open file " << filename << endl;
        exit(1);
    }

    string line;
    for (int i = 0; i < CHAIN_LENGTH;) {
        getline(infile, line);
        if (line.empty() || line[0] == '#') continue;
        stringstream ss(line);
        ss >> chain1[i].x >> chain1[i].y >> chain1[i].z ;
        i++;
    }
    for (int i = 0; i < CHAIN_LENGTH;) {
        getline(infile, line);
        if (line.empty() || line[0] == '#') continue;
        stringstream ss(line);
        ss >> chain2[i].x >> chain2[i].y >> chain2[i].z ;
        i++;
    }
    infile.close();
}

// Softened potentials
double chain_chain_energy(double dist_sq) {
    const double sigma_sq = 4.0;  // WALL_DISTANCE_SQUARED
    const double barrier_sq = 10.0; // MAX_BOND_LENGTH_SQUARED
    if (dist_sq > barrier_sq) return 0.0; // Zero-energy cutoff
    return -1.0 * (1.0 - (dist_sq - sigma_sq) / (barrier_sq - sigma_sq)); // Quadratic decay
}

double chain_wall_energy(double z) {
    const double z_peak = 2.0;
    const double z_cutoff = 3.0; // Zero-energy cutoff
    if (z > z_cutoff) return 0.0;
    return -1.0 * (1.0 - pow((z - z_peak) / (z_cutoff - z_peak), 2)); // Quadratic decay
}

// Function to calculate squared distance between two monomers in 3D
double distance_squared(const Monomer& m1, const Monomer& m2) {
    double mx=m1.x-m2.x;
    double my=m1.y-m2.y;
    double mz=m1.z-m2.z;
    return (mx*mx) + (my*my) + (mz*mz);
}

// Function to check if a bond length is allowed
bool is_valid_bond_length(double dist_sq) {
    for (double allowed_dist_sq : allowed_bond_lengths_squared) {
        if (fabs(dist_sq - allowed_dist_sq) < 1e-6) { // Check if distance is one of the allowed values, 1e-6 for overflow control
            return true;
        }
    }
    return false;
}

// Function to check if a proposed move is valid (self-avoidance and bond length constraints)
bool is_valid_move(const vector<Monomer>& chain1, const vector<Monomer>& chain2, int chain_index, int monomer_index, const Monomer& new_pos) {
    const vector<Monomer>& current_chain = (chain_index == 0) ? chain1 : chain2;
    const vector<Monomer>& other_chain = (chain_index == 0) ? chain2 : chain1;

    if (new_pos.z < Z_MIN- 1e-6) {
        return false;  // z cannot be less than Z_MIN
    }

    // Check bond length constraints within the current chain
    if ((monomer_index > 0) && !is_valid_bond_length(distance_squared(new_pos, current_chain[monomer_index - 1]))) return false;
    if ((monomer_index < (CHAIN_LENGTH - 1)) && !is_valid_bond_length(distance_squared(new_pos, current_chain[monomer_index + 1]))) return false;
    
    // Check self-avoidance within the current chain
    for (int i = 0; i < CHAIN_LENGTH; i++) {
        if (i != monomer_index && distance_squared(new_pos, current_chain[i]) < MIN_BOND_LENGTH_SQUARED) {
            return false;  // Overlap detected
        }
    }

    // Check self avoidance with other chain
    for (int i = 0; i < CHAIN_LENGTH; i++) {
        if (distance_squared(new_pos, other_chain[i]) < MIN_BOND_LENGTH_SQUARED) {
            return false;  // Overlap detected
        }
    }

    return true;
}

// Function to calculate the energy change for a monomer move
double calculate_energy_change(const Monomer& current, const Monomer& other, const Monomer& new_pos) {
    double dE = 0.0;
    dE += chain_wall_energy(new_pos.z) - chain_wall_energy(current.z);
    dE += chain_chain_energy(distance_squared(new_pos, other)) 
          - chain_chain_energy(distance_squared(current, other));
    return dE;
}

double calculate_total_energy(const vector<Monomer>& chain1, const vector<Monomer>& chain2) {
    double energy = 0.0;
    for (int i = 0; i < CHAIN_LENGTH; ++i) {
        energy += chain_wall_energy(chain1[i].z) + chain_wall_energy(chain2[i].z);
        energy += chain_chain_energy(distance_squared(chain1[i], chain2[i]));
    }
    return energy;
}

// Function to count the number of contacts between the two chains
double count_contacts(const vector<Monomer>& chain1, const vector<Monomer>& chain2) {
    double contacts = 0;
    for (int i = 0; i < CHAIN_LENGTH; ++i) {
        double dist_sq = distance_squared(chain1[i], chain2[i]);
        if (dist_sq < 10.0000) {
            contacts+=1.00;
        }
    }
    return contacts;
}

double count_contacts_wall(const vector<Monomer>& chain) {
    double contacts = 0;
    for (int i = 0; i < CHAIN_LENGTH; ++i) {
        if (chain[i].z < 3.00) {
            contacts+=1.00;
        }
    }
    return contacts;
}

// average energy and contact calculations
void average(vector<Monomer>& chain1, vector<Monomer>& chain2,double& total_energy,double& E,double& Esquared,double& E4,double& c1,double& c1squared,double& c2,double& c2squared,double& c3,double& c3squared,double& BETA, long long int& sweep,long long int& M) {
    // Select random chain and monomer (excluding last monomer)
    int chain_index = uniform_int_distribution<int>(0, 1)(rng);  // 0 for chain1, 1 for chain2
    int monomer_index = uniform_int_distribution<int>(0, CHAIN_LENGTH - 2)(rng);  // Exclude last monomer

    // Get references to the selected chain and monomer
    vector<Monomer>& selected_chain = (chain_index == 0) ? chain1 : chain2;
    Monomer current_monomer = selected_chain[monomer_index];
    Monomer other_chain_monomer = (chain_index == 0) ? chain2[monomer_index] : chain1[monomer_index];

    // Attempt to move the monomer in a random direction by 1 unit
    Monomer trial_pos = current_monomer;  // Initialize trial position
    int direction = uniform_int_distribution<int>(0, 3)(rng); // Random direction (0: +x, 1: -x, 2: +y, 3: -y, 4: +z, 5: -z)
    
    switch (direction) {
        //   case 0: trial_pos.x += 1; break;  // +x
        //   case 1: trial_pos.x -= 1; break;  // -x
          case 0: trial_pos.y += 1; break;  // +y
          case 1: trial_pos.y -= 1; break;  // -y
          case 2: trial_pos.z += 1; break;  // +z
          case 3: trial_pos.z -= 1; break;  // -z
    }
    
    // Check if the move satisfies SAW and bond length conditions
    if (is_valid_move(chain1, chain2, chain_index, monomer_index, trial_pos)) {
        double delta_energy = calculate_energy_change(current_monomer, other_chain_monomer, trial_pos);

        // Metropolis acceptance criterion
        if (delta_energy <= 1e-6) {
            selected_chain[monomer_index] = trial_pos;
            total_energy += delta_energy;
        } else {
            double xi = dist(rng);
            double p = exp(-BETA * delta_energy);
            if (xi < p) {
                selected_chain[monomer_index] = trial_pos;
                total_energy += delta_energy;
            }
        }
    }
     double co1=0.00,co2=0.00,co3=0.00;
    if (sweep % 1000 == 0) {
                E += total_energy;
                Esquared += (total_energy * total_energy);
                E4 += (total_energy * total_energy * total_energy * total_energy);
                co1 = count_contacts_wall(chain1);
                c1+=co1;
                c1squared+=(co1*co1);
                co2 = count_contacts_wall(chain2);
                c2+=co2;
                c2squared+=(co2*co2);
                co3 = count_contacts(chain1, chain2);
                c3+=co3;
                c3squared+=(co3*co3);
                ++M;
          }
           
           
}

// Function to save results to a CSV file
void save_results_to_csv(const vector<double>& T_values, const vector<double>& E_, const vector<double>& Esquared_,const vector<double>& Binder, const vector<double>& c_1,const vector<double>& c_1_fluct,const vector<double>& c_2,const vector<double>& c_2_fluct,const vector<double>& c_3,const vector<double>& c_3_fluct, const vector<double>& Cv_, const string& filename) {
    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }

    file << "Temperature,E,Esquared,Binder,C1,C1fluc,C2,C2fluc,C3,C3fluc,Cv\n";
    for (size_t i = 0; i < T_values.size(); ++i) {
        file << T_values[i] << "," << E_[i] << "," << Esquared_[i] << "," << Binder[i] << "," << c_1[i] << "," << c_1_fluct[i] << "," << c_2[i] << "," << c_2_fluct[i] << "," << c_3[i] << "," << c_3_fluct[i] << "," << Cv_[i] << "\n";
    }
    file.close();
}


// Main simulation loop
int main() {
    vector<double> E_, Esquared_,Binder, c_1,c_1_fluct,c_2,c_2_fluct,c_3,c_3_fluct, Cv_, T_values;
    vector<Monomer> chain1(CHAIN_LENGTH);
    vector<Monomer> chain2(CHAIN_LENGTH);
    string input_file = "polymer_input20.txt";
    initialize_chains_from_file(chain1, chain2, input_file);

    double total_energy = calculate_total_energy(chain1, chain2);
    cout << "Initial Energy (Before Simulation): " << total_energy << endl;

    double E = 0.00, Esquared = 0.00,E4=0.00, c1 = 0.00,c1squa = 0.00,c2 = 0.00,c2squa = 0.00 ,c3 = 0.00,c3squa = 0.00,Cv=0.00,c1fluc=0.00,c2fluc=0.00,c3fluc=0.00,B=0.00;

    cout<<"Initial contacts : "<< " c1 : "<< count_contacts_wall(chain1)<<" c2 :"<< count_contacts_wall(chain2) << " c3 :" << count_contacts(chain1,chain2)<<endl;

   for (double T = T_START; T <= T_END + 0.0001; T += T_STEP) {
        total_energy = calculate_total_energy(chain1, chain2);
        double BETA = 1.0 / (k * T);
        T_values.push_back(T);
        long long int monte_carlo_step = 2*CHAIN_LENGTH;
        long long int monte_carlo_sweep1 = 1e6 ;
        long long int M=0;

            for (long long int sweep1 = 1; sweep1 <= monte_carlo_sweep1; sweep1++) {
                for (int step = 1; step <= monte_carlo_step; step++) {
                    average(chain1,chain2,total_energy,E,Esquared,E4,c1,c1squa,c2,c2squa,c3,c3squa,BETA,sweep1,M);
                }
                
            }
       
        E = 0.00, Esquared = 0.00,E4 = 0.00, c1 = 0.00,c1squa = 0.00,c2 = 0.00,c2squa = 0.00 ,c3 = 0.00,c3squa = 0.00,Cv = 0.00,c1fluc = 0.00,c2fluc = 0.00,c3fluc = 0.00;
        M=0;

        for (long long int sweep2 = 1; sweep2 <= monte_carlo_sweep1; sweep2++){
            for (int step = 1; step <= monte_carlo_step; step++) {
                average(chain1,chain2,total_energy,E,Esquared,E4,c1,c1squa,c2,c2squa,c3,c3squa,BETA,sweep2,M);
            }
        }

        cout << "\nRunning Monte Carlo at T=" << T << " with BETA=" << BETA << endl;
        double steps2 = static_cast<double>(M);

        E = E/steps2;
        E_.push_back(E);

        Esquared = Esquared/steps2;
        Esquared_.push_back(Esquared);

        E4= E4/steps2;
        B= 1.0 - (E4/(3.0 * Esquared * Esquared));
        Binder.push_back(B);

        c1 = c1/steps2;
        c_1.push_back(c1);

        c1squa = c1squa / steps2;
        c1fluc= (c1squa - (c1 * c1));
        c_1_fluct.push_back(c1fluc);

        c2 = c2 / steps2;
        c_2.push_back(c2);

        c2squa = c2squa / steps2;
        c2fluc= (c2squa - (c2 * c2));
        c_2_fluct.push_back(c2fluc);

        c3 = c3 / steps2;
        c_3.push_back(c3);

        c3squa = c3squa / steps2;
        c3fluc= (c3squa - (c3 * c3));
        c_3_fluct.push_back(c3fluc);

        Cv = (Esquared - (E * E)) / (k * T * T );
        Cv_.push_back(Cv);

        cout << "E: " << E << endl;
        cout << "Esquared: " << Esquared << endl;
        cout << "Contacts chain1-wall: " << c1 << endl;
        cout << "Contacts chain1-wall fluctiations: " << c1fluc << endl;
        cout << "Contacts chain2-wall: " << c2 << endl;
        cout << "Contacts chain2-wall fluctiations: " << c2fluc << endl;
        cout << "Contacts chain1-chain2: " << c3 << endl;
        cout << "Contacts chain1-chain2 fluctiations: " << c3fluc << endl;
        cout << "Cv : " << Cv << endl;
        
        E = 0.00, Esquared = 0.00,E4=0.00, c1 = 0.00,c1squa = 0.00,c2 = 0.00,c2squa = 0.00 ,c3 = 0.00,c3squa = 0.00,Cv=0.00,c1fluc=0.00,c2fluc=0.00,c3fluc=0.00;
    }

    save_results_to_csv(T_values, E_, Esquared_,Binder, c_1,c_1_fluct, c_2,c_2_fluct, c_3,c_3_fluct, Cv_,"results20quadraticpotentialcc1cw120_2D.csv");
    return 0;
}