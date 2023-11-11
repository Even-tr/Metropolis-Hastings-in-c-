/*
A program which runs a symemtric Metropolis Hastings algorithm on the 
on the Correlates of War data set. It assumes a exponential model
*/

#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <stdexcept>


// Defining generators and distributions.
const double prop_stdev_theta = 0.3; // standard deviation for the proposition distribution of theta
const double prop_stdev_year = 2; // standard deviation for the proposition distribution of a and b
std::default_random_engine generator;
std::normal_distribution<double> normal_theta(0, prop_stdev_theta);
std::normal_distribution<double> normal_year(0, prop_stdev_year);
std::uniform_real_distribution<double> uniform_01(0, 1);
std::uniform_int_distribution<int> RandomShift(0, 8);


double prop_year(){
    double sample = normal_year(generator);
    return (double) std::round(sample);

}

double gamma_pdf(double x, double alpha, double beta){
    return pow(beta, alpha)*pow(x, alpha -1) * exp(-beta*x)/std::tgamma(alpha);
}

int min_a = 10;
int max_b = 40;

// Struct which store the parameters.
// Implements som quality of life methods like addition and 
// string representation.
struct Parameters
{
    double thetaL;
    double thetaR;
    double a;
    double b;

    // Methods
    Parameters operator+(const Parameters& other) const {
        Parameters result;
        result.thetaL = this->thetaL + other.thetaL;
        result.thetaR = this->thetaR + other.thetaR;
        result.a = this->a + other.a;
        result.b = this->b + other.b;
        return result;
    };

    std::string toStringBracket() const {
        return "{" + toString() + "}"; 
    };
    std::string toString() const {
        return std::to_string(thetaL) + ", " + std::to_string(thetaR) +", "+ std::to_string(a)+ ", "+ std::to_string(b); 
    };
    std::string parameterNames()const{
        return "thetaL, thetaR, a, b";
    }
    
    bool outsideBoundaries(){
        return (this->thetaL < 0) or 
        (this->thetaR   < 0) or 
        (this->a        < min_a) or 
        (this->b        > max_b) or 
        (this->a        >= this-> b);
        // (this-> thetaL  < this->thetaR);
    };

    std::string verboseOutsideBoundaries(){
        if (this->thetaL < 0){
            return "thetaL too low";
        };
        if (this->thetaR < 0) {
            return "thetaR too low";
        };
        if (this->a < min_a) {
            return "a too low";
        };
        if (this->b > max_b) {
            return "b too high";
        };
        if (this->a >= this-> b){
            return "a greather than or equal to b";
        };
        if (this-> thetaL < this->thetaR){
            return "thetaR greather thetaL";
        }
        
    return "";
    };
};

// A struct for a single observation. Edit to fit your needs
struct Observation
{
    double year;
    double y;

    std::string toString() const {
        return  std::to_string(year) + ", " + std::to_string(y); 
    };

    std::string colnames() const {
        return  "x \ty"; 
    };

};

// struct which stores the observation data s.t. it can be referenced with a pointer
struct Data
{
    int size;
    Observation* rows;

    // Displays the table of observations.
    void display(){
        std::cout << rows[0].colnames() << std::endl;
        for (size_t i = 0; i < size; i++)
        {
            std::cout << (rows[i]).toString() << std::endl;
        }
    }
};

double exponential(double x, double theta){
    // theta is the scale parameter. This means E[x] = theta
    return 1/theta*exp(-1/theta*x);
}

// Function which returns the prior for a set of parameters. 
// Currently constant, i.e. all parameters equally likely a prior.
double lprior(Parameters P){

    // return log((double) (P.b-P.a) == 1); // Changepoint
    // return 1;
    // Exponential decay prior for window length
    double dprior = exp(-std::abs(P.b - P.a));
    double thetaLprior = 1/P.thetaL;
    double thetaRprior = 1/P.thetaR;
    // double thetaLprior = gamma_pdf(P.thetaL, 1, 1);
    // double thetaRprior = gamma_pdf(P.thetaR, 1, 1);
    return log(dprior) + log(thetaLprior) + log(thetaRprior);
}

// Defines the model for parameter change, wich is linear from thetaR to thetaL
double theta_linear(Parameters P, double year){
    // return P.thetaL;
    if (year < P.a){
        return P.thetaL;
    } else if (year > P.b){
        return P.thetaR;
    } else{
        return P.thetaL + (year - P.a)/(P.b - P.a) * (P.thetaR - P.thetaL);
    }
} 

// Function which returns a proposition parameter
Parameters prop(Parameters P){
    Parameters shift;
    shift.thetaL = normal_theta(generator);
    shift.thetaR = normal_theta(generator);
    // Continous shift
    // shift.a = normal_year(generator);
    // shift.b = normal_year(generator);
    // int index = RandomShift(generator);
    // shift.a = shiftArray[index][0];
    // shift.b = shiftArray[index][1];
    shift.a = prop_year();
    shift.b = prop_year();
    return P + shift;
};

// Calculates likelihood of a single observation y given P parameters
double model(Parameters P,  Observation y){
    return exponential(y.y, theta_linear(P, y.year));
}

// calculates log likelihood of a given combination of parameters and observations
double llik(Parameters P, Data d){
    double l = 0;
    for (size_t i = 0; i < d.size; i++){
        l = l + log(model(P, d.rows[i]));
    }
    return l;
}

// Reads a file and fills a Data struct with its contents.
void readData(Data* d, std::string filename){
    std::string line;
    std::ifstream file(filename);
    std::vector<std::string> xs;
    std::vector<std::string> ys;
    
    int length = 0; // to be overwritten
    // Read file into column vectors and count length
    if (!file) {
        std::cerr << "Failed to open the file." << std::endl;
        return;
    }

    int col = 0; // keeps track of which columns to fill
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string word;
        while (iss >> word) {
            switch (col%2)
            {
            case 0:
                xs.push_back(word);
                break;
            case 1:
                ys.push_back(word);
                break;
            default:
            std::cout<< "ERROR in reading " + filename << std::endl;
                break;
            }
            col ++;            
        }
        length++;
    }
    file.close();

    // Initialize data array
    d->size = length;
    d->rows = new Observation[length];

    // fill data array with column vectors, converting to floats.
    for (size_t i = 0; i < length; i++){
        Observation O = {std::stof(xs[i]), std::stof(ys[i])};
        d->rows[i] = O;
    }
}


// Formats the chain to a nice string for printing or saving
std::string chainToString(Parameters * chain, int chain_length){
    std::string res = "";
    for (size_t i = 0; i < chain_length; i++){
        res = res + chain[i].toString()  + "\n";
    }
    return res;
    
}

int main(){
    // Read british mining data set
    Data data;
    std::string data_path = "data/COW2.txt";
    // std::string data_path = "synthetic.txt";
    std::string out_path = "results/COW2_chain.txt";
    std::cout << "Reading data from "<< data_path << std::endl;
    readData(&data, data_path);
    // data.display(); std::cout<< "\n\n";

    // Initialize the chain
    const int chain_length = 10000;
    int accepted = 0;
    Parameters* chain = new Parameters[chain_length];
    // Parameters chain [chain_length];
    Parameters thetaZero = {2, 1, 10, 15};
    chain[0] = thetaZero;

    Parameters theta;
    Parameters theta_prop;
    
    // The MCMC must start at a legal vaule, or else it might diverge
    if (thetaZero.outsideBoundaries()){
        std::cout << "Illegal start: " << thetaZero.toStringBracket() << std::endl;
        return 1;
    }

    // Run chain
    for (size_t i = 0; i < chain_length - 1; i++){
        theta = chain[i];
        theta_prop = prop(theta);

        // disregard impossible values
        if (theta_prop.outsideBoundaries()){
            chain[i+1] = theta;
        } 

        else{
            // Calculate transition probabilities min(1, p(x') / p(x)).
            // where p(x) = exp(log(prior(x)) + log(likelihood(y | x)))
            // with x being the parameters, while y is the observation
            double d = std::min(1.0, exp(llik(theta_prop, data) + lprior(theta_prop) - llik(theta, data) - lprior(theta)));
            if (uniform_01(generator) < d){
                accepted ++;
                chain[i+1] = theta_prop;
            }
            else{
                 chain[i+1] = theta;
            }
        }        
    }

    std::cout << "Acceptance rate " << (float) accepted / (float) chain_length << std::endl;



    std::cout << "Saving chain to " << out_path << std::endl;
    // Write the chain to file 
    std::ofstream myfile;
    myfile.open(out_path);
    myfile << theta_prop.parameterNames() << std::endl;
    myfile << chainToString(chain, chain_length);
    myfile.close();
    std::cout << "Done!" << std::endl;

    // Free memory
    delete[] chain;
    delete[] data.rows; 
}