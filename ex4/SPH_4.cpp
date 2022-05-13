/* ============================================================
 * Theoretical and Computational Star Formation
 * Assignment sheet 3
 * Clara Kretzschmar
 * Tom Carron
 * Timon Danowski
 * ============================================================ */
 
#include <iostream>
#include <cmath>
#include <fstream> 
#include <vector>

// Create a data structure to store the important quantities
// mass, positions, velocities and accelerations are member variables of the class Particle
class Particle
{
  public:
    double position;
    double mass;
    double internalEnergy;
    double density;
    double pressure;
    double velocity;
    double acceleration;
    double changeInInternalEnergy;
    double timestep;
    std::vector<int> neighbours; 
    
    Particle()
    {
        // member variables are set to 0
        position = 0.0;
        mass = 0.0;
        internalEnergy = 0.0;
        density = 0.0;
        pressure = 0.0;
        acceleration = 0.0;
        velocity = 0.0;
        changeInInternalEnergy = 0.0;
        timestep = 0.0;
    }
};

// function declaration
void initParticles(std::vector<Particle>& particles, uint N_sph);
double getSmoothingLength(std::vector<Particle>& particles, double eta, uint N_sph);
void getNeighbours(std::vector<Particle>& particles, double h, uint N_sph);
void getDensity(std::vector<Particle>& particles, double h, uint N_sph);
void getPressure(std::vector<Particle>& particles, double gamma, uint N_sph);
void getAcceleration(std::vector<Particle>& particles, double h, uint N_sph);
void getChangeInInternalEnergy(std::vector<Particle>& particles, double h, uint N_sph);
double getTimestep(std::vector<Particle>& particles, double h, double gamma, double gamma_CFL, uint N_sph);
void timeIntegration(std::vector<Particle>& particles, double timestep, uint N_sph);

void writeToFile(std::vector<Particle>& particles, uint N_sph, double time, uint * createFile);

double kernel_function(double r_ij, double h);
double deriv_kernel_function(double r_ij, double h);

int main(){
    // set initial variables
    // number of particles
    int N_sph = 100;
    // parameter eta for neighbours, between 2 and 10
    double eta = 3.0;
    // Courant-Friedrichs-Lew constant between 0 and 1
    double gamma_CFL = 0.5;
    // ratio of specific heats
    double gamma = 5.0/3.0;
    
    // Stopping time
    double stop_t=0.309;
    // set initial time
    double init_t=0.0;
    
    // Initialise system
    // create vector of N particles
    std::vector<Particle> particles; 
    particles.resize(N_sph);
    initParticles(particles, N_sph);
    
    // define variables with arbitrary values
    double SmoothingLength = 1.0;
    double globalTimestep = 1.0;
  
    // file does not exist yet
    uint createFile = 0;
    // first time to write to the file: t=0.0
    double writeToFile_time = 0.0;
    
    // current time
    double t=init_t;
    // time loop
    while (t<=stop_t)   
    {
        // print current time
        std::cout << "________________________________________________" << std::endl;
        std::cout << "t = " << t << std::endl;
        
        // call the functions
        SmoothingLength = getSmoothingLength(particles, eta, N_sph);
        std::cout << "SmoothingLength = " << SmoothingLength << std::endl;
    
        getNeighbours(particles, SmoothingLength, N_sph);
        getDensity(particles, SmoothingLength, N_sph);
    
        getPressure(particles, gamma, N_sph);
        getAcceleration(particles, SmoothingLength, N_sph);
        getChangeInInternalEnergy(particles, SmoothingLength, N_sph);
       
        // calculate timestep
        globalTimestep = getTimestep(particles, SmoothingLength, gamma, gamma_CFL, N_sph);
        std::cout << "globalTimestep = " << globalTimestep << std::endl;
        
        // output
        std::cout << "particles[0].position = " << particles[0].position << std::endl;
        std::cout << "particles[0].velocity = " << particles[0].velocity << std::endl;
        std::cout << "particles[0].density = " << particles[0].density << std::endl;
        std::cout << "particles[0].internalEnergy = " << particles[0].internalEnergy << std::endl;
        
        // write to file every 0.01 time
        if (t >= writeToFile_time)
        {
            // call the function to write postion, density, .. to a file
            writeToFile(particles, N_sph, t, &createFile);
            // write to file every 0.01 time units
            writeToFile_time += 0.01;
        }
        
        // time integration
        timeIntegration(particles, globalTimestep, N_sph);
    
        t += globalTimestep;
        
    }
    
    return 0;
}


/* ============================================================
 * functions
 * ============================================================ */

// A method to initialise the problem
void initParticles(std::vector<Particle>& particles, uint N_sph)
{
    // compute the distance dx between the evenly spaced particles
    double dx = 1.0 / (N_sph);
    // loop over all particles
    for(uint i=0; i<N_sph; i++)
    {
        // particles evenly spaced between 0 and 1
        particles[i].position = 0.005 + i * dx; //(i + 1.0) * dx;
        // mass = 1 / number of particles
        particles[i].mass = 1.0 / N_sph;
        // all particles have the same internal energy 1
        particles[i].internalEnergy = 1.0;
    }
}

// Function to compute the smoothing length
double getSmoothingLength(std::vector<Particle>& particles, double eta, uint N_sph)
{
    // calculate the mean distance to the closest neighbour, first it is set to 0
    double mean_dis = 0.0;
    // create array that stores the distance to the closest neighbour for each particle
    double min_r[N_sph];
    // loop over all particles
    for(uint i=0; i<N_sph; i++)
    {
        // distance to the closest neighbour is set to the maximum distance 1
        min_r[i] = 100.0;
        for(uint j=0; j<N_sph; j++)
        {
            if (i!=j)
            {
                // compute the distance between particle i and j
                double r_ij_abs = fabs(particles[i].position - particles[j].position);
                // if distance to particle j is smaller than the former minimum distance: replace it
                min_r[i] = fmin(min_r[i], r_ij_abs);
            }
        }
        // Add minimum distance of particle i 
        mean_dis += min_r[i];
            
    }
    // calculate the mean distance by dividing by the number of particles
    mean_dis /= N_sph;
    
    // smoothing length = eta * mean distance to the closest neighbour
    // return the smoothing length
    return eta * mean_dis;
}

// Function to determine the neighbours of each particle
void getNeighbours(std::vector<Particle>& particles, double h, uint N_sph)
{
    // loop over all particles
    for(uint i=0; i<N_sph; i++)
    {
        // remove all elements of particles[i].neighbours
        particles[i].neighbours.clear();
        
        for(uint j=0; j<N_sph; j++)
        {
            // if a particle is closer than 2*h, j is a neighbour of i: 
            // (i is also a neighbour of i)
            // append index of j to neighbours array of particle i
            if (fabs(particles[i].position - particles[j].position) < 2.0 * h)
            {
                particles[i].neighbours.push_back(j);
            }
        }
    }
}

// Function to calculate the density for each particle
void getDensity(std::vector<Particle>& particles, double h, uint N_sph)
{
    // loop over all particles
    for(uint i=0; i<N_sph; i++)
    {
        // density is set to zero
        particles[i].density = 0.0;
        // Calculate the number of neighbours of particle i, i.e. the length of the neighbours array
        uint N_neighbours = particles[i].neighbours.size();
        // summation over all neighbours of particle i
        for(uint j=0; j<N_neighbours; j++)
        {
            // index of particle j in the neighbours array of particle i
            uint index = particles[i].neighbours[j];
            // r_ij_abs is the (absolute value of the) distance between i and j
            double r_ij_abs = fabs(particles[i].position - particles[index].position);
            // add the contribution of neighbours j to the density. In the beginning the density is set to zero.
            particles[i].density += particles[index].mass * kernel_function(r_ij_abs, h);
        }
    }
}

// Function to calculate the pressure for each particle
void getPressure(std::vector<Particle>& particles, double gamma, uint N_sph)
{
    // loop over all particles
    for(uint i=0; i<N_sph; i++)
    {
        // Calculate the pressure with equation (3)
        particles[i].pressure = (gamma - 1.0) * particles[i].density * particles[i].internalEnergy;
    }
}

// Function to calculate the acceleration for each particle
void getAcceleration(std::vector<Particle>& particles, double h, uint N_sph)
{
    // loop over all particles
    for(uint i=0; i<N_sph; i++)
    {
        //acceleration is set to zero
        particles[i].acceleration = 0.0;
        // Calculate the number of neighbours of particle i, i.e. the length of the neighbours array
        uint N_neighbours = particles[i].neighbours.size();
        // summation over all neighbours of particle i
        for(uint j=0; j<N_neighbours; j++)
        {
            // index of particle j in the neighbours array of particle i
            uint index = particles[i].neighbours[j];
            // r_ij is the distance between i and j
            double r_ij = particles[i].position - particles[index].position;
            // r_ij_abs is the absolute value of r_ij
            double r_ij_abs = fabs(r_ij);
            // sign of r_ij
            double sign_r_ij = 0.0;
            if(r_ij >= 0.0)
            { sign_r_ij = 1.0; }
            else if (r_ij < 0.0)
            { sign_r_ij = -1.0; }
            else
            { std::cout << "ERROR: r_ij is not real (getAcceleration)" << r_ij << std::endl; }
            
            // add the contribution of neighbours j to the acceleraton. In the beginning the acceleration is set to zero.
            // Calculate the acceleration with equation (1)
            
            double temp = 1.0;
            if ((particles[i].density != 0.0) && particles[index].density != 0.0)
            { double temp1 = (particles[i].pressure / particles[i].density / particles[i].density);
              double temp2 = (particles[index].pressure / particles[index].density / particles[index].density);
              temp = temp1 + temp2; }
            else 
            { std::cout << "ERROR: density is zero (getAcceleration)" << std::endl; }
            particles[i].acceleration += particles[index].mass * temp * deriv_kernel_function(r_ij_abs, h) * sign_r_ij;
        }
        // change the sign of the acceleration
        particles[i].acceleration *= -1.0;
    }
}

// Function to calculate the change in internal energy for each particle
void getChangeInInternalEnergy(std::vector<Particle>& particles, double h, uint N_sph)
{
    // loop over all particles
    for(uint i=0; i<N_sph; i++)
    {
        // Change in internal energy is set to zero
        particles[i].changeInInternalEnergy = 0.0;
        // Calculate the number of neighbours of particle i, i.e. the length of the neighbours array
        uint N_neighbours = particles[i].neighbours.size();
        // summation over all neighbours of particle i
        for(uint j=0; j<N_neighbours; j++)
        {
            // index of particle j in the neighbours array of particle i
            uint index = particles[i].neighbours[j];
            // r_ij is the distance between i and j
            double r_ij = particles[i].position - particles[index].position;
            // r_ij_abs is the absolute value of r_ij
            double r_ij_abs = fabs(r_ij);
            // sign of r_ij
            double sign_r_ij = 0.0;
            if(r_ij >= 0.0)
            { sign_r_ij = 1.0; }
            else if (r_ij < 0.0)
            { sign_r_ij = -1.0; }
            else
            { std::cout << "ERROR: r_ij is not real (getChangeInInternalEnergy)" << std::endl; }
            
            // add the contribution of neighbours j to the acceleraton. In the beginning the density is set to zero.
            // Calculate the change in internal energy with equation (2)
            particles[i].changeInInternalEnergy += particles[index].mass * (particles[i].velocity - particles[index].velocity) * deriv_kernel_function(r_ij_abs, h) * sign_r_ij;
            
        // Mulitply the change with the prefactor
        if(particles[index].density!=0.0)
        { particles[i].changeInInternalEnergy *= particles[i].pressure / particles[i].density / particles[i].density; }
        else 
        { std::cout << "ERROR: density is zero" << std::endl; }
        
        }
    }
}

// Function to calculate the global timestep 
double getTimestep(std::vector<Particle>& particles, double h, double gamma, double gamma_CFL, uint N_sph)
{
    double epsilon = 1e-50;
    // loop over all particles
    for(uint i=0; i<N_sph; i++)
    {
        // Calculate the timestep for each particle with equation (8)
        
        // Calculate the sound speed with equation (4)
        double sound_speed = sqrt(gamma * (gamma - 1.0) * particles[i].internalEnergy);
        
        // calculate the divergence of the velocity (Exercise sheet 2 equation (6))
        double grad_v = 0.0;
        // Calculate the number of neighbours of particle i, i.e. the length of the neighbours array
        uint N_neighbours = particles[i].neighbours.size();
        // summation over all neighbours of particle i
        for(uint j=0; j<N_neighbours; j++)
        {
            // index of particle j in the neighbours array of particle i
            uint index = particles[i].neighbours[j];
            // r_ij is the distance between i and j
            double r_ij = particles[i].position - particles[index].position;
            // r_ij_abs is the absolute value of r_ij
            double r_ij_abs = fabs(r_ij);
            // sign of r_ij
            double sign_r_ij = 0.0;
            if(r_ij >= 0.0)
            { sign_r_ij = 1.0; }
            else if (r_ij < 0.0)
            { sign_r_ij = -1.0; }
            else
            { std::cout << "ERROR: r_ij is not real (getTimestep)" << std::endl; }
            
            // add the contribution of neighbours j to the divergence. In the beginning the density is set to zero.
            if(particles[index].density!=0.0)
            { grad_v += particles[index].mass / particles[index].density * particles[index].velocity * deriv_kernel_function(r_ij_abs, h) * sign_r_ij; }
            else
            { std::cout << "ERROR: density is zero (getTimestep)" << std::endl; }
        }
        
        // set temp1 to an arbitrary high number
        double temp1 = 100.0;
        // if the denominator is non-zero: calculate temp1, otherwise print warning
        if((h * fabs(grad_v) + sound_speed) != 0.0)
        { temp1 = h / (h * fabs(grad_v) + sound_speed); }
        else 
        { std::cout << "ERROR: divided by zero (getTimestep)" << std::endl; }
       
        double temp2 = sqrt(h / (fabs(particles[i].acceleration) + epsilon));
    
        particles[i].timestep = gamma_CFL * fmin(temp1, temp2);
    }
    // determine the minimum of the timestep of all particles (equation (9))
    // set min_timestep to 0.001
    double min_timestep=0.01;
    // loop over all particles
    for(uint i=0; i<N_sph; i++)
    {
        min_timestep = fmin(min_timestep, particles[i].timestep);
    }
    // return the global timestep
    return min_timestep;
}

void timeIntegration(std::vector<Particle>& particles, double timestep, uint N_sph)
{
    // loop over all particles
    for(uint i=0; i<N_sph; i++)
    {
        // update velocity, position and internal energy (overwrite the member variables)
        // with equation (6) and (7) (first-order Euler time integration)
        particles[i].internalEnergy += particles[i].changeInInternalEnergy * timestep;
        particles[i].position += particles[i].velocity * timestep;
        particles[i].velocity += particles[i].acceleration * timestep;
    }
}


// Function to write each particles positions and velocities to a file
void writeToFile(std::vector<Particle>& particles, uint N_sph, double time, uint * createFile)
{
    // open file
    std::ofstream file01("tcsf_SPH_sheet3.txt", std::ofstream::app);
        
    if (*createFile == 0)
    {
        // create head line
        file01 << "time" << "\t" << "position" << "\t" << "velocity" << "\t" << "internal energy" << std::endl;
        *createFile = 1;
    }
    
    // first column: current time
    file01 << time << "\t";
    
    // loop over all particles
    for(uint i=0; i<N_sph; i++)
    {
        // write particle positions to file
        file01 << particles[i].position << "\t";
    }
    // seperator 1 between position and density
    //file01 << "\t" << "\t";
    
    // loop over all particles again
    for(uint i=0; i<N_sph; i++)
    {
        // write particle densities to file
        file01 << particles[i].density << "\t";
    }
    // seperator 2 between density and velocity
    //file01 << "\t" << "\t";
     // loop over all particles again
    for(uint i=0; i<N_sph; i++)
    {
        // write particle velocities to file
        file01 << particles[i].velocity << "\t";
    }
    
    // seperator 2 between velocity and internal energy
    //file01 << "\t" << "\t";
     // loop over all particles again
    for(uint i=0; i<N_sph; i++)
    {
        // write particle internal energy to file
        file01 << particles[i].internalEnergy << "\t";
    }

    // end line
    file01 << std::endl;
     
    // close file
    file01.close();
}

/* ============================================================
 * kernel function and its derivative
 * ============================================================ */

// Compute kernel function for given distance r_ij and smoothing length h
double kernel_function(double r_ij, double h)
{
    // print a warning if the smoothing length is zero
    if (h == 0)
    {
        std::cout << "WARNING: smoothing length is zero" << std::endl;
    }
    // Define the constant s
    double s = r_ij / h;
    // normalisation constant in 1D
    double W = 2.0 / 3.0 / h;
    // consider cases for s, then calculate W
    if (s >= 0.0 && s <= 1.0)
    {
        W *= (1.0 - 1.5 * s * s + 0.75 * pow(s,3));
    }
    else if (s >= 1.0 && s <= 2.0)
    {
        W *= 0.25 * pow((2.0 - s),3.0);
    }
    else if (s > 2.0)
    {
        W = 0.0;
    }
    else 
    {
        std::cout << "WARNING: kernel function cannot be computed" << std::endl;
    }
    // return the value of the kernel function
    return W;
}

// Derivative of the kernel function for given distance r_ij and smoothing length h
double deriv_kernel_function(double r_ij, double h)
{
    // print a warning if the smoothing length is zero
    if (h == 0.0)
    {
        std::cout << "WARNING: smoothing length is zero" << std::endl;
    }
    // Define the constant s
    double s = r_ij / h;
    // normalisation constant in 1D
    double dW = 2.0 / 3.0 / h / h;
    // consider cases for s, then calculate dW
    if (s >= 0.0 && s < 1.0)
    {
        dW *= -3.0 * s + 2.25  * s * s;
    }
    else if (s >= 1.0 && s < 2.0)
    {
        dW *= -0.75 * pow((2.0 - s),2.0);
    }
    else if (s >= 2.0)
    {
        dW = 0.0;
    }
    else 
    {
        std::cout << "WARNING: Derivative of the kernel function cannot be computed" << std::endl;
    }
    // return the value of the derivative of the kernel function
    return dW;
}
