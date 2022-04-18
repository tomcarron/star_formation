#include<iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
using std::abs;
using std::pow;
using std::setw;
using std::vector;
using namespace std;
using std::string;
using std::setprecision;

/* Code for exercise 2 of tcsf, basic functions for Smoothed Particles Hydrodynamics (SPH) */

int n=100;                          //number of particles in system
float mass=1.0;                     //mass of all the particles is the same and set to 1.
vector<float> etas={2,4,6,8,10};    //eta for smoothing length. Is of order unity and between 2 and 10.

//Function to print the values of a vector of floats
void printVec(vector<float> vec){
    for (int i = 0; i < vec.size(); i++){
        cout << vec.at(i) << "\n";
    }
}
//Function to print the values of a vector of integers
void printVecInt(vector<int> vec){
    for (int i = 0; i < vec.size(); i++){
        cout << vec.at(i) << "\n";
    }
}
// Function to take a vector of floats and save to a txt file named "name.txt" for plotting in python.
void output(vector<float> vec, string name){
    fstream file;
    file.open(name,ios_base::out);

    for(int i=0;i<vec.size();i++)
    {
        file<<vec.at(i)<<endl;
    }

    file.close();


}
//Now write a function to edit the values of the vector to be evenly spaced between 0 and n but not exactly at n.
vector<float> setPositions(vector<float> vec){
    for ( int i = 0; i < vec.size(); i++ ) {
        float i2=i;
        vec.at(i) = ((i2 / (vec.size()))+(1.0/(vec.size()+1))); // set element at location i to i /n
    }
    return vec;
}
//Function to calculate the smoothing length
float SmoothLength(float eta,vector<float> vec){
    float sum=0.0;  //dummy value to sum all the distances to nearest neighbour
    for (int i = 0; i < (vec.size()-1); i++){
        float dist=vec.at(i+1)-vec.at(i);   //distance between i and i+1 positions
        sum=sum+dist;                          //summing distances
    }
    float mean=sum/(vec.size()-1);            //mean of the distances, dividing sum of distances by number of distances
    float h=eta*(mean);
    return h;
}

//Function to find the positions of neighbours of particle at i. Return positions as a vector.
vector<int> neighbours(int pos, vector<float> vec, float h){
    vector<int> newvec;
    for (int i=0; i < vec.size(); i++){
        if ((abs(vec.at(pos)-vec.at(i)) < 2.0*h )){ //& (i != pos)
            newvec.push_back(i);
        }
    }
    return newvec;
}
//Function to evaluate the M4 spline kernel for a given smoothing length and neigbour distance
float M4kernel(float pos1,float pos2,float h){
    float norm=(2.0/(3.0*h));
    float s=(abs(pos1-pos2) / h);
    float w=0.0;
    if (s>=0.0 && s<=1.0){
        w=(norm*(1.0-((3.0/2.0)*(pow(s,2.0)))+((3.0/4.0)*(pow(s,3.0)))));
    }
    else if (s>1.0 && s<=2.0){
        w=(norm*((1.0/4.0)*(pow((2-s),3.0))));
    }
    else if (s>2.0){
        w=0.0;
    }
    else {
        cout << "error in kernel evaluation bozo";
    }
    return w;
}

//Now a function to find the neighbours of each particle using a brute force approach.
//Particle j is defined as a neighbour of i if |r_ij|<=2h
//then we calculate the density at each position and return it as a vector.

vector<float> density(float mass,float h, vector<float> vec){
    vector<float> density;                                                       //intialising vector
    for ( int i = 0; i < (vec.size()); i++ ){                                   //loop through each position
        vector<int> nn=neighbours(i,vec,h);                                     //creating a vector of the nearest neighbours of particle i.
        float temp=0.0;
        for (int j =0; j < nn.size(); j++){                                     //loop through each neighbour
                float w = (mass*(M4kernel(vec.at(i),vec.at(nn.at(j)),h)));      //contribution of neighbour a j to density at i.
                temp =temp + w;
        }
        density.push_back(temp);                                               //adding density to end of vector

    }
    return density;                                                             //returning vector of densities
}


int main () {
    vector<float> positions(n,0);                               //vector of size n, each element is 0.
    vector<float> newpos=setPositions(positions);               // set the positions to n evenly spaced values between 0 and 1, not including 0 and 1.
    for (int i=0; i < etas.size(); i++){                        // loop through values of vector of eta values.
        float eta=etas[i];                                      //current eta
        string eta_string=to_string(llround(eta));              //eta as string for naming
        float h=SmoothLength(eta,newpos);                       //calculate smoothing length
        vector<float> densities=density(mass,h,newpos);         // compute densities at each position
        string filename=eta_string+"_densities.txt";            //filename to save densities vector
        output(densities,filename);                             //output densities to .txt
        string file2=eta_string+"_positions.txt";               //filename to save positions vector
        output(newpos,file2);                                   //output densities to .txt

    }
}
