#include<iostream>
#include <vector>
using namespace std;

#include <iomanip>
#include <cmath>
using std::abs;
using std::pow;
using std::setw;
using std::vector;

/* Code for exercise 2 of tcsf, basic functions for Smoothed Particles Hydrodynamics (SPH) */

int n=10;  //number of particles in system
float n2=n;
float mass=1; //mass of all the particles is the same and set to 1.
float eta=2;  //eta for smoothing length. Is of order unity and between 2 and 10.

//Function to print the values of a vector
void printVec(vector<float> vec){
    for (int i = 0; i < vec.size(); i++){
        cout << vec.at(i) << "\n";
    }
}
//Now write a function to edit the values of the vector to be evenly spaced between 0 and n but not exactly at n.
vector<float> setPositions(vector<float> vec){
    for ( int i = 0; i < vec.size(); i++ ) {
        float i2=i;
        vec.at(i) = (i2 / vec.size())+0.1; // set element at location i to i /n
    }
    vec.pop_back();  //removes back value
    printVec(vec);
    return vec;
}
//Function to calculate the smoothing length
float SmoothLength(float eta,vector<float> vec){
    float sum = 0;  //dummy value to sum all the distances to nearest neighbour
    for (int i = 0; i < (vec.size()-1); i++){
        float dist=vec.at(i+1)-vec.at(i);   //distance between i and i+1 positions
        sum=sum+dist;                          //summing distances
    }
    //cout << "sum -> " << sum <<"\n";
    float mean=sum/(vec.size()-1);            //mean of the distances, dividing sum of distances by number of distances
    //cout << "mean -> " << mean <<"\n";
    //cout << "n2 -> " << n2 <<"\n";
    float h=eta*(mean);
    //cout << "h -> " << h <<"\n";
    return h;
}

//Function to find the positions of neighbours of particel at i. Return positions as a vector.
vector<float> neighbours(int pos, float h){
    vector<float> neighbours(0,0);
    return neighbours;
}
//Function to evaluate the M4 spline kernel for a given smoothing length and neigbour distance
float M4kernel(float pos1,float pos2,float h){
    float s=abs(pos1-pos2) / h;
    float norm=2/(3*h);
    float w=0;
    if (s>=0 && s<1){
        float w=norm*(1-((3/2)*(pow(s,2.0)))+((3/4)*(pow(s,3.0))));
    } else if (s>=1 && s<=2){
        float w=norm*(1/4)*(pow((2-s),3.0));
    } else if (s>2){
        float w=0.0;
    } else {
        cout << "error in kernel evaluation bozo";
    }
    return w;


}

//Now a function to find the neighbours of each particle using a brute force approach.
//Particle j is defined as a neighbour of i if |r_ij|<=2h
//then we calculate the density at each position and return it as a vector.
/*
vector<float> density(float mass,float h, vector<float> vec){
    for ( int i = 0; i < (vec.size()); i++ ){
        rho=mass*M4kernel(vec.at(i+1),vec.at(i),h);
    }
    vector<float>density(n,0);
    return density;
}
*/

int main () {
    vector<float> positions(n+1,0);  //vector of size n+1, as last element will need to be removed. each element is 0.
    //printVec(positions);
    vector<float> newpos=setPositions(positions);
    cout <<"smoothing length -> " << SmoothLength(eta,newpos) << "\n";

}
