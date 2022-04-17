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

/* Code for exercise 2 of tcsf, basic functions for Smoothed Particles Hydrodynamics (SPH) */

int n=100;  //number of particles in system
//float n2=n;
float mass=1.0; //mass of all the particles is the same and set to 1.
float eta=8;  //eta for smoothing length. Is of order unity and between 2 and 10.

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
void output(vector<float> vec, char name[]){
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
    //cout << "vecsize= " << vec.size() <<"\n";
    for ( int i = 0; i < vec.size(); i++ ) {
        float i2=i;
        vec.at(i) = ((i2 / (vec.size()))+(1.0/(vec.size()+1))); // set element at location i to i /n
    }
   // vec.pop_back();  //removes back value
    //vec.pop_back();  //removes back value
    //printVec(vec);
    return vec;
}
//Function to calculate the smoothing length
float SmoothLength(float eta,vector<float> vec){
    float sum=0.0;  //dummy value to sum all the distances to nearest neighbour
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

//Function to find the positions of neighbours of particle at i. Return positions as a vector.
vector<int> neighbours(int pos, vector<float> vec, float h){
    vector<int> newvec;
    for (int i=0; i < vec.size(); i++){
        if ((abs(vec.at(pos)-vec.at(i)) < 2.0*h )){ //& (i != pos)
            newvec.push_back(i);
        }
    }
    //printVec(newvec);
    return newvec;
}
//Function to evaluate the M4 spline kernel for a given smoothing length and neigbour distance
float M4kernel(float pos1,float pos2,float h){
    float norm=(2.0/(3.0*h));
    float s=(abs(pos1-pos2) / h);
    float w=0.0;
    //cout << "pos1,pos2= " << pos1 << " " << pos2 << "\n";
    //cout << "s= " << s << "\n";
    //cout << "s^2= " << pow(s,2.0) << "\n";
    //cout << "(2-s)^3= " << pow((2-s),3.0) << "\n";
    if (s>=0.0 && s<=1.0){
        w=(norm*(1.0-((3.0/2.0)*(pow(s,2.0)))+((3.0/4.0)*(pow(s,3.0)))));
        //cout << "s= " << s << "\n";
        //cout << "w= " << w << "\n";
    }
    else if (s>1.0 && s<=2.0){
        w=(norm*((1.0/4.0)*(pow((2-s),3.0))));
        //cout << "s= " << s << "\n";
        //cout << "w= " << w << "\n";
    }
    else if (s>2.0){
        w=0.0;
        //cout << "s= " << s << "\n";
        //cout << "w= " << w << "\n";
    }
    else {
        cout << "error in kernel evaluation bozo";
    }
    //cout << "w= " << w << "\n";
    //cout << "norm= " << norm << "\n";
    return w;
}

//Now a function to find the neighbours of each particle using a brute force approach.
//Particle j is defined as a neighbour of i if |r_ij|<=2h
//then we calculate the density at each position and return it as a vector.

vector<float> density(float mass,float h, vector<float> vec){
    vector<float> density;
    for ( int i = 0; i < (vec.size()); i++ ){
        vector<int> nn=neighbours(i,vec,h);     //creating a vector of the nearest neighbours of particle i.
        //printVecInt(nn);
        //cout << i;
        float temp=0.0;
        //cout << temp << "=temp\n";
        for (int j =0; j < nn.size(); j++){
                //cout << "i= " << i <<"j= " << j << "\n";
                //cout << "i -> " << i << "j -> " << (nn.at(j)) << "\n";
                float w = (mass*(M4kernel(vec.at(i),vec.at(nn.at(j)),h)));
                //cout << "M4kernel= " << (M4kernel(vec.at(i),vec.at(nn.at(j)),h)) << "\n";
                //cout << "vec.at(i)= " << vec.at(i) << " vec.at(nn.at(j))= " << vec.at(nn.at(j)) << "\n";
                //cout << "w= " << w << "\n";
                temp =temp + w;
                //cout << temp << w << "temp, w\n";
        }
        density.push_back(temp);

    }
    return density;
}


int main () {
    vector<float> positions(n,0);  //vector of size n, each element is 0.
    //printVec(positions);
    //cout << positions.size() << "positions size\n";
    vector<float> newpos=setPositions(positions);
    //cout << newpos.size() << "newpos size\n";
    printVec(newpos);
    float h=SmoothLength(eta,newpos);
    cout <<"smoothing length -> " << h << "\n";
    //vector<int> test=neighbours(2,newpos,h);
    //printVecInt(test);
    vector<float> densities=density(mass,h,newpos);
    printVec(densities);
    cout << densities.size() << " <- densities size\n";
    //vector<int> nn=neighbours(12,newpos,h);
    //printVecInt(nn);
    char filename[]="densities.txt";
    output(densities,filename);
    char file2[]="positions.txt";
    output(newpos,file2);

}
