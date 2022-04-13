#include<iostream>
/* Code for exercise 2 of tcsf, basic functions for Smoothed Particles Hydrodynamics (SPH) */
using namespace std;

#include <iomanip>
using std::setw;

int n=100;  //number of particles in system
float n2=100; //same as n but as float for calulations
int main () {

   float sys[ n ]; // sys is an array of n floats

   // initialize elements of array n to 0
   for ( int i = 0; i < n; i++ ) {
      sys[ i ] = i / n2; // set element at location i to i /n
   }
   cout << "Element" << setw( 13 ) << "Value" << endl;

   // output each array element's value
   for ( int j = 0; j < n; j++ ) {
      cout << setw( 7 )<< j << setw( 13 ) << sys[ j ] << endl;
   }

   return 0;
}
