#include <iostream>
#include <math.h>
#include <vector>
#include <cassert>
#include <map>
#include <stdlib.h> 
using namespace std;


//////Define a function to make dot product
double dot(vector<double> a, vector<double> b){
double d=0;
assert(a.size()==b.size());
for(int i=0; i<a.size();i++){
d+=a[i]*b[i];
}
return d;
};////End of dot




int main(){

///define some variables	

////Here insert the minimum radius of all your dust particles
map<string, double> q_map; ///this is a vector containing charge on the particles. 
map<string, double> d_map; ///
vector<double> s,p,lambda,temp_vec;
bool bool_value;
double e=1.602e-19;///charge of ion and electron
double k_B=1.38044e-23;///Boltzman constant
double T_ion=298;///temperature of ion
double T_elec=26400;//temperature of electron
double n_ion=1.6e15;///number density of ions
double m_ion=40*1.67e-27;///mass of ion
double f_ion,f_elec,J_ion,I_ion,q_total_ion,q_patch_ion;
double q_total_AGG=0;
double pi=3.14;
double v0;
double R_min=4.5e-6;
double del_alpha = pi/8;
double del_beta = pi/8;
double del_v = 500;
double del_t=1e-8;///del_time
double del_theta, del_phi;
int N; ///nember of intervals for a particle of radius R.
map<string , vector<double> > AGG;///this is a map between the particles name in string and a vector with the following structure. The zeroth, first, second, and third components are x, y, z, and raduis of the dust particle. 


/////Define your dust particles here 
AGG["par0"].push_back(0);AGG["par0"].push_back(0);AGG["par0"].push_back(0);AGG["par0"].push_back(4.5e-6);
AGG["par1"].push_back(9e-6);AGG["par1"].push_back(0);AGG["par1"].push_back(0);AGG["par1"].push_back(4.5e-6);


q_map.clear();///to make sure the vector is empty.
for(map<string, vector<double> >::iterator it=AGG.begin(); it!=AGG.end(); it++){///it has two methods ->first and ->second
N=abs(4/R_min*it->second[3]);
del_theta=pi/N;
del_phi=pi/N;

q_total_ion=0;////for each patch make the total charge zero. This is to avoid adding diff. particles charges.
for(double theta=pi/(2*N);theta<pi;theta+=del_theta){///this is theta of each patch
for(double phi=pi/(2*N);phi<2*pi;phi+=del_phi){///this is phi of each patch on the sphere, " dust".
J_ion=0;
d_map.clear();

if(theta > pi/2)v0=0;
else v0=7000;

for(double v=0; v< 30000; v+=del_v){///this is loop over v of incoming particles              ????????????????????????????????????
for(double alpha=0; alpha <= pi/2;alpha+=del_alpha){///this is loop over alpha
for(double beta=0;beta<2*pi;beta+=del_beta){//this is loop over beta of each patch on the dust sphere

	
for(map<string, vector<double> >::iterator itt=AGG.begin(); itt!=AGG.end(); itt++){///itt is the second sum over particles and is used to calculate the "line of sight"
if(itt->first==it->first)continue;///line of sigh of particle is calculated with respect ot other particles not itself.
s.clear();///vector s is from the center of the patch to the center of other particles
p.clear();
lambda.clear();
temp_vec.clear();

////to understand these vectors please see the documentation.
p.push_back(it->second[3]*sin(theta)*cos(phi));
p.push_back(it->second[3]*sin(theta)*sin(phi));
p.push_back(it->second[3]*cos(theta));
s.push_back(itt->second[0]-it->second[0]-p[0]);
s.push_back(itt->second[1]-it->second[1]-p[1]);
s.push_back(itt->second[2]-it->second[2]-p[2]);
lambda.push_back(sin(alpha)*cos(beta));//this is the unit vector
lambda.push_back(sin(alpha)*sin(beta));
lambda.push_back(cos(alpha));
temp_vec.push_back(s[0]-dot(s,lambda)*lambda[0]);
temp_vec.push_back(s[1]-dot(s,lambda)*lambda[1]);
temp_vec.push_back(s[2]-dot(s,lambda)*lambda[2]);
d_map[itt->first]=sqrt(dot(temp_vec,temp_vec));

}////this is end of loop over itt. The second loop over particles.	


/////////////////determine bool_value
bool_value = true;
for(map<string, vector<double> >::iterator itt=AGG.begin(); itt!=AGG.end(); itt++){///itt is the second sum over particles and is used to calculate the "line of sight"
if(itt->first==it->first)continue;///line of sigh of particle is calculated with respect ot other particles not itself.
if(d_map[itt->first] < itt->second[3]) bool_value=false;
}///end of loop over itt the third loop over particles
//////end of bool_value determination


if(bool_value==true){///this is condition on the line of sight
f_ion=e*n_ion*pow((m_ion/(2*pi*k_B*T_ion)),3/2)*pow(v,3)*cos(alpha)*sin(alpha)*exp(m_ion*(pow(v,2)*pow(v0,2)-2*v*v0*(cos(alpha)*cos(theta)+sin(alpha)*cos(beta)*sin(theta)))/(2*k_B*T_ion));
J_ion=f_ion*del_v*del_alpha*del_beta+J_ion;
}///end of condition on the line of sight


}///end of loop over beta
}///end of loop over alpha	
}///end of loop over v of incoming particles, ions and electrons.


////here J_ion of a patch is calculated. Now we need to calculate I and q_patch and q_total
I_ion=J_ion*pow(it->second[3],2)*sin(theta)*del_theta*del_phi;///J*del_area is current
q_patch_ion=I_ion*del_t;
q_total_ion+=q_patch_ion;


/////save q_patch_ion and q_total_ion 
//??????????????????????????????????????????
//??????????????????????????????????
//???????????????????????????????????????????


}///en of loop over phi of patch
}///end of loop over theta of patch

cout << " par: " << it->first << "\n q_total_ion: " << q_total_ion << endl;

q_total_AGG+=q_total_ion;///this is the sum of charges over all the particles 
cout << "q_total_AGG: " << q_total_AGG << endl;

}////end of loop over par "dust particles"	







}
