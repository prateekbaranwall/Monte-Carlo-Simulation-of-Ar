#include <bits/stdc++.h>
#define ed "\n"
#define sp " "
using namespace std;


void LJ_E_Change(vector<vector<double>>& coords, vector<double>& rTrial, double part, double L, double& dE, bool& cutoff)
{
  dE=0;

  // number of particle
  double num_particle = coords[2].size();

  for(double otherPart=0;otherPart<num_particle;otherPart++)
  {
    // to ignore the self interaction
    if(otherPart==part)
      continue;

    double hL = L/2.0;
    // getting the new distance b/w the particle
    vector<double> drNew(3);
    for(double i=0;i<3;i++)
    {
      drNew[i] = coords[i][otherPart] - rTrial[i];
      // Boundary condition
      if(drNew[i]>hL)
        drNew[i] = drNew[i]-L;
      else if(drNew[i]<-hL)
        drNew[i] = drNew[i]+L;
    }

    // getting the old distance b/w the particle
    vector<double> drOld(3);
    for(double i=0;i<3;i++)
    {
      drOld[i] = coords[i][otherPart] - coords[i][part];
      // Boundary condition
      if(drOld[i]>hL)
        drOld[i] = drOld[i]-L;
      else if(drOld[i]<-hL)
        drOld[i] = drOld[i]+L;
    }

    // getting the distance squared
    double dr2_New=0, dr2_Old=0;
      for(double i=0;i<3;i++)
      {
        dr2_New+=drNew[i]*drNew[i];
        dr2_Old+=drOld[i]*drOld[i];
      }
      if(cutoff==false || (cutoff==true && 2.5*2.5<=dr2_New&&dr2_New<=3.0*3.0 && 2.5*2.5<=dr2_Old&&dr2_Old<=3.0*3.0))
      {
          // LJ potential formula
          // U(r) = 4 * epsilon * [(sigma/r)^12 - (sigma/r)^6]
          // setting epislon = sigme = 1
          // U(r) = 4 * [(1/r)^12 - (1/r)^6]
          // dr2 is the square of the distance b/w the particle
          // invDR6 is the inverse to the power of 6 of the distance b/w the particle
          double invDR6_New = 1.0/(pow(dr2_New,3));
          double invDR6_Old = 1.0/(pow(dr2_Old,3));

          // calculating potential energy
          double eNew = (invDR6_New*(invDR6_New-1));
          double eOld = (invDR6_Old*(invDR6_Old-1));

          dE += eNew - eOld; 
      }
  }
  dE*=4;
}

void LJ_E(vector<vector<double>>& coords, double& L, double& energy, bool& cutoff)
{
  // number of particle
  double num_particle = coords[2].size();

  // taking the pair of distinct particle
  for(double partA=0; partA<num_particle-1; partA++)
  {
    for(double partB=partA+1; partB<num_particle; partB++)
    {
      double hL = L/2.0;
      vector<double> dr(3);
      // calculating the distance between two particle
      for(double i=0;i<3;i++) 
      {
          dr[i] = coords[i][partA]-coords[i][partB];
          // Boundary condition
          if(dr[i]>hL)
            dr[i] = dr[i]-L;
          else if(dr[i]<-hL)
            dr[i] = dr[i]+L;
      }

      // get the distance squared
      double dr2=0;
      for(double i=0;i<3;i++)
        dr2+=dr[i]*dr[i];
     if(cutoff == false || (cutoff==true && 2.5*2.5<=dr2&&dr2<=3.0*3.0)) // cut-off condition
     {
          // LJ potential formula
          // U(r) = 4 * epsilon * [(sigma/r)^12 - (sigma/r)^6]
          // taking epislon = sigma = 1
          // U(r) = 4 * [(1/r)^12 - (1/r)^6]
          // dr2 is the square of the distance b/w the particle
          // invDR6 is the inverse to the power of 6 of the distance b/w the particle
          double invDR6 = 1.0/(pow(dr2,3));

          // calculating potential energy
          energy = energy + (invDR6 * (invDR6-1));
     }
    }
  }

  energy*=4;
}
 
void initialize(double num_particle, double density, vector<vector<double>>& coords, double& L)
{
  // length of the cube
  L = pow((num_particle/density),(1.0)/3);

  // Finding the lowest perfect cube greater than or equal to the number of particle so that total particle fit in the cube
  double nCube = 2;
  while(nCube*nCube*nCube<num_particle)
    nCube++;

  // for counting the vacant spot
  vector<double> index {0,0,0};

  for(double part=0;part<num_particle;part++)
  {
    for(double i=0;i<3;i++)
      // set coordinates
    coords[i][part] = (int32_t)((index[i]+0.5)*(L/nCube));

      // Advancing the index
    index[0] = index[0]+1;
    if(index[0] == nCube)
    {
      index[0] = 0;
      index[1] = index[1]+1;
      if(index[1]==nCube)
      {
        index[1]=0;
        index[2] = index[2]+1;
      }
    }
  }
}


void Monte_carlo()
{
  double num_particle = 700; // total number of particle
  double density = 0.7; // Density of the system
  bool cutoff = false; // setting the cut-off

  double Temp=1.0, beta, maxDr=0.1; // max displacement
  beta = 1.0/Temp; // inverse of temperature

  double nSteps = 10000; // number of steps
  double L=0;
  vector<vector<double>> coords(3,vector<double>(num_particle,0));

  // inittializing the coordinates of the each particle
  initialize(num_particle, density, coords,L);


  // Calculate the initial potential energy
  double energy=0;
  LJ_E(coords,L,energy,cutoff);


  srand(time(0));
  for(double step=1;step<=nSteps;step++)
  {
    
    if(((float)rand()/RAND_MAX)*(num_particle+1) + 1 < num_particle) // partial displacement
    {
      for(double i=0;i<num_particle;i++)
      {
        // trail move for a partical i
        vector<double> rTrial(3);
        for(double j=0;j<3;j++)
        {
          rTrial[j] = coords[j][i] + maxDr*(((float)rand()/RAND_MAX)-0.5);

          // boundary condition
          if(rTrial[j]>L) 
            rTrial[j] = rTrial[j]-L;
          else if(rTrial[j]<0)
            rTrial[j] = rTrial[j]+L;
        }

         // Calculate the change in energy due to this trial move
        double dE=0;
         LJ_E_Change(coords,rTrial,i,L,dE,cutoff);

         if(((float)rand()/RAND_MAX)< exp(-beta*dE)) // condition to accept displacement move
         {
           // updating the position and energy
           for(double j=0;j<3;j++)
           {
             coords[j][i] = rTrial[j];
             energy = energy + dE;
           }
         }

      }
    }
    if((int32_t)step%10==0)
    cout<<energy<<ed;
  }
}

 
 signed main()
{
  ios_base::sync_with_stdio(false);cin.tie(NULL);cout.tie(NULL);
  Monte_carlo();
return 0;
}