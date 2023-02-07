#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cassert>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;
///----------------------------------------------------------------------------------------------------------------------------------
///----------------------------------------------------------------------------------------------------------------------------------
const bool plot_vtk = true;
const int nx = 150, ny = nx/2, np = 9;
vector<const int> cx = {0, 1, 0, -1,  0, 1, -1, -1,  1},
									cy = {0, 0, 1,  0, -1, 1,  1, -1, -1},
								 opp = {0, 3, 4, 1, 2, 7, 8, 5, 6};
const double cs2 = 1./3., rho0 = 1., Pr = 0.71, Ra = 2000;
const double T_h = 1.0, T_l = 0.0, T_0 = (T_h+T_l)*0.5, DT = T_h-T_l;
const double U_ref = 0.1, gbeta = U_ref*U_ref/(DT*(ny-1)), ni = U_ref*(ny-1)*sqrt(Pr/Ra), tau = ni/cs2+0.5, omega = 1./tau, omega1 = 1.-omega, T_ref = (ny-1)/U_ref;
const int nsteps = (int)(10000*T_ref)+1, n_out = (int)(0.1*T_ref);
vector<const double> wf = {4/9., 1/9., 1/9., 1/9., 1/9., 1/36., 1/36., 1/36., 1/36.};
vector<double> f1(nx*ny*np,0.), f2(nx*ny*np,0.), rho(nx*ny,0.), u(nx*ny,0.), v(nx*ny,0), Feq(np,0.);
double R, U, V, Fy, T, gradx_T, Nu, ftemp;
int newx, newy, id, idn;
// Temperature field
const double alpha = ni/Pr, tauT = alpha*3.+0.5, omegaT = 1./tauT, omegaT1 = 1.-omegaT;
vector<double> g1(nx*ny*np,0.), g2(nx*ny*np,0.), temperature(nx*ny,0.);
///----------------------------------------------------------------------------------------------------------------------------------
///----------------------------------------------------------------------------------------------------------------------------------
void write_fluid_vtk(int time)
{
	const int nz = 1;
	stringstream output_filename;
	output_filename << "vtk_fluid/fluid_t" << time << ".vtk";
	ofstream output_file;
	output_file.open(output_filename.str().c_str());

	output_file << "# vtk DataFile Version 3.0\n";
	output_file << "fluid_state\n";
	output_file << "ASCII\n";
	output_file << "DATASET RECTILINEAR_GRID\n";
	output_file << "DIMENSIONS " << nx << " " << ny << " " << nz << " " << "\n";
	output_file << "X_COORDINATES " << nx << " float\n";
	for(int i = 0; i < nx; ++i)
		output_file << i << " ";
	output_file << "\n";
	output_file << "Y_COORDINATES " << ny  << " float\n";
	for(int j = 0; j < ny ; ++j)
		output_file << j  << " ";
	output_file << "\n";
	output_file << "Z_COORDINATES " << nz  << " float\n";
	for(int z = 0; z < nz ; ++z)
		output_file << z  << " ";
	output_file << "\n";
	output_file << "POINT_DATA " << (nx) * (ny) * (nz)  << "\n";

	output_file << "SCALARS Temperature float 1\n";
	output_file << "LOOKUP_TABLE default\n";
	for(int Z = 0; Z < nz; ++Z)
	for(int Y = 0; Y < ny ; ++Y)
		for(int X = 0; X < nx; ++X)
		{
			id = (X*ny+Y)*nz+Z;
			if(fabs(temperature[id])<1e-16)
				temperature[id] = 0.;
			output_file << temperature[id] << "\n";
		}

	output_file << "VECTORS velocity_vector double\n";
	for(int Z = 0; Z < nz; ++Z)
		for(int Y = 0; Y < ny ; ++Y)
			for(int X = 0; X < nx; ++X)
			{
				id = (X*ny+Y)*nz+Z;
				if(fabs(u[id])<1e-16)
					u[id] = 0.;
				if(fabs(v[id])<1e-16)
					v[id] = 0.;
				output_file << u[id] << " " << v[id] << " 0\n";
			}
	output_file.close();
}
///----------------------------------------------------------------------------------------------------------------------------------
void initial_state()
{
	double NX1 = (double)(nx-1);
	for(int x=0; x<nx; x++)
    for(int y=0; y<ny; y++)
		{
			id = x*ny+y;
			R = rho[id] = rho0;
			U = u[id] = 0.;
			V = v[id] = 0.;
			T = temperature[id] = T_l;
			if(y==0)
				T = temperature[id] = T_h;
			//if(y==ny/2)
				R = rho[id] = rho0*(1.+0.01*cos(2*M_PI*x/NX1));
			for(int k=0; k<np; k++)
			{
				f1[id*np+k] = f2[id*np+k] = wf[k]*R*(1. + 1./cs2*(U*cx[k]+V*cy[k]) + 0.5/cs2/cs2*pow(U*cx[k]+V*cy[k],2) - 0.5/cs2*(U*U+V*V));
				g1[id*np+k] = g2[id*np+k] = wf[k]*T*(1. + 1./cs2*(U*cx[k]+V*cy[k]) + 0.5/cs2/cs2*pow(U*cx[k]+V*cy[k],2) - 0.5/cs2*(U*U+V*V));
			}
		}
}
///----------------------------------------------------------------------------------------------------------------------------------
int algo_LB()
{
	int hh = 0;
	Nu = 0;
	for(int x=0; x<nx; x++)
    for(int y=0; y<ny; y++)
			{
				id = x*ny+y;
				U = V = R = T = 0.;
				for(int k=0; k<np; k++)
				{
					ftemp = f1[id*np+k];
	        R += ftemp;
	        U += ftemp*cx[k];
	        V += ftemp*cy[k];
					T += g1[id*np+k];
				}
				U /= R;
				V /= R;
				temperature[id] = T;
				Fy = R*gbeta*(T-T_0)/DT;
				V += 0.5*Fy/R;
				Nu += V*T;
				u[id] = U;
				v[id] = V;
				rho[id] = R;
				for(int k=0; k<np; k++)
				{
					f1[id*np+k] = omega1*f1[id*np+k] + omega*wf[k]*R*(1. + 1./cs2*(U*cx[k]+V*cy[k]) + 0.5/cs2/cs2*pow(U*cx[k]+V*cy[k],2) - 0.5/cs2*(U*U+V*V));
					f1[id*np+k] += wf[k]*(1.-0.5*omega)*Fy*(3.*(cy[k]-V)+9.*V*cy[k]*cy[k]);
					g1[id*np+k] = omegaT1*g1[id*np+k] + omegaT*wf[k]*T*(1. + 1./cs2*(U*cx[k]+V*cy[k]) + 0.5/cs2/cs2*pow(U*cx[k]+V*cy[k],2) - 0.5/cs2*(U*U+V*V));;
					newx = x+cx[k];
					newy = y+cy[k];
					if(x==0 || x==nx-1)
						newx = (newx+nx)%nx;
					if(y==0 || y==ny-1)
						newy = (newy+ny)%ny;
					idn = newx*ny+newy;
					f2[idn*np+k] = f1[id*np+k];
					g2[idn*np+k] = g1[id*np+k];
				}
			}
	Nu = 1.+Nu/(alpha*DT)/(nx-1);
  return hh;
}
///----------------------------------------------------------------------------------------------------------------------------------
void boundary_conditions()
{
	// south & north
	for(int x=0; x<nx; x++)
		for(int k=0; k<np; k++)
		{
			id = x*ny+0;
			f2[id*np+k] = f1[id*np+opp[k]];
			g2[id*np+k] = wf[k]*T_h;

			id = x*ny+ny-1;
			f2[id*np+k] = f1[id*np+opp[k]];
			g2[id*np+k] = wf[k]*T_l;
		}
}
///----------------------------------------------------------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  system("mkdir vtk_fluid");
  FILE *data = fopen("data.txt","wt");
	initial_state();
	int check_mach = 0;
	clock_t c_start = clock();
	for(int i=0; i<nsteps; i++)
  {
    check_mach = algo_LB();
		boundary_conditions();
		f1 = f2;
		g1 = g2;
    if(check_mach==1)
      goto labelA;
		if(plot_vtk==true && i%n_out==0)
			write_fluid_vtk(i);
   	if(i%100==0)
      printf("Iteration %lf of %lf. Nu=%e\n", (double)(i/T_ref), (double)(nsteps/T_ref), Nu);
  }
	  labelA:
  clock_t c_end = clock();
	double time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;
	cout << "CPU time used: " << time_elapsed_ms << " ms\n";

  return 0;
}
///----------------------------------------------------------------------------------------------------------------------------------
///----------------------------------------------------------------------------------------------------------------------------------
