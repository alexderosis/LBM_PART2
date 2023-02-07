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
///----------------------------------------------------------------------------------
const bool plot_vtk = true;
const int nx = 64, ny = nx, nz = nx, np = 19;
vector<const int> cx = {0, 1, -1, 0,  0, 0,  0, 1, -1,  1, -1, 1, -1,  1, -1, 0,  0,  0,  0},
									cy = {0, 0,  0, 1, -1, 0,  0, 1, -1, -1,  1, 0,  0,  0,  0, 1, -1,  1, -1},
									cz = {0, 0,  0, 0,  0, 1, -1, 0,  0,  0,  0, 1, -1, -1,  1, 1, -1, -1,  1};
const double cs = 1./sqrt(3.), cs2 = 1./3., rho0 = 1.;
const double Ma = 0.1, v0 = Ma*cs, Reynolds = 1600., ni = v0*(nx-1)/Reynolds, tau = ni/cs2+0.5, omega = 1./tau, omega1 = 1.-omega;
const double T_ref = ((double)nx-1)/(v0);
const int nsteps = (int)(10*T_ref)+1, n_out = (int)(T_ref/1)+1;
vector<const double> wf = {1/3., 1/18., 1/18., 1/18., 1/18., 1/18., 1/18., 1/36., 1/36., 1/36., 1/36., 1/36., 1/36., 1/36., 1/36., 1/36., 1/36., 1/36., 1/36.};
vector<double> f1(nx*ny*nz*np,0.), f2(nx*ny*nz*np,0.), u(nx*ny*nz,0.), v(nx*ny*nz,0.), w(nx*ny*nz,0.);
double R, U, V, W, ftemp, feq, kinetic_energy, kinetic_energy_0;
int newx, newy, newz, check, id, idn;
///----------------------------------------------------------------------------------------------------------------------------------
///----------------------------------------------------------------------------------------------------------------------------------
void write_fluid_vtk(int time)
{
	/// Create filename
	stringstream output_filename;
	output_filename << "vtk_fluid/fluid_t" << time << ".vtk";
	ofstream output_file;

	/// Open file
	output_file.open(output_filename.str().c_str());

	/// Write VTK header
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

	/// Write velocity
	output_file << "VECTORS velocity_vector float\n";
	for(int x=0; x<nx; ++x)
		for(int y=0; y<ny; ++y)
			for(int z=0; z<nz; ++z)
			{
				id = (x*ny+y)*nz+z;
				output_file << u[id] << " " << v[id] << " " << w[id] << "\n";
			}
	/// Close file
	output_file.close();
}
///---------------------------------------------------------------------------------------------------------------------------------
void initial_state()
{
	kinetic_energy_0 = 0.;
	double II, JJ, ZZ;
	for(int x=0; x<nx; ++x)
		for(int y=0; y<ny; ++y)
			for(int z=0; z<nz; ++z)
			{
				id = (x*ny+y)*nz+z;
				II = 2*M_PI*(double)x/(double)(nx-1);
      	JJ = 2*M_PI*(double)y/(double)(ny-1);
      	ZZ = 2*M_PI*(double)z/(double)(nz-1);
				R = rho0;
				U = u[id] = v0*cos(II)*sin(JJ)*sin(ZZ);
        V = v[id] = -v0*sin(II)*cos(JJ)*sin(ZZ)*0.5;
        W = w[id] = -v0*sin(II)*sin(JJ)*cos(ZZ)*0.5;
				for(int k=0; k<np;k++)
					f1[id*np+k] = wf[k] * R * (1. + 1./cs2*(U*cx[k]+V*cy[k]+W*cz[k]) + 0.5/cs2/cs2*pow(U*cx[k]+V*cy[k]+W*cz[k],2) - 0.5/cs2*(U*U+V*V+W*W));
				kinetic_energy_0 += R*(U*U+V*V+W*W);
			}
}
///----------------------------------------------------------------------------------------------------------------------------------
int algo_LB()
{
	check = 0.;
	kinetic_energy = 0.;
	for(int x=0; x<nx; x++)
		for(int y=0; y<ny; y++)
			for(int z=0; z<nz; z++)
			{
				id = (x*ny+y)*nz+z;
				U = V = W = R = 0.;
				for(int k=0; k<np; k++)
				{
					ftemp = f1[id*np+k];
					R += ftemp;
					U += ftemp*cx[k];
					V += ftemp*cy[k];
					W += ftemp*cz[k];
				}
				U /= R;
				V /= R;
				W /= R;
				u[id] = U;
				v[id] = V;
				w[id] = W;
				for(int k=0; k<np; k++)
				{
					feq = wf[k]*R*(1. + 1./cs2*(U*cx[k]+V*cy[k]+W*cz[k]) + 0.5/cs2/cs2*pow(U*cx[k]+V*cy[k]+W*cz[k],2) - 0.5/cs2*(U*U+V*V+W*W));
					f1[id*np+k] = f1[id*np+k] + omega*(feq-f1[id*np+k]);
					newx = x+cx[k];
					newy = y+cy[k];
					newz = z+cz[k];
					if(x==0 || x==nx-1)
						newx = (newx+nx)%nx;
					if(y==0 || y==ny-1)
						newy = (newy+ny)%ny;
					if(z==0 || z==nz-1)
						newz = (newz+nz)%nz;
					idn = (newx*ny+newy)*nz+newz;
					f2[idn*np+k] = f1[id*np+k];
				}
				kinetic_energy += R*(U*U+V*V+W*W);
				if(isnan(U) || isinf(U))
					check = 1;
			}
	return check;
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
    f1 = f2;
    if(check_mach==1)
      goto labelA;
		 if(plot_vtk==true && i%n_out==0)
		 	write_fluid_vtk(i);
		fprintf(data,"%lf    %e\n", (double)i/T_ref, kinetic_energy/kinetic_energy_0);
   	 if(i%1==0)
       printf("Time %lf of %lf. Energy=%e\n", (double)i/T_ref, (double)nsteps/T_ref, kinetic_energy/kinetic_energy_0 );
  }
  labelA:
	clock_t c_end = clock();
	double time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;
	cout << "CPU time used: " << time_elapsed_ms << " ms\n";
  fclose(data);
  return 0;
}
///----------------------------------------------------------------------------------------------------------------------------------
///----------------------------------------------------------------------------------------------------------------------------------
