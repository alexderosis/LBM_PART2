///----------------------------------------------------------------------------------------------------------------------------------
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cassert>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <algorithm>
///----------------------------------------------------------------------------------------------------------------------------------
using namespace std;
const bool plot_vtk = true;
/// Flow quantities
const bool cms_hydro = true, cms_phase = true;
const int D_wedge = 500, radius = D_wedge/2, nx = 4*D_wedge, ny = 2*D_wedge, np = 9, yc = 3*ny/4;
const double beta_angle = 4.*M_PI/180., cs2 = 1./3., cs4 = cs2*cs2, cs6 = cs4*cs2, cs8 = cs6*cs2;
const double xi = 4., v0 = -0.01, Scale_velocity = 1./(-v0), Scale_length = 1/((double)D_wedge), Scale_time = Scale_length/Scale_velocity, Scale_acc = Scale_length/(Scale_time*Scale_time),
						gravity = 0./Scale_acc, nu_phys = 1E-6, Scale_nu = Scale_length*Scale_length/Scale_time, nu = nu_phys/Scale_nu, Reynolds = -v0*D_wedge/nu, T = 0.0061/Scale_time;
const int nsteps = (int)(T+1), n_out = (int)(0.001/Scale_time);
const double rhoL = 1., rhoH = 100., niH = nu, tauH = niH/cs2, niL = nu*15., tauL = niL/cs2, sigma = 1E-5, Scale_rho = 1000./rhoH, Scale_PressRes = Scale_rho*pow(Scale_velocity,2)*Scale_length;
// HYDRO            		0   1   2   3   4   5  6  7  8
const vector<int> cx = {0, 1, 0, -1,  0, 1, -1, -1,  1},
									cy = {0, 0, 1,  0, -1, 1,  1, -1, -1},
         				 opp = {0, 3, 4, 1, 2, 7, 8, 5, 6};
const vector<double> wf = {4./9., 1./9., 1./9., 1./9., 1./9., 1./36., 1./36., 1./36., 1./36.};
double A, B, C, R, U, V, ftemp, CX, CY, U2, V2, U2V2, UV, U2V, UV2, U0, V0, Press0;
double Press, tau, omega, omega1, feq, ni;
vector<double> temp_pop(np, 0.), f1(nx*ny*np, 0.), f2(nx*ny*np, 0.), press(nx*ny, 0.), u(nx*ny, 0.), v(nx*ny, 0.), rho(nx*ny, 0.), u_old(nx*ny, 0.), v_old(nx*ny, 0.);
double Fpx, Fpy, Fsx, Fsy, Fm, Fmx, Fmy, mu;
double gradx_u, grady_u, gradx_v, grady_v;
int newi, newj, check;
double k1, k2, k3, k4, k5, k6, k7, k8;
double r1, r2, r3, r4, r5, r6, r7, r8;
double third_order, fourth_order;
// PHASE
const double PhiH = 1., PhiL = 0., Phi0 = 0.5*(PhiH+PhiL);
const double M = 0.25, tau_phase = M/cs2+0.5, omega_phase = 1./tau_phase, omega_phase1 = 1.-omega_phase;
const double beta = 12.*sigma/xi, kappa = 3.*sigma*xi/2.;
double grad_phix, grad_phiy, laplPhi;
double Phi, Phi_prev, geq;
vector<double> g1(nx*ny*np, 0.), g2(nx*ny*np, 0.), temp_pop_phase(np, 0.), phase(nx*ny, 0.), phase_old(nx*ny, 0.);
double Fx, Fy, gtemp, P, Fx_phase, Fy_phase;
double k1_g, k2_g, k3_g, k4_g, k5_g, k6_g, k7_g, k8_g;
double Gamma, Nx, Ny;
double Phii, thisPhi;
int newii, newjj, id, idn;
double gradx_of_u, gradx_of_v, grady_of_u, grady_of_v;
vector<double> collision_operator(np, 0.), collision_operator_phase(np, 0.);
// IBM
const double ds = 0.5;
const int N = (int)(2*D_wedge/ds+1);
vector<double> X(N, 0.), Y(N, 0.), X0(N, 0.), Y0(N, 0.), G1(2*N, 0.), gx(nx*ny, 0.), gy(nx*ny, 0.);
double Force_body_x, Force_body_y;
vector<double> f_star(np, 0.), df(N*np, 0.);
double wet_area, wet_arean, buoyancy_force, buoyancy_forcen, F0y, F0yn, UY[nsteps], DY, Y_center;
const double rho_cyl = 0.62*rhoH, mass = rho_cyl*M_PI*radius*radius, gravity_force = mass*gravity;
double press_wedge[N], Pressure_resultant;
const double rho0 = rhoL;
int wet_surface;
///---------------------------------------------------------------------------------
///---------------------------------------------------------------------------------
void bilinear_interpolation(int index)
{
	double x, y, x1, x2, y1, y2, R1, R2;
	x = X[index];
	y = Y[index]-1.5;
	//id = (int)(x*ny+(y-1));
	//if(phase[id]>0.8)
	{
		x1 = floor(x);
		y1 = floor(y);
		x2 = x1+1;
		y2 = y1+1;
		R1 = (x2-x)/1*(press[(int)x1*ny+(int)y1])*rhoH*cs2 + (x-x1)*(press[(int)x2*ny+(int)y1])*rhoH*cs2;
		R2 = (x2-x)/1*(press[(int)x1*ny+(int)y2])*rhoH*cs2 + (x-x1)*(press[(int)x2*ny+(int)y2])*rhoH*cs2;
		press_wedge[index] = (y2-y)*R1+(y-y1)*R2;

	press_wedge[index] = press[(int)x*ny+(int)y]*rho[(int)x*ny+(int)y];
	 }
	 //else
	// press_wedge[index] = 0.;
	Pressure_resultant += press_wedge[index];
}
///---------------------------------------------------------------------------------
void write_pressure(int time)
{
	int xx, yy;
	stringstream output_filename;
	output_filename << "presst" << time << ".txt";
	ofstream output_file;
	output_file.open(output_filename.str().c_str());
	for(int k=0; k<N; k++)
	{
		//bilinear_interpolation(k);
		output_file << k << " " << press_wedge[k] << "\n";
	}
	output_file.close();
}
///----------------------------------------------------------------------------------------------------------------------------------
///----------------------------------------------------------------------------------------------------------------------------------
void diagnostics(int l)
{
	double wet_surface_vk = -DY/tan(beta_angle),
	wet_surface_w = wet_surface_vk*M_PI*0.5,
	vk = rhoH*wet_surface_vk*M_PI*0.5*pow(v0,2)/tan(beta_angle),
	wagner = rhoH*wet_surface_w*M_PI*M_PI*0.25*pow(v0,2)/tan(beta_angle),
	err_vk = fabs(vk-Pressure_resultant)/fabs(vk)*100.,
	err_w = fabs(wagner-Pressure_resultant)/fabs(wagner)*100.;
    if(l%1==0)
    {
    	printf("t %lf of %lf. %e %e %e\n", l*Scale_time, nsteps*Scale_time, Pressure_resultant*Scale_PressRes, wagner*Scale_PressRes, err_w);
			//printf("t %lf of %lf. %e %e %e\n", l*Scale_time, nsteps*Scale_time, DY, tan(beta_angle), wet_surface_vk);
    }
}
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
	output_file << "DIMENSIONS " << nx << " " << ny << " 1" << "\n";
	output_file << "X_COORDINATES " << nx << " double\n";
	for(int i = 0; i < nx; ++i)
		output_file << i << " ";
	output_file << "\n";
	output_file << "Y_COORDINATES " << ny  << " double\n";
	for(int j = 0; j < ny ; ++j)
		output_file << j  << " ";
	output_file << "\n";
	output_file << "Z_COORDINATES " << 1 << " double\n";
	output_file << 0 << "\n";
	output_file << "POINT_DATA " << (nx) * (ny) << "\n";

	output_file << "SCALARS phase double 1\n";
	output_file << "LOOKUP_TABLE default\n";
	for(int Y = 0; Y < ny ; ++Y)
		for(int X = 0; X < nx; ++X)
		{
			if(phase[X*ny+Y]>1E-10)
				output_file << phase[X*ny+Y]<< "\n";
			else
				output_file << 0 << "\n";
			}

  output_file << "SCALARS density double 1\n";
  output_file << "LOOKUP_TABLE default\n";
  for(int Y = 0; Y < ny ; ++Y)
    for(int X = 0; X < nx; ++X)
    	output_file << rho[X*ny+Y]<< "\n";

			output_file << "SCALARS pressure double 1\n";
		  output_file << "LOOKUP_TABLE default\n";
		  for(int Y = 0; Y < ny ; ++Y)
		    for(int X = 0; X < nx; ++X)
		    	output_file << press[X*ny+Y]<< "\n";

  output_file << "VECTORS velocity_vector double\n";
  for(int Y = 0; Y < ny ; ++Y)
  	for(int X = 0; X < nx; ++X)
  		output_file << u[X*ny+Y] << " " << v[X*ny+Y] << " 0\n";

	/// Close file
	output_file.close();
}
///---------------------------------------------------------------------------------
void write_cylinder_vtk(int time)
{
	/// Create filename
	stringstream output_filename;
	output_filename << "vtk_solid/solid_t" << time << ".vtk";
	ofstream output_file;

	/// Open file
	output_file.open(output_filename.str().c_str());

	/// Write VTK header
	output_file << "# vtk DataFile Version 3.0\n";
	output_file << "particle_state\n";
	output_file << "ASCII\n";
	output_file << "DATASET POLYDATA\n";

	/// Write node positions
	output_file << "POINTS " << N << " double\n";
	for(int n = 0; n < N; ++n)
	    output_file << X[n] << " " << Y[n] << " 0\n";

	/// Write lines between neighboring nodes
	output_file << "LINES " << N-1 << " " << 3 * (N-1) << "\n";
	for(int n = 0; n < N-1; ++n)
	    output_file << "2 " << n << " " << (n + 1) % N << "\n";

	/// Write vertices
	output_file << "VERTICES 1 " << N + 1 << "\n";
	output_file << N << " ";
	for(int n = 0; n < N; ++n)
	    output_file << n << " ";

	/// Close file
	output_file.close();
}
///----------------------------------------------------------------------------------------------------------------------------------
void update_position(int tempo)
{
	tempo++;
// 	wet_arean = wet_area;
// 	double depth;
// 	if(*min_element(Y.begin(),Y.end())>ny/2)
// 		wet_area = 0.;
// 	else
// 	{
// 		depth = (double)(ny/2-*min_element(Y.begin(),Y.end()));
// 		wet_area = pow(radius,2) * acos((radius-depth)/radius) - (radius-depth)*sqrt(2.*radius*depth-pow(depth,2));
// 	}
//
// 	buoyancy_force = rhoH*gravity*wet_area,
//   buoyancy_forcen = rhoH*gravity*wet_arean;
// //buoyancy_force = buoyancy_forcen = 0.;
	// F0y += -72.5*Scale_length/(1000*pow(Scale_length,3));
	// F0yn += -72.5*Scale_length/(1000*pow(Scale_length,3));
	UY[tempo] = v0;
	DY += UY[tempo];
	//UY[tempo] = UY[tempo-1]+(1.5*F0y-.5*F0yn)/mass;
	for(int k=0; k<N; k++)
		Y[k] = Y0[k]+DY;
	Y_center = 0.5*(*min_element(Y.begin(),Y.end()) + *max_element(Y.begin(),Y.end()));
}
///----------------------------------------------------------------------------------------------------------------------------------
void define_obstacle()
{
	double K;
	for(int k=0; k<N; k++)
	{
		X0[k] = X[k] = nx/2-D_wedge+k*0.5;
		Y0[k] = Y[k] = yc+fabs((X0[k]-nx/2)*tan(beta_angle));
	}
	Y_center = 0.5*(*min_element(Y.begin(),Y.end()) + *max_element(Y.begin(),Y.end()));
}
//----------------------------------------------------------------------------------------------------------------------------------
double wIB(double r)
{
	double wib;
	r = fabs(r);
	if(r<=1.)
		wib = 1./8.*(3.-2.*r+sqrt(1.+4.*r-4.*r*r));
	else if(1.<r && r<=2.)
		wib = 1./8.*(5.-2.*r-sqrt(-7.+12.*r-4.*r*r));
	else
		wib = 0.;
	return wib;
}
//----------------------------------------------------------------------------------------------------------------------------------
int immersed_boundary(int mn)
{
	double eps = 0, dot, toll = pow(10.,-4), minX, minY, maxX, maxY, tempx, tempy, temp, tempo = (double)mn, Ux, Uy;
  int l = 0;
	Ux = 0.;
	Uy = UY[mn];
	std::fill(G1.begin(), G1.end(), 0.);
	for(int k=0; k<N; k++)
	{
		minX = X[k]-3;
		maxX = X[k]+3;
		minY = Y[k]-3;
		maxY = Y[k]+3;
		tempx = tempy = 0.;
		for(int i=minX; i<maxX; i++)
			for(int j=minY; j<maxY; j++)
			//if((i-nx/2)*(i-nx/2) + (j-Y_center)*(j-Y_center)>radius*radius)
			{
				id = i*ny+j;
				temp = wIB(i-X[k])*wIB(j-Y[k]);
				tempx += u[id]*temp;
				tempy += v[id]*temp;
			}
		G1[k] = Ux-tempx;
		G1[k+N] = Uy-tempy;
		eps += Ux*Ux+Uy*Uy;
	}
	eps = sqrt(eps);
	std::fill(gx.begin(), gx.end(), 0.);
	std::fill(gy.begin(), gy.end(), 0.);
	do
	{
		l++;
		dot = Force_body_x = Force_body_y = 0.;
		minX = *min_element(X.begin(), X.end())-3;
		maxX = *max_element(X.begin(), X.end())+3;
		minY = *min_element(Y.begin(), Y.end())-3;
		maxY = *max_element(Y.begin(), Y.end())+3;
		for(int i=minX; i<maxX; i++)
			for(int j=minY; j<maxY; j++)
			//if((i-nx/2)*(i-nx/2) + (j-Y_center)*(j-Y_center)>radius*radius)
			{
				tempx = tempy = 0.;
				for(int k=0; k<N; k++)
				{
					temp = wIB(i-X[k])*wIB(j-Y[k])*ds;
					tempx += G1[k]*temp;
					tempy += G1[k+N]*temp;
				}
				id = i*ny+j;
				gx[id] = tempx;
				gy[id] = tempy;
				Force_body_x -= rho[id]*tempx;
				Force_body_y -= rho[id]*tempy;
			}

		for(int k=0; k<N; k++)
		{
			minX = X[k]-3;
			maxX = X[k]+3;
			minY = Y[k]-3;
			maxY = Y[k]+3;
			tempx = tempy = 0.;
      for(int i=minX; i<maxX; i++)
        for(int j=minY; j<maxY; j++)
				//if((i-nx/2)*(i-nx/2) + (j-Y_center)*(j-Y_center)>radius*radius)
				{
					id = i*ny+j;
					temp = wIB(i-X[k])*wIB(j-Y[k]);
					tempx += (u[id]+gx[id])*temp;
					tempy += (v[id]+gy[id])*temp;
				}
			G1[k] += Ux-tempx;
			G1[k+N] += Uy-tempy;
			dot += pow(Ux-tempx,2)+pow(Uy-tempy,2);
		}
    dot = sqrt(dot);
	}while(dot>=eps*toll+1e-4 && l<1);
	return l;
}
//----------------------------------------------------------------------------------------------------------------------------------
void immersed_boundary2(int mn)
{
	double minX, minY, maxX, maxY, temp, coeff, kernel;
	Force_body_x = Force_body_y = 0.;
	for(int k=0; k<N; k++)
	{
		minX = X[k]-3;
		maxX = X[k]+3;
		minY = Y[k]-3;
		maxY = Y[k]+3;
		coeff = 0.;
		std::fill(f_star.begin(), f_star.end(), 0.);
		for(int i=minX; i<maxX; i++)
			for(int j=minY; j<maxY; j++)
			{
				id = i*ny+j;
				kernel = wIB(i-X[k])*wIB(j-Y[k]);
				coeff += 2.*ds*rho[id]*kernel*kernel;
				for(int n=0; n<np; n++)
					f_star[n] += f1[id*np+n]*kernel;
			}
		coeff = 1./coeff;
		id = (int)(X[k]*ny+Y[k]);
		for(int n=0; n<np; n++)
		{
			df[k*np+n] = coeff*(f_star[opp[n]]+6.*wf[n]*cy[n]*UY[mn]-f_star[n]);
			Force_body_x -= df[k*np+n]*cx[n];
			Force_body_y -= df[k*np+n]*cy[n];
		}
	}
	minX = *min_element(X.begin(), X.end())-3;
	maxX = *max_element(X.begin(), X.end())+3;
	minY = *min_element(Y.begin(), Y.end())-3;
	maxY = *max_element(Y.begin(), Y.end())+3;
	for(int n=0; n<np; n++)
		for(int i=minX; i<maxX; i++)
			for(int j=minY; j<maxY; j++)
			{
				id = i*ny+j;
				temp = 0;
				for(int k=0; k<N; k++)
					temp += df[k*np+n]*wIB(i-X[k])*wIB(j-Y[k])*ds;
				f1[id*np+n] += temp;
			}
}
///----------------------------------------------------------------------------------------------------------------------------------
void compute_grad_phi(int i, int j)
{
	grad_phix = grad_phiy = laplPhi = 0.;
	gradx_of_u = gradx_of_v = grady_of_u = grady_of_v = 0.;
  thisPhi = phase_old[i*ny+j];
	for(int k=1; k<np; k++)
	{
		newii = i+cx[k];
		newjj = j+cy[k];
		if(i==0 || i==nx-1)
			newii = (newii+nx)%nx;
		if(j==0 || j==ny-1)
			newjj = (newjj+ny)%ny;
		Phii = phase_old[newii*ny+newjj];
		grad_phix += Phii*wf[k]*cx[k];
		grad_phiy += Phii*wf[k]*cy[k];
    laplPhi += wf[k]*(Phii-thisPhi);

		gradx_of_u += u[newii*ny+newjj]*wf[k]*cx[k];
		gradx_of_v += v[newii*ny+newjj]*wf[k]*cx[k];
		grady_of_u += u[newii*ny+newjj]*wf[k]*cy[k];
		grady_of_v += v[newii*ny+newjj]*wf[k]*cy[k];
	}
	grad_phix /= cs2;
	grad_phiy /= cs2;
	gradx_of_u /= cs2;
	gradx_of_v /= cs2;
	grady_of_u /= cs2;
	grady_of_v /= cs2;
  laplPhi *= 2./cs2;

	//if(i==0 || i==nx-1)
		//grad_phix = grad_phiy = laplPhi = 0.;
	if(j==0 || j==ny-1)
		grad_phix = grad_phiy = laplPhi = 0.;
	//if((i-nx/2)*(i-nx/2) + (j-Y_center)*(j-Y_center)<radius*radius)
		//grad_phix = grad_phiy = laplPhi = 0.;
}
///----------------------------------------------------------------------------------------------------------------------------------
void initial_state()
{
  double extra_term;
	for(int i=0; i<nx; i++)
    for(int j=0; j<ny; j++)
		{
			id = i*ny+j;
			if(j<yc)
				phase[id] = 1.;
      else
				phase[id] = 0.;
			gx[id] = gy[id] = 0.;
      rho[id] = rhoL + phase[id]*(rhoH-rhoL);
		}
	for(int i=0; i<nx; i++)
    for(int j=0; j<ny; j++)
    {
			id = i*ny+j;
			u[id] = U = 0.;
			v[id] = V = 0.;
      R = rho[id];
      Phi = phase[id];
		  compute_grad_phi(i,j);
      Nx  = grad_phix/(sqrt(pow(grad_phix,2)+pow(grad_phiy,2))+1E-12);
      Ny  = grad_phiy/(sqrt(pow(grad_phix,2)+pow(grad_phiy,2))+1E-12);

      Fpx = -Press*cs2*(rhoH-rhoL)/(PhiH-PhiL)*grad_phix;
      Fpy = -Press*cs2*(rhoH-rhoL)/(PhiH-PhiL)*grad_phiy;
      mu = 4.*beta*(Phi-PhiL)*(Phi-PhiH)*(Phi-Phi0)-kappa*laplPhi;
      Fsx = mu*grad_phix;
      Fsy = mu*grad_phiy;
      Fx = Fpx+Fsx;
      Fy = Fpy+Fsy;
      for(int k=0; k<np;k++)
			{
        A = U*cx[k]+V*cy[k];
        third_order = 1./(2.*cs6)*((cx[k]*cx[k]-cs2)*cy[k]*U*U*V+(cy[k]*cy[k]-cs2)*cx[k]*U*V*V);
        fourth_order = 1./(4.*cs8)*((cx[k]*cx[k]-cs2)*(cy[k]*cy[k]-cs2)*U*U*V*V);
        extra_term = 4.*Phi*(1.-Phi)/xi*wf[k]*(cx[k]*Nx+cy[k]*Ny);
        geq = Phi*wf[k]*(1.+3.*A+4.5*((cx[k]*cx[k]-cs2)*U*U+(cy[k]*cy[k]-cs2)*V*V+
                2.*cx[k]*cy[k]*U*V)+third_order+fourth_order) -
              0.5*extra_term;
        g1[id*np+k] = g2[id*np+k] = geq+extra_term;

        feq = wf[k]*(Press+3.*A+4.5*((cx[k]*cx[k]-cs2)*U*U+(cy[k]*cy[k]-cs2)*V*V+
                   2.*cx[k]*cy[k]*U*V)+third_order+fourth_order)-
                   0.5*wf[k]/(R*cs2)*(cx[k]*Fx+cy[k]*Fy);
        f1[id*np+k] = f2[id*np+k] = feq+wf[k]/(R*cs2)*(cx[k]*Fx+cy[k]*Fy);
			}
		}
}
///----------------------------------------------------------------------------------------------------------------------------------
int algorithm_lattice_boltzmann(int tempo)
{
  double third_order, fourth_order, extra_term;
	int check = 0;
	u_old = u;
	v_old = v;
	phase_old = phase;
	for(int i=0; i<nx; i++)
		for(int j=0; j<ny; j++)
		{
			P = 0.;
			id = i*ny+j;
			for(int k=0; k<np; k++)
				P += g1[id*np+k];
			phase[id] = P;
		}
	for(int i=0; i<nx; i++)
		for(int j=0; j<ny; j++)
		{
			id = i*ny+j;
			compute_grad_phi(i,j);
			U0 = u[id];
			V0 = v[id];
			Phi_prev = phase[id];
			Press0 = press[id];
			U = V = Press = Phi = 0.;
			for(int k=0; k<np; k++)
      {
        temp_pop[k] = ftemp = f1[id*np+k];
        temp_pop_phase[k] = g1[id*np+k];
        Press += ftemp;
        U += ftemp*cx[k];
        V += ftemp*cy[k];
      }
			Phi = phase[id];
      press[id] = Press;
      Fpx = -Press*cs2*(rhoH-rhoL)*grad_phix;
      Fpy = -Press*cs2*(rhoH-rhoL)*grad_phiy;
      mu = 1.5*sigma*(32.*Phi*(Phi-1.)*(Phi-0.5)/xi-xi*laplPhi);
			//mu = 4.*beta*(Phi-PhiL)*(Phi-PhiH)*(Phi-Phi0)-kappa*laplPhi;
      Fsx = mu*grad_phix;
      Fsy = mu*grad_phiy;
      Fmx = Fmy = 0.;
      // for(int k=0; k<np; k++)
      // {
			// 	A = U0*cx[k]+V0*cy[k];
      //   feq = wf[k]*(Press0+3.*A+4.5*((cx[k]*cx[k]-cs2)*U0*U0+(cy[k]*cy[k]-cs2)*V0*V0+
      //            2.*cx[k]*cy[k]*U0*V0));
      //   Fmx += cx[k]*cy[k]*(ftemp-feq);
      //   Fmy += cy[k]*cx[k]*(ftemp-feq);
      // }
			phase[id] = Phi;
      press[id] = Press;
      Fpx = -Press0*cs2*(rhoH-rhoL)/(PhiH-PhiL)*grad_phix;
      Fpy = -Press0*cs2*(rhoH-rhoL)/(PhiH-PhiL)*grad_phiy;
      mu = 4.*beta*(Phi_prev-PhiL)*(Phi_prev-PhiH)*(Phi_prev-Phi0)-kappa*laplPhi;
      Fsx = mu*grad_phix;
      Fsy = mu*grad_phiy;
      tau = tauL+(Phi_prev-PhiL)/(PhiH-PhiL)*(tauH-tauL);
      ni = tau*cs2;
      omega = 1./(tau+0.5);
      omega1 = 1.-omega;
      //Fmx *= -ni/((tau+0.5)*cs2)*(rhoH-rhoL)/(PhiH-PhiL)*grad_phix;
      //Fmy *= -ni/((tau+0.5)*cs2)*(rhoH-rhoL)/(PhiH-PhiL)*grad_phiy;
			Fmx = ni* (gradx_of_u*(rhoH-rhoL)/(PhiH-PhiL)*grad_phix + gradx_of_v*(rhoH-rhoL)/(PhiH-PhiL)*grad_phiy);
			Fmy = ni* (grady_of_u*(rhoH-rhoL)/(PhiH-PhiL)*grad_phix + grady_of_v*(rhoH-rhoL)/(PhiH-PhiL)*grad_phiy);
			rho[id] = R = rhoL+(Phi_prev-PhiL)/(PhiH-PhiL)*(rhoH-rhoL);
      Fx = Fpx+Fsx+Fmx;// + gx[id]*rho[id];
      Fy = Fpy+Fsy+Fmy;// + gy[id]*rho[id];
			if(phase[id]>0.5)
			{
				Fx += gx[id]*rho[id];
	      Fy += gy[id]*rho[id];
			}
      //Fy += -(R - 0.5*(rhoH+rhoL))*gravity;
			Fy += -R*gravity;
      U += 0.5*Fx/R;
      V += 0.5*Fy/R;
			// if((i-nx/2)*(i-nx/2) + (j-Y_center)*(j-Y_center)<radius*radius)
			// {
			// 	U = Fx = Fy = 0.;
			// 	V = UY[tempo];
			// }
    	u[id] = U;
			v[id] = V;
      Nx  = grad_phix/(sqrt(pow(grad_phix,2)+pow(grad_phiy,2))+1E-12);
      Ny  = grad_phiy/(sqrt(pow(grad_phix,2)+pow(grad_phiy,2))+1E-12);
			U2 = U*U;
			V2 = V*V;
			UV = U*V;
			U2V = U2*V;
			UV2 = U*V2;
			U2V2 = U2*V2;
      if(cms_hydro)
      {
  			k1 = k2 = k3 = k4 = k5 = k6 = k7 = k8 = 0.;
  			for(int k=0; k<np; k++)
  			{
  				CX = cx[k]-U;
  				CY = cy[k]-V;
  				ftemp = temp_pop[k];
  				k4 += ftemp*(CX*CX-CY*CY);
  				k5 += ftemp*CX*CY;
  			}

        k1 = Fx/(2.*R) - U*(Press - 1.);
        k2 = Fy/(2.*R) - V*(Press - 1.);
        k3 = (2.*Press)/3. + Press*U2 + Press*V2 - U2 - V2;
        k4 = omega1*k4+omega*(U2 - V2)*(Press - 1.);
        k5 = omega1*k5+omega*U*V*(Press - 1.);
        k6 = Fy/(6.*R) - (V*(3.*U2 + 1.)*(Press - 1.))/3.;
        k7 = Fx/(6.*R) - (U*(3.*V2 + 1.)*(Press - 1.))/3.;
        k8 = Press/9. + (Press*U2)/3. + (Press*V2)/3. - U2/3. - V2/3. - U2*V2 + Press*U2*V2;

        /// reconstruct post-collision populations
        collision_operator[0] = (U2V2-U2-V2+1.)*Press+(2.*U*V2-2.*U)*k1+(2.*V*U2-2.*V)*k2+(U2/2.+V2/2.-1.)*k3+(V2/2.-U2/2.)*k4+4.*U*V*k5+2.*V*k6+2.*U*k7+k8;
        collision_operator[1] = (U2/2.-(U2*V2)/2.-(U*V2)/2.+U/2.)*Press+(U-U*V2-V2/2.+1/2.)*k1+(-V*U2-V*U)*k2+(-U2/4.-U/4.-V2/4.+1/4.)*k3+(U2/4.+U/4.-V2/4.+1/4.)*k4+(-V-2.*UV)*k5+(-V)*k6+(-U-1/2.)*k7-k8/2.;
        collision_operator[2] = (V2/2.-(U2*V)/2.-(U2*V2)/2.+V/2.)*Press+(-U*V2-UV)*k1+(V-U2*V-U2/2.+1/2.)*k2+(-U2/4.-V2/4.-V/4.+1/4.)*k3+(U2/4.-V2/4.-V/4.-1/4.)*k4+(-U-2.*UV)*k5+(-V-1/2.)*k6+(-U)*k7-k8/2.;
        collision_operator[3] = (-(U2V2)/2.+U2/2.+(U*V2)/2.-U/2.)*Press+(U-U*V2+V2/2.-1/2.)*k1+(-V*U2+V*U)*k2+(-U2/4.+U/4.-V2/4.+1/4.)*k3+(U2/4.-U/4.-V2/4.+1/4.)*k4+(V-2.*UV)*k5+(-V)*k6+(1/2.-U)*k7-k8/2.;
        collision_operator[4] = (-(U2V2)/2.+(U2*V)/2.+V2/2.-V/2.)*Press+(-U*V2+UV)*k1+(V-U2*V+U2/2.-1/2.)*k2+(-U2/4.-V2/4.+V/4.+1/4.)*k3+(U2/4.-V2/4.+V/4.-1/4.)*k4+(U-2.*UV)*k5+(1/2.-V)*k6+(-U)*k7-k8/2.;
        collision_operator[5] = ((U2V2)/4.+(U2*V)/4.+(U*V2)/4.+(UV)/4.)*Press+(V/4.+(UV)/2.+(U*V2)/2.+V2/4.)*k1+(U/4.+(UV)/2.+(U2*V)/2.+U2/4.)*k2+(U2/8.+U/8.+V2/8.+V/8.)*k3+(-U2/8.-U/8.+V2/8.+V/8.)*k4+(U/2.+V/2.+U*V+1/4.)*k5+(V/2.+1/4.)*k6+(U/2.+1/4.)*k7+k8/4.;
        collision_operator[6] = ((U2V2)/4.+(U2*V)/4.-(U*V2)/4.-(UV)/4.)*Press+((UV)/2.-V/4.+(U*V2)/2.-V2/4.)*k1+((U2*V)/2.-(UV)/2.-U/4.+U2/4.)*k2+(U2/8.-U/8.+V2/8.+V/8.)*k3+(-U2/8.+U/8.+V2/8.+V/8.)*k4+(U/2.-V/2.+U*V-1/4.)*k5+(V/2.+1/4.)*k6+(U/2.-1/4.)*k7+k8/4.;
        collision_operator[7] = ((U2V2)/4.-(U2*V)/4.-(U*V2)/4.+(UV)/4.)*Press+(V/4.-(UV)/2.+(U*V2)/2.-V2/4.)*k1+(U/4.-(UV)/2.+(U2*V)/2.-U2/4.)*k2+(U2/8.-U/8.+V2/8.-V/8.)*k3+(-U2/8.+U/8.+V2/8.-V/8.)*k4+(U*V-V/2.-U/2.+1/4.)*k5+(V/2.-1/4.)*k6+(U/2.-1/4.)*k7+k8/4.;
        collision_operator[8] = ((U2V2)/4.-(U2*V)/4.+(U*V2)/4.-(UV)/4.)*Press+((U*V2)/2.-(UV)/2.-V/4.+V2/4.)*k1+((UV)/2.-U/4.+(U2*V)/2.-U2/4.)*k2+(U2/8.+U/8.+V2/8.-V/8.)*k3+(-U2/8.-U/8.+V2/8.-V/8.)*k4+(V/2.-U/2.+U*V-1/4.)*k5+(V/2.-1/4.)*k6+(U/2.+1/4.)*k7+k8/4.;
      }
      if(cms_phase)
      {
  			k1 = k2 = k3 = k4 = k5 = k6 = k7 = k8 = 0.;
  			for(int k=0; k<np; k++)
  			{
  				CX = cx[k]-U;
  				CY = cy[k]-V;
  				ftemp = temp_pop_phase[k];
          k1 += ftemp*CX;
          k2 += ftemp*CY;
  			}
        /// post-collision state
        k1 = omega_phase1*k1;
        k2 = omega_phase1*k2;
  			k3 = 2.*Phi/3.;
  			k4 = 0.;
  			k5 = 0.;
  			k6 = 0.;
  			k7 = 0.;
  			k8 = Phi/9.;
        Fx_phase = cs2*(1.-4.*pow((Phi-Phi0),2))/xi*Nx;
        Fy_phase = cs2*(1.-4.*pow((Phi-Phi0),2))/xi*Ny;
        k1 += 0.5*Fx_phase;
        k2 += 0.5*Fy_phase;
        k6 += 0.5*cs2*Fy_phase;
        k7 += 0.5*cs2*Fx_phase;
        /// reconstruct post-collision populations
        collision_operator_phase[0] = (U2V2-U2-V2+1.)*Phi+(2.*U*V2-2.*U)*k1+(2.*V*U2-2.*V)*k2+(U2/2.+V2/2.-1.)*k3+(V2/2.-U2/2.)*k4+4.*U*V*k5+2.*V*k6+2.*U*k7+k8;
        collision_operator_phase[1] = (U2/2.-(U2*V2)/2.-(U*V2)/2.+U/2.)*Phi+(U-U*V2-V2/2.+1/2.)*k1+(-V*U2-V*U)*k2+(-U2/4.-U/4.-V2/4.+1/4.)*k3+(U2/4.+U/4.-V2/4.+1/4.)*k4+(-V-2.*UV)*k5+(-V)*k6+(-U-1/2.)*k7-k8/2.;
        collision_operator_phase[2] = (V2/2.-(U2*V)/2.-(U2*V2)/2.+V/2.)*Phi+(-U*V2-UV)*k1+(V-U2*V-U2/2.+1/2.)*k2+(-U2/4.-V2/4.-V/4.+1/4.)*k3+(U2/4.-V2/4.-V/4.-1/4.)*k4+(-U-2.*UV)*k5+(-V-1/2.)*k6+(-U)*k7-k8/2.;
        collision_operator_phase[3] = (-(U2V2)/2.+U2/2.+(U*V2)/2.-U/2.)*Phi+(U-U*V2+V2/2.-1/2.)*k1+(-V*U2+V*U)*k2+(-U2/4.+U/4.-V2/4.+1/4.)*k3+(U2/4.-U/4.-V2/4.+1/4.)*k4+(V-2.*UV)*k5+(-V)*k6+(1/2.-U)*k7-k8/2.;
        collision_operator_phase[4] = (-(U2V2)/2.+(U2*V)/2.+V2/2.-V/2.)*Phi+(-U*V2+UV)*k1+(V-U2*V+U2/2.-1/2.)*k2+(-U2/4.-V2/4.+V/4.+1/4.)*k3+(U2/4.-V2/4.+V/4.-1/4.)*k4+(U-2.*UV)*k5+(1/2.-V)*k6+(-U)*k7-k8/2.;
        collision_operator_phase[5] = ((U2V2)/4.+(U2*V)/4.+(U*V2)/4.+(UV)/4.)*Phi+(V/4.+(UV)/2.+(U*V2)/2.+V2/4.)*k1+(U/4.+(UV)/2.+(U2*V)/2.+U2/4.)*k2+(U2/8.+U/8.+V2/8.+V/8.)*k3+(-U2/8.-U/8.+V2/8.+V/8.)*k4+(U/2.+V/2.+U*V+1/4.)*k5+(V/2.+1/4.)*k6+(U/2.+1/4.)*k7+k8/4.;
        collision_operator_phase[6] = ((U2V2)/4.+(U2*V)/4.-(U*V2)/4.-(UV)/4.)*Phi+((UV)/2.-V/4.+(U*V2)/2.-V2/4.)*k1+((U2*V)/2.-(UV)/2.-U/4.+U2/4.)*k2+(U2/8.-U/8.+V2/8.+V/8.)*k3+(-U2/8.+U/8.+V2/8.+V/8.)*k4+(U/2.-V/2.+U*V-1/4.)*k5+(V/2.+1/4.)*k6+(U/2.-1/4.)*k7+k8/4.;
        collision_operator_phase[7] = ((U2V2)/4.-(U2*V)/4.-(U*V2)/4.+(UV)/4.)*Phi+(V/4.-(UV)/2.+(U*V2)/2.-V2/4.)*k1+(U/4.-(UV)/2.+(U2*V)/2.-U2/4.)*k2+(U2/8.-U/8.+V2/8.-V/8.)*k3+(-U2/8.+U/8.+V2/8.-V/8.)*k4+(U*V-V/2.-U/2.+1/4.)*k5+(V/2.-1/4.)*k6+(U/2.-1/4.)*k7+k8/4.;
        collision_operator_phase[8] = ((U2V2)/4.-(U2*V)/4.+(U*V2)/4.-(UV)/4.)*Phi+((U*V2)/2.-(UV)/2.-V/4.+V2/4.)*k1+((UV)/2.-U/4.+(U2*V)/2.-U2/4.)*k2+(U2/8.+U/8.+V2/8.-V/8.)*k3+(-U2/8.-U/8.+V2/8.-V/8.)*k4+(V/2.-U/2.+U*V-1/4.)*k5+(V/2.-1/4.)*k6+(U/2.+1/4.)*k7+k8/4.;
      }
      for(int k=0; k<np; k++)
			{
        A = U*cx[k]+V*cy[k];
        third_order = 1./(2.*cs6)*((cx[k]*cx[k]-cs2)*cy[k]*U*U*V+(cy[k]*cy[k]-cs2)*cx[k]*U*V*V);
        fourth_order = 1./(4.*cs8)*((cx[k]*cx[k]-cs2)*(cy[k]*cy[k]-cs2)*U*U*V*V);
        if(cms_hydro)
          f1[id*np+k] = collision_operator[k];
        else
        {
          feq = wf[k]*(Press+3.*A+4.5*((cx[k]*cx[k]-cs2)*U*U+(cy[k]*cy[k]-cs2)*V*V+
                         2.*cx[k]*cy[k]*U*V)+third_order+fourth_order)-
                         0.5*wf[k]/(R*cs2)*(cx[k]*Fx+cy[k]*Fy);
          f1[id*np+k] = omega1*f1[id*np+k]+omega*feq+wf[k]/(R*cs2)*(cx[k]*Fx+cy[k]*Fy);
        }
        if(cms_phase)
          g1[id*np+k] = collision_operator_phase[k];
        else
        {
          extra_term = 4.*Phi*(1.-Phi)/xi*wf[k]*(cx[k]*Nx+cy[k]*Ny);
          geq = Phi*wf[k]*(1.+3.*A+4.5*((cx[k]*cx[k]-cs2)*U*U+(cy[k]*cy[k]-cs2)*V*V+
                          2.*cx[k]*cy[k]*U*V)+third_order+fourth_order) -
                0.5*extra_term;
          g1[id*np+k] = omega_phase1*g1[id*np+k] + omega_phase*geq + extra_term;
        }
        newi = i+cx[k];
				newj = j+cy[k];
				if(i==0 || i==nx-1)
					newi = (newi+nx)%nx;
				if(j==0 || j==ny-1)
					newj = (newj+ny)%ny;
				idn = newi*ny+newj;
        f2[idn*np+k] = f1[id*np+k];
				g2[idn*np+k] = g1[id*np+k];
			}
			if(fabs(U)>1. || isnan(U))
			  check = 1;
		}
    return check;
}
///----------------------------------------------------------------------------------------------------------------------------------
void boundary_conditions()
{
  for(int i=0; i<nx; i++)
    for(int k=0; k<np; k++)
    {
			if(cy[k]>0)
			{
				id = i*ny+0;
				idn = i*ny+1;
	      f2[id*np+k] = f2[idn*np+k];
	      g2[id*np+k] = g2[idn*np+k];
			}
			if(cy[k]<0)
			{
				id = i*ny+ny-1;
				idn = i*ny+ny-2;
				f2[id*np+k] = f2[idn*np+k];
				g2[id*np+k] = g2[idn*np+k];
			}
    }
	// for(int j=0; j<ny; j++)
	// {
	// 	id = 0*ny+j;
	// 	idn = 1*ny+j;
	// 	f2[id*np+1] = f2[idn*np+1];
	// 	g2[id*np+1] = g2[idn*np+1];
	// 	f2[id*np+5] = f2[idn*np+5];
	// 	g2[id*np+5] = g2[idn*np+5];
	// 	f2[id*np+8] = f2[idn*np+8];
	// 	g2[id*np+8] = g2[idn*np+8];
	//
	// 	id = (nx-1)*ny+j;
	// 	idn = (nx-2)*ny+j;
	// 	f2[id*np+3] = f2[idn*np+3];
	// 	g2[id*np+3] = g2[idn*np+3];
	// 	f2[id*np+6] = f2[idn*np+6];
	// 	g2[id*np+6] = g2[idn*np+6];
	// 	f2[id*np+7] = f2[idn*np+7];
	// 	g2[id*np+7] = g2[idn*np+7];
	// }
}
///---------------------------------------------------------------------------------------------------------------------------------
void bb(int tempo)
{
	int minX = (int)(*min_element(X.begin(), X.end())-1),
			maxX = (int)(*max_element(X.begin(), X.end())+1),
			minY = (int)(*min_element(Y.begin(), Y.end())-1),
			maxY = (int)(*max_element(Y.begin(), Y.end())+1);
	for(int i=minX; i<maxX; i++)
		for(int j=minY; j<maxY; j++)
			if((i-nx/2)*(i-nx/2) + (j-Y_center)*(j-Y_center)>(radius-1)*(radius-1) &&
			 		(i-nx/2)*(i-nx/2) + (j-Y_center)*(j-Y_center)<=(radius)*(radius))
			{
				id = i*ny+j;
				for(int k=0; k<np; k++)
					f2[id*np+k] = f1[id*np+opp[k]]+6.*wf[k]*cy[k]*v0*rho[id];
				}
}
///----------------------------------------------------------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
	FILE *data = fopen("data.txt","wt");
	system("mkdir vtk_fluid");
	system("mkdir vtk_solid");
	int check_mach = 0, t;
	initial_state();
	define_obstacle();
	wet_arean = wet_area = 0.;
  //printf("rhoH = %lf, rhoL = %lf\n", rhoH, rhoL);
  //printf("gravity = %e, mobility = %lf, St = %lf\n", gravity, M, Scale_time);
	UY[0] = v0;
	F0yn = F0y = 0.;
	double wet_surface_w, wagner, err_w, wet_surface_vk, vk, err_vk;
  for(t=0; t<nsteps; t++)
  {
		immersed_boundary(t);
		//immersed_boundary2(t);
		F0yn = F0y;
		F0y = Force_body_y;
    check_mach = algorithm_lattice_boltzmann(t);
    boundary_conditions();
		//bb(t);

		swap(f1,f2);
		swap(g1,g2);
		if(plot_vtk==true && t%n_out==0)
		{
			write_fluid_vtk(t);
			write_cylinder_vtk(t);
			write_pressure(t);
		}

		Pressure_resultant = 0.;
		for(int k=0; k<N/2+1; k++)
			bilinear_interpolation(k);
		if(t%n_out==0)
			diagnostics(t);
		wet_surface_vk = -DY/tan(beta_angle);
		wet_surface_w = wet_surface_vk*M_PI*0.5;
		vk = rhoH*wet_surface_vk*M_PI*0.5*pow(v0,2)/tan(beta_angle);
		wagner = rhoH*wet_surface_w*M_PI*M_PI*0.25*pow(v0,2)/tan(beta_angle);
		err_vk = fabs(vk-Pressure_resultant)/fabs(vk)*100.;
		err_w = fabs(wagner-Pressure_resultant)/fabs(wagner)*100.;
		fprintf(data,"%lf %e  %e %e\n", t*Scale_time, Pressure_resultant*Scale_PressRes, wagner*Scale_PressRes, vk*Scale_PressRes);

				update_position(t);
		if(check_mach==1) // check the Mach number...if too high, it exits!
      goto labelA;
  }
  labelA:
	  fclose(data);
  return 0;
}
///----------------------------------------------------------------------------------------------------------------------------------
///----------------------------------------------------------------------------------------------------------------------------------
