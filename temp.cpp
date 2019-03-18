#include <math.h>
#include <stdio.h>
#include <iostream>
#include <armadillo>
using namespace std;
using namespace arma;
static float sqrarg;
static float cubearg;
static float pow4arg;
#define pi 3.14159
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)
#define CUBE(a) ((cubearg=(a)) == 0.0 ? 0.0 : cubearg*cubearg*cubearg)
#define POW4(a) ((pow4arg=(a)) == 0.0 ? 0.0 :  pow4arg*pow4arg*pow4arg*pow4arg)

//Global Variables
float	omhh,		/* Omega_matter*h^2 */
	obhh,		/* Omega_baryon*h^2 */
	theta_cmb,	/* Tcmb in units of 2.7 K */
	z_equality,	/* Redshift of matter-radiation equality, really 1+z */
	k_equality,	/* Scale of equality, in Mpc^-1 */
	z_drag,		/* Redshift of drag epoch */
	R_drag,		/* Photon-baryon ratio at drag epoch */
	R_equality,	/* Photon-baryon ratio at equality epoch */
	sound_horizon,	/* Sound horizon at drag epoch, in Mpc */
	k_silk,		/* Silk damping scale, in Mpc^-1 */
	alpha_c,	/* CDM suppression */
	beta_c,		/* CDM log shift */
	alpha_b,	/* Baryon suppression */
	beta_b,		/* Baryon envelope shift */
	beta_node,	/* Sound horizon shift */
	k_peak,		/* Fit to wavenumber of first peak, in Mpc^-1 */
	sound_horizon_fit,	/* Fit to sound horizon, in Mpc */
	alpha_gamma;	/* Gamma suppression in approximate TF */


void TFset_parameters(float omega0hh, float f_baryon, float Tcmb)
/* Set all the scalars quantities for Eisenstein & Hu 1997 fitting formula */
/* Input: omega0hh -- The density of CDM and baryons, in units of critical dens, multiplied by the square of the Hubble constant, in units of 100 km/s/Mpc */
/* 	  f_baryon -- The fraction of baryons to CDM */
/*        Tcmb -- The temperature of the CMB in Kelvin.  Tcmb<=0 forces use of the COBE value of  2.728 K. */
/* Output: Nothing, but set many global variables used in TFfit_onek(). 
/* Note: Units are always Mpc, never h^-1 Mpc. */
{
	float z_drag_b1, z_drag_b2;
    	float alpha_c_a1, alpha_c_a2, beta_c_b1, beta_c_b2, alpha_b_G, y;

    	if (f_baryon<=0.0 || omega0hh<=0.0) 
	{
		cout<<"TFset_parameters(): Illegal input.\n";
		return;
	}
    	omhh = omega0hh;
    	obhh = omhh*f_baryon;
    	if (Tcmb<=0.0) 
		Tcmb=2.728;	/* COBE FIRAS */
    	theta_cmb = Tcmb/2.7;

    	z_equality = 2.50e4*omhh/POW4(theta_cmb);  /* Really 1+z */
    	k_equality = 0.0746*omhh/SQR(theta_cmb);

    	z_drag_b1 = 0.313*pow(omhh,-0.419)*(1+0.607*pow(omhh,0.674));
    	z_drag_b2 = 0.238*pow(omhh,0.223);
    	z_drag = 1291*pow(omhh,0.251)/(1+0.659*pow(omhh,0.828))*(1+z_drag_b1*pow(obhh,z_drag_b2));
    
    	R_drag = 31.5*obhh/POW4(theta_cmb)*(1000/(1+z_drag));
    	R_equality = 31.5*obhh/POW4(theta_cmb)*(1000/z_equality);

    	sound_horizon = 2./3./k_equality*sqrt(6./R_equality)*log((sqrt(1+R_drag)+sqrt(R_drag+R_equality))/(1+sqrt(R_equality)));

    	k_silk = 1.6*pow(obhh,0.52)*pow(omhh,0.73)*(1+pow(10.4*omhh,-0.95));

    	alpha_c_a1 = pow(46.9*omhh,0.670)*(1+pow(32.1*omhh,-0.532));
    	alpha_c_a2 = pow(12.0*omhh,0.424)*(1+pow(45.0*omhh,-0.582));
    	alpha_c = pow(alpha_c_a1,-f_baryon)*pow(alpha_c_a2,-CUBE(f_baryon));
    
    	beta_c_b1 = 0.944/(1+pow(458*omhh,-0.708));
    	beta_c_b2 = pow(0.395*omhh, -0.0266);
    	beta_c = 1.0/(1+beta_c_b1*(pow(1-f_baryon, beta_c_b2)-1));

    	y = z_equality/(1+z_drag);
    	alpha_b_G = y*(-6.*sqrt(1+y)+(2.+3.*y)*log((sqrt(1+y)+1)/(sqrt(1+y)-1)));
    	alpha_b = 2.07*k_equality*sound_horizon*pow(1+R_drag,-0.75)*alpha_b_G;

    	beta_node = 8.41*pow(omhh, 0.435);
    	beta_b = 0.5+f_baryon+(3.-2.*f_baryon)*sqrt(pow(17.2*omhh,2.0)+1);

    	k_peak = 2.5*3.14159*(1+0.217*omhh)/sound_horizon;
    	sound_horizon_fit = 44.5*log(9.83/omhh)/sqrt(1+10.0*pow(obhh,0.75));

    	alpha_gamma = 1-0.328*log(431.0*omhh)*f_baryon + 0.38*log(22.3*omhh)*SQR(f_baryon);

    	return;
}


float TFfit_onek(float k, float *tf_baryon, float *tf_cdm)
/* Input: 	k -- Wavenumber at which to calculate transfer function, in Mpc^-1.
	  	*tf_baryon, *tf_cdm -- Input value not used; replaced on output if the input was not NULL. */
/* Output: 	Returns the value of the full transfer function fitting formula.
		*tf_baryon -- The baryonic contribution to the full fit.
	  	*tf_cdm -- The CDM contribution to the full fit. */
/* Notes: 	Units are Mpc, not h^-1 Mpc. */
{
	float T_c_ln_beta, T_c_ln_nobeta, T_c_C_alpha, T_c_C_noalpha;
    	float q, xx, xx_tilde, q_eff;
    	float T_c_f, T_c, s_tilde, T_b_T0, T_b, f_baryon, T_full;
    	float T_0_L0, T_0_C0, T_0, gamma_eff; 
    	float T_nowiggles_L0, T_nowiggles_C0, T_nowiggles;

    	k = fabs(k);	/* Just define negative k as positive */
    	if (k==0.0) 
	{
		if (tf_baryon!=NULL) *tf_baryon = 1.0;
		if (tf_cdm!=NULL) *tf_cdm = 1.0;
		return 1.0;
	}

    	q = k/13.41/k_equality;
    	xx = k*sound_horizon;

    	T_c_ln_beta = log(2.718282+1.8*beta_c*q);
    	T_c_ln_nobeta = log(2.718282+1.8*q);
    	T_c_C_alpha = 14.2/alpha_c + 386.0/(1+69.9*pow(q,1.08));
    	T_c_C_noalpha = 14.2 + 386.0/(1+69.9*pow(q,1.08));

    	T_c_f = 1.0/(1.0+POW4(xx/5.4));
    	T_c = T_c_f*T_c_ln_beta/(T_c_ln_beta+T_c_C_noalpha*SQR(q)) + (1-T_c_f)*T_c_ln_beta/(T_c_ln_beta+T_c_C_alpha*SQR(q));
    
	s_tilde = sound_horizon*pow(1+CUBE(beta_node/xx),-1./3.);
    	xx_tilde = k*s_tilde;

    	T_b_T0 = T_c_ln_nobeta/(T_c_ln_nobeta+T_c_C_noalpha*SQR(q));
    	T_b = sin(xx_tilde)/(xx_tilde)*(T_b_T0/(1+SQR(xx/5.2)) + alpha_b/(1+CUBE(beta_b/xx))*exp(-pow(k/k_silk,1.4)));
    
    	f_baryon = obhh/omhh;
    	T_full = f_baryon*T_b + (1-f_baryon)*T_c;

    	if (tf_baryon!=NULL) 
		*tf_baryon = T_b;
    	if (tf_cdm!=NULL) 
		*tf_cdm = T_c;
   	return T_full;
}

float wFs(float x)
{
	float y = (3*(sin(x) - x*cos(x)))/(pow(x, 3));
	return pow(y, 2);
}

float calcSigmaIntegrand(float y)
{
	float h = 0.7;	
	float ey = exp(y);
	float z = TFfit_onek(ey*h, NULL, NULL);
	float result = z*z;
	result *= wFs(ey * 8/h);
	result *= exp(4*y);
	return result;

}

float simpson(int type, float upper, float lower)
{
 	double h = (upper-lower)/(double)pow(10, 4);
	double x = lower+h;
	double sum = 0.0;
	while(x+2*h <= upper)
	{
		if (type == 1)
		{
			sum += calcSigmaIntegrand(x);
			sum += 4*calcSigmaIntegrand(x+h);
			sum += calcSigmaIntegrand(x+2*h);
		}
		x+=2*h;
		//cout<<sum<<endl;	
	}
	sum *= (h/3);
	return sum;
}

int main()
{
	float H0;
	float omega[4];
	float z;
	float w[4] = {0.3333, 0.0, 0.0, -1.0};
	int N;
	int def;
	float kmin, kmax; 
	cout<<"Enter 1. for Default values (H = 69.6, DM = 0.3, Lambda = 0.7, z = 0, N = 1000)\n";
	cout<<"Enter 2. to set values \n";
	cout<<"Enter choice : ";
	cin>>def;
	cout<<endl;

	if (def == 1)
	{
		H0 = 70;
		omega[0] = 0.0;
		omega[1] = 0.02;
		omega[2] = 0.28;
		omega[3] = 0.7; 
		z = 2.0;
		N = 100;
		kmin = 0.001;
		kmax = 10;
	}
	else if (def == 2)
	{
		cout<<"Enter H0 : ";
		cin>>H0;
		cout<<"Enter omegaR : ";
		cin>>omega[0];
		cout<<"Enter omegaB : ";
		cin>>omega[1];
		cout<<"Enter omegaDM : ";
		cin>>omega[2];
		cout<<"Enter omegaLambda : ";
		cin>>omega[3];
		cout<<"Enter z : ";
		cin>>z;
		cout<<"Enter N : ";
		cin>>N;
		cout<<"Enter Kmin : ";
		cin>>kmin;
		cout<<"Enter Kmax : ";
		cin>>kmax;
	}
	
	float h = H0/100;

	TFset_parameters(h*h*(omega[1]+omega[2]), omega[1]/omega[2], -1);

	mat values = mat(N, 2);
	float k;
	float alpha = ((float)log10(kmax/kmin))/(N-1);
	float A = kmin/pow(10, alpha); 
	cout<<"Alpha : "<<alpha<<" A : "<<A<<endl;

	double sigma2 = (1/(2*pi*pi))*simpson(1, 5, -10);
	double sigma = pow(sigma2, 0.5);
	cout<<"Sigma2 value : "<<sigma2<<endl;
	cout<<"Sigma value : "<<sigma<<endl;

	double rscaling = 0.82/sigma;
	cout<<"Scaling = "<<rscaling<<endl;	
	for (int i=1; i<=N; i++)
	{
		k = A*pow(10, i*alpha);
		values(i-1, 0) = k;
		float z = TFfit_onek(k*h, NULL, NULL);
		values(i-1, 1) = rscaling*rscaling*k*z*z;
	}
	string str = "values.txt";
	values.save(str, raw_ascii);
	
}

