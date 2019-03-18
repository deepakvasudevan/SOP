#include<iostream>
#include<math.h>
//#include<armadillo>
#include<stdio.h>
using namespace std;
//using namespace arma;
static float sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

/*Global variables*/
/* The following are set in TFmdm_set_cosm() */
float   alpha_gamma,	/* sqrt(alpha_nu) */
	alpha_nu,	/* The small-scale suppression */
	beta_c,		/* The correction to the log in the small-scale */
	num_degen_hdm,	/* Number of degenerate massive neutrino species */
	f_baryon,	/* Baryon fraction */
	f_bnu,		/* Baryon + Massive Neutrino fraction */
	f_cb,		/* Baryon + CDM fraction */
	f_cdm,		/* CDM fraction */
	f_hdm,		/* Massive Neutrino fraction */
	growth_k0,	/* D_1(z) -- the growth function as k->0 */
	growth_to_z0,	/* D_1(z)/D_1(0) -- the growth relative to z=0 */
	hhubble,	/* Need to pass Hubble constant to TFmdm_onek_hmpc() */
	k_equality,	/* The comoving wave number of the horizon at equality*/
	obhh,		/* Omega_baryon * hubble^2 */
	omega_curv,	/* = 1 - omega_matter - omega_lambda */
	omega_lambda_z, /* Omega_lambda at the given redshift */
	omega_matter_z,	/* Omega_matter at the given redshift */
	omhh,		/* Omega_matter * hubble^2 */
	onhh,		/* Omega_hdm * hubble^2 */
	p_c,		/* The correction to the exponent before drag epoch */
	p_cb,		/* The correction to the exponent after drag epoch */
	sound_horizon_fit,  /* The sound horizon at the drag epoch */
	theta_cmb,	/* The temperature of the CMB, in units of 2.7 K */
	y_drag,		/* Ratio of z_equality to z_drag */
	z_drag,		/* Redshift of the drag epoch */
	z_equality;	/* Redshift of matter-radiation equality */

/* The following are set in TFmdm_onek_mpc() */
float	gamma_eff,	/* Effective \Gamma */
	growth_cb,	/* Growth factor for CDM+Baryon perturbations */
	growth_cbnu,	/* Growth factor for CDM+Baryon+Neutrino pert. */
	max_fs_correction,  /* Correction near maximal free streaming */
	qq,		/* Wavenumber rescaled by \Gamma */
	qq_eff,		/* Wavenumber rescaled by effective Gamma */
	qq_nu,		/* Wavenumber compared to maximal free streaming */
	tf_master,	/* Master TF */
	tf_sup,		/* Suppressed TF */
	y_freestream; 	/* The epoch of free-streaming for a given scale */


/* Finally, TFmdm_onek_mpc() and TFmdm_onek_hmpc() give their answers as */
float   tf_cb,		/* The transfer function for density-weighted 
			CDM + Baryon perturbations. */
	tf_cbnu;	/* The transfer function for density-weighted
			CDM + Baryon + Massive Neutrino perturbations. */

int TFmdm_set_cosm(float omega_matter, float omega_baryon, float omega_hdm, int degen_hdm, float omega_lambda, float hubble, float redshift)
{
	/* This routine takes cosmological parameters and a redshift and sets up
	all the internal scalar quantities needed to compute the transfer function. */
	/* INPUT: omega_matter -- Density of CDM, baryons, and massive neutrinos,
				in units of the critical density. */
	/* 	  omega_baryon -- Density of baryons, in units of critical. */
	/* 	  omega_hdm    -- Density of massive neutrinos, in units of critical */
	/* 	  degen_hdm    -- (Int) Number of degenerate massive neutrino species */
	/*        omega_lambda -- Cosmological constant */
	/* 	  hubble       -- Hubble constant, in units of 100 km/s/Mpc */
	/*        redshift     -- The redshift at which to evaluate */
	/* OUTPUT: Returns 0 if all is well, 1 if a warning was issued.  Otherwise,
		sets many global variables for use in TFmdm_onek_mpc() */
	float z_drag_b1, z_drag_b2, omega_denom;
	int qwarn;
	qwarn = 0;

	theta_cmb = 2.728/2.7;	/* Assuming T_cmb = 2.728 K */

    /* Look for strange input */
	if (omega_baryon<0.0) 
	{
		cout<<"TFmdm_set_cosm(): Negative omega_baryon set to trace amount.\n";
		qwarn = 1;
    	}
    	if (omega_hdm<0.0) 
	{
		cout<<"TFmdm_set_cosm(): Negative omega_hdm set to trace amount.\n";
		qwarn = 1;
    	}	
    	
	if (hubble<=0.0) 
	{
		cout<<"TFmdm_set_cosm(): Negative Hubble constant illegal.\n";
		return 1;
	} 
	else if (hubble>2.0) 
	{
		cout<<"TFmdm_set_cosm(): Hubble constant should be in units of 100 km/s/Mpc.\n";
		qwarn = 1;
	}
    
	if (redshift<=-1.0) 
	{
		cout<<"TFmdm_set_cosm(): Redshift < -1 is illegal.\n";
		return 1;
	}
    	else if (redshift>99.0) 
	{
		cout<<"TFmdm_set_cosm(): Large redshift entered.  TF may be inaccurate.\n";
		qwarn = 1;
    	}	
    
	if (degen_hdm<1) 
		degen_hdm=1;
	num_degen_hdm = (float) degen_hdm;	
	/* Have to save this for TFmdm_onek_mpc() */
    	/* This routine would crash if baryons or neutrinos were zero, so don't allow that */
	if (omega_baryon<=0) omega_baryon=1e-5;
    	if (omega_hdm<=0) omega_hdm=1e-5;

    	omega_curv = 1.0-omega_matter-omega_lambda;
    	omhh = omega_matter*SQR(hubble);
    	obhh = omega_baryon*SQR(hubble);
    	onhh = omega_hdm*SQR(hubble);
    	f_baryon = omega_baryon/omega_matter;
    	f_hdm = omega_hdm/omega_matter;
    	f_cdm = 1.0-f_baryon-f_hdm;
    	f_cb = f_cdm+f_baryon;
    	f_bnu = f_baryon+f_hdm;

    	/* Compute the equality scale. */
    	z_equality = 25000.0*omhh/SQR(SQR(theta_cmb));		/* Actually 1+z_eq */
    	k_equality = 0.0746*omhh/SQR(theta_cmb);

    	/* Compute the drag epoch and sound horizon */
    	z_drag_b1 = 0.313*pow(omhh,-0.419)*(1+0.607*pow(omhh,0.674));
    	z_drag_b2 = 0.238*pow(omhh,0.223);
    	z_drag = 1291*pow(omhh,0.251)/(1.0+0.659*pow(omhh,0.828))*(1.0+z_drag_b1*pow(obhh,z_drag_b2));
	y_drag = z_equality/(1.0+z_drag);

	sound_horizon_fit = 44.5*log(9.83/omhh)/sqrt(1.0+10.0*pow(obhh,0.75));

    	/* Set up for the free-streaming & infall growth function */
    	p_c = 0.25*(5.0-sqrt(1+24.0*f_cdm));
    	p_cb = 0.25*(5.0-sqrt(1+24.0*f_cb));

    	omega_denom = omega_lambda+SQR(1.0+redshift)*(omega_curv+omega_matter*(1.0+redshift));
    	omega_lambda_z = omega_lambda/omega_denom;
    	omega_matter_z = omega_matter*SQR(1.0+redshift)*(1.0+redshift)/omega_denom;
    	growth_k0 = z_equality/(1.0+redshift)*2.5*omega_matter_z/(pow(omega_matter_z,4.0/7.0)-omega_lambda_z+(1.0+omega_matter_z/2.0)*(1.0+omega_lambda_z/70.0));
    	growth_to_z0 = z_equality*2.5*omega_matter/(pow(omega_matter,4.0/7.0)-omega_lambda + (1.0+omega_matter/2.0)*(1.0+omega_lambda/70.0));
    	growth_to_z0 = growth_k0/growth_to_z0;	
    
    	/* Compute small-scale suppression */
    	alpha_nu = f_cdm/f_cb*(5.0-2.*(p_c+p_cb))/(5.-4.*p_cb)*pow(1+y_drag,p_cb-p_c)*(1+f_bnu*(-0.553+0.126*f_bnu*f_bnu))/(1-0.193*sqrt(f_hdm*num_degen_hdm)+0.169*f_hdm*pow(num_degen_hdm,0.2))*(1+(p_c-p_cb)/2*(1+1/(3.-4.*p_c)/(7.-4.*p_cb))/(1+y_drag));
    	alpha_gamma = sqrt(alpha_nu);
    	beta_c = 1/(1-0.949*f_bnu);
    	/* Done setting scalar variables */
    	hhubble = hubble;	/* Need to pass Hubble constant to TFmdm_onek_hmpc() */
    	return qwarn;

}


float TFmdm_onek_mpc(float kk)
{
	/* Given a wavenumber in Mpc^-1, return the transfer function for the
	cosmology held in the global variables. */
	/* Input: kk -- Wavenumber in Mpc^-1 */
	/* Output: The following are set as global variables:
			growth_cb -- the transfer function for density-weighted CDM + Baryon perturbations. 
		 	growth_cbnu -- the transfer function for density-weighted CDM + Baryon + Massive Neutrino perturbations. */
	/* The function returns growth_cb */
	float tf_sup_L, tf_sup_C;
    	float temp1, temp2;

    	qq = kk/omhh*SQR(theta_cmb);

    	/* Compute the scale-dependent growth functions */
    	y_freestream = 17.2*f_hdm*(1+0.488*pow(f_hdm,-7.0/6.0))*SQR(num_degen_hdm*qq/f_hdm);
	temp1 = pow(growth_k0, 1.0-p_cb);
    	temp2 = pow(growth_k0/(1+y_freestream),0.7);
    	growth_cb = pow(1.0+temp2, p_cb/0.7)*temp1;
    	growth_cbnu = pow(pow(f_cb,0.7/p_cb)+temp2, p_cb/0.7)*temp1;

    	/* Compute the master function */
    	gamma_eff =omhh*(alpha_gamma+(1-alpha_gamma)/(1+SQR(SQR(kk*sound_horizon_fit*0.43))));
	qq_eff = qq*omhh/gamma_eff;

    	tf_sup_L = log(2.71828+1.84*beta_c*alpha_gamma*qq_eff);
   	tf_sup_C = 14.4+325/(1+60.5*pow(qq_eff,1.11));
    	tf_sup = tf_sup_L/(tf_sup_L+tf_sup_C*SQR(qq_eff));

    	qq_nu = 3.92*qq*sqrt(num_degen_hdm/f_hdm);
    	max_fs_correction = 1+1.2*pow(f_hdm,0.64)*pow(num_degen_hdm,0.3+0.6*f_hdm)/(pow(qq_nu,-1.6)+pow(qq_nu,0.8));
   	tf_master = tf_sup*max_fs_correction;

    	/* Now compute the CDM+HDM+baryon transfer functions */
    	tf_cb = tf_master*growth_cb/growth_k0;
    	tf_cbnu = tf_master*growth_cbnu/growth_k0;
    	return tf_cb;		
}



float TFmdm_onek_hmpc(float kk)
{
	/* Given a wavenumber in h Mpc^-1, return the transfer function for the cosmology held in the global variables. */
	/* Input: kk -- Wavenumber in h Mpc^-1 */
	/* Output: The following are set as global variables:
			growth_cb -- the transfer function for density-weighted CDM + Baryon perturbations. 
		 	growth_cbnu -- the transfer function for density-weighted CDM + Baryon + Massive Neutrino perturbations. */
	/* The function returns growth_cb */
	return TFmdm_onek_mpc(kk*hhubble);
}


double convertAgetoYears(double x)
{
	double sec = x*3.086*pow(10, 19);
	sec = sec/(60.0*60.0*24.0*365.25);
	return sec;
}
double convertSectoMPc(double x)
{
	return x*300000;
}
double calcAgeIntegrand(double x, double H0, double omega[], double w[])
{
	double sum = 0.0, temp = 0.0;
	for (int i=0; i<4; i++)
		temp += omega[i]*pow(x, -3*(1+w[i]));
	//cout<<temp<<endl;
	sum = x*H0*pow(temp, 0.5);
	return (1.0/sum);
		
}
double calcRadIntegrand(double x, double H0, double omega[], double w[])
{
	double sum = 0.0, temp = 0.0;
	for (int i=0; i<4; i++)
		temp += omega[i]*pow(x, -3*(1+w[i]));
	//cout<<temp<<endl;
	sum = x*H0*pow(temp, 0.5);
	sum = (1.0/(sum*x));
	return sum;
}

double calcRad(double H0, double omega[], double w[], double endz)
{
	double upper, lower;	
	lower = 0.0;
	upper = 1.0/(1+endz);

	double h = (upper-lower)/pow(10, 5);
	double x = lower+h;
	double sum = 0.0;
	while(x+2*h <= upper)
	{
		sum += calcRadIntegrand(x, H0, omega, w);
		sum += 4*calcRadIntegrand(x+h, H0, omega, w);
		sum += calcRadIntegrand(x+2*h, H0, omega, w);
		x+=2*h;	
	}
	sum *= (h/3);
	return sum;
}
double calcAge(double H0, double omega[], double w[], double z)
{
	double upper, lower;
	lower = 0.0;
	upper = 1.0/(1+z);

	double h = (upper-lower)/pow(10, 5);
	double x = lower+h;
	double sum = 0.0;
	while(x+2*h <= upper)
	{
		sum += calcAgeIntegrand(x, H0, omega, w);
		sum += 4*calcAgeIntegrand(x+h, H0, omega, w);
		sum += calcAgeIntegrand(x+2*h, H0, omega, w);
		x+=2*h;	
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
		H0 = 69.6;
		omega[0] = 0.0;
		omega[1] = 0.0;
		omega[2] = 0.3;
		omega[3] = 0.7; 
		z = 0.0;
		N = 1000;
		kmin = 0.001;
		kmax = 10;
	}
	else
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
	int error = TFmdm_set_cosm(omega[2], omega[1], 0.0, 1, omega[3], H0/100.0, z);
	if (error == 1)
	{
		cout<<"Invalid Parameters.\n";
		return 0;
	}
	mat values = mat(N, 2);
	float k;
	float alpha = ((float)log10(kmax/kmin))/(N-1);
	float A = kmin/pow(10, alpha); 
	cout<<"Alpha : "<<alpha<<" A : "<<A<<endl;
	for (int i=1; i<=N; i++)
	{
		k = A*pow(10, i*alpha);
		values(i-1, 0) = k;
		values(i-1, 1) = TFmdm_onek_hmpc(k);
	}
	string str = "values.txt";
	values.save(str, raw_ascii);
	
}
