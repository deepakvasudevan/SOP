`#include <math.h>
#include <stdio.h>
#include <iostream>
#include <armadillo>
using namespace std;
using namespace arma;

double h;
double H0, omegaM0;
int inp;
int points = 5000;
double k;
double w_0, w_a;

double y_new;
double y_old;
double x_new;
double x_old = log(0.001);
double z_new;
double z_old; 

double k0, k1, k2, k3, l0, l1, l2, l3;

double calcZ(double x)
{
	double ans = exp(-1.0*x) - 1.0;
	return ans;
}	


double W_7CPL(double z)
{
	double ans;
	//ans = w_0 + w_a * pow((z/(1+z)), 1);
	ans = w_0 + w_a * pow((z/(1+z)), 7);
	return ans;
}
double calcHRadIntegrand(double z)
{
	double ans;
	ans = 1.0 + W_7CPL(z);
	ans /= (1.0 + z);
	return ans;
}

double calcHRad(double endz)
{
	double upper, lower;	
	lower = 0.0;
	upper = endz;
	if (upper == lower)
		return 0.0;
	double h = (upper-lower)/(5.0*pow(10, 3));
	double x = lower;
	double sum = 0.0;
	while(x+2*h <= upper)
	{
		sum += calcHRadIntegrand(x);
		sum += 4*calcHRadIntegrand(x+h);
		sum += calcHRadIntegrand(x+2*h);
		x+=2*h;	
	}
	sum *= (h/3);
	return sum;
}


double calcH(double z)
{
	double ans;
	ans = omegaM0 * pow((1+z), 3);
	ans += (1-omegaM0) * exp(3*calcHRad(z));
	ans = pow(ans, 0.5);
	ans *= H0;
	return ans;

}

double calcHLCDM(double z)
{
	double ans;
	ans = omegaM0 * pow((1+z), 3);
	ans += (1-omegaM0);
	ans = pow(ans, 0.5);
	ans *= H0;
	return ans;
}

double omegaM(double z)
{
	double ans;
	ans = pow((1+z), 3);
	ans *= omegaM0;
	if (inp == 2)
		ans *= pow((H0/calcHLCDM(z)), 2);
	else
		ans *= pow((H0/calcH(z)), 2);
	return ans;
}

double omegaQ(double z)
{
	double ans;
	ans = (1-omegaM0) * exp(3*calcHRad(z));
	return ans;
}


double f(double x, double y, double z)
{
	return z;
}

double g(double x, double y, double z)
{
	double ans;
	double redshift = calcZ(x);
	double temp = omegaM(redshift);
	ans = 3 * temp * y - 1.0 + 3 * W_7CPL(redshift) * (1-temp);
	return ans/2.0;
}

/*double g(double x, double y, double z)
{
	return 2.0;
}*/

int main()
{
	double a = log(0.001);
	double b = 0.0;
	h = (b-a)/points;
	
	cout<<"Input 1. for normal 2. for LCDM : ";
	cin>>inp;
	if (inp == 1)
	{
		w_0 = -1.1;
		w_a = 0.1;
	}
	else if (inp == 2)
	{
		w_0 = -1.0;
		w_a = 0.0;
	}
	else
	{
		w_0 = -0.9;
		w_a = -0.1;
	}
	cout<<"Input scaling : ";
	cin>>k;
	z_old = k * 0.001;
	y_old = k * 0.001;
	H0 = 69.0;
	omegaM0 = 0.313;
	mat values = mat(points+1, 2);
	for(int i=0;i<=points;i++)
	{ 
		k0 = h*f(x_old, y_old, z_old);
      		l0 = h*g(x_old, y_old, z_old);
		
	      	k1 = h*f(x_old + 0.5*h, y_old + 0.5*k0, z_old + 0.5*l0);
	     	l1 = h*g(x_old + 0.5*h, y_old + 0.5*k0, z_old + 0.5*l0);

	      	k2 = h*f(x_old + 0.5*h, y_old + 0.5*k1, z_old + 0.5*l1);
	      	l2 = h*g(x_old + 0.5*h, y_old + 0.5*k1, z_old + 0.5*l1);

	      	k3 = h*f(x_old + h, y_old + k2, z_old + l2);
	      	l3 = h*g(x_old + h, y_old + k2, z_old + l2);

		x_new=x_old+h;	      	
		y_new = y_old+ (k0+2*k1+2*k2+k3)/6;
	      	z_new = z_old+ (l0+2*l1+2*l2+l3)/6;

   		values(i, 0) = calcZ(x_old);
		values(i, 1) = y_old;

		x_old=x_new;
      		y_old=y_new;
      		z_old=z_new;
		if (i%10 == 0)
			cout<<((double)i/(double)points)*100<<"% done\n";
	}
	string str;
	if (inp == 1)
		str = "upper.txt";
	else if (inp == 2)
		str = "LCDM.txt";
	else 
		str = "lower.txt";
	values.save(str, raw_ascii);
}
