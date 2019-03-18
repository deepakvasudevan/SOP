#include <math.h>
#include <stdio.h>
#include <iostream>
#include <armadillo>
using namespace std;
using namespace arma;

double omegaM0 = 0.3;
double H0 = 70;
double w_0, w_a;

double W_7CPL(double z)
{
	double ans;
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
	double h = (upper-lower)/100;
	//cout<<h<<endl;
	double x = lower + h;
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

double calcIntegrand(double x)
{
	double ans = 1.0;
	ans /= pow(x, 3);
	ans /= pow(calcH(1.0/x - 1.0), 3);
	return ans;

}
double calcIntegral(double a)
{
	double upper, lower;	
	lower = 0.0;
	upper = a;
	if (upper == lower)
		return 0.0;
	double h = (upper-lower)/100;
	//cout<<h<<endl;
	double x = lower + h;
	double sum = 0.0;
	while(x+2*h <= upper)
	{
		sum += calcIntegrand(x);
		sum += 4*calcIntegrand(x+h);
		sum += calcIntegrand(x+2*h);
		x+=2*h;	
	}
	sum *= (h/3);
	return sum;

}

double getAfromZ(double z)
{
	double ans;
	ans = 1.0/(1.0 + z);
	return ans;
}

double calcdelta(double z)
{	
	double a = getAfromZ(z);
	double ans = 5.0/2.0;
	ans *= omegaM0;
	ans *= calcH(z);
	ans *= calcIntegral(a);
	return ans;
}

int main()
{
	int N = 1000;	
	mat deltaG = mat(N+1, 2);
	double za = 0.0;
	double zb = 10.0;
	double temp = za;
	double h = (zb - za)/N;
	int inp;
	cout<<"Enter 1. for Normal 2. for LCDM\n";
	cin>>inp;
	if (inp == 1)
	{
		cout<<"Enter w_0 : ";
		cin>>w_0;
		cout<<"Enter w_a : ";
		cin>>w_a;
	}
	else
	{
		w_0 = -1.0;
		w_a = 0.0;
	}
	for (int i=0; i<=N; i++)
	{
		deltaG(i, 0) = temp;
		if (i%10 == 0)
			cout<<(float)i*100.0/(float)N<<"% done\n";
		deltaG(i, 1) = calcdelta(temp);
		temp = temp + h;
		
	}
	string str;
	if (inp == 1)
		str = "newRkData.txt";
	else
		str = "newRkDataLCDM.txt";
	deltaG.save(str, raw_ascii);
	
}
