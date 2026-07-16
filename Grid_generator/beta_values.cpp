#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
using namespace std;

void cl(vector<double> *x, double dxi,double Lx, double beta, int Nx, double x0)
{
	int i;
	double alpha = 0.0;
	double temp1 = (1+beta)/(beta-1);
	for(i=0;i<Nx;i++)
	{
		double xi_i = i*dxi;
		double temp2 = (xi_i-alpha)/(1-alpha);
		double exp = pow(temp1,temp2);

		(*x)[i] = x0 + Lx*((beta+2*alpha)*exp - beta + 2*alpha)/((1+2*alpha)*(1+exp));
	}
}

void c0(vector<double> *x, double dxi,double Lx, double beta, int Nx, double x0)
{
	int i;
	double temp1 = (1+beta)/(beta-1);
	for(i=0;i<Nx;i++)
	{
		double xi_i = i*dxi;
		double temp2 = 1-xi_i;
		double exp = pow(temp1,temp2);

		(*x)[i] = x0 + Lx*((beta +1)-(beta-1)*exp)/(1+exp);
	}
}


void cs(vector<double> *x, double dxi,double Lx, double  beta, int Nx, double x0)
{
	int i;
	double alpha = 0.5;
	double temp1 = (1+beta)/(beta-1);
	for(i=0;i<Nx;i++)
	{
		double xi_i = i*dxi;
		double temp2 = (xi_i-alpha)/(1-alpha);
		double exp = pow(temp1,temp2);

		(*x)[i] = x0 + Lx*((beta+2*alpha)*exp - beta + 2*alpha)/((1+2*alpha)*(exp+1));
	}

}

void writeVectorToCSV(const vector<double>& data, const string& filename)
{
    ofstream file(filename);

    if (!file.is_open()) {
        cerr << "Error opening file!" << endl;
        return;
    }

    for (size_t i = 0; i < data.size(); i++) {
        file << data[i]<<"\n";

    }

    file << endl;

    file.close();
}

double error(vector <double> x2, vector <double> x3)
{
  double dx2 = abs(x2.back() - x2[x2.size()-2]);
  double dx23 = abs(x3[0] - x2.back());

  double e23 = abs((dx23 - dx2)/dx2);
  return e23;
}

int main()
{
vector<double> betasx,Lx;
vector <int> Nx;

ifstream filem("grid_2.txt");
double bx,nx,lx;
while(filem>>bx>>nx>>lx)
{
	betasx.push_back(bx);
	Nx.push_back(nx);
	Lx.push_back(lx);
}
filem.close(); 

vector <double> x1(Nx[0]),x2(Nx[1]),x3(Nx[2]),x4(Nx[3]),x5(Nx[4]);

double dxi[5];
for(int i =0;i<5;i++)
{
  dxi[i]=(double)1/(Nx[i]-1);
}

cl(&x1,dxi[0],Lx[0],betasx[0],Nx[0], 0);
cs(&x2,dxi[1],Lx[1],betasx[1],Nx[1], Lx[0]);
cs(&x3,dxi[2],Lx[2],betasx[2],Nx[2], Lx[0]+Lx[1]);
cs(&x4,dxi[3],Lx[3],betasx[3],Nx[3], Lx[0]+Lx[1]+Lx[2]);
c0(&x5,dxi[4],Lx[4],betasx[4],Nx[4], Lx[0]+Lx[1]+Lx[2]+Lx[3]);

double beta1 = 1.001;
double beta2 = 150.0;

double target = 0.003;
double step = 0.0001;

double dx0 = abs(x1.back()-x1[x1.size()-2]);
while(true)
{

    cl(&x1, dxi[0], Lx[0], beta1, Nx[0], 0);
    cs(&x2, dxi[1], Lx[1], beta2, Nx[1], Lx[0]);

    
    vector<double> xtemp = x2;
    xtemp.erase(xtemp.begin());

    
    double dx1  = abs(x1.back() - x1[x1.size()-2]);
    double dx12 = abs(xtemp[0] - x1.back());

   
    double e12 = abs((dx12 - dx1)/dx1);

   /* cout << beta1 << " "
         << beta2 << " "
         << e12 << endl;
*/
    if(e12 < target)
    {
        cout << "\nconverged\n";

        cout << "beta1 = " << beta1 << endl;
        cout << "beta2 = " << beta2 << endl;
        cout << "e12 = " << e12 << endl;
        
        break;
    }


    beta1 +=beta1*step;
    beta2 -=beta2*step;

    

    if(beta2 <= 1.01 || beta1>=2)
    {
        cout << "no convergence\n";
        break;
    }
}

double beta3 = 1.001;
double beta2_2 = 150.0;

double target3 = 0.002;
double step3 = 0.0001;

double dx20 = abs(x2.back()-x2[x2.size()-2]);
while(true)
{

    cs(&x2, dxi[1], Lx[1], beta2_2, Nx[1], Lx[0]);
    cs(&x3, dxi[2], Lx[2], beta3, Nx[2], Lx[0]+Lx[1]);

    
    vector<double> xtemp = x3;
    xtemp.erase(xtemp.begin());

    
    double dx2  = abs(x2.back() - x2[x2.size()-2]);
    double dx23 = abs(xtemp[0] - x2.back());

   
    double e23 = abs((dx23 - dx2)/dx2);
/*
    cout << beta2 << " "
         << beta3 << " "
         << e23 << endl;
*/
    if(e23 < target3)
    {
        cout << "\nconverged\n";

        cout << "beta2_2 = " << beta2 << endl;
        cout << "beta3 = " << beta3 << endl;
        cout << "e23 = " << e23 << endl;
        
        break;
    }


    beta2_2 -= beta2_2*step3;

    beta3 += beta3*step3;

    

    if(beta2_2 <= 1.01 || beta3 >= 50)
    {
        cout << "no convergence\n";
        break;
    }
}


double beta5 = 1.001;
double beta4 = 150;

double target2 = 0.003;

double step2 = 0.0001;
while(true)
{

    c0(&x5, dxi[4], Lx[4], beta5, 41, Lx[0]+Lx[1]+Lx[2]+Lx[3]);
    cs(&x4, dxi[3], Lx[3], beta4,91, Lx[0]+Lx[1]+Lx[2]);

    
    vector<double> xtemp = x5;
    xtemp.erase(xtemp.begin());

    
    double dx4  = abs(x4.back() - x4[x4.size()-2]);
    double dx45 = abs(xtemp[0] - x4.back());

   
    double e45 = abs((dx45 - dx4)/dx4);

    /*cout << beta4 << " "
         << beta5 << " "
         << e45 << endl;
*/
    if(e45 < target2)
    {
        cout << "\nconverged\n";

        cout << "beta4 = " << beta4 << endl;
        cout << "beta5 = " << beta5 << endl;
        cout << "e45 = " << e45 << endl;
        
        break;
    }


    beta5 += beta5*step2;

    beta4 -= beta4*step2;

    

    if(beta5 >=2 || beta4<=90)
    {
        cout << "no convergence\n";
        break;
    }
}

vector<double> betasy,Ly;
vector <int> Ny;

ifstream filey("gridy.txt");
double by,ny,ly;
while(filey>>by>>ny>>ly)
{
	betasy.push_back(by);
	Ny.push_back(ny);
	Ly.push_back(ly);
}
filey.close(); 
vector <double> y1(Ny[0]),y2(Ny[1]);

double detai[2];
for(int i =0;i<2;i++)
{
    detai[i]=(double)1/(Ny[i]-1);
}

cs(&y1,detai[0],Ly[0],betasy[0],Ny[0],0);
c0(&y2,detai[1],Ly[1],betasy[1],Ny[1],Ly[0]);

double beta1y = 1.01;
double beta2y = 1.10;

double targety = 0.005 ;
double stepy = 0.0001;

double dy0 = abs(y1.back()-y1[y1.size()-2]);
while(true)
{

    cs(&y1, detai[0], Ly[0], beta1y, Ny[0], 0);
    c0(&y2, detai[1], Ly[1], beta2y, Ny[1], Ly[0]);

    
    vector<double> ytemp = y2;
    ytemp.erase(ytemp.begin());

    
    double dy1  = abs(y1.back() - y1[y1.size()-2]);
    double dy12 = abs(ytemp[0] - y1.back());

   
    double e12 = abs((dy12 - dy1)/dy1);

    /*cout << beta1y << " "
         << beta2y << " "
         << e12 << endl;
*/
    if(e12 < targety)
    {
        cout << "\nconverged\n";

        cout << "beta1y = " << beta1y << endl;
        cout << "beta2y = " << beta2y << endl;
        cout << "e12y = " << e12 << endl;
        
        break;
    }


    beta1y += beta1y*stepy;

    beta2y -= beta2y*stepy;

    

    if(beta2y <= 1.0 || beta1y>=2)
    {
        cout << "no convergence\n";
        break;
    }
  }

	return 0;
}
