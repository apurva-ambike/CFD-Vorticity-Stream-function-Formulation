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

x2.erase(x2.begin());
x3.erase(x3.begin());
x4.erase(x4.begin());
x5.erase(x5.begin());
  
vector<vector<double>> all = {x1, x2, x3, x4, x5};
vector<double> resultx;

for (const auto& v : all) 
{
    resultx.insert(resultx.end(), v.begin(), v.end());
}
  
vector < double> delx(resultx.size());

for(int i =0; i<resultx.size()-1;i++)
{
	delx[i] = resultx[i+1] -resultx[i];
}
cout<<error(x1,x2)*100<<endl;
cout<<error(x2,x3)*100<<endl;
cout<<error(x3,x4)*100<<endl;
cout<<error(x4,x5)*100<<endl;

writeVectorToCSV(resultx,"xcoords_3.csv");
writeVectorToCSV(delx,"deltax_3.csv");

y2.erase(y2.begin());

vector<double> resulty;
resulty.insert(resulty.end(), y1.begin(), y1.end());
resulty.insert(resulty.end(), y2.begin(), y2.end());

vector < double> dely(resulty.size());

for(int i =0; i<resulty.size()-1;i++)
{
	dely[i] = resulty[i+1] -resulty[i];
}

writeVectorToCSV(resulty,"ycoords.csv");
writeVectorToCSV(dely,"deltay.csv");

cout<<error(y1,y2)*100;

	return 0;
}
