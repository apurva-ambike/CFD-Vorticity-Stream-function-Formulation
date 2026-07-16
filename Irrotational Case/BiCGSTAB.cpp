#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>

using namespace std;

long double as(int i, int j, vector<vector<double>>& Y)
{   if(j>0)
    return 1.0/((Y[i][j]-Y[i][j-1])*(Y[i][j+1]-Y[i][j-1]));

    else if(j==0)
    return 0.0;

    else{
    cerr<<"index out of bounds";
    return 0.0;}
}

long double an(int i, int j, vector<vector<double>>& Y)
{   if(j>0)
    return 1.0/((Y[i][j+1]-Y[i][j])*(Y[i][j+1]-Y[i][j-1]));

    else if(j==0)
    return 0.0;

    else{
    cerr<<"index out of bounds";
    return 0.0;
    }
}

long double aw(int i, int j, vector<vector<double>>& X)
{
    if(i>0)
    return 1.0/((X[i][j]-X[i-1][j])*(X[i+1][j]-X[i-1][j]));

    else if(i==0)
    return 0.0;

    else{
    cerr<<"index out of bounds";
    return 0.0;}
}  

long double ae(int i, int j, vector<vector<double>>& X)
{
    if(i>0)
    return 1.0/((X[i+1][j]-X[i][j])*(X[i+1][j]-X[i-1][j]));

    else if(i==0)
    return 0.0;

    else{
    cerr<<"index out of bounds";
    return 0.0;}
}

long double ap(int i, int j, vector<vector<double>>& X, vector <vector <double>>& Y)
{
    double f1 = -1.0/((X[i+1][j]-X[i][j])*(X[i][j]-X[i-1][j]));

    double f2 = -1.0/((Y[i][j+1]-Y[i][j])*(Y[i][j]-Y[i][j-1]));

    return f1 + f2;
}

vector <vector<double>> createXcoords(vector <double>& xc, vector <double>& yc)
{
    int Nx = xc.size();
    int Ny = yc.size();

    vector <vector<double>> X(Nx, vector <double> (Ny));
    for(int i =0;i<Nx;i++)
    {
        for(int j =0;j<Ny;j++)
        {
            X[i][j] = xc[i];
        }
    }
    return X;
}

vector <vector<double>> createYcoords(vector <double>& xc, vector <double>& yc)
{
    int Nx = xc.size();
    int Ny = yc.size();
    vector <vector<double>> Y(Nx, vector <double> (Ny));
    for(int i =0;i<Nx;i++)
    {
        for(int j =0;j<Ny;j++)
        {
            Y[i][j] = yc[j];
        }
    }
    return Y;
}


long double VectorDot(const vector <vector<long double>>& x, const vector <vector<long double>> & y)
{   long double result = 0;
    if( x.size() == y.size() && x[0].size() == y[0].size())
    {
        for(int i=0; i< x.size(); i++)
        {
            for(int j =0; j<x[0].size();j++)
            {
                result += x[i][j]*y[i][j];
            }
        }

    return result;
    }

    else
    {
        cerr<<"vectors of unequal size cannot undergo dot product"<<endl;
        return 0.0;
    }
}

int isInsideBlock(double l1, double l2, double y1, int i, int j, vector<vector<double>>& X, vector<vector<double>>& Y)
{
    if(X[i][j]/l1>=1 && X[i][j]/l2 <= 1 && Y[i][j]/y1 <= 1)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

vector <vector<long double>> VectorSubtraction(const vector<vector<long double>>& X, const vector<vector<long double>>& Y)
{
    vector <vector<long double>> res (X.size(),vector<long double>(X[0].size(),0.0));

    if(X.size() == Y.size() && X[0].size() == Y[0].size())
    {
        for(int i =0;i<X.size();i++)
        {
            for(int j =0; j<X[0].size();j++)
            {
                res[i][j] = X[i][j] - Y[i][j];
            }
        }

        return res;
    }
    
    else
    {
        cerr<<"vectors have unequal sizes";
        return {};
    }
}


vector <vector<long double>> ScalarMultiplication(const vector<vector<long double>> &X, const long double alpha)
{
    vector <vector<long double>> res(X.size(),vector<long double>(X[0].size(),0.0));
    for(int i =0;i<X.size();i++)
        {
            for(int j =0; j<X[0].size();j++)
            {
                res[i][j] = alpha*X[i][j];
            }
        }
    return res;
}


class SparseMatrix
{public:
    int nx;
    int ny;
    vector <vector<vector<long double > > > a;

SparseMatrix(int r, int c)
{   nx = r;
    ny = c;
    a.resize(5,vector<vector<long double>>(nx,vector<long double>(ny,0.0)));
}

void setValues(vector <vector<double>>& X, vector<vector<double>>& Y)
{
   int Nx = X.size();
   int Ny = X[0].size();

   for(int i =1; i<Nx-1; i++)
   {
    for(int j = 1; j<Ny-1;j++)
    {
        a[0][i][j] = ap(i,j,X,Y);
        a[1][i][j] = ae(i,j,X);
        a[2][i][j] = as(i,j,Y);
        a[3][i][j] = aw(i,j,X);
        a[4][i][j] = an(i,j,Y);
    }
   }

   for(int k =0;k<Nx;k++)
   {
    a[0][k][0] = 1.0;
    a[0][k][Ny-1] = 1.0;
   }

   for(int l =0;l<Ny;l++)
   {
    a[0][0][l] =1.0;
    a[0][Nx-1][l] = 1.0;
   }
}

void setBoundaryConditions(double l1, double l2,double l3, double l4,double y1, vector<vector<double>>& X, vector<vector<double>>& Y)
{
    int Nx = X.size();
    int Ny = X[0].size();

//finding indices of boundary points:
int i_1=-1,i_2=-1,i_3=-1,i_4=-1, j_1=-1;
int i,j;

for(i =1; i<Nx-1;i++)
{
    if(abs(X[i][1]/l1 -1.0) < 1e-8)
    {
        i_1 = i;
    }

    if(abs(X[i][1]/l2 -1.0) < 1e-8)
    {
        i_2 = i;
    }

    if(abs(X[i][1]/l3 - 1.0) < 1e-8)
    {
        i_3 = i;
    }

    if(abs(X[i][1]/l4 -1.0) < 1e-8)
    {
        i_4 = i;
        break;
    }
}

for(int k =1;k<Ny-1;k++)
{
    if(abs(Y[1][k]/y1 - 1.0) < 1e-8)
    {
        j_1 = k;
        break;
    }
}


//std::cout<<i_1<<" "<<i_2<<" "<<i_3<<" "<<i_4<<" "<<j_1<<endl;
if(i_1 == -1 || i_2 == -1 || i_3 == -1 || i_4 == -1 || j_1 == -1)
{
    cerr << "Boundary index not found\n";
    return;
}

    for(int p =1; p<Ny-1; p++)
    {
        a[3][1][p] = 0;
        a[0][Nx-2][p] += a[1][Nx-2][p];
        a[1][Nx-2][p] = 0;
    }

    for(int q =1; q<Nx-1;q++)
    {
        a[4][q][Ny-2] = 0;
    }

    for(int l =1; l<Nx-1; l++)
    {
        for(int m =1; m<Ny-1; m++)
        {
            if(isInsideBlock(l1,l2,y1,l,m,X,Y) == 1 || isInsideBlock(l3,l4,y1,l,m,X,Y)== 1)
            {
                a[0][l][m] = 1.0;
                a[1][l][m]=0;
                a[2][l][m]=0;
                a[3][l][m]=0;
                a[4][l][m] = 0;
            }
       }    
    }

for(i =1;i<i_1;i++)
{
    a[2][i][1] = 0;
}

for(j = 1; j<j_1+1;j++)
{
    a[1][i_1-1][j] = 0;
}

for(i= i_1+1;i<i_2+1;i++)
{
    a[2][i][j_1+1] = 0;
}

for(j = 1; j<j_1+1;j++)
{
    a[3][i_2+1][j] = 0;
}

for(i =i_2+1;i<i_3;i++)
{
    a[2][i][1] = 0;
}

for(j = 1; j<j_1+1;j++)
{
    a[1][i_3-1][j] = 0;
}

for(i = i_3+1;i<i_4+1;i++)
{
    a[2][i][j_1+1] = 0;
}

for(j = 1; j<j_1+1;j++)
{
    a[3][i_4+1][j] = 0;
}

for(i = i_4+1; i<Nx-1;i++)
{
    a[2][i][1] = 0;
}

}

vector <vector <long double>> MatrixVectorMultiplication(vector <vector<long double>>& phi)
{       int Nx = nx;
        int Ny = ny;

    if(phi.size() == nx && phi[0].size() == ny)
    {   
        vector <vector<long double>> result (Nx, vector<long double>(Ny,0.0));
        
        for(int i =1; i<Nx-1;i++)
        {
            for(int j =1; j<Ny-1;j++)
            {    
                result[i][j] = a[0][i][j]*phi[i][j] + a[1][i][j]*phi[i+1][j] +a[2][i][j]*phi[i][j-1] +a[3][i][j]*phi[i-1][j] +a[4][i][j]*phi[i][j+1] ;
            }
        }

        vector <vector<long double>> interior(Nx-2, vector<long double>(Ny-2));
        for(int k =1;k<Nx-1;k++)
        {
            for(int l =1;l<Ny-1;l++)
            {
                interior[k-1][l-1] = result[k][l];
            }
        }
        return interior;
    }

    if(phi.size()==nx-2 && phi[0].size()==ny-2)
    {
        vector <vector<long double>> result2 (Nx-2, vector<long double>(Ny-2,0.0));
       for(int i = 1; i < nx-1; i++)
    {
        for(int j = 1; j < ny-1; j++)
        {
            int I = i - 1;
            int J = j - 1;

            result2[I][J] += a[0][i][j] * phi[I][J];

            if(i < nx-2)
                result2[I][J] += a[1][i][j] * phi[I+1][J]; 

            if(i > 1)
                result2[I][J] += a[3][i][j] * phi[I-1][J]; 

            if(j > 1)
                result2[I][J] += a[2][i][j] * phi[I][J-1]; 

            if(j < ny-2)
                result2[I][J] += a[4][i][j] * phi[I][J+1]; 
            }
        }

    return result2;
    }

    else
    {
        cerr<<"matrix and vector of these sizes cannot be multiplied"<<endl;
        return {};
    }
}

double Residual(vector<vector<long double>>& b,vector<vector<long double>>& psi )
{
    vector <vector<long double>> psi_new = MatrixVectorMultiplication(psi);
    vector <vector<long double>> temp = VectorSubtraction(b,psi_new);
    double temp2 = VectorDot(temp,temp);
    double res = sqrtl(temp2);
    return res;
}

};

int main()
{  
    vector <double> xc,yc;

    ifstream filex("xcoords_3.csv");
    double x;

    while(filex>>x)
    {
        xc.push_back(x);
    }
    filex.close();

    ifstream filey("ycoords.csv");
    double y;

    while(filey>>y)
    {
        yc.push_back(y);
    }
    filey.close();

    vector <vector<double> > X = createXcoords(xc,yc);
    vector <vector<double> > Y = createYcoords(xc,yc);

    int Nx = xc.size();
    int Ny = yc.size();

    SparseMatrix A(Nx,Ny);
    A.setValues(X,Y);
double u_inf = 1.0;

vector <vector<long double>> psi(Nx, vector<long double>(Ny,0.0));

for(int i = 0; i < Nx; i++) //bottom
{
     psi[i][0] = 0.0;
}

for(int j = 0; j < Ny; j++)  //left
{
    psi[0][j] = u_inf * Y[0][j];
}

for(int i = 0; i < Nx; i++) 
{
    psi[i][Ny-1] = psi[0][Ny-1]; //top
}

vector <vector<long double>> b(Nx-2, vector<long double>(Ny-2,0.0));

  for(int l = 1; l<Ny-1;l++)
{
    b[0][l-1] += 0 - aw(1,l,X)*psi[0][l];
}

for(int m =1; m<Nx-1;m++)
{
    b[m-1][Ny-3] += 0 - an(m,Ny-2,Y)*psi[m][Ny-1];
}

//std::cout<<"b vector created"<<endl;
A.setBoundaryConditions(5.0,6.0,8.0,9.0,1.0,X,Y);

double tol =1e-8;
int count =0;
int count_max = 1e4;

vector<vector<long double>> r = VectorSubtraction(b,A.MatrixVectorMultiplication(psi));
vector<vector<long double>> r_0 = r;
long double alpha =0, beta =0,ksi =0;
long double res = A.Residual(b,psi);
vector <vector<long double>> p = r_0;

ofstream file2("residualBiCGSTAB.csv");

while(res>tol && count<count_max)
{
    vector <vector<long double>> temp  = VectorSubtraction(p,ScalarMultiplication(A.MatrixVectorMultiplication(p),ksi));
    for(int i =0;i<p.size();i++)
    {
        for(int j =0;j<p[0].size();j++)
        {
            p[i][j] = r[i][j] + beta*temp[i][j];
        }
    }

    alpha = (VectorDot(r,r_0))/(VectorDot(A.MatrixVectorMultiplication(p),r_0));

    vector<vector<long double>> z = VectorSubtraction(r,ScalarMultiplication(A.MatrixVectorMultiplication(p),alpha));

    ksi = (VectorDot(A.MatrixVectorMultiplication(z),z))/(VectorDot(A.MatrixVectorMultiplication(z),A.MatrixVectorMultiplication(z)));

    for(int i = 1;i<Nx-1;i++)
    {
        for(int j =1;j<Ny-1;j++)
        {
            psi[i][j] = psi[i][j] + alpha*p[i-1][j-1] + ksi*z[i-1][j-1];
        }
    }

    double rho = (VectorDot(r_0,r));

    r = VectorSubtraction(z,ScalarMultiplication(A.MatrixVectorMultiplication(z),ksi));

    beta = (alpha*(VectorDot(r_0,r)))/(ksi*rho);

    count ++;
    res = A.Residual(b,psi);

    std::cout<<count<<" "<<res<<endl;
    file2<<log10(count)<<","<<log10(res)<<","<<count<<","<<res<<endl;

    for(int j = 1; j < Ny-1; j++)
    {
        psi[Nx-1][j] = psi[Nx-2][j];
    }
    
}

if(res<= tol)
{
    std::cout<<"converged in "<<count<<" steps"<<endl;
}

else
{
    std::cout <<"no convergence"<<endl;
}
ofstream filer("resultBiCGSTAB.csv");

for(int i = 0; i < Nx; i++)
{
    for(int j = 0; j < Ny; j++)
    {
        filer << psi[i][j];
        if(j < Ny-1) filer << ",";
    }
    filer << "\n";
}

filer.close();
return 0;
}