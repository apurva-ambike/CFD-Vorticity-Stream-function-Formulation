#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <chrono>
using namespace std;
using namespace std::chrono;

void WriteFieldCSV(const string& filename,
                   const vector<vector<double>>& phi)
{
    ofstream file(filename);

    for(int i = 0; i < phi.size(); i++)
    {
        for(int j = 0; j < phi[0].size(); j++)
        {
            file << phi[i][j];

            if(j != phi[0].size()-1)
                file << ",";
        }
        file << "\n";
    }

    file.close();
}

vector <vector<double>> VectorSubtraction(const vector<vector<double>>& X, const vector<vector<double>>& Y)
{
    vector <vector<double>> res (X.size(),vector<double>(X[0].size(),0.0));

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

vector <vector<double>> ScalarMultiplication(const vector<vector<double>> &X, const double alpha)
{
    vector <vector<double>> res (X.size(),vector<double>(X[0].size(),0.0));
    for(int i =0;i<X.size();i++)
        {
            for(int j =0; j<X[0].size();j++)
            {
                res[i][j] = alpha*X[i][j];
            }
        }
    return res;
}

double VectorDot(const vector <vector<double>>& x, const vector <vector<double>> & y)
{   double result = 0;
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

vector <int> boundarypoints(const vector <vector<double>> &X, const vector<vector<double>> &Y,double l1, double l2,double l3, double l4,double y1)
{
    int Nx = X.size();
    int Ny = X[0].size();

    vector<int> I(5,-1);

//finding indices of boundary points:
int i,j;

for(i = 1; i<Nx-1;i++)
{
    if(abs(X[i][1]/l1 -1.0) < 1e-8)
    {
        I[0] = i;
    }

    if(abs(X[i][1]/l2 -1.0) < 1e-8)
    {
        I[1] = i;
    }

    if(abs(X[i][1]/l3 - 1.0) < 1e-8)
    {
        I[2] = i;
    }

    if(abs(X[i][1]/l4 -1.0) < 1e-8)
    {
        I[3] = i;
        break;
    }
}

for(int k =1;k<Ny-1;k++)
{
    if(abs(Y[1][k]/y1 - 1.0) < 1e-8)
    {
        I[4] = k;
        break;
    }
}

if(I[0] == -1 || I[1] == -1 || I[2] == -1 || I[3] == -1 || I[4] == -1)
{
    cerr << "Boundary index not found\n";
    return {};
}

else{
    return I;
    }
}



class SparseMatrix
{public:
    int nx;
    int ny;
    vector <vector<vector<double > > > a;

SparseMatrix(int r, int c)
{   nx = r;
    ny = c;
    a.resize(5,vector<vector<double>>(nx,vector<double>(ny,0.0)));
}

void setValuesPsi(vector <vector<double>>& X, vector<vector<double>>& Y)
{
   int Nx = X.size();
   int Ny = X[0].size();

   for(int i =1; i<Nx-1; i++)
   {
    for(int j = 1; j<Ny-1;j++)
    {
        a[0][i][j] = -(((Y[i][j+1] - Y[i][j-1])*(X[i+1][j] - X[i-1][j]))/((X[i][j]-X[i-1][j])*(X[i+1][j]-X[i][j])) + ((Y[i][j+1] -Y[i][j-1])*(X[i+1][j] - X[i-1][j]))/((Y[i][j+1]-Y[i][j])*(Y[i][j]-Y[i][j-1])));
        a[1][i][j] = (Y[i][j+1]-Y[i][j-1])/(X[i+1][j] - X[i][j]);
        a[2][i][j] = (X[i+1][j]-X[i-1][j])/(Y[i][j] - Y[i][j-1]);
        a[3][i][j] = (Y[i][j+1]-Y[i][j-1])/(X[i][j] - X[i-1][j]);
        a[4][i][j] = (X[i+1][j]-X[i-1][j])/(Y[i][j+1]-Y[i][j]);
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

void setValuesOmega(vector <vector<double>>& X, vector<vector<double>>& Y, double delta_t, double Re, vector<vector<double>> &u,vector<vector<double>> &v )
{
   int Nx = X.size();
   int Ny = X[0].size();

   for(int i =1; i<Nx-1; i++)
   {
    for(int j = 1; j<Ny-1;j++)
    {
        a[0][i][j] = (1.0+(delta_t*u[i][j]*(X[i+1][j] + X[i-1][j] - 2.0*X[i][j]))/((X[i+1][j]- X[i][j])*(X[i][j]-X[i-1][j])) + (delta_t*v[i][j]*(Y[i][j-1] + Y[i][j+1] - 2.0*Y[i][j]))/((Y[i][j+1]- Y[i][j])*(Y[i][j]-Y[i][j-1])) + (2.0*delta_t/Re)*(1.0/((X[i][j]-X[i-1][j])*(X[i+1][j]-X[i][j]))+ 1.0/((Y[i][j]-Y[i][j-1])*(Y[i][j+1]-Y[i][j]))));

        a[1][i][j] = (delta_t/((X[i+1][j]-X[i-1][j])*(X[i+1][j]-X[i][j])))*(u[i][j]*(X[i][j]-X[i-1][j]) - (2.0/Re));

        a[2][i][j] = (-delta_t/((Y[i][j+1]-Y[i][j-1])*(Y[i][j]-Y[i][j-1])))*(v[i][j]*(Y[i][j+1]-Y[i][j]) + (2.0/Re));

        a[3][i][j] = (-delta_t/((X[i+1][j]-X[i-1][j])*(X[i][j] - X[i-1][j])))*(u[i][j]*(X[i+1][j]-X[i][j]) + (2.0/Re));

        a[4][i][j] = (delta_t/((Y[i][j+1]-Y[i][j-1])*(Y[i][j+1]-Y[i][j])))*(v[i][j]*(Y[i][j]-Y[i][j-1]) - (2.0/Re));
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

    vector <int> I = boundarypoints(X,Y,l1,l2,l3,l4,y1);
    int i_1 = I[0], i_2 = I[1], i_3 = I[2], i_4 = I[3], j_1 = I[4];
    int i,j;

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

for(i= i_1;i<=i_2;i++)
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

for(i = i_3;i<=i_4;i++)
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

vector <vector <double>> MatrixVectorMultiplication(vector <vector<double>>& phi)
{       int Nx = nx;
        int Ny = ny;

    if(phi.size() == nx && phi[0].size() == ny)
    {   
        vector <vector<double>> result (Nx, vector<double>(Ny,0.0));
        
        for(int i =1; i<Nx-1;i++)
        {
            for(int j =1; j<Ny-1;j++)
            {    
                result[i][j] = a[0][i][j]*phi[i][j] + a[1][i][j]*phi[i+1][j] +a[2][i][j]*phi[i][j-1] +a[3][i][j]*phi[i-1][j] +a[4][i][j]*phi[i][j+1] ;
            }
        }

        vector <vector<double>> interior(Nx-2, vector<double>(Ny-2));
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
        vector <vector<double>> result2 (Nx-2, vector<double>(Ny-2,0.0));
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

double Residual(vector<vector<double>>& b,vector<vector<double>>& psi )
{
    vector <vector<double>> psi_new = MatrixVectorMultiplication(psi);
    vector <vector<double>> temp = VectorSubtraction(b,psi_new);
    double temp2 = VectorDot(temp,temp);
    double res = sqrt(temp2);
    return res;
}

};

class StoneCoeffs
{
 private:
 int nx;
 int ny;
 double alpha;
 vector <vector<vector<double>>> L;
 vector <vector<vector<double>>> U;
 public:
 StoneCoeffs(int r, int c, double a)
 {
    nx =r;
    ny =c;
    alpha =a;
    L.resize(3,vector<vector<double>>(nx,vector<double>(ny,0.0)));
    U.resize(2,vector<vector<double>>(nx,vector<double>(ny,0.0)));
 }

 void SetCoeffs(vector<vector<double>> &X, vector<vector<double>> &Y, SparseMatrix &A)
 {
    int Nx =X.size();
    int Ny = X[0].size();

    for(int i =1;i<Nx-1;i++)
    {
        for(int j =1;j<Ny-1;j++)
        {
            L[2][i][j] =A.a[2][i][j]/(1+alpha*(U[0][i][j-1])); //Ls

            L[1][i][j] =  A.a[3][i][j]/(1+alpha*(U[1][i-1][j])); //Lw
            
            L[0][i][j] = A.a[0][i][j] +alpha*L[1][i][j]*U[1][i-1][j] - L[1][i][j]*U[0][i-1][j] + alpha*L[2][i][j]*U[0][i][j-1] - L[2][i][j]*U[1][i][j-1]; //Lp

            U[0][i][j] = (A.a[1][i][j] - alpha*(L[2][i][j]*U[0][i][j-1]))/L[0][i][j]; //Ue

            U[1][i][j] = (A.a[4][i][j] - alpha*(L[1][i][j]*U[1][i-1][j]))/L[0][i][j]; //Un
        }
    }
 }

vector<vector<double>> M_Inverse( const vector<vector<double>> &p)
{
    vector<vector<double>> q (p.size(),vector<double>(p[0].size(),0.0));
    vector<vector<double>> y (p.size(),vector<double>(p[0].size(),0.0));

    //Forward Sweep
    double sum =0;
    int i=0,j=0;
    for( i=0;i<p.size();i++)
    {
        for(j=0;j<p[0].size();j++)
        {   int I = i+1;
            int J =j+1;
            sum = 0;
            if(j>0)
            sum += L[2][I][J]*q[i][j-1];

            if(i>0)
            sum += L[1][I][J]*q[i-1][j];

            q[i][j] = (p[i][j] -sum)/L[0][I][J];
        }
    }

    //Backward sweep
    sum = 0;
    for(i = p.size()-1;i >=0;i--)
    {
        for(j = p[0].size()-1;j>=0;j--)
        {   int I =i+1;
            int J =j+1;
            sum =0;
            if(i<p.size()-1)
            sum += U[0][I][J]*y[i+1][j];

            if(j<p[0].size()-1)
            sum += U[1][I][J]*y[i][j+1];

            y[i][j] = (q[i][j] - sum);
        }
    }
    
    return y;
}
};


void UpdateOmegaAndb_omegaAtBoundaries(vector<vector<double>> &omega, vector<vector<double>> &b_omega, vector<vector<double>> &psi, SparseMatrix &A_omega_raw, vector<vector<double>> &X, vector<vector<double>> &Y, double u_inf)
{

    vector<int> id = boundarypoints(X,Y,5.0,6.0,8.0,9.0,1.0);
    int Nx = X.size();
    int Ny = X[0].size();

    int i,j;

       for( j = 1; j < Ny-1; j++)  //left boundary
    {
        omega[0][j] = 0.0;
        b_omega[0][j-1] -=   A_omega_raw.a[3][1][j]*omega[0][j];
    }   
    
    omega[0][Ny-1] = 0.0;

      for( i = 1; i < Nx-1; i++) 
    {
        omega[i][Ny-1] = 0.0; //top boundary
        b_omega[i-1][Ny-3] -= A_omega_raw.a[4][i][Ny-2]*omega[i][Ny-1];
    }

        for(i =1;i<id[0];i++)
    {
        omega[i][0] = (-2.0/((Y[i][1]-Y[i][0])*(Y[i][1]-Y[i][0])))*(psi[i][1]-psi[i][0]);
        b_omega[i-1][0] -= A_omega_raw.a[2][i][1]*omega[i][0];
    }

    omega[0][0] =0.0; 

    for(j = 1; j<id[4];j++)
    {
        omega[id[0]][j] = (-2.0/((X[id[0]][j]-X[id[0]-1][j])*(X[id[0]][j]-X[id[0]-1][j])))*(psi[id[0]-1][j]-psi[id[0]][j]);
        b_omega[id[0]-2][j-1] -= A_omega_raw.a[1][id[0]-1][j] * omega[id[0]][j];

    }

    omega[id[0]][0] = 0.5*(((-2.0/((Y[id[0]][1]-Y[id[0]][0])*(Y[id[0]][1]-Y[id[0]][0])))*(psi[id[0]][1]-psi[id[0]][0]))+((-2.0/((X[id[0]][0]-X[id[0]-1][0])*(X[id[0]][0]-X[id[0]-1][0])))*(psi[id[0]-1][0]-psi[id[0]][0])));
    
    omega[id[0]][id[4]] = 0.5*(((-2.0/((X[id[0]][id[4]]-X[id[0]-1][id[4]])*(X[id[0]][id[4]]-X[id[0]-1][id[4]])))*(psi[id[0]-1][id[4]]-psi[id[0]][id[4]]))+((-2.0/((Y[id[0]][id[4]+1]-Y[id[0]][id[4]])*(Y[id[0]][id[4]+1]-Y[id[0]][id[4]])))*(psi[id[0]][id[4]+1]-psi[id[0]][id[4]])));
    b_omega[id[0]-2][id[4]-1] -= A_omega_raw.a[1][id[0]-1][id[4]]*omega[id[0]][id[4]];
    b_omega[id[0]-1][id[4]]  -=  A_omega_raw.a[2][id[0]][id[4]+1]*omega[id[0]][id[4]];

    for(i= id[0]+1;i<id[1];i++)
    {
        omega[i][id[4]] = (-2.0/((Y[i][id[4]+1]-Y[i][id[4]])*(Y[i][id[4]+1]-Y[i][id[4]])))*(psi[i][id[4]+1]-psi[i][id[4]]);
        b_omega[i-1][id[4]] -= A_omega_raw.a[2][i][id[4]+1] * omega[i][id[4]];
    }

    for(j = 1; j<id[4];j++)
    {
        omega[id[1]][j] = (-2.0/((X[id[1]][j]-X[id[1]+1][j])*(X[id[1]][j]-X[id[1]+1][j])))*(psi[id[1]+1][j]-psi[id[1]][j]);
        b_omega[id[1]][j-1] -= A_omega_raw.a[3][id[1]+1][j] * omega[id[1]][j];
    }

    omega[id[1]][id[4]] = 0.5*(((-2.0/((Y[id[1]][id[4]+1]-Y[id[1]][id[4]])*(Y[id[1]][id[4]+1]-Y[id[1]][id[4]])))*(psi[id[1]][id[4]+1]-psi[id[1]][id[4]]))+((-2.0/((X[id[1]][id[4]]-X[id[1]+1][id[4]])*(X[id[1]][id[4]]-X[id[1]+1][id[4]])))*(psi[id[1]+1][id[4]]-psi[id[1]][id[4]])));
    
    omega[id[1]][0] = 0.5*(((-2.0/((X[id[1]][0]-X[id[1]+1][0])*(X[id[1]][0]-X[id[1]+1][0])))*(psi[id[1]+1][0]-psi[id[1]][0]))+((-2.0/((Y[id[1]][1]-Y[id[1]][0])*(Y[id[1]][1]-Y[id[1]][0])))*(psi[id[1]][1]-psi[id[1]][0])));

    b_omega[id[1]][id[4]-1] -= A_omega_raw.a[3][id[1]+1][id[4]]*omega[id[1]][id[4]];
    b_omega[id[1]-1][id[4]] -=A_omega_raw.a[2][id[1]][id[4]+1]*omega[id[1]][id[4]];

    for(i =id[1]+1;i<id[2];i++)
    {
        omega[i][0] = (-2.0/((Y[i][1]-Y[i][0])*(Y[i][1]-Y[i][0])))*(psi[i][1]-psi[i][0]);
        b_omega[i-1][0] -= A_omega_raw.a[2][i][1] * omega[i][0];
    }

    for(j = 1; j<id[4];j++)
    {
        omega[id[2]][j] =  (-2.0/((X[id[2]][j]-X[id[2]-1][j])*(X[id[2]][j]-X[id[2]-1][j])))*(psi[id[2]-1][j]-psi[id[2]][j]);
        b_omega[id[2]-2][j-1] -= A_omega_raw.a[1][id[2]-1][j] * omega[id[2]][j];
    }

    omega[id[2]][0] = 0.5*(((-2.0/((Y[id[2]][1]-Y[id[2]][0])*(Y[id[2]][1]-Y[id[2]][0])))*(psi[id[2]][1]-psi[id[2]][0]))+((-2.0/((X[id[2]][0]-X[id[2]-1][0])*(X[id[2]][0]-X[id[2]-1][0])))*(psi[id[2]-1][0]-psi[id[2]][0])));

    omega[id[2]][id[4]] = 0.5*(((-2.0/((X[id[2]][id[4]]-X[id[2]-1][id[4]])*(X[id[2]][id[4]]-X[id[2]-1][id[4]])))*(psi[id[2]-1][id[4]]-psi[id[2]][id[4]]))+((-2.0/((Y[id[2]][id[4]+1]-Y[id[2]][id[4]])*(Y[id[2]][id[4]+1]-Y[id[2]][id[4]])))*(psi[id[2]][id[4]+1]-psi[id[2]][id[4]])));

    b_omega[id[2]-2][id[4]-1] -= A_omega_raw.a[1][id[2]-1][id[4]]*omega[id[2]][id[4]];
    b_omega[id[2]-1][id[4]] -= A_omega_raw.a[2][id[2]][id[4]+1]*omega[id[2]][id[4]];

    for(i = id[2]+1;i<id[3];i++)
    {
        omega[i][id[4]] = (-2.0/((Y[i][id[4]+1]-Y[i][id[4]])*(Y[i][id[4]+1]-Y[i][id[4]])))*(psi[i][id[4]+1]-psi[i][id[4]]);
        b_omega[i-1][id[4]] -= A_omega_raw.a[2][i][id[4]+1] * omega[i][id[4]];
    }

    for(j = 1; j<id[4];j++)
    {
         omega[id[3]][j] = (-2.0/((X[id[3]][j]-X[id[3]+1][j])*(X[id[3]][j]-X[id[3]+1][j])))*(psi[id[3]+1][j]-psi[id[3]][j]);
         b_omega[id[3]][j-1] -= A_omega_raw.a[3][id[3]+1][j] * omega[id[3]][j];

    }

    omega[id[3]][id[4]] = 0.5*(((-2.0/((Y[id[3]][id[4]+1]-Y[id[3]][id[4]])*(Y[id[3]][id[4]+1]-Y[id[3]][id[4]])))*(psi[id[3]][id[4]+1]-psi[id[3]][id[4]]))+((-2.0/((X[id[3]][id[4]]-X[id[3]+1][id[4]])*(X[id[3]][id[4]]-X[id[3]+1][id[4]])))*(psi[id[3]+1][id[4]]-psi[id[3]][id[4]])));

    omega[id[3]][0] = 0.5*(((-2.0/((X[id[3]][0]-X[id[3]+1][0])*(X[id[3]][0]-X[id[3]+1][0])))*(psi[id[3]+1][0]-psi[id[3]][0]))+((-2.0/((Y[id[3]][1]-Y[id[3]][0])*(Y[id[3]][1]-Y[id[3]][0])))*(psi[id[3]][1]-psi[id[3]][0])));

    b_omega[id[3]][id[4]-1] -= A_omega_raw.a[3][id[3]+1][id[4]]*omega[id[3]][id[4]];
    b_omega[id[3]-1][id[4]] -= A_omega_raw.a[2][id[3]][id[4]+1]*omega[id[3]][id[4]];

    for(i = id[3]+1; i<Nx-1;i++)
    {
        omega[i][0] = (-2.0/((Y[i][1]-Y[i][0])*(Y[i][1]-Y[i][0])))*(psi[i][1]-psi[i][0]);
        b_omega[i-1][0] -= A_omega_raw.a[2][i][1] * omega[i][0];

    }



b_omega[id[0]-2][0] -= A_omega_raw.a[1][id[0]-1][1]*omega[id[0]][0];  // left neighbor
b_omega[id[0]-1][0] -= A_omega_raw.a[2][id[0]][1]*omega[id[0]][0];    // above neighbor

b_omega[id[1]][0]   -= A_omega_raw.a[3][id[1]+1][1]*omega[id[1]][0];  // right neighbor
b_omega[id[1]-1][0] -= A_omega_raw.a[2][id[1]][1]*omega[id[1]][0];    // above neighbor


b_omega[id[2]-2][0] -= A_omega_raw.a[1][id[2]-1][1]*omega[id[2]][0];  // left neighbor
b_omega[id[2]-1][0] -= A_omega_raw.a[2][id[2]][1]*omega[id[2]][0];    // above neighbor


b_omega[id[3]][0]   -= A_omega_raw.a[3][id[3]+1][1]*omega[id[3]][0];  // right neighbor
b_omega[id[3]-1][0] -= A_omega_raw.a[2][id[3]][1]*omega[id[3]][0];    // above neighbor

}


void Update_u_v(vector<vector<double>> &u, vector<vector<double>> &v, vector <vector<double>> &X, vector<vector<double>> &Y, vector<vector<double>> &psi)
{
    int Nx = X.size();
    int Ny = X[0].size();

    for(int i =1 ;i<Nx-1;i++)
    {
        for(int j =1; j<Ny-1;j++)
        {
            u[i][j] = (psi[i][j]*(Y[i][j+1]+ Y[i][j-1] -2.0*Y[i][j])/((Y[i][j+1]-Y[i][j])*(Y[i][j]-Y[i][j-1]))) + (psi[i][j+1]*(Y[i][j]-Y[i][j-1])/(Y[i][j+1]-Y[i][j]) - psi[i][j-1]*(Y[i][j+1]-Y[i][j])/(Y[i][j]-Y[i][j-1]))/(Y[i][j+1]-Y[i][j-1]);
            v[i][j] = -((psi[i][j]*(X[i+1][j] + X[i-1][j] - 2.0*X[i][j])/((X[i+1][j]-X[i][j])*(X[i][j]-X[i-1][j]))) +  (psi[i+1][j]*(X[i][j]-X[i-1][j])/(X[i+1][j]-X[i][j]) - psi[i-1][j]*(X[i+1][j]-X[i][j])/(X[i][j]-X[i-1][j]))/(X[i+1][j]- X[i-1][j]));
        }
    }

}

void ApplyVelocityBC(vector<vector<double>>& u,
                     vector<vector<double>>& v,
                     vector<vector<double>>& X,
                     vector<vector<double>>& Y,
                     double u_inf)
{
    int Nx = X.size();
    int Ny = X[0].size();

    vector<int> id = boundarypoints(X, Y, 5.0, 6.0, 8.0, 9.0, 1.0);

    // left inlet
    for(int j = 0; j < Ny; j++)
    {
        u[0][j] = u_inf;
        v[0][j] = 0.0;
    }

    // top boundary
    for(int i = 0; i < Nx; i++)
    {
        u[i][Ny-1] = u_inf;
        v[i][Ny-1] = 0.0;
    }

    // bottom wall
    for(int i = 0; i < Nx; i++)
    {
        u[i][0] = 0.0;
        v[i][0] = 0.0;
    }

    // block 1: x = id[0] to id[1], y = 0 to id[4]
    for(int i = id[0]; i <= id[1]; i++)
    {
        for(int j = 0; j <= id[4]; j++)
        {
            u[i][j] = 0.0;
            v[i][j] = 0.0;
        }
    }

    // block 2: x = id[2] to id[3], y = 0 to id[4]
    for(int i = id[2]; i <= id[3]; i++)
    {
        for(int j = 0; j <= id[4]; j++)
        {
            u[i][j] = 0.0;
            v[i][j] = 0.0;
        }
    }
}

void ApplyPsiSolidBC(vector<vector<double>>& psi,
                     vector<vector<double>>& X,
                     vector<vector<double>>& Y)
{
    int Nx = X.size();
    int Ny = X[0].size();

    vector<int> id = boundarypoints(X, Y, 5.0, 6.0, 8.0, 9.0, 1.0);

    // bottom wall
    for(int i = 0; i < Nx; i++)
    {
        psi[i][0] = 0.0;
    }

    // block 1
    for(int i = id[0]; i <= id[1]; i++)
    {
        for(int j = 0; j <= id[4]; j++)
        {
            psi[i][j] = 0.0;
        }
    }

    // block 2
    for(int i = id[2]; i <= id[3]; i++)
    {
        for(int j = 0; j <= id[4]; j++)
        {
            psi[i][j] = 0.0;
        }
    }
}

int SIP_BiCGSTAB(vector<vector<double>> &x, vector<vector<double>> &b, SparseMatrix &A, StoneCoeffs &S)
{
double tol =1e-6;
int count =0;
int count_max = 1e4;
int Nx = x.size();
int Ny = x[0].size();

vector<vector<double>> r_true = VectorSubtraction(b,A.MatrixVectorMultiplication(x));
vector<vector<double>> r = S.M_Inverse(r_true);
vector<vector<double>> r_0 = r;
double alpha =0, beta =0,ksi =0;
double res = A.Residual(b,x);
vector <vector<double>> p = r_0;
vector <vector<double>> temp(r.size(),vector<double>(r[0].size(),0.0));
vector<vector<double>> v(r.size(),vector<double>(r[0].size(),0.0));
vector<vector<double>> z(r.size(),vector<double>(r[0].size(),0.0)),t(r.size(),vector<double>(r[0].size(),0.0));


while(res>tol && count<count_max)
{  
    temp  = VectorSubtraction(p,ScalarMultiplication(v,ksi));
    for(int i =0;i<p.size();i++)
    {
        for(int j =0;j<p[0].size();j++)
        {
            p[i][j] = r[i][j] + beta*temp[i][j];
        }
    }
    v = S.M_Inverse(A.MatrixVectorMultiplication(p));
    alpha = (VectorDot(r,r_0))/(VectorDot(v,r_0));

     z = VectorSubtraction(r,ScalarMultiplication(v,alpha));
     t = S.M_Inverse(A.MatrixVectorMultiplication(z));
    
    ksi = (VectorDot(t,z))/VectorDot(t,t);
    
    for(int i = 0; i < Nx-2; i++)
    {
        for(int j = 0; j < Ny-2; j++)
        {
            x[i+1][j+1] = x[i+1][j+1] + alpha*p[i][j] + ksi*z[i][j];
        }
    }
    double rho = (VectorDot(r_0,r));

    r = VectorSubtraction(z,ScalarMultiplication(t,ksi));

    beta = (alpha*(VectorDot(r_0,r)))/(ksi*rho);

    count ++;
    res = A.Residual(b,x);

    for(int j = 1; j < Ny-1; j++)
    {
        x[Nx-1][j] = x[Nx-2][j];
    }

}

if(res<= tol)
{
    return count;
}

else
{
   return -1;
}

}

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

    double delta_t = 1e-3;
    double Re = 100;
    int Nx = xc.size();
    int Ny = yc.size();

    SparseMatrix A_psi(Nx,Ny);
    SparseMatrix A_omega_raw(Nx,Ny);
    SparseMatrix A_omega(Nx,Ny);
    vector <vector<double>> u (Nx,vector<double>(Ny,0.0));
    vector <vector<double>> v (Nx,vector<double>(Ny,0.0));
    vector <vector<double>> psi (Nx, vector<double>(Ny,0.0));
    vector <vector<double>> omega (Nx, vector<double>(Ny,0.0));

    A_psi.setValuesPsi(X,Y);
    A_omega_raw.setValuesOmega(X,Y,delta_t,Re,u,v);
    A_omega.setValuesOmega(X,Y,delta_t,Re,u,v);
    A_psi.setBoundaryConditions(5.0,6.0,8.0,9.0,1.0,X,Y);
    A_omega.setBoundaryConditions(5.0,6.0,8.0,9.0,1.0,X,Y);

    StoneCoeffs S_psi(Nx,Ny,0.9);
    S_psi.SetCoeffs(X,Y,A_psi);


    double u_inf = 1.0;

    vector <vector<double>> b_psi(Nx-2, vector<double>(Ny-2,0.0));
    vector <vector<double>> b_omega(Nx-2, vector<double>(Ny-2,0.0));

    //Boundary values:
    //psi

    for(int j = 0; j < Ny; j++)  //left
    {
        psi[0][j] = u_inf * Y[0][j];
        u[0][j] = u_inf;
    }

    for(int i = 0; i < Nx; i++) 
    {
        psi[i][Ny-1] = psi[0][Ny-1]; //top
        u[i][Ny-1] = u_inf;
    }

    double t = 0.0;
    double t_max = 60;
    vector<vector<double>> omega_0(Nx, vector<double>(Ny,1.0));
    vector<vector<double>> p = VectorSubtraction(omega,omega_0);
    double rem = sqrt((VectorDot(p,p))/(Nx*Ny));
    double crit = 1e-6;
    int step = 0;


ofstream file_log("venv/VTE_data_1.csv");
ofstream file_time("venv/CPU_time.csv");

auto block_start = high_resolution_clock::now();
    while(t<t_max )
    {
    int i, j;

    for(i =0;i<Nx;i++)
    {
        for(j=0;j<Ny;j++)
        {
            omega_0[i][j] = omega[i][j];
        }
    }

    Update_u_v(u,v,X,Y,psi);
    ApplyVelocityBC(u,v,X,Y,u_inf);
    A_omega.setValuesOmega(X,Y,delta_t,Re,u,v);
    A_omega_raw.setValuesOmega(X,Y,delta_t,Re,u,v);
    A_omega.setBoundaryConditions(5.0,6.0,8.0,9.0,1.0,X,Y);
    StoneCoeffs S_omega(Nx,Ny,0.9);
    S_omega.SetCoeffs(X,Y,A_omega);        

   for(int i = 0; i < Nx-2; i++)
    {   
        for(int j = 0; j < Ny-2; j++)
        {   
            b_omega[i][j] = omega_0[i+1][j+1];
        }
    }

    UpdateOmegaAndb_omegaAtBoundaries(omega,b_omega, psi, A_omega_raw, X, Y, u_inf);
    int iter_omega = SIP_BiCGSTAB(omega,b_omega,A_omega,S_omega);

  
        for(int i = 0; i<Nx-2;i++)
    {
        for(int j =0; j<Ny-2;j++)
        {
            int I = i+1;
            int J = j+1;
            if(isInsideBlock(5.0,6.0,1.0,I,J,X,Y) ||isInsideBlock(8.0,9.0,1.0,I,J,X,Y))
            {
                b_psi[i][j] = 0.0;   // solid block psi value
            }
            else
            {
              b_psi[i][j] = (-omega[I][J]/2.0)*(X[I+1][J]-X[I-1][J])*(Y[I][J+1]-Y[I][J-1]);
            }
        }
    }
         for(int l = 0; l<Ny-2;l++)
    {   int L = l+1;
        b_psi[0][l] +=   -(Y[1][L+1]-Y[1][L-1])/(X[1][L] - X[0][L])*psi[0][L];
    }

    for(int m =0; m<Nx-2;m++)
    {   int M = m+1;
        b_psi[m][Ny-3] +=  -(X[M+1][Ny-2]-X[M-1][Ny-2])/(Y[M][Ny-1]-Y[M][Ny-2])*psi[M][Ny-1] ;
    }
    ApplyPsiSolidBC(psi,X,Y);
    int iter_psi = SIP_BiCGSTAB(psi,b_psi,A_psi,S_psi);
    t = t+delta_t;
    p = VectorSubtraction(omega,omega_0);
    rem = sqrt((VectorDot(p,p))/(Nx*Ny));

    step++;
if(step%100 ==0)
{   
  file_log << t << ","<<log10(t)<<","<<iter_omega<<","<<log10(iter_omega)<<","<<iter_psi<<","<<log10(iter_psi)<<","<<rem<<","<<log10(rem)<<"\n";
  std::cout << t << "," << iter_omega << "," << iter_psi << ","  << rem << "\n";
}

if (step % 1000 == 0)
{
    double cpu_time =
        duration<double>(high_resolution_clock::now() - block_start).count();

    file_time << step << "," << t << "," << cpu_time << "\n";
    block_start = high_resolution_clock::now();
}

if(step > 500 && rem < crit)
{
    cout << "Steady state reached at t = " << t << endl;
    break;
}
    if(step%500 ==0)
    {
        WriteFieldCSV("venv/psi_t_1.csv", psi);
        WriteFieldCSV("venv/omega_t_1.csv", omega);

    }
}

file_log.close();
file_time.close();
WriteFieldCSV("venv/psi_final_1.csv", psi);
WriteFieldCSV("venv/omega_final_1.csv", omega);
WriteFieldCSV("venv/u_final_1.csv",u);
WriteFieldCSV("venv/v_final_1.csv",v);

return 0;
}


