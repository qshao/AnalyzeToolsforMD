/*
Hydration Dynamics
*/
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <cstdio>
#include <cstdlib>
using namespace std;

typedef struct {
double *x ;  
double *y ;
double *z ;
double CoorNum;  // coordination number
} solute;

typedef struct {
int id;
double xo,yo,zo;
double xh1,yh1,zh1;
double xh2,yh2,zh2;
double xw,yw,zw;
bool vacant;
} water;

double Distance (double &, double &, double &, double &, double &, double &, double &, double &, double &); // calculate the distance
bool IsHydration(solute , water , double &, double &, double &, const double &, const int &); // IsHydration?
void Display(); // Display information
void doAnalysis(solute *, water *, double *, const int &,const int &, const int &, const int &, double &, double &, double &, const double &); // Do analysis per frame
void myReadFile(FILE *, FILE *, solute *, water *, const int &, const int &, const int &, double &, double &, double &); // read file 
void doOutput(double *, double *, int , double , FILE *);
int doList(solute *, water *, int *, const int &, const int &, double &, double &, double &, const double &, const int &);
double doCompare(int *, int *, int , int );


double Distance(double &x1, double &y1, double &z1,
			   double &x2, double &y2, double &z2,
			   double &px, double &py, double &pz) // calculate the distance
{
double r1=0, r1x=0, r1y=0, r1z=0;
//cout<<"px: "<<px<<endl;
r1x=abs(x1-x2);
if (r1x>px/2) r1x=px-r1x;  //pbc correction
r1y=abs(y1-y2);
if (r1y>py/2) r1y=py-r1y;   // pbc correction
r1z=abs(z1-z2);
if (r1z>pz/2) r1z=pz-r1z;  // pbc correction
r1=r1x*r1x+r1y*r1y+r1z*r1z;
r1= sqrt(r1);
return r1;
}			   

int doList(solute zw[], water sol[], int *idStore,
           const int &nmol, const int &nsol,
		   double &px, double &py, double &pz, const double &d_hydra, const int &natom) // create the frame list
{
int n=0;
for (int imol=0;imol<nmol;imol++)
{
	for (int isol=0;isol<nsol;isol++)
	{
		bool Hydr=false;
		Hydr=IsHydration(zw[imol], sol[isol], px, py, pz,d_hydra, natom);
		if (Hydr && sol[isol].vacant)
		{
		idStore[n]=sol[isol].id; n++;
//		 cout<<sol[isol].id<<endl;
		 sol[isol].vacant=false;}
	}
}
return n;
}

double doCompare(int *idStoreNew, int *idStoreOld, int nNew, int nOld) // return remain percentage
{
double nSame = 0.0; 
for (int i=0;i<nOld;i++)
{
	for (int j=0;j<nNew;j++)
	{
	if (idStoreOld[i] == idStoreNew[j])
	nSame ++ ; 
	}
}
double ratio= double(nSame/nOld);
//cout << "nOld: "<< nOld;
//cout<< ratio<<endl;
return ratio;
}

bool IsHydration(solute zw, water sol,
				double &px, double &py, double &pz,
				const double &d, const int &natom) // in hydration layer?
{
	for (int iatom = 0; iatom<natom; iatom++)
	{
	double r=Distance(zw.x[iatom],zw.y[iatom],zw.z[iatom],sol.xo,sol.yo,sol.zo,px,py,pz);
	if (r<=d) 
	{
//	cout<<"r: "<<r<<endl;
	return true;}
	}
	return false;
}



void myReadFile(FILE *pfile_sol, FILE *pfile_zw, 
                solute zw[], water sol[], 
				const int &nmol, const int &nsol, const int &natom, 
				double &px, double &py, double &pz)
{
	// read the sol file
	int n,n1,n2;
	char char_name[80];
	char char1[5], char2[5];
	double x=0,y=0,z=0;
	fgets (char_name,80,pfile_sol);
	fscanf (pfile_sol, "%d\n", &n); // read the amount of sol molecules
	for (int isol=0;isol<nsol;isol++)
	{
// read the coord for OW  	
    fscanf(pfile_sol,"%5d%5s%5s%5d%lf%lf%lf\n", &n1, char1, char2, &n2, &x, &y, &z);
	sol[isol].xo=x; sol[isol].yo=y; sol[isol].zo=z;
	sol[isol].id=n1;
// read the coord for HW1  
	fscanf(pfile_sol,"%5d%5s%5s%5d%lf%lf%lf\n", &n1, char1, char2, &n2, &x, &y, &z);
//	cout<<n1<<" "<<char1<<" "<<char2<<" "<<n2<<" "<<x<<" "<<y<<" "<<z<<endl;
	sol[isol].xh1=x; sol[isol].yh1=y; sol[isol].zh1=z;
// read the coord for HW2  
	fscanf(pfile_sol,"%5d%5s%5s%5d%lf%lf%lf\n", &n1, char1, char2, &n2, &x, &y, &z );
//	cout<<n1<<" "<<char1<<" "<<char2<<" "<<n2<<" "<<x<<" "<<y<<" "<<z<<endl;
	sol[isol].xh2=x; sol[isol].yh2=y; sol[isol].zh2=z;
// read the coord for MW  
	fscanf(pfile_sol,"%5d%5s%5s%5d%lf%lf%lf\n", &n1, char1, char2, &n2, &x, &y, &z );
	sol[isol].xw=x; sol[isol].yw=y; sol[isol].zw=z;
	}
// read the periodic condition
 fscanf (pfile_sol,"%lf%lf%lf\n" , &px, &py, &pz) ; //  
// end of sol file

// read the zw file  -- modified this function
	fgets (char_name,80,pfile_zw);
	fscanf (pfile_zw, "%d\n", &n); // read the amount of sol molecules
  for (int imol=0;imol<nmol;imol++)
  {
// read the coord for O-zw
		for (int iatom=0;iatom<natom;iatom++)
		{
		   fscanf(pfile_zw,"%5d%5s%5s%5d%lf%lf%lf\n", &n1, char1, char2, &n2, &x, &y, &z);
		   //cout<<n1<<" "<<char1<<" "<<char2<<" "<<n2<<" "<<x<<" "<<y<<" "<<z<<endl;
		   zw[imol].x[iatom]=x; zw[imol].y[iatom]=y; zw[imol].z[iatom]=z;
		}
 
  }
// read the periodic condition
 fscanf (pfile_zw,"%lf%lf%lf\n" , &px, &py, &pz) ; // 	
// end of zw file
}

void doOutput(double arr1[], double arr2[], int nsize, double scale, FILE *pfile_output)// Output an array
{
for (int i=0; i<nsize;i++)
{ double d=i*scale;
  arr1[i]=arr1[i]/arr2[i];
  fprintf (pfile_output, "%8.4lf%12.6lf\n", d,arr1[i]); }
}

void Display() // display notation
{
cout<<" This program calculate the lifetime of Hydration Layer"<<endl;
cout<<" HydDyn nmol nsol natom nframe d_hydra "<<endl;
cout<<" nmol: number of solute"<<endl;
cout<<" nsol: number of sol"<<endl;
cout<<" natom: number of atoms in solute"<<endl;
cout <<" nframe: number of frames "<<endl;
cout<<" d_hydra: hydration layer size"<<endl;
}


int main (int argc, char *argv[])
{
	const int num = 20;
	const int tframe = 500;
	Display(); 
// arguments
	const int nmol=atoi(argv[1]); // number of zw
	const int nsol=atoi(argv[2]); // number of water
	const int natom=atoi(argv[3]); // number of O
	const int nframe=atoi(argv[4]);
	const double d_hydra=atof(argv[5]); // distance for hydration layer
	cout<<"nmol: "<<nmol<<endl;
	cout<<"nsol: "<<nsol<<endl;
	cout<<"natom: "<<natom<<endl;
	cout<<"nframe: "<<nframe<<endl;
	cout<<"d_hydra: "<<d_hydra<<endl;
//  memory allocation
	solute *zw;
	zw=(solute *) malloc(nmol*sizeof(solute));
// array initial
	int *n= new int[nframe]; // number of water per frame
	int **idStore = new int* [nframe]; // id of water in layer
	for (int i=0;i<nframe;i++)
	idStore[i]=new int[num];
	
	double C[nframe];
	double total[nframe];
	
	for (int imol=0;imol<nmol;imol++)
	{
	zw[imol].x=(double *) malloc(natom* sizeof(double));
	zw[imol].y=(double *) malloc(natom* sizeof(double));
	zw[imol].z=(double *) malloc(natom* sizeof(double));
	}
	water *sol;
	sol = (water *) malloc(nsol*sizeof(water));
	double px,py,pz; // pbc size	
// open file pointers	
	FILE *pfile_zw, *pfile_sol, *pfile_lifetime;
    pfile_sol=fopen("sol.gro","r"); // water gro file
    pfile_zw=fopen("O.gro","r"); // O-zw gro file
	pfile_lifetime=fopen("Dynamic_zw.dat","w+");// Distribution of Coord Num

for (int iframe=0;iframe<nframe; iframe++)	
{
// initialize
	double px=0.0;
	double py=0.0;
	double pz=0.0;
	for (int isol=0;isol<nsol;isol++)
	{
	sol[isol].xo=0;sol[isol].yo=0;sol[isol].zo=0;
	sol[isol].xh1=0;sol[isol].yh1=0;sol[isol].zh1=0;
	sol[isol].xh2=0;sol[isol].yh2=0;sol[isol].zh2=0;
	sol[isol].xw=0;sol[isol].yw=0;sol[isol].zw=0;
	sol[isol].vacant=true;
	}
	for (int imol=0;imol<nmol;imol++)
	{
	for (int iatom=0;iatom<natom;iatom++)
	{ zw[imol].x[iatom]=0;zw[imol].y[iatom]=0;zw[imol].z[iatom]=0;}
	zw[imol].CoorNum=0;
	}
	myReadFile(pfile_sol, pfile_zw, zw, sol, nmol, nsol, natom, px,py,pz);
	n[iframe]=doList(zw,sol,idStore[iframe],nmol, nsol,px,py,pz,d_hydra, natom); // create the list	
//	cout<<n[iframe]<<endl;

}
for (int iframe=0;iframe<nframe;iframe++)
{
	int nframe2=iframe+tframe;
	if (nframe2>nframe) 
	nframe2=nframe;	
	for (int iframe2=iframe;iframe2<nframe2;iframe2++)
	{
	int dframe=iframe2-iframe;
	if (n[iframe]>0)
	{
	double dc = doCompare(idStore[iframe2], idStore[iframe], n[iframe2], n[iframe]);
	C[dframe]=C[dframe]+ dc;
//	cout<<"dc: "<<dc<<endl;
	total[dframe]++;
	}
	}
}
// output
	const double delt=1;
	doOutput(C,total,tframe,delt,pfile_lifetime);
// free the memory
 return 0;
}
