/*
Hydration Structure
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
double xo,yo,zo;
double xh1,yh1,zh1;
double xh2,yh2,zh2;
double xw,yw,zw;
bool vacant;
} water;

double Distance (double &, double &, double &, double &, double &, double &, double &, double &, double &); // calculate the distance
double Orient(double &, double &, double &, double &, double &, double &, water & ); // Dipole Orientation
bool IsHydration(double , double , double , water *, const double &); // IsHydration?
void Display(); // Display information
void doAnalysis(solute *, water *, double *, const int &,const int &, const int &, const int &, double &, double &, double &, const double &); // Do analysis per frame
void myReadFile(FILE *, FILE *, solute *, water *, const int &, const int &, const int &, double &, double &, double &); // read file 
void doOutput(double *, int , double , FILE *);




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

double Orient(double &x, double &y, double &z,
			double &px, double &py, double &pz,
			water &sol) // Dipole angle
{
double theta=0.0;
double xh=0.0,yh=0.0,zh=0.0;
//cout<<"x: "<<sol.xo<<" y: "<<sol.yo<<" zo: "<<sol.zo<<endl;
xh=(sol.xh1+sol.xh2)/2.0;
yh=(sol.yh1+sol.yh2)/2.0;
zh=(sol.zh1+sol.zh2)/2.0;
double r1=0,r2=0,r3=0;
r1=Distance(x,y,z,sol.xo,sol.yo,sol.zo,px,py,pz);
r2=Distance(sol.xo,sol.yo,sol.zo,xh,yh,zh,px,py,pz);
r3=Distance(x,y,z,xh,yh,zh,px,py,pz);
//cout<<"r1 "<<r1<<endl;
//cout<<"r2 "<<r2<<endl;
//cout<<"r3 "<<r3<<endl;
theta=(r1*r1+r2*r2-r3*r3)/2.0/r1/r2;
//cout<<r1<<" "<<r2<<" "<<r3<<endl;
//cout<<"theta: "<<theta<<endl;
return theta; 
} 

bool IsHydration(double x, double y, double z, 
				water &sol,
				double &px, double &py, double &pz,
				const double &d) // in hydration layer?
{
 double r=Distance(x,y,z,sol.xo,sol.yo,sol.zo,px,py,pz);
 if (r<=d && r>0.15)
 {  return true; }
 else 
 {return false;}
}

void doAnalysis(solute zw[], water sol[], 
				double Dipole[], const int &nDipole,      // Dipole Distribution 
				const int &nmol, const int &nsol, const int &natom,   
				double &px, double &py, double &pz, const double &d_hydra) // calculate the hydration structure per frame
{
//cout<<"px: "<<px<<endl;
	double Angle=0.0;
	const double dAngle= 2.0/nDipole;
//	cout<<"nDipole: "<<nDipole<<endl;
	for (int imol=0;imol<nmol;imol++)
	{
	for (int iatom=0;iatom<natom;iatom++)
	{
	for (int isol=0;isol<nsol;isol++)
	{
		bool Hydra=false;
		Hydra=IsHydration(zw[imol].x[iatom],zw[imol].y[iatom],zw[imol].z[iatom],
				 sol[isol], px, py, pz, d_hydra);
		if (Hydra && sol[isol].vacant)  
		{
			// set the Coordination
			zw[imol].CoorNum++;
			sol[isol].vacant=false; }
		if (Hydra)
		{
			// calculate the orient
		Angle=Orient(zw[imol].x[iatom], zw[imol].y[iatom],zw[imol].z[iatom], 
							px,py,pz, 
						sol[isol]);
//			cout<<"Angle: "<<Angle<<endl;
			if (Angle<1 && Angle>-1)
			{int binAngle=0;
			binAngle=int((Angle+1.0)/dAngle);
//			cout<<"binAngle: "<<binAngle<<endl;
			Dipole[binAngle]++;}
		}
	}
	}
	}
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

void doOutput(double arr[], int nsize, double scale, FILE *pfile_output)// Output an array
{
  double sum=0.0;
  for (int i=0;i<nsize;i++)
  sum=sum+arr[i];
for (int i=0; i<nsize;i++)
{ double d=i*scale;
  fprintf (pfile_output, "%8.4lf%12.6lf\n", d,arr[i]/sum); }
}

void Display() // display notation
{
cout<<" This program analyze the hydration number and orientation distribution"<<endl;
cout<<" HydStruct nmol nsol natom nframe d_hydra "<<endl;
cout<<" nmol: number of solute"<<endl;
cout<<" nsol: number of sol"<<endl;
cout<<" natom: number of atoms in solute"<<endl;
cout <<" nframe: number of frames "<<endl;
cout<<" d_hydra: hydration layer size"<<endl;
}


int main (int argc, char *argv[])
{
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
	for (int imol=0;imol<nmol;imol++)
	{
	zw[imol].x=(double *) malloc(natom* sizeof(double));
	zw[imol].y=(double *) malloc(natom* sizeof(double));
	zw[imol].z=(double *) malloc(natom* sizeof(double));
	}
	water *sol;
	sol = (water *) malloc(nsol*sizeof(water));
	const int nCoorDis = 10;
	const int nDipole = 100;
	double CoorDis[nCoorDis] = {0}; // Distribution of CoordNumer from 1-10
	double Dipole[nDipole] = {0}; // Dipole Distribution from -1 to 1, 100 bins
	double px,py,pz; // pbc size	
// open file pointers	
	FILE *pfile_zw, *pfile_sol, *pfile_CoorNum, *pfile_Dis;
    pfile_sol=fopen("sol.gro","r"); // water gro file
    pfile_zw=fopen("O.gro","r"); // O-zw gro file
	pfile_CoorNum=fopen("CoordNum_zw.dat","w+");// Distribution of Coord Num
	pfile_Dis=fopen("Orient_Dist_zw.dat","w+"); // Orient Distribution 

for (int iframe=0;iframe<nframe; iframe++)	
{
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
//	cout<<" End of "<<iframe<<" frame reading"<<endl;
//	cout<<"zw-O1: "<<sol[0].xo<<" "<<sol[0].yo<<" "<<sol[0].zo<<endl;
// 	cout<<" the pbc size is "<<px<<" "<<py<<" "<<pz<<endl;
	// initial the variables

// main loop
	doAnalysis(zw, sol, Dipole, nDipole,
			nmol, nsol, natom, px, py, pz, d_hydra);
	
	for (int imol=0; imol<nmol;imol++)
		{
//		cout<<" The CoorNum is: "<<zw[imol].CoorNum<<endl;
		int n=int(zw[imol].CoorNum);
		CoorDis[n]++;
		}
}
//	cout<<Dipole[50]<<endl;
	doOutput(CoorDis, nCoorDis,1.0,pfile_CoorNum);
	doOutput(Dipole, nDipole,2.0/nDipole,pfile_Dis);
 return 0;
}
