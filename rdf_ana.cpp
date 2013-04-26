/*
analyze the rdf of two types of atoms
*/
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <cstdio>
#include <cstdlib>
using namespace std;

typedef struct {
double x; double y; double z;	
int label;       
} atom;

double distance(atom *, atom *, double &, double &, double &);// calculate a-b distance
void doRdf(atom *, atom *, const int &, const int &, double &, double &, double &, double *, const int &, const double &); 
void doReadFile(FILE *, atom *, int ,  double &, double &, double &); 
void doOutput(double *, const int &, double &, const int &); 

double distance(atom *atom1, atom *atom2, double &px, double &py, double &pz)
{
  double x,y,z;
  double r;
  x=abs(atom1->x-atom2->x);
  if (x>px/2.0) x=px-x;
  y=abs(atom1->y-atom2->y);
  if (y>py/2.0) y=py-y;
  z=abs(atom1->z-atom2->z);
  if (z>pz/2.0) z=pz-z;
  r=x*x+y*y+z*z;
  return pow(r,0.5);
}

void doRdf(atom *atom1, atom *atom2, const int &natom1, const int &natom2, 
		   double &px, double &py, double &pz,
		   double rdf[], const int &maxbin,const double &rmax)
{
  double r; // distance
  double dv, s; // volume and area
  int nbin; // array bins
  double *cn1;
  double *cn2;
  double *para1;
  double *para2;
  const double PI=3.1415926536;
  cn1=(double *) malloc(maxbin*sizeof(double));
  para1=(double *) malloc(maxbin*sizeof(double));
  cn2=(double *) malloc(maxbin*sizeof(double));
  para2=(double *) malloc(maxbin*sizeof(double));
  
  double dr;
  dr=rmax/maxbin;
  double dens1,dens2;
  dens1=natom1/(px*py*pz);
  dens2=natom2/(px*py*pz);
  // calculate s/dv
  for (int bin=0;bin<maxbin;bin++)
  {
  cn1[bin]=0.0;
  cn2[bin]=0.0;
  r=(bin+1)*dr;
  para1[bin]=4*PI*pow(r,2)*dr;
  para2[bin]=4*PI*pow(r,2)*dr;
  }
  // the amount
  for (int iatom1=0;iatom1<natom1;iatom1++)
    {
      for (int iatom2=0;iatom2<natom2;iatom2++)
	{ 
	  if ( atom1[iatom1].label != atom2[iatom2].label)
	  {
		r=distance(&atom1[iatom1],&atom2[iatom2],px,py,pz);// a-b distance
		if (r<rmax-0.1)
		{
//	  	cout<<"r: "<<r<<endl;
		nbin=int(r/dr);
//	  cout<<"nbin: "<<nbin<<endl;
		cn2[nbin]++;	  
		}	  	  
	  }
	  
	}
    }
  // number*para
  for (int bin=0;bin<maxbin;bin++)
  {
  rdf[bin]=rdf[bin]+cn2[bin]/para2[bin]/dens2/natom1;
  }
}

void doReadFile(FILE *pfile,  
                atom atom1[], int natom,
				double &px, double &py, double &pz)
{
	// read the same atoms from the file
	int n,n1,n2;
	char char_name[80];
	char char1[5], char2[5];
	double x=0,y=0,z=0;
	fgets (char_name,80,pfile);
	fscanf (pfile, "%d\n", &n); // read the amount of sol molecules
	for (int iatom=0;iatom<natom;iatom++)
	{
// read the coord for OW  	
    fscanf(pfile,"%5d%5s%5s%5d%lf%lf%lf\n", &n1, char1, char2, &n2, &x, &y, &z);
	atom1[iatom].x=x; atom1[iatom].y=y; atom1[iatom].z=z;atom1[iatom].label=n1;
	}
// read the periodic condition
 fscanf (pfile,"%lf%lf%lf\n" , &px, &py, &pz) ; //  
// end of atom file
}

void doOutput(double rdf[], const int &maxbin, double &dr, const int &nframe)
{
FILE *pfile_output;
pfile_output=fopen("rdf.dat","w+");
for (int bin=0;bin<maxbin;bin++)
{
fprintf (pfile_output, "%8.4lf%12.6lf\n", (bin+1)*dr,rdf[bin]/nframe); 
}
}

void Display() // display notation
{
cout<<" This program analyze the rdf of two types of atoms"<<endl;
cout<<" rdf_ana natom1 natom2 rmax maxbin nframe  "<<endl;
cout<<" natom1: amount of atom1"<<endl;
cout<<" natom2: amount of atom2"<<endl;
cout<<" rmax: rmax"<<endl;
cout <<" maxbin: maxbin "<<endl;
cout<<" n_frame: nframe"<<endl;
}

int main(int argc, char *argv[])
{
   Display();
// arguments
	const int natom1=atoi(argv[1]); // amount of atom1
	const int natom2=atoi(argv[2]); // atom2 amount
	const double rmax=atof(argv[3]); // max r
	const int maxbin=atoi(argv[4]); // bin amount 
	const int nframe=atoi(argv[5]); // bin amount 
	cout<<"natom1: "<<natom1<<endl;
	cout<<"natom2: "<<natom2<<endl;
	cout<<"rmax: "<<rmax<<endl;
	cout<<"maxbin: "<<maxbin<<endl;
	cout<<"n_frame: "<<nframe<<endl;
	atom *atom1;
	atom *atom2;
	atom1=(atom *) malloc(natom1*sizeof(atom));
	atom2=(atom *) malloc(natom2*sizeof(atom));
	double *rdf;
	rdf=(double*) malloc(maxbin*sizeof(double));
	for (int bin=0;bin<maxbin;bin++)
	{rdf[bin]=0;}
	
	double dr;
	dr=rmax/maxbin;	
	double px,py,pz;
	px=py=pz=0;
// open files
	FILE *pfile_atom1, *pfile_atom2;
	pfile_atom1=fopen("O.gro","r");
	pfile_atom2=fopen("N.gro","r");
  //read the file
  for (int iframe=0;iframe<nframe;iframe++)
    {
	doReadFile( pfile_atom1, atom1, natom1, px, py, pz);
	doReadFile( pfile_atom2, atom2, natom2, px, py, pz);	
    doRdf(atom1,atom2,natom1,natom2,px,py,pz,rdf,maxbin,rmax); //analyze the rdf
    }
  // output data
  doOutput(rdf,maxbin,dr,nframe);
  return 0;
}
