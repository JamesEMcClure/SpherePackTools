// Created by James McClure
// Copyright 2008-2011

// Read in sphere pack (pack.out) and write out a PM file (zeros and ones)
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>

using namespace std;

int main (int argc, char * const argv[]) {
	// Use of non-constant grid spacing will give ellipsoid packings!
	cout << "Digitizing Media" << endl;
	// number of nodes, in x, y, z direction
	int dimX,dimY,dimZ;
	// startX,startY,startZ - start location inside of pack
	// lenX,lenY,lenZ - real domain length
	double startX,startY,startZ,lenX,lenY,lenZ; 
	string Pack,PMname;
	
	// Read input variables
	ifstream infile ("Digitize.in");
	infile >> Pack;
	cout << "Packing file is: " << Pack << endl;
	infile >> PMname;
	cout << "Digitized media file is: " << PMname << endl;
	infile >> dimX;
	infile >> dimY;
	infile >> dimZ;
	cout << "Mesh size is : " << dimX << "x" << dimY << "x"<< dimZ << endl;
	infile >> startX;
	infile >> startY;
	infile >> startZ;
	infile >> lenX;
	infile >> lenY;
	infile >> lenZ;
	
	double dx =  lenX / (dimX-1);
	double dy =  lenY / (dimY-1);
	double dz =  lenZ / (dimZ-1);
	
	// Read in the full sphere pack
	int count,nspheres,nBC,i,j,k,p,q,n;  
	double x,y,z,r,num,hx,hy,hz;
	double maxX,maxY,maxZ;
	maxX = startX+lenX;
	maxY = startY+lenY;
	maxZ = startZ+lenZ;
	
	// declare dyanmaic arrays to store sphere info
	vector<double> cx;
	vector<double> cy;
	vector<double> cz;
	vector<double> rad;
	// spheres that cross a boundary stored here (centroids outside of the domain)
	vector<double> cxBC;
	vector<double> cyBC;
	vector<double> czBC;
	vector<double> radBC;
	
	//...... READ IN THE SPHERES...................................
	char * trsh;
	trsh = new char[100];
	cout << "Reading the packing file..." << endl;
	ifstream pack (Pack.c_str());
	//.........Trash the header lines..........
	pack.getline(trsh, 100);
	pack.getline(trsh, 100);
	pack.getline(trsh, 100);
	pack.getline(trsh, 100);
	pack.getline(trsh, 100);
	//........read the spheres..................
	nspheres = 0;
	pack >> x;
	pack >> y;
	pack >> z;
	pack >> r;
	while (! pack.eof()){
		cx.push_back(x);
		cy.push_back(y);
		cz.push_back(z);
		rad.push_back(r);
		pack >> x;
		pack >> y;
		pack >> z;
		pack >> r;
		nspheres++;
	//	cout << "successfully read sphere " << nspheres << endl;
	}
	pack.close();
	cout << "Number of spheres extracted is: " << nspheres << endl;
	// .............................................................
	// ..... DETERMINE SPHERES THAT INTERSECT A BOUNDARY ...........
	nBC = 0;
	for (i=0;i<nspheres;i++){
		// Upper Boundary -faces
		if (cx[i]+rad[i] > startX+lenX) {
			cxBC.push_back(cx[i]-lenX);
			cyBC.push_back(cy[i]);
			czBC.push_back(cz[i]);
			radBC.push_back(rad[i]);
			nBC++;
		}
		if (cy[i]+rad[i] > startY+lenY) {
			cxBC.push_back(cx[i]);
			cyBC.push_back(cy[i]-lenY);
			czBC.push_back(cz[i]);
			radBC.push_back(rad[i]);
			nBC++;
		}
		if (cz[i]+rad[i] > startZ+lenZ) {
			cxBC.push_back(cx[i]);
			cyBC.push_back(cy[i]);
			czBC.push_back(cz[i]-lenZ);
			radBC.push_back(rad[i]);
			nBC++;
		}
		// Lower boundary -faces
		if (cx[i]-rad[i] < startX) {
			cxBC.push_back(cx[i]+lenX);
			cyBC.push_back(cy[i]);
			czBC.push_back(cz[i]);
			radBC.push_back(rad[i]);
			nBC++;
		}
		if (cy[i]-rad[i] < startY) {
			cxBC.push_back(cx[i]);
			cyBC.push_back(cy[i]+lenY);
			czBC.push_back(cz[i]);
			radBC.push_back(rad[i]);
			nBC++;
		}
		if (cz[i]-rad[i] < startZ) {
			cxBC.push_back(cx[i]);
			cyBC.push_back(cy[i]);
			czBC.push_back(cz[i]+lenZ);
			radBC.push_back(rad[i]);
			nBC++;
		}
		// Edges
		if (cx[i]-rad[i] < startX && cy[i]+rad[i] > startY+lenY) {
			cxBC.push_back(cx[i]+lenX);
			cyBC.push_back(cy[i]-lenY);
			czBC.push_back(cz[i]);
			radBC.push_back(rad[i]);
			nBC++;
		}
		if (cx[i]-rad[i] < startX && cz[i]+rad[i] > startZ+lenZ) {
			cxBC.push_back(cx[i]+lenX);
			cyBC.push_back(cy[i]);
			czBC.push_back(cz[i]-lenZ);
			radBC.push_back(rad[i]);
			nBC++;
		}
		if (cx[i]-rad[i] < startX && cy[i]-rad[i] < startY) {
			cxBC.push_back(cx[i]+lenX);
			cyBC.push_back(cy[i]+lenY);
			czBC.push_back(cz[i]);
			radBC.push_back(rad[i]);
			nBC++;
		}
		if (cx[i]-rad[i] < startX && cz[i]-rad[i] < startZ) {
			cxBC.push_back(cx[i]+lenX);
			cyBC.push_back(cy[i]);
			czBC.push_back(cz[i]+lenZ);
			radBC.push_back(rad[i]);
			nBC++;
		}
		// 
		if (cx[i]+rad[i] > startX+lenX && cy[i]+rad[i] > startY+lenY) {
			cxBC.push_back(cx[i]-lenX);
			cyBC.push_back(cy[i]-lenY);
			czBC.push_back(cz[i]);
			radBC.push_back(rad[i]);
			nBC++;
		}
		if (cx[i]+rad[i] > startX+lenX && cz[i]+rad[i] > startZ+lenZ) {
			cxBC.push_back(cx[i]-lenX);
			cyBC.push_back(cy[i]);
			czBC.push_back(cz[i]-lenZ);
			radBC.push_back(rad[i]);
			nBC++;
		}
		if (cx[i]+rad[i] > startX+lenX && cy[i]-rad[i] < startY) {
			cxBC.push_back(cx[i]-lenX);
			cyBC.push_back(cy[i]+lenY);
			czBC.push_back(cz[i]);
			radBC.push_back(rad[i]);
			nBC++;
		}
		if (cx[i]+rad[i] > startX+lenX && cz[i]-rad[i] < startZ) {
			cxBC.push_back(cx[i]-lenX);
			cyBC.push_back(cy[i]);
			czBC.push_back(cz[i]+lenZ);
			radBC.push_back(rad[i]);
			nBC++;
		}
		//
		if (cy[i]-rad[i] < startY && cz[i]-rad[i] < startZ) {
			cxBC.push_back(cx[i]);
			cyBC.push_back(cy[i]+lenY);
			czBC.push_back(cz[i]+lenZ);
			radBC.push_back(rad[i]);
			nBC++;
		}
		if (cy[i]+rad[i] > startY+lenY && cz[i]-rad[i] < startZ) {
			cxBC.push_back(cx[i]);
			cyBC.push_back(cy[i]-lenY);
			czBC.push_back(cz[i]+lenZ);
			radBC.push_back(rad[i]);
			nBC++;
		}
		if (cy[i]-rad[i] < startY && cz[i]+rad[i] > startZ+lenZ) {
			cxBC.push_back(cx[i]);
			cyBC.push_back(cy[i]+lenY);
			czBC.push_back(cz[i]-lenZ);
			radBC.push_back(rad[i]);
			nBC++;
		}
		if (cy[i]+rad[i] > startY+lenY && cz[i]+rad[i] > startZ+lenZ) {
			cxBC.push_back(cx[i]);
			cyBC.push_back(cy[i]-lenY);
			czBC.push_back(cz[i]-lenZ);
			radBC.push_back(rad[i]);
			nBC++;
		}
		// Corners
		if (cx[i]+rad[i]> startX+lenX && cy[i]-rad[i] < startY && cz[i]-rad[i] < startZ) {
			cxBC.push_back(cx[i]-lenX);
			cyBC.push_back(cy[i]+lenY);
			czBC.push_back(cz[i]+lenZ);
			radBC.push_back(rad[i]);
			nBC++;
		}
		if (cx[i]+rad[i]> startX+lenX && cy[i]+rad[i] > startY+lenY && cz[i]-rad[i] < startZ) {
			cxBC.push_back(cx[i]-lenX);
			cyBC.push_back(cy[i]-lenY);
			czBC.push_back(cz[i]+lenZ);
			radBC.push_back(rad[i]);
			nBC++;
		}
		if (cx[i]+rad[i]> startX+lenX && cy[i]-rad[i] < startY && cz[i]+rad[i] > startZ+lenZ) {
			cxBC.push_back(cx[i]-lenX);
			cyBC.push_back(cy[i]+lenY);
			czBC.push_back(cz[i]-lenZ);
			radBC.push_back(rad[i]);
			nBC++;
		}
		if ( cx[i]+rad[i]> startX+lenX && cy[i]+rad[i] > startY+lenY && cz[i]+rad[i] > startZ+lenZ) {
			cxBC.push_back(cx[i]-lenX);
			cyBC.push_back(cy[i]-lenY);
			czBC.push_back(cz[i]-lenZ);
			radBC.push_back(rad[i]);
			nBC++;
		}
		//
		if (cx[i]-rad[i] < startX && cy[i]-rad[i] < startY && cz[i]-rad[i] < startZ) {
			cxBC.push_back(cx[i]+lenX);
			cyBC.push_back(cy[i]+lenY);
			czBC.push_back(cz[i]+lenZ);
			radBC.push_back(rad[i]);
			nBC++;
		}
		if (cx[i]-rad[i] < startX && cy[i]+rad[i] > startY+lenY && cz[i]-rad[i] < startZ) {
			cxBC.push_back(cx[i]+lenX);
			cyBC.push_back(cy[i]-lenY);
			czBC.push_back(cz[i]+lenZ);
			radBC.push_back(rad[i]);
			nBC++;
		}
		if (cx[i]-rad[i] < startX && cy[i]-rad[i] < startY && cz[i]+rad[i] > startZ+lenZ) {
			cxBC.push_back(cx[i]+lenX);
			cyBC.push_back(cy[i]+lenY);
			czBC.push_back(cz[i]-lenZ);
			radBC.push_back(rad[i]);
			nBC++;
		}
		if ( cx[i]-rad[i] < startX && cy[i]+rad[i] > startY+lenY && cz[i]+rad[i] > startZ+lenZ) {
			cxBC.push_back(cx[i]+lenX);
			cyBC.push_back(cy[i]-lenY);
			czBC.push_back(cz[i]-lenZ);
			radBC.push_back(rad[i]);
			nBC++;
		}
		// End of Assigning the boundary spheres
	}
	cout << "Number of boundary spheres: " << nBC << endl;
	// .............................................................	
	// Use sphere lists to determine which nodes are in porespace
	// Write out binary file for nodes
	char value;
	int size = dimX*dimY*dimZ;
	int nodno;
	int nodecount=0;
	hx = lenX/(dimX-1);
	hy = lenY/(dimY-1);
	hz = lenZ/(dimZ-1);	
	vector<char> nodes(size);
	int imin,imax,jmin,jmax,kmin,kmax;
	char * indicator;
	indicator = new char [size];
	// initialize indicator to one
	for (i=0;i<size;i++){
		indicator[i]=1;
	}
	// ............ REGULAR SPHERES .................
	for (p=0;p<nspheres;p++){
		imin = int ((cx[p]-rad[p])/hx)-1;
		imax = int ((cx[p]+rad[p])/hx)+1;
		jmin = int ((cy[p]-rad[p])/hy)-1;
		jmax = int ((cy[p]+rad[p])/hy)+1;
		kmin = int ((cz[p]-rad[p])/hz)-1;
		kmax = int ((cz[p]+rad[p])/hz)+1;
		if (imin<0){ imin = 0;}
		if (imax>dimX){imax = dimX;}
		if (jmin<0){ jmin = 0;}
		if (jmax>dimY){jmax = dimY;}
		if (kmin<0){ kmin = 0;}
		if (kmax>dimZ){kmax = dimZ;}
		for (i=imin;i<imax;i++){
			for (j=jmin;j<jmax;j++){
				for (k=kmin;k<kmax;k++){
					// initialize to porespace value of one
					x = hx*i;
					y = hy*j;
					z = hz*k;	
					value = 1;
					// if inside sphere, set to zero
					if ( (cx[p]-x)*(cx[p]-x)+(cy[p]-y)*(cy[p]-y)
						+(cz[p]-z)*(cz[p]-z) < rad[p]*rad[p]){
						value=0;
					}
					
					// get the position in the list
					nodno = i+j*dimX+k*dimX*dimY;
					if ( indicator[nodno] != 0 ){
						indicator[nodno] = value;
					}
					
				}
			}
		}
	}
	// ............ BOUNDARY SPHERES .................
	for (p=0;p<nBC;p++){
		imin = int ((cxBC[p]-radBC[p])/hx)-1;
		imax = int ((cxBC[p]+radBC[p])/hx)+1;
		jmin = int ((cyBC[p]-radBC[p])/hy)-1;
		jmax = int ((cyBC[p]+radBC[p])/hy)+1;
		kmin = int ((czBC[p]-radBC[p])/hz)-1;
		kmax = int ((czBC[p]+radBC[p])/hz)+1;
		if (imin<0){ imin = 0;}
		if (imax>dimX){imax = dimX;}
		if (jmin<0){ jmin = 0;}
		if (jmax>dimY){jmax = dimY;}
		if (kmin<0){ kmin = 0;}
		if (kmax>dimZ){kmax = dimZ;}
		for (i=imin;i<imax;i++){
			for (j=jmin;j<jmax;j++){
				for (k=kmin;k<kmax;k++){
					// initialize to porespace value of one
					x = hx*i;
					y = hy*j;
					z = hz*k;	
					value = 1;
					// if inside sphere, set to zero
					if ( (cxBC[p]-x)*(cxBC[p]-x)+(cyBC[p]-y)*(cyBC[p]-y)
						+(czBC[p]-z)*(czBC[p]-z) < radBC[p]*radBC[p]){
						value=0;
					}
					
					// get the position in the list
					nodno = i+j*dimX+k*dimX*dimY;
					if ( indicator[nodno] != 0 ){
						indicator[nodno] = value;
					}
					
				}
			}
		}
	}
	// ............... INDICATOR FUNCTION HAS BEEN DETERMINED ..............
	
	// Write output file
	char writeval;
	ofstream PM(PMname.c_str(),ios::binary);
	for (i=0;i<size;i++){
		writeval = indicator[i];
		if (writeval == 1){
			nodecount++;
		}
		PM.write((char *) (&writeval), sizeof(writeval));
	}
	PM.close();
		
	double porosity = double(nodecount) / double(size);
	
	cout <<  "Porosity: " << porosity << endl;
	
	// write sphere packing
	Pack.append("BC");
	ofstream cpack (Pack.c_str());
	// Write boundary of spheres
	for (p=0;p<nBC;p++){
		cpack << cxBC[p] << " " << cyBC[p] << " " << czBC[p] << " " << radBC[p] << endl;
	}
	cpack.close();
	
	return 0;
}