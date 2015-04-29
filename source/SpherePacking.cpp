// Created by James McClure
// Copyright 2008-2011
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>

//#define DEBUG

using namespace std;

double PI = 3.141592653589793;
//******************************************************
// STRUCTURES TO STORE SPHERES
//******************************************************

// ........................................................
//   CELL STORAGE STRUCTURE FOR PARTICLES
//		- cell (icx,icy,icz) :
//			Lower Boundary = icx*lencx,icy*lency,icz*lencz
//			Upper Boundary = (icx+1)*lencx,(icy+1)*lency,(icz+1)*lencz
//		- particle is inside of cell if cell contains centroid
//		- particle data stored externally (centroid, geometry)
// ........................................................
struct CellStorage {
	CellStorage(int &nx, int &ny, int &nz, double &lenx, double &leny, double &lenz);
	~CellStorage();
	// Set up cells to store all spheres
	
	int ncells,ncx,ncy,ncz;		// # cells
	int maxcell;				// max # particles per cell
	double lencx,lency,lencz;	// cell lengths
	
	int * DATA;					// store particles
	int * COUNT;				// actual # particles per cell

	int & CellCount(int icx, int icy, int icz) {
		return COUNT[icz*ncx*ncy+icy*ncx+icx]; 
	}
	
	int & CellEntry(int icx, int icy, int icz, int index) {
		i = icz*ncx*ncy+icy*ncx+icx;
		return DATA[maxcell*i+index];
	}
	
	void Reset();
	
private:
	int i;
};
//........ constructor .....................................
CellStorage::CellStorage(int &nx, int &ny, int &nz, double &lenx, double &leny, double &lenz)
{
	ncx = nx;					// # cells in x direction
	ncy = ny;					// # cells in y direction
	ncz = nz;					// # cells in z direction
	ncells = ncx*ncy*ncz;		// total # cells
	
	lencx = lenx / ncx;			// cell length in x direction
	lency = leny / ncy;			// cell length in y direction
	lencz = lenz / ncz;			// cell length in z direction
	
	maxcell = 10000;
	
	COUNT = new int [ncells];	
	DATA = new int [maxcell*ncells];	
}
//........ destructor .....................................
CellStorage::~CellStorage()
{
	delete COUNT;
	delete DATA;
}

// ...... reset the cell count ............................
void CellStorage::Reset()
{
	for (i=0;i<ncells;i++) COUNT[i] = 0;
}

// ........................................................
//	STORAGE FOR A COLLECTION OF SPHERES
//		- store particle centroid, radius	(cx,cy,cz,Radius)
//		- store the distance to move each particle (dx,dy,dz)
// ........................................................
struct SphereCollection {
	SphereCollection(int &n);
	~SphereCollection();
	int N;					// number of particles
	//.............................................
	// | centroid | radius | 
	double * data;		
	//.......... Access ..................
	// center of mass / centroid for particle i
	double & cx(int i) { return data[7*i];}
	double & cy(int i) { return data[7*i+1];}
	double & cz(int i){ return data[7*i+2];}
	// distance to move particle i
	double & dx(int i) {return data[7*i+3];}
	double & dy(int i){ return data[7*i+4];}
	double & dz(int i){ return data[7*i+5];}
	// cylinder length
	double & Radius(int i){ return data[7*i+6];}	
//	double & Length(int &i){ return data[8*i+7];}
	
	// ..... Initialize a system of spheres 
	void Initialize(double &Lx, double &Ly, double &Lz, double &porosity, 
					double &mu, double &sig);
	
//	int GetContacts(int & i);
//	int & Contact(int &i) {	return work[i];}
private:
	int i;
};

//........ constructor..........................
SphereCollection::SphereCollection(int &n)
{
	N = n;
	data = new double [7*N];
}

//........ destructor ...........................
SphereCollection::~SphereCollection()
{
	delete data;
}

// ...... initialization .........................
void SphereCollection::Initialize(double &Lx, double &Ly, double &Lz, double &porosity, 
				double &mu, double &sig)
{
	// double mu;
	double cxi,cyi,czi,r;
	srand((unsigned)time(0));
	// Expected volume for the 
	//  r = CUBE_ROOT( 3*VOLUME/4/PI) 
	mu = 0.333333333333*log(0.75*(1.0-porosity)*Lx*Ly*Lz/PI/N) - 1.5*sig;

	// r = exp( 0.333333333333333*log((1.0-porosity)*0.75*Lx*Ly*Lz/N/PI));
	for (i=0;i<N;i++){
//		cout << "INITIALIZE SPHERE " << 1.0*rand()/RAND_MAX << endl;
		// Uniformly distribute the centroids
		cxi = rand()*Lx/RAND_MAX;
		cyi = rand()*Ly/RAND_MAX;
		czi = rand()*Lz/RAND_MAX;
		cx(i) = cxi;
		cy(i) = cyi;
		cz(i) = czi;
		// Generate radii from lognormal distribution, truncate lower & upper  1.5*sig
		r = exp(sqrt(-2*log(1.0*rand()/RAND_MAX))*cos(2.0*PI*rand()/RAND_MAX)*sqrt(sig)+mu);
		while (r < exp(mu)/pow(exp(1.5*sqrt(sig)),2) || r > exp(mu)*pow(exp(1.5*sqrt(sig)),2) ){
			r = exp(sqrt(-2*log(1.0*rand()/RAND_MAX))*cos(2.0*PI*rand()/RAND_MAX)*sqrt(sig)+mu);
		}
		Radius(i) = r;
//		cout << "	Radius value: " << r << endl; 
	}
}

// *********** FUNCTION DECLARATIONS *******************
inline void Random_Displacement(SphereCollection &Particles, double size);
inline bool CHECK_OVERLAPS(SphereCollection &Particles, CellStorage &Storage, double tol);
inline double CoordinationNumber(SphereCollection &Particles, double size,
								 double Lx, double Ly, double Lz);

//******************************************************
//					MAIN CODE 
//******************************************************
int main (int argc, char * const argv[]) {

	// ........ Input variables ..............
	int N,iteration_cutoff;
	int ncx,ncy,ncz;
	double porosity_target, porosity_initial,Lx,Ly,Lz;
	double M,Mt,S,factor,tol;
//	double mean_radius, mean_length;
//	double stddev_radius, stddev_length;
//	double shift_factor;
	//........................................	
	
	// ..................... READ INPUT VARIABLES ...........,...............
	ifstream input("pack.in");
	input >> N;
	cout << "Number of particles is: " << N << endl;
	//input >> M;
	input >> S;
	//cout << "standard deviation for the radii is: " << stddev << endl;
	input >> porosity_initial;
	cout << "Initial porosity is: " << porosity_initial<< endl;
	input >> porosity_target;
	cout << "Target porosity is: " << porosity_target<< endl;
	input >> Lx >> Ly >> Lz;
	cout << "Domain size is: "<< Lx << "," << Ly << "," << Lz << endl;
	input >> ncx >> ncy >> ncz;
	cout << "Number of cells is: "<< ncx << "," << ncy << "," << ncz << endl;
	input >> iteration_cutoff;
	cout << "Maximum number of iterations is: " << iteration_cutoff << endl;
	input >> factor;
	cout << "Radius scaling factor: " << factor << endl;
	input >> tol;
	cout << "Error tolerance: " << tol << endl;	
	input.close();
	// ......................................................................
	
	// internal variables
	bool time_to_stop; //, random_shift;
	int icx,icy,icz,ix,iy,iz,inx,iny,inz;
	int i,j,ii,jj,count,iter,overlaps;
	int index,mainCellCount,failCount;
	double ax,ay,az,a,kx,ky,kz,k,delta,d;
	double porosity,V,Mprev,Sprev;
	double bcx,bcy,bcz;	
	double dxmin,dymin,dzmin;
	double err,max_overlap;

	SphereCollection Particles(N),Save(N);
	CellStorage Storage(ncx,ncy,ncz,Lx,Ly,Lz);
	// Generate a system of particles so that log(r) ~ N(M,S)
	Particles.Initialize(Lx,Ly,Lz,porosity_initial,M,S);
	
	
	// value of mu which corresponds with the target porosity (approx)
	Mt = 0.333333333333*log(0.75*(1.0-porosity_target)*Lx*Ly*Lz/PI/N) - 1.5*S;


	// ............. CHECK INITIAL POROSITY .............................
	V = 0;
	for (i=0;i<Particles.N;i++)	V += 4*PI*pow(Particles.Radius(i),3)/3;
	porosity = 1.0 - V/Lx/Ly/Lz;
	cout << "Initial porosity (actual) " << porosity << endl;
	cout << "Initial value for mu: " << M << endl;
	cout << "Initial value for sigma (input): " << sqrt(S) << endl;
	// ..................................................................
	
	double input_factor = factor;
	factor = 1.025;
	
	//cin  >> time_to_stop;
	// min_radius = max_radius = mean_radius;
	cout << "BEGIN ITERATIONS "<< endl;
	time_to_stop = false;
	failCount = 0;
	while (time_to_stop == false){
		// ..................................................................
		// .............. INCREASE SIZE OF ALL RADII ........................

		cout << "	Increasing size of radii... " << endl;
		V = 0.0;
		// Pre-compute the new porosity:
		for (i=0;i<Particles.N;i++){
			V += 4*PI*pow(factor*Particles.Radius(i),3)/3;
		}
		porosity = 1.0-V/Lx/Ly/Lz;
		if (porosity < porosity_target){
			cout << "********************************************" << endl;
			cout << "Preparing for final iteration... " << endl;
			time_to_stop = true;
			// current sphere volume
			V = 0.0;
			for (i=0;i<Particles.N;i++){
				V += 4*PI*pow(Particles.Radius(i),3)/3;
			}
			
			factor = pow(Lx*Ly*Lz*(1-porosity_target)/V,1.0/3.0);
			cout << "	Radii will be rescaled by: " << factor << endl;
		}
		
		
		for (i=0;i<Particles.N;i++) Particles.Radius(i) = factor*Particles.Radius(i);
		// update the mean & variance 
		Mprev = M;
		Sprev = S;
		M = Mprev+log(factor);
		cout << "	log(r) now normally distributed with mean " << M 
				<< " and variance " << S << endl;
		
		// ..................................................................
		// ... ITERATE UNTIL OVERLAPS ARE ELIMINATED OR CUTOFF IS REACHED ...
		iter = 0;
		err = 1.0;
		while (iter < iteration_cutoff && err > tol ){
			err = 0.0;
			iter++;
			overlaps = 0;
			max_overlap = 0.0;
			// ..............................................................
			// ...... ASSIGN EACH PARTICLE TO THE CORRECT CELL ..............
			Storage.Reset();
			for (i=0;i<Particles.N;i++){
				// Enforce periodic BC
				if ( Particles.cx(i) < 0 ) Particles.cx(i) += Lx;
				if ( Particles.cy(i) < 0 ) Particles.cy(i) += Ly;
				if ( Particles.cz(i) < 0 ) Particles.cz(i) += Lz;
				if ( !(Particles.cx(i) < Lx) ) Particles.cx(i) -= Lx;
				if ( !(Particles.cy(i) < Ly) ) Particles.cy(i) -= Ly;
				if ( !(Particles.cz(i) < Lz) ) Particles.cz(i) -= Lz;
				// !!!!!! Check this guy !!!!!!!!
				icx = int(floor(Particles.cx(i)/Storage.lencx));
				icy = int(floor(Particles.cy(i)/Storage.lency));
				icz = int(floor(Particles.cz(i)/Storage.lencz));
//				cout << "particle to" << icx << ","<< icy << ","<< icz << ","<<endl;
				index = Storage.CellCount(icx,icy,icz)++;
				Storage.CellEntry(icx,icy,icz,index) = i;
			}
			// ..............................................................
			// .....  COMPUTE DISTANCE TO MOVE EACH PARTICLE ................ 
			// .....  Loop over all cells  ..................................
			for (icz=0;icz<ncz;icz++){
				for (icy=0;icy<ncy;icy++){
					for (icx=0;icx<ncx;icx++){
						// Get the number of partices in cell (icx,icy,icz)
						mainCellCount = Storage.CellCount(icx,icy,icz);
						// Loop over all particles in cell icx,icy,icz
						for (i=0;i<mainCellCount;i++){
							ii = Storage.CellEntry(icx,icy,icz,i);
#ifdef DEBUG
							cout << "PARTICLE " << ii << endl;
#endif
							// .............................................................
							//	COMPUTE ALL OVERLAPS FOR PARTICLE ii
							//	Loop over all neighboring cells (inx,iny,inz)
							d = Lx+Ly+Lz;
							ax = ay = az = 0;
							dxmin = dymin = dzmin = Lx+Ly+Lz;
							for (iz=-1;iz<2;iz++){
								for (iy=-1;iy<2;iy++){
									for (ix=-1;ix<2;ix++){
										bcx = bcy = bcz = 0.0;
										// Determine the neighbor cell (inx,iny,inz)
										inx = icx+ix;
										iny = icy+iy;
										inz = icz+iz;
										if (inx < 0){		// periodic BC
											inx += ncx;		// neighbor
											bcx = -Lx;		// shift 
										}
										if (iny < 0){		// periodic BC
											iny += ncy;		// neighbor
											bcy = -Ly;		// shift
										}
										if (inz < 0){		// periodic BC
											inz += ncz;		// neighbor
											bcz = -Lz;		// shift
										}
										if (!(inx < ncx)){	// periodic BC
											inx -= ncx;		// neighbor
											bcx = Lx;		// shift 
										}
										if (!(iny < ncy)){	// periodic BC
											iny -= ncy;		// neighbor
											bcy = Ly;		// shift
										}
										if (!(inz < ncz)){	// periodic BC
											inz -= ncz;		// neighbor
											bcz = Lz;		// shift
										}
										// Number of spheres in neighbor cell
										count = Storage.CellCount(inx,iny,inz);
										// Loop over all possible overlaps
										for (j=0;j<count;j++){
											jj = Storage.CellEntry(inx,iny,inz,j);
											// Don't compute overlap with ii
											if (ix == iy == iz == 0 && jj == ii) ;
											else{
												// Determine overlap between particles ii, jj
												// Compute the vector k - points toward ii
												kx = Particles.cx(ii) - bcx - Particles.cx(jj);
												ky = Particles.cy(ii) - bcy - Particles.cy(jj);
												kz = Particles.cz(ii) - bcz - Particles.cz(jj);
												// Compute the norm of k
												k = sqrt(kx*kx+ky*ky+kz*kz);
												// Determine magnitude of overlap: delta = diameter-k
												delta = Particles.Radius(jj)+Particles.Radius(ii)-k;
												if (delta > 0 ) {	// Particles ii, jj overlap
													overlaps++;
#ifdef DEBUG
													cout << "	overlap " << jj << ": "
													<<  Particles.cx(jj) + bcx << "," 
													<< Particles.cy(jj)+bcy << " " 
													<< Particles.cz(jj) + bcz << ", size "
													<< delta << endl;
#endif
													// Vector contribution to direction to move
													ax += delta*kx/k;
													ay += delta*ky/k;
													az += delta*kz/k;
													// Determine minimum overlap size 'd' (local to ii)
													if (delta < d) d = delta;	
													// Determine maximum overlap size (global)
													if (delta > max_overlap) max_overlap = delta;
													// determine minimum overlap vector
													if (delta*kx/k < dxmin) dxmin = delta*kx/k;
													if (delta*ky/k < dymin) dymin = delta*ky/k;
													if (delta*kz/k < dzmin) dzmin = delta*kz/k;
													if (delta > err) err = delta;
												}
											}
										}
									}
								}
							}	// End Loop over Neighbor cells
							// .............................................................
							// Determine the distance to move particle ii
							a = sqrt(ax*ax+ay*ay+az*az);	// norm of the direction vector
							// Choose the distance d to ensure that the smallest overlap is eliminated
							//d = 0.501*max(dxmin/ax,max(dymin/ay,dzmin/az));
							d = 0.601*a;					// distance to move particle ii
							// cout << "Particle " << ii  << ", size " << a << endl;
							if ( a > 0 ){
								Particles.dx(ii) = d*ax/a;		// store the displacement vector
								Particles.dy(ii) = d*ay/a;
								Particles.dz(ii) = d*az/a;
//								cout << "dx,dy,dz " << ii <<": " << Particles.dx(ii) << "" << Particles.dy(ii)
//								<<"" << Particles.dz(ii) << endl;
							}
							else {
								Particles.dx(ii) = Particles.dy(ii) = Particles.dz(ii) = 0.0;
							}
							// .............................................................
						}
						// .... displacement vector (dx,dy,dz) updated for particle i ......
					}
				}
			}
//			cout << "	 Iteration: "<< iter << ", error: " << err << " index " << jj << endl;
			
			// ..................................................................
			// ............. MOVE ALL OF THE PARTICLES ..........................
//			cout << "Radius: " << Particles.Radius(0) << endl;
			for (i=0;i<Particles.N;i++){
				Particles.cx(i) += Particles.dx(i);
				Particles.cy(i) += Particles.dy(i);
				Particles.cz(i) += Particles.dz(i);
			}

			// ..................................................................
		}
		cout << "Finished a loop in "<< iter <<" with error = " << err << endl;
		// ..................................................................
		// .................... COMPUTE THE POROSITY ........................
		V = 0;
		for (i=0;i<Particles.N;i++){
			V += 4*PI*pow(Particles.Radius(i),3)/3;
		}
		porosity = 1.0-V/Lx/Ly/Lz;
		cout << "	Current porosity is: " << porosity << endl;
		if (porosity < 0.4){
			factor = input_factor;
			cout << "Specified rescaling factor now being used: " << factor << endl;
		}
		if (iter == iteration_cutoff){
			// RESCALE (DECREASE) THE RADII TO ELIMINATE OVERLAPS
			cout << "Maxed out iterations... reverting to radii from previous timestep." << endl;
			M = Mprev;
			S = Sprev;
			for (i=0;i<Particles.N;i++){
				Particles.cx(i) = Save.cx(i);
				Particles.cy(i) = Save.cy(i);
				Particles.cz(i) = Save.cz(i);
				Particles.Radius(i) = Save.Radius(i);
			}
			failCount++;
			cout << "Failure number " << failCount << endl;
			if (failCount > 3){
				cout << "Failed 4 times, EXIT" << endl;
				time_to_stop = true;
			}
			else {
				cout << "Decreasing the radius factor from " << factor;
				factor = 1.0 + 0.5*(factor-1.0);
				cout << " to " << factor << endl;
				// If the radius factor gets too close to tolerance, give up
				if (factor-1.0 < .001){
					cout << "EXIT: factor became too small." << endl;
					time_to_stop = true;
				}
			}
		}
		else failCount = 0;
		// else	random_shift = false;
		if (porosity <= porosity_target)  time_to_stop = true;
		// ..................................................................
		// ............ SAVE THE CURRENT SPHERE ARRANGEMENT .................
		if (iter != iteration_cutoff) {
			for (i=0;i<Particles.N;i++) {
				Save.cx(i) = Particles.cx(i);
				Save.cy(i) = Particles.cy(i);
				Save.cz(i) = Particles.cz(i);
				Save.Radius(i) = Particles.Radius(i);
			}
		}
		// ..................................................................
	}
	// ...............  END OF ITERATIONS ...................................
	cout << "SIMULATION COMPLETE " << endl;
	V = 0;
	for (i=0;i<Save.N;i++){
		V += 4*PI*pow(Save.Radius(i),3)/3;
	}
	porosity = 1.0-V/Lx/Ly/Lz;
	// Compute the mean coordination number
	double CoordNo = CoordinationNumber(Particles,tol*10,Lx,Ly,Lz);
	
	cout << "Target Porosity was " << porosity_target << endl;
	cout << "Actual Porosity is " << porosity << endl;
	cout << "Mean coordination No. " << CoordNo << endl;
	cout << "FINAL DISTRIBUTION PARAMETERS: " << endl;
	cout << "	log(r) normally distributed with mean " << M << " and variance " << S << endl;
	
//	if (CHECK_OVERLAPS(Save, Storage, tol*100) == true){
//		cout << "WARNING: overlaps exist greater than tolerance!" << endl;
//	}
	
	double max = 0.0;
	for (i=0; i<Particles.N; i++){
		if (Particles.Radius(i) > max)	max = Particles.Radius(i);
	}
	if ( max > 0.5*Storage.lencx || max > 0.5*Storage.lency || max > 0.5*Storage.lencz){
		cout << "WARNING: maximum radius exceeds bin width!" << endl;
	}

	// ........  WRITE THE PACKING TO OUTPUT FILE ...........................
	ofstream pack ("pack.out");
	pack << "Number of Spheres: " << Particles.N << endl;
	pack << "Domain Length (x,y,z):  " << Lx <<", "<< Ly << ", "<< Lz << endl; 
	pack << "Media porosity: " << porosity << endl;
	pack << "log(r) Normal with mean " << M << ", variance " << S << endl;
	pack << "Mean coordination No. " << CoordNo << endl;

	for (i=0; i<Particles.N; i++){
//#ifdef DEBUG
//		pack << Particles.cx(i) << " " << Particles.cy(i) << " " 
//			<< Particles.cz(i) << " " << Particles.Radius(i) << endl;
//#else
		pack << Save.cx(i) << " " << Save.cy(i) << " " 
		<< Save.cz(i) << " " << Save.Radius(i) << endl;
//#endif
	}
}
// **************** END OF THE MAIN CODE *********************************************
inline double CoordinationNumber(SphereCollection &Particles, double size,	
								 double Lx, double Ly, double Lz)
{
	double dist;
	int count;
	count = 0;
	for (int i=0;i<Particles.N;i++){
		// Compute the number of contacts for particle i
		for (int j=0;j<Particles.N;j++){
			// compute distance
			dist = sqrt(pow(Particles.cx(i)-Particles.cx(j),2)
					+ pow(Particles.cy(i)-Particles.cy(j),2)
					+ pow(Particles.cz(i)-Particles.cz(j),2))
					- Particles.Radius(i)-Particles.Radius(j);
			if ( dist < size ){
				count++;
			}
			// compute distance ( -x Boundary Condition )
			dist = sqrt(pow(Particles.cx(i)-Particles.cx(j)-Lx,2)
						+ pow(Particles.cy(i)-Particles.cy(j),2)
						+ pow(Particles.cz(i)-Particles.cz(j),2))
			- Particles.Radius(i)-Particles.Radius(j);
			if ( dist < size ){
				count++;
			}	
			// compute distance ( +x Boundary Condition )
			dist = sqrt(pow(Particles.cx(i)-Particles.cx(j)+Lx,2)
						+ pow(Particles.cy(i)-Particles.cy(j),2)
						+ pow(Particles.cz(i)-Particles.cz(j),2))
			- Particles.Radius(i)-Particles.Radius(j);
			if ( dist < size ){
				count++;
			}
			// compute distance ( -y Boundary Condition )
			dist = sqrt(pow(Particles.cx(i)-Particles.cx(j),2)
						+ pow(Particles.cy(i)-Particles.cy(j)-Ly,2)
						+ pow(Particles.cz(i)-Particles.cz(j),2))
			- Particles.Radius(i)-Particles.Radius(j);
			if ( dist < size ){
				count++;
			}
			// compute distance ( +y Boundary Condition )
			dist = sqrt(pow(Particles.cx(i)-Particles.cx(j),2)
						+ pow(Particles.cy(i)-Particles.cy(j)+Ly,2)
						+ pow(Particles.cz(i)-Particles.cz(j),2))
			- Particles.Radius(i)-Particles.Radius(j);
			if ( dist < size ){
				count++;
			}	
			// compute distance ( -z Boundary Condition )
			dist = sqrt(pow(Particles.cx(i)-Particles.cx(j),2)
						+ pow(Particles.cy(i)-Particles.cy(j),2)
						+ pow(Particles.cz(i)-Particles.cz(j)-Lz,2))
			- Particles.Radius(i)-Particles.Radius(j);
			if ( dist < size ){
				count++;
			}	
			// compute distance ( -z Boundary Condition )
			dist = sqrt(pow(Particles.cx(i)-Particles.cx(j),2)
						+ pow(Particles.cy(i)-Particles.cy(j),2)
						+ pow(Particles.cz(i)-Particles.cz(j)+Lz,2))
			- Particles.Radius(i)-Particles.Radius(j);
			if ( dist < size ){
				count++;
			}
			dist = sqrt(pow(Particles.cx(i)-Particles.cx(j),2)
						+ pow(Particles.cy(i)-Particles.cy(j),2)
						+ pow(Particles.cz(i)-Particles.cz(j),2))
			- Particles.Radius(i)-Particles.Radius(j);
			if ( dist < size ){
				count++;
			}
			// compute distance ( -x Boundary Condition )
			dist = sqrt(pow(Particles.cx(i)-Particles.cx(j)-Lx,2)
						+ pow(Particles.cy(i)-Particles.cy(j)-Ly,2)
						+ pow(Particles.cz(i)-Particles.cz(j),2))
			- Particles.Radius(i)-Particles.Radius(j);
			if ( dist < size ){
				count++;
			}	
			// compute distance ( +x Boundary Condition )
			dist = sqrt(pow(Particles.cx(i)-Particles.cx(j)+Lx,2)
						+ pow(Particles.cy(i)-Particles.cy(j)-Ly,2)
						+ pow(Particles.cz(i)-Particles.cz(j),2))
			- Particles.Radius(i)-Particles.Radius(j);
			if ( dist < size ){
				count++;
			}
			// compute distance ( -y Boundary Condition )
			dist = sqrt(pow(Particles.cx(i)-Particles.cx(j)+Lx,2)
						+ pow(Particles.cy(i)-Particles.cy(j)-Ly,2)
						+ pow(Particles.cz(i)-Particles.cz(j),2))
			- Particles.Radius(i)-Particles.Radius(j);
			if ( dist < size ){
				count++;
			}
			// compute distance ( +y Boundary Condition )
			dist = sqrt(pow(Particles.cx(i)-Particles.cx(j)+Lx,2)
						+ pow(Particles.cy(i)-Particles.cy(j)+Ly,2)
						+ pow(Particles.cz(i)-Particles.cz(j),2))
			- Particles.Radius(i)-Particles.Radius(j);
			if ( dist < size ){
				count++;
			}	
			// compute distance ( -z Boundary Condition )
			dist = sqrt(pow(Particles.cx(i)-Particles.cx(j)-Lx,2)
						+ pow(Particles.cy(i)-Particles.cy(j),2)
						+ pow(Particles.cz(i)-Particles.cz(j)-Lz,2))
			- Particles.Radius(i)-Particles.Radius(j);
			if ( dist < size ){
				count++;
			}	
			// compute distance ( -z Boundary Condition )
			dist = sqrt(pow(Particles.cx(i)-Particles.cx(j)+Lx,2)
						+ pow(Particles.cy(i)-Particles.cy(j),2)
						+ pow(Particles.cz(i)-Particles.cz(j)+Lz,2))
			- Particles.Radius(i)-Particles.Radius(j);
			if ( dist < size ){
				count++;
			}	
			// compute distance ( -z Boundary Condition )
			dist = sqrt(pow(Particles.cx(i)-Particles.cx(j)+Lx,2)
						+ pow(Particles.cy(i)-Particles.cy(j),2)
						+ pow(Particles.cz(i)-Particles.cz(j)-Lz,2))
			- Particles.Radius(i)-Particles.Radius(j);
			if ( dist < size ){
				count++;
			}	
			// compute distance ( -z Boundary Condition )
			dist = sqrt(pow(Particles.cx(i)-Particles.cx(j)-Lx,2)
						+ pow(Particles.cy(i)-Particles.cy(j),2)
						+ pow(Particles.cz(i)-Particles.cz(j)+Lz,2))
			- Particles.Radius(i)-Particles.Radius(j);
			if ( dist < size ){
				count++;
			}	
			// compute distance ( -z Boundary Condition )
			dist = sqrt(pow(Particles.cx(i)-Particles.cx(j),2)
						+ pow(Particles.cy(i)-Particles.cy(j)-Ly,2)
						+ pow(Particles.cz(i)-Particles.cz(j)-Lz,2))
			- Particles.Radius(i)-Particles.Radius(j);
			if ( dist < size ){
				count++;
			}	
			// compute distance ( -z Boundary Condition )
			dist = sqrt(pow(Particles.cx(i)-Particles.cx(j),2)
						+ pow(Particles.cy(i)-Particles.cy(j)+Ly,2)
						+ pow(Particles.cz(i)-Particles.cz(j)+Lz,2))
			- Particles.Radius(i)-Particles.Radius(j);
			if ( dist < size ){
				count++;
			}	
			// compute distance ( -z Boundary Condition )
			dist = sqrt(pow(Particles.cx(i)-Particles.cx(j),2)
						+ pow(Particles.cy(i)-Particles.cy(j)+Ly,2)
						+ pow(Particles.cz(i)-Particles.cz(j)-Lz,2))
			- Particles.Radius(i)-Particles.Radius(j);
			if ( dist < size ){
				count++;
			}	
			// compute distance ( -z Boundary Condition )
			dist = sqrt(pow(Particles.cx(i)-Particles.cx(j),2)
						+ pow(Particles.cy(i)-Particles.cy(j)-Ly,2)
						+ pow(Particles.cz(i)-Particles.cz(j)+Lz,2))
			- Particles.Radius(i)-Particles.Radius(j);
			if ( dist < size ){
				count++;
			}	
			// compute distance ( -z Boundary Condition )
			dist = sqrt(pow(Particles.cx(i)-Particles.cx(j)-Lx,2)
						+ pow(Particles.cy(i)-Particles.cy(j)-Ly,2)
						+ pow(Particles.cz(i)-Particles.cz(j)-Lz,2))
			- Particles.Radius(i)-Particles.Radius(j);
			if ( dist < size ){
				count++;
			}	
			// compute distance ( -z Boundary Condition )
			dist = sqrt(pow(Particles.cx(i)-Particles.cx(j)-Lx,2)
						+ pow(Particles.cy(i)-Particles.cy(j)+Ly,2)
						+ pow(Particles.cz(i)-Particles.cz(j)+Lz,2))
			- Particles.Radius(i)-Particles.Radius(j);
			if ( dist < size ){
				count++;
			}	
			// compute distance ( -z Boundary Condition )
			dist = sqrt(pow(Particles.cx(i)-Particles.cx(j)-Lx,2)
						+ pow(Particles.cy(i)-Particles.cy(j)+Ly,2)
						+ pow(Particles.cz(i)-Particles.cz(j)-Lz,2))
			- Particles.Radius(i)-Particles.Radius(j);
			if ( dist < size ){
				count++;
			}	
			// compute distance ( -z Boundary Condition )
			dist = sqrt(pow(Particles.cx(i)-Particles.cx(j)-Lx,2)
						+ pow(Particles.cy(i)-Particles.cy(j)-Ly,2)
						+ pow(Particles.cz(i)-Particles.cz(j)+Lz,2))
			- Particles.Radius(i)-Particles.Radius(j);
			if ( dist < size ){
				count++;
			}
			// compute distance ( -z Boundary Condition )
			dist = sqrt(pow(Particles.cx(i)-Particles.cx(j)+Lx,2)
						+ pow(Particles.cy(i)-Particles.cy(j)-Ly,2)
						+ pow(Particles.cz(i)-Particles.cz(j)-Lz,2))
			- Particles.Radius(i)-Particles.Radius(j);
			if ( dist < size ){
				count++;
			}	
			// compute distance ( -z Boundary Condition )
			dist = sqrt(pow(Particles.cx(i)-Particles.cx(j)+Lx,2)
						+ pow(Particles.cy(i)-Particles.cy(j)+Ly,2)
						+ pow(Particles.cz(i)-Particles.cz(j)+Lz,2))
			- Particles.Radius(i)-Particles.Radius(j);
			if ( dist < size ){
				count++;
			}	
			// compute distance ( -z Boundary Condition )
			dist = sqrt(pow(Particles.cx(i)-Particles.cx(j)+Lx,2)
						+ pow(Particles.cy(i)-Particles.cy(j)+Ly,2)
						+ pow(Particles.cz(i)-Particles.cz(j)-Lz,2))
			- Particles.Radius(i)-Particles.Radius(j);
			if ( dist < size ){
				count++;
			}	
			// compute distance ( -z Boundary Condition )
			dist = sqrt(pow(Particles.cx(i)-Particles.cx(j)+Lx,2)
						+ pow(Particles.cy(i)-Particles.cy(j)-Ly,2)
						+ pow(Particles.cz(i)-Particles.cz(j)+Lz,2))
			- Particles.Radius(i)-Particles.Radius(j);
			if ( dist < size ){
				count++;
			}
		}
	}
	// Compute the averge coordination number 
	double cn = double (count) / double(Particles.N);
	return cn;
}
inline void Random_Displacement(SphereCollection &Particles, double size)
{
	int i;
	for (i=0;i<Particles.N;i++){
		Particles.cx(i) + size*(2*rand()/RAND_MAX - 1)*Particles.Radius(i);
		Particles.cy(i) + size*(2*rand()/RAND_MAX - 1)*Particles.Radius(i);
		Particles.cz(i) + size*(2*rand()/RAND_MAX - 1)*Particles.Radius(i);
	}
}

inline bool CHECK_OVERLAPS(SphereCollection &Particles, CellStorage &Storage, double tol)
{
	bool toReturn = false;
	int icx,icy,icz,jcx,jcy,jcz,number,count,i,j;
	int ii,jj;
	// Check for any overlaps
	for (icz=0; icz<Storage.ncz; icz++){
		for (icy=0; icy<Storage.ncy; icy++){
			for (icx=0; icx<Storage.ncx; icx++){
				// How many particles in this cell 
				number = Storage.CellCount(icx,icy,icz);
				for (i=0; i<number; i++){
					ii = Storage.CellEntry(icx,icy,icz,i);
					// Go over all the other cells 
					for (jcz=0; jcz<Storage.ncz; jcz++){
						for (jcy=0; jcy<Storage.ncy; jcy++){
							for (jcx=0; jcx<Storage.ncx; jcx++){
								count = Storage.CellCount(jcx,jcy,jcz);
								for (j=0; j<count; j++){
									jj = Storage.CellEntry(jcx,jcy,jcz,j);
									// Check for overlap between ii & jj
									if ( icx == jcx && icy == jcy && icz == jcz && ii==jj){
										// Same particle - don't count overlap with self
									}
									else {
										if ( sqrt( pow(Particles.cx(ii)-Particles.cx(jj),2)
												  +pow(Particles.cy(ii)-Particles.cy(jj),2)
												  +pow(Particles.cz(ii)-Particles.cz(jj),2) )
											- Particles.Radius(ii) - Particles.Radius(jj) > tol){
#ifdef DEBUG
											cout << "ERROR: Overlap between" << ii << " and " << jj << endl; 
#endif
											toReturn = true;
											
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	if ( toReturn == true){
		cout << "WARNING: Overlaps exist,  consider decreasing the number of cells!" << endl;
	}
	return toReturn;
}

