// classes example
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <ctime>
#include <unistd.h>
#include <cstdlib>
#include <sstream>
#include <cmath>
#include <ctime>

using namespace std;


///
/// This program calculates the distribution of: 
///
/// (1) Number of first neigbours
/// (2) First neighbours distances
///
/// Based on the connectivity criterion imposed by topological clustering and NOT
/// on the distance based clustering
///



//**********************************************************************************
// ATOM class and methods
//**********************************************************************************


// defining Atom class interface
//    -> functions and properties of the class Atom 
class Atom {
	int index ; /// index = index position in the list as given by the file
	int type ;  /// type = type of atom
	double xpos, ypos, zpos ;
   public:
	void SetType(int) ;
	int GetType() ;
	void SetIndex(int) ;
	int GetIndex() ;
	
	void SetX(double) ;
	double GetX() ;
	void SetY(double) ;
	double GetY() ;
	void SetZ(double) ;
	double GetZ() ;
};


void Atom::SetType (int tt) {
	type = tt ;
}

int Atom::GetType () {
	return type ;
}


void Atom::SetIndex (int ii) {
	index = ii ;
}


int Atom::GetIndex () {
	return index ;
}

// Atom class method
//   method to save/set the Atom X position
void Atom::SetX (double xx) {
	xpos = xx ;
}

// Atom class method
//   method to get the atom X position
double Atom::GetX () {
	return xpos ;
}

// Atom class method
//   method to save/set the Atom X position
void Atom::SetY (double yy) {
	ypos = yy ;
}

// Atom class method
//   method to get the atom X position
double Atom::GetY () {
	return ypos ;
}

// Atom class method
//   method to save/set the Atom X position
void Atom::SetZ (double zz) {
	zpos = zz ;
}

// Atom class method
//   method to get the atom X position
double Atom::GetZ () {
	return zpos ;
}



//**********************************************************************************
// AtomSet class and methods
//**********************************************************************************
class AtomSet
{
	private:
		size_t AtomSetCard ;
		std::vector<Atom> AtomList ;
	public:
		AtomSet(size_t dim_AtomSetCard) ; // Constructor
		Atom GetAtom(size_t dim_AtomSetCard);
};

// CONSTRUCTOR Building atom objects and placing them into the AtomSet object (Atom list vector)
AtomSet::AtomSet(size_t dim_AtomSetCard) : AtomSetCard(dim_AtomSetCard), AtomList()
{
	AtomList.resize(AtomSetCard);
	
	std::ifstream file( "traj_100000.xyz" ) ; // input file
	std::string line ;    // variable string to read the line
	unsigned int element ;
	double x, y, z ;
	
	
	for ( size_t k = 0; k < AtomSetCard; ++k ) {
		std::getline(file, line) ;
		std::stringstream aa(line) ;
		aa >> element >> x >> y >> z ;
		// creating an Atom object that will contain the type and position
		Atom Atom_object = Atom() ;
		// Assingning the atom type to the Atom object
		Atom_object.SetType(element) ;
		Atom_object.SetIndex(k) ;
		Atom_object.SetX(x) ;
		Atom_object.SetY(y) ;
		Atom_object.SetZ(z) ;
		// asigning this atom object to the i position of the atom list
		AtomList[k] = Atom_object ;
	}
}

// defining AtomSet class methods
//   method get the Atom in a location of the AtomList
Atom AtomSet::GetAtom(size_t dim_AtomSetCard)
{
	return AtomList[dim_AtomSetCard];
}












	

//**********************************************************************************
// main function
//**********************************************************************************
int main()
{
	// Parameters:
	//----------------------------------------------------------------------------------------------------------------------------
	//----------------------------------------------------------------------------------------------------------------------------	
	unsigned int NodeType = 1 ; // Tells the program what type of atoms will be clustered all other will be consider as sea atoms.
	size_t Natoms = 100000 ; // Total number of atoms in the file (Node atoms + Sea atoms).
	double R1 = 2.5 ; // Radius of type 1 atoms.
	double R2 = 1.7 ; // Radius of type 2 atoms.
	//----------------------------------------------------------------------------------------------------------------------------
	//----------------------------------------------------------------------------------------------------------------------------		
	
	
	
	const double pi = 3.1415926535897 ;
	double RT ;
	
	double ijx, ijy, ijz, ikx, iky, ikz, modij, Vol, ikpar, modik, ikver, den ;
	unsigned int Atom_count = 0 ;
	unsigned int memo = 0 ;
	unsigned int count = 0 ;
	unsigned int ifriends = 0 ;
	
	// Count execution time
	int start_s=clock();
	
	ofstream Frequencies ;
	Frequencies.open ("Frequencies.dat", ios::out | ios::trunc); // app=append
	ofstream Distances ;
	Distances.open ("Distances.dat", ios::out | ios::trunc); // app=append
	
	// Define the object set of atoms called AtomSet_object with Natoms
	// It is automatically initialized with the data in the file that we 
	// defined in the constructor
	
	AtomSet AtomSet_object(Natoms);
	
	for ( size_t i = 0; i < Natoms; ++i ) {
		if ( AtomSet_object.GetAtom(i).GetType() == NodeType ) {
			ifriends = 0 ;
			for ( size_t j = 0; j < Natoms; ++j ) {
				if ( AtomSet_object.GetAtom(j).GetType() == NodeType ) { 
					ijx = AtomSet_object.GetAtom(j).GetX() - AtomSet_object.GetAtom(i).GetX() ;
					ijy = AtomSet_object.GetAtom(j).GetY() - AtomSet_object.GetAtom(i).GetY() ;
					ijz = AtomSet_object.GetAtom(j).GetZ() - AtomSet_object.GetAtom(i).GetZ() ;
					modij = sqrt (ijx*ijx + ijy*ijy + ijz*ijz) ;
					Vol = pi*RT*RT*modij ;
					ijx = ijx/modij ;
					ijy = ijy/modij ;
					ijz = ijz/modij ;
				
					Atom_count = 0 ;
					
					for ( size_t k = 0; k < Natoms; ++k ) {
						if ( (k != i) && (k != j) ) {						
							ikx = AtomSet_object.GetAtom(k).GetX() - AtomSet_object.GetAtom(i).GetX() ;
							iky = AtomSet_object.GetAtom(k).GetY() - AtomSet_object.GetAtom(i).GetY() ;	
							ikz = AtomSet_object.GetAtom(k).GetZ() - AtomSet_object.GetAtom(i).GetZ() ;	
							modik = sqrt (ikx*ikx + iky*iky + ikz*ikz) ;						
							ikpar = ijx*ikx + ijy*iky + ijz*ikz ;
							ikver = sqrt (modik*modik - ikpar*ikpar) ;						
							// k atoms inside the cylinder
							unsigned int typei = AtomSet_object.GetAtom(i).GetType() ;
							unsigned int typek = AtomSet_object.GetAtom(k).GetType() ;
							
							if ( (typek ==  typei) ) { 
								if ( typei == 1) {
									RT = R1 + R1 ;
								}
								else {
									RT = R2 + R2 ;
								}
							}
							else {
								RT = R1 + R2 ;
							}				
							
							
							if ( (ikpar > 0) && (ikpar < modij) && (ikver < RT)) {
								++Atom_count ;
							}// If atom inside cylinder i-j then count one atom to Atom_count
						}// If k is different from i or j			
					}// Sum over all atoms
					
					// This If can be removed to speed up just kept for strict consistency
					if ( i != j ) { 
						// Conectivity criterion
						den = Atom_count/Vol ;
						//cout << i << " " << j << " " << den << endl ;
						// computing i-friends
						if ( den < 0.000000001 ) {
						++ifriends ;				
						//cout << " density = " << i+1 << " " << j+1 << " " << ifriends << " " << den << endl ;	
						Distances << modij << endl ; // Store in file Distances.dat
						}// If density is lower than 0.000...
					}
			
			
			
				}// i and j are nodes (type=1)
			}// End j
			Frequencies << ifriends << endl ; // Store in file Frequencies.dat	
			cout << "Atom " << i+1 << " has " << ifriends << " first neighours" << endl ;		
		}	
	}// End i
	
	Frequencies.close() ;
	Distances.close() ;

	
	
	//cout << "\033[1m Results: \033[0m" << endl ;
	cout << endl;
	//cout << "Mean number of first nighbours = \033[1m" << "\033[0m" << endl ;
	cout << endl;
	//cout << "Standard deviation of the number of first nighbours = \033[1m" << "\033[0m" << endl ;		
	
	// End counting time 	
 	int stop_s=clock();
 	
 	cout << endl;
	cout << "Execution time: \033[1m" << ( stop_s - start_s )/double(CLOCKS_PER_SEC) << "\033[0m seconds" << endl;
	cout << endl;
	
}




// Add construction of Histograms?

// count:
////////http://www.cplusplus.com/reference/algorithm/count/
//////// count algorithm example
//////#include <iostream>     // std::cout
//////#include <algorithm>    // std::count
//////#include <vector>       // std::vector

//////int main () {
//////  // counting elements in array:
//////  int myints[] = {10,20,30,30,20,10,10,20};   // 8 elements
//////  int mycount = std::count (myints, myints+8, 10);
//////  std::cout << "10 appears " << mycount << " times.\n";

//////  // counting elements in container:
//////  std::vector<int> myvector (myints, myints+8);
//////  mycount = std::count (myvector.begin(), myvector.end(), 20);
//////  std::cout << "20 appears " << mycount  << " times.\n";

//////  return 0;
//////}

