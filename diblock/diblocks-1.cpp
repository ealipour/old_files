#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <cstdlib>
#include <sstream>
#include <string>
#include <stdint.h>

using namespace std;

const double PI = 4*atan(1.);

const int maxAtoms = 500000;
const int maxBonds = 400000;
const int nUnitCellx = 57;
const int nUnitCelly = 33;

int main(int argc, char* argv[])
{
  cout << "Start" << endl;


  // Memory

  vector <double> xPos; xPos.resize(maxAtoms);
  vector <double> yPos; yPos.resize(maxAtoms);
  vector <double> zPos; zPos.resize(maxAtoms);
  vector <double> type; type.resize(maxAtoms);
  vector <double> mole; mole.resize(maxAtoms);

  vector <double> bondType; bondType.resize(maxBonds);
  vector <double> bondAtom1; bondAtom1.resize(maxBonds);
  vector <double> bondAtom2; bondAtom2.resize(maxBonds);

  int counter = 0;
  int index;
  string line;
  double lengthBox[3];
  double hBox[3];
  int bondCounter = 0;
  int molCounter = 0;

  lengthBox[0] = 43;
  lengthBox[1] = 43;
  lengthBox[2] = 43;

  hBox[0] = lengthBox[0]/2;
  hBox[1] = lengthBox[1]/2;
  hBox[2] = lengthBox[2]/2;


  // Insert the diblocks (SC lattice)

  int nGridx = 2;
  int nGridy = 40;
  int nGridz = 40;
  int nBead = 20;
  double gridLengthX = lengthBox[0]/nGridx;
  double gridLengthY = lengthBox[1]/nGridy;
  double gridLengthZ = lengthBox[2]/nGridz;

  for (int iGridx = 0; iGridx < nGridx; iGridx++)
  {
    for (int iGridy = 0; iGridy < nGridy; iGridy++)
    {
      for (int iGridz = 0; iGridz < nGridz; iGridz++)
      {
        for (int iBead = 0; iBead < nBead; iBead++)
        {
          // The Atom positions 
          xPos[counter] = iGridx*gridLengthX + iBead;
          yPos[counter] = iGridy*gridLengthY;
          zPos[counter] = iGridz*gridLengthZ;
          mole[counter] = molCounter;
          type[counter] = 1;
          if (iBead<nBead/2)
	  {
	    type[counter]=2;
	  }
          if (iBead > 0)
          {
	    // The bonds lower slab
	    bondType[bondCounter] = 1;
	    bondAtom1[bondCounter] = counter;
	    bondAtom2[bondCounter] = counter + 1;
            bondCounter++;
          }
          // Counters up
          counter++;
        }
        molCounter++;
      }
    }
  }

  int nAtom = counter;
  int nBonds = bondCounter;

  // Initialise and write the outputfile in LAMMPS format

  ofstream myfileLAMMPS;
  myfileLAMMPS.open ("diblocks.pos");
  myfileLAMMPS << "# The init config for the diblocks" << endl;
  myfileLAMMPS << nAtom << "  atoms" << endl;
  myfileLAMMPS << nBonds << "  bonds" << endl;
  int nType = 2;
  int nBondType = 1;
  myfileLAMMPS << nType << "  atom types" << endl;
  myfileLAMMPS << nBondType << "  bond types" << endl;
  myfileLAMMPS << endl;
  double xlo = 0.0, ylo = 0.0, zlo = 0.0;
  myfileLAMMPS << xlo << " " << lengthBox[0] << "  xlo xhi" << endl;
  myfileLAMMPS << ylo << " " << lengthBox[1] << "  ylo yhi" << endl;
  myfileLAMMPS << zlo << " " << lengthBox[2] << "  zlo zhi" << endl;
  myfileLAMMPS << endl;
  double mass = 1.0;
  myfileLAMMPS << " Masses" << endl;
  myfileLAMMPS << endl;
  for (int iType = 0; iType < nType; iType++)
  {
    myfileLAMMPS << (iType+1) << " " << mass << endl;
  }
  myfileLAMMPS << endl;
  myfileLAMMPS << " Atoms" << endl;
  myfileLAMMPS << endl;
  for (int iAtom = 0; iAtom < nAtom; iAtom++)
  {
    myfileLAMMPS << (iAtom+1) << " " << mole[iAtom] << " " << type[iAtom] << " " <<
         xPos[iAtom] << " " << yPos[iAtom] << " " << zPos[iAtom] << endl;
  }
  myfileLAMMPS << endl;
  myfileLAMMPS << " Bonds" << endl;
  myfileLAMMPS << endl;
  for (int iBond = 0; iBond < nBonds; iBond++)
  {
    myfileLAMMPS << (iBond+1) << " " << bondType[iBond] << " " << 
         bondAtom1[iBond] << " " << bondAtom2[iBond] << endl;
  }
  myfileLAMMPS << endl;

  myfileLAMMPS.close();

  
  cout << "Klaar!!" << endl;
  
  return (0);

}
