#ifndef myHistogram_h
#define myHistogram_h 1

#include <string>
#include <fstream>
#include "G4Threading.hh"

// Written by Grant Berland
// Header-only histogramming class to record energy deposition per altitude bin

class myHistogram
{
public:

  myHistogram();
  myHistogram(double, double, int);
  
  ~myHistogram();

  // Tells program whether or not to write to histogram if other 
  // data collections methods are selected
  // void InitializeHistogram(){ initializedFlag = 1; };
  
  // Method to fill histogram with data
  void AddCountToBin(unsigned int, double);

  void AddCountTo2DHistogram(unsigned int, double);

  void WriteHistogramToFile(std::string); // Writes vector to file name provided

  void Write2DHistogram(std::string);
  
  // Always base 10
  void GenerateLogspaceBins(double, double, int, double[]);

private:
  // Need to run InitializeHistogram() method to write particle energy depostion to histogram
  // int initializedFlag = 0; // TODO delete?
  
  // Array initialized to zeros (with fixed resolution)
  double histogramArray[1000] = {};
  double binEdges[101] = {};

  double twoDhistogramArray[1000][100] = {}; // TODO delete?
};

// Inline constructor and destructor methods
inline myHistogram::myHistogram()
  : histogramArray()
{}

inline myHistogram::myHistogram(double binStart, double binStop, int nBins) 
  : binEdges(),
    twoDhistogramArray()
{
  GenerateLogspaceBins(binStart, binStop, nBins, binEdges);
}

inline myHistogram::~myHistogram(){}

inline void myHistogram::AddCountToBin(unsigned int binAddress, double amountToAdd)
{
  histogramArray[binAddress] += amountToAdd;
}

inline void myHistogram::WriteHistogramToFile(std::string filename)
{
  // Open file
  std::ofstream outputFile;
  outputFile.open(filename, std::ios_base::app);

  // Write data
  outputFile << "altitude_km,energy_deposition_kev" << "\n"; // Header
  for(unsigned int i=0; i<1000; i++)
  {
    outputFile << i << "," << histogramArray[i] << "\n";
  }
  // Close file
  outputFile.close();
}

inline void myHistogram::AddCountTo2DHistogram(unsigned int address1, double value)
{
  for(unsigned int i=0; i<101; i++)
  {   
    if(binEdges[i] > value)
    { 
      twoDhistogramArray[address1][i] += 1;
      break; 
    }
  }
}

inline void myHistogram::Write2DHistogram(std::string filename) // TODO delete?
{
  // Open file
  std::ofstream outputFile;
  outputFile.open(filename, std::ios_base::app);

  // Write data
  for(unsigned int i=0; i<1000; i++)
  {
    for (unsigned int j=0; j<100-1; j++)
    {
        outputFile << twoDhistogramArray[i][j] << ",";
    }
    outputFile << twoDhistogramArray[i][100-1] << "\n";
  }
  
  // Close file
  outputFile.close();
}

inline void myHistogram::GenerateLogspaceBins(G4double start, G4double end, G4int nBins, G4double binArray[])
{
  // step size 
  double c = (end - start) / (nBins - 1);
  
  // fill vector 
  for (int i = 0; i < nBins-1; ++i)
  {
    binArray[i] = std::pow(10., start + i * c);
  }

  // fix last entry to 10^b 
  binArray[nBins-1] = std::pow(10., end);

}

#endif