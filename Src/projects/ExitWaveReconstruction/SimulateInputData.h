#ifndef EWR_SIMULATIONBASE_H
#define EWR_SIMULATIONBASE_H

/// \cond
#include "Physics.h"
#include "Image.h"
#include "Parameters.h"
/// \endcond
#include "OpenCLKernelData.h"
#include "SimulateImage.h"
/// \cond
#include "Utilities.h"

#include <omp.h>

#include <boost/math/special_functions/bessel.hpp>
/// \endcond

/**
 * Calculates the projected potential of a given list of point charges with respect to the
 * parameters in Param in the unit Volt * nanometer.
 * 
 * The parameters pointChargeX and pointChargeY are expected to be in nanometers and the
 * vector charge contains the multiplier for the electric charge of the point charges as
 * multiples of the elementary charge. (Ex.: if charge[0] = -4 then the first point charge
 * in the list has a (negative) charge of -4 * Constants::e, where Constants::e is the
 * elementary charge as defined in Physics.h.)
 */
template <typename RealType>
void CalculateProjectedPointChargePotential ( const Parameters<RealType>& Param,
                                              qc::ScalarArray<RealType, qc::QC_2D>& ProjectedPotential,
                                              const vector<int>& charge,
                                              const vector<RealType>& chargeX,
                                              const vector<RealType>& chargeY,
                                              bool verbose = false ) {
  // Check vector sizes
  if ( charge.size ( ) != chargeX.size ( ) || charge.size ( ) != chargeY.size ( ) )
    throw aol::Exception ( "CalculateProjectedPointChargePotential: invalid parameters!", __FILE__, __LINE__ );
  
  // Calculate the projected potential of the point charge lattice
  if ( verbose )
    cerr << "\tCalculating the projected potential... " << setw ( 3 ) << 0 << "%";
  
  for ( int y = 0; y < Param.Y ; y++ ) {
    for ( int x = 0; x < Param.X ; x++ ) {
      const RealType posX = ( x - Param.X / 2 ) * Param.lenX / Param.X;
      const RealType posY = ( y - Param.Y / 2 ) * Param.lenY / Param.Y;
      
      RealType arg = 1;
      for ( unsigned int j = 0; j < charge.size ( ) ; j++ ) {
        // Calculate the squared distance of the pixel (x,y) to the j-th point charge in nm^2
        RealType distSqr = aol::Sqr<RealType> ( posX - chargeX[j] ) + aol::Sqr<RealType> ( posY - chargeY[j] );
        
        // Avoid overflow errors if the distance is too small, because the projected potential goes to +-infinity close the the point charges
        if ( distSqr < 1e-5 )
          distSqr = 1e-5;
        
        arg *= pow ( distSqr, -charge[j] );
      }
      
      // Calculate the projected potential
      ProjectedPotential.set ( x, y, 1 / ( 4 * aol::NumberTrait<RealType>::pi * Constants::eps0 ) * Constants::e * log ( arg ) );
    }
    
    if ( verbose )
      cerr << "\r\tCalculating the projected potential... " << setw ( 3 ) << static_cast<int> ( ( ( y + 1 ) / ( Param.Y ) ) * 100 ) << "%";
  }
  
  if ( verbose )
    cerr << " done." << endl;
}

/**
 * \brief Reads atomic potential coefficients from a parameter file
 * 
 * These coefficients are used for the calculation of the projected potential of atoms.
 * The coefficients are given in E. J. Kirkland, Advanced Computing in Electron Microscopy
 * (Appendix D).
 */
template <typename RealType>
void GetPotentialCoefficients ( const string& path, const vector<int>& atomicNumber, vector<aol::MultiVector<RealType>>& coeff ) {
  aol::ParameterParser parser ( path );
  
  coeff.resize ( atomicNumber.size ( ) );
  for ( unsigned int i = 0; i < atomicNumber.size ( ) ; i++ ) {
    // Read coefficients for each atom from the file
    string variable = string ( "z" ) + valToString ( atomicNumber[i] );
    
    if ( !parser.hasVariable ( variable.c_str ( ) ) )
      throw aol::Exception ( string ( "Potential coefficients for the atomic number " ) + valToString ( atomicNumber[i] ) + " not found!", __FILE__, __LINE__ );
    
    parser.getRealMultiVec<RealType> ( variable.c_str ( ), coeff[i] );
    
    // Scale coefficients to the correct magnitude (from Angstrom to nanometer resp. from 1 / Angstrom to 1 / nanometer)
    coeff[i][0] *= 10;
    coeff[i][1] *= 100;
    coeff[i][2] /= 10;
    coeff[i][3] /= 100;
  }
}

/**
 * Calculates the projected atomic potential of the specimen given by the parameters in Param
 * in the unit Volt * nanometer. The parameters atomPosX and atomPosY are expected to be in
 * nanometers.
 * 
 * Formula taken from E. J. Kirkland, Advanced Computing in Electron Microscopy.
 */
template <typename RealType>
static void CalculateProjectedAtomicPotential ( const Parameters<RealType>& Param,
                                                qc::ScalarArray<RealType, qc::QC_2D>& ProjectedPotential,
                                                const vector<int>& atomNumber,
                                                const vector<RealType>& atomPosX,
                                                const vector<RealType>& atomPosY,
                                                bool verbose = false ) {
  // Check vector sizes
  if ( atomNumber.size ( ) != atomPosX.size ( ) || atomNumber.size ( ) != atomPosY.size ( ) )
    throw aol::Exception ( "CalculateProjectedAtomicPotential: invalid parameters!", __FILE__, __LINE__ );
  
  // Read potential coefficients from file
  vector<aol::MultiVector<RealType>> coeff;
  GetPotentialCoefficients ( Param.specimenPotentialCoefficientsFile, atomNumber, coeff );
  
  // Calculate the projected potential
  const RealType e = 1.44;
  
  if ( verbose )
    cerr << "\tCalculating the projected potential... " << setw ( 3 ) << 0 << "%";
  
  for ( unsigned int j = 0; j < coeff.size ( ) ; j++ ) {
    for ( int y = 0; y < Param.Y ; y++ ) for ( int x = 0; x < Param.X ; x++ ) {
      const RealType posX = ( x - Param.X / 2 ) * Param.lenX / Param.X;
      const RealType posY = ( y - Param.Y / 2 ) * Param.lenY / Param.Y;
      
      RealType sum1 = 0;
      RealType sum2 = 0;
      
      RealType distSqr = aol::Sqr<RealType> ( posX - atomPosX[j] ) + aol::Sqr<RealType> ( posY - atomPosY[j] );
      
      if ( distSqr < 1e-5 )
        distSqr = 1e-5;
      
      for ( int i = 0; i < 3 ; i++ ) {
        sum1 += coeff[j][0][i] * boost::math::cyl_bessel_k<RealType, RealType> ( 0, 2 * aol::NumberTrait<RealType>::pi * sqrt ( distSqr ) * sqrt ( coeff[j][1][i] ) );
        sum2 += coeff[j][2][i] / coeff[j][3][i] * exp ( -aol::Sqr<RealType> ( aol::NumberTrait<RealType>::pi ) * distSqr / coeff[j][3][i] );
      }
      
      ProjectedPotential.set ( x, y, ProjectedPotential.get ( x, y ) + 2 * aol::Sqr<RealType> ( aol::NumberTrait<RealType>::pi ) * Constants::a0 * e * ( 2 * sum1 + sum2 ) );
    }
    
    if ( verbose )
      cerr << "\r\tCalculating the projected potential... " << setw ( 3 ) << static_cast<int> ( ( ( j + 1 ) / static_cast<double> ( coeff.size ( ) ) ) * 100 ) << "%";
  }
  
  if ( verbose )
    cerr << " done." << endl;
}

/**
 * \brief Simulates an exit wave in real space using the phase-object approximation
 */
template <typename RealType>
void SimulateExitWave ( Image<RealType>& ExitWave,
                        const Parameters<RealType>& Param,
                        const bool verbose = false ) {
  // Adjust image dimension parameters for the exit wave simulation
  Parameters<RealType> inputSimulationParam = Param;
  
  inputSimulationParam.X = Param.X + 2 * Param.inputDataBufferZone;
  inputSimulationParam.Y = Param.Y + 2 * Param.inputDataBufferZone;
  
  inputSimulationParam.lenX *= static_cast<RealType> ( inputSimulationParam.X ) / Param.X;
  inputSimulationParam.lenY *= static_cast<RealType> ( inputSimulationParam.Y ) / Param.Y;
  
  // Allocate memory and set all entries to zero
  ExitWave.reallocate ( inputSimulationParam.X, inputSimulationParam.Y );
  
  // Calculate the projected potential of the specimen
  switch ( Param.specimenFlag ) {
    case 0: {  // Lattice of alternating positive and negative point charges
      // Calculate the point charges' positions in nanometer
      const int numChargesX = 2 * Param.specimenPointsX;
      const int numChargesY = 2 * Param.specimenPointsY;
      
      vector<int> charge ( numChargesX * numChargesY );
      vector<RealType> chargeX ( numChargesX * numChargesY );
      vector<RealType> chargeY ( numChargesX * numChargesY );
      
      for ( int x = 0; x < numChargesX ; x++ ) for ( int y = 0; y < numChargesY ; y++ ) {
        charge[ y + numChargesY * x ] = 2 * ( ( y + x ) % 2 ) - 1;
        chargeX[ y + numChargesY * x ] = Param.specimenOffsetX + Param.lenX / Param.specimenPointsX * ( x - numChargesX / 2 + 0.5 );
        chargeY[ y + numChargesY * x ] = Param.specimenOffsetY + Param.lenY / Param.specimenPointsY * ( y - numChargesY / 2 + 0.5 );
      }
      
      CalculateProjectedPointChargePotential ( inputSimulationParam, ExitWave[0], charge, chargeX, chargeY, verbose );
      }
      break;
    case 1: {  // Single atom
      vector<int> atomNumber = { Param.specimenAtomicNumber };
      vector<RealType> atomPosX = { Param.specimenOffsetX };
      vector<RealType> atomPosY = { Param.specimenOffsetY };
      
      CalculateProjectedAtomicPotential ( inputSimulationParam, ExitWave[0], atomNumber, atomPosX, atomPosY, verbose );
      }
      break;
    case 2: {  // Standard lattice of atoms
      // Calculate the atoms' positions in nanometers and set the atomic numbers
      const int numAtomsX = 2 * Param.specimenPointsX;
      const int numAtomsY = 2 * Param.specimenPointsY;
      
      vector<int> atomNumber ( numAtomsX * numAtomsY, Param.specimenAtomicNumber );
      vector<RealType> atomPosX ( numAtomsX * numAtomsY, Param.specimenOffsetX );
      vector<RealType> atomPosY ( numAtomsX * numAtomsY, Param.specimenOffsetY );
      
      for ( int x = 0; x < numAtomsX ; x++ ) for ( int y = 0; y < numAtomsY ; y++ ) {
        const int index = y + x * numAtomsY;
        
        atomPosX[ index ] += Param.lenX / Param.specimenPointsX * ( x - numAtomsX / 2 + 0.5 );
        atomPosY[ index ] += Param.lenY / Param.specimenPointsY * ( y - numAtomsY / 2 + 0.5 );
      }
      
      CalculateProjectedAtomicPotential ( inputSimulationParam, ExitWave[0], atomNumber, atomPosX, atomPosY, verbose );
      
      }
      break;
    case 3: {  // Five atoms on a line (carbon, silicon, copper, gold, uranium)
      vector<int> atomNumber = { 6, 14, 29, 79, 92 };
      vector<RealType> atomPosX ( 5, Param.specimenOffsetX );
      vector<RealType> atomPosY ( 5, Param.specimenOffsetY );
      
      for ( int i = 0; i < 5 ; i++ )
        atomPosX[i] += ( i - 2 ) * Param.lenX / 5;
      
      CalculateProjectedAtomicPotential ( inputSimulationParam, ExitWave[0], atomNumber, atomPosX, atomPosY, verbose );
      }
      break;
    case 4: {  // Single atom near the image border
      vector<int> atomNumber = { Param.specimenAtomicNumber };
      vector<RealType> atomPosX = { 9 * Param.lenX / 20 };
      vector<RealType> atomPosY = { 0 };
      
      CalculateProjectedAtomicPotential ( inputSimulationParam, ExitWave[0], atomNumber, atomPosX, atomPosY, verbose );
      }
      break;
    case 5: {  // Atoms arranged in a honeycomb structure
      vector<RealType> atomPosX, atomPosY;
      
      const RealType cos30 = cos ( aol::NumberTrait<RealType>::pi / 6 );
      const RealType sin30 = sin ( aol::NumberTrait<RealType>::pi / 6 );
      
      const int nx = static_cast<int> ( ceil ( Param.lenX / ( cos30 * Param.specimenHoneycombParam ) ) );
      const int ny = static_cast<int> ( ceil ( Param.lenX / ( sin30 * Param.specimenHoneycombParam ) ) );
      
      RealType posY[2] = { Param.specimenOffsetY, Param.specimenOffsetY - ( 1 + sin30 ) * Param.specimenHoneycombParam };
      for ( int y = 0 ; y < ny ; y++ ) {
        RealType posX[2] = { Param.specimenOffsetX, Param.specimenOffsetX };
        posX[(y+1)%2] = Param.specimenOffsetX + cos30 * Param.specimenHoneycombParam;
        
        for ( int i = 0; i < 2 ; i++ ) {
          atomPosX.push_back ( posX[i] );
          atomPosY.push_back ( posY[i] );
        }
        
        for ( int x = 1; x < nx ; x++ ) {
          for ( int i = 0; i < 2 ; i++ ) {
            atomPosX.push_back ( posX[i] + x * cos30 * Param.specimenHoneycombParam );
            atomPosY.push_back ( posY[i] + (x % 2) * sin30 * Param.specimenHoneycombParam );
            atomPosX.push_back ( posX[i] - x * cos30 * Param.specimenHoneycombParam );
            atomPosY.push_back ( posY[i] + (x % 2) * sin30 * Param.specimenHoneycombParam );
          }
        }
        
        posY[0] += ( 1 + sin30 ) * Param.specimenHoneycombParam;
        posY[1] -= ( 1 + sin30 ) * Param.specimenHoneycombParam;
      }
      
      vector<int> atomNumber ( atomPosX.size ( ), Param.specimenAtomicNumber );
      
      CalculateProjectedAtomicPotential ( inputSimulationParam, ExitWave[0], atomNumber, atomPosX, atomPosY, verbose );
      
      }
      break;
    default:
      throw aol::Exception ( "Invalid specimen flag!", __FILE__, __LINE__ );
      break;
  }
  
  // Apply a bandwidth limit to the projected potential
  ExitWave.LowpassFilter ( Param );
  
  // Calculate the interaction parameter
  const RealType sigma = getInteractionParameter<RealType> ( 1000 * Param.AcceleratingVoltage );
  
  // Generate exit wave from the projected potential using the phase-object approximation
  ExitWave[0] *= sigma;
  
  ExitWave[1] = ExitWave[0];
  
  ExitWave[0].apply ( cos );
  ExitWave[1].apply ( sin );
  
  // Rotate the exit wave and apply a bandwidth limit
  ExitWave.Rotate ( Param.specimenRotation, true );
  ExitWave.LowpassFilter ( Param );
}

/**
 * \brief Simulates the focus image series that is used as input data to the algorithm
 * 
 * The images are simulated with the simulation mode given by Param.inputSimulationMode.
 */
template <typename RealType>
void SimulateInputImageSeries ( vector<Image<RealType>>& ImageSeries,
                                const aol::Vector<RealType>& Focus,
                                const qc::ScalarArray<RealType, qc::QC_2D>& Translation,
                                const Image<RealType>& RealSpaceExitWave,
                                const Parameters<RealType>& Param,
                                const bool verbose = false ) {
  if ( !RealSpaceExitWave.correctFormat ( Space::RealSpace, false ) )
    throw aol::Exception ( "Invalid exit wave format!", __FILE__, __LINE__ );
  
  if ( verbose )
    cerr << "\tSimulating the input image series... " << setw ( 3 ) << 0 << "%";
  
  // Calculate the size of the buffer zone that was already removed
  const int removedBufferX = ( Param.X + 2 * Param.inputDataBufferZone - RealSpaceExitWave.getNumX ( ) ) / 2;
  const int removedBufferY = ( Param.Y + 2 * Param.inputDataBufferZone - RealSpaceExitWave.getNumY ( ) ) / 2;
  
  // Adjust image dimension parameters for the input series simulation
  Parameters<RealType> inputSimulationParam = Param;
  
  inputSimulationParam.lenX *= static_cast<RealType> ( RealSpaceExitWave.getNumX ( ) ) / Param.X;
  inputSimulationParam.lenY *= static_cast<RealType> ( RealSpaceExitWave.getNumY ( ) ) / Param.Y;
  
  inputSimulationParam.X = RealSpaceExitWave.getNumX ( );
  inputSimulationParam.Y = RealSpaceExitWave.getNumY ( );
  
  // Prepare the OpenCL image simulation kernel
  KernelData<RealType> kData ( inputSimulationParam, true );
  
  // Calculate the Fourier space exit wave
  Image<RealType> ExitWave ( inputSimulationParam.X, inputSimulationParam.Y, Space::FourierSpace );
  RealSpaceExitWave.FourierTransformTo ( ExitWave );
  
  // Initialize the image series
  ImageSeries.clear ( );
  ImageSeries.resize ( Param.N, Image<RealType> ( Param.subsection_width, Param.subsection_height, Space::FourierSpace ) );
  
  int count = 0;
  const int nThreads = ( Param.inputSimulationMode == 1 ? static_cast<int> ( kData.getNumGPUs ( ) ) : Param.numThreads );
  #pragma omp parallel for num_threads ( nThreads )
  for ( int i = 0; i < Param.N ; i++ ) {
    Image<RealType> ImageWithBufferZone ( inputSimulationParam.X, inputSimulationParam.Y );
    Image<RealType> ImageWithBufferZoneAndSpecimenDrift ( inputSimulationParam.X, inputSimulationParam.Y );
    Image<RealType> RealSpaceImage ( Param.subsection_width, Param.subsection_height );
    
    // Simulate the real space TEM image
    SimulateImage ( ImageWithBufferZone, ExitWave, Focus[i], inputSimulationParam, kData, omp_get_thread_num ( ), true, false );
    
    // Check if the rectangular image subsection given by the four subsection_*
    // parameters is still completely contained within RealSpaceExitWave after
    // the image shift below
    const RealType TranslationInPixel[2] = { Translation.get ( i, 0 ) * Param.X / Param.lenX,
                                             Translation.get ( i, 1 ) * Param.Y / Param.lenY };
    
    if ( Param.subsection_x + Param.inputDataBufferZone - removedBufferX - TranslationInPixel[0] < 0 ||
         Param.subsection_x + Param.subsection_width + Param.inputDataBufferZone + removedBufferX - TranslationInPixel[1] > Param.X + 2 * Param.inputDataBufferZone ||
         Param.subsection_y + Param.inputDataBufferZone - removedBufferY - TranslationInPixel[1] < 0 ||
         Param.subsection_y + Param.subsection_height + Param.inputDataBufferZone + removedBufferY - TranslationInPixel[1] > Param.Y + 2 * Param.inputDataBufferZone )
      throw aol::Exception ( "The size of the buffer zone is too small for the given specimen drift! (Param: inputDataBufferZone, specimenDriftX, specimenDriftY)", __FILE__, __LINE__ );
    
    // Apply the specimen drift to the simulated image
    ImageWithBufferZone.Shift ( ImageWithBufferZoneAndSpecimenDrift, Translation.get ( i, 0 ), Translation.get ( i, 1 ), inputSimulationParam );
    
    // Remove the buffer zone
    ImageWithBufferZoneAndSpecimenDrift[0].copyBlockTo ( Param.subsection_x + Param.inputDataBufferZone - removedBufferX,
                                                         Param.subsection_y + Param.inputDataBufferZone - removedBufferY,
                                                         RealSpaceImage[0] );
    RealSpaceImage.FourierTransformTo ( ImageSeries[i] );
    
    if ( verbose ) {
      #pragma omp critical ( InputSeriesSimulationOutput )
      {
        ++count;
        cerr << "\r\tSimulating the input image series... " << setw ( 3 ) << static_cast<int> ( count / static_cast<double> ( Param.N ) * 100 ) << "%";
      }
    }
  }
  
  if ( verbose )
    cerr << " done." << endl;
}

#endif  // EWR_SIMULATIONBASE_H
