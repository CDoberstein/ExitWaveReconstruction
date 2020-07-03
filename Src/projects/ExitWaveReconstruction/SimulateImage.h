#ifndef EWR_SIMULATION_H
#define EWR_SIMULATION_H

#include "OpenCLKernelData.h"
/// \cond
#include "Parameters.h"
#include "Image.h"
#include "Physics.h"

#include <thread>
#include <ctime>
#include <ratio>
#include <chrono>
/// \endcond

/**
 * \brief Calculates the pure phase transfer function multiplied by the aperture function
 * 
 * \param [out] PhaseTransfer the resulting phase transfer function as a complex-valued image
 * 
 * \param [in] Param the microscope parameters
 * 
 * \param [in] Focus the Focus parameter (required for the evaluation of the aberration function)
 */
void CalculatePhaseTransferFunction ( qc::MultiArray<RealType, 2, 2>& PhaseTransfer,
                                      const Parameters<RealType>& Param,
                                      const RealType Focus ) {
  if ( PhaseTransfer.getNumX ( ) != Param.X || PhaseTransfer.getNumY ( ) != Param.Y )
    throw aol::Exception ( "Invalid size of the PhaseTransfer MultiArray!", __FILE__, __LINE__ );
  
  const RealType invLenX = 1 / Param.lenX;
  const RealType invLenY = 1 / Param.lenY;
  
  const RealType Xhalf = static_cast<RealType> ( Param.X ) / 2;
  const RealType Yhalf = static_cast<RealType> ( Param.Y ) / 2;
  
  for ( int y = 0; y < Param.Y ; y++ ) for ( int x = 0; x < Param.X ; x++ ) {
    const RealType freqX = invLenX * ( x - Xhalf );
    const RealType freqY = invLenY * ( y - Yhalf );
    
    if ( ApertureFunction ( freqX, freqY, Param ) != 0 ) {
      const RealType phaseTransferExp = -2 * aol::NumberTrait<RealType>::pi * AberrationFunction ( freqX, freqY, Focus, Param );
      
      PhaseTransfer[0].set ( x, y, cos ( phaseTransferExp ) );
      PhaseTransfer[1].set ( x, y, sin ( phaseTransferExp ) );
    } else {
      PhaseTransfer[0].set ( x, y, 0 );
      PhaseTransfer[1].set ( x, y, 0 );
    }
  }
}

/**
 * \brief Simulates a TEM image using the focal integration approximation and periodic continuation
 * 
 * Simulates a TEM image in Fourier space or real space from a given exit wave and
 * transfer function, where the transfer function is given by the focal integration
 * approximation of Ishizuka's TCC. The input exit wave is expected to be in Fourier
 * space coordinates with the lowest frequency located at the image center.
 * 
 * The space of the resulting simulated image is determined by the image space of Img
 * when calling this function.
 */
template <typename RealType>
void SimulateImageFocalInt ( Image<RealType>& Img,
                             const Image<RealType>& ExitWave,
                             const RealType Focus,
                             const Parameters<RealType>& Param,
                             const bool verbose = false ) {
  if ( verbose )
    cerr << "Simulating a TEM image with the focal integration approximation" << endl
         << "and periodic continuation... ";
  
  // Reset the result
  Img[0].setZero ( );
  Img[1].setZero ( );
  
  // Intermediate results
  qc::MultiArray<RealType, 2, 2> tmp0 ( Param.X, Param.Y );
  qc::MultiArray<RealType, 2, 2> tmp1 ( Param.X, Param.Y );
  qc::MultiArray<RealType, 2, 2> tmp2 ( Param.X, Param.Y );
  qc::MultiArray<RealType, 2, 2> tmp3 ( Param.X, Param.Y );
  
  // Calculate Psi * t_{j,k} partially (everything except the temporal coherence factor, i.e. everything that does not depend on k)
  const RealType invLenX = 1 / Param.lenX;
  const RealType invLenY = 1 / Param.lenY;
  
  for ( int y = 0; y < Param.Y ; ++y ) for ( int x = 0; x < Param.X ; ++x ) {
    const RealType freqX = invLenX * ( x - Param.X / 2 );
    const RealType freqY = invLenY * ( y - Param.Y / 2 );
    
    if ( ApertureFunction ( freqX, freqY, Param ) == 0 )
      continue;
    
    const RealType phaseTransferExp = -2 * aol::NumberTrait<RealType>::pi * AberrationFunction ( freqX, freqY, Focus, Param );
    
    const RealType sin_phaseTransfer = sin ( phaseTransferExp );
    const RealType cos_phaseTransfer = cos ( phaseTransferExp );
    
    const RealType spatialCoherence = SpatialCoherenceEnvelope ( freqX, freqY, Focus, Param );
    
    tmp0[0].set ( x, y, spatialCoherence * ( ExitWave[0].get ( x, y ) * cos_phaseTransfer - ExitWave[1].get ( x, y ) * sin_phaseTransfer ) );
    tmp0[1].set ( x, y, spatialCoherence * ( ExitWave[0].get ( x, y ) * sin_phaseTransfer + ExitWave[1].get ( x, y ) * cos_phaseTransfer ) );
  }
  
  // Calculate the summands
  for ( int k = -Param.FocalIntegration_M; k <= Param.FocalIntegration_M ; k++ ) {
    // 1. step: calculate Psi * t_{j,k}
    for ( int y = 0; y < Param.Y ; ++y ) for ( int x = 0; x < Param.X ; ++x ) {
      const RealType freqX = invLenX * ( x - Param.X / 2 );
      const RealType freqY = invLenY * ( y - Param.Y / 2 );
      
      const RealType temporalCoherenceExp = -aol::NumberTrait<RealType>::pi * Param.lambda * ( k * Param.FocalIntegration_delta ) * ( freqX * freqX + freqY * freqY );
      
      const RealType cos_temporalCoherence = cos ( temporalCoherenceExp );
      const RealType sin_temporalCoherence = sin ( temporalCoherenceExp );
      
      tmp1[0].set ( x, y, tmp0[0].get ( x, y ) * cos_temporalCoherence - tmp0[1].get ( x, y ) * sin_temporalCoherence );
      tmp1[1].set ( x, y, tmp0[1].get ( x, y ) * cos_temporalCoherence + tmp0[0].get ( x, y ) * sin_temporalCoherence );
    }
    
    // 2. step: calculate the inverse Fourier transform of Psi * t_{j,k}
    /* Note: Theoretically we should perform a frequency shift of *tmp1
     *       here, because the lowest frequency of *tmp1 is located in the
     *       center. However, this only changes the entries of *tmp2 after
     *       the following inverse Fourier transform by a multiplication
     *       with (-1)^(r+c), where r and c are the row and column numbers
     *       of the entries. Therefore a fourier shift has no effect
     *       at all after the subsequent calculation of the pointwise
     *       squared absolute value.
     */
    qc::FourierTransform ( tmp1, tmp2, qc::FTBackward );
    
    // 3. step: calculate the (pointwise) squared absolute value
    for ( int y = 0; y < Param.Y ; ++y ) for ( int x = 0; x < Param.X ; ++x )
      tmp1[0].set ( x, y, aol::Sqr<RealType> ( tmp2[0].get ( x, y ) ) + aol::Sqr<RealType> ( tmp2[1].get ( x, y ) ) );
    
    // 4. step: calculate the oefficient c_k
    RealType c_k = Param.FocalIntegration_delta / ( sqrt ( 2 * aol::NumberTrait<RealType>::pi ) * Param.FocalSpread )
                   * exp ( -aol::Sqr<RealType> ( k * Param.FocalIntegration_delta ) / ( 2 * aol::Sqr<RealType> ( Param.FocalSpread ) ) );
    
    // 5. step: add up the partial images
    if ( Img.getImageSpace ( ) == Space::FourierSpace )
      tmp3[0].addMultiple ( tmp1[0], c_k );
    else
      Img[0].addMultiple ( tmp1[0], c_k );
    
    if ( verbose )
      cerr << "\rand periodic continuation... " << setw ( 3 ) << static_cast<int> ( static_cast<RealType> ( k + 1 + Param.FocalIntegration_M ) / ( 2 * Param.FocalIntegration_M + 1 ) * 100 ) << "%";
  }
  
  // Apply the forward Fourier transform if necessary
  if ( Img.getImageSpace ( ) == Space::FourierSpace ) {
    Image<RealType> tmp3_flatcopy ( tmp3, Space::RealSpace, false, aol::FLAT_COPY );
    tmp3_flatcopy.FourierTransformTo ( Img );
  }
}

/**
 * \brief Simulates a TEM image using the full TCC on the GPU and zero continuation
 * 
 * Simulates a TEM image in Fourier space or real space from a given exit wave using
 * Ishizuka's TCC. The input exit wave is expected to be in Fourier space coordinates
 * with the lowest frequency located at the image center.
 * 
 * The space of the resulting simulated image is determined by the image space of Img
 * when calling this function.
 * 
 * \param [in] kData OpenCL kernel data for the computation on the GPU
 * 
 * \param [in] gpu_id index of the GPU that is used for the computation
 */
void SimulateImageFullInt ( Image<RealType>& Img,
                            const Image<RealType>& ExitWave,
                            const RealType Focus,
                            const Parameters<RealType>& Param,
                            KernelData<RealType>& kData,
                            const int gpu_id,
                            const bool verbose = false ) {
  chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now ( );
  
  SimulateImage_KernelData<RealType>* skData = kData.SimulateImageData[gpu_id];
  
  // Check kernel data parameters
  if ( skData->Param.X != Param.X || skData->Param.Y != Param.Y ||
       skData->Param.lenX != Param.lenX || skData->Param.lenY != Param.lenY )
    throw aol::Exception ( "Invalid parameters!", __FILE__, __LINE__ );
  
  // Calculate phase transfer function
  CalculatePhaseTransferFunction ( skData->PhaseTransfer, Param, Focus );
  
  // Set input arguments
  for ( int i = 0; i < Param.X * Param.Y ; i++ ) {
    skData->interleaved_input[ 4 * i + 0 ] = ExitWave[0][i];
    skData->interleaved_input[ 4 * i + 1 ] = ExitWave[1][i];
    skData->interleaved_input[ 4 * i + 2 ] = skData->PhaseTransfer[0][i];
    skData->interleaved_input[ 4 * i + 3 ] = skData->PhaseTransfer[1][i];
  }
  
  cl::size_t<3> origin;
  cl::size_t<3> region;
  
  origin[0] = origin[1] = origin[2] = 0;
  
  region[0] = Param.X;
  region[1] = Param.Y;
  region[2] = 1;
  
  try {
    if ( skData->useKernelV1 ) {
      skData->ocl.queue.enqueueWriteImage ( skData->imInput, CL_TRUE, origin, region, 0, 0, skData->interleaved_input.data ( ) );
      skData->ocl.kernel.setArg ( 2, Focus );
    } else {
      skData->ocl.queue.enqueueWriteBuffer ( skData->bufInput, CL_TRUE, 0, 4 * Param.X * Param.Y * sizeof ( RealType ), skData->interleaved_input.data ( ) );
      skData->ocl.kernel.setArg ( 2, Focus );
    }
  } catch ( cl::Error *err ) {
    throw aol::Exception ( err->what ( ), __FILE__, __LINE__ );
  }
  
  // Perform image simulation on the GPU
  chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now ( );
  
  try {
    skData->ocl.queue.finish ( );
    skData->ocl.queue.enqueueNDRangeKernel ( skData->ocl.kernel, cl::NullRange, cl::NDRange ( Param.X, Param.Y ), cl::NullRange );
    skData->ocl.queue.finish ( );
  } catch ( cl::Error *err ) {
    throw aol::Exception ( err->what ( ), __FILE__, __LINE__ );
  }
  
  chrono::high_resolution_clock::time_point t3 = chrono::high_resolution_clock::now ( );
  
  // Read the results
  try {
    if ( skData->useKernelV1 )
      skData->ocl.queue.enqueueReadImage ( skData->imSimulatedImage, CL_TRUE, origin, region, 0, 0, skData->interleaved_output.data ( ) );
    else
      skData->ocl.queue.enqueueReadBuffer ( skData->bufSimulatedImage, CL_TRUE, 0, 2 * Param.X * Param.Y * sizeof ( RealType ), skData->interleaved_output.data ( ) );
  } catch ( cl::Error *err ) {
    throw aol::Exception ( err->what ( ), __FILE__, __LINE__ );
  }
  
  // Transform image to real space if necessary
  if ( Img.getImageSpace ( ) == Space::RealSpace ) {
    Image<RealType> tmp ( Img.getNumX ( ), Img.getNumY ( ), Space::FourierSpace );
    for ( int i = 0; i < Param.X * Param.Y ; i++ ) {
      tmp[0][i] = skData->interleaved_output[ 2 * i + 0 ];
      tmp[1][i] = skData->interleaved_output[ 2 * i + 1 ];
    }
    tmp.FourierTransformTo ( Img );
  } else {
    for ( int i = 0; i < Param.X * Param.Y ; i++ ) {
      Img[0][i] = skData->interleaved_output[ 2 * i + 0 ];
      Img[1][i] = skData->interleaved_output[ 2 * i + 1 ];
    }
  }
  
  // Print computation time if verbose is true
  if ( verbose ) {
    chrono::duration<RealType> time_span12 = chrono::duration_cast<chrono::duration<RealType>> ( t2 - t1 );
    chrono::duration<RealType> time_span23 = chrono::duration_cast<chrono::duration<RealType>> ( t3 - t2 );
    
    cerr << "Simulating the image on the GPU using the full TCC took " << time_span23.count ( ) << " seconds. (Preparation time: "
         << time_span12.count ( ) << " seconds)" << endl;
  }
}

/**
 * \brief Simulates a TEM image in Fourier space or real space from a given exit wave and focus value
 * 
 * This function calculates the weighted autocorrelation of the exit wave, where the
 * weight is determined by the simulationMode parameter (if inputDataSimulation == false)
 * or the inputSimulationMode parameter (if inputDataSimulation == true ).
 * 
 * The space of the resulting simulated image is determined by the image space of Img
 * when calling this function.
 * 
 * \param [in] kData OpenCL kernel data  (accessed only if the computation is to be performed on the GPU, i.e. if Param.simulationMode == 1)
 * 
 * \param [in] gpu_id index of the GPU that is used for the computation (if Param.simulationMode == 1)
 */
template <typename RealType>
void SimulateImage ( Image<RealType>& Img,
                     const Image<RealType>& ExitWave,
                     const RealType Focus,
                     const Parameters<RealType>& Param,
                     KernelData<RealType>& kData,
                     const int gpu_id,
                     const bool inputDataSimulation = false,
                     const bool verbose = false ) {
  // Check the exit wave's size and format
  if ( !(ExitWave.getNumX ( ) == Param.X && ExitWave.getNumY ( ) == Param.Y &&
         ExitWave.correctFormat ( Space::FourierSpace, true ) ) )
    throw aol::Exception ( "Invalid exit wave size or format!", __FILE__, __LINE__ );
  
  // Check the image's size and format
  if ( !(Img.getNumX ( ) == Param.X && Img.getNumY ( ) == Param.Y &&
         ( Img.correctFormat ( Space::FourierSpace, true ) || Img.correctFormat ( Space::RealSpace, false ) ) ) )
    throw aol::Exception ( "Invalid image size or format!", __FILE__, __LINE__ );
  
  const int simulationMode = ( inputDataSimulation ? Param.inputSimulationMode : Param.simulationMode );
  switch ( simulationMode ) {
    case 0:
      // focal integration approximation with periodic continuation (CPU)
      SimulateImageFocalInt ( Img, ExitWave, Focus, Param, verbose );
      break;
    case 1:
      // full computation of the weighted autocorrelation with zero continuation (GPU)
      SimulateImageFullInt ( Img, ExitWave, Focus, Param, kData, gpu_id, verbose );
      break;
    default:
      throw aol::Exception ( "Invalid simulation mode!", __FILE__, __LINE__ );
      break;
  }
}

#endif  // EWR_SIMULATION_H
