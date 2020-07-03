#ifndef EWR_FOCUSDERIVATIVEIMAGE_H
#define EWR_FOCUSDERIVATIVEIMAGE_H

#include "ReconstructionBase.h"

/**
 * \brief Calculates the autocorrelation of the exit wave weighted by the Focus derivative of the TCC using the focal integration approximation with periodic continuation
 * 
 * \param [out] FocusDerivativeImage The resulting image corresponding to the focus value given as the third parameter
 * 
 * \param [in] ExitWave The current estimate of the exit wave in Fourier space
 * 
 * \param [in] Focus The focus value
 * 
 * \param [in] Param Microscope parameters
 */
template <typename RealType>
void CalculateFocusDerivativeImageFocalInt ( Image<RealType>& FocusDerivativeImage,
                                             const Image<RealType>& ExitWave,
                                             const RealType Focus,
                                             const Parameters<RealType>& Param ) {
  // Intermediate results
  Image<RealType> tmp0 ( Param.X, Param.Y, Space::FourierSpace );
  Image<RealType> tmp00 ( Param.X, Param.Y, Space::FourierSpace );
  Image<RealType> tmp1 ( Param.X, Param.Y, Space::FourierSpace );
  Image<RealType> tmp11 ( Param.X, Param.Y, Space::FourierSpace );
  Image<RealType> tmp2 ( Param.X, Param.Y, Space::RealSpace );
  Image<RealType> tmp22 ( Param.X, Param.Y, Space::RealSpace );
  Image<RealType> tmp3 ( Param.X, Param.Y, Space::RealSpace );
  
  // Calculate tmp0 = Psi * t_{j,k} and tmp00 = Psi * t_{j,k} * g partially, where
  //   g = (-Pi*i*lambda) * abs(freq)^2 - 2*(Pi*alpha)^2 * (Z*abs(freq)^2 + C_s*lambda^2*abs(freq)^4).
  // (we calculate everything except the temporal coherence factor, i.e. everything
  // that does not depend on k)
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
    
    const RealType FreqAbsSqr = freqX * freqX + freqY * freqY;
    
    const RealType derivativeCoeffRealPart = -2 * aol::Sqr<RealType> ( aol::NumberTrait<RealType>::pi * Param.alpha )
                                            * ( Focus * FreqAbsSqr + Param.SphericalAberration * Param.lambda * Param.lambda * FreqAbsSqr * FreqAbsSqr );
    const RealType derivativeCoeffImagPart = -aol::NumberTrait<RealType>::pi * Param.lambda * FreqAbsSqr;
    
    tmp00[0].set ( x, y, derivativeCoeffRealPart * tmp0[0].get ( x, y ) - derivativeCoeffImagPart * tmp0[1].get ( x, y ) );
    tmp00[1].set ( x, y, derivativeCoeffRealPart * tmp0[1].get ( x, y ) + derivativeCoeffImagPart * tmp0[0].get ( x, y ) );
  }
  
  // Calculate the summands
  for ( int k = -Param.FocalIntegration_M; k <= Param.FocalIntegration_M ; k++ ) {
    // 1. step: finish the calculation of Psi * t_{j,k} and Psi * t_{j,k} * g
    for ( int y = 0; y < Param.Y ; ++y ) for ( int x = 0; x < Param.X ; ++x ) {
      const RealType freqX = invLenX * ( x - Param.X / 2 );
      const RealType freqY = invLenY * ( y - Param.Y / 2 );
      
      const RealType temporalCoherenceExp = -aol::NumberTrait<RealType>::pi * Param.lambda * ( k * Param.FocalIntegration_delta ) * ( freqX * freqX + freqY * freqY );
      
      const RealType cos_temporalCoherence = cos ( temporalCoherenceExp );
      const RealType sin_temporalCoherence = sin ( temporalCoherenceExp );
      
      tmp1[0].set ( x, y, tmp0[0].get ( x, y ) * cos_temporalCoherence - tmp0[1].get ( x, y ) * sin_temporalCoherence );
      tmp1[1].set ( x, y, tmp0[1].get ( x, y ) * cos_temporalCoherence + tmp0[0].get ( x, y ) * sin_temporalCoherence );
      
      tmp11[0].set ( x, y, tmp00[0].get ( x, y ) * cos_temporalCoherence - tmp00[1].get ( x, y ) * sin_temporalCoherence );
      tmp11[1].set ( x, y, tmp00[1].get ( x, y ) * cos_temporalCoherence + tmp00[0].get ( x, y ) * sin_temporalCoherence );
    }
    
    // 2. step: calculate the inverse Fourier transforms
    tmp1.FourierTransformTo ( tmp2 );
    tmp11.FourierTransformTo ( tmp22 );
    
    // 3. step: calculate 2 * Re ( conjugate(tmp2) * tmp22 )
    for ( int y = 0; y < Param.Y ; ++y ) for ( int x = 0; x < Param.X ; ++x )
      tmp1[0].set ( x, y, 2 * ( tmp2[0].get ( x, y ) * tmp22[0].get ( x, y ) + tmp2[1].get ( x, y ) * tmp22[1].get ( x, y ) ) );
    
    // 4. step: calculate the oefficient c_k
    RealType c_k = Param.FocalIntegration_delta / ( sqrt ( 2 * aol::NumberTrait<RealType>::pi ) * Param.FocalSpread )
                   * exp ( -aol::Sqr<RealType> ( k * Param.FocalIntegration_delta ) / ( 2 * aol::Sqr<RealType> ( Param.FocalSpread ) ) );
    
    // 5. step: add up the partial images
    tmp3[0].addMultiple ( tmp1[0], c_k );
  }
  
  // Transform the result to Fourier space
  tmp3.FourierTransformTo ( FocusDerivativeImage );
}

/**
 * \brief Calculates the autocorrelation of ExitWave weighted by the Focus derivative of the TCC using the full TCC
 * 
 * \param [out] FocusDerivativeImage The resulting image corresponding to the focus value given as the third parameter
 * 
 * \param [in] ExitWave The current estimate of the exit wave in Fourier space
 * 
 * \param [in] Focus The focus value
 * 
 * \param [in] Param Microscope parameters
 * 
 * \param [in] kData OpenCL kernel data for the computation on the GPU
 * 
 * \param [in] gpu_id Index of the GPU that is used for the computation
 */
template <typename RealType>
void CalculateFocusDerivativeImageFullInt ( Image<RealType>& FocusDerivativeImage,
                                            const Image<RealType>& ExitWave,
                                            const RealType Focus,
                                            const Parameters<RealType>& Param,
                                            KernelData<RealType>& kData,
                                            const int gpu_id ) {
  SimulateImageFocusDerivative_KernelData<RealType>* skData = kData.SimulateImageFocusDerivativeData[gpu_id];
  
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
  try {
    skData->ocl.queue.finish ( );
    skData->ocl.queue.enqueueNDRangeKernel ( skData->ocl.kernel, cl::NullRange, cl::NDRange ( Param.X, Param.Y ), cl::NullRange );
    skData->ocl.queue.finish ( );
  } catch ( cl::Error *err ) {
    throw aol::Exception ( err->what ( ), __FILE__, __LINE__ );
  }
  
  // Read the results
  try {
    if ( skData->useKernelV1 )
      skData->ocl.queue.enqueueReadImage ( skData->imSimulatedImage, CL_TRUE, origin, region, 0, 0, skData->interleaved_output.data ( ) );
    else
      skData->ocl.queue.enqueueReadBuffer ( skData->bufSimulatedImage, CL_TRUE, 0, 2 * Param.X * Param.Y * sizeof ( RealType ), skData->interleaved_output.data ( ) );
  } catch ( cl::Error *err ) {
    throw aol::Exception ( err->what ( ), __FILE__, __LINE__ );
  }
  
  for ( int i = 0; i < Param.X * Param.Y ; i++ ) {
    FocusDerivativeImage[0][i] = skData->interleaved_output[ 2 * i + 0 ];
    FocusDerivativeImage[1][i] = skData->interleaved_output[ 2 * i + 1 ];
  }
}

/**
 * \brief Calculates the autocorrelation of ExitWave weighted by the Focus derivative of the TCC
 * 
 * \param [out] FocusDerivativeImage The resulting image corresponding to the focus value given as the third parameter
 * 
 * \param [in] ExitWave The current estimate of the exit wave in Fourier space
 * 
 * \param [in] Focus The focus value
 * 
 * \param [in] Param Microscope parameters
 * 
 * \param [in] kData OpenCL kernel data  (accessed only if the computation is to be performed on the GPU, i.e. if Param.simulationMode == 1)
 * 
 * \param [in] gpu_id Index of the GPU that is used for the computation (if Param.simulationMode == 1)
 */
template <typename RealType>
void CalculateFocusDerivativeImage ( Image<RealType>& FocusDerivativeImage,
                                     const Image<RealType>& ExitWave,
                                     const RealType Focus,
                                     const Parameters<RealType>& Param,
                                     KernelData<RealType>& kData,
                                     const int gpu_id ) {
  // Check the format of FocusDerivativeImage and ExitWave
  if ( !FocusDerivativeImage.correctFormat ( Space::FourierSpace, true ) )
    throw aol::Exception ( "Invalid format of FocusDerivativeImage!", __FILE__, __LINE__ );
  
  if ( !ExitWave.correctFormat ( Space::FourierSpace, true ) )
    throw aol::Exception ( "Invalid exit wave format!", __FILE__, __LINE__ );
  
  switch ( Param.simulationMode ) {
    case 0:
      // focal integration approximation with periodic continuation (CPU)
      CalculateFocusDerivativeImageFocalInt ( FocusDerivativeImage, ExitWave, Focus, Param );
      break;
    case 1:
      // full computation of the weighted autocorrelation with zero continuation (GPU)
      CalculateFocusDerivativeImageFullInt ( FocusDerivativeImage, ExitWave, Focus, Param, kData, gpu_id );
      break;
    default:
      throw aol::Exception ( "Invalid simulation mode!", __FILE__, __LINE__ );
      break;
  }
}
        

#endif  // EWR_FOCUSDERIVATIVEIMAGE_H
