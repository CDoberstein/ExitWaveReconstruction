#ifndef EWR_EXITWAVEDERIVATIVE_H
#define EWR_EXITWAVEDERIVATIVE_H

#include "ReconstructionBase.h"

/**
 * \brief Calculates one summand of the exit wave derivative with the focal integration approximation and periodic continuation
 * 
 * \param [out] Summand The resulting summand corresponding to the given Focus value
 * 
 * \param [in] ExitWave The current estimate of the exit wave in Fourier space
 * 
 * \param [in] DiffImage The inner part of the functional's data term (corresponding to the given Focus value!)
 * 
 * \param [in] Focus The Focus value
 * 
 * \param [in] Param Microscope parameters
 */
template <typename RealType>
void CalculateExitWaveDerivativeSummandFocalInt ( Image<RealType>& Summand,
                                                  const Image<RealType>& ExitWave,
                                                  const Image<RealType>& DiffImage,
                                                  const RealType Focus,
                                                  const Parameters<RealType>& Param ) {
  // Initialize Summand
  Summand.setZero ( );
  
  // Transform DiffImage to real space
  Image<RealType> RealSpaceDiffImage ( Param.X, Param.Y );
  DiffImage.FourierTransformTo ( RealSpaceDiffImage );
  
  // Intermediate results
  Image<RealType> tmp0 ( Param.X, Param.Y, Space::FourierSpace );
  Image<RealType> tmp1 ( Param.X, Param.Y, Space::FourierSpace );
  Image<RealType> tmp2 ( Param.X, Param.Y, Space::RealSpace );
  
  // Calculate t_{j,k} partially (everything except the temporal coherence factor, i.e. everything that does not depend on k)
  Image<RealType> t_jk ( Param.X, Param.Y, Space::FourierSpace );
  
  const RealType invLenX = 1 / Param.lenX;
  const RealType invLenY = 1 / Param.lenY;
  
  for ( int y = 0; y < Param.Y ; ++y ) for ( int x = 0; x < Param.X ; ++x ) {
    const RealType freqX = invLenX * ( x - Param.X / 2 );
    const RealType freqY = invLenY * ( y - Param.Y / 2 );
    
    if ( ApertureFunction ( freqX, freqY, Param ) == 0 ) {
      t_jk[0].set ( x, y, 0 );
      t_jk[1].set ( x, y, 0 );
      continue;
    }
    
    const RealType phaseTransferExp = -2 * aol::NumberTrait<RealType>::pi * AberrationFunction ( freqX, freqY, Focus, Param );
    const RealType spatialCoherence = SpatialCoherenceEnvelope ( freqX, freqY, Focus, Param );
    
    t_jk[0].set ( x, y, spatialCoherence * cos ( phaseTransferExp ) );
    t_jk[1].set ( x, y, spatialCoherence * sin ( phaseTransferExp ) );
  }
  
  // Sum over all k from the focal integration approximation
  for ( int k = -Param.FocalIntegration_M; k <= Param.FocalIntegration_M ; k++ ) {
    // 1. step: finish the calculation of t_{j,k}
    for ( int y = 0; y < Param.Y ; ++y ) for ( int x = 0; x < Param.X ; ++x ) {
      const RealType freqX = invLenX * ( x - Param.X / 2 );
      const RealType freqY = invLenY * ( y - Param.Y / 2 );
      
      const RealType temporalCoherenceExp = -aol::NumberTrait<RealType>::pi * Param.lambda * ( k * Param.FocalIntegration_delta ) * ( freqX * freqX + freqY * freqY );
      
      const RealType cos_temporalCoherence = cos ( temporalCoherenceExp );
      const RealType sin_temporalCoherence = sin ( temporalCoherenceExp );
      
      tmp0[0].set ( x, y, t_jk[0].get ( x, y ) * cos_temporalCoherence - t_jk[1].get ( x, y ) * sin_temporalCoherence );
      tmp0[1].set ( x, y, t_jk[1].get ( x, y ) * cos_temporalCoherence + t_jk[0].get ( x, y ) * sin_temporalCoherence );
    }
    
    // 2. step: calculate ExitWave * t_{j,k}
    for ( int y = 0; y < Param.Y ; y++ ) for ( int x = 0; x < Param.X ; x++ ) {
      tmp1[0].set ( x, y, ExitWave[0].get ( x, y ) * tmp0[0].get ( x, y ) - ExitWave[1].get ( x, y ) * tmp0[1].get ( x, y ) );
      tmp1[1].set ( x, y, ExitWave[0].get ( x, y ) * tmp0[1].get ( x, y ) + ExitWave[1].get ( x, y ) * tmp0[0].get ( x, y ) );
    }
    
    // 3. step: calculate the inverse Fourier transform of ExitWave * t_{j,k}
    tmp1.FourierTransformTo ( tmp2 );
    
    // 4. step: multipy the result with RealSpaceDiffImage
    // Note: it follows from the definition of DiffImage that RealSpaceDiffImage
    //       is real-valued
    tmp2[0] *= RealSpaceDiffImage[0];
    tmp2[1] *= RealSpaceDiffImage[0];
    
    // 5. step: calculate the Fourier transform of the product
    tmp2.FourierTransformTo ( tmp1 );
    
    // 6. step: Multipy the result with c_k * conjugate(t_{j,k})
    RealType c_k = Param.FocalIntegration_delta / ( sqrt ( 2 * aol::NumberTrait<RealType>::pi ) * Param.FocalSpread )
                  * exp ( - ( k * Param.FocalIntegration_delta ) * ( k * Param.FocalIntegration_delta ) /
                          ( 2 * Param.FocalSpread * Param.FocalSpread ) );
    
    for ( int y = 0; y < Param.Y ; y++ ) for ( int x = 0; x < Param.X ; x++ ) {
      Summand[0].add ( x, y, c_k * ( tmp1[0].get ( x, y ) * tmp0[0].get ( x, y ) + tmp1[1].get ( x, y ) * tmp0[1].get ( x, y ) ) );
      Summand[1].add ( x, y, c_k * ( tmp1[1].get ( x, y ) * tmp0[0].get ( x, y ) - tmp1[0].get ( x, y ) * tmp0[1].get ( x, y ) ) );
    }
  }
}

/**
 * \brief Calculates one summand of the exit wave derivative with the full TCC on the GPU
 * 
 * \param [out] Summand The resulting summand corresponding to the given Focus value
 * 
 * \param [in] ExitWave The current estimate of the exit wave in Fourier space
 * 
 * \param [in] DiffImage The inner part of the functional's data term (corresponding to the given Focus value!)
 * 
 * \param [in] Focus The Focus value
 * 
 * \param [in] Param Microscope parameters
 * 
 * \param [in] kData OpenCL kernel data for the computation on the gpu
 * 
 * \param [in] gpu_id Index of the GPU that is used for the computation
 */
template <typename RealType>
void CalculateExitWaveDerivativeSummandFullInt ( Image<RealType>& Summand,
                                                 const Image<RealType>& ExitWave,
                                                 const Image<RealType>& DiffImage,
                                                 const RealType Focus,
                                                 const Parameters<RealType>& Param,
                                                 KernelData<RealType>& kData,
                                                 const int gpu_id ) {
  ExitWaveDerivative_KernelData<RealType>* skData = kData.ExitWaveDerivativeData[gpu_id];
  
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
  
  for ( int i = 0; i < Param.X * Param.Y ; i++ ) {
    skData->interleaved_DiffImage[ 2 * i + 0 ] = DiffImage[0][i];
    skData->interleaved_DiffImage[ 2 * i + 1 ] = DiffImage[1][i];
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
      skData->ocl.queue.enqueueWriteImage ( skData->imDiffImage, CL_TRUE, origin, region, 0, 0, skData->interleaved_DiffImage.data ( ) );
      skData->ocl.kernel.setArg ( 3, Focus );
    } else {
      skData->ocl.queue.enqueueWriteBuffer ( skData->bufInput, CL_TRUE, 0, 4 * Param.X * Param.Y * sizeof ( RealType ), skData->interleaved_input.data ( ) );
      skData->ocl.queue.enqueueWriteBuffer ( skData->bufDiffImage, CL_TRUE, 0, 2 * Param.X * Param.Y * sizeof ( RealType ), skData->interleaved_DiffImage.data ( ) );
      skData->ocl.kernel.setArg ( 3, Focus );
    }
  } catch ( cl::Error *err ) {
    throw aol::Exception ( err->what ( ), __FILE__, __LINE__ );
  }
  
  // Execute kernel
  try {
    skData->ocl.queue.finish ( );
    skData->ocl.queue.enqueueNDRangeKernel ( skData->ocl.kernel, cl::NullRange, cl::NDRange ( Param.X, Param.Y ), cl::NullRange );
    skData->ocl.queue.finish ( );
  } catch ( cl::Error *err ) {
    throw aol::Exception ( err->what ( ), __FILE__, __LINE__ );
  }
  
  // Read results
  try {
    if ( skData->useKernelV1 )
      skData->ocl.queue.enqueueReadImage ( skData->imDerivativePart, CL_TRUE, origin, region, 0, 0, skData->interleaved_DerivativePart.data ( ) );
    else
      skData->ocl.queue.enqueueReadBuffer ( skData->bufDerivativePart, CL_TRUE, 0, 2 * Param.X * Param.Y * sizeof ( RealType ), skData->interleaved_DerivativePart.data ( ) );
  } catch ( cl::Error *err ) {
    throw aol::Exception ( err->what ( ), __FILE__, __LINE__ );
  }
  
  // Write the interleaved data to Summand (TODO: note the minus sign in the derivative formula)
  for ( int i = 0; i < Param.X * Param.Y ; i++ ) {
    Summand[0][i] = skData->interleaved_DerivativePart[ 2 * i + 0 ];
    Summand[1][i] = -skData->interleaved_DerivativePart[ 2 * i + 1 ];
  }
}

/**
 * \brief Calculates one summand of the exit wave derivative
 * 
 * This function calculates the L2 scalar product of DiffImage and WCC(x,y)
 * for all pixels (x,y) of the exit wave and stores the results in the
 * corresponding pixels of Summand. Here, WCC(x,y) is the cross correlation
 * of ExitWave and delta_(x,y) (resp. i*delta_(x,y) for the imaginary part
 * of the summand) weighted by the TCC given by Param.simulationMode, where
 * delta_(x,y) is the image that is 1 at (x,y) and zero everywhere else.
 * 
 * \warning The focus value given as the fourth parameter must be equal to
 *          the focus value used for the simulation of DiffImage in order
 *          for this function to calculate a correct summand of the exit
 *          wave derivative.
 * 
 * \param [out] Summand The resulting summand corresponding to the given Focus value
 * 
 * \param [in] ExitWave The current estimate of the exit wave in Fourier space
 * 
 * \param [in] DiffImage The inner part of the functional's data term (corresponding to the given Focus value!)
 * 
 * \param [in] Focus The Focus value
 * 
 * \param [in] Param Microscope parameters
 * 
 * \param [in] kData OpenCL kernel data (accessed only if the computation is to be performed on the GPU, i.e. if Param.simulationMode == 1)
 * 
 * \param [in] gpu_id Index of the GPU that is used for the computation (if Param.simulationMode == 1)
 */
template <typename RealType>
void CalculateExitWaveDerivativeSummand ( Image<RealType>& Summand,
                                          const Image<RealType>& ExitWave,
                                          const Image<RealType>& DiffImage,
                                          const RealType Focus,
                                          const Parameters<RealType>& Param,
                                          KernelData<RealType>& kData,
                                          const int gpu_id ) {
  // Check the size and format of Summand, ExitWave and DiffImage
  if ( !(Summand.getNumX ( ) == Param.X && Summand.getNumY ( ) == Param.Y &&
         Summand.correctFormat ( Space::FourierSpace, true ) ) )
    throw aol::Exception ( "Invalid size or format of Summand!", __FILE__, __LINE__ );
  
  if ( !(ExitWave.getNumX ( ) == Param.X && ExitWave.getNumY ( ) == Param.Y &&
         ExitWave.correctFormat ( Space::FourierSpace, true ) ) )
    throw aol::Exception ( "Invalid exit wave size or format!", __FILE__, __LINE__ );
  
  if ( !(DiffImage.getNumX ( ) == Param.X && DiffImage.getNumY ( ) == Param.Y &&
         DiffImage.correctFormat ( Space::FourierSpace, true ) ) )
    throw aol::Exception ( "Invalid size or format of DiffImage!", __FILE__, __LINE__ );
  
  switch ( Param.simulationMode ) {
    case 0:
      // focal integration approximation with periodic continuation (CPU)
      CalculateExitWaveDerivativeSummandFocalInt ( Summand, ExitWave, DiffImage, Focus, Param );
      break;
    case 1:
      // full computation of the weighted autocorrelation with zero continuation (GPU)
      CalculateExitWaveDerivativeSummandFullInt ( Summand, ExitWave, DiffImage, Focus, Param, kData, gpu_id );
      break;
    default:
      throw aol::Exception ( "Invalid simulation mode!", __FILE__, __LINE__ );
      break;
  }
}

#endif  // EWR_EXITWAVEDERIVATIVE_H
