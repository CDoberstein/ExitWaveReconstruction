#ifndef EWR_FUNCTIONALDERIVATIVE_H
#define EWR_FUNCTIONALDERIVATIVE_H

#include "ReconstructionBase.h"
#include "IntegrationDomains.h"
#include "ExitWaveDerivative.h"
#include "FocusDerivativeImage.h"

/**
 * \brief Implements the numerical calculation of the objective functional's derivative as an operator derived from aol::Op
 */
template <typename RealType>
class FunctionalDerivative : public aol::Op<Arguments<RealType>> {
private:
  const InputData<RealType> input;  //!< Constant input data for the functional. (See ReconstructionBase.h)
  IntegrationDomains domains;       //!< Position and size of the integration domains. (See IntegrationDomains.h)
  Arguments<RealType> ScaleMask;    //!< Scale mask for the arguments that is applied prior to the evaluation of the functional
  KernelData<RealType>& kData;      //!< OpenCL kernel data for computations on the GPU(s)
  
  int current_phase;                //!< Current phase id for the alternating minimization of the functional
  
  mutable Arguments<RealType> Arg;  //!< Current estimate of the exit wave, translation and focus values without the scaling from the scale mask
  
public:
  FunctionalDerivative ( const InputData<RealType> _input,
                         const IntegrationDomains& _domains,
                         const Arguments<RealType>& _ScaleMask,
                         KernelData<RealType>& _kData )
    : input ( _input ),
      domains ( _domains ),
      ScaleMask ( _ScaleMask ),
      kData ( _kData ),
      current_phase ( 0 ),
      Arg ( _ScaleMask, aol::STRUCT_COPY ) {
    input.checkValidity ( );
  }
  
  void apply ( const Arguments<RealType>& ScaledArg, Arguments<RealType>& Dest ) const {
    // Remove the scaling from the argument
    Arg = ScaledArg;
    Arg *= ScaleMask;
    
    // Check the validity of Arg
    if ( !input.isCompatibleTo ( Arg ) )
      throw aol::Exception ( "Incompatible input data and argument!", __FILE__, __LINE__ );
    
    Arg.checkValidity ( input.Param.Focus[0] );
    
    // Reset Dest
    Dest.setZero ( );
    
    vector<Image<RealType>> ExitWaveSummands ( input.Param.N, Image<RealType> ( input.Param.X, input.Param.Y, Space::FourierSpace ) );
    
    const int nThreads = ( input.Param.simulationMode == 1 ? static_cast<int> ( kData.getNumGPUs ( ) ) : input.Param.numThreads );
    #pragma omp parallel for num_threads ( nThreads )
    for ( int i = 0; i < input.Param.N ; i++ ) {
      // Simulate a TEM image with the current estimate of the exit wave and focus value
      Image<RealType> SimulatedImage ( input.Param.X, input.Param.Y, Space::FourierSpace );
      SimulateImage ( SimulatedImage, Arg.ExitWave, Arg.Focus[i], input.Param, kData, omp_get_thread_num ( ) );
      
      // Modulate the i-th input image
      Image<RealType> ModulatedImage ( input.Param.X, input.Param.Y, Space::FourierSpace );
      input.ImageSeries[i].Modulate ( ModulatedImage, Arg.Translation.get ( i, 0 ), Arg.Translation.get ( i, 1 ), input.Param );
      
      // Calculate the difference and set all pixels outside of the integration domain
      // to zero
      Image<RealType> DiffImage ( SimulatedImage );
      DiffImage -= ModulatedImage;
      
      Image<RealType> RealSpaceDiff ( input.Param.X, input.Param.Y, Space::RealSpace );
      DiffImage.FourierTransformTo ( RealSpaceDiff ); // note that RealSpaceDiff is real-valued
      
      Image<RealType> RealSpaceDiffSubsection ( input.Param.X, input.Param.Y, Space::RealSpace );
      for ( int y = 0; y < domains[i].height ; y++ ) for ( int x = 0; x < domains[i].width ; x++ ) {
        const int xCoord = x + input.Param.reconstructionBufferZone + domains[i].x;
        const int yCoord = y + input.Param.reconstructionBufferZone + domains[i].y;
        
        RealSpaceDiffSubsection[0].set ( xCoord, yCoord, RealSpaceDiff[0].get ( xCoord, yCoord ) );
      }
      
      RealSpaceDiffSubsection.FourierTransformTo ( DiffImage );
      
      DiffImage /= domains[i].width * domains[i].height;
      
      // Undo the normalization from the Fourier transforms of DiffImage
      DiffImage *= input.Param.X * input.Param.Y;
      
      // Derivative: exit wave
      // Calculates the i-th summand of the exit wave derivative
      // Note: the partial results are not added to Dest.ExitWave here, since the order
      //       of the additions would be essentially random if more than 1 thread is used.
      if ( input.Param.b_opt_exitwave[ current_phase ] ) {
        CalculateExitWaveDerivativeSummand ( ExitWaveSummands[i], Arg.ExitWave, DiffImage, Arg.Focus[i], input.Param, kData, omp_get_thread_num ( ) );
        ExitWaveSummands[i] *= static_cast<RealType> ( 4 ) / input.Param.N;
      }
      
      // Derivative: translation
      // Calculates the derivatives in x and y direction for the i-th image (only for i > 0
      // since the translation of the first image is fixed)
      if ( input.Param.b_opt_translation[ current_phase ] && i != 0 ) {
        RealType resX = 0;
        RealType resY = 0;
        
        const RealType invLenX = 1 / input.Param.lenX;
        const RealType invLenY = 1 / input.Param.lenY;
        
        for ( int y = 0; y < input.Param.Y ; y++ ) for ( int x = 0; x < input.Param.X ; x++ ) {
          const RealType freqX = invLenX * ( x - input.Param.X / 2 );
          const RealType freqY = invLenY * ( y - input.Param.Y / 2 );
          
          // Calculate the real part of DiffImage(x,y) * Conjugate(ModulatedImage(x,y)*i)
          const RealType realPart = -DiffImage[0].get ( x, y ) * ModulatedImage[1].get ( x, y ) +
                                     DiffImage[1].get ( x, y ) * ModulatedImage[0].get ( x, y );
          
          resX += 2 * aol::NumberTrait<RealType>::pi * freqX * realPart;
          resY += 2 * aol::NumberTrait<RealType>::pi * freqY * realPart;
        }
        
        Dest.Translation.set ( i, 0, -static_cast<RealType> ( 2 ) / input.Param.N * resX );
        Dest.Translation.set ( i, 1, -static_cast<RealType> ( 2 ) / input.Param.N * resY );
      }
      
      // Derivative: focus
      // Calculates the focus derivative for the i-th image (only for i > 0 since the
      // focus parameter of the first image is fixed)
      if ( input.Param.b_opt_focus[ current_phase ] && i != 0 ) {
        Image<RealType> FocusDerivativeImage ( input.Param.X, input.Param.Y, Space::FourierSpace );
        CalculateFocusDerivativeImage ( FocusDerivativeImage, Arg.ExitWave, Arg.Focus[i], input.Param, kData, omp_get_thread_num ( ) );
        
        Dest.Focus[i] = DiffImage * FocusDerivativeImage;
        Dest.Focus[i] *= static_cast<RealType> ( 2 ) / input.Param.N;
      }
    }
    
    // Addition of floating point numbers is not associative. In order to ensure the
    // reproducibility of the results, the exit wave summands are added to Dest.ExitWave
    // in the fixed order 0, 1, ..., N-1.
    for ( int i = 0; i < input.Param.N ; i++ )
      Dest.ExitWave += ExitWaveSummands[i];
    
    // Scale the result with the scale mask (chain rule)
    Dest *= ScaleMask;
  }
  
  void applyAdd ( const Arguments<RealType>& Arg, Arguments<RealType>& Dest ) const {
    throw aol::UnimplementedCodeException ( "Not implemented!", __FILE__, __LINE__ );
    
    aol::doNothingWithArgumentToPreventUnusedParameterWarning ( Arg );
    aol::doNothingWithArgumentToPreventUnusedParameterWarning ( Dest );
  }
  
  /**
   * \brief Adjusts the integration domains
   */
  void setIntegrationDomains ( const IntegrationDomains& new_domains ) {
    domains = new_domains;
  }
  
  /**
   * \brief Set a new scale mask
   */
  void setScaleMask ( const Arguments<RealType>& new_ScaleMask ) {
    ScaleMask = new_ScaleMask;
  }
  
  /**
   * \brief Set the current phase id for the alternating minimization of the functional
   */
  void setOptimizationPhase ( const int phase ) {
    current_phase = phase;
  }
};

#endif  // EWR_FUNCTIONALDERIVATIVE_H
