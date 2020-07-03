#ifndef EWR_FUNCTIONAL_H
#define EWR_FUNCTIONAL_H

#include "ReconstructionBase.h"
#include "IntegrationDomains.h"

/**
 * \brief Calculates the inner part of the dataterm for the i-th summand
 * 
 * This function calculates I_sim(i) - I_modexp(i), where
 * 
 *   -I_sim(i) is the TEM image simulated from Arg for the i-th focus value,
 *   -I_modexp(i) is the i-th input image modulated by the i-th translation vector in Arg.
 * 
 * The result is calculated in Fourier space and stored in FourierSpaceDiffImage.
 * 
 * \param [in] kData OpenCL kernel data  (accessed only if the computation is to be performed on the GPU)
 * 
 * \param [in] gpu_id Index of the GPU that is used for the computation
 */
void Residual ( Image<RealType>& FourierSpaceDiffImage,
                const Arguments<RealType>& Arg,
                const InputData<RealType>& input,
                const int i,
                KernelData<RealType>& kData,
                const int gpu_id ) {
  // Check the size and format of FourierSpaceDiffImage
  if ( FourierSpaceDiffImage.getNumX ( ) != input.Param.X || FourierSpaceDiffImage.getNumY ( ) != input.Param.Y ||
       !FourierSpaceDiffImage.correctFormat ( Space::FourierSpace, true ) )
    throw aol::Exception ( "Invalid size or format of DiffImage!", __FILE__, __LINE__ );
  
  // Simulate a TEM image with the current estimate of the exit wave and focus value
  SimulateImage ( FourierSpaceDiffImage, Arg.ExitWave, Arg.Focus[i], input.Param, kData, gpu_id );
  
  // Modulate the i-th input image
  Image<RealType> ModulatedImage ( input.Param.X, input.Param.Y, Space::FourierSpace );
  input.ImageSeries[i].Modulate ( ModulatedImage, Arg.Translation.get ( i, 0 ), Arg.Translation.get ( i, 1 ), input.Param );
  
  // Calculate the difference of the simulated image and the modulated input image
  FourierSpaceDiffImage -= ModulatedImage;
}

/**
 * \brief Implements the numerical evaluation of the objective functional as an operator derived from aol::Op
 */
template <typename RealType>
class Functional : public aol::Op<Arguments<RealType>, aol::Scalar<RealType>> {
private:
  const InputData<RealType> input;  //!< Constant input data for the functional. (See ReconstructionBase.h)
  IntegrationDomains domains;       //!< Position and size of the integration domains. (See IntegrationDomains.h)
  Arguments<RealType> ScaleMask;    //!< Scale mask for the arguments that is applied prior to the evaluation of the functional
  KernelData<RealType>& kData;      //!< OpenCL kernel data for computations on the GPU(s)
  
  mutable Arguments<RealType> Arg;  //!< Current estimate of the exit wave, translation and focus values without the scaling from the scale mask
  
public:
  Functional ( const InputData<RealType> _input,
               const IntegrationDomains& _domains,
               const Arguments<RealType>& _ScaleMask,
               KernelData<RealType>& _kData )
    : input ( _input ),
      domains ( _domains ),
      ScaleMask ( _ScaleMask ),
      kData ( _kData ),
      Arg ( _ScaleMask, aol::STRUCT_COPY ) {
    input.checkValidity ( );
  }
  
  /**
   * \brief Evaluates the data term of the objective functional for the given arguments scaled by the scale mask
   * 
   * \param [in,out] RealSpaceDiffImages Optional parameter that can be used to store the inner part of the L2-norms of the data term in real space coordinates for later use.
   */
  RealType DataTerm ( const Arguments<RealType>& ScaledArg,
                      vector<Image<RealType>> *RealSpaceDiffImages = nullptr ) const {
    // RealSpaceDiffImages can be used to store the simulated images for later computations
    // unrelated to the energy computation here
    if ( RealSpaceDiffImages != nullptr )
      if ( static_cast<int> ( RealSpaceDiffImages->size ( ) ) != input.Param.N )
        throw aol::Exception ( "Invalid size of RealSpaceDiffImages!", __FILE__, __LINE__ );
    
    // Remove the scaling from the arguments
    Arg = ScaledArg;
    Arg *= ScaleMask;
    
    // Check the validity of Arg
    if ( !input.isCompatibleTo ( Arg ) )
      throw aol::Exception ( "Incompatible input data and argument!", __FILE__, __LINE__ );
    
    Arg.checkValidity ( input.Param.Focus[0] );
    
    // Calculate the value of the data term
    RealType result = 0;
    
    vector<RealType> partial_norms ( input.Param.N, 0 );
    
    const int nThreads = ( input.Param.simulationMode == 1 ? static_cast<int> ( kData.getNumGPUs ( ) ) : input.Param.numThreads );
    #pragma omp parallel for num_threads ( nThreads )
    for ( int i = 0; i < input.Param.N ; i++ ) {
      // Print a warning if the real-space domain of the modulated i-th input image does
      // not cover the entire integration domain
      const RealType tx_pixel = Arg.Translation.get ( i, 0 ) * input.Param.X / input.Param.lenX;
      const RealType ty_pixel = Arg.Translation.get ( i, 1 ) * input.Param.Y / input.Param.lenY;
      
      if ( -tx_pixel > domains[i].x ||
           -ty_pixel > domains[i].y ||
           -tx_pixel + input.Param.X - 2 * input.Param.reconstructionBufferZone < domains[i].x + domains[i].width ||
           -ty_pixel + input.Param.Y - 2 * input.Param.reconstructionBufferZone < domains[i].y + domains[i].height ) {
        #pragma omp critical ( DataTermResidualWarning )
        {
          cerr << endl << "Warning (i = " << i << "): the integration domain "
               << "\'(" << domains[i].x << ", " << domains[i].y << "), " << domains[i].width << " x " << domains[i].height << "\'"
               << " is not covered entirely by the image domain shifted by (" << tx_pixel << ", " << ty_pixel
               << ") pixel." << endl;
        }
      }
      
      // Calculate the squared 2-norm of the i-th residual on the given integration domain
      Image<RealType> FourierSpaceDiffImage ( input.Param.X, input.Param.Y, Space::FourierSpace );
      Residual ( FourierSpaceDiffImage, Arg, input, i, kData, omp_get_thread_num ( ) );
      
      Image<RealType> RealSpaceDiffImage ( input.Param.X, input.Param.Y, Space::RealSpace );
      FourierSpaceDiffImage.FourierTransformTo ( RealSpaceDiffImage ); // note that RealSpaceDiffImage is real-valued
      
      RealType norm = 0;
      for ( int y = 0; y < domains[i].height ; y++ ) for ( int x = 0; x < domains[i].width ; x++ )
        norm += aol::Sqr<RealType> ( RealSpaceDiffImage[0].get ( input.Param.reconstructionBufferZone + domains[i].x + x,
                                                                 input.Param.reconstructionBufferZone + domains[i].y + y ) );
      
      norm /= domains[i].width * domains[i].height;
      
      partial_norms[i] = norm;
      
      if ( RealSpaceDiffImages != nullptr )
        (*RealSpaceDiffImages)[i] = RealSpaceDiffImage;
    }
    
    for ( int i = 0; i < input.Param.N ; i++ )
      result += partial_norms[i];
    
    return result / input.Param.N;
  }
  
  /**
   * \brief Calculates the generalized non-linear Tikhonov regularizer for the given arguments scaled by the scale mask
   */
  RealType Regularizer ( const Arguments<RealType>& ScaledArg ) const {
    // Remove the scaling from the argument
    Arg = ScaledArg;
    Arg *= ScaleMask;
    
    // Check the validity of Arg
    if ( !input.isCompatibleTo ( Arg ) )
      throw aol::Exception ( "Incompatible input data and argument!", __FILE__, __LINE__ );
    
    Arg.checkValidity ( input.Param.Focus[0] );
    
    // Calculate the regularizer
    return input.Param.tikhonov_coeff * Arg.ExitWave.normSqr ( ) / ( static_cast<RealType> ( input.Param.X ) * input.Param.Y );
  }
  
  /**
   * \brief The objective functional's value is calculated as the sum of the data term and the regularizer
   */
  void apply ( const Arguments<RealType>& ScaledArg, aol::Scalar<RealType>& Dest ) const {
    Dest = DataTerm ( ScaledArg ) + Regularizer ( ScaledArg );
  }
  
  void applyAdd ( const Arguments<RealType>& ScaledArg, aol::Scalar<RealType>& Dest ) const {
    aol::Scalar<RealType> tmp;
    apply ( ScaledArg, tmp );
    Dest += tmp;
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
};

#endif  // EWR_FUNCTIONAL_H
