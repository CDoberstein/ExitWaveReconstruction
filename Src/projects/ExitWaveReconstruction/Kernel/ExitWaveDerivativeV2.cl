/**
 * \brief OpenCL kernel for the calculation of the exit wave derivative of the Fourier space functionals
 * 
 * This kernel only computes one summand of the exit wave derivative for
 * a given j. That is, in order to compute the complete derivative, this
 * kernel needs to be executed N times, the results have to be added and
 * finally multiplied by 4/N, where N is the number of images in the TEM
 * image series.
 * 
 * Each work item computes exactly one pixel (u, v) of the given part of
 * the derivative, where u = get_global_id ( 0 ) and v = get_global_id ( 1 ).
 * 
 * This kernel is equivalent to ExitWaveDerivativeV1, with the exception
 * that buffer objects instead of images are used to transmit data
 * between the host program and the GPU. Thus this kernel can also
 * be run on GPUs without OpenCL image support and / or with double
 * precision (if the GPU supports the double datatype).
 * 
 * 
 * Input: buffer object containing the real and imaginary parts of
 *   the exit wave and the phase transfer function.
 * 
 * A: buffer object containing the real and imaginary parts of the
 *   variable A in the formula for the derivative.
 * 
 * DerivativePart: buffer object containing the real and imaginary
 *   parts of the corresponding part of the derivative (the result).
 * 
 * focus: focus value that was used to compute the pure phase transfer
 *   function contained in Input.
 * 
 * 
 * The following variables are replaced with appropriate values when
 * building the kernel:
 *   - RealType, RealTypeN (float or double)
 *   - DIM_X, DIM_Y (image dimensions in pixel)
 *   - DIM_X_HALF, DIM_Y_HALF ( = DIM_X / 2, DIM_Y / 2)
 *   - INV_LEN_X, INV_LEN_Y (inverse of the image sizes in 1/nm)
 *   - WAVELENGTH (electron wavelength in nm)
 *   - WAVELENGTH_CUB (cubed electron wavelength in nm^3)
 *   - SPHERICAL_ABERRATION (spherical aberration in nm)
 *   - SPATIAL_COHERENCE_COEFF (equal to -(Pi*alpha/lambda)^2)
 *   - TEMPORAL_COHERENCE_COEFF (equal to -0.5*(Pi*Delta*lambda)^2)
 */
kernel void ExitWaveDerivativeV2 ( global RealType4 *Input,
                                   global RealType2 *A,
                                   global RealType2 *DerivativePart,
                                   RealType focus ) {
  // Temporary variables
  RealType2 tmp1, tmp2;
  
  #define tmp                  tmp2
  
  #define grad_chi_X           tmp1
  #define grad_chi_Y           tmp2
  
  #define pure_phase_transfer  tmp1
  #define partial_prod         tmp2
  
  // Current pixel
  const int u = get_global_id ( 0 );
  const int v = get_global_id ( 1 );
  
  // Resulting pixel value (complex valued)
  RealType2 integral = (RealType2)( 0, 0 );
  
  // Min and max values for x and y in the loops
  const int y_start = max ( v - DIM_Y_HALF, 0 );
  const int x_start = max ( u - DIM_X_HALF, 0 );
  const int y_end = min ( v + DIM_Y_HALF, DIM_Y );
  const int x_end = min ( u + DIM_X_HALF, DIM_X );
  
  // Initial value for the y frequency
  RealType2 freqY = (RealType2)( INV_LEN_Y * ( y_start - DIM_Y_HALF ), INV_LEN_Y * ( v - DIM_Y_HALF ) );
  
  for ( int y = y_start; y < y_end ; ++y, freqY.s0 += INV_LEN_Y ) {
    // Initial value for the x frequency
    RealType2 freqX = (RealType2)( INV_LEN_X * ( x_start - DIM_X_HALF ), INV_LEN_X * ( u - DIM_X_HALF ) );
    
    for ( int x = x_start; x < x_end ; ++x, freqX.s0 += INV_LEN_X ) {
      // Compute spatial coherence envelope
      tmp = focus * WAVELENGTH + SPHERICAL_ABERRATION * WAVELENGTH_CUB * ( freqX * freqX + freqY * freqY );
      
      grad_chi_X = freqX * tmp;
      grad_chi_Y = freqY * tmp;
      
      RealType spatial_coh = ( grad_chi_X.s0 - grad_chi_X.s1 ) * ( grad_chi_X.s0 - grad_chi_X.s1 ) + ( grad_chi_Y.s0 - grad_chi_Y.s1 ) * ( grad_chi_Y.s0 - grad_chi_Y.s1 );
      spatial_coh = exp ( SPATIAL_COHERENCE_COEFF * spatial_coh );
      
      // Compute temporal coherence envelope
      RealType temporal_coh = freqX.s0 * freqX.s0 + freqY.s0 * freqY.s0 - freqX.s1 * freqX.s1 - freqY.s1 * freqY.s1;
      temporal_coh = exp ( TEMPORAL_COHERENCE_COEFF * temporal_coh * temporal_coh );
      
      // Compute product of phase transfer functions (complex valued)
      RealType4 p1 = Input[ u + v * DIM_X ];
      RealType4 p2 = Input[ x + y * DIM_X ];
      
      pure_phase_transfer = (RealType2)( p1.s2 * p2.s2 + p1.s3 * p2.s3, -p1.s2 * p2.s3 + p1.s3 * p2.s2 );
      
      // Compute complete transfer function
      pure_phase_transfer *= spatial_coh * temporal_coh;
      
      // Compute product of the exit wave with the transfer function
      partial_prod = (RealType2)( p2.s0 * pure_phase_transfer.s0 + p2.s1 * pure_phase_transfer.s1,
                               p2.s0 * pure_phase_transfer.s1 - p2.s1 * pure_phase_transfer.s0 );
      
      // Compute final value of the current summand and add it to the result
      p1.s01 = A[ ( x - u + DIM_X_HALF ) + ( y - v + DIM_Y_HALF ) * DIM_X ];
      
      integral.s0 += p1.s0 * partial_prod.s0 - p1.s1 * partial_prod.s1;
      integral.s1 += p1.s0 * partial_prod.s1 + p1.s1 * partial_prod.s0;
    }
  }
  
  DerivativePart[ u + v * DIM_X ] = integral;
}

