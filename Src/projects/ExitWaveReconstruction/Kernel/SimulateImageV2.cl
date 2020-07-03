/**
 * \brief OpenCL kernel for TEM image simulation
 * 
 * Each work item computes exactly one pixel (u, v) of the simulated
 * image, where u = get_global_id ( 0 ) and v = get_global_id ( 1 ).
 * 
 * This kernel is equivalent to SimulateImageV1, with the exception
 * that buffer objects instead of images are used to transmit data
 * between the host program and the GPU. Thus this kernel can also
 * be run on GPUs without OpenCL image support and / or with double
 * precision (if the GPU supports the double datatype).
 * 
 * 
 * Input: buffer object containing the real and imaginary parts of
 *   the exit wave and the phase transfer function.
 * 
 * SimulatedImage: buffer object containing the real and imaginary
 *   parts of the result.
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
 * 
 */
kernel void SimulateImageV2 ( global RealType4 *Input,
                              global RealType2 *SimulatedImage,
                              RealType focus ) {
  // Temporary variables
  RealType2 tmp1, tmp2;
  
  #define tmp                tmp2
  
  #define grad_chi_X         tmp1
  #define grad_chi_Y         tmp2
  
  #define exit_wave_product  tmp1
  #define transfer_function  tmp2
  
  // Current pixel
  const int u = get_global_id ( 0 );
  const int v = get_global_id ( 1 );
  
  // Resulting pixel value (complex valued)
  RealType2 integral = (RealType2)( 0, 0 );
  
  // Min and max values for x and y in the loops
  const int y_start = max ( DIM_Y_HALF - v, 0 );
  const int x_start = max ( DIM_X_HALF - u, 0 );
  const int y_end = min ( DIM_Y + DIM_Y_HALF - v, DIM_Y );
  const int x_end = min ( DIM_X + DIM_X_HALF - u, DIM_X );
  
  // Initial value for the y frequency
  RealType2 freqY = (RealType2)( INV_LEN_Y * ( y_start - DIM_Y_HALF ), INV_LEN_Y * ( y_start + v - DIM_Y ) );
  
  for ( int y = y_start; y < y_end ; ++y, freqY += (RealType2)( INV_LEN_Y, INV_LEN_Y ) ) {
    // Initial value for the x frequency
    RealType2 freqX = (RealType2)( INV_LEN_X * ( x_start - DIM_X_HALF ), INV_LEN_X * ( x_start + u - DIM_X ) );
    
    for ( int x = x_start; x < x_end ; ++x, freqX += (RealType2)( INV_LEN_X, INV_LEN_X ) ) {
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
      RealType4 p1 = Input[ ( x + u - DIM_X_HALF ) + ( y + v - DIM_Y_HALF ) * DIM_X ];
      RealType4 p2 = Input[ x + y * DIM_X ];
      
      RealType2 pure_phase_transfer = (RealType2)( p1.s2 * p2.s2 + p1.s3 * p2.s3, -p1.s2 * p2.s3 + p1.s3 * p2.s2 );
      
      // Compute product of the exit wave entries and the complete transfer function
      exit_wave_product = (RealType2)( p1.s0 * p2.s0 + p1.s1 * p2.s1, -p1.s0 * p2.s1 + p1.s1 * p2.s0 );
      transfer_function = spatial_coh * temporal_coh * pure_phase_transfer;
      
      // Compute final value of the current summand and add it to the result
      integral.s0 += exit_wave_product.s0 * transfer_function.s0 - exit_wave_product.s1 * transfer_function.s1;
      integral.s1 += exit_wave_product.s0 * transfer_function.s1 + exit_wave_product.s1 * transfer_function.s0;
    }
  }
  
  SimulatedImage[ u + v * DIM_X ] = integral;
}
