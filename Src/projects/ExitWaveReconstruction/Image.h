#ifndef EWR_IMAGE_H
#define EWR_IMAGE_H

/// \cond
#include "Parameters.h"

#include <convolution.h>
#include <multiArray.h>
/// \endcond

enum class Space { RealSpace, FourierSpace };

/**
 * \brief Utility class for complex-valued images that keeps track of the image space (i.e. real space or Fourier space) and the frequency shift
 */
template <typename RealType>
class Image : public qc::MultiArray<RealType, 2, 2> {
private:
  Space image_space;  //!< Either Space::RealSpace or Space::FourierSpace
  bool freq_shift;    //!< true if a frequency shift was applied to the image, false otherwise
  
public:
  /**
   * \brief Default constructor for the construction of an empty image in real space
   */
  Image ( const int sizeX = 0,
          const int sizeY = 0,
          const Space _image_space = Space::RealSpace )
    : qc::MultiArray<RealType, 2, 2> ( sizeX, sizeY ),
      image_space ( _image_space ),
      freq_shift ( _image_space == Space::RealSpace ? false : true ) { }
  
  /**
   * \brief Constructs an image of the given size and format
   */
  Image ( const int sizeX,
          const int sizeY,
          const Space _image_space,
          const bool _freq_shift )
    : qc::MultiArray<RealType, 2, 2> ( sizeX, sizeY ),
      image_space ( _image_space ),
      freq_shift ( _freq_shift ) { }
  
  /**
   * \brief Constructs an image from a qc::MultiArray with the given format
   */
  Image ( const qc::MultiArray<RealType, 2, 2>& img,
          const Space _image_space,
          const bool _freq_shift,
          aol::CopyFlag copytype = aol::DEEP_COPY )
    : qc::MultiArray<RealType, 2, 2> ( img, copytype ),
      image_space ( _image_space ),
      freq_shift ( _freq_shift ) { }
  
  /**
   * \returns the current image space
   */
  Space getImageSpace ( ) const { return image_space; }
  
  /**
   * \returns the current frequency shift boolean (true if a frequency shift was applied to the image, false otherwise)
   */
  bool getFreqShift ( ) const { return freq_shift; }
  
  /**
   * \brief Checks if the image format coincides with the format given by the arguments to this function
   * 
   * \returns true if the formats coincide, false otherwise
   */
  bool correctFormat ( const Space correct_image_space, const bool correct_freq_shift ) const {
    return ( correct_image_space == image_space && correct_freq_shift == freq_shift );
  }
  
  /**
   * \brief Transforms the image to the given space and performs a frequency shift if necessary
   */
  void setFormat ( const Space target_image_space, const bool target_freq_shift ) {
    if ( image_space != target_image_space ) {
      Image<RealType> tmp ( *this, image_space, freq_shift );
      tmp.FourierTransformTo ( *this );
    }
    
    if ( freq_shift != target_freq_shift )
      FrequencyShift ( );
  }
  
  /**
   * \brief Performs a 2D Fourier transformation of *this using FFTW
   * 
   * If *this is in real space coordinates, a forward Fourier transformation is
   * performed, otherwise a inverse Fourier transformation.
   * 
   * \note In addition to the Fourier transformation this function automatically
   *       applies frequency shifts to ensure that a resulting real space image has
   *       no frequency shift whereas a resulting Fourier space image is returned with
   *       a frequency shift.
   * 
   * \param [in,out] transformed_img the resulting transformed image; transformed_img is
   *        expected to be of the same size as *this.
   */
  void FourierTransformTo ( Image<RealType>& transformed_img ) const {
    if ( this->getNumX ( ) % 2 != 0 || this->getNumY ( ) % 2 != 0 )
      throw aol::Exception ( "Image::FourierTransformTo() is not implemented for odd image sizes!", __FILE__, __LINE__ );
    
    if ( transformed_img.getNumX ( ) != this->getNumX ( ) || transformed_img.getNumY ( ) != this->getNumY ( ) )
      throw aol::Exception ( "Invalid image size!", __FILE__, __LINE__ );
    
    if ( image_space == Space::RealSpace ) {
      // Forward transformation
      qc::FourierTransform ( *this, transformed_img, qc::FTForward );
      
      // Normalize the result
      transformed_img /= static_cast<RealType> ( this->getNumX ( ) ) * this->getNumY ( );
      
      // Set the member variables appropriately
      transformed_img.image_space = Space::FourierSpace;
      transformed_img.freq_shift = freq_shift;
      
      // Perform a frequency shift if necessary
      if ( !transformed_img.freq_shift )
        transformed_img.FrequencyShift ( );
      
    } else {
      // Inverse transformation
      qc::FourierTransform ( *this, transformed_img, qc::FTBackward );
      
      // Set member variables appropriately
      transformed_img.image_space = Space::RealSpace;
      transformed_img.freq_shift = freq_shift;
      
      // Perform a frequency shift if necessary
      if ( transformed_img.freq_shift )
        transformed_img.FrequencyShift ( );
    }
  }
  
  /**
   * \brief Performs a 2D frequency shift
   * 
   * In Fourier space this function divides the image in four evenly sized rectangles
   * and interchanges the top left with the bottom right rectangle and the top right
   * with the bottom left rectangle. (This is equivalent to applying a translation by
   * half the image size in each direction.)
   * 
   * In real space the image is modulated by a modulation that corresponds exactly to
   * the translation (frequency shift) in Fourier space.
   */
  void FrequencyShift ( ) {
    if ( this->getNumX ( ) % 2 != 0 || this->getNumY ( ) % 2 != 0 )
      throw aol::Exception ( "Image::FrequencyShift() is not implemented for odd image sizes!", __FILE__, __LINE__ );
    
    if ( image_space == Space::RealSpace ) {
      // Modulate image
      for ( int i = 0; i < 2 ; i++ ) {
        int x_start = 1;
        for ( int y = 0; y < this->getNumY ( ) ; y++, x_start = 1 - x_start )
          for ( int x = x_start; x < this->getNumX ( ) ; x += 2 )
            (*this)[i].getReference ( x, y ) *= -1;
      }
    } else {
      // Translate image (frequency shift)
      const int Xh = this->getNumX ( ) / 2;
      const int Yh = this->getNumY ( ) / 2;
      
      for ( int i = 0; i < 2 ; i++ )
        for ( int y = 0; y < Yh ; y++ ) for ( int x = 0; x < Xh ; x++ ) {
          std::swap ( (*this)[i].getReference ( x, y ), (*this)[i].getReference ( x + Xh, y + Yh ) );
          std::swap ( (*this)[i].getReference ( x + Xh, y ), (*this)[i].getReference ( x, y + Yh ) );
        }
    }
    
    freq_shift = !freq_shift;
  }
  
  /**
   * \brief Applies a low-pass filter with the given cutoff frequency (in nm^-1) to *this
   * 
   * If cutoff_frequency == -1, the cutoff frequency is set maximal such that the
   * elliptical domain of the lowpass filter is still completely contained within the
   * image domain.
   * 
   * \param [in] Param Microscope parameters (required for the conversion from pixels to nm)
   * 
   * \returns the number of pixels that were not set to zero.
   */
  int LowpassFilter ( const Parameters<RealType>& Param,
                      const RealType cutoff_frequency = -1 ) {
    const int X = this->getNumX ( );
    const int Y = this->getNumY ( );
    
    const RealType lenX = ( Param.lenX * X ) / Param.X;
    const RealType lenY = ( Param.lenY * Y ) / Param.Y;
    
    // Calculate the squared cutoff_frequency
    RealType cutoff_frequency_sqr = cutoff_frequency * cutoff_frequency;
    if ( cutoff_frequency == -1 ) {
      cutoff_frequency_sqr = min ( X / ( 2 * lenX ), Y / ( 2 * lenY ) );
      cutoff_frequency_sqr *= cutoff_frequency_sqr;
    }
    
    const bool orig_freq_shift = freq_shift;
    
    const RealType invLenX = 1 / lenX;
    const RealType invLenY = 1 / lenY;
    
    int count = X * Y;
    
    if ( image_space == Space::RealSpace ) {
      // Convert the image to Fourier space
      Image<RealType> fimg ( X, Y );
      FourierTransformTo ( fimg );
      
      // Apply the low-pass filter
      for ( int y = 0; y < Y ; y++ ) for ( int x = 0; x < X ; x++ ) {
        const RealType freqX = invLenX * ( x - X / 2 );
        const RealType freqY = invLenY * ( y - Y / 2 );
        if ( aol::Sqr<RealType> ( freqX ) + aol::Sqr<RealType> ( freqY ) > cutoff_frequency_sqr ) {
          fimg[0].set ( x, y, 0 );
          fimg[1].set ( x, y, 0 );
          --count;
        }
      }
      
      // Convert result back to real space
      fimg.FourierTransformTo ( *this );
    } else {
      // Make sure that the lowest frequency is located at the image's center
      if ( !orig_freq_shift )
        FrequencyShift ( );
      
      // Apply the low-pass filter
      for ( int y = 0; y < Y ; y++ ) for ( int x = 0; x < X ; x++ ) {
        const RealType freqX = invLenX * ( x - X / 2 );
        const RealType freqY = invLenY * ( y - Y / 2 );
        if ( aol::Sqr<RealType> ( freqX ) + aol::Sqr<RealType> ( freqY ) > cutoff_frequency_sqr ) {
          (*this)[0].set ( x, y, 0 );
          (*this)[1].set ( x, y, 0 );
          --count;
        }
      }
    }
    
    // Reset the frequency shift to the original state
    if ( freq_shift != orig_freq_shift )
      FrequencyShift ( );
    
    return count;
  }
  
  /**
   * \brief Rotate the image by angle degrees using bilinear interpolation and constant continuation
   * 
   * If adjust_size is true, the rotated image's size is reduced to the largest axis
   * aligned rectangle with the same aspect ratio as the original image such that the
   * result is contained entirely within the rotated domain of the original image.
   * 
   * \param [in] fillPx value that is used for pixels beyond the image borders when performing the bilinear interpolation
   */
  void Rotate ( RealType angle,
                bool adjust_size = false,
                RealType fillPx = 0 ) {
    // No rotation is necessary if the angle is equal to zero
    if ( angle == 0 )
      return;
    
    // Only allow rotation of images with a format where the rotation is mathematically meaningful
    if ( !( correctFormat ( Space::RealSpace, false ) || correctFormat ( Space::FourierSpace, true ) ) )
      throw aol::Exception ( "Unexpected image format!", __FILE__, __LINE__ );
    
    // Perform the rotation using bilinear interpolation
    qc::MultiArray<RealType, 2, 2> copy ( *this );
    
    const RealType cos_angle = cos ( angle / 180 * aol::NumberTrait<RealType>::pi );
    const RealType sin_angle = sin ( angle / 180 * aol::NumberTrait<RealType>::pi );
    
    const int X = this->getNumX ( );
    const int Y = this->getNumY ( );
    
    const RealType x_half = static_cast<RealType> ( X ) / 2;
    const RealType y_half = static_cast<RealType> ( Y ) / 2;
    for ( int i = 0; i < 2 ; i++ )
      for ( int y = 0; y < Y ; y++ ) for ( int x = 0; x < X ; x++ ) {
        RealType v[2] = { x - x_half, y - y_half };
        RealType v_rot[2] = { cos_angle * v[0] + sin_angle * v[1],
                             -sin_angle * v[0] + cos_angle * v[1] };
        
        v_rot[0] += x_half;
        v_rot[1] += y_half;
        
        int xint = static_cast<int> ( v_rot[0] );
        int yint = static_cast<int> ( v_rot[1] );
        
        RealType xdiff = v_rot[0] - xint;
        RealType ydiff = v_rot[1] - yint;
        
        RealType val = 0;
        if ( xint < 0 || yint < 0 || xint >= X || yint >= Y )
          val += ( 1 - xdiff ) * ( 1 - ydiff ) * fillPx;
        else 
          val += ( 1 - xdiff ) * ( 1 - ydiff ) * copy[i].get ( xint, yint );
          
        if ( xint+1 < 0 || yint < 0 || xint+1 >= X || yint >= Y )
          val += xdiff * ( 1 - ydiff ) * fillPx;
        else
          val += xdiff * ( 1 - ydiff ) * copy[i].get ( xint+1, yint );
          
        if ( xint < 0 || yint+1 < 0 || xint >= X || yint+1 >= Y ) 
          val += ( 1 - xdiff ) * ydiff * fillPx;
        else
          val += ( 1 - xdiff ) * ydiff * copy[i].get ( xint, yint+1 );
          
        if ( xint+1 < 0 || yint+1 < 0 || xint+1 >= X || yint+1 >= Y )
          val += xdiff * ydiff * fillPx;
        else
          val += xdiff * ydiff * copy[i].get ( xint+1, yint+1 );
        
        (*this)[i].set ( x, y, val );
      }
    
    if ( adjust_size ) {
      // Calculate the coordinates of one corner of the new rectangular domain
      int angle_mod = static_cast<int> ( angle ) % 180;
      if ( angle_mod < 0 )
        angle_mod += 180;
      const RealType sign = ( angle_mod >= 90 ? -1 : 1 );
      
      const int x = static_cast<int> ( abs ( X * X / ( 2 * cos_angle * X + 2 * sign * sin_angle * Y ) ) );
      const int y = static_cast<int> ( abs ( X * Y / ( 2 * cos_angle * X + 2 * sign * sin_angle * Y ) ) );
      
      // Crop the rotated image
      Image<RealType> copy ( *this );
      
      this->reallocate ( 2 * x, 2 * y );
      
      copy[0].copyBlockTo ( X / 2 - x, Y / 2 - y, (*this)[0] );
      copy[1].copyBlockTo ( X / 2 - x, Y / 2 - y, (*this)[1] );
    }
  }
  
  /**
   * \brief Modulate the image by (dx, dy) nanometers
   * 
   * \param [out] result the modulated image; result is expected to be of the same size as *this
   * 
   * \param [in] Param Microscope parameters (required for the conversion from pixels to nm)
   */
  void Modulate ( Image<RealType>& result,
                  const RealType dx,
                  const RealType dy,
                  const Parameters<RealType>& Param ) const {
    // Only allow the modulation for Fourier space images with the highest frequency in the center of the image
    if ( !(correctFormat ( Space::FourierSpace, true ) && result.correctFormat ( Space::FourierSpace, true ) ) )
      throw aol::Exception ( "Invalid image format!", __FILE__, __LINE__ );
    
    if ( this->getNumX ( ) != result.getNumX ( ) || this->getNumY ( ) != result.getNumY ( ) )
      throw aol::Exception ( "Invalid image size!", __FILE__, __LINE__ );
    
    const int X = this->getNumX ( );
    const int Y = this->getNumY ( );
    
    const RealType invLenX = 1 / ( Param.lenX * X / static_cast<RealType> ( Param.X ) );
    const RealType invLenY = 1 / ( Param.lenY * Y / static_cast<RealType> ( Param.Y ) );
    
    const RealType Xhalf = static_cast<RealType> ( X ) / 2;
    const RealType Yhalf = static_cast<RealType> ( Y ) / 2;
    
    for ( int y = 0; y < Y ; ++y ) for ( int x = 0; x < X ; ++x ) {
      const RealType freqX = invLenX * ( x - Xhalf );
      const RealType freqY = invLenY * ( y - Yhalf );
      
      const RealType arg = 2 * aol::NumberTrait<RealType>::pi * ( freqX * dx + freqY * dy );
      
      const RealType cos_arg = cos ( arg );
      const RealType sin_arg = sin ( arg );
      
      const RealType img_Re = (*this)[0].get ( x, y );
      const RealType img_Im = (*this)[1].get ( x, y );
      
      result[0].set ( x, y, img_Re * cos_arg - img_Im * sin_arg );
      result[1].set ( x, y, img_Re * sin_arg + img_Im * cos_arg );
    }
  }
  
  /**
   * \brief Shift the image by (dx, dy) nanometers using bilinear interpolation with periodic continuation
   * 
   * \param [out] result the shifted image; result is expected to be of the same size as *this
   * 
   * \param [in] Param Microscope parameters (required for the conversion from pixels to nm)
   */
  void Shift ( Image<RealType>& result,
               const RealType dx,
               const RealType dy,
               const Parameters<RealType>& Param ) const {
    // Only allow an image shift for real space images without a frequency shift
    if ( !(correctFormat ( Space::RealSpace, false ) && result.correctFormat ( Space::RealSpace, false ) ) )
      throw aol::Exception ( "Invalid image format!", __FILE__, __LINE__ );
    
    if ( this->getNumX ( ) != result.getNumX ( ) || this->getNumY ( ) != result.getNumY ( ) )
      throw aol::Exception ( "Invalid image size!", __FILE__, __LINE__ );
    
    const int X = this->getNumX ( );
    const int Y = this->getNumY ( );
    
    const RealType dx_pixel = dx * Param.X / Param.lenX;
    const RealType dy_pixel = dy * Param.Y / Param.lenY;
    
    const int dx_int = static_cast<int> ( floor ( dx_pixel ) ) % X;
    const int dy_int = static_cast<int> ( floor ( dy_pixel ) ) % Y;
    
    const RealType dx_diff = dx_pixel - dx_int;
    const RealType dy_diff = dy_pixel - dy_int;
    
    for ( int i = 0; i < 2 ; i++ ) {
      for ( int y = 0; y < Y ; ++y ) for ( int x = 0; x < X ; ++x ) {
        int ul[2] = { ( x + dx_int + X ) % X, ( y + dy_int + Y ) % Y };
        int lr[2] = { ( x + dx_int + 1 + X ) % X, ( y + dy_int + 1 + Y ) % Y };
        
        result[i].set ( x, y, (*this)[i].get ( ul[0], ul[1] ) * ( 1 - dx_diff ) * ( 1 - dy_diff ) +
                              (*this)[i].get ( ul[0], lr[1] ) * ( 1 - dx_diff ) * dy_diff +
                              (*this)[i].get ( lr[0], ul[1] ) * dx_diff * ( 1 - dy_diff ) +
                              (*this)[i].get ( lr[0], lr[1] ) * dx_diff * dy_diff );
      }
    }
  }
  
  /**
   * \brief Calculate the amplitude and phase of the image
   * 
   * The amplitude is written to the first component of result and the phase
   * is written to the second component.
   */
  void CalculateAmplitudeAndPhase ( qc::MultiArray<RealType, 2, 2>& result ) const {
    if ( this->getNumX ( ) != result.getNumX ( ) || this->getNumY ( ) != result.getNumY ( ) )
      throw aol::Exception ( "Invalid image size!", __FILE__, __LINE__ );
    
    for ( int i = 0; i < (*this)[0].size ( ) ; i++ ) {
      const RealType real = (*this)[0][i];
      const RealType imag = (*this)[1][i];
      
      // Amplitude
      result[0][i] = sqrt ( real * real + imag * imag );
      
      // Phase
      if ( real == 0 && imag == 0 )
        result[1][i] = 0;
      else
        result[1][i] = atan2 ( imag, real );
    }
  }
  
  /**
   * \brief Checks if all pixels of the image are equal to zero
   * 
   * \returns true if all pixels are equal to zero, false otherwise
   */
  bool isZero ( ) const {
    for ( int i = 0; i < 2; i++ )
      for ( int y = 0; y < this->getNumY ( ) ; y++ ) for ( int x = 0; x < this->getNumX ( ) ; x++ )
        if ( (*this)[i].get ( x, y ) != static_cast<RealType> ( 0 ) )
          return false;
    return true;
  }
  
  /**
   * \brief Removes a buffer zone of the size buffer_zone_size in real space from *this and stores the resulting smaller image in result
   */
  void removeRealSpaceBufferZone ( Image<RealType>& result,
                                   const int buffer_zone_size ) const {
    if ( result.getNumX ( ) != this->getNumX ( ) - 2 * buffer_zone_size ||
         result.getNumY ( ) != this->getNumY ( ) - 2 * buffer_zone_size )
      throw aol::Exception ( "Invalid image size!", __FILE__, __LINE__ );
    
    Image<RealType> RealSpaceImage ( *this );
    RealSpaceImage.setFormat ( Space::RealSpace, false );
    
    const bool target_freq_shift = result.freq_shift;
    if ( result.image_space == Space::RealSpace ) {
      result.freq_shift = false;
      
      RealSpaceImage[0].copyBlockTo ( buffer_zone_size, buffer_zone_size, result[0] );
      RealSpaceImage[1].copyBlockTo ( buffer_zone_size, buffer_zone_size, result[1] );
    } else {
      Image<RealType> RealSpaceImageNoBuf ( this->getNumX ( ) - 2 * buffer_zone_size, this->getNumY ( ) - 2 * buffer_zone_size, Space::RealSpace );
      
      RealSpaceImage[0].copyBlockTo ( buffer_zone_size, buffer_zone_size, RealSpaceImageNoBuf[0] );
      RealSpaceImage[1].copyBlockTo ( buffer_zone_size, buffer_zone_size, RealSpaceImageNoBuf[1] );
      
      RealSpaceImageNoBuf.FourierTransformTo ( result );
    }
    
    if ( result.freq_shift != target_freq_shift )
      result.FrequencyShift ( );
  }
};

/**
 * \brief Save an image in the Quocmesh .q2bz format or in the .tiff format
 * 
 * \param [in] image the image to be saved
 * 
 * \param [in] real_valued if true, only the first component of image is saved. Otherwise both the real and the imaginary component are saved in separate files
 * 
 * \param [in] path the path and base file name of the file without extension
 * 
 * \param [in] APimage if true, then "Amplitude" and "Phase" are used as file name suffices instead of "Re" and "Im"
 * 
 * \param [in] output_image_format integer describing the output image format according to outputImageFormat from the parameter file
 */
void SaveImage ( const qc::MultiArray<RealType, qc::QC_2D>& image,
                 const bool real_valued,
                 const string& path,
                 const int output_image_format,
                 const bool APimage = false ) {
  switch ( output_image_format ) {
    case 0: // .q2bz
      if ( real_valued ) {
        image[0].saveToFile ( ( path + ".q2bz" ).c_str ( ) );
      } else {
        image[0].saveToFile ( ( path + ( APimage ? "_Amplitude.q2bz" : "_Re.q2bz" ) ).c_str ( ) );
        image[1].saveToFile ( ( path + ( APimage ? "_Phase.q2bz" : "_Im.q2bz" ) ).c_str ( ) );
      }
      break;
    case 1: // .tiff
      if ( real_valued ) {
        image[0].saveTIFF ( ( path + ".tiff" ).c_str ( ) );
      } else {
        image[0].saveTIFF ( ( path + ( APimage ? "_Amplitude.tiff" : "_Re.tiff" ) ).c_str ( ) );
        image[1].saveTIFF ( ( path + ( APimage ? "_Phase.tiff" : "_Im.tiff" ) ).c_str ( ) );
      }
      break;
    default:
      throw aol::Exception ( "Invalid output image format!", __FILE__, __LINE__ );
      break;
  }
}

#endif  // EWR_IMAGE_H
