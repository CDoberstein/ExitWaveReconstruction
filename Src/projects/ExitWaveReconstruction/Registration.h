#ifndef EWR_REGISTRATION_H
#define EWR_REGISTRATION_H

/// \cond
#include "Image.h"
/// \endcond

/**
 * \brief Register two images with the cross-correlation
 * 
 * This function determines a shift (x,y) with abs(x-xhint) <= maxAbsShift and
 * abs(y-yhint) <= maxAbsShift such that
 *   ReferenceImage(a+x,b+y) is approximately equal to RegImage(a,b)
 * for all a,b.
 */
template <typename RealType>
void CrossCorrelationRegistration ( const Image<RealType>& ReferenceImage,
                                    const Image<RealType>& RegImage,
                                    int maxAbsShift,
                                    RealType& x,
                                    RealType& y,
                                    int xhint = 0,
                                    int yhint = 0) {
  // Check the image sizes
  if ( ReferenceImage.getNumX ( ) != RegImage.getNumX ( ) ||
       ReferenceImage.getNumY ( ) != RegImage.getNumY ( ) )
    throw aol::Exception ( "Images have different sizes!", __FILE__, __LINE__ );
  
  const int X = ReferenceImage.getNumX ( );
  const int Y = ReferenceImage.getNumY ( );
  
  maxAbsShift = aol::Clamp<int> ( maxAbsShift, 0, min ( X / 2, Y / 2 ) );
  
  // Transform the images to Fourier space and remove the frequency shift if necessary
  Image<RealType> ReferenceImageF ( ReferenceImage );
  Image<RealType> RegImageF ( RegImage );
  
  ReferenceImageF.setFormat ( Space::FourierSpace, false );
  RegImageF.setFormat ( Space::FourierSpace, false );
  
  // Calculate the pointwise product of ReferenceImageF and the complex conjugate of
  // RegImageF and store the result in ReferenceImageF
  for ( int i = 0; i < ReferenceImageF[0].size ( ) ; i++ ) {
    const RealType a = ReferenceImageF[0][i];
    const RealType b = ReferenceImageF[1][i];
    
    const RealType c = RegImageF[0][i];
    const RealType d = -RegImageF[1][i];
    
    ReferenceImageF[0][i] = a * c - b * d;
    ReferenceImageF[1][i] = a * d + c * b;
  }
  
  // Subtract the mean
  ReferenceImageF[0].set ( 0, 0, 0 );
  ReferenceImageF[1].set ( 0, 0, 0 );
  
  // Calculate the inverse Fourier transform and store the result in RegImageF
  ReferenceImageF.FourierTransformTo ( RegImageF );
  
  // Determine the (x,y)-coordinates (within the window given by maxAbsShift) where
  // RegImageF has the largest value
  RealType max_value = 0;
  int xint = 0;
  int yint = 0;
  
  for ( int yy = -maxAbsShift + yhint; yy < maxAbsShift + yhint ; yy++ ) for ( int xx = -maxAbsShift + xhint; xx < maxAbsShift + xhint; xx++ ) {
    const int pos_x = ( xx + X ) % X;
    const int pos_y = ( yy + Y ) % Y;
    
    if ( RegImageF[0].get ( pos_x, pos_y ) > max_value ) {
      max_value = RegImageF[0].get ( pos_x, pos_y );
      xint = pos_x;
      yint = pos_y;
    }
  }
  
  if ( xint >= X / 2 ) xint -= X;
  if ( yint >= Y / 2 ) yint -= Y;
  
  // Consider the 3x3 window around (xint,yint) and take the average of the pixel
  // coordinates weighted by the pixel values to get a (hopefully) better result
  // than the integral shift (xint,yint)
  RealType sum = 0;
  x = 0;
  y = 0;
  for ( int dy = -1; dy <= 1 ; dy++ ) for ( int dx = -1; dx <= 1 ; dx++ ) {
    const int pos_x = ( xint + dx + X ) % X;
    const int pos_y = ( yint + dy + Y ) % Y;
    
    const RealType pxValue = RegImageF[0].get ( pos_x, pos_y );
    
    x += pxValue * ( xint + dx );
    y += pxValue * ( yint + dy );
    
    sum += pxValue;
  }
  
  x /= sum;
  y /= sum;
  
  x = aol::Clamp<RealType> ( x, -maxAbsShift, maxAbsShift );
  y = aol::Clamp<RealType> ( y, -maxAbsShift, maxAbsShift );
}

#endif  // EWR_REGISTRATION_H
