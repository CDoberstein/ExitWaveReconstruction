#ifndef EWR_PHYSICS_H
#define EWR_PHYSICS_H

#include "Parameters.h"

namespace Constants {
  const RealType c = 299792458;           // Speed of light in meters per second
  const RealType e = 1.60217662e-19;      // Proton charge in Coulomb
  const RealType me = 9.10938356e-31;     // Electron rest mass in kg
  const RealType h = 6.626070040e-34;     // Planck constant in Joule times seconds
  const RealType eps0 = 8.854187817e-21;  // Electric constant in Farads per nanometer
  const RealType a0 = 5.2917721067e-2;    // Bohr radius in nanometer
}

/**
 * Calculates the electron wavelength from the accelerating voltage U
 * in Volt and returns the result in nanometers.
 */
template <typename RealType>
RealType getElectronWavelength ( const RealType U ) {
  using namespace Constants;
  
  // The calculation is divided into three steps in order to avoid
  // errors with single precision floating point values
  const RealType v1 = h / sqrt ( 2 * me * U );
  const RealType v2 = 1 / sqrt ( e );
  const RealType v3 = 1 / sqrt ( 1 + e * U / ( 2 * me * c * c ) );
  
  return v1 * v2 * v3 * 1e9;
}

/**
 * Calculates the interaction parameter from the accelerating voltage U and returns the
 * result in radians / ( volt * nanometer ). The accelerating voltage unit is volts.
 * 
 * Formula taken from E. J. Kirkland, Advanced Computing in Electron Microscopy, page 65.
 */
template <typename RealType>
RealType getInteractionParameter ( const RealType U ) {
  using namespace Constants;
  
  // eU and mec2 are approximately of the same order of magnitude
  const RealType eU = e * U;
  const RealType mec2 = me * c * c;
  
  const RealType lambda = getElectronWavelength ( U );
  
  return 2 * aol::NumberTrait<RealType>::pi / ( lambda * U ) * ( ( mec2 + eU ) / ( 2 * mec2 + eU ) );
}

/**
 * Calculates the aberration function for a given frequency and aberration parameters.
 * The focus is expected to be given in nm and the frequencies are expected to be given
 * in nm^-1.
 */
template <typename RealType>
RealType AberrationFunction ( const RealType freqX, const RealType freqY, const RealType Focus, const Parameters<RealType>& Param ) {
  const RealType FreqAbsSqr = freqX * freqX + freqY * freqY;
  
  const RealType FocusTerm = 0.5 * Focus * Param.lambda * FreqAbsSqr;
  const RealType SphericalAberrationTerm = 0.25 * Param.SphericalAberration * aol::Cub<RealType> ( Param.lambda ) * aol::Sqr<RealType> ( FreqAbsSqr );
  
  return FocusTerm + SphericalAberrationTerm;
}

/**
 * Calculates the value of the aperture function for a given frequency and semiangle
 * alpha_max. The frequencies are assumed to be given in nm^-1.
 */
template <typename RealType>
int ApertureFunction ( const RealType freqX, const RealType freqY, const Parameters<RealType>& Param ) {
  const RealType freqAbs = sqrt ( freqX * freqX + freqY * freqY );
  
  return ( Param.lambda * freqAbs < Param.alpha_max ? 1 : 0 );
}

/**
 * Calculates the value of the spatial coherence envelope E_s(v,0) for a given
 * frequency v = (freqX, freqY).
 * 
 * The focus parameter is expected to be given in nm and the frequencies are expected to
 * be given in nm^-1.
 */
template <typename RealType>
RealType SpatialCoherenceEnvelope ( const RealType freqX, const RealType freqY, const RealType Focus, const Parameters<RealType>& Param ) {
  const RealType c1 = Param.lambda * freqX;
  const RealType c2 = Param.lambda * freqY;
  
  const RealType grad1 = Focus * c1 + Param.SphericalAberration * ( c1 * c1 * c1 + c1 * c2 * c2 );
  const RealType grad2 = Focus * c2 + Param.SphericalAberration * ( c2 * c2 * c2 + c1 * c1 * c2 );
  
  return exp ( -aol::Sqr<RealType> ( aol::NumberTrait<RealType>::pi * Param.alpha / Param.lambda ) * ( grad1 * grad1 + grad2 * grad2 ) );
}

/**
 * Interpretes the aperture function as a lowpass filter and calculates its radius in
 * pixels.
 * 
 * The aperture function can only be interpreted as a circular lowpass filter, if X and Y
 * are equal as well as lenX and lenY. (Otherwise its domain has the shape of an ellipse.)
 */
template <typename RealType>
RealType ApertureRadius ( const Parameters<RealType>& Param ) {
  if ( Param.X != Param.Y || Param.lenX != Param.lenY )
    return -1;
  
  return Param.alpha_max / Param.lambda * Param.lenX;
}
#endif  // EWR_PHYSICS_H
