#ifndef EWR_PRECONDITIONER_H
#define EWR_PRECONDITIONER_H

#include "Functional.h"
#include "FunctionalDerivative.h"

/**
 * \brief Auxiliary function to find the optimal power-of-2 stepsize in the given direction starting with hint as the initial exponent
 */
int FindStepsize ( const Arguments<RealType>& ScaledPos,
                   const Arguments<RealType>& Dir,
                   const Functional<RealType>& energy,
                   const RealType pos_energy,
                   const int hint = 0 ) {
  const int max_abs_exp = 20;
  
  int exp = hint;
  RealType tau = pow ( 2.0, hint );
  
  Arguments<RealType> tmp ( Dir );
  
  aol::Scalar<RealType> cur_energy;
  tmp = Dir; tmp *= tau; tmp += ScaledPos;
  energy.apply ( tmp, cur_energy );
  
  if ( cur_energy < pos_energy ) {
    ++exp;
    tau *= 2;
    aol::Scalar<RealType> min_energy;
    tmp = Dir; tmp *= tau; tmp += ScaledPos;
    energy.apply ( tmp, min_energy );
    if ( min_energy < cur_energy ) {
      do {
        ++exp;
        tau *= 2;
        tmp = Dir; tmp *= tau; tmp += ScaledPos;
        energy.apply ( tmp, cur_energy );
        if ( cur_energy >= min_energy ) return exp - 1;
        min_energy = cur_energy;
      } while ( abs(exp) < max_abs_exp );
      return exp;
    }
  } else {
    do {
      --exp;
      tau /= 2;
      tmp = Dir; tmp *= tau; tmp += ScaledPos;
      energy.apply ( tmp, cur_energy );
    } while ( cur_energy >= pos_energy && abs(exp) < max_abs_exp );
  
    if ( abs(exp) >= max_abs_exp ) return exp;
  }
  
  aol::Scalar<RealType> min_energy = cur_energy;
  do {
    --exp;
    tau /= 2;
    tmp = Dir; tmp *= tau; tmp += ScaledPos;
    energy.apply ( tmp, cur_energy );
    if ( cur_energy >= min_energy ) return exp + 1;
    min_energy = cur_energy;
  } while ( abs(exp) < max_abs_exp );
  
  return exp;
}

/**
 * \brief Calculate a circular exit wave scale mask
 * 
 * \param [out] ExitWaveScaleMask the resulting scale mask
 * 
 * \param [in] ScaledPos current estimate for the arguments (i.e. exit wave, translations and focus values)
 * 
 * \param [in] energy energy operator for the evaluation of the objective functional
 * 
 * \param [in] derivative derivative operator for the evaluation of the objective functional's derivative
 * 
 * \param [in] Param microscope parameters
 */
void CreateCircularExitWaveScaleMask ( Image<RealType>& ExitWaveScaleMask,
                                       const Arguments<RealType>& ScaledPos,
                                       const Functional<RealType>& energy,
                                       const FunctionalDerivative<RealType>& derivative,
                                       const Parameters<RealType>& Param ) {
  const int X = ScaledPos.ExitWave.getNumX ( );
  
  if ( X != ScaledPos.ExitWave.getNumY ( ) || Param.lenX != Param.lenY )
    throw aol::Exception ( "The computation of a circular exit wave scale mask is only implemented for square images!", __FILE__, __LINE__ );
  
  if ( ExitWaveScaleMask.getNumX ( ) != X || ExitWaveScaleMask.getNumY ( ) != X )
    throw aol::Exception ( "Invalid format of ExitWaveScaleMask!", __FILE__, __LINE__ );
  
  cerr << "Updating the exit wave scale mask (circular shape)... ";
  
  // Calculate the energy at the current position
  aol::Scalar<RealType> pos_energy;
  energy.apply ( ScaledPos, pos_energy );
  
  // Calculate the derivative
  Arguments<RealType> negDir ( ScaledPos, aol::STRUCT_COPY );
  derivative.apply ( ScaledPos, negDir );
  
  // Calculate the optimal power-of-2 stepsizes for the pixels on the diagonal (d,d)
  Arguments<RealType> pxDir ( negDir, aol::STRUCT_COPY );
  vector<int> stepsize_exp[2] = { vector<int> ( X / 2 + 1 ), vector<int> ( X / 2 + 1 ) };
  
  int stepsize_hint = 0;
  for ( int j = 0; j < 2 ; j++ )
    for ( int d = 0; d <= X / 2 ; d++ ) {
      cerr << "\rUpdating the exit wave scale mask (circular shape)... " << setw ( 5 ) << setprecision ( 2 )
           << 100 * static_cast<RealType> ( d + j * ( X / 2 + 1 ) ) / ( 2 * ( X / 2 + 1 ) ) << "%";
      
      pxDir.ExitWave[j].set ( d, d, -negDir.ExitWave[j].get ( d, d ) );
      /*stepsize_hint =*/ stepsize_exp[j][ X / 2 - d ] = FindStepsize ( ScaledPos, pxDir, energy, pos_energy, stepsize_hint );
      pxDir.ExitWave[j].set ( d, d, 0 );
    }
  
  // Smooth the stepsize vectors with the mean value kernel of size 2 * kernel_size + 1
  // and constant continuation
  vector<int> smooth_stepsize_exp[2] = { vector<int> ( X / 2 + 1 ), vector<int> ( X / 2 + 1 ) };
  const int kernel_size = 2;
  
  for ( int j = 0; j < 2 ; j++ )
    for ( int d = 0; d <= X / 2 ; d++ ) {
      int sum = 0;
      for ( int k = -kernel_size ; k <= kernel_size ; k++ )
        sum += stepsize_exp[j][ aol::Clamp ( d + k, 0, X / 2 ) ];
      
      smooth_stepsize_exp[j][d] = sum / ( 2 * kernel_size + 1 );
    }
  
  // Calculate the scale mask
  for ( int j = 0; j < 2 ; j++ )
    for ( int y = 0; y < X ; y++ )
      for ( int x = 0; x < X ; x++ ) {
        // Calculate the distance of (x,y) to (X/2,X/2)
        RealType scaled_norm = sqrt ( aol::Sqr<RealType> ( X / 2 - x ) + aol::Sqr<RealType> ( X / 2 - y ) );
        
        // Get the stepsize exponents of the two pixels on the diagonal whose distance to (X/2,X/2) is closest to scaled_norm
        // and calculate the stepsize exponent at (x,y) as the linear interpolation of the two exponents
        scaled_norm /= sqrt ( 2.0 );
        RealType int_part = static_cast<int> ( scaled_norm );
        RealType frac_part = scaled_norm - int_part;
        
        int exp;
        if ( int_part < X / 2 )
          exp = static_cast<int> ( ( 1 - frac_part ) * smooth_stepsize_exp[j][int_part] + frac_part * smooth_stepsize_exp[j][int_part+1] );
        else
          exp = smooth_stepsize_exp[j][int_part];
        
        // Set the scale mask entry
        ExitWaveScaleMask[j].set ( x, y, Param.scaleMask_exitwave * pow ( 2.0, exp / 6.0 ) );
      }
  
  cerr << "\rUpdating the exit wave scale mask (circular shape)... done.               " << endl;
}

#endif  // EWR_PRECONDITIONER_H
