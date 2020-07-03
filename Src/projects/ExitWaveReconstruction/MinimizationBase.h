#ifndef EWR_MINIMIZATIONBASE_H
#define EWR_MINIMIZATIONBASE_H

#include "ReconstructionBase.h"
#include "Functional.h"
#include "FunctionalDerivative.h"
/// \cond
#include "StepSaver.h"

#include <gradientDescent.h>
#include <Newton.h>
/// \endcond

/**
 * \brief Calculates the current phase id and remaining iterations of the current phase from a given total iteration count
 * 
 * \param [out] phase_id the current phase id
 * 
 * \param [out] phase_iterations the remaining iterations of the current phase
 * 
 * \param [in] current_iteration the current total iteration count
 * 
 * \param [in] Param microscope parameters and minimization settings
 */
template <typename RealType>
void getPhaseAndIterations ( int& phase_id,
                             int& phase_iterations,
                             const int current_iteration,
                             const Parameters<RealType>& Param ) {
  if ( Param.AlternatingMinimizationSteps.size ( ) == 1 ) {
    // No alternating minimization
    phase_id = 0;
    if ( Param.enableIntegrationDomainFitting )
      phase_iterations = Param.IntegrationDomainFittingSteps;
    else if ( Param.enableEWCircularScaleMask )
      phase_iterations = Param.EWCircularScaleMaskUpdateSteps;
    else
      phase_iterations = Param.max_iterations;
    
    const int current_iteration_mod = current_iteration % phase_iterations;
    phase_iterations -= current_iteration_mod;
  } else {
    // Alternating minimization
    int sum = 0;
    for ( unsigned int i = 0; i < Param.AlternatingMinimizationSteps.size ( ) ; i++ )
      sum += Param.AlternatingMinimizationSteps[i];
    const int current_iteration_mod = current_iteration % sum;
    
    int it_count = 0;
    for ( unsigned int i = 0; i < Param.AlternatingMinimizationSteps.size ( ) ; i++ ) {
      it_count += Param.AlternatingMinimizationSteps[i];
      if ( current_iteration_mod < it_count ) {
        phase_id = i;
        phase_iterations = it_count - current_iteration_mod;
        break;
      }
    }
  }
  
  if ( current_iteration + phase_iterations > Param.max_iterations )
    phase_iterations = Param.max_iterations - current_iteration;
}

/**
 * \brief Prints some information about the current optimization phase to cerr
 * 
 * \param [in] phase_id the current phase id
 * 
 * \param [in] phase_iterations number of iterations during the current phase
 * 
 * \param [in] Param microscope parameters and minimization settings
 */
template <typename RealType>
void printPhaseInfo ( const int phase_id,
                      const int phase_iterations,
                      const Parameters<RealType>& Param ) {
  vector<string> optArg;
  if ( Param.b_opt_exitwave[ phase_id ] ) optArg.push_back ( "exit wave" );
  if ( Param.b_opt_translation[ phase_id ] ) optArg.push_back ( "translations" );
  if ( Param.b_opt_focus[ phase_id ] ) optArg.push_back ( "focus values" );
  
  cerr << "Phase " << phase_id << " (" << phase_iterations << " iterations): optimizing ";
  
  if ( optArg.size ( ) == 0 ) {
    cerr << "---" << endl << endl;
    return;
  }
  
  if ( optArg.size ( ) == 1 ) {
    cerr << optArg[0] << "..." << endl << endl;
    return;
  }
  
  cerr << optArg[0];
  for ( unsigned int i = 1; i < optArg.size ( ) - 1 ; i++ )
    cerr << ", " << optArg[i];
  cerr << " and " << optArg[ optArg.size ( ) - 1 ] << "..." << endl << endl;
}

/**
 * \brief Convenience function for the creation of a minimization algorithm operator according to the setting Param.minAlg
 * 
 * \returns a pointer to the minimization algorithm operator
 */
template <typename RealType>
aol::Op<Arguments<RealType>, Arguments<RealType>>* CreateNewMinimizationAlgorithm (
          const Functional<RealType>& energy,
          const FunctionalDerivative<RealType>& derivative,
          const Parameters<RealType>& Param,
          const StepSaver<RealType>& stepSaver,
          const int maxIterations ) {
  switch ( Param.minAlg ) {
    case 0: { // Gradient descent
      auto alg = new aol::GridlessGradientDescent<RealType, Arguments<RealType>> ( energy, derivative, maxIterations, 1, Param.stop_epsilon );
      alg->setStepSaverReference ( stepSaver );
      return alg;
    } break;
    case 1: { // Conjugate gradient method
      auto alg = new aol::GridlessGradientDescent<RealType, Arguments<RealType>> ( energy, derivative, maxIterations, 1, Param.stop_epsilon );
      alg->setConfigurationFlags ( aol::GridlessGradientDescent<RealType, Arguments<RealType>>::USE_NONLINEAR_CG );
      alg->setStepSaverReference ( stepSaver );
      return alg;
    }  break;
    case 2: { // Quasi-Newton (BFGS)
      auto alg = new aol::QuasiNewtonBFGS<RealType, Arguments<RealType>, Arguments<RealType>> ( energy, derivative, maxIterations, Param.stop_epsilon );
      alg->setStepSaverReference ( stepSaver );
      return alg;
    } break;
    case 3: { // Gauss-Newton (LSMR)
      throw aol::Exception ( "Not implemented yet!", __FILE__, __LINE__ );
    } break;
    case 4: { // Gauss-Newton (QR)
      throw aol::Exception ( "Not implemented yet!", __FILE__, __LINE__ );
    } break;
    default:
      throw aol::Exception ( "Invalid minimization algorithm! (Param.minAlg)", __FILE__, __LINE__ );
      break;
  }
  
  return nullptr;
}

/**
 * \param [in] MinimizationAlgorithm the minimization algorithm operator
 * 
 * \param [in] Param microscope parameters and minization settings
 * 
 * \returns The number of iterations that have been performed with the last evaluation of the minimization operator
 */
template <typename RealType>
int numIterationsPerformed ( const aol::Op<Arguments<RealType>, Arguments<RealType>>& MinimizationAlgorithm,
                             const Parameters<RealType>& Param ) {
  int result = 0;
  
  switch ( Param.minAlg ) {
    case 0: // Gradient descent
    case 1: // Conjugate gradient method
      result = dynamic_cast< const aol::GridlessGradientDescent<RealType, Arguments<RealType>>& > ( MinimizationAlgorithm ).getIterations ( );
      break;
    case 2: // Quasi-Newton (BFGS)
      result = dynamic_cast< const aol::QuasiNewtonBFGS<RealType, Arguments<RealType>, Arguments<RealType>>& > ( MinimizationAlgorithm ).getNewtonInfo ( ).getIterationCount ( );
      break;
    case 3: // Gauss-Newton (LSMR)
      //
      break;
    case 4: // Gauss-Newton (QR)
      //
      break;
    default:
      throw aol::Exception ( "Invalid minimization algorithm! (Param.minAlg)", __FILE__, __LINE__ );
      break;
  }
  
  return result;
}

#endif  // EWR_MINIMIZATIONBASE_H
