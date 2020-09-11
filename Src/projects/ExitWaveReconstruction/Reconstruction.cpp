// Floating point data type that is used for the reconstruction
// (currently only float and double are supported)
typedef double RealType;

/// \cond
#include "MinimizationBase.h"
/// \endcond
#include "Functional.h"
#include "FunctionalDerivative.h"
#include "Preconditioner.h"

int main ( int argc, char **argv ) {
  if ( argc != 2 && argc != 3 ) {
    PrintProgramUsage ( cerr, argv[0] );
    return 0;
  }
  
  UserPrompts default_action = UserPrompts::Ask;
  if ( argc == 3 ) {
    if ( string ( argv[2] ) == "-y" )
      default_action = UserPrompts::Accept;
    else if ( string ( argv[2] ) == "-n" )
      default_action = UserPrompts::Decline;
    else {
      PrintProgramUsage ( cerr, argv[0] );
      return 0;
    }
  }
  
  // Print time and floating point data type
  PrintTime ( cerr );
  PrintRealType<RealType> ( cerr );
  
  /********************/
  /** Initialization **/
  /********************/
  cerr << endl
       << "********************" << endl
       << "** Initialization **" << endl
       << "********************" << endl
       << endl;
  
  // Load the parameter file
  Parameters<RealType> Param ( argv[1] );
  
  // Create the reconstruction directory path (all output is written to reconstructionDir)
  auto pos = string ( argv[1] ).find_last_of ( '/' );
  const string path = ( pos == string::npos ? string ( ) : string ( argv[1] ).substr ( 0, pos + 1 ) );
  
  auto pos2 = string ( argv[1] ).find_last_of ( '.' );
  auto substr_len = ( ( pos2 < ( pos == string::npos ? 0 : pos ) || pos2 == string::npos ) ? ( string::npos ) : ( pos2 - ( pos == string::npos ? 0 : ( pos + 1 ) ) ) );
  const string local_reconstructionDir = "Reconstruction_" + string ( argv[1] ).substr ( ( pos == string::npos ? 0 : ( pos + 1 ) ),  substr_len );
  
  const string reconstructionDir = path + local_reconstructionDir + '/';
  
  // Check for the existence of files from a previous reconstruction
  bool resumeReconstruction = false;
  if ( aol::directoryExists ( reconstructionDir ) ) {
    cerr << "Found data from a previous reconstruction (" << local_reconstructionDir << "). Resume the reconstruction? [y/n]";
    
    char answer = ' ';
    if ( default_action == UserPrompts::Ask )
      cin >> answer;
    if ( default_action == UserPrompts::Accept || answer == 'y' || answer == 'Y' ) {
      cerr << "Resuming the reconstruction..." << endl << endl;
      resumeReconstruction = true;
      Param = Parameters<RealType> ( reconstructionDir + "InputData/ParameterFile.param" );
    } else {
      cerr << "Program stopped." << endl;
      return 0;
    }
  }
  
  // Create the scale mask
  Arguments<RealType> ScaleMask ( Param.subsection_width + 2 * Param.reconstructionBufferZone,
                                  Param.subsection_height + 2 * Param.reconstructionBufferZone,
                                  Param.N );
  CreateScaleMask ( ScaleMask, Param, resumeReconstruction, reconstructionDir );
  
  // Generate and save the input data for the reconstruction, "input", which consists of
  // the image series, the aberration coefficients and several other parameters
  //
  // If requested in the parameter file, the image series is simulated from an exit wave
  // that is itself simulated or load from file. In this case a cropped version of this
  // exit wave of the same size as the images in the input image series is stored in
  // the "inputArg" variable (along with the other arguments, i.e. the correct
  // translations and focus values).
  InputData<RealType> input;
  Arguments<RealType> inputArg ( Param.subsection_width, Param.subsection_height, Param.N );
  
  CreateInputData ( input, inputArg, Param, resumeReconstruction, reconstructionDir, default_action, true );
  SaveInputData ( input, inputArg, reconstructionDir );
  Param = input.Param;
  
  // Create and save the initial guess
  Arguments<RealType> InitialGuess ( Param.X, Param.Y, Param.N );
  
  CreateInitialGuess ( InitialGuess, input, inputArg, ScaleMask, resumeReconstruction, reconstructionDir );
  SaveInitialGuess ( InitialGuess, ScaleMask, reconstructionDir );
  
  // Change the format of inputArg so that the functional can be evaluated at inputArg
  AddReconstructionBufferZone ( inputArg, Param );
  inputArg /= ScaleMask;
  inputArg.Translation *= -1;
  
  SaveScaledInputArg ( inputArg, ScaleMask, reconstructionDir );
  
  // Add the reconstruction buffer zone to the focus series and adapt the parameters
  // X, Y, lenX and lenY accordingly
  AddReconstructionBufferZone ( input );
  Param = input.Param;
  
  // Prepare OpenCL kernels for computations on the GPU(s)
  KernelData<RealType> kData ( Param );
  
  // Create the integration domains
  IntegrationDomains domains;
  GetIntegrationDomains ( domains, InitialGuess, ScaleMask, Param, resumeReconstruction, reconstructionDir );
  
  if ( !resumeReconstruction )
    domains.save ( reconstructionDir + "InitialGuess/.integration_domains" );
  
  // Create the energy and derivative operators
  Functional<RealType> energy ( input, domains, ScaleMask, kData );
  FunctionalDerivative<RealType> derivative ( input, domains, ScaleMask, kData );
  
  // Create energy and derivative operators for the scale mask updates
  InputData<RealType> inputPreconditioner ( input );
  inputPreconditioner.Param.simulationMode = inputPreconditioner.Param.EWCircularScaleMaskComputationMode;
  
  Arguments<RealType> ScaleMaskPreconditioner ( ScaleMask, aol::STRUCT_COPY );
  CreateScaleMask ( ScaleMaskPreconditioner, Param, false, "" );
  
  Functional<RealType> energyPreconditioner ( inputPreconditioner, domains, ScaleMaskPreconditioner, kData );
  FunctionalDerivative<RealType> derivativePreconditioner ( inputPreconditioner, domains, ScaleMaskPreconditioner, kData );
  
  // Create the step saver
  StepSaver<RealType> stepSaver ( energy, Param, inputArg, ScaleMask, resumeReconstruction, reconstructionDir, domains );
  
  /*******************************/
  /** Print general information **/
  /*******************************/
  cerr << endl
       << "*************************" << endl
       << "** General information **" << endl
       << "*************************" << endl
       << endl;
  
  // Print all parameters
  Param.print ( cerr, true );
  
  // Print the aperture radius in pixel
  RealType aperture_radius = ApertureRadius ( Param );
  if ( aperture_radius != -1 )
    cerr << "Aperture radius: " << aperture_radius << " (pixel)" << endl << endl;
  
  // Print the estimated translations and compare them with the correct translations (if available)
  if ( !inputArg.ExitWave.isZero ( ) )
    PrintTranslations ( InitialGuess.Translation, ScaleMask.Translation, Param, inputArg.Translation );
  else
    PrintTranslations ( InitialGuess.Translation, ScaleMask.Translation, Param );
  
  // Print the estimated focus values and compare them with the correct values (if available)
  if ( !inputArg.ExitWave.isZero ( ) )
    PrintFocusValues ( InitialGuess.Focus, ScaleMask.Focus, inputArg.Focus );
  else
    PrintFocusValues ( InitialGuess.Focus, ScaleMask.Focus );
  
  // Print the integration domains
  domains.print ( cerr );
  
  // Print the energy of the current estimate
  if ( resumeReconstruction )
    cerr << "Calculating the energy for the current estimate... ";
  else
    cerr << "Calculating the energy for the initial guess... ";
  
  aol::Scalar<RealType> res;
  energy.apply ( InitialGuess, res );
  cerr << res << endl;
  
  // Print the energy of inputArg (if available)
  if ( !inputArg.ExitWave.isZero ( ) ) {
    cerr << "Calculating the energy for the (cropped) arguments used to simulate the input data... ";
    
    energy.apply ( inputArg, res );
    cerr << res << endl;
  }
  
  /******************/
  /** Minimization **/
  /******************/
  cerr << endl
       << "******************" << endl
       << "** Minimization **" << endl
       << "******************" << endl
       << endl;
  
  // Perform the minimization
  aol::Op<Arguments<RealType>, Arguments<RealType>> *MinimizationAlgorithm = nullptr;
  
  int current_iteration = stepSaver.getIterationOffset ( );
  vector<bool> min_stopped ( Param.max_phase_id + 1, false );
  while ( current_iteration < Param.max_iterations ) {
    // Update the scale mask if necessary
    if ( Param.enableEWCircularScaleMask && current_iteration % Param.EWCircularScaleMaskUpdateSteps == 0 ) {
      Arguments<RealType> newScaleMask ( ScaleMask );
      
      Arguments<RealType> CurrentPos ( InitialGuess );
      CurrentPos *= ScaleMask;
      CurrentPos /= ScaleMaskPreconditioner;
      
      CreateCircularExitWaveScaleMask ( newScaleMask.ExitWave, CurrentPos, energyPreconditioner, derivativePreconditioner, Param );
    
      InitialGuess *= ScaleMask;
      InitialGuess /= newScaleMask;
      ScaleMask = newScaleMask;
      
      energy.setScaleMask ( ScaleMask );
      derivative.setScaleMask ( ScaleMask );
      stepSaver.setScaleMask ( ScaleMask );
    }
    
    // Calculate the current phase id and the remaining number of iterations of the
    // current phase
    int phase_id = 0;
    int phase_iterations = 0;
    getPhaseAndIterations ( phase_id, phase_iterations, current_iteration, Param );
    
    printPhaseInfo ( phase_id, phase_iterations, Param );
    
    derivative.setOptimizationPhase ( phase_id );
    MinimizationAlgorithm = CreateNewMinimizationAlgorithm ( energy, derivative, Param, stepSaver, phase_iterations );
    
    // Perform phase_iterations steps of the minimization algorithm
    MinimizationAlgorithm->applySingle ( InitialGuess );
    
    // Stop if the last update of all phases was smaller than stop_epsilon
    min_stopped[ phase_id ] = ( numIterationsPerformed ( *MinimizationAlgorithm, Param ) < phase_iterations );
      
    bool finished = true;
    for ( unsigned int i = 0; i < min_stopped.size ( ) ; i++ )
      if ( !min_stopped[i] ) {
        finished = false;
        break;
      }
    
    delete MinimizationAlgorithm;
    
    if ( finished ) {
      cerr << endl << "The updates of the energy are smaller than Param.stop_epsilon = "
           << Param.stop_epsilon << " for all phases. Stopped the minimization." << endl << endl;
      break;
    }
    
    current_iteration += phase_iterations;
    
    stepSaver.addToIterationOffset ( phase_iterations );
    stepSaver.setResumeReconstruction ( );
    
    // Update the integration domains if necessary
    if ( Param.enableIntegrationDomainFitting && current_iteration % Param.IntegrationDomainFittingSteps == 0 ) {
      vector<qc::ScalarArray<RealType, qc::QC_2D>> prevTranslations;
      stepSaver.getAndResetCachedTranslations ( prevTranslations );
      
      for ( unsigned int i = 0; i < prevTranslations.size ( ) ; i++ )
        prevTranslations[i] *= ScaleMask.Translation;
    
      UpdateIntegrationDomains ( domains, prevTranslations, Param, true );
      
      energy.setIntegrationDomains ( domains );
      derivative.setIntegrationDomains ( domains );
      stepSaver.setIntegrationDomains ( domains );
      
      energyPreconditioner.setIntegrationDomains ( domains );
      derivativePreconditioner.setIntegrationDomains ( domains );
    }
    
    if ( phase_id == Param.max_phase_id )
      cerr << endl << "--------------------------------------------------------------" << endl
           << "Iteration count: " << current_iteration << endl;
    
    cerr << endl << endl;
  }
  
  // Save the result
  stepSaver.writeLastTimeStep ( );
  
  return 0;
}
