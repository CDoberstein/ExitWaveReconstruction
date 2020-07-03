#ifndef EWR_PARAMETERCLASS_H
#define EWR_PARAMETERCLASS_H

#include <parameterParser.h>

template <typename RealType> class Parameters;
template <typename RealType> RealType getElectronWavelength ( const RealType AcceleratingVoltage );
template <typename RealType> string AlgorithmName ( const Parameters<RealType>& Param );

/**
 * \brief Container class that holds all parameters from the parameter file
 */
template <typename RealType>
class Parameters {
public:
  /******************************/
  /** General program settings **/
  /******************************/
  //! Maximum number of threads and GPUs
  int numThreads;
  
  //! Format of the output images
  int outputImageFormat;
  
  //! Format of the output image numbering
  string outputNumberingFormat;
  
  //! TEM image simulation mode
  int simulationMode;
  
  //! Input image series simulation mode
  int inputSimulationMode;
  
  //! Focal integration parameter
  int FocalIntegration_M;
  //! Focal integration parameter
  RealType FocalIntegration_delta;
  
  //! Minimization algorithm index
  int minAlg;
  
  //! Input data source (0 = simulated, 1 = read from file)
  int inputDataSource;
  
  //! Input data file path (if inputDataSource == 1)
  string inputDataFile;
  
  //! Iterations for the computation
  int collectDataIterations;
  //! Iterations for the saving of intermediate results
  int saveDataIterations;
  
  /***************************/
  /** Instrument parameters **/
  /***************************/
  //! Accelerating voltage of the electrons in kV
  RealType AcceleratingVoltage;
  
  //! Maximum semiangle allowed by the objective aperture in radians
  RealType alpha_max;
  
  //! Half angle of beam convergence in radians
  RealType alpha;
  
  //! Focal spread parameter (Delta) in nanometers
  RealType FocalSpread;
  
  /*****************************/
  /** Aberration coefficients **/
  /*****************************/
  //! Spherical aberration coefficient Cs in nanometers
  RealType SphericalAberration;
  
  /********************/
  /** Image settings **/
  /********************/
  //! Number of images
  int N;
  
  //! Width and height of the images in pixels
  int X, Y;
  
  //! Width and height of the image section in nanometers
  RealType lenX, lenY;
  
  //! Focus values
  vector<RealType> Focus;
  
  //! Subsection of the images used for the reconstruction
  int subsection_x, subsection_y;
  //! Subsection of the images used for the reconstruction
  int subsection_width, subsection_height;
  
  /************************************/
  /** Input data simulation settings **/
  /************************************/
  //! Specimen drift mode
  int specimenDriftMode;
  
  //! Specimen drift between successive images in nanometers (if specimenDriftMode == 0)
  RealType specimenDriftX, specimenDriftY;
  
  //! Specimen position for each image (if specimenDriftMode == 1)
  vector<RealType> specimenPosX, specimenPosY;
  
  //! Specimen rotation in degrees
  RealType specimenRotation;
  
  //! Size of the buffer zone around the simulated input images in pixels
  int inputDataBufferZone;
  
  /***********************/
  /** Specimen settings **/
  /***********************/
  //! Specimen flag
  int specimenFlag;
  
  //! Number of atoms or point charges in the visible area in each direction
  int specimenPointsX, specimenPointsY;
  
  //! Parameters for the honeycomb structure
  RealType specimenHoneycombParam;
  
  //! Offset of the specimen in nanometers
  RealType specimenOffsetX, specimenOffsetY;
  
  //! Atomic number of the simulated atoms
  int specimenAtomicNumber;
  
  //! Path to the potential coefficients parameter file
  string specimenPotentialCoefficientsFile;
  
  /*****************************/
  /** Reconstruction settings **/
  /*****************************/
  //! Maximum number of iterations
  int max_iterations;
  
  //! Stopping criterion epsilon
  RealType stop_epsilon;
  
  //! Scalar coefficient of the Tikhonov regularizer for the exit wave
  RealType tikhonov_coeff;
  
  //! Scale mask factors
  RealType scaleMask_exitwave, scaleMask_translation, scaleMask_focus;
  
  //! Parameter for the circular shaped exit wave scale mask
  bool enableEWCircularScaleMask;
  //! Parameter for the circular shaped exit wave scale mask
  int EWCircularScaleMaskUpdateSteps;
  //! Parameter for the circular shaped exit wave scale mask
  int EWCircularScaleMaskComputationMode;
  
  //! Size of the reconstruction buffer zone in pixel
  int reconstructionBufferZone;
  
  //! Parameter for the automatical adjustment of the integration domains
  bool enableIntegrationDomainFitting;
  //! Parameter for the automatical adjustment of the integration domains
  int IntegrationDomainFittingSteps;
  
  //! Parameters that determine which arguments are optimized at what time
  vector<int> opt_exitwave, opt_translation, opt_focus;
  vector<int> AlternatingMinimizationSteps;
  
  /*******************/
  /** Initial guess **/
  /*******************/
  //! Exit wave intial guess
  int initialGuess_ExitWave;
  
  //! Translation initial guess
  int initialGuess_Translation;
  //! Translation initial guess
  RealType initialGuess_TranslationMaxShift;
  //! Translation initial guess
  vector<int> initialGuess_TranslationHint_Img;
  //! Translation initial guess
  vector<int> initialGuess_TranslationHint_ShiftX;
  //! Translation initial guess
  vector<int> initialGuess_TranslationHint_ShiftY;
  //! Translation initial guess
  RealType initialGuess_TranslationStdDev;
  //! Translation initial guess
  RealType initialGuess_TranslationErrorBound;
  
  //! Base focus value
  int initialGuess_Focus;
  //! Base focus value
  RealType initialGuess_FocusStdDev;
  
  /***************************/
  /** Additional parameters **/
  /***************************/
  //! Electron wavelength in nm
  RealType lambda;
  
  //! Modified version of the parameters that determine which arguments are optimized
  int max_phase_id;
  //! Modified version of the parameters that determine which arguments are optimized
  vector<bool> b_opt_exitwave, b_opt_translation, b_opt_focus;
  
  //! Path to the parameter file
  string parameterfile;
  
public:
  Parameters ( ) = default;
  
  /**
   * \brief Reads the parameters from a file located at path using aol::ParameterParser
   */
  Parameters ( const string& path ) : parameterfile ( path ) {
    aol::ParameterParser parser ( path );
    
    // Read data from parameter file
    numThreads = parser.getInt ( "numThreads" );
    outputImageFormat = parser.getInt ( "outputImageFormat" );
    outputNumberingFormat = parser.getString ( "outputNumberingFormat" );
    simulationMode = parser.getInt ( "simulationMode" );
    inputSimulationMode = parser.getInt ( "inputSimulationMode" );
    FocalIntegration_M = parser.getInt ( "FocalIntegration_M" );
    FocalIntegration_delta = parser.getReal<RealType> ( "FocalIntegration_delta" );
    minAlg = parser.getInt ( "minAlg" );
    inputDataSource = parser.getInt ( "inputDataSource" );
    inputDataFile = aol::expandTildeOrEnvVar ( parser.getString ( "inputDataFile" ) );
    collectDataIterations = parser.getInt ( "collectDataIterations" );
    saveDataIterations = parser.getInt ( "saveDataIterations" );
    
    AcceleratingVoltage = parser.getReal<RealType> ( "AcceleratingVoltage" );
    alpha_max = parser.getReal<RealType> ( "alpha_max" );
    alpha = parser.getReal<RealType> ( "alpha" );
    FocalSpread = parser.getReal<RealType> ( "FocalSpread" );
    
    SphericalAberration = parser.getReal<RealType> ( "SphericalAberration" );
    
    N = parser.getInt ( "N" );
    X = parser.getInt ( "X" );
    Y = parser.getInt ( "Y" );
    lenX = parser.getReal<RealType> ( "lenX" );
    lenY = parser.getReal<RealType> ( "lenY" );
    aol::Vector<RealType> aol_Focus;
    parser.getRealVec ( "Focus", aol_Focus );
    subsection_x = parser.getInt ( "subsection_x" );
    subsection_y = parser.getInt ( "subsection_y" );
    subsection_width = parser.getInt ( "subsection_width" );
    subsection_height = parser.getInt ( "subsection_height" );
    
    specimenDriftMode = parser.getInt ( "specimenDriftMode" );
    specimenDriftX = parser.getReal<RealType> ( "specimenDriftX" );
    specimenDriftY = parser.getReal<RealType> ( "specimenDriftY" );
    aol::Vector<RealType> aol_specimenPosX, aol_specimenPosY;
    parser.getRealVec( "specimenPosX", aol_specimenPosX );
    parser.getRealVec( "specimenPosY", aol_specimenPosY );
    specimenRotation = parser.getReal<RealType> ( "specimenRotation" );
    inputDataBufferZone = parser.getInt ( "inputDataBufferZone" );
    
    specimenFlag = parser.getInt ( "specimenFlag" );
    specimenPointsX = parser.getInt ( "specimenPointsX" );
    specimenPointsY = parser.getInt ( "specimenPointsY" );
    specimenHoneycombParam = parser.getReal<RealType> ( "specimenHoneycombParam" );
    specimenOffsetX = parser.getReal<RealType> ( "specimenOffsetX" );
    specimenOffsetY = parser.getReal<RealType> ( "specimenOffsetY" );
    specimenAtomicNumber = parser.getInt ( "specimenAtomicNumber" );
    specimenPotentialCoefficientsFile = aol::expandTildeOrEnvVar ( parser.getString ( "specimenPotentialCoefficientsFile" ) );
    
    max_iterations = parser.getInt ( "max_iterations" );
    stop_epsilon = parser.getReal<RealType> ( "stop_epsilon" );
    tikhonov_coeff = parser.getReal<RealType> ( "tikhonov_coeff" );
    scaleMask_exitwave = parser.getReal<RealType> ( "scaleMask_exitwave" );
    scaleMask_translation = parser.getReal<RealType> ( "scaleMask_translation" );
    scaleMask_focus = parser.getReal<RealType> ( "scaleMask_focus" );
    enableEWCircularScaleMask = parser.getBool ( "enableEWCircularScaleMask" );
    EWCircularScaleMaskUpdateSteps = parser.getInt ( "EWCircularScaleMaskUpdateSteps" );
    EWCircularScaleMaskComputationMode = parser.getInt ( "EWCircularScaleMaskComputationMode" );
    reconstructionBufferZone = parser.getInt ( "reconstructionBufferZone" );
    enableIntegrationDomainFitting = parser.getBool ( "enableIntegrationDomainFitting" );
    IntegrationDomainFittingSteps = parser.getInt ( "IntegrationDomainFittingSteps" );
    aol::Vector<int> aol_opt_exitwave, aol_opt_translation, aol_opt_focus;
    aol::Vector<int> aol_AlternatingMinimizationSteps;
    parser.getIntVec ( "opt_exitwave", aol_opt_exitwave );
    parser.getIntVec ( "opt_translation", aol_opt_translation );
    parser.getIntVec ( "opt_focus", aol_opt_focus );
    parser.getIntVec ( "AlternatingMinimizationSteps", aol_AlternatingMinimizationSteps );
    
    initialGuess_ExitWave = parser.getInt ( "initialGuess_ExitWave" );
    initialGuess_Translation = parser.getInt ( "initialGuess_Translation" );
    initialGuess_TranslationMaxShift = parser.getReal<RealType> ( "initialGuess_TranslationMaxShift" );
    aol::Vector<int> aol_initialGuess_TranslationHint_Img, aol_initialGuess_TranslationHint_ShiftX, aol_initialGuess_TranslationHint_ShiftY;
    parser.getIntVec ( "initialGuess_TranslationHint_Img", aol_initialGuess_TranslationHint_Img );
    parser.getIntVec ( "initialGuess_TranslationHint_ShiftX", aol_initialGuess_TranslationHint_ShiftX );
    parser.getIntVec ( "initialGuess_TranslationHint_ShiftY", aol_initialGuess_TranslationHint_ShiftY );
    initialGuess_TranslationStdDev = parser.getReal<RealType> ( "initialGuess_TranslationStdDev" );
    initialGuess_TranslationErrorBound = parser.getReal<RealType> ( "initialGuess_TranslationErrorBound" );
    initialGuess_Focus = parser.getInt ( "initialGuess_Focus" );
    initialGuess_FocusStdDev = parser.getReal<RealType> ( "initialGuess_FocusStdDev" );
    
    // Calculate the electron wavelength
    lambda = getElectronWavelength ( 1000 * AcceleratingVoltage );
    
    // Limit the number of threads by the number of images
    numThreads = min ( numThreads, N );
    
    // Adapt the path of the input data file if it is a relative path
    auto pos = path.find_last_of ( '/' );
    if ( inputDataFile[0] != '/' && pos != string::npos )
      inputDataFile = path.substr ( 0, pos + 1 ) + inputDataFile;
    
    // Adapt the path of the potential coefficients parameter file if it is a relative path
    if ( specimenPotentialCoefficientsFile[0] != '/' && pos != string::npos )
      specimenPotentialCoefficientsFile = path.substr ( 0, pos + 1 ) + specimenPotentialCoefficientsFile;
    
    // Check the value of outputImageFormat
    if ( outputImageFormat != 0 && outputImageFormat != 1 )
      throw aol::Exception ( "Invalid format of the output images specified! (Param.outputImageFormat)", __FILE__, __LINE__ );
    
    // Convert the translation hint vectors
    if ( aol_initialGuess_TranslationHint_Img.size ( ) != aol_initialGuess_TranslationHint_ShiftX.size ( ) ||
         aol_initialGuess_TranslationHint_Img.size ( ) != aol_initialGuess_TranslationHint_ShiftY.size ( ) )
      throw aol::Exception ( "Invalid translation hint vector provided! (Param.initialGuess_TranslationHint_*)", __FILE__, __LINE__ );
    
    const int num_hints = aol_initialGuess_TranslationHint_Img.size ( );
    initialGuess_TranslationHint_Img.resize ( num_hints );
    initialGuess_TranslationHint_ShiftX.resize ( num_hints );
    initialGuess_TranslationHint_ShiftY.resize ( num_hints );
    
    for ( int i = 0; i < num_hints ; i++ ) {
      initialGuess_TranslationHint_Img[i] = aol_initialGuess_TranslationHint_Img[i];
      initialGuess_TranslationHint_ShiftX[i] = aol_initialGuess_TranslationHint_ShiftX[i];
      initialGuess_TranslationHint_ShiftY[i] = aol_initialGuess_TranslationHint_ShiftY[i];
    }

    // Convert (and possibly extend) the focus vector
    if ( aol_Focus.size ( ) < 2 && aol_Focus.size ( ) != N )
      throw aol::Exception ( "Insufficient number of focus values provided!", __FILE__, __LINE__ );
    
    Focus.resize ( N );
    for ( int i = 0; i < min ( N, aol_Focus.size ( ) ) ; i++ )
      Focus[i] = aol_Focus[i];
    
    if ( aol_Focus.size ( ) < N ) {
      RealType diff = aol_Focus[ aol_Focus.size ( ) - 1 ] - aol_Focus[ aol_Focus.size ( ) - 2 ];
      for ( int i = aol_Focus.size ( ) ; i < N ; i++ )
        Focus[i] = Focus[i-1] + diff;
    }
    
    // Convert the aol_specimenPosX and aol_specimenPosY vectors and check their sizes if necessary
    if ( specimenDriftMode == 1 ) {
      if ( aol_specimenPosX.size ( ) != N || aol_specimenPosY.size ( ) != N )
        throw aol::Exception ( "Invalid number of specimen positions with specimenPosX or specimenPosY provided!", __FILE__, __LINE__ );
      
      specimenPosX.resize ( N );
      specimenPosY.resize ( N );
      for ( int i = 0; i < N ; i++ ) {
        specimenPosX[i] = aol_specimenPosX[i] - aol_specimenPosX[0];
        specimenPosY[i] = aol_specimenPosY[i] - aol_specimenPosY[0];
      }
    }
    
    // Bound the subsection width and height
    if ( subsection_width == -1 || subsection_width + subsection_x > X ) {
      subsection_width = X - subsection_x;
      if ( subsection_width % 2 == 1 )
        --subsection_width;
    }
    
    if ( subsection_height == -1 || subsection_height + subsection_y > Y ) {
      subsection_height = Y - subsection_y;
      if ( subsection_height % 2 == 1 )
        --subsection_height;
    }
    
    // Check if the image dimensions are even integers
    if ( X % 2 == 1 || Y % 2 == 1 || subsection_width % 2 == 1 || subsection_height % 2 == 1 )
      throw aol::Exception ( "The parameters X, Y, subsection_width and subsection_height must be even integers!", __FILE__, __LINE__ );
    
    // Process the parameters that determine which arguments are optimized
    vector<aol::Vector<int>*> aol_opt = { &aol_opt_exitwave, &aol_opt_translation, &aol_opt_focus };
    vector<vector<int>*>          opt = {     &opt_exitwave,     &opt_translation,     &opt_focus };
    vector<vector<bool>*>       b_opt = {   &b_opt_exitwave,   &b_opt_translation,   &b_opt_focus };
    
    bool any_opt = false;
    for ( unsigned int i = 0; i < aol_opt.size ( ) ; i++ )
      if ( aol_opt[i]->size ( ) > 0 ) {
        any_opt = true;
        break;
      }
    if ( !any_opt )
      throw aol::Exception ( "No optimization arguments specified! (Param.opt_*)", __FILE__, __LINE__ );
    
    for ( unsigned int i = 0; i < opt.size ( ) ; i++ ) {
      opt[i]->resize ( aol_opt[i]->size ( ) );
      for ( int j = 0; j < aol_opt[i]->size ( ) ; j++ )
        (*opt[i])[j] = (*aol_opt[i])[j];
    }
    
    AlternatingMinimizationSteps.resize ( aol_AlternatingMinimizationSteps.size ( ) );
    for ( int j = 0; j < aol_AlternatingMinimizationSteps.size ( ) ; j++ )
      AlternatingMinimizationSteps[j] = aol_AlternatingMinimizationSteps[j];
    
    max_phase_id = 0;
    for ( unsigned int i = 0; i < aol_opt.size ( ) ; i++ )
      if ( aol_opt[i]->size ( ) > 0 ) {
        max_phase_id = max ( max_phase_id, aol_opt[i]->getMaxValue ( ) );
        if ( aol_opt[i]->getMinValue ( ) < 0 )
          throw aol::Exception ( "The elements of the opt_* vectors must be nonnegative!", __FILE__, __LINE__ );
      }
    
    if ( max_phase_id == 0 ) {
      AlternatingMinimizationSteps.resize ( 1 );
      AlternatingMinimizationSteps[0] = max_iterations;
    }
    
    if ( static_cast<int> ( AlternatingMinimizationSteps.size ( ) ) != max_phase_id + 1 )
      throw aol::Exception ( ( "Invalid length of the AlternatingMinimizationSteps vector! (expected " + to_string ( max_phase_id + 1 ) + " elements)" ).c_str ( ), __FILE__, __LINE__ );
    
    for ( unsigned int i = 0; i < opt.size ( ) ; i++ ) {
      *b_opt[i] = vector<bool> ( max_phase_id + 1, false );
      for ( unsigned int j = 0; j < opt[i]->size ( ) ; j++ )
        (*b_opt[i])[ (*opt[i])[j] ] = true;
    }
    
    // Adjust the integration domain fitting parameters if necessary
    if ( IntegrationDomainFittingSteps < 5 )
      IntegrationDomainFittingSteps = 5;
    
    if ( max_phase_id > 0 ) {
      // Round up IntegrationDomainFittingSteps to the next multiple of the sum of the
      // elements of AlternatingMinimizationSteps
      int sum = 0;
      for ( unsigned int i = 0; i < AlternatingMinimizationSteps.size ( ) ; i++ )
        sum += AlternatingMinimizationSteps[i];
      
      IntegrationDomainFittingSteps = ( ( IntegrationDomainFittingSteps - 1 ) / sum + 1 ) * sum;
    }
    
    // Disable automatic integration domain fitting if the translations are not optimized
    if ( opt_translation.empty ( ) )
      enableIntegrationDomainFitting = false;
    
    // Adjust the circular exit wave scale mask parameters if necessary
    if ( max_phase_id > 0 ) {      
      // Round up EWCircularScaleMaskUpdateSteps to the next multiple of the sum of the
      // elements of AlternatingMinimizationSteps
      int sum = 0;
      for ( unsigned int i = 0; i < AlternatingMinimizationSteps.size ( ) ; i++ )
        sum += AlternatingMinimizationSteps[i];
      
      EWCircularScaleMaskUpdateSteps = ( ( EWCircularScaleMaskUpdateSteps - 1 ) / sum + 1 ) * sum;
    } else if ( enableIntegrationDomainFitting ) {
      // Round up EWCircularScaleMaskUpdateSteps to the next multiple of
      // IntegrationDomainFittingSteps
      EWCircularScaleMaskUpdateSteps = ( ( EWCircularScaleMaskUpdateSteps - 1 ) / IntegrationDomainFittingSteps + 1 ) * IntegrationDomainFittingSteps;
    }
    
    // Disable the circular exit wave scale mask if the exit wave is not optimized
    if ( opt_exitwave.empty ( ) )
      enableEWCircularScaleMask = false;
    
    // Since X, Y and N are implemented as ints, the results will likely be wrong if
    // their product exceeds the 32 bit int range, so we stop the program in this case
    if ( static_cast<long> ( X + 2 * inputDataBufferZone ) * static_cast<long> ( Y + 2 * inputDataBufferZone ) * static_cast<long> ( N ) > ( 1l << 30 ) ||
         static_cast<long> ( X + 2 * reconstructionBufferZone ) * static_cast<long> ( Y + 2 * reconstructionBufferZone ) * static_cast<long> ( N ) > ( 1l << 30 ) )
      throw aol::Exception ( "The images are too large!", __FILE__, __LINE__ );
  }
  
  /**
   * \brief Save the parameters to a file
   * 
   * The new parameter file created this way is only read if an existing
   * reconstruction is resumed.
   * 
   * \note the paths inputDataFile and specimenPotentialCoefficientsFile will likely be
   *       wrong in the new parameter file if they are relative paths to the original
   *       parameter file. However, since these two paths are only needed for the input
   *       data generation, they are not adapted here. (At this point it is expected that
   *       the input data generation is already finished. Therefore these two paths are
   *       not needed anymore.)
   */
  void save ( const string& filename ) const {
    ofstream paramfile ( filename );
    if ( !paramfile.is_open ( ) )
      throw aol::Exception ( ( "Unable to open \"" + filename + "\"!" ).c_str ( ), __FILE__, __LINE__ );
    
    paramfile << "## Automatically generated Parameter file. Do not edit ##" << endl << endl;
    
    paramfile << "# General program settings" << endl
              << "numThreads " << numThreads << endl
              << "outputImageFormat " << outputImageFormat << endl
              << "outputNumberingFormat \"" << outputNumberingFormat << "\"" << endl
              << "simulationMode " << simulationMode << endl
              << "inputSimulationMode " << inputSimulationMode << endl
              << "FocalIntegration_M " << FocalIntegration_M << endl
              << "FocalIntegration_delta " << FocalIntegration_delta << endl
              << "minAlg " << minAlg << endl
              << "inputDataSource " << inputDataSource << endl
              << "inputDataFile \"" << inputDataFile << "\"" << endl
              << "collectDataIterations " << collectDataIterations << endl
              << "saveDataIterations " << saveDataIterations << endl << endl;
    
    paramfile << "# Instrument parameters" << endl
              << "AcceleratingVoltage " << AcceleratingVoltage << endl
              << "alpha_max " << alpha_max << endl
              << "alpha " << alpha << endl
              << "FocalSpread " << FocalSpread << endl << endl;
    
    paramfile << "# Aberration coefficients" << endl
              << "SphericalAberration " << SphericalAberration << endl << endl;
    
    paramfile << "# Image settings" << endl
              << "N " << N << endl
              << "X " << X << endl
              << "Y " << Y << endl
              << "lenX " << lenX << endl
              << "lenY " << lenY << endl
              << "Focus "; printVec ( Focus, paramfile ); paramfile << endl
              << "subsection_x " << subsection_x << endl
              << "subsection_y " << subsection_y << endl
              << "subsection_width " << subsection_width << endl
              << "subsection_height " << subsection_height << endl << endl;
    
    paramfile << "# Input data simulation settings" << endl
              << "specimenDriftMode " << specimenDriftMode << endl
              << "specimenDriftX " << specimenDriftX << endl
              << "specimenDriftY " << specimenDriftY << endl
              << "specimenPosX "; printVec ( specimenPosX, paramfile ); paramfile << endl
              << "specimenPosY "; printVec ( specimenPosY, paramfile ); paramfile << endl
              << "specimenRotation " << specimenRotation << endl
              << "inputDataBufferZone " << inputDataBufferZone << endl << endl;
    
    paramfile << "# Specimen settings" << endl
              << "specimenFlag " << specimenFlag << endl
              << "specimenPointsX " << specimenPointsX << endl
              << "specimenPointsY " << specimenPointsY << endl
              << "specimenHoneycombParam " << specimenHoneycombParam << endl
              << "specimenOffsetX " << specimenOffsetX << endl
              << "specimenOffsetY " << specimenOffsetY << endl
              << "specimenAtomicNumber " << specimenAtomicNumber << endl
              << "specimenPotentialCoefficientsFile \"" << specimenPotentialCoefficientsFile << "\"" << endl << endl;
    
    paramfile << "# Reconstruction settings" << endl
              << "max_iterations " << max_iterations << endl
              << "stop_epsilon " << stop_epsilon << endl
              << "tikhonov_coeff " << tikhonov_coeff << endl
              << "scaleMask_exitwave " << scaleMask_exitwave << endl
              << "scaleMask_translation " << scaleMask_translation << endl
              << "scaleMask_focus " << scaleMask_focus << endl
              << "enableEWCircularScaleMask " << enableEWCircularScaleMask << endl
              << "EWCircularScaleMaskUpdateSteps " << EWCircularScaleMaskUpdateSteps << endl
              << "EWCircularScaleMaskComputationMode " << EWCircularScaleMaskComputationMode << endl
              << "reconstructionBufferZone " << reconstructionBufferZone << endl
              << "enableIntegrationDomainFitting " << enableIntegrationDomainFitting << endl
              << "IntegrationDomainFittingSteps " << IntegrationDomainFittingSteps << endl
              << "opt_exitwave "; printVec ( opt_exitwave, paramfile ); paramfile << endl
              << "opt_translation "; printVec ( opt_translation, paramfile ); paramfile << endl
              << "opt_focus "; printVec ( opt_focus, paramfile ); paramfile << endl
              << "AlternatingMinimizationSteps "; printVec ( AlternatingMinimizationSteps, paramfile ); paramfile << endl << endl;
    
    paramfile << "# Initial guess" << endl
              << "initialGuess_ExitWave " << initialGuess_ExitWave << endl
              << "initialGuess_Translation " << initialGuess_Translation << endl
              << "initialGuess_TranslationMaxShift " << initialGuess_TranslationMaxShift << endl
              << "initialGuess_TranslationHint_Img "; printVec ( initialGuess_TranslationHint_Img, paramfile ); paramfile << endl
              << "initialGuess_TranslationHint_ShiftX "; printVec ( initialGuess_TranslationHint_ShiftX, paramfile ); paramfile << endl
              << "initialGuess_TranslationHint_ShiftY "; printVec ( initialGuess_TranslationHint_ShiftY, paramfile ); paramfile << endl
              << "initialGuess_TranslationStdDev " << initialGuess_TranslationStdDev << endl
              << "initialGuess_TranslationErrorBound " << initialGuess_TranslationErrorBound << endl
              << "initialGuess_Focus " << initialGuess_Focus << endl
              << "initialGuess_FocusStdDev " << initialGuess_FocusStdDev << endl;
  }
  
  /**
   * \brief Prints the parameters to out
   */
  void print ( ostream& out,
               const bool reconstruction_buffer_zone_added = false ) const {
    out << "Parameter file path" << endl
        << "\t\"" << parameterfile << "\"" << endl << endl;
    
    out << "General program settings" << endl
        << "\tMaximum number of threads (and GPUs): " << numThreads << endl
        << "\tOutput image format: " << ( outputImageFormat == 0 ? "q2bz" : "tiff" ) << endl
        << "\tImage simulation mode: " << simulationMode << endl
        << "\tInput image series simulation mode: " << inputSimulationMode << endl
        << "\t\tFocal integration parameters: " << FocalIntegration_M << ", " << FocalIntegration_delta << " nm" << endl
        << "\tMinimization algorithm: " << minAlg << " (" << AlgorithmName ( *this ) << ")" << endl
        << "\tInput data source: " << ( inputDataSource == 0 ? "simulated" : "read from file" ) << endl
        << "\t\tInput data file path: \"" << inputDataFile << "\"" << endl
        << "\tIntermediate results are computed every " << collectDataIterations << " iterations and saved every " << saveDataIterations << " iterations" << endl << endl;
    
    out << "Instrument parameters" << endl
        << "\tAccelerating voltage: " << AcceleratingVoltage << " kV" << endl
        << "\tElectron wavelength: " << lambda << " nm" << endl
        << "\tObjective aperture semiangle: " << alpha_max << " rad" << endl
        << "\tHalf angle of beam convergence: " << alpha << " rad" << endl
        << "\tFocal spread: " << FocalSpread << " nm" << endl << endl;
    
    out << "Aberration coefficients" << endl
        << "\tSpherical aberration: " << SphericalAberration << " nm" << endl << endl;
    
    const int origX = ( reconstruction_buffer_zone_added ? X - 2 * reconstructionBufferZone : X );
    const int origY = ( reconstruction_buffer_zone_added ? Y - 2 * reconstructionBufferZone : Y );
    const RealType origLenX = ( lenX * origX ) / X;
    const RealType origLenY = ( lenY * origY ) / Y;
    
    out << "Image settings" << endl
        << "\tNumber of images: " << N << endl
        << "\tWidth and height: " << origX << " x " << origY << " (pixel)" << endl
        << "\tWidth and height: " << origLenX << " x " << origLenY << " (nm)" << endl
        << "\tFocus vector: "; printVec ( Focus, out, ',' ); out << " (nm)" << endl
        << "\tSubsection used for the reconstruction: (" << subsection_x << ", " << subsection_y << "), "
          << subsection_width << " x " << subsection_height << " (pixel)" << endl << endl;
    
    out << "Input data simulation settings" << endl
        << "\tSpecimen drift mode: " << specimenDriftMode << endl
        << "\tSpecimen drift: " << specimenDriftX << ", " << specimenDriftY << " (nm)" << endl
        << "\tSpecimen pos (X): "; printVec ( specimenPosX, out, ',' ); out << " (nm)" << endl
        << "\tSpecimen pos (Y): "; printVec ( specimenPosY, out, ',' ); out << " (nm)" << endl
        << "\tSpecimen rotation: " << specimenRotation << " (deg)" << endl
        << "\tBuffer zone size for input data simulation: " << inputDataBufferZone << " (pixel)" << endl << endl;
    
    if ( inputDataSource == 0 )
      out << "Specimen settings" << endl
          << "\tSpecimen flag: " << specimenFlag << endl
          << "\tNumber of atoms or point charges in each direction: " << specimenPointsX << ", " << specimenPointsY << endl
          << "\tHoneycomb structure parameter: " << specimenHoneycombParam << endl
          << "\tSpecimen offset: " << specimenOffsetX << ", " << specimenOffsetY << " (nm)" << endl
          << "\tAtomic number of simulated atoms: " << specimenAtomicNumber << endl
          << "\tPotential coefficients parameter file: \"" << specimenPotentialCoefficientsFile << "\"" << endl << endl;
    
    out << "Reconstruction settings" << endl
        << "\tMaximum number of iterations: " << max_iterations << endl
        << "\tStopping criterion epsilon: " << stop_epsilon << endl
        << "\tCoefficient of the Tikhonov regularizer: " << tikhonov_coeff << endl
        << "\tScale mask factors: " << scaleMask_exitwave << ", " << scaleMask_translation << ", " << scaleMask_focus << endl
        << "\tCircular exit wave scale mask: " << ( enableEWCircularScaleMask ? "yes" : "no" ) << " (updated every " << EWCircularScaleMaskUpdateSteps << " iterations)" << endl
        << "\tCircular exit wave scale mask computation mode: " << EWCircularScaleMaskComputationMode << endl
        << "\tBuffer zone size for the reconstruction: " << reconstructionBufferZone << " (pixel)" << endl
        << "\tAutomatic integration domain fitting: " << ( enableIntegrationDomainFitting ? "enabled" : "disabled" ) << " (every " << IntegrationDomainFittingSteps << " iterations)" << endl
        << "\tOptimization of the exit wave: "; printVec ( opt_exitwave, out, ',' ); out << endl
        << "\tOptimization of the translations: "; printVec ( opt_translation, out, ',' ); out << endl
        << "\tOptimization of the focus values: "; printVec ( opt_focus, out, ',' ); out << endl
        << "\tMinimization steps in each phase: "; printVec ( AlternatingMinimizationSteps, out, ',' ); out << endl << endl;
    
    out << "Initial guess" << endl
        << "\tExit wave: " << initialGuess_ExitWave << endl
        << "\tTranslation: " << initialGuess_Translation << " (Max. shift: " << initialGuess_TranslationMaxShift << " nm; std. dev.: " << initialGuess_TranslationStdDev << " nm)" << endl
        << "\tTranslation hints: "; printVec ( initialGuess_TranslationHint_Img, out, ',' ); out << endl
        << "\t                   "; printVec ( initialGuess_TranslationHint_ShiftX, out, ',' ); out << endl
        << "\t                   "; printVec ( initialGuess_TranslationHint_ShiftY, out, ',' ); out << endl
        << "\tUpper bound on the error of the translation initial guess: " << initialGuess_TranslationErrorBound << " (pixel)" << endl
        << "\tFocus: " << initialGuess_Focus << " (Std. dev.: " << initialGuess_FocusStdDev << " nm)" << endl << endl;
  }
  
  /**
   * \brief Auxiliary function to print a vector in a format that is compatible with aol::ParameterParser
   */
  template <typename T>
  static void printVec ( const vector<T>& v, ostream& out, const char value_sep = ' ' ) {
    if ( v.empty ( ) ) {
      out << "{ }";
      return;
    }
    
    out << "{ " << v[0];
    for ( unsigned int i = 1; i < v.size ( ) ; i++ )
      out << value_sep << " " << v[i];
    out << " }";
  }
};

#include "Physics.h"

#endif  // EWR_PARAMETERCLASS_H
