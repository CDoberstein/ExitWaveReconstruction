#ifndef EWR_RECONSTRUCTIONBASE_H
#define EWR_RECONSTRUCTIONBASE_H

/// \cond
#include "Image.h"
#include "Parameters.h"
/// \endcond
#include "SimulateInputData.h"
/// \cond
#include "Utilities.h"
#include "LoadUtilities.h"
/// \endcond
#include "OpenCLKernelData.h"
#include "Registration.h"

/// \cond
#include <randomGenerator.h>
/// \endcond

/**
 * \brief Wrapper class for the functional's arguments
 * 
 * The currently implemented optimization parameters are the exit wave, the translations
 * and the focus values.
 * 
 * \note the dimensions of the underlying MultiVector should not be changed.
 */
template <typename RealType>
class Arguments : public aol::MultiVector<RealType> {
public:
  struct DimRange {
    int X, Y, N;
    
    DimRange ( int _X, int _Y, int _N ) : X ( _X ), Y ( _Y ), N ( _N ) { }
  };

private:
  DimRange P;

public:
  //! Argument 1: Fourier space exit wave
  Image<RealType> ExitWave;
  
  //! Argument 2: translations in nm as a ScalarArray of size N x 2
  qc::ScalarArray<RealType, qc::QC_2D> Translation;
  
  //! Argument 3: Focus values in nm
  aol::Vector<RealType> Focus;
  
  /**
   * \brief Structure initialization constructor for an exit wave of size X x Y and N images
   */
  Arguments ( const int X, const int Y, const int N )
    : P ( X, Y, N ), ExitWave ( X, Y, Space::FourierSpace ), Translation ( N, 2 ), Focus( N ) {
    this->appendReference ( ExitWave[0] );
    this->appendReference ( ExitWave[1] );
    this->appendReference ( Translation );
    this->appendReference ( Focus );
  }
  
  Arguments ( const DimRange& d )
    : Arguments ( d.X, d.Y, d.N ) { }
  
  /**
   * \brief Constructs an Arguments class from an aol::MultiVector with an exit wave of size X x Y and N images
   */
  Arguments ( const aol::MultiVector<RealType>& mvec, const int X, const int Y, const int N )
    : Arguments ( X, Y, N ) {
    for ( int i = 0; i < mvec.numComponents ( ) ; i++ )
      (*this)[i] = mvec[i];
  }
  
  /**
   * \brief Copy constructor
   */
  Arguments ( const Arguments& Arg, aol::CopyFlag copyFlag = aol::DEEP_COPY )
    : Arguments ( Arg.P.X, Arg.P.Y, Arg.P.N ) {
    switch ( copyFlag ) {
      case aol::FLAT_COPY:
        throw aol::UnimplementedCodeException ( "Not implemented!", __FILE__, __LINE__ );
        break;
      case aol::DEEP_COPY:
        this->operator= ( Arg );
        break;
      case aol::STRUCT_COPY:
        // Nothing to do here
        break;
      case aol::STRUCT_COPY_UNINIT:
        throw aol::UnimplementedCodeException ( "Not implemented!", __FILE__, __LINE__ );
        break;
      default:
        throw aol::Exception ( "Invalid copy flag!", __FILE__, __LINE__ );
        break;
    }
  }
  
  /**
   * \brief Checks the validity of the arguments structure
   */
  void checkValidity ( const RealType Im0Focus ) const {
    // Check exit wave space and format
    if ( !ExitWave.correctFormat ( Space::FourierSpace, true ) )
      throw aol::Exception ( "Invalid exit wave format!", __FILE__, __LINE__ );
    
    // Check if the sizes of the ExitWave, Translation and Focus are
    // consistent with the underlying MultiVector
    const int num_components = 4;
    int components_size[4] = { ExitWave.getNumX ( ) * ExitWave.getNumY ( ),
                               ExitWave.getNumX ( ) * ExitWave.getNumY ( ),
                               Translation.getNumX ( ) * Translation.getNumY ( ),
                               Focus.size ( ) };
    
    if ( this->numComponents ( ) != num_components )
      throw aol::Exception ( "Invalid number of components!", __FILE__, __LINE__ );
    
    for ( int i = 0; i < num_components ; i++ )
      if ( this->operator[] ( i ).size ( ) != components_size[i] )
        throw aol::Exception ( "Invalid component size!", __FILE__, __LINE__ );
    
    // Check if the first translation is equal to zero
    if ( Translation.get ( 0, 0 ) != 0 || Translation.get ( 0, 1 ) != 0 )
      throw aol::Exception ( "The translation of the first image is not zero!", __FILE__, __LINE__ );
    
    // Check if the first focus value is kept constant
    if ( Focus[0] != Im0Focus )
      throw aol::Exception ( "The focus value of the first image was changed!", __FILE__, __LINE__ );
  }
  
  /**
   * \brief Saves the arguments in various human-readable formats
   */
  void saveArguments ( const string& dir,
                       const int output_image_format,
                       const Arguments<RealType>& subArg = Arguments<RealType> ( 0, 0, 0 )) const {
    // Exit wave
    Arguments<RealType> Arg ( *this );
    if ( subArg.ExitWave.getTotalSize ( ) != 0 )
      Arg -= subArg;
    
    SaveImage ( Arg.ExitWave, false, dir + "EW_FourierSpace", output_image_format  );
    
    Image<RealType> RealSpaceExitWave ( P.X, P.Y, Space::RealSpace );
    Arg.ExitWave.FourierTransformTo ( RealSpaceExitWave );
    SaveImage ( RealSpaceExitWave, false, dir + "EW_RealSpace", output_image_format );
    
    if ( subArg.ExitWave.getTotalSize ( ) == 0 ) {
      // No subtrahend: calculate the power spectrum and the amplitude and phase directly
      qc::MultiArray<RealType, 2, 2> APImage ( P.X, P.Y );
      Arg.ExitWave.CalculateAmplitudeAndPhase ( APImage );
      SaveImage ( APImage, true, dir + "EW_PowerSpectrum", output_image_format );
      
      RealSpaceExitWave.CalculateAmplitudeAndPhase ( APImage );
      SaveImage ( APImage, false, dir + "EW_RealSpace", output_image_format, true );
    } else {
      // Difference of two exit waves: calculate the power spectra, amplitudes and phases
      // first and the difference afterwards
      qc::MultiArray<RealType, 2, 2> APImage1 ( P.X, P.Y );
      qc::MultiArray<RealType, 2, 2> APImage2 ( P.X, P.Y );
      
      this->ExitWave.CalculateAmplitudeAndPhase ( APImage1 );
      subArg.ExitWave.CalculateAmplitudeAndPhase ( APImage2 );
      
      APImage1 -= APImage2;
      SaveImage ( APImage1, true, dir + "EW_PowerSpectrum", output_image_format );
      
      Image<RealType> RealSpaceEW1 ( P.X, P.Y, Space::RealSpace );
      Image<RealType> RealSpaceEW2 ( P.X, P.Y, Space::RealSpace );
      
      this->ExitWave.FourierTransformTo ( RealSpaceEW1 );
      subArg.ExitWave.FourierTransformTo ( RealSpaceEW2 );
      
      RealSpaceEW1.CalculateAmplitudeAndPhase ( APImage1 );
      RealSpaceEW2.CalculateAmplitudeAndPhase ( APImage2 );
      
      APImage1 -= APImage2;
      SaveImage ( APImage1, false, dir + "EW_RealSpace", output_image_format, true );
    }
    
    // Translation
    ofstream translationfile ( dir + "Translation.txt" );
    if ( translationfile.is_open ( ) ) {
      for ( int i = 0; i < Translation.getNumX ( ) ; i++ )
        translationfile << setw ( 3 ) << i << ": "
                        << setw ( 14 ) << setprecision ( 8 ) << Translation.get ( i, 0 ) << " nm, "
                        << setw ( 14 ) << setprecision ( 8 ) << Translation.get ( i, 1 ) << " nm" << endl;
    }
    
    // Focus values
    ofstream focusfile ( dir + "Focus.txt" );
    if ( focusfile.is_open ( ) ) {
      for ( int i = 0; i < Focus.size ( ) ; i++ )
        focusfile << setw ( 3 ) << i << ": "
                  << setw ( 14 ) << setprecision ( 8 ) << Focus[i] << " nm" << endl;
    }
  }
  
  /**
   * Overwrites the MultiVector::loadFromFile method since it does not work with FLAT_COPY
   * vectors (due to appendReference in the Arguments constructor)
   * 
   * \note *this must already have the correct size when calling loadFromFile
   */
  void loadFromFile ( const string& file ) {
    aol::MultiVector<RealType> tmp;
    tmp.loadFromFile ( file.c_str ( ) );
    *this = Arguments<RealType> ( tmp, P.X, P.Y, P.N );
  }
  
  void reallocateExitWave ( const int X, const int Y ) {
    ExitWave.reallocate ( X, Y );
    
    P.X = X;
    P.Y = Y;
  }
  
  /**
   * Removes a buffer zone of the size buffer_zone_size from this->ExitWave and stores
   * the result in Arg
   */
  void removeBufferZone ( Arguments<RealType>& Arg, const int buffer_zone_size ) {
    this->ExitWave.removeRealSpaceBufferZone ( Arg.ExitWave, buffer_zone_size );
    Arg.Translation = this->Translation;
    Arg.Focus = this->Focus;
  }
};

/**
 * \brief Constant input parameters for the functional
 * 
 * This container holds all the constant input data to the
 * functional. This includes the focus image series, aberration
 * coefficients and other TEM parameters.
 */
template <typename RealType>
class InputData {
public:
  //! Focus image series
  vector<Image<RealType>> ImageSeries;
  
  //! Instrument settings and other aberration coefficients
  Parameters<RealType> Param;
  
  /**
   * \brief Checks if the input data struct satisfies all constraints of the implementation
   */
  void checkValidity ( ) const {
    // Check if at least one image is given
    if ( ImageSeries.size ( ) == 0 )
      throw aol::Exception ( "Empty image series!", __FILE__, __LINE__ );
    
    // Check if all input images have the correct size
    for ( unsigned int i = 0; i < ImageSeries.size ( ) ; i++ )
      if ( ImageSeries[i].getNumX ( ) != Param.X || ImageSeries[i].getNumY ( ) != Param.Y )
        throw aol::Exception ( "Invalid image sizes!", __FILE__, __LINE__ );
    
    // Check if the sizes are smaller than or equal to 2^14 = 16384 to prevent overflows
    if ( Param.X > 16384 || Param.Y > 16384 )
      throw aol::Exception ( "Images are too large!", __FILE__, __LINE__ );
    
    // Check if all images are given in Fourier space with a frequency shift
    for ( unsigned int i = 1; i < ImageSeries.size ( ) ; i++ )
      if ( ImageSeries[i].getImageSpace ( ) != Space::FourierSpace || !ImageSeries[i].getFreqShift ( ) )
        throw aol::Exception ( "Invalid image format!", __FILE__, __LINE__ );
  }
  
  /**
   * \brief Checks compatibility with an Arguments structure
   */
  bool isCompatibleTo ( const Arguments<RealType>& Arg ) const {
    if ( Arg.ExitWave.getNumX ( ) != ImageSeries[0].getNumX ( ) || Arg.ExitWave.getNumY ( ) != ImageSeries[0].getNumY ( ) )
      return false;
    
    if ( Arg.Translation.getNumX ( ) != static_cast<int> ( ImageSeries.size ( ) ) )
      return false;
    
    if ( Arg.Focus.size ( ) != static_cast<int> ( ImageSeries.size ( ) ) )
      return false;
    
    return true;
  }
};

/**
 * \brief Initializes the scale mask with the Param.scaleMask settings
 * 
 * \param [in] resumeReconstruction true, if a previous reconstruction is resumed, and false otherwise
 * 
 * \param [in] reconstructionDir path to the base directory for the reconstruction
 */
template <typename RealType>
void CreateScaleMask ( Arguments<RealType>& ScaleMask,
                       const Parameters<RealType>& Param,
                       const bool resumeReconstruction,
                       const string& reconstructionDir,
                       const bool verbose = true ) {
  if ( resumeReconstruction ) {
    // An earlier reconstruction is resumed: load the most recent scale mask from file
    if ( verbose )
      cerr << "Loading the most recent scale mask from the previous reconstruction ..." << endl;
    
    ifstream fs ( reconstructionDir + "IntermediateResults/.lastresult" );
    if ( !fs.is_open ( ) ) {
      // The file likely could not be opened because it does not exist, which means that
      // less than Param.saveDataIterations minimization steps have been performed in the
      // last run; in this case the scale mask from the initial guess is loaded instead
      cerr << "Warning: unable to find the last estimate from the previous reconstruction. Loading the scale mask from the initial guess instead..." << endl;
      
      ScaleMask.loadFromFile ( reconstructionDir + "InitialGuess/.ScaleMask" );
    } else {
      string currentEstimateDir;
      getline ( fs, currentEstimateDir );
      ScaleMask.loadFromFile ( reconstructionDir + currentEstimateDir + ".ScaleMask" );
    }
  } else {
    ScaleMask.ExitWave.setAll ( Param.scaleMask_exitwave );
    ScaleMask.Translation.setAll ( Param.scaleMask_translation );
    ScaleMask.Focus.setAll ( Param.scaleMask_focus );
  }
}

/**
 * \brief Simulate input data or read input data from file if a previous reconstruction is resumed
 * 
 * \param [out] input the input data
 * 
 * \param [out] inputArg the arguments used to simulate the input data (if available)
 * 
 * \param [in] Param microscope settings
 * 
 * \param [in] resumeReconstruction true, if a previous reconstruction is resumed, and false otherwise
 * 
 * \param [in] reconstructionDir path to the base directory for the reconstruction
 */
template <typename RealType>
void CreateInputData ( InputData<RealType>& input,
                       Arguments<RealType>& inputArg,
                       const Parameters<RealType>& Param,
                       const bool resumeReconstruction,
                       const string& reconstructionDir,
                       bool verbose = true ) {
  if ( verbose )
    cerr << "Creating the input data..." << endl;
  
  // input: copy Param
  input.Param = Param;
  
  if ( resumeReconstruction ) {
    // Load the input data from the previous reconstruction...
    // ... read the focus image series from file
    LoadInputImageSeries ( input.ImageSeries, reconstructionDir + "InputData.dm3", Param, verbose );
    
    // ... also load inputArg if available
    if ( aol::fileExists ( reconstructionDir + "InputData/CroppedArg/.Arg" ) )
      inputArg.loadFromFile ( reconstructionDir + "InputData/CroppedArg/.Arg" );
  } else {
    string inputDataFileExt = Param.inputDataFile;
    if ( inputDataFileExt.length ( ) > 4 )
      inputDataFileExt = inputDataFileExt.substr ( inputDataFileExt.length ( ) - 4 );
    
    if ( Param.inputDataSource == 0 || inputDataFileExt == ".wav" ) {
      // Simulate input data...
      // ... large exit wave that initially has the size
      //    (Param.X+2*Param.inputDataBufferZone) x (Param.Y+2*Param.inputDataBufferZone),
      //   which is reduced by the rotation of the exit wave in SimulateExitWave and
      //   LoadExitWave
      Image<RealType> RealSpaceExitWave;
      if ( Param.inputDataSource == 0 )
        SimulateExitWave ( RealSpaceExitWave, Param, verbose );
      else
        LoadExitWave ( RealSpaceExitWave, Param, verbose );
      
      // Check if the rectangular image subsection given by the four subsection_*
      // parameters is still completely contained within RealSpaceExitWave after the
      // rotation in SimulateExitWave or LoadExitWave
      const int removedBufferX = ( Param.X + 2 * Param.inputDataBufferZone - RealSpaceExitWave.getNumX ( ) ) / 2;
      const int removedBufferY = ( Param.Y + 2 * Param.inputDataBufferZone - RealSpaceExitWave.getNumY ( ) ) / 2;
      
      if ( Param.subsection_x + Param.inputDataBufferZone - removedBufferX < 0 ||
           Param.subsection_x + Param.subsection_width + Param.inputDataBufferZone + removedBufferX > Param.X + 2 * Param.inputDataBufferZone ||
           Param.subsection_y + Param.inputDataBufferZone - removedBufferY < 0 ||
           Param.subsection_y + Param.subsection_height + Param.inputDataBufferZone + removedBufferY > Param.Y + 2 * Param.inputDataBufferZone )
        throw aol::Exception ( "The size of the buffer zone is too small for the given rotation angle! (Param: inputDataBufferZone, specimenRotation)", __FILE__, __LINE__ );
      
      // ... inputArg: exit wave of size Param.subsection_width x Param.subsection_height
      Image<RealType> croppedRealSpaceExitWave ( Param.subsection_width, Param.subsection_height );
      
      const int x = Param.subsection_x + Param.inputDataBufferZone - removedBufferX;
      const int y = Param.subsection_y + Param.inputDataBufferZone - removedBufferY;
      
      RealSpaceExitWave[0].copyBlockTo ( x, y, croppedRealSpaceExitWave[0] );
      RealSpaceExitWave[1].copyBlockTo ( x, y, croppedRealSpaceExitWave[1] );
      
      croppedRealSpaceExitWave.FourierTransformTo ( inputArg.ExitWave );
      
      // ... inputArg: translation matrix
      if ( verbose )
        cerr << "\tCreating the translation matrix..." << endl;
      
      switch ( Param.specimenDriftMode ) {
        case 0:
          for ( int i = 0; i < Param.N ; i++ ) {
            inputArg.Translation.set ( i, 0, i * Param.specimenDriftX );
            inputArg.Translation.set ( i, 1, i * Param.specimenDriftY );
          }
          break;
        case 1:
          for ( int i = 0; i < Param.N ; i++ ) {
            inputArg.Translation.set ( i, 0, Param.specimenPosX[i] );
            inputArg.Translation.set ( i, 1, Param.specimenPosY[i] );
          }
          break;
        default:
          throw aol::Exception ( "Invalid specimenDriftMode parameter!", __FILE__, __LINE__ );
      }
      
      // ... inputArg: focus values
      for ( int i = 0; i < Param.N ; i++ )
        inputArg.Focus[i] = Param.Focus[i];
      
      // ... input: focus image series
      SimulateInputImageSeries ( input.ImageSeries, inputArg.Focus, inputArg.Translation, RealSpaceExitWave, Param, verbose );
    } else {
      // Read the focus image series from file...
      LoadInputImageSeries ( input.ImageSeries, Param.inputDataFile, Param, verbose );
    }
    
    // input: adjust Param
    input.Param.lenX *= static_cast<RealType> ( Param.subsection_width ) / Param.X;
    input.Param.lenY *= static_cast<RealType> ( Param.subsection_height ) / Param.Y;
    
    input.Param.X = Param.subsection_width;
    input.Param.Y = Param.subsection_height;
    
    input.Param.subsection_x = 0;
    input.Param.subsection_y = 0;
    
    inputArg.ExitWave.LowpassFilter ( input.Param );
  }
  
  if ( verbose )
    cerr << endl;
}

/**
 * \brief Saves the input data and inputArg to a subdirectory at reconstructionDir
 */
template <typename RealType>
void SaveInputData ( const InputData<RealType>& input,
                     const Arguments<RealType>& inputArg,
                     const string& reconstructionDir,
                     const bool verbose = true ) {
  // input data
  if ( !aol::directoryExists ( reconstructionDir + "InputData" ) ) {
    if ( verbose )
      cerr << "Saving the input data..." << endl;
    
    MakeDirectories ( reconstructionDir + "InputData/" );
    
    // input: ImageSeries
    Image<RealType> RealSpaceImage ( input.Param.X, input.Param.Y );
    
    // The image series is always saved in the Quocmesh format. If, however,
    // the format requested with outputImageFormat is different from the
    // Quocmesh format, then a subdirectory with a copy of the input series
    // in the requested format is created as well.
    const string subdir_name = ( input.Param.outputImageFormat == 0 ? "" : "tiff/" );
    if ( input.Param.outputImageFormat != 0 )
      MakeDirectories ( reconstructionDir + "InputData/" + subdir_name );
    
    for ( int i = 0; i < input.Param.N ; i++ ) {
      input.ImageSeries[i].FourierTransformTo ( RealSpaceImage );
      
      char id[1024];
      sprintf ( id, input.Param.outputNumberingFormat.c_str ( ), i );
      
      RealSpaceImage[0].saveToFile ( ( reconstructionDir + "InputData/Image" + id + ".q2bz" ).c_str ( ) );
      
      if ( input.Param.outputImageFormat != 0 )
        SaveImage ( RealSpaceImage, true, reconstructionDir + "InputData/" + subdir_name + "Image" + id, input.Param.outputImageFormat );
      
      if ( verbose )
        cerr << "\r\tSaving the image series... " << setw ( 3 ) << i+1 << " / " << input.Param.N;
    }
    ofstream hiddenfile ( reconstructionDir + "InputData/.inpainting_done" );
    if ( hiddenfile.is_open ( ) )
      hiddenfile.close ( );
    
    if ( verbose )
      cerr << "\r\tSaving the image series... done.          " << endl;
    
    // input: Param
    input.Param.save ( reconstructionDir + "InputData/ParameterFile.param" );
  }
  
  // inputArg
  if ( !aol::directoryExists ( reconstructionDir + "InputData/CroppedArg" ) && !inputArg.ExitWave.isZero ( ) ) {
    MakeDirectories ( reconstructionDir + "InputData/CroppedArg/" );
    inputArg.saveToFile ( ( reconstructionDir + "InputData/CroppedArg/.Arg" ).c_str ( ) );
  }
  
  if ( verbose )
    cerr << endl;
}

/**
 * \brief Saves the scaled input arguments to a subdirectory of reconstructionDir
 */
template <typename RealType>
void SaveScaledInputArg ( const Arguments<RealType>& ScaledInputArg,
                          const Arguments<RealType>& ScaleMask,
                          const string& reconstructionDir ) {
  if ( !aol::fileExists ( ( reconstructionDir + "InputData/CroppedArg/.ScaledArg" ).c_str ( ) ) &&
       !ScaledInputArg.ExitWave.isZero ( ) ) {
    ScaledInputArg.saveToFile ( ( reconstructionDir + "InputData/CroppedArg/.ScaledArg" ).c_str ( ) );
    ScaleMask.saveToFile ( ( reconstructionDir + "InputData/CroppedArg/.ScaleMask" ).c_str ( ) );
    
    // Create shell scripts for quick image generation, derivative tests and step size tests
    WriteShellScript ( reconstructionDir + "InputData/CroppedArg/GenerateImages.sh", ".ScaledArg ../../" );
    //WriteShellScript ( reconstructionDir + "InputData/CroppedArg/DerivativeTest.sh", ".ScaledArg ../../ 20" );
    //WriteShellScript ( reconstructionDir + "InputData/CroppedArg/StepsizeTest.sh", ".ScaledArg ../../" );
  }
}

/**
 * \brief Creates the initial guess or reads it from file if a reconstruction is resumed
 */
template <typename RealType>
void CreateInitialGuess ( Arguments<RealType>& InitialGuess,
                          const InputData<RealType>& input,
                          const Arguments<RealType>& inputArg,
                          const Arguments<RealType>& ScaleMask,
                          const bool resumeReconstruction,
                          const string& reconstructionDir,
                          const bool verbose = true ) {
  // Check the size and format of the initial guess
  if ( InitialGuess.ExitWave.getNumX ( ) != input.Param.X ||
       InitialGuess.ExitWave.getNumY ( ) != input.Param.Y ||
       InitialGuess.Translation.getNumX ( ) != input.Param.N ||
       InitialGuess.Translation.getNumY ( ) != 2 ||
       InitialGuess.Focus.size ( ) != input.Param.N )
    throw aol::Exception ( "Invalid initial guess size!", __FILE__, __LINE__ );
  
  if ( !InitialGuess.ExitWave.correctFormat ( Space::FourierSpace, true ) )
    throw aol::Exception ( "Invalid exit wave format!", __FILE__, __LINE__ );
  
  if ( resumeReconstruction ) {
    // An earlier reconstruction is resumed: read the "initial guess" from file
    if ( verbose )
      cerr << "Loading the last estimate from the previous reconstruction as the initial guess..." << endl;
    
    InitialGuess.reallocateExitWave ( input.Param.X + 2 * input.Param.reconstructionBufferZone, input.Param.Y + 2 * input.Param.reconstructionBufferZone );
    
    ifstream fs ( reconstructionDir + "IntermediateResults/.lastresult" );
    if ( !fs.is_open ( ) ) {
      // The file likely could not be opened because it does not exist, which means that
      // less than Param.saveDataIterations minimization steps have been performed in the
      // last run; in this case the initial guess is loaded instead
      cerr << "Warning: unable to find the last estimate from the previous reconstruction. Loading the initial guess instead..." << endl;
      
      InitialGuess.loadFromFile ( reconstructionDir + "InitialGuess/.ScaledArg" );
    } else {
      string currentEstimateDir;
      getline ( fs, currentEstimateDir );
      InitialGuess.loadFromFile ( reconstructionDir + currentEstimateDir + ".ScaledArg" );
    }
  } else {
    // Generate the initial guess
    if ( verbose )
      cerr << "Generating the initial guess..." << endl;
    
    InitialGuess.setZero ( );
    
    // Initial guess for the exit wave
    switch ( input.Param.initialGuess_ExitWave ) {
      case 0:   // zero
        break;
      case 1: { // square root of the mean intensity of the focus series
        RealType mean = 0;
        for ( unsigned int i = 0; i < input.ImageSeries.size ( ) ; i++ )
          if ( input.ImageSeries[i].correctFormat ( Space::RealSpace, false ) )
            mean += input.ImageSeries[i][0].getMeanValue ( );
          else if ( input.ImageSeries[i].correctFormat ( Space::FourierSpace, true ) )
            mean += input.ImageSeries[i][0].get ( input.Param.X / 2, input.Param.Y / 2 );
          else
            throw aol::Exception ( "Invalid image format of the focus series!", __FILE__, __LINE__ );
        
        InitialGuess.ExitWave[0].set ( input.Param.X / 2, input.Param.Y / 2, sqrt ( mean / input.ImageSeries.size ( ) ) );
      } break;
      case 2:   // simulated exit wave (if available)
        if ( inputArg.ExitWave.isZero ( ) )
          throw aol::Exception ( "The simulated exit wave is not available for the initial guess! (Param.initialGuess_ExitWave)", __FILE__, __LINE__ );
        InitialGuess.ExitWave = inputArg.ExitWave;
        break;
      default:
        throw aol::Exception ( "Invalid initial guess for the exit wave! (Param.initialGuess_ExitWave)", __FILE__, __LINE__ );
        break;
    }
    
    switch ( input.Param.initialGuess_Translation ) {
      case 0:   // zero
        break;
      case 1: { // cross-correlation of successive images in the series
        const int maxAbsShift = static_cast<int> ( ceil ( input.Param.initialGuess_TranslationMaxShift * max ( input.Param.X / input.Param.lenX, input.Param.Y / input.Param.lenY ) ) );
        
        RealType dx = 0, dy = 0;
        for ( unsigned int i = 0; i < input.ImageSeries.size ( ) - 1 ; i++ ) {
          if ( verbose )
            cerr << "\r\tRegistering successive images with the cross-correlation... " << setw ( 3 ) << i+1 << " / " << input.ImageSeries.size ( ) - 1;
          
          const int numHints = input.Param.initialGuess_TranslationHint_Img.size ( );
          
          int j = 0;
          for ( ; j < numHints ; j++ )
            if ( i == static_cast<unsigned int> ( input.Param.initialGuess_TranslationHint_Img[j] ) )
              break;
          
          RealType x, y;
          if ( j < numHints ) {
            // Use the manually set estimate as a hint for the cross-correlation registration
            CrossCorrelationRegistration ( input.ImageSeries[i],
                                           input.ImageSeries[i+1],
                                           maxAbsShift,
                                           x,
                                           y,
                                           -input.Param.initialGuess_TranslationHint_ShiftX[j],
                                           -input.Param.initialGuess_TranslationHint_ShiftY[j] );
          } else {
            // Register the images i and i+1 with the cross-correlation
            CrossCorrelationRegistration ( input.ImageSeries[i], input.ImageSeries[i+1], maxAbsShift, x, y );
          }
          
          dx += x;
          dy += y;
          
          InitialGuess.Translation.set ( i + 1, 0, -dx * input.Param.lenX / input.Param.X );
          InitialGuess.Translation.set ( i + 1, 1, -dy * input.Param.lenY / input.Param.Y );
        }
        if ( verbose )
          cerr << " done." << endl;
        
      } break;
      case 2:   // correct translation values (if available)
        if ( inputArg.ExitWave.isZero ( ) )
          throw aol::Exception ( "The correct translation values are not available for the initial guess! (Param.initialGuess_Translation)", __FILE__, __LINE__ );
        
        InitialGuess.Translation = inputArg.Translation;
        InitialGuess.Translation *= -1;
        break;
      case 3: { // add Gaussian noise to the correct translations (if available)
        if ( inputArg.ExitWave.isZero ( ) )
          throw aol::Exception ( "The correct translation values are not available for the initial guess! (Param.initialGuess_Translation)", __FILE__, __LINE__ );
        
        aol::RandomGenerator rndGen;
        rndGen.randomize ( );
        
        for ( int j = 0; j < 2 ; j++ ) for ( int i = 1; i < input.Param.N ; i++ )
          InitialGuess.Translation.set ( i, j, rndGen.normalrReal<RealType> ( -inputArg.Translation.get ( i, j ), input.Param.initialGuess_TranslationStdDev ) );
      } break;
      default:
        throw aol::Exception ( "Invalid initial guess for the translations! (Param.initialGuess_Translation)", __FILE__, __LINE__ );
        break;
    }
    
    // Initial guess for the base focus value
    switch ( input.Param.initialGuess_Focus ) {
      case 0:   // zero (except for the first focus value)
        InitialGuess.Focus[0] = input.Param.Focus[0];
        break;
      case 1:   // the values given by the Focus vector
        for ( int i = 0; i < input.Param.N ; i++ )
          InitialGuess.Focus[i] = input.Param.Focus[i];
        break;
      case 2: { // the values given by the Focus vector distorted by Gaussian noise
        aol::RandomGenerator rndGen;
        rndGen.randomize ( );
        
        InitialGuess.Focus[0] = input.Param.Focus[0];
        for ( int i = 1; i < input.Param.N ; i++ )
          InitialGuess.Focus[i] = rndGen.normalrReal<RealType> ( input.Param.Focus[i], input.Param.initialGuess_FocusStdDev );
      } break;
      default:
        throw aol::Exception ( "Invalid initial guess for the focus values! (Param.initialGuess_Focus)", __FILE__, __LINE__ );
        break;
    }
    
    // Add the reconstruction buffer zone to the exit wave and scale the initial guess
    AddReconstructionBufferZone ( InitialGuess, input.Param );
    InitialGuess /= ScaleMask;
  }
}

/**
 * \brief Saves the initial guess to the reconstructionDir/InitialGuess directory
 */
template <typename RealType>
void SaveInitialGuess ( const Arguments<RealType>& InitialGuess,
                        const Arguments<RealType>& ScaleMask,
                        const string& reconstructionDir ) {
  if ( !aol::directoryExists ( reconstructionDir + "InitialGuess" ) ) {
    MakeDirectories ( reconstructionDir + "InitialGuess/" );
    
    InitialGuess.saveToFile ( ( reconstructionDir + "InitialGuess/.ScaledArg" ).c_str ( ) );
    ScaleMask.saveToFile ( ( reconstructionDir + "InitialGuess/.ScaleMask" ).c_str ( ) );
    
    // Create shell scripts in the same directory for quick image generation, derivative tests and stepsize tests
    WriteShellScript ( reconstructionDir + "InitialGuess/GenerateImages.sh", ".ScaledArg ../" );
    //WriteShellScript ( reconstructionDir + "InitialGuess/DerivativeTest.sh", ".ScaledArg ../ 20" );
    //WriteShellScript ( reconstructionDir + "InitialGuess/StepsizeTest.sh", ".ScaledArg ../" );
  }
}

/**
 * Adds the reconstruction buffer zone to the focus series using periodic continuation and
 * adjusts the parameters X, Y, lenX and lenY of input.Param accordingly.
 */
template <typename RealType>
void AddReconstructionBufferZone ( InputData<RealType>& input,
                                   const bool verbose = true ) {
  if ( input.Param.reconstructionBufferZone == 0 )
    return;
  
  // Add the reconstruction buffer zone to the image series
  if ( verbose )
    cerr << "Adding the reconstruction buffer zone to the focus series... " << setw ( 3 ) << 0 << " / " << input.ImageSeries.size ( );
  
  for ( unsigned int i = 0; i < input.ImageSeries.size ( ) ; i++ ) {
    // Transform the i-th image to real space
    Image<RealType> RealSpaceImage ( input.Param.X, input.Param.Y );
    input.ImageSeries[i].FourierTransformTo ( RealSpaceImage );
    
    // Add the buffer zone (recall that RealSpaceImage is real-valued)
    Image<RealType> RealSpaceImageBuf ( input.Param.X + 2 * input.Param.reconstructionBufferZone, input.Param.Y + 2 * input.Param.reconstructionBufferZone );
    for ( int y = 0; y < RealSpaceImageBuf.getNumY ( ) ; y++ ) for ( int x = 0; x < RealSpaceImageBuf.getNumX ( ) ; x++ ) {
      int Img_x = ( x - input.Param.reconstructionBufferZone ) % input.Param.X;
      int Img_y = ( y - input.Param.reconstructionBufferZone ) % input.Param.Y;
      
      if ( Img_x < 0 ) Img_x += input.Param.X;
      if ( Img_y < 0 ) Img_y += input.Param.Y;
      
      RealSpaceImageBuf[0].set ( x, y, RealSpaceImage[0].get ( Img_x, Img_y ) );
    }
    // Transform the image back to Fourier space
    input.ImageSeries[i].reallocate ( RealSpaceImageBuf.getNumX ( ), RealSpaceImageBuf.getNumY ( ) );
    RealSpaceImageBuf.FourierTransformTo ( input.ImageSeries[i] );
    
    if ( verbose )
      cerr << "\rAdding the reconstruction buffer zone to the focus series... " << setw ( 3 ) << i + 1 << " / " << input.ImageSeries.size ( );
  }
  
  if ( verbose )
    cerr << endl;
  
  // Adjust the parameters X, Y, lenX and lenY
  input.Param.lenX *= static_cast<RealType> ( input.Param.X + 2 * input.Param.reconstructionBufferZone ) / input.Param.X;
  input.Param.lenY *= static_cast<RealType> ( input.Param.Y + 2 * input.Param.reconstructionBufferZone ) / input.Param.Y;
  
  input.Param.X += 2 * input.Param.reconstructionBufferZone;
  input.Param.Y += 2 * input.Param.reconstructionBufferZone;
}

/**
 * Adds the reconstruction buffer zone to an exit wave given by Arg using constant continuation
 * with the mean value
 */
template <typename RealType>
void AddReconstructionBufferZone ( Arguments<RealType>& Arg,
                                   const Parameters<RealType>& Param ) {
  // Transform the exit wave to real space
  Image<RealType> RealSpaceExitWave ( Param.X, Param.Y );
  Arg.ExitWave.FourierTransformTo ( RealSpaceExitWave );
  
  // Add the buffer zone
  RealType mean[2] = { RealSpaceExitWave[0].getMeanValue ( ), RealSpaceExitWave[1].getMeanValue ( ) };
  
  Image<RealType> RealSpaceExitWaveBuf ( Param.X + 2 * Param.reconstructionBufferZone, Param.Y + 2 * Param.reconstructionBufferZone );
  
  for ( int j = 0; j < 2 ; j++ )
    for ( int y = 0; y < RealSpaceExitWaveBuf.getNumY ( ) ; y++ ) for ( int x = 0; x < RealSpaceExitWaveBuf.getNumX ( ) ; x++ )
      RealSpaceExitWaveBuf[j].set ( x, y, mean[j] );
  
  for ( int j = 0; j < 2 ; j++ )
    for ( int y = 0; y < Param.Y ; y++ ) for ( int x = 0; x < Param.X ; x++ )
      RealSpaceExitWaveBuf[j].set ( x + Param.reconstructionBufferZone, y + Param.reconstructionBufferZone, RealSpaceExitWave[j].get ( x, y ) );
  
  // Transform the exit wave back to Fourier space
  Arg.reallocateExitWave ( RealSpaceExitWaveBuf.getNumX ( ), RealSpaceExitWaveBuf.getNumY ( ) );
  RealSpaceExitWaveBuf.FourierTransformTo ( Arg.ExitWave );
}

#endif  // EWR_RECONSTRUCTIONBASE_H
