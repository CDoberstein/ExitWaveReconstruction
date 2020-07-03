#ifndef EWR_GENERATEIMAGES_H
#define EWR_GENERATEIMAGES_H

#include "../ReconstructionBase.h"
#include "../IntegrationDomains.h"

// Print information on the command line usage
void PrintGenerateImageUsage ( ostream& out, const string& argv0 ) {
  const string bold_modifier = "\033[1m";
  const string reset_modifier = "\033[0m";

  out << bold_modifier << "Usage" << reset_modifier << endl
      << "\t" << argv0 << " <ScaledArgPath> <reconstructionDir>" << endl
      << endl
      << bold_modifier << "Arguments" << reset_modifier << endl
      << "\tScaledArgPath     Path to the scaled arguments in the .q2bz format." << endl
      << "\reconstructionDir  Path to the base output directory." << endl;
}

void GenerateImages ( const string& ScaledArgPath, const string& reconstructionDir ) {
  // Load the parameter file
  Parameters<RealType> Param ( reconstructionDir + "InputData/ParameterFile.param" );
  
  // Adjust the parameters X, Y, lenX and lenY
  Param.lenX *= static_cast<RealType> ( Param.X + 2 * Param.reconstructionBufferZone ) / Param.X;
  Param.lenY *= static_cast<RealType> ( Param.Y + 2 * Param.reconstructionBufferZone ) / Param.Y;
  
  Param.X += 2 * Param.reconstructionBufferZone;
  Param.Y += 2 * Param.reconstructionBufferZone;
  
  // Get the directory where the arguments are saved
  string ArgDir = "";
  auto pos = ScaledArgPath.find_last_of ( '/' );
  if ( pos != string::npos )
    ArgDir = ScaledArgPath.substr ( 0, pos + 1 );
  
  // Initialize the scale mask
  Arguments<RealType> ScaleMask ( Param.X, Param.Y, Param.N );
  ScaleMask.loadFromFile ( ArgDir + ".ScaleMask" );
  
  // Load the integration domains
  IntegrationDomains domains;
  
  if ( aol::fileExists ( ( ArgDir + ".integration_domains" ).c_str ( )  ) )
    domains.load ( ArgDir + ".integration_domains" );
  else
    domains.load ( reconstructionDir + "InitialGuess/.integration_domains" );
  
  // Prepare the OpenCL image simulation kernel
  KernelData<RealType> kData ( Param, true );
  
  // Load the arguments and undo the scaling
  Arguments<RealType> Arg ( Param.X, Param.Y, Param.N );
  Arg.loadFromFile ( ScaledArgPath );
  Arg *= ScaleMask;
  
  // Calculate and save the exit wave images
  cerr << "[CurrentEstimate/] Generating the exit wave images..." << endl;
  MakeDirectories ( ArgDir + "CurrentEstimate/NoReconstructionBuffer/" );
  
  Arg.saveArguments ( ArgDir + "CurrentEstimate/", Param.outputImageFormat );
  
  Arguments<RealType> ArgNoBuffer ( Param.X - 2 * Param.reconstructionBufferZone, Param.Y - 2 * Param.reconstructionBufferZone, Param.N );
  Arg.removeBufferZone ( ArgNoBuffer, Param.reconstructionBufferZone );
  ArgNoBuffer.saveArguments ( ArgDir + "CurrentEstimate/NoReconstructionBuffer/", Param.outputImageFormat );
  
  // Draw the integration domains on top of the real space exit wave images
  MakeDirectories ( ArgDir + "CurrentEstimate/IntegrationDomains/" );
  
  Image<RealType> RealSpaceExitWave ( Param.X, Param.Y );
  Arg.ExitWave.FourierTransformTo ( RealSpaceExitWave );
  
  qc::MultiArray<RealType, 2, 2> RealSpaceExitWaveAP ( Param.X, Param.Y );
  RealSpaceExitWave.CalculateAmplitudeAndPhase ( RealSpaceExitWaveAP );
  
  DrawIntegrationDomains ( RealSpaceExitWave[0], domains, Param, true );
  DrawIntegrationDomains ( RealSpaceExitWave[1], domains, Param, true );
  DrawIntegrationDomains ( RealSpaceExitWaveAP[0], domains, Param, true );
  DrawIntegrationDomains ( RealSpaceExitWaveAP[1], domains, Param, true );
  
  SaveImage ( RealSpaceExitWave, false, ArgDir + "CurrentEstimate/IntegrationDomains/EW_RealSpace", Param.outputImageFormat );
  SaveImage ( RealSpaceExitWaveAP, false, ArgDir + "CurrentEstimate/IntegrationDomains/EW_RealSpace", Param.outputImageFormat, true );
  
  // Calculate and save the focus image series
  cerr << "[FocusSeries/] Simulating the focus image series... " << setw ( 3 ) << 0 << " / " << Param.N;
  MakeDirectories ( ArgDir + "FocusSeries/NoReconstructionBuffer/" );
  
  vector<Image<RealType>> RealSpaceImages ( Param.N, Image<RealType> ( Param.X - 2 * Param.reconstructionBufferZone, Param.Y - 2 * Param.reconstructionBufferZone ) );
  int count = 0;
  const int nThreads = ( Param.simulationMode == 1 ? static_cast<int> ( kData.getNumGPUs ( ) ) : Param.numThreads );
  #pragma omp parallel for num_threads ( nThreads )
  for ( int i = 0; i < Param.N ; i++ ) {
    // Simulate the i-th image
    Image<RealType> SimulatedImage ( Param.X, Param.Y, Space::FourierSpace );
    SimulateImage ( SimulatedImage, Arg.ExitWave, Arg.Focus[i], Param, kData, omp_get_thread_num ( ) );
    
    // Modulate the image with the inverse translation for an easier comparison with the
    // input focus series
    Image<RealType> ModulatedImage ( Param.X, Param.Y, Space::FourierSpace );
    SimulatedImage.Modulate ( ModulatedImage, -Arg.Translation.get ( i, 0 ), -Arg.Translation.get ( i, 1 ), Param );
    
    // Transform the image to real space
    Image<RealType> RealSpaceImageBuf ( Param.X, Param.Y );
    ModulatedImage.FourierTransformTo ( RealSpaceImageBuf );
    
    RealSpaceImageBuf[0].copyBlockTo ( Param.reconstructionBufferZone, Param.reconstructionBufferZone, RealSpaceImages[i][0] );
    
    // Save the images
    char id[1024];
    sprintf ( id, Param.outputNumberingFormat.c_str ( ), i );
    
    SaveImage ( RealSpaceImageBuf, true, ArgDir + "FocusSeries/Image" + id, Param.outputImageFormat );
    SaveImage ( RealSpaceImages[i], true, ArgDir + "FocusSeries/NoReconstructionBuffer/Image" + id, Param.outputImageFormat );
    
    #pragma omp critical ( FocusSeriesSimulationOutput )
    {
      ++count;
      cerr << "\r[FocusSeries/] Simulating the focus image series... " << setw ( 3 ) << count << " / " << Param.N;
    }
  }
  
  cerr << endl;
  
  // Calculate and save the absolute value of the differences of the simulated
  // images and the input images on the respective integration domains
  cerr << "[FocusSeriesDiff/] Calculating the differences of simulated and input images..." << endl;
  MakeDirectories ( ArgDir + "FocusSeriesDiff/" );
  
  for ( int i = 0; i < Param.N ; i++ ) {
    char id[1024];
    sprintf ( id, Param.outputNumberingFormat.c_str ( ), i );
    
    qc::ScalarArray<RealType, qc::QC_2D> AbsDiffImage;
    AbsDiffImage.loadFromFile ( ( reconstructionDir + "InputData/Image" + id + ".q2bz" ).c_str ( ) );
    AbsDiffImage -= RealSpaceImages[i][0];
    AbsDiffImage.apply ( aol::Abs<RealType> );
    
    qc::ScalarArray<RealType, qc::QC_2D> AbsDiffImageID ( domains[i].width, domains[i].height );
    const RealType tx = Arg.Translation.get ( i, 0 ) * Param.X / Param.lenX;
    const RealType ty = Arg.Translation.get ( i, 1 ) * Param.Y / Param.lenY;
    AbsDiffImage.copyBlockTo ( static_cast<int> ( floor ( tx + domains[i].x ) ),
                               static_cast<int> ( floor ( ty + domains[i].y ) ),
                               AbsDiffImageID );
    
    qc::MultiArray<RealType, qc::QC_2D> tmp_AbsDiffImageID ( AbsDiffImageID.getNumX ( ), AbsDiffImageID.getNumY ( ) );
    tmp_AbsDiffImageID[0] = AbsDiffImageID;
    SaveImage ( tmp_AbsDiffImageID, true, ArgDir + "FocusSeriesDiff/Image" + id, Param.outputImageFormat );
  }
  
  // Calculate and save the difference of Arg and inputArg (if inputArg is available)
  if ( aol::fileExists ( reconstructionDir + "InputData/CroppedArg/.ScaledArg" ) ) {
    cerr << "[CurrentEstimateDiff/] Calculating the difference of the current estimate and inputArg..." << endl;
    MakeDirectories ( ArgDir + "CurrentEstimateDiff/" );
    
    Arguments<RealType> inputArg ( Param.X, Param.Y, Param.N );
    inputArg.loadFromFile ( reconstructionDir + "InputData/CroppedArg/.ScaledArg" );
    inputArg *= ScaleMask;
    
    Arg.saveArguments ( ArgDir + "CurrentEstimateDiff/", Param.outputImageFormat, inputArg );
  }
}

#endif  // EWR_GENERATEIMAGES_H
