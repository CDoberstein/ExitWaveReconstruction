#ifndef EWR_OPENCLKERNELDATA_H
#define EWR_OPENCLKERNELDATA_H

/// \cond
#include "OpenCLWrapper.h"
#include "Parameters.h"
#include "Image.h"
#include "Utilities.h"
/// \endcond

template <typename RealType>
class SimulateImage_KernelData {
public:
  // OpenCL environment
  OpenCLWrapper ocl;
  
  // Kernel arguments (buffer or image depending on the kernel)
  cl::Buffer bufInput;
  cl::Image2D imInput;
  
  cl::Buffer bufSimulatedImage;
  cl::Image2D imSimulatedImage;
  
  // Parameters used in the kernel
  const Parameters<RealType> Param;
  
  // Boolean that determines whether the V1 or the V2 kernel is used
  bool useKernelV1;
  
  // Temporary variables for the image simulation
  qc::MultiArray<RealType, 2, 2> PhaseTransfer;
  vector<RealType> interleaved_input;
  vector<RealType> interleaved_output;
  
  /**
   * \brief Constructor
   * 
   * Initializes the OpenCL environment, builds the kernel and associates
   * buffer or image objects with GPU memory.
   * 
   * If float precision is used for RealType and the GPU has OpenCL
   * image support the kernel SimulateImageV1 is chosen due to its
   * performance benefit over SimulateImageV2. Otherwise the kernel
   * SimulateImageV2 is chosen.
   * 
   * select_default: If this is true, the first OpenCL driver that is
   *   found is selected. Otherwise the user is asked to select a
   *   driver if multiple drivers are installed.
   * 
   * gpu_number: integral value specifying the number of the GPU that
   *   is used (starting from 0).
   */
   SimulateImage_KernelData ( const Parameters<RealType>& _Param,
                              const bool select_default = true,
                              const int gpu_number = 0 )
      : Param ( _Param ),
        PhaseTransfer ( _Param.X, _Param.Y ),
        interleaved_input ( 4 * _Param.X * _Param.Y ),
        interleaved_output ( 2 * _Param.X * _Param.Y ) {
    // Initialize the OpenCL context
    ocl.initContext ( select_default, gpu_number );
    
    // Choose an appropriate kernel
    string kernel_file;
    SelectKernel<RealType> ( ocl, "SimulateImage", useKernelV1, kernel_file );
    
    string path = __FILE__;
    path = path.substr ( 0, path.find_last_of ( '/' ) + 1 );
    
    kernel_file = path + "Kernel/" + kernel_file;
    
    // Build a vector with compiler definitions
    string datatype = ( is_same<RealType, float>::value ? "float" : "double" );
    string datatype_suffix = ( is_same<RealType, float>::value ? "f" : "" );
    
    vector<string> constPrimitiveValues = {
        "RealType=" + datatype,
        "RealType2=" + datatype + '2',
        "RealType3=" + datatype + '3',
        "RealType4=" + datatype + '4',
        "DIM_X=" + valToString ( Param.X ),
        "DIM_Y=" + valToString ( Param.Y ),
        "DIM_X_HALF=" + valToString ( Param.X / 2 ),
        "DIM_Y_HALF=" + valToString ( Param.Y / 2 ),
        "INV_LEN_X=" + valToString ( 1 / Param.lenX ) + datatype_suffix,
        "INV_LEN_Y=" + valToString ( 1 / Param.lenY ) + datatype_suffix,
        "WAVELENGTH=" + valToString ( Param.lambda ) + datatype_suffix,
        "WAVELENGTH_CUB=" + valToString ( aol::Cub<RealType> ( Param.lambda ) ) + datatype_suffix,
        "SPHERICAL_ABERRATION=" + valToString ( Param.SphericalAberration ),
        "SPATIAL_COHERENCE_COEFF=" + valToString ( -aol::Sqr<RealType> ( aol::NumberTrait<RealType>::pi * Param.alpha / Param.lambda ) ) + datatype_suffix,
        "TEMPORAL_COHERENCE_COEFF=" + valToString ( -0.5 * aol::Sqr<RealType> ( aol::NumberTrait<RealType>::pi * Param.FocalSpread * Param.lambda ) ) + datatype_suffix
      };
    
    // Initialize the OpenCL kernel
    ocl.initKernel ( kernel_file, constPrimitiveValues );
    
    // Assign the kernel arguments
    try {
      if ( useKernelV1 ) {
        imInput          = cl::Image2D ( ocl.context, CL_MEM_READ_ONLY, cl::ImageFormat ( CL_RGBA, CL_FLOAT ), Param.X, Param.Y );
        imSimulatedImage = cl::Image2D ( ocl.context, CL_MEM_WRITE_ONLY, cl::ImageFormat ( CL_RG, CL_FLOAT ), Param.X, Param.Y );
        
        ocl.kernel.setArg ( 0, imInput );
        ocl.kernel.setArg ( 1, imSimulatedImage );
      } else {
        bufInput          = cl::Buffer ( ocl.context, CL_MEM_READ_WRITE | CL_MEM_HOST_WRITE_ONLY, 4 * Param.X * Param.Y * sizeof ( RealType ) );
        bufSimulatedImage = cl::Buffer ( ocl.context, CL_MEM_WRITE_ONLY | CL_MEM_HOST_READ_ONLY, 2 * Param.X * Param.Y * sizeof ( RealType ) );
        
        ocl.kernel.setArg ( 0, bufInput );
        ocl.kernel.setArg ( 1, bufSimulatedImage );
      }
    } catch ( cl::Error *err ) {
      throw aol::Exception ( err->what ( ), __FILE__, __LINE__ );
    }
  }
};

template <typename RealType>
class SimulateImageFocusDerivative_KernelData {
public:
  // OpenCL environment
  OpenCLWrapper ocl;
  
  // Kernel arguments (buffer or image depending on the kernel)
  cl::Buffer bufInput;
  cl::Image2D imInput;
  
  cl::Buffer bufSimulatedImage;
  cl::Image2D imSimulatedImage;
  
  // Parameters used in the kernel
  const Parameters<RealType> Param;
  
  // Boolean that determines whether the V1 or the V2 kernel is used
  bool useKernelV1;
  
  // Temporary variables for the image simulation
  qc::MultiArray<RealType, 2, 2> PhaseTransfer;
  vector<RealType> interleaved_input;
  vector<RealType> interleaved_output;
  
  /**
   * \brief Constructor
   * 
   * Initializes OpenCL environment, builds the kernel and associates
   * buffer or image objects with GPU memory.
   * 
   * If float precision is used for RealType and the GPU has OpenCL
   * image support the kernel SimulateImageFocusDerivativeV1 is chosen
   * due to its performance benefit over SimulateImageFocusDerivativeV2.
   * Otherwise the kernel SimulateImageFocusDerivativeV2 is chosen.
   * 
   * select_default: If this is true, the first OpenCL driver that is
   *   found is selected. Otherwise the user is asked to select a
   *   driver if multiple drivers are installed.
   * 
   * gpu_number: integral value specifying the number of the GPU that
   *   is used (starting from 0).
   */
  SimulateImageFocusDerivative_KernelData ( const Parameters<RealType>& _Param,
                                      const bool select_default = true,
                                      const int gpu_number = 0 )
      : Param ( _Param ),
        PhaseTransfer ( _Param.X, _Param.Y ),
        interleaved_input ( 4 * _Param.X * _Param.Y ),
        interleaved_output ( 2 * _Param.X * _Param.Y ) {
    // Initialize the OpenCL context
    ocl.initContext ( select_default, gpu_number );
    
    // Choose an appropriate kernel
    string kernel_file;
    SelectKernel<RealType> ( ocl, "SimulateImageFocusDerivative", useKernelV1, kernel_file );
    
    string path = __FILE__;
    path = path.substr ( 0, path.find_last_of ( '/' ) + 1 );
    
    kernel_file = path + "Kernel/" + kernel_file;
    
    // Build a vector with compiler definitions
    string datatype = ( is_same<RealType, float>::value ? "float" : "double" );
    string datatype_suffix = ( is_same<RealType, float>::value ? "f" : "" );
    
    vector<string> constPrimitiveValues = {
        "RealType=" + datatype,
        "RealType2=" + datatype + '2',
        "RealType3=" + datatype + '3',
        "RealType4=" + datatype + '4',
        "DIM_X=" + valToString ( Param.X ),
        "DIM_Y=" + valToString ( Param.Y ),
        "DIM_X_HALF=" + valToString ( Param.X / 2 ),
        "DIM_Y_HALF=" + valToString ( Param.Y / 2 ),
        "INV_LEN_X=" + valToString ( 1 / Param.lenX ) + datatype_suffix,
        "INV_LEN_Y=" + valToString ( 1 / Param.lenY ) + datatype_suffix,
        "WAVELENGTH=" + valToString ( Param.lambda ) + datatype_suffix,
        "WAVELENGTH_SQR=" + valToString ( aol::Sqr<RealType> ( Param.lambda ) ) + datatype_suffix,
        "WAVELENGTH_CUB=" + valToString ( aol::Cub<RealType> ( Param.lambda ) ) + datatype_suffix,
        "WAVELENGTH_QRT=" + valToString ( aol::Sqr<RealType> ( aol::Sqr<RealType> ( Param.lambda ) ) ) + datatype_suffix,
        "SPHERICAL_ABERRATION=" + valToString ( Param.SphericalAberration ),
        "SPATIAL_COHERENCE_COEFF=" + valToString ( -aol::Sqr<RealType> ( aol::NumberTrait<RealType>::pi * Param.alpha / Param.lambda ) ) + datatype_suffix,
        "TEMPORAL_COHERENCE_COEFF=" + valToString ( -0.5 * aol::Sqr<RealType> ( aol::NumberTrait<RealType>::pi * Param.FocalSpread * Param.lambda ) ) + datatype_suffix,
        "PI=3.14159265358979324" + datatype_suffix
      };
    
    // Initialize the OpenCL kernel
    ocl.initKernel ( kernel_file, constPrimitiveValues );
    
    // Assign the kernel arguments
    try {
      if ( useKernelV1 ) {
        imInput          = cl::Image2D ( ocl.context, CL_MEM_READ_ONLY, cl::ImageFormat ( CL_RGBA, CL_FLOAT ), Param.X, Param.Y );
        imSimulatedImage = cl::Image2D ( ocl.context, CL_MEM_WRITE_ONLY, cl::ImageFormat ( CL_RG, CL_FLOAT ), Param.X, Param.Y );
        
        ocl.kernel.setArg ( 0, imInput );
        ocl.kernel.setArg ( 1, imSimulatedImage );
      } else {
        bufInput          = cl::Buffer ( ocl.context, CL_MEM_READ_WRITE | CL_MEM_HOST_WRITE_ONLY, 4 * Param.X * Param.Y * sizeof ( RealType ) );
        bufSimulatedImage = cl::Buffer ( ocl.context, CL_MEM_WRITE_ONLY | CL_MEM_HOST_READ_ONLY, 2 * Param.X * Param.Y * sizeof ( RealType ) );
        
        ocl.kernel.setArg ( 0, bufInput );
        ocl.kernel.setArg ( 1, bufSimulatedImage );
      }
    } catch ( cl::Error *err ) {
      throw aol::Exception ( err->what ( ), __FILE__, __LINE__ );
    }
  }
};

template <typename RealType>
class ExitWaveDerivative_KernelData {
public:
  // OpenCL environment
  OpenCLWrapper ocl;
  
  // Kernel arguments (buffer or image depending on the kernel)
  cl::Buffer bufInput;
  cl::Image2D imInput;
  
  cl::Buffer bufDiffImage;
  cl::Image2D imDiffImage;
  
  cl::Buffer bufDerivativePart;
  cl::Image2D imDerivativePart;
  
  // Parameters used in the kernel
  const Parameters<RealType> Param;
  
  // Boolean that determines whether the V1 or the V2 kernel is used
  bool useKernelV1;
  
  // Temporary variables for the derivative calculation
  qc::MultiArray<RealType, 2, 2> PhaseTransfer;
  vector<RealType> interleaved_input;
  vector<RealType> interleaved_DiffImage;
  vector<RealType> interleaved_DerivativePart;
  
  /**
   * \brief Constructor
   * 
   * Initializes OpenCL environment, builds the kernel and associates
   * buffer or image objects with GPU memory.
   * 
   * If float precision is used for RealType and the GPU has OpenCL
   * image support the kernel ExitWaveDerivativeV1 is chosen due to
   * its performance benefit over ExitWaveDerivativeV2. Otherwise the
   * kernel ExitWaveDerivativeV2 is chosen.
   * 
   * select_default: If this is true, the first OpenCL driver that is
   *   found is selected. Otherwise the user is asked to select a
   *   driver if multiple drivers are installed.
   * 
   * gpu_number: integral value specifying the number of the GPU that
   *   is used (starting from 0).
   */
  ExitWaveDerivative_KernelData ( const Parameters<RealType>& _Param,
                                  const bool select_default = true,
                                  const int gpu_number = 0 )
      : Param ( _Param ),
        PhaseTransfer ( _Param.X, _Param.Y ),
        interleaved_input ( 4 * _Param.X * _Param.Y ),
        interleaved_DiffImage ( 2 * _Param.X * _Param.Y ),
        interleaved_DerivativePart ( 2 * _Param.X * _Param.Y ) {
    // Initialize the OpenCL context
    ocl.initContext ( select_default, gpu_number );
    
    // Choose an appropriate kernel
    string kernel_file;
    SelectKernel<RealType> ( ocl, "ExitWaveDerivative", useKernelV1, kernel_file );
    
    string path = __FILE__;
    path = path.substr ( 0, path.find_last_of ( '/' ) + 1 );
    
    kernel_file = path + "Kernel/" + kernel_file;
    
    // Build a vector with compiler definitions
    string datatype = ( is_same<RealType, float>::value ? "float" : "double" );
    string datatype_suffix = ( is_same<RealType, float>::value ? "f" : "" );
    
    vector<string> constPrimitiveValues = {
        "RealType=" + datatype,
        "RealType2=" + datatype + '2',
        "RealType3=" + datatype + '3',
        "RealType4=" + datatype + '4',
        "DIM_X=" + valToString ( Param.X ),
        "DIM_Y=" + valToString ( Param.Y ),
        "DIM_X_HALF=" + valToString ( Param.X / 2 ),
        "DIM_Y_HALF=" + valToString ( Param.Y / 2 ),
        "INV_LEN_X=" + valToString ( 1 / Param.lenX ) + datatype_suffix,
        "INV_LEN_Y=" + valToString ( 1 / Param.lenY ) + datatype_suffix,
        "WAVELENGTH=" + valToString ( Param.lambda ) + datatype_suffix,
        "WAVELENGTH_CUB=" + valToString ( aol::Cub<RealType> ( Param.lambda ) ) + datatype_suffix,
        "SPHERICAL_ABERRATION=" + valToString ( Param.SphericalAberration ),
        "SPATIAL_COHERENCE_COEFF=" + valToString ( -aol::Sqr<RealType> ( aol::NumberTrait<RealType>::pi * Param.alpha / Param.lambda ) ) + datatype_suffix,
        "TEMPORAL_COHERENCE_COEFF=" + valToString ( -0.5 * aol::Sqr<RealType> ( aol::NumberTrait<RealType>::pi * Param.FocalSpread * Param.lambda ) ) + datatype_suffix
      };
    
    // Initialize the OpenCL kernel
    ocl.initKernel ( kernel_file, constPrimitiveValues );
    
    // Assign the kernel arguments
    try {
      if ( useKernelV1 ) {
        imInput          = cl::Image2D ( ocl.context, CL_MEM_READ_ONLY, cl::ImageFormat ( CL_RGBA, CL_FLOAT ), Param.X, Param.Y );
        imDiffImage      = cl::Image2D ( ocl.context, CL_MEM_READ_ONLY, cl::ImageFormat ( CL_RG, CL_FLOAT ), Param.X, Param.Y );
        imDerivativePart = cl::Image2D ( ocl.context, CL_MEM_WRITE_ONLY, cl::ImageFormat ( CL_RG, CL_FLOAT ), Param.X, Param.Y );
        
        ocl.kernel.setArg ( 0, imInput );
        ocl.kernel.setArg ( 1, imDiffImage );
        ocl.kernel.setArg ( 2, imDerivativePart );
      } else {
        bufInput          = cl::Buffer ( ocl.context, CL_MEM_READ_WRITE | CL_MEM_HOST_WRITE_ONLY, 4 * Param.X * Param.Y * sizeof ( RealType ) );
        bufDiffImage      = cl::Buffer ( ocl.context, CL_MEM_READ_WRITE | CL_MEM_HOST_WRITE_ONLY, 2 * Param.X * Param.Y * sizeof ( RealType ) );
        bufDerivativePart = cl::Buffer ( ocl.context, CL_MEM_WRITE_ONLY | CL_MEM_HOST_READ_ONLY, 2 * Param.X * Param.Y * sizeof ( RealType ) );
        
        ocl.kernel.setArg ( 0, bufInput );
        ocl.kernel.setArg ( 1, bufDiffImage );
        ocl.kernel.setArg ( 2, bufDerivativePart );
      }
    } catch ( cl::Error *err ) {
      throw aol::Exception ( err->what ( ), __FILE__, __LINE__ );
    }
  }
};

template <typename RealType> class KernelData;

template <typename RealType>
class KernelData {
public:
  vector<SimulateImage_KernelData<RealType>*> SimulateImageData;
  vector<SimulateImageFocusDerivative_KernelData<RealType>*> SimulateImageFocusDerivativeData;
  vector<ExitWaveDerivative_KernelData<RealType>*> ExitWaveDerivativeData;
  
  KernelData ( const Parameters<RealType>& Param, bool imageSimulationOnly = false ) {
    // Don't create kernel data if the image simulation is performed on the CPU
    string inputDataFileExt = Param.inputDataFile;
    if ( inputDataFileExt.length ( ) > 4 )
      inputDataFileExt = inputDataFileExt.substr ( inputDataFileExt.length ( ) - 4 );
    
    if ( Param.simulationMode == 0 && ( Param.inputSimulationMode == 0 ||
         ( Param.inputDataSource == 1 && inputDataFileExt != ".wav" ) ) )
      return;
    
    // Determine the number of GPUs that will be used
    int nGPUs = aol::Clamp ( Param.numThreads, 1, Param.N );
    
    int maxNumGPUs = static_cast<int> ( ::getNumGPUs ( ) );
    if ( maxNumGPUs == 0 )
      throw aol::Exception ( "No GPUs found on the default platform!", __FILE__, __LINE__ );
    
    nGPUs = aol::Clamp ( nGPUs, 1, maxNumGPUs );
    
    // Prepare the OpenCL environment and kernels
    SimulateImageData.resize ( nGPUs );
    if ( !imageSimulationOnly ) {
      SimulateImageFocusDerivativeData.resize ( nGPUs );
      ExitWaveDerivativeData.resize ( nGPUs );
    }
    
    for ( int i = 0; i < nGPUs ; i++ ) {
      SimulateImageData[i] = new SimulateImage_KernelData<RealType> ( Param, true, i );
      if ( !imageSimulationOnly ) {
        SimulateImageFocusDerivativeData[i] = new SimulateImageFocusDerivative_KernelData<RealType> ( Param, true, i );
        ExitWaveDerivativeData[i] = new ExitWaveDerivative_KernelData<RealType> ( Param, true, i );
      }
    }
  }
  
  ~KernelData ( ) {
    for ( unsigned int i = 0; i < SimulateImageData.size ( ) ; i++ )
      delete SimulateImageData[i];
    for ( unsigned int i = 0; i < SimulateImageFocusDerivativeData.size ( ) ; i++ )
      delete SimulateImageFocusDerivativeData[i];
    for ( unsigned int i = 0; i < ExitWaveDerivativeData.size ( ) ; i++ )
      delete ExitWaveDerivativeData[i];
  }
  
  unsigned int getNumGPUs ( ) const {
    return SimulateImageData.size ( );
  }
};

#endif  // EWR_OPENCLKERNELDATA_H
