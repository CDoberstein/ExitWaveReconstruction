#ifndef EWR_OPENCLWRAPPER_H
#define EWR_OPENCLWRAPPER_H

#define CL_TARGET_OPENCL_VERSION 120

// Enable OpenCL exceptions (instead of C-style error handling)
#define __CL_ENABLE_EXCEPTIONS

#ifdef __APPLE__
#include "cl.hpp"
#else
#include <CL/cl.hpp>
#endif

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>

/**
 * \brief Initializes an OpenCL context and kernel
 */
class OpenCLWrapper {
public:
  // Platform (= driver)
  cl::Platform platform;
  
  // Device (= GPU)
  cl::Device device;
  
  // OpenCL context
  cl::Context context;
  
  // OpenCL C language program
  cl::Program program;
  
  // Device command queue
  cl::CommandQueue queue;
  
  // Kernel
  cl::Kernel kernel;

private:
  //! Convert value to string using std::stringstream
  template <typename T>
  std::string valToString ( const T& value ) {
    std::stringstream s;
    s << std::setprecision ( 20 ) << value;
    return s.str ( );
  }

  //! Print platform (driver) information
  void _printPlatformInfo ( const cl::Platform& p ) const {
    std::cerr << "       Name: " << p.getInfo<CL_PLATFORM_NAME> ( ) << std::endl
              << "       Vendor: " << p.getInfo<CL_PLATFORM_VENDOR> ( ) << std::endl
              << "       Profile: " << p.getInfo<CL_PLATFORM_PROFILE> ( ) << std::endl
              << "       Version: " << p.getInfo<CL_PLATFORM_VERSION> ( ) << std::endl
              << "       Extensions: " << p.getInfo<CL_PLATFORM_EXTENSIONS> ( ) << std::endl;
  }
  
  //! Print device (gpu) information
  void _printDeviceInfo ( const cl::Device& d ) const {
    std::cerr << "       Name: " << d.getInfo<CL_DEVICE_NAME> ( ) << std::endl
              << "       Vendor: " << d.getInfo<CL_DEVICE_VENDOR> ( ) << std::endl
              << "       Profile: " << d.getInfo<CL_DEVICE_PROFILE> ( ) << std::endl
              << "       OpenCL C version: " << d.getInfo<CL_DEVICE_OPENCL_C_VERSION> ( ) << std::endl
              << "       Version: " << d.getInfo<CL_DEVICE_VERSION> ( ) << std::endl
              << "       Built in kernels: " << d.getInfo<CL_DEVICE_BUILT_IN_KERNELS> ( ) << std::endl
              << "       Vendor ID: " << d.getInfo<CL_DEVICE_VENDOR_ID> ( ) << std::endl
              << "       Max. number of compute units: " << d.getInfo<CL_DEVICE_MAX_COMPUTE_UNITS> ( ) << std::endl
              << "       Max. work item dimensions: " << d.getInfo<CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS> ( ) << std::endl;
    
    std::vector<size_t> work_item_sizes = d.getInfo<CL_DEVICE_MAX_WORK_ITEM_SIZES> ( );
    std::cerr << "       Max. local work item sizes: " << work_item_sizes[0];
    for ( unsigned int j = 1; j < work_item_sizes.size ( ) ; j++ )
      std::cerr << ", " << work_item_sizes[j];
    std::cerr << std::endl;
    
    std::cerr << "       Max. work group size: " << d.getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE> ( ) << std::endl
              << "       Preferred vector width (char): " << d.getInfo<CL_DEVICE_PREFERRED_VECTOR_WIDTH_CHAR> ( ) << std::endl
              << "       Preferred vector width (short): " << d.getInfo<CL_DEVICE_PREFERRED_VECTOR_WIDTH_SHORT> ( ) << std::endl
              << "       Preferred vector width (int): " << d.getInfo<CL_DEVICE_PREFERRED_VECTOR_WIDTH_INT> ( ) << std::endl
              << "       Preferred vector width (half): " << d.getInfo<CL_DEVICE_PREFERRED_VECTOR_WIDTH_HALF> ( ) << std::endl
              << "       Preferred vector width (float): " << d.getInfo<CL_DEVICE_PREFERRED_VECTOR_WIDTH_FLOAT> ( ) << std::endl
              << "       Preferred vector width (double): " << d.getInfo<CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE> ( ) << std::endl
              << "       Max. clock frequency (Mhz): " << d.getInfo<CL_DEVICE_MAX_CLOCK_FREQUENCY> ( ) << std::endl
              << "       Address space size: " << d.getInfo<CL_DEVICE_ADDRESS_BITS> ( ) << std::endl
              << "       Max. size of memory allocation: " << d.getInfo<CL_DEVICE_MAX_MEM_ALLOC_SIZE> ( ) << std::endl
              << "       Image support: " << d.getInfo<CL_DEVICE_IMAGE_SUPPORT> ( ) << std::endl
              << "       Global memory size (bytes): " << d.getInfo<CL_DEVICE_GLOBAL_MEM_SIZE> ( ) << std::endl
              << "       Local memory size (bytes): " << d.getInfo<CL_DEVICE_LOCAL_MEM_SIZE> ( ) << std::endl
              << "       Max. size of constant buffer allocation (bytes): " << d.getInfo<CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE> ( ) << std::endl
              << "       Max. number of constant arguments: " << d.getInfo<CL_DEVICE_MAX_CONSTANT_ARGS> ( ) << std::endl
              << "       Built-in kernels: " << d.getInfo<CL_DEVICE_BUILT_IN_KERNELS> ( ) << std::endl
              << "       Double datatype support: " << d.getInfo<CL_DEVICE_DOUBLE_FP_CONFIG> ( ) << " (zero means no support)" << std::endl
              << "       Profiling timer resolution (ns): " << d.getInfo<CL_DEVICE_PROFILING_TIMER_RESOLUTION> ( ) << std::endl
              << "       Local memory type: " << ( d.getInfo<CL_DEVICE_LOCAL_MEM_TYPE> ( ) == CL_LOCAL ? "local" : "global" ) << std::endl
              << "       Extensions: " << d.getInfo<CL_DEVICE_EXTENSIONS> ( ) << std::endl;
  }
  
  //! Select a platform
  void selectPlatform ( bool select_default ) {
    // Get a vector of available platforms
    std::vector<cl::Platform> all_platforms;
    cl::Platform::get ( &all_platforms );
    
    if ( all_platforms.empty ( ) ) {
      std::cerr << "No platforms found. Check the OpenCL installation!" << std::endl;
      throw;
    }
    
    // Automatically choose platform if only one is available or select_default is true,
    // otherwise let the user select a platform
    if ( all_platforms.size ( ) == 1 || select_default ) {
      platform = all_platforms[0];
      return;
    }
    
    for ( unsigned int i = 0; i < all_platforms.size ( ) ; i++ ) {
      std::cerr << " - Platform " << i << ":" << std::endl;
      _printPlatformInfo ( all_platforms[i] );
    }
    
    unsigned int platform_ind;
    do {
      std::cerr << "Enter platform number (0-" << all_platforms.size ( ) - 1 << "): ";
      std::cin >> platform_ind;
    } while ( platform_ind > all_platforms.size ( ) - 1 );
    std::cerr << std::endl;
    
    platform = all_platforms[platform_ind];
  }
  
  //! Select device
  void selectDevice ( int gpu_number ) {
    // Get a vector of all GPUs on the current platform
    std::vector<cl::Device> all_gpu_devices;
    platform.getDevices ( CL_DEVICE_TYPE_GPU, &all_gpu_devices );
    
    if ( all_gpu_devices.empty ( ) ) {
      std::cerr << "No GPUs found. Check the OpenCL installation and the selected platform!" << std::endl;
      throw;
    }
    
    // Check if gpu_number is a valid GPU number
    if ( gpu_number < 0 || gpu_number >= static_cast<int> ( all_gpu_devices.size ( ) ) ) {
      std::cerr << "Invalid GPU number! (" << gpu_number << ")" << std::endl << std::endl;
      
      std::cerr << "Available GPUs:" << std::endl;
      for ( unsigned int i = 0; i < all_gpu_devices.size ( ) ; i++ ) {
        std::cerr << " - GPU " << i << ":" << std::endl;
        _printDeviceInfo ( all_gpu_devices[i] );
      }
      std::cerr << std::endl;
      
      do {
        std::cerr << "Enter GPU number (0-" << all_gpu_devices.size ( ) - 1 << "): ";
        std::cin >> gpu_number;
      } while ( gpu_number < 0 || gpu_number >= static_cast<int> ( all_gpu_devices.size ( ) ) );
      std::cerr << std::endl;
    }
    
    // Select GPU
    device = all_gpu_devices[gpu_number];
  }
  
  //! Initialize context
  void _initContext ( ) {
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wold-style-cast"
    cl_context_properties properties[3] = { CL_CONTEXT_PLATFORM, (cl_context_properties)(platform)(), 0 };
    #pragma GCC diagnostic pop
    
    cl_int err;
    context = cl::Context ( { device }, properties, NULL, NULL, &err );
    if ( err != CL_SUCCESS ) {
      std::cerr << "Creating context failed! (error code = " + valToString ( err ) + ")" << std::endl;
      throw;
    }
  }
  
  //! Build program
  void buildProgram ( const std::string& kernel_file, const std::vector<std::string>& constPrimitiveValues ) {
    // Read kernel sourcecode from file
    std::ifstream kFile ( kernel_file );
    if ( !kFile.good ( ) ) {
      std::cerr << "Kernel file not found!" << std::endl;
      throw;
    }
    
    std::stringstream buffer;
    buffer << kFile.rdbuf ( );
    std::string kernel_code = buffer.str ( );

    // Build kernel
    cl::Program::Sources sources;
    sources.push_back ( { kernel_code.c_str ( ), kernel_code.length ( ) } );
    
    program = cl::Program ( context, sources );
    
    std::string constPrimitiveValuesDefines = "";
    for ( unsigned int i = 0; i < constPrimitiveValues.size ( ) ; i++ )
      constPrimitiveValuesDefines += " -D" + constPrimitiveValues[i];
    try {
      program.build ( { device }, constPrimitiveValuesDefines.c_str ( ) );
    } catch ( ... ) {
      std::cerr << program.getBuildInfo<CL_PROGRAM_BUILD_LOG> ( device ) << std::endl;
      std::cerr << "Building failed!" << std::endl;
      throw;
    }
    
  }
  
  //! Create command queue
  void createCommandQueue ( bool enable_profiling ) {
    cl_int err;
    queue = cl::CommandQueue ( context, device, ( enable_profiling ? CL_QUEUE_PROFILING_ENABLE : 0 ), &err );
    if ( err != CL_SUCCESS ) {
      std::cerr << "Creating command queue failed! (error code = " + valToString ( err ) + ")" << std::endl;
      throw;
    }
  }
  
  //! Create kernel
  void createKernel ( const std::string& kernel_file ) {
    std::string kernel_name = kernel_file.substr ( kernel_file.find_last_of ( '/' ) + 1 );
    kernel_name = kernel_name.substr ( 0, kernel_name.find_last_of ( '.' ) );
    
    cl_int err;
    kernel = cl::Kernel ( program, kernel_name.c_str ( ), &err );
    if ( err != CL_SUCCESS ) {
      std::cerr << "Creating kernel failed! (error code = " + valToString ( err ) + ")" << std::endl;
      throw;
    }
  }

public:
  /**
   * \brief Initialize OpenCL environment
   * 
   * kernel_file: Path to the kernel program file that is executed on the
   *   GPU. The name of the kernel has to be the same as the plain file-
   *   name (i.e. kernel_file without extension and path).
   * 
   * constPrimitiveValues: vector of strings of the format "VAR=VALUE",
   *   where the variable VAR in the kernel code is replaced with VALUE
   *   at compile time.
   * 
   * select_default: if this is true, the first available platform is
   *   selected. Otherwise the user is asked to select a platform if
   *   multiple platforms are found.
   * 
   * gpu_number: integral value specifying the number of the GPU that is
   *   used (starting from 0).
   */
  void init ( const std::string& kernel_file,
              const std::vector<std::string>& constPrimitiveValues = std::vector<std::string> ( ),
              const bool select_default = true,
              const int gpu_number = 0 ) {
    initContext ( select_default, gpu_number );
    initKernel ( kernel_file, constPrimitiveValues );
  }

  /**
   * \brief Initialize OpenCL environment (context)
   * 
   * This function should be called first.
   * 
   * select_default: if this is true, the first available platform is
   *   selected. Otherwise the user is asked to select a platform if
   *   multiple platforms are found.
   * 
   * gpu_number: integral value specifying the number of the GPU that is
   *   used (starting from 0).
   */
  void initContext ( const bool select_default = true,
                     const int gpu_number = 0 ) {
    selectPlatform ( select_default );
    selectDevice ( gpu_number );
    _initContext ( );
  }
  
  /**
   * \brief Initialize OpenCL environment (kernel)
   * 
   * This function should be called second.
   * 
   * kernel_file: Path to the kernel program file that is executed on the
   *   GPU. The name of the kernel has to be the same as the plain file-
   *   name (i.e. kernel_file without extension and path).
   * 
   * constPrimitiveValues: vector of strings of the format "VAR=VALUE",
   *   where the variable VAR in the kernel code is replaced with VALUE
   *   at compile time.
   */
  void initKernel ( const std::string& kernel_file,
                    const std::vector<std::string>& constPrimitiveValues = std::vector<std::string> ( ),
                    bool enable_profiling = false ) {
    buildProgram ( kernel_file, constPrimitiveValues );
    createCommandQueue ( enable_profiling );
    createKernel ( kernel_file );
  }
                     

private:
  // Convert cl_channel_order bitfield to a human-readable string
  static std::string oclChannelLayoutToString ( const cl_channel_order& image_channel_order ) {
    std::string channel_layout;
    switch ( image_channel_order ) {
      case CL_R: channel_layout = "R"; break;
      case CL_Rx: channel_layout = "Rx"; break;
      case CL_A: channel_layout = "A"; break;
      case CL_INTENSITY: channel_layout = "INTENSITY"; break;
      case CL_LUMINANCE: channel_layout = "LUMINANCE"; break;
      case CL_DEPTH: channel_layout = "DEPTH"; break;
      case CL_RG: channel_layout = "RG"; break;
      case CL_RGx: channel_layout = "RGx"; break;
      case CL_RA: channel_layout = "RA"; break;
      case CL_RGB: channel_layout = "RGB"; break;
      case CL_RGBx: channel_layout = "RGBx"; break;
      case CL_RGBA: channel_layout = "RGBA"; break;
      case CL_ARGB: channel_layout = "ARGB"; break;
      case CL_BGRA: channel_layout = "BGRA"; break;
      default: channel_layout = "unkown"; break;
    }
    return channel_layout;
  }
  
  // Convert cl_channel_type to a human-readable string
  static std::string oclChannelTypeToString ( const cl_channel_type& image_channel_data_type ) {
    std::string channel_data_type;
    switch ( image_channel_data_type ) {
      case CL_SNORM_INT8: channel_data_type = "SNORM_INT8"; break;
      case CL_SNORM_INT16: channel_data_type = "SNORM_INT16"; break;
      case CL_UNORM_INT8: channel_data_type = "UNORM_INT8"; break;
      case CL_UNORM_INT16: channel_data_type = "UNORM_INT16"; break;
      case CL_UNORM_SHORT_565: channel_data_type = "UNORM_SHORT_565"; break;
      case CL_UNORM_SHORT_555: channel_data_type = "UNORM_SHORT_555"; break;
      case CL_UNORM_INT_101010: channel_data_type = "UNORM_INT_101010"; break;
      case CL_SIGNED_INT8: channel_data_type = "SIGNED_INT8"; break;
      case CL_SIGNED_INT16: channel_data_type = "SIGNED_INT16"; break;
      case CL_SIGNED_INT32: channel_data_type = "SIGNED_INT32"; break;
      case CL_UNSIGNED_INT8: channel_data_type = "UNSIGNED_INT8"; break;
      case CL_UNSIGNED_INT16: channel_data_type = "UNSIGNED_INT16"; break;
      case CL_UNSIGNED_INT32: channel_data_type = "UNSIGNED_INT32"; break;
      case CL_HALF_FLOAT: channel_data_type = "HALF_FLOAT"; break;
      case CL_FLOAT: channel_data_type = "FLOAT"; break;
      default: channel_data_type = "unknown"; break;
    }
    return channel_data_type;
  }

public:
  //! Print information on the used platform
  void printPlatformInfo ( ) const {
    std::cerr << "Platform:" << std::endl;
    _printPlatformInfo ( platform );
  }
  
  //! Print information on the used GPU
  void printDeviceInfo ( ) const {
    std::cerr << "GPU:" << std::endl;
    _printDeviceInfo ( device );
  }
  
  //! Print information on the supported image formats by the context
  void printSupportedImageFormats ( ) const {
    std::vector<cl::ImageFormat> image_formats;
    context.getSupportedImageFormats ( CL_MEM_READ_WRITE, CL_MEM_OBJECT_IMAGE2D, &image_formats );
    std::cerr << "Supported image formats (read-write, 2D): " << std::endl
              << "       format number: channel layout, channel data type" << std::endl;
    for ( unsigned int i = 0; i < image_formats.size ( ) ; i++ )
      std::cerr << "       " << i << ": " << oclChannelLayoutToString ( image_formats[i].image_channel_order ) << ", "
                                          << oclChannelTypeToString ( image_formats[i].image_channel_data_type ) << std::endl;
  }
  
  //! Print information on the kernel work groups
  void printWorkGroupInfo ( ) const {
    std::cerr << "Kernel work group:" << std::endl
              << "       Max. work group size: " << kernel.getWorkGroupInfo<CL_KERNEL_WORK_GROUP_SIZE> ( device ) << std::endl
              << "       Preferred work group size multiple: " << kernel.getWorkGroupInfo<CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE> ( device ) << std::endl
              << "       Kernel local memory usage (bytes): " << kernel.getWorkGroupInfo<CL_KERNEL_LOCAL_MEM_SIZE> ( device ) << std::endl
              << "       Kernel private memory usage (bytes): " << kernel.getWorkGroupInfo<CL_KERNEL_PRIVATE_MEM_SIZE> ( device ) << std::endl;
  }
};

/**
 * \brief Select a kernel based on the GPUs capabilities and RealType
 */
template <typename RealType>
void SelectKernel ( const OpenCLWrapper& ocl,
                    const std::string& kernel_name,
                    bool& useKernelV1,
                    string& kernel_file ) {
  // Check type of RealType
  if ( !std::is_same<RealType, float>::value && !std::is_same<RealType, double>::value ) {
    std::cerr << "Invalid RealType!" << std::endl;
    throw 0;
  }
  
  // Choose kernel
  if ( ocl.device.getInfo<CL_DEVICE_IMAGE_SUPPORT> ( ) == CL_TRUE && std::is_same<RealType, float>::value )
    useKernelV1 = true;
  else
    useKernelV1 = false;
  
  // Check GPU double datatype support if RealType == double
  if ( std::is_same<RealType, double>::value && ocl.device.getInfo<CL_DEVICE_DOUBLE_FP_CONFIG> ( ) == 0 ) {
    std::cerr << "GPU does not support the double datatype; RealType needs to be changed to float!" << std::endl;
    throw 0;
  }
  
  // Get path to the kernel file
  kernel_file = kernel_name + ( useKernelV1 ? "V1.cl" : "V2.cl" );
}

/**
 * \brief Get number of available GPUs on the default platform
 */
unsigned int getNumGPUs ( ) {
  std::vector<cl::Platform> all_platforms;
  std::vector<cl::Device> all_gpu_devices;
  
  try {
    cl::Platform::get ( &all_platforms );
    
    if ( all_platforms.empty ( ) ) {
      std::cerr << "No platforms found. Check the OpenCL driver!" << std::endl;
      throw;
    }
    
    all_platforms[0].getDevices ( CL_DEVICE_TYPE_GPU, &all_gpu_devices );
  } catch ( cl::Error *err ) {
    std::cerr << err->what ( ) << std::endl;
    throw err;
  }
  
  return all_gpu_devices.size ( );
}

#endif  // EWR_OPENCLWRAPPER_H
