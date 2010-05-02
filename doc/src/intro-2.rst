============
Installation
============

Linux Users
-----------

Documentation
~~~~~~~~~~~~~

Hermes2D Sphinx documentation (tutorial, benchmarks, examples) can be found in
the directory 'doc/'. Type 'make html' there to build it. The documentation is
available online at http://hpfem.org/hermes2d/doc/index.html.

Developer documentation can be compiled by running 'doxygen' in 'src/'.


Compilation
~~~~~~~~~~~

If you are using a Debian-based system, install the required libraries first:

:: 

    apt-get install cmake g++ gfortran libjudy-dev freeglut3-dev libsuitesparse-dev libglew1.5-dev python-dev python-numpy python-scipy cython python-matplotlib

(Note: cmake has to be at least version 2.6 or later, matplotlib has to be at
least 0.98.5.2 or higher.)

Configure:

::
    
    cmake .

Build:

::

    make

If you have more than one CPU, you can use "make -jN" where N is
the number of CPUs of your computer.

Tests
~~~~~

To execute all tests, do:

::

    make test

Note that some of the tests take a long time to finish. To just execute the
short running tests, do:

::

    make test-quick

More options
~~~~~~~~~~~~

You can turn on and off various components to build, just create the CMake.vars
file and add the following:

::

    set(WITH_EXAMPLES NO)
    set(WITH_PYTHON YES)

(and any other option that you would like to change, see CMakeLists.txt for the
whole list).

You can also easily generate it from a script (e.g. a debian/rules file) by:

::

    python -c 'print "set(H2D_COMPLEX no)\nset(WITH_EXAMPLES no)\nset(WITH_TUTORIAL no)\nset(WITH_PYTHON yes)\nset(WITH_GLUT no)\nset(WITH_UTIL no)"' > CMake.vars

If you are on OS X, you have to disable GLUT as the glut library is not easily
installable on OS X. To do so, just put the following line into your
CMake.vars:

::

    set(WITH_GLUT NO)


For development, it is good to say (in global CMake.vars):

::

    set(DEBUG YES) to compile debug versions
    set(RELEASE YES) to compile release versions

Then type:

::
 
    make debug    (to build debug versions)
    make release  (to build release versions)

Debugging with Eclipse
~~~~~~~~~~~~~~~~~~~~~~

To use eclipse as debugger, in the root folder of the project:
mkdir eclipse_build
cd eclipse_build
cmake -G"Eclipse CDT4 - Unix Makefiles" -D CMAKE_BUILD_TYPE=Debug ../

In Eclipse:
Import project using Menu File->Import
Select General->Existing projects into workspace:
Browse where your build tree is and select the root build tree directory. Keep "Copy projects into workspace" unchecked.


Install Hermes2D
~~~~~~~~~~~~~~~~

::

    cmake -DCMAKE_INSTALL_PREFIX=~/usr .
    make
    make install



Mac OS X Users
--------------

Documentation
~~~~~~~~~~~~~

Hermes2D Sphinx documentation (tutorial, benchmarks, examples) can be found in
the directory 'doc/'. Type 'make html' there to build it. The documentation is
available online at http://hpfem.org/hermes2d/doc/index.html.

Developer documentation can be compiled by running 'doxygen' in 'src/'.

Compilation
~~~~~~~~~~~

Download and install Xcode. You have to register as Apple Developer but this is quick. Do not pay any fees. Xcode contains make and other things that are required to build Hermes2D.

Also install cmake, cython, judy, and git. 

(Note: cmake has to be at least version 2.6 or later, matplotlib has to be at
least 0.98.5.2 or higher.)

Configure:

::

    cmake .

Build:

::

    make

If you have more than one CPU, you can use "make -j N" where N is
the number of CPUs of your computer.

Tests
~~~~~

To execute all tests, do:

::
 
    make test

Note that some of the tests take a long time to finish. To just execute the
short running tests, do:

::

    make test-quick


More options
~~~~~~~~~~~~

You can turn on and off various components to build, just create the CMake.vars
file and add the following:

::

    set(WITH_EXAMPLES NO)
    set(WITH_PYTHON YES)

(and any other option that you would like to change, see CMakeLists.txt for the
whole list).

You can also easily generate it from a script (e.g. a debian/rules file) by:

::

    python -c 'print "set(H2D_COMPLEX no)\nset(WITH_EXAMPLES no)\nset(WITH_TUTORIAL no)\nset(WITH_PYTHON yes)\nset(WITH_GLUT no)\nset(WITH_UTIL no)"' > CMake.vars

If you are on OS X, you have to disable GLUT as the glut library is not easily
installable on OS X. To do so, just put the following line into your
CMake.vars:

::

    set(WITH_GLUT NO)


For development, it is good to say (in global CMake.vars):

::

    set(DEBUG YES) to compile debug versions
    set(RELEASE YES) to compile release versions

Then type:

::

    make debug    (to build debug versions)
    make release  (to build release versions)

Install Hermes2D
~~~~~~~~~~~~~~~~

::

    cmake -DCMAKE_INSTALL_PREFIX=~/usr .
    make
    make install


Windows Cygwin Users
--------------------

Download and install the Linux emulator Cygwin from `here <http://www.cygwin.com/>`_ (the small icon in the top-right corner). While running setup.exe, you need to install 

cmake, gcc4, gfortran, git, gitk, libX11-devel, libXext-devel, libXt-devel, libXt, libXext, make, m4, openssl-devel, perl, 
python, wget, xextproto.

Then download, unpack, and build FEMhub as in Linux:

::

    git clone http://hpfem.org/git/hermes2d.git
    cd hermes2d
    cmake .
    make

For more details go to the Linux section above.

Windows MSVC Users
------------------

This section describes how to build and use Hermes2D in Microsoft Visual C++ 2008 (Express Edition). 
These instructions should probably work even for older versions of MS Visual C++ up to version 2003.

Known limitations and issues
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 - Stand-alone viewers do not compile.
 - Trilinos not supported.
 - Python not supported.

Building Hermes2D
~~~~~~~~~~~~~~~~~

 In order to build the library and examples, you need to:

 - Prepare dependecy libraries, see 'Dependency Check-list'
 - Copy a file 'my_hermes_root\MSVC2008\CMake.vars' to 'my_hermes_root'. The file contains settings for the projekct.
 - Modify the file 'my_hermes_root\CMake.vars'. Only the first line has to be modified, the rest of lines depends on project settings. Please, follow comments in the file.
 - In the directory 'my_hermes_root', run CMAKE using: 'cmake . -G "Visual Studio 9 2008". This will create project files.
 - Open a SLN file 'my_hermes_root\hermes2d.sln' and build Hermes2D. Actually, this step is not necessary if you plan to use Hermes2D in your projects rather than explore tutorials.

Configuration options
~~~~~~~~~~~~~~~~~~~~~

 Hermes2D is configured through preprocessor directives. Directives are generated by CMAKE and your settings might be overriden by CMAKE. The directives are:

  - NOGLUT : excludes GLUT-dependant parts. This replaces viewers with an empty implementation that does nothing if invoked. If used, a library 'freeglut.lib' does not need to be linked. 

  - ENABLE_VIEWER_GUI : enables GUI for viewers. Currently, only 'ScalarView' support GUI because this is an experimental feature. This directive is mutually exclusive with NOGLUT. If used, a library 'AntTweakBar.lib' does not need to be linked.

Using Hermes2D
~~~~~~~~~~~~~~
 
In order to used Hermes2D in your project, you need to do following steps. Steps has 5, 6, and 7 to be repeated for every configuration, i.e., Debug, Release. Except the step 7b, this can be done easily by setting the drop-down Configuration to 'All configurations' in the Project Property dialog.

  - Prepare Hermes2D to be buildable by MSVC, see 'Building Hermes2D'
  - Create your project in MSVC. Set the project to be empty Win32 console project.
  - Add either 'my_hermes2d_root\src\hermes2d-real.vcproj' or 'my_hermes2d_root\src\hermes2d-cplx.vcproj' project to your solution (<right click on solution>\Add\Existing Project...)
  - Set that your project depends on hermes2d-* project (<right click on your project>\Project Dependences...)
  - Add directories 'my_hermes2d_directory\src' and 'dependencies\include' to additional include directories (<right click on your project>\Properties\Configuration Properties\C/C++\Additional Include Directories)
  - Add directories 'dependencies\lib' to additional library directories (<right click on your project>\Properties\Configuration Properties\Linker\Additional Library Directories)
  - Deny warnings that are not indicating anything dangerous:
    - Avoid warnings about STL in DLL by denying a warning 4251 (<right click on your project>\Properties\Configuration Properties\C/C++\Advanced\Disable Specific Warnings, enter 4251)
    - Avoid warnings about standard functions that are not safe (<right click on your project>\Properties\Configuration Properties\C/C++\Preprocessor\Preprocessor Definitions, add _CRT_SECURE_NO_WARNINGS)
 
Dependency check-list
~~~~~~~~~~~~~~~~~~~~~

This list works for 32-bit version of Hermes2D. If you intend to cross-compile 64-bit version, you have to cross-compile all libraries. Asthe first step, create a  directory structure
	
  - in order to create the structure, execute 'prepare_dep_dir.bat'. Be sure to include a directory 'dependecies\bin' into 'PATH' environment variable.
  - all Hermes2D project files assumes that dependency libraries are available in a fixed directory structure. The root of this structure has to have the same parent as does Hermes2D director, i.e., if 'C:\my_work\hermes2d\' is a root of the Hermes2D directory, then 'C:\my_work\dependecies\' is a root of the dependency directory. Subdirectories are:    
    > dependencies\include: Header files (*.h) of dependency libraries.
    > dependencies\lib: Library files (*.lib) of dependency libraries.   
    > dependencies\bin: Binary modules (*.dll) of dependency libraries. Be sure to include a directory 'dependecies\bin' into 'PATH' environment variable.
	
  - JUDY
    - download judy (http://sourceforge.net/projects/judy/) and upack it 
    - open a command promt with MSVC variables set up: either use a command prompt in MSVC start menu or start a command prompt and execute VCVARS.BAT from the MSVC directory
    - switch to a directory containing JUDY sources, e.g., 'my_judy_root/src'
    - compile JUDY with 'build.bat': this creates Judy.dll and Judy.lib
    - copy 'Judy.dll', 'Judy.h', and 'Judy.lib' to 'bin', 'include', and 'lib' dependecy directories respectively

  - PTHREAD
    - download pthread binaries version 2.8.0 (ftp://sourceware.org/pub/pthreads-win32/)
    - copy 'lib\pthreadVCE2.dll', 'include\*.h' and 'lib\pthreadVCE2.lib' to 'bin', 'include', and 'lib' dependecy directories respectively.

  - UMFPACK
    - download UMFPACK source file package (http://www.cise.ufl.edu/research/sparse/umfpack/current/)
    - unpack source file into a directory that has the same parent as a directory where you unpacked UFconfig
    - copy the file 'my_hermes2d_root\UMFPACK.nmake' to the 'my_umfpack_root/Lib' directory
    - run MSVC command prompt similar as in the case of Judy and switch to 'my_umfpack_root\Lib'
    - compile UMFPACK using 'nmake -f UMFPACK.nmake'. Linking might take some time, please, be patient.
    - copy 'libumfpack.dll', all include files, and 'libumfpack.lib' to 'bin', 'include', and 'lib' dependecy directories, respectively.

    - UFConfig:

      - download UFconfig source file package (http://www.cise.ufl.edu/research/sparse/UFconfig/)
      - unpack it
      - copy UFconfig.h to 'include' dependecy directory
    - AMD:

      - download AMD source file package (http://www.cise.ufl.edu/research/sparse/amd/)
      - unpack source file into a directory that has the same parent as a directory where you unpacked UFconfig
      - copy the file 'my_hermes2d_root\MSVC2008\AMD.nmake' to a directory 'my_amd_directory\Lib'
      - run MSVC command prompt similar as in the case of Judy and switch to 'my_amd_directory\Lib'
      - compile AMD using 'nmake -f AMD.nmake'
      - copy 'amd.h', 'amd_internal.h', and 'libamd.lib' to 'include', and 'lib' dependecy directories respectively

  - CMAKE

    - download CMAKE (http://www.cmake.org/cmake/resources/software.html) version 2.6.X at minimum
    - install CMAKE such that it is accessible from every location

  - OpenGL support (optional)
    - if a directive NOGLUT is used, this step and all its substeps can be skipped

    - FREEGLUT 
      - download freeglut 2.6.0 (http://freeglut.sourceforge.net/) and unpack it
        - open a DSW or DSP file in MSVC, MSVC will convert file into a newer format, i.e., SLN or VCPROJ
        - compile either Debug or Release version. Debug version is recommended in a case of debugging.
        - copy 'freeglut.dll', 'freeglut.h', and 'freeglut.lib' to 'bin', 'include', and 'lib' dependency directories, respectively/
  
  - GLEW
    - download glew 1.5.2 (http://glew.sourceforge.net/) and unpack it
      - open a DSW file 'my_glew_root/builds/vc6' and let MSVC to convert it 
      - switch to 'Release' version
      - build a project 'glew_shared': this will create DLL file
      - copy 'my_glew_root/bin/glew32.dll', 'my_glew_root/include/GL/*.h', and 'my_glew_root/bin/glew32.lib' to 'bin', 'include/GL', and 'lib' dependency directories respectively
 	
  - AntTweakBar (optional)
    - if a directive ENABLE_VIEWER_GUI is *not* used, this step can be skipped
      - download a modified version 1.1.3 of AntTweakView (http://hpfem.org/hermes2d/)
      - unpack it
      - open SLN file in MSVC and compile it
      - copy 'AntTweakBar.dll', 'AntTweakBar.h', and 'AntTweakBar.lib' to 'bin', 'include', and 'lib' dependency directories respectively
	
  - ExodusII (optional)
    - if a directive WITH_EXODUSII is *not* used, this step including all sub-steps can be skipped.
	
    - Zlib
      - download sources of version 1.2.3 (http://www.zlib.net/) and unpack them
	- open 'my_zlib_root/projects/visualc6/zlib.dsw' (Visual C++ 6 Solution File) in MSVC and let MSVC to convert it
	- switch a configuration to 'Release DLL'
	- build project 'zlib': this will create DLL/LIB files in 'my_zlib_root/projects/visual6/Win32_DLL_Release'
	- copy 'zlib1.dll', 'zlib.h', and 'zlib1.lib' to 'bin', 'include', and 'lib' dependency directories respectively
 
    - HDF5
      - download sources of version 1.8 (ftp://ftp.hdfgroup.org/HDF5/current/src/) and unpack them
	- since SLIB is not used, comment out a line '#define H5_HAVE_FILTER_SZIP 1' in the header file 'my_hdf5_root/windows/src/H5pubconf.h'
	- copy the file 'my_hdf5_root/windows/src/H5pubconf.h' to the directory 'my_hdf5_root/src/'
	- run MSVC Command Prompt and switch to a directory 'my_hdf5_root/windows/proj'
	- set variable HDF5_EXT_ZLIB to 'my_dependencies\lib\zlib1.lib', e.g.,
	  set HDF5_EXT_ZLIB="C:\fem\dependencies\lib\zlib1.lib"
	- if SLIB is used, set variable HDF5_EXT_SLIB similarly, .e.g,
	  set HDF5_EXT_SLIB="C:\fem\dependencies\lib\slib.lib"
	- open SLN file in MSVC executing 'VCExpress.exe all/all.sln' in the command prompt and let MSVC to convert files
	- switch a configuration to 'Release'
	- build project 'hdf5_hldll': this will create DLL/LIB files in 'my_hdf5_root/proj/hdf5_hldll/Release/' and 'my_hdf5_root/proj/hdf5dll/Release/'
	- copy 'hdf5dll.dll' and 'hdf5dll.lib' to 'bin' and 'lib' dependency directories respectively
	- copy 'hdf5_hldll.dll' and 'hdf5_hldll.lib' to 'bin' and 'lib' dependency directories respectively
 
    - NetCDF
      - download sources of version 4.0.1 (http://www.unidata.ucar.edu/downloads/netcdf/netcdf-4_0_1/index.jsp) and unpack them
      - open a SLN file 'my_netcfd_root/win32/NET/netcdf.sln'
      - switch to 'Release' version
      - in properties of the project 'netcdf'
        > add paths 'my_hdf5_root/src/' and 'my_hdf5_root/hl/src' to 'C/C++ \ Additional Include Directories'
        > add a path 'dependencies/lib/' to 'Linker \ Additional Library Directories'
      - build project 'netcdf': this will create DLL/LIB files in 'my_netcdf_root/win32/NET/Release'
      - copy 'netcdf.dll' and 'netcdf.lib' to 'bin' and 'lib' dependency directories respectively
      - copy 'my_netcdf_root/libsrc4/netcdf.h' to 'include' dependency directory

    - ExodusII
      - download sources (http://sourceforge.net/projects/exodusii/) and unpack 'exodusii'
      - add a line 'set(NETCDF_INCLUDE_DIR "my_netcdf_root/libsrc4")' to the file 'my_exodusii_root/CMakeLists.txt' just after the line 'PROJECT(Exodusii)', .e.g.,
      - set(NETCDF_INCLUDE_DIR "C:/fem/dependencies/src/netcdf-4.0.1/libsrc4") # Be sure to use a slash '/' instead of a backslash '\'
      - generate MSVC project files using CMAKE, i.e., in command prompt run
        cmake . -G "Visual Studio 9 2008"
      - open a SLN file 'my_exodusii_root/ExodusII.sln' in MSVC
      - switch to 'Release' version
      - build a project 'exoIIv2c': this will create a LIB file in 'my_exodusii_root/cbind/Release'
      - copy 'exoIIv2c.lib' to 'lib' dependency directory structure
      - copy 'my_exodusii_root/cbind/include/exodusII.h' and 'my_netcdf_root/libsrc4/exodusII_ext.h' to 'include' dependency directory
	