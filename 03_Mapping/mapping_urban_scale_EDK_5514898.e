slurmstepd-cobalt1028: error: task/cgroup: task[0] unable to set taskset '0x0'
+ SCRIPT_PID=23559
+ set +x
+ /bin/bash -x /tmp/jobstart.23526
+ [[ hxB =~ i ]]
+ export -f module
+ export -f switchml
+ ENV=/usr/share/modules-tcl/init/profile.sh
+ export ENV
+ BASH_ENV=/usr/share/modules-tcl/init/bash
+ export BASH_ENV
+ [[ ! :/usr/share/modules-tcl/init/ksh-functions: =~ :/usr/share/modules-tcl/init/ksh-functions: ]]
+ '[' 4 -ge 3 ']'
+ [[ hxB =~ i ]]
+ MODULESHOME=/usr/share/modules-tcl
+ export MODULESHOME
+ [[ ! :/ccc/products/geos-3.4.2/intel--17.0.4.196/default/bin:/ccc/products/cdo-1.7.2rc6/intel--17.0.4.196/default/bin:/ccc/products/openmpi-2.0.4/intel--17.0.6.256/default/bin:/ccc/products/icc-17.0.6.256/system/default/17.0.6.256/bin/intel64:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:/opt/ibutils/bin:/ccc/products/ccc_users_env/bin:/ccc/products/ifort-17.0.6.256/system/default/17.0.6.256/bin/intel64:/ccc/products/nco-4.6.0/intel--17.0.4.196/default/bin:/ccc/products/proj-4.9.1/system/default/bin:/ccc/products/gdal-2.1.0/intel--17.0.4.196/default/bin:/ccc/products/r-3.4.4/intel--17.0.4.196__openmpi--2.0.2/default/bin:/ccc/cont004/home/ineris/gressena/bin: =~ :/usr/bin: ]]
++ manpath
+ manpath=/ccc/products/openmpi-2.0.4/intel--17.0.6.256/default/share/man:/ccc/products/icc-17.0.6.256/system/default/17.0.6.256/man/common:/ccc/products/share/man:/usr/share/man/en:/usr/share/man:/ccc/products/ccc_users_env/man/en:/ccc/products/ifort-17.0.6.256/system/default/17.0.6.256/man/common:/ccc/products/nco-4.6.0/intel--17.0.4.196/default/share/man:/ccc/products/r-3.4.4/intel--17.0.4.196__openmpi--2.0.2/default/share/man
+ [[ ! :/ccc/products/openmpi-2.0.4/intel--17.0.6.256/default/share/man:/ccc/products/icc-17.0.6.256/system/default/17.0.6.256/man/common:/ccc/products/share/man:/usr/share/man/en:/usr/share/man:/ccc/products/ccc_users_env/man/en:/ccc/products/ifort-17.0.6.256/system/default/17.0.6.256/man/common:/ccc/products/nco-4.6.0/intel--17.0.4.196/default/share/man:/ccc/products/r-3.4.4/intel--17.0.4.196__openmpi--2.0.2/default/share/man: =~ :/usr/share/man: ]]
+ '[' /usr/share/modulefiles/applications:/usr/share/modulefiles/environment:/usr/share/modulefiles/configurations:/usr/share/modulefiles/tools:/usr/share/modulefiles/graphics:/usr/share/modulefiles/parallel:/usr/share/modulefiles/libraries:/usr/share/modulefiles/compilers = '' ']'
+ '[' ccc/1.0:datadir/ineris:datadir/own:dfldatadir/own:flavor/buildtarget/x86_64:feature/mkl/single_node:flavor/buildmpi/openmpi/2.0:flavor/openmpi/standard:feature/openmpi/net/auto:.tuning/openmpi/2.0.4:flavor/buildcompiler/intel/17:feature/mkl/lp64:feature/mkl/sequential:mkl/17.0.6.256:licsrv/intel:c/intel/17.0.6.256:c++/intel/17.0.6.256:fortran/intel/17.0.6.256:intel/17.0.6.256:mpi/openmpi/2.0.4:nco/4.6.0:cdo/1.7.2rc6:geos/3.4.2:proj/4.9.1:gdal/2.1.0:r/3.4.4 = '' ']'
+ '[' -r /usr/share/modules-tcl/init/modulerc -a /usr/share/modulefiles/applications:/usr/share/modulefiles/environment:/usr/share/modulefiles/configurations:/usr/share/modulefiles/tools:/usr/share/modulefiles/graphics:/usr/share/modulefiles/parallel:/usr/share/modulefiles/libraries:/usr/share/modulefiles/compilers = '' -a ccc/1.0:datadir/ineris:datadir/own:dfldatadir/own:flavor/buildtarget/x86_64:feature/mkl/single_node:flavor/buildmpi/openmpi/2.0:flavor/openmpi/standard:feature/openmpi/net/auto:.tuning/openmpi/2.0.4:flavor/buildcompiler/intel/17:feature/mkl/lp64:feature/mkl/sequential:mkl/17.0.6.256:licsrv/intel:c/intel/17.0.6.256:c++/intel/17.0.6.256:fortran/intel/17.0.6.256:intel/17.0.6.256:mpi/openmpi/2.0.4:nco/4.6.0:cdo/1.7.2rc6:geos/3.4.2:proj/4.9.1:gdal/2.1.0:r/3.4.4 = '' ']'
+ Rlog=/ccc/work/cont004/ineris/gressena/microcapteurs_carto/MAPPING_SENSOR_DIR/03_Mapping/mapping_urban_scale_EDK.Rlog
+ module load gnu
+++ [[ hxB =~ i ]]
+++ echo 1
++ export MODULEINTERACTIVE=1
++ MODULEINTERACTIVE=1
++ /usr/bin/tclsh /usr/share/modules-tcl/libexec/modulecmd.tcl bash load gnu
load module c/gnu/7.3.0 (GNU C compiler)
load module c++/gnu/7.3.0 (GNU C++ compiler)
load module fortran/gnu/7.3.0 (GNU Fortran compiler)
load module gnu/7.3.0 (GNU Compiler Collection)
unload module r/3.4.4 (R)
unload module gdal/2.1.0 (GDAL)
unload module proj/4.9.1 (PROJ)
unload module geos/3.4.2 (GEOS)
unload module mpi/openmpi/2.0.4 (Open MPI)
unload module feature/openmpi/net/auto (MPI Network backend feature)
unload module flavor/openmpi/standard (Open MPI flavor)
unload module flavor/buildmpi/openmpi/2.0 (MPI build flavor)
load module flavor/buildmpi/openmpi/2.0 (MPI build flavor)
load module feature/openmpi/mpi_compiler/intel (MPI Compiler feature)
load module flavor/openmpi/standard (Open MPI flavor)
load module feature/openmpi/net/auto (MPI Network backend feature)
load module mpi/openmpi/2.0.4 (Open MPI)
load module geos/3.4.2 (GEOS)
load module proj/4.9.1 (PROJ)
load module gdal/2.1.0 (GDAL)
load module r/3.4.4 (R)
++ unset MODULEINTERACTIVE
+ eval 'GDAL_EXEDIR=/ccc/products/gdal-2.1.0/intel--17.0.4.196/default/bin;' export 'GDAL_EXEDIR;' 'PORTALS4_CQ_MAXPIDS=2;' export 'PORTALS4_CQ_MAXPIDS;' 'MPI_ROOT=/ccc/products/openmpi-2.0.4/intel--17.0.6.256/default;' export 'MPI_ROOT;' 'FEATURE_OPENMPI_MPI_COMPILER_TOOLCHAIN=;' export 'FEATURE_OPENMPI_MPI_COMPILER_TOOLCHAIN;' 'CCC_CFLAGS_modshare=-I/ccc/products/icc-17.0.6.256/system/default/17.0.6.256/include:2:-I/ccc/products/mkl-17.0.6.256/intel--17.0.6.256__openmpi--2.0.4/default/17.0.6.256/mkl/include:1:-I/ccc/products/geos-3.4.2/intel--17.0.4.196/default/include:1:-I/ccc/products/mkl-17.0.6.256/intel--17.0.6.256__openmpi--2.0.4/default/17.0.6.256/mkl/include/intel64/lp64:1;' export 'CCC_CFLAGS_modshare;' 'OPENMPI_ROOT=/ccc/products/openmpi-2.0.4/intel--17.0.6.256/default;' export 'OPENMPI_ROOT;' 'CCC_F77=gfortran;' export 'CCC_F77;' 'OPENMPI_LIBDIR=/ccc/products/openmpi-2.0.4/intel--17.0.6.256/default/lib;' export 'OPENMPI_LIBDIR;' 'CCC_MPIFC=mpifort;' export 'CCC_MPIFC;' 'OPENMPI_INCDIR=/ccc/products/openmpi-2.0.4/intel--17.0.6.256/default/include;' export 'OPENMPI_INCDIR;' 'CCC_CFLAGS=-I/ccc/products/mkl-17.0.6.256/intel--17.0.6.256__openmpi--2.0.4/default/17.0.6.256/mkl/include\' '-I/ccc/products/mkl-17.0.6.256/intel--17.0.6.256__openmpi--2.0.4/default/17.0.6.256/mkl/include/intel64/lp64\' '-I/ccc/products/icc-17.0.6.256/system/default/17.0.6.256/include\' '-I/ccc/products/geos-3.4.2/intel--17.0.4.196/default/include;' export 'CCC_CFLAGS;' 'CCC_MPIF77=mpif77;' export 'CCC_MPIF77;' 'GEOS_LDFLAGS=-L/ccc/products/geos-3.4.2/intel--17.0.4.196/default/lib\' '-lgeos;' export 'GEOS_LDFLAGS;' 'OMPI_MCA_fs_lustre_stripe_size=1048576;' export 'OMPI_MCA_fs_lustre_stripe_size;' 'GEOS_LIBDIR=/ccc/products/geos-3.4.2/intel--17.0.4.196/default/lib;' export 'GEOS_LIBDIR;' 'GEOS_INCDIR=/ccc/products/geos-3.4.2/intel--17.0.4.196/default/include;' export 'GEOS_INCDIR;' '_LMSTICKY_=ccc/1.0:dfldatadir/own:flavor/buildtarget/x86_64:nco/4.6.0:cdo/1.7.2rc6:gnu/7.3.0:mpi/openmpi/2.0.4:r/3.4.4;' export '_LMSTICKY_;' 'CXX_GNU_ROOT=/ccc/products/gcc-7.3.0/system/default;' export 'CXX_GNU_ROOT;' 'CCC_F90=gfortran;' export 'CCC_F90;' 'OPENMPI_TOOLCHAIN=intel--17.0.6.256;' export 'OPENMPI_TOOLCHAIN;' 'FORTRAN_GNU_ROOT=/ccc/products/gcc-7.3.0/system/default;' export 'FORTRAN_GNU_ROOT;' 'MANPATH=/ccc/products/openmpi-2.0.4/intel--17.0.6.256/default/share/man:/ccc/products/gcc-7.3.0/system/default/share/man:/ccc/products/icc-17.0.6.256/system/default/17.0.6.256/man/common:/ccc/products/share/man:/usr/share/man/en:/usr/share/man:/ccc/products/ccc_users_env/man/en:/ccc/products/ifort-17.0.6.256/system/default/17.0.6.256/man/common:/ccc/products/nco-4.6.0/intel--17.0.4.196/default/share/man:/ccc/products/r-3.4.4/intel--17.0.4.196__openmpi--2.0.2/default/share/man;' export 'MANPATH;' 'OMPI_CC=icc;' export 'OMPI_CC;' 'GNU_TOOLCHAIN=;' export 'GNU_TOOLCHAIN;' 'MANPATH_modshare=/ccc/products/r-3.4.4/intel--17.0.4.196__openmpi--2.0.2/default/share/man:1:/ccc/products/openmpi-2.0.4/intel--17.0.6.256/default/share/man:1:/usr/share/man/en:1:/ccc/products/share/man:1:/ccc/products/nco-4.6.0/intel--17.0.4.196/default/share/man:1:/ccc/products/ccc_users_env/man/en:1:/ccc/products/icc-17.0.6.256/system/default/17.0.6.256/man/common:2:/usr/share/man:1:/ccc/products/gcc-7.3.0/system/default/share/man:3:/ccc/products/ifort-17.0.6.256/system/default/17.0.6.256/man/common:1;' export 'MANPATH_modshare;' 'CXX_GNU_INCDIR=/ccc/products/gcc-7.3.0/system/default/include/c++/7.3.0;' export 'CXX_GNU_INCDIR;' 'CCC_MPIF90=mpif90;' export 'CCC_MPIF90;' '_LMPREREQ__modshare=r/3.4.4\&mpi\&mkl\&geos\&gdal:1:mkl/17.0.6.256\&feature/mkl/lp64\|feature/mkl/ilp64\&feature/mkl/sequential\|feature/mkl/multi-threaded\&feature/mkl/single_node\|feature/mkl/mpi:1:datadir/ineris\&ccc:1:gnu/7.3.0\&c/gnu/7.3.0\&c++/gnu/7.3.0\&fortran/gnu/7.3.0:1:c++/gnu/7.3.0\&flavor/buildcompiler/gcc/7\|flavor/buildcompiler:1:datadir/own\&ccc\&datadir/ineris:1:dfldatadir/own\&datadir/own\&ccc\&datadir/ineris:1:geos/3.4.2\&intel:1:c/gnu/7.3.0\&flavor/buildcompiler/gcc/7\|flavor/buildcompiler:1:gdal/2.1.0\&proj:1:fortran/intel/17.0.6.256\&licsrv/intel:1:cdo/1.7.2rc6\&intel:1:mpi/openmpi/2.0.4\&feature/mkl/single_node\|feature/mkl/mpi\&flavor/buildmpi/openmpi/2.0\&feature/openmpi/mpi_compiler/intel\|feature/openmpi/mpi_compiler\&flavor/openmpi\&.tuning/openmpi/2.0.4\&intel\|gnu\|pgi:1:intel/17.0.6.256\&flavor/buildcompiler/intel/17\|flavor/buildcompiler\&mkl/17.0.6.256\|mkl\&c/intel/17.0.6.256\&c++/intel/17.0.6.256\&fortran/intel/17.0.6.256:1:c++/intel/17.0.6.256\&licsrv/intel:1:.tuning/openmpi/2.0.4\&feature/openmpi/net\&c/intel\&c++/intel\&fortran/intel:1:fortran/gnu/7.3.0\&flavor/buildcompiler/gcc/7\|flavor/buildcompiler:1:c/intel/17.0.6.256\&licsrv/intel:1:nco/4.6.0\&intel:1;' export '_LMPREREQ__modshare;' 'PROJ_LIBDIR=/ccc/products/proj-4.9.1/system/default/lib;' export 'PROJ_LIBDIR;' 'LD_LIBRARY_PATH_modshare=/ccc/products/gcc-7.3.0/system/default/lib:3:/ccc/products/mkl-17.0.6.256/intel--17.0.6.256__openmpi--2.0.4/default/17.0.6.256/mkl/lib/intel64:1:/ccc/products/proj-4.9.1/system/default/lib:1:/ccc/products/gdal-2.1.0/intel--17.0.4.196/default/lib:1:/ccc/products/openmpi-2.0.4/intel--17.0.6.256/default/lib/openmpi:1:/opt/intel/ifort-17.0.6.256/system/default/lib/intel64:1:/ccc/products/openmpi-2.0.4/intel--17.0.6.256/default/lib:1:/opt/intel/icc-17.0.6.256/system/default/lib/intel64:2:/ccc/products/gcc-7.3.0/system/default/lib64:3:/ccc/products/nco-4.6.0/intel--17.0.4.196/default/lib:1:/ccc/products/r-3.4.4/intel--17.0.4.196__openmpi--2.0.2/default/lib64/R/lib:1:/ccc/products/geos-3.4.2/intel--17.0.4.196/default/lib:1;' export 'LD_LIBRARY_PATH_modshare;' 'R_TOOLCHAIN=intel--17.0.4.196__openmpi--2.0.2;' export 'R_TOOLCHAIN;' 'CCC_CXX=g++;' export 'CCC_CXX;' 'PROJ_INCDIR=/ccc/products/proj-4.9.1/system/default/include;' export 'PROJ_INCDIR;' 'TUNING_OPENMPI_TOOLCHAIN=;' export 'TUNING_OPENMPI_TOOLCHAIN;' 'MXM_TLS=self,shm,ud,rc;' export 'MXM_TLS;' 'C_GNU_TOOLCHAIN=system;' export 'C_GNU_TOOLCHAIN;' 'GEOS_CXXFLAGS=-I/ccc/products/geos-3.4.2/intel--17.0.4.196/default/include;' export 'GEOS_CXXFLAGS;' 'PATH=/ccc/products/geos-3.4.2/intel--17.0.4.196/default/bin:/ccc/products/openmpi-2.0.4/intel--17.0.6.256/default/bin:/ccc/products/gcc-7.3.0/system/default/bin:/ccc/products/cdo-1.7.2rc6/intel--17.0.4.196/default/bin:/ccc/products/icc-17.0.6.256/system/default/17.0.6.256/bin/intel64:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:/opt/ibutils/bin:/ccc/products/ccc_users_env/bin:/ccc/products/ifort-17.0.6.256/system/default/17.0.6.256/bin/intel64:/ccc/products/nco-4.6.0/intel--17.0.4.196/default/bin:/ccc/cont004/home/ineris/gressena/bin:/ccc/products/proj-4.9.1/system/default/bin:/ccc/products/gdal-2.1.0/intel--17.0.4.196/default/bin:/ccc/products/r-3.4.4/intel--17.0.4.196__openmpi--2.0.2/default/bin;' export 'PATH;' 'CCC_LDFLAGS_modshare=-L/ccc/products/geos-3.4.2/intel--17.0.4.196/default/lib:1::1:-lmkl_intel_lp64:1:-lmkl_core:1:-lpthread:1:-lm:1:-lgeos:1:-L/ccc/products/mkl-17.0.6.256/intel--17.0.6.256__openmpi--2.0.4/default/17.0.6.256/mkl/lib/intel64:1:-lmkl_sequential:1;' export 'CCC_LDFLAGS_modshare;' 'R_LIBDIR=/ccc/products/r-3.4.4/intel--17.0.4.196__openmpi--2.0.2/default/lib64/R/lib;' export 'R_LIBDIR;' 'OMPI_F77=ifort;' export 'OMPI_F77;' 'OPENMPI_EXEDIR=/ccc/products/openmpi-2.0.4/intel--17.0.6.256/default/bin;' export 'OPENMPI_EXEDIR;' 'LD_LIBRARY_PATH=/ccc/products/geos-3.4.2/intel--17.0.4.196/default/lib:/ccc/products/openmpi-2.0.4/intel--17.0.6.256/default/lib/openmpi:/ccc/products/openmpi-2.0.4/intel--17.0.6.256/default/lib:/ccc/products/gcc-7.3.0/system/default/lib64:/ccc/products/gcc-7.3.0/system/default/lib:/opt/intel/ifort-17.0.6.256/system/default/lib/intel64:/opt/intel/icc-17.0.6.256/system/default/lib/intel64:/ccc/products/mkl-17.0.6.256/intel--17.0.6.256__openmpi--2.0.4/default/17.0.6.256/mkl/lib/intel64:/ccc/products/nco-4.6.0/intel--17.0.4.196/default/lib:/ccc/products/proj-4.9.1/system/default/lib:/ccc/products/gdal-2.1.0/intel--17.0.4.196/default/lib:/ccc/products/r-3.4.4/intel--17.0.4.196__openmpi--2.0.2/default/lib64/R/lib;' export 'LD_LIBRARY_PATH;' 'CCC_MPICXX=mpicxx;' export 'CCC_MPICXX;' 'CCC_MPIC=mpicc;' export 'CCC_MPIC;' 'MXM_LOG_LEVEL=error;' export 'MXM_LOG_LEVEL;' 'CCC_CXXFLAGS_modshare=-I/ccc/products/icc-17.0.6.256/system/default/17.0.6.256/include:2:-I/ccc/products/mkl-17.0.6.256/intel--17.0.6.256__openmpi--2.0.4/default/17.0.6.256/mkl/include:1:-I/ccc/products/geos-3.4.2/intel--17.0.4.196/default/include:1:-I/ccc/products/mkl-17.0.6.256/intel--17.0.6.256__openmpi--2.0.4/default/17.0.6.256/mkl/include/intel64/lp64:1;' export 'CCC_CXXFLAGS_modshare;' 'LOADEDMODULES=ccc/1.0:datadir/ineris:datadir/own:dfldatadir/own:flavor/buildtarget/x86_64:feature/mkl/single_node:flavor/buildcompiler/intel/17:feature/mkl/lp64:feature/mkl/sequential:mkl/17.0.6.256:licsrv/intel:c/intel/17.0.6.256:c++/intel/17.0.6.256:fortran/intel/17.0.6.256:intel/17.0.6.256:nco/4.6.0:cdo/1.7.2rc6:c/gnu/7.3.0:c++/gnu/7.3.0:fortran/gnu/7.3.0:gnu/7.3.0:flavor/buildmpi/openmpi/2.0:feature/openmpi/mpi_compiler/intel:flavor/openmpi/standard:feature/openmpi/net/auto:.tuning/openmpi/2.0.4:mpi/openmpi/2.0.4:geos/3.4.2:proj/4.9.1:gdal/2.1.0:r/3.4.4;' export 'LOADEDMODULES;' 'MPI_LIBDIR=/ccc/products/openmpi-2.0.4/intel--17.0.6.256/default/lib;' export 'MPI_LIBDIR;' 'GDAL_ROOT=/ccc/products/gdal-2.1.0/intel--17.0.4.196/default;' export 'GDAL_ROOT;' 'FORTRAN_GNU_TOOLCHAIN=system;' export 'FORTRAN_GNU_TOOLCHAIN;' 'GDAL_LIBDIR=/ccc/products/gdal-2.1.0/intel--17.0.4.196/default/lib;' export 'GDAL_LIBDIR;' 'MPI_INCDIR=/ccc/products/openmpi-2.0.4/intel--17.0.6.256/default/include;' export 'MPI_INCDIR;' 'GEOS_TOOLCHAIN=intel--17.0.4.196;' export 'GEOS_TOOLCHAIN;' 'GEOS_EXEDIR=/ccc/products/geos-3.4.2/intel--17.0.4.196/default/bin;' export 'GEOS_EXEDIR;' 'PROJ_ROOT=/ccc/products/proj-4.9.1/system/default;' export 'PROJ_ROOT;' 'CCC_CXXFLAGS=-I/ccc/products/mkl-17.0.6.256/intel--17.0.6.256__openmpi--2.0.4/default/17.0.6.256/mkl/include\' '-I/ccc/products/mkl-17.0.6.256/intel--17.0.6.256__openmpi--2.0.4/default/17.0.6.256/mkl/include/intel64/lp64\' '-I/ccc/products/icc-17.0.6.256/system/default/17.0.6.256/include\' '-I/ccc/products/geos-3.4.2/intel--17.0.4.196/default/include;' export 'CCC_CXXFLAGS;' 'GDAL_INCDIR=/ccc/products/gdal-2.1.0/intel--17.0.4.196/default/include;' export 'GDAL_INCDIR;' 'GEOS_ROOT=/ccc/products/geos-3.4.2/intel--17.0.4.196/default;' export 'GEOS_ROOT;' 'MAN_PATH_modshare=/ccc/products/proj-4.9.1/system/default/share/man:1;' export 'MAN_PATH_modshare;' 'PROJ_TOOLCHAIN=system;' export 'PROJ_TOOLCHAIN;' 'CCC_MPICC=mpicc;' export 'CCC_MPICC;' 'LOADEDMODULES_modshare=r/3.4.4:1:c++/intel/17.0.6.256:1:mkl/17.0.6.256:1:feature/openmpi/net/auto:1:flavor/openmpi/standard:1:c/gnu/7.3.0:1:c/intel/17.0.6.256:1:fortran/intel/17.0.6.256:1:proj/4.9.1:1:datadir/own:1:feature/mkl/lp64:1:feature/openmpi/mpi_compiler/intel:1:nco/4.6.0:1:flavor/buildtarget/x86_64:1:flavor/buildcompiler/intel/17:1:licsrv/intel:1:gnu/7.3.0:1:ccc/1.0:1:intel/17.0.6.256:1:feature/mkl/sequential:1:geos/3.4.2:1:fortran/gnu/7.3.0:1:datadir/ineris:1:dfldatadir/own:1:.tuning/openmpi/2.0.4:1:flavor/buildmpi/openmpi/2.0:1:feature/mkl/single_node:1:gdal/2.1.0:1:mpi/openmpi/2.0.4:1:cdo/1.7.2rc6:1:c++/gnu/7.3.0:1;' export 'LOADEDMODULES_modshare;' 'OMPI_FC=ifort;' export 'OMPI_FC;' 'GEOS_CFLAGS=-I/ccc/products/geos-3.4.2/intel--17.0.4.196/default/include;' export 'GEOS_CFLAGS;' '_LMPREREQ_=datadir/ineris\&ccc:datadir/own\&ccc\&datadir/ineris:dfldatadir/own\&datadir/own\&ccc\&datadir/ineris:mkl/17.0.6.256\&feature/mkl/lp64\|feature/mkl/ilp64\&feature/mkl/sequential\|feature/mkl/multi-threaded\&feature/mkl/single_node\|feature/mkl/mpi:c/intel/17.0.6.256\&licsrv/intel:c++/intel/17.0.6.256\&licsrv/intel:fortran/intel/17.0.6.256\&licsrv/intel:intel/17.0.6.256\&flavor/buildcompiler/intel/17\|flavor/buildcompiler\&mkl/17.0.6.256\|mkl\&c/intel/17.0.6.256\&c++/intel/17.0.6.256\&fortran/intel/17.0.6.256:nco/4.6.0\&intel:cdo/1.7.2rc6\&intel:c/gnu/7.3.0\&flavor/buildcompiler/gcc/7\|flavor/buildcompiler:c++/gnu/7.3.0\&flavor/buildcompiler/gcc/7\|flavor/buildcompiler:fortran/gnu/7.3.0\&flavor/buildcompiler/gcc/7\|flavor/buildcompiler:gnu/7.3.0\&c/gnu/7.3.0\&c++/gnu/7.3.0\&fortran/gnu/7.3.0:.tuning/openmpi/2.0.4\&feature/openmpi/net\&c/intel\&c++/intel\&fortran/intel:mpi/openmpi/2.0.4\&feature/mkl/single_node\|feature/mkl/mpi\&flavor/buildmpi/openmpi/2.0\&feature/openmpi/mpi_compiler/intel\|feature/openmpi/mpi_compiler\&flavor/openmpi\&.tuning/openmpi/2.0.4\&intel\|gnu\|pgi:geos/3.4.2\&intel:gdal/2.1.0\&proj:r/3.4.4\&mpi\&mkl\&geos\&gdal;' export '_LMPREREQ_;' 'C_GNU_ROOT=/ccc/products/gcc-7.3.0/system/default;' export 'C_GNU_ROOT;' 'GDAL_TOOLCHAIN=intel--17.0.4.196;' export 'GDAL_TOOLCHAIN;' 'MAN_PATH=/ccc/products/proj-4.9.1/system/default/share/man;' export 'MAN_PATH;' 'FLAVOR_OPENMPI_TOOLCHAIN=;' export 'FLAVOR_OPENMPI_TOOLCHAIN;' 'FEATURE_OPENMPI_NET_TOOLCHAIN=;' export 'FEATURE_OPENMPI_NET_TOOLCHAIN;' 'LIBRARY_PATH=/ccc/products/gcc-7.3.0/system/default/lib64;' export 'LIBRARY_PATH;' 'INFOPATH=/ccc/products/gcc-7.3.0/system/default/share/info;' export 'INFOPATH;' 'CCC_FC=gfortran;' export 'CCC_FC;' '_LMFILES_=/usr/share/modulefiles/environment/ccc/1.0:/usr/share/modulefiles/environment/datadir/ineris:/usr/share/modulefiles/environment/datadir/own:/usr/share/modulefiles/environment/dfldatadir/own:/usr/share/modulefiles/configurations/flavor/buildtarget/x86_64:/usr/share/modulefiles/configurations/feature/mkl/single_node:/usr/share/modulefiles/configurations/flavor/buildcompiler/intel/17:/usr/share/modulefiles/configurations/feature/mkl/lp64:/usr/share/modulefiles/configurations/feature/mkl/sequential:/usr/share/modulefiles/libraries/mkl/17.0.6.256:/usr/share/modulefiles/environment/licsrv/intel:/usr/share/modulefiles/compilers/c/intel/17.0.6.256:/usr/share/modulefiles/compilers/c++/intel/17.0.6.256:/usr/share/modulefiles/compilers/fortran/intel/17.0.6.256:/usr/share/modulefiles/compilers/intel/17.0.6.256:/usr/share/modulefiles/tools/nco/4.6.0:/usr/share/modulefiles/tools/cdo/1.7.2rc6:/usr/share/modulefiles/compilers/c/gnu/7.3.0:/usr/share/modulefiles/compilers/c++/gnu/7.3.0:/usr/share/modulefiles/compilers/fortran/gnu/7.3.0:/usr/share/modulefiles/compilers/gnu/7.3.0:/usr/share/modulefiles/configurations/flavor/buildmpi/openmpi/2.0:/usr/share/modulefiles/configurations/feature/openmpi/mpi_compiler/intel:/usr/share/modulefiles/configurations/flavor/openmpi/standard:/usr/share/modulefiles/configurations/feature/openmpi/net/auto:/usr/share/modulefiles/configurations/.tuning/openmpi/2.0.4:/usr/share/modulefiles/parallel/mpi/openmpi/2.0.4:/usr/share/modulefiles/libraries/geos/3.4.2:/usr/share/modulefiles/libraries/proj/4.9.1:/usr/share/modulefiles/libraries/gdal/2.1.0:/usr/share/modulefiles/tools/r/3.4.4;' export '_LMFILES_;' 'LIBRARY_PATH_modshare=/ccc/products/gcc-7.3.0/system/default/lib64:3;' export 'LIBRARY_PATH_modshare;' 'INFOPATH_modshare=/ccc/products/gcc-7.3.0/system/default/share/info:3;' export 'INFOPATH_modshare;' 'PROJ_EXEDIR=/ccc/products/proj-4.9.1/system/default/bin;' export 'PROJ_EXEDIR;' 'OMPI_CXX=icpc;' export 'OMPI_CXX;' '_LMSTICKY__modshare=r/3.4.4:1:nco/4.6.0:1:mpi/openmpi/2.0.4:1:cdo/1.7.2rc6:1:flavor/buildtarget/x86_64:1:dfldatadir/own:1:gnu/7.3.0:1:ccc/1.0:1;' export '_LMSTICKY__modshare;' 'OMPI_MCA_mca_base_param_files=/opt/ccc_etc/openmpi-2.0.4/openmpi-mca-params.conf;' export 'OMPI_MCA_mca_base_param_files;' '_LMCONFLICT_=dfldatadir/own\&dfldatadir:flavor/buildtarget/x86_64\&flavor/buildtarget:feature/mkl/single_node\&feature/mkl/mpi:flavor/buildcompiler/intel/17\&flavor/buildcompiler:mkl/17.0.6.256\&mkl:c/intel/17.0.6.256\&c/intel:c++/intel/17.0.6.256\&c++/intel:fortran/intel/17.0.6.256\&fortran/intel:intel/17.0.6.256\&intel:nco/4.6.0\&nco:cdo/1.7.2rc6\&cdo:c/gnu/7.3.0\&c/gnu:c++/gnu/7.3.0\&c++/gnu:fortran/gnu/7.3.0\&fortran/gnu:gnu/7.3.0\&gnu:flavor/buildmpi/openmpi/2.0\&flavor/buildmpi:feature/openmpi/mpi_compiler/intel\&feature/openmpi/mpi_compiler\&mpi/intelmpi\&mpi/mvapich2:feature/openmpi/net/auto\&mpi/intelmpi\&mpi/mvapich2\&feature/openmpi/net:.tuning/openmpi/2.0.4\&.tuning/openmpi:mpi/openmpi/2.0.4\&mpi/openmpi\&mpi/intelmpi:geos/3.4.2\&geos:proj/4.9.1\&proj:gdal/2.1.0\&gdal:r/3.4.4\&r;' export '_LMCONFLICT_;' 'MXM_LOG_FILE=/dev/null;' export 'MXM_LOG_FILE;' 'R_EXEDIR=/ccc/products/r-3.4.4/intel--17.0.4.196__openmpi--2.0.2/default/bin;' export 'R_EXEDIR;' 'ROMIO_HINTS=/opt/ccc_etc/openmpi-2.0.4/romio-hints;' export 'ROMIO_HINTS;' 'CCC_C=gcc;' export 'CCC_C;' '_LMCONFLICT__modshare=flavor/buildcompiler/intel/17\&flavor/buildcompiler:1:fortran/intel/17.0.6.256\&fortran/intel:1:mkl/17.0.6.256\&mkl:1:c++/intel/17.0.6.256\&c++/intel:1:intel/17.0.6.256\&intel:1:.tuning/openmpi/2.0.4\&.tuning/openmpi:1:nco/4.6.0\&nco:1:flavor/buildtarget/x86_64\&flavor/buildtarget:1:mpi/openmpi/2.0.4\&mpi/openmpi\&mpi/intelmpi:1:c/intel/17.0.6.256\&c/intel:1:feature/mkl/single_node\&feature/mkl/mpi:1:c++/gnu/7.3.0\&c++/gnu:1:feature/openmpi/mpi_compiler/intel\&feature/openmpi/mpi_compiler\&mpi/intelmpi\&mpi/mvapich2:1:gnu/7.3.0\&gnu:1:dfldatadir/own\&dfldatadir:1:r/3.4.4\&r:1:feature/openmpi/net/auto\&mpi/intelmpi\&mpi/mvapich2\&feature/openmpi/net:1:flavor/buildmpi/openmpi/2.0\&flavor/buildmpi:1:cdo/1.7.2rc6\&cdo:1:c/gnu/7.3.0\&c/gnu:1:gdal/2.1.0\&gdal:1:proj/4.9.1\&proj:1:fortran/gnu/7.3.0\&fortran/gnu:1:geos/3.4.2\&geos:1;' export '_LMCONFLICT__modshare;' '_LMFILES__modshare=/usr/share/modulefiles/libraries/geos/3.4.2:1:/usr/share/modulefiles/compilers/c/gnu/7.3.0:1:/usr/share/modulefiles/tools/cdo/1.7.2rc6:1:/usr/share/modulefiles/environment/datadir/ineris:1:/usr/share/modulefiles/environment/dfldatadir/own:1:/usr/share/modulefiles/configurations/feature/mkl/lp64:1:/usr/share/modulefiles/configurations/feature/openmpi/mpi_compiler/intel:1:/usr/share/modulefiles/tools/r/3.4.4:1:/usr/share/modulefiles/configurations/flavor/buildtarget/x86_64:1:/usr/share/modulefiles/configurations/flavor/buildcompiler/intel/17:1:/usr/share/modulefiles/libraries/gdal/2.1.0:1:/usr/share/modulefiles/configurations/feature/mkl/sequential:1:/usr/share/modulefiles/compilers/gnu/7.3.0:1:/usr/share/modulefiles/libraries/mkl/17.0.6.256:1:/usr/share/modulefiles/configurations/.tuning/openmpi/2.0.4:1:/usr/share/modulefiles/tools/nco/4.6.0:1:/usr/share/modulefiles/compilers/c++/intel/17.0.6.256:1:/usr/share/modulefiles/configurations/flavor/buildmpi/openmpi/2.0:1:/usr/share/modulefiles/compilers/fortran/gnu/7.3.0:1:/usr/share/modulefiles/configurations/feature/mkl/single_node:1:/usr/share/modulefiles/compilers/c/intel/17.0.6.256:1:/usr/share/modulefiles/libraries/proj/4.9.1:1:/usr/share/modulefiles/compilers/fortran/intel/17.0.6.256:1:/usr/share/modulefiles/environment/datadir/own:1:/usr/share/modulefiles/parallel/mpi/openmpi/2.0.4:1:/usr/share/modulefiles/environment/licsrv/intel:1:/usr/share/modulefiles/configurations/feature/openmpi/net/auto:1:/usr/share/modulefiles/configurations/flavor/openmpi/standard:1:/usr/share/modulefiles/environment/ccc/1.0:1:/usr/share/modulefiles/compilers/c++/gnu/7.3.0:1:/usr/share/modulefiles/compilers/intel/17.0.6.256:1;' export '_LMFILES__modshare;' 'CXX_GNU_TOOLCHAIN=system;' export 'CXX_GNU_TOOLCHAIN;' 'CCC_LDFLAGS=-L/ccc/products/mkl-17.0.6.256/intel--17.0.6.256__openmpi--2.0.4/default/17.0.6.256/mkl/lib/intel64\' '-lmkl_intel_lp64\' '-lmkl_core\' '-lmkl_sequential\' '-lpthread\' '-lm\' '\' '-L/ccc/products/geos-3.4.2/intel--17.0.4.196/default/lib\' '-lgeos;' export 'CCC_LDFLAGS;' 'R_ROOT=/ccc/products/r-3.4.4/intel--17.0.4.196__openmpi--2.0.2/default;' export 'R_ROOT;' 'FLAVOR_BUILDMPI_OPENMPI_TOOLCHAIN=;' export 'FLAVOR_BUILDMPI_OPENMPI_TOOLCHAIN;' 'MPI_TOOLCHAIN=intel--17.0.6.256;' export 'MPI_TOOLCHAIN;' 'MPI_EXEDIR=/ccc/products/openmpi-2.0.4/intel--17.0.6.256/default/bin;' export 'MPI_EXEDIR;' 'PATH_modshare=/ccc/products/geos-3.4.2/intel--17.0.4.196/default/bin:1:/usr/bin:1:/ccc/products/ifort-17.0.6.256/system/default/17.0.6.256/bin/intel64:1:/ccc/products/gcc-7.3.0/system/default/bin:3:/usr/local/bin:1:/ccc/products/proj-4.9.1/system/default/bin:1:/ccc/products/gdal-2.1.0/intel--17.0.4.196/default/bin:1:/opt/ibutils/bin:1:/ccc/products/cdo-1.7.2rc6/intel--17.0.4.196/default/bin:1:/ccc/products/r-3.4.4/intel--17.0.4.196__openmpi--2.0.2/default/bin:1:/ccc/products/openmpi-2.0.4/intel--17.0.6.256/default/bin:1:/ccc/cont004/home/ineris/gressena/bin:1:/sbin:1:/ccc/products/nco-4.6.0/intel--17.0.4.196/default/bin:1:/ccc/products/icc-17.0.6.256/system/default/17.0.6.256/bin/intel64:2:/ccc/products/ccc_users_env/bin:1:/usr/sbin:1:/usr/local/sbin:1;' export 'PATH_modshare;'
++ GDAL_EXEDIR=/ccc/products/gdal-2.1.0/intel--17.0.4.196/default/bin
++ export GDAL_EXEDIR
++ PORTALS4_CQ_MAXPIDS=2
++ export PORTALS4_CQ_MAXPIDS
++ MPI_ROOT=/ccc/products/openmpi-2.0.4/intel--17.0.6.256/default
++ export MPI_ROOT
++ FEATURE_OPENMPI_MPI_COMPILER_TOOLCHAIN=
++ export FEATURE_OPENMPI_MPI_COMPILER_TOOLCHAIN
++ CCC_CFLAGS_modshare=-I/ccc/products/icc-17.0.6.256/system/default/17.0.6.256/include:2:-I/ccc/products/mkl-17.0.6.256/intel--17.0.6.256__openmpi--2.0.4/default/17.0.6.256/mkl/include:1:-I/ccc/products/geos-3.4.2/intel--17.0.4.196/default/include:1:-I/ccc/products/mkl-17.0.6.256/intel--17.0.6.256__openmpi--2.0.4/default/17.0.6.256/mkl/include/intel64/lp64:1
++ export CCC_CFLAGS_modshare
++ OPENMPI_ROOT=/ccc/products/openmpi-2.0.4/intel--17.0.6.256/default
++ export OPENMPI_ROOT
++ CCC_F77=gfortran
++ export CCC_F77
++ OPENMPI_LIBDIR=/ccc/products/openmpi-2.0.4/intel--17.0.6.256/default/lib
++ export OPENMPI_LIBDIR
++ CCC_MPIFC=mpifort
++ export CCC_MPIFC
++ OPENMPI_INCDIR=/ccc/products/openmpi-2.0.4/intel--17.0.6.256/default/include
++ export OPENMPI_INCDIR
++ CCC_CFLAGS='-I/ccc/products/mkl-17.0.6.256/intel--17.0.6.256__openmpi--2.0.4/default/17.0.6.256/mkl/include -I/ccc/products/mkl-17.0.6.256/intel--17.0.6.256__openmpi--2.0.4/default/17.0.6.256/mkl/include/intel64/lp64 -I/ccc/products/icc-17.0.6.256/system/default/17.0.6.256/include -I/ccc/products/geos-3.4.2/intel--17.0.4.196/default/include'
++ export CCC_CFLAGS
++ CCC_MPIF77=mpif77
++ export CCC_MPIF77
++ GEOS_LDFLAGS='-L/ccc/products/geos-3.4.2/intel--17.0.4.196/default/lib -lgeos'
++ export GEOS_LDFLAGS
++ OMPI_MCA_fs_lustre_stripe_size=1048576
++ export OMPI_MCA_fs_lustre_stripe_size
++ GEOS_LIBDIR=/ccc/products/geos-3.4.2/intel--17.0.4.196/default/lib
++ export GEOS_LIBDIR
++ GEOS_INCDIR=/ccc/products/geos-3.4.2/intel--17.0.4.196/default/include
++ export GEOS_INCDIR
++ _LMSTICKY_=ccc/1.0:dfldatadir/own:flavor/buildtarget/x86_64:nco/4.6.0:cdo/1.7.2rc6:gnu/7.3.0:mpi/openmpi/2.0.4:r/3.4.4
++ export _LMSTICKY_
++ CXX_GNU_ROOT=/ccc/products/gcc-7.3.0/system/default
++ export CXX_GNU_ROOT
++ CCC_F90=gfortran
++ export CCC_F90
++ OPENMPI_TOOLCHAIN=intel--17.0.6.256
++ export OPENMPI_TOOLCHAIN
++ FORTRAN_GNU_ROOT=/ccc/products/gcc-7.3.0/system/default
++ export FORTRAN_GNU_ROOT
++ MANPATH=/ccc/products/openmpi-2.0.4/intel--17.0.6.256/default/share/man:/ccc/products/gcc-7.3.0/system/default/share/man:/ccc/products/icc-17.0.6.256/system/default/17.0.6.256/man/common:/ccc/products/share/man:/usr/share/man/en:/usr/share/man:/ccc/products/ccc_users_env/man/en:/ccc/products/ifort-17.0.6.256/system/default/17.0.6.256/man/common:/ccc/products/nco-4.6.0/intel--17.0.4.196/default/share/man:/ccc/products/r-3.4.4/intel--17.0.4.196__openmpi--2.0.2/default/share/man
++ export MANPATH
++ OMPI_CC=icc
++ export OMPI_CC
++ GNU_TOOLCHAIN=
++ export GNU_TOOLCHAIN
++ MANPATH_modshare=/ccc/products/r-3.4.4/intel--17.0.4.196__openmpi--2.0.2/default/share/man:1:/ccc/products/openmpi-2.0.4/intel--17.0.6.256/default/share/man:1:/usr/share/man/en:1:/ccc/products/share/man:1:/ccc/products/nco-4.6.0/intel--17.0.4.196/default/share/man:1:/ccc/products/ccc_users_env/man/en:1:/ccc/products/icc-17.0.6.256/system/default/17.0.6.256/man/common:2:/usr/share/man:1:/ccc/products/gcc-7.3.0/system/default/share/man:3:/ccc/products/ifort-17.0.6.256/system/default/17.0.6.256/man/common:1
++ export MANPATH_modshare
++ CXX_GNU_INCDIR=/ccc/products/gcc-7.3.0/system/default/include/c++/7.3.0
++ export CXX_GNU_INCDIR
++ CCC_MPIF90=mpif90
++ export CCC_MPIF90
++ _LMPREREQ__modshare='r/3.4.4&mpi&mkl&geos&gdal:1:mkl/17.0.6.256&feature/mkl/lp64|feature/mkl/ilp64&feature/mkl/sequential|feature/mkl/multi-threaded&feature/mkl/single_node|feature/mkl/mpi:1:datadir/ineris&ccc:1:gnu/7.3.0&c/gnu/7.3.0&c++/gnu/7.3.0&fortran/gnu/7.3.0:1:c++/gnu/7.3.0&flavor/buildcompiler/gcc/7|flavor/buildcompiler:1:datadir/own&ccc&datadir/ineris:1:dfldatadir/own&datadir/own&ccc&datadir/ineris:1:geos/3.4.2&intel:1:c/gnu/7.3.0&flavor/buildcompiler/gcc/7|flavor/buildcompiler:1:gdal/2.1.0&proj:1:fortran/intel/17.0.6.256&licsrv/intel:1:cdo/1.7.2rc6&intel:1:mpi/openmpi/2.0.4&feature/mkl/single_node|feature/mkl/mpi&flavor/buildmpi/openmpi/2.0&feature/openmpi/mpi_compiler/intel|feature/openmpi/mpi_compiler&flavor/openmpi&.tuning/openmpi/2.0.4&intel|gnu|pgi:1:intel/17.0.6.256&flavor/buildcompiler/intel/17|flavor/buildcompiler&mkl/17.0.6.256|mkl&c/intel/17.0.6.256&c++/intel/17.0.6.256&fortran/intel/17.0.6.256:1:c++/intel/17.0.6.256&licsrv/intel:1:.tuning/openmpi/2.0.4&feature/openmpi/net&c/intel&c++/intel&fortran/intel:1:fortran/gnu/7.3.0&flavor/buildcompiler/gcc/7|flavor/buildcompiler:1:c/intel/17.0.6.256&licsrv/intel:1:nco/4.6.0&intel:1'
++ export _LMPREREQ__modshare
++ PROJ_LIBDIR=/ccc/products/proj-4.9.1/system/default/lib
++ export PROJ_LIBDIR
++ LD_LIBRARY_PATH_modshare=/ccc/products/gcc-7.3.0/system/default/lib:3:/ccc/products/mkl-17.0.6.256/intel--17.0.6.256__openmpi--2.0.4/default/17.0.6.256/mkl/lib/intel64:1:/ccc/products/proj-4.9.1/system/default/lib:1:/ccc/products/gdal-2.1.0/intel--17.0.4.196/default/lib:1:/ccc/products/openmpi-2.0.4/intel--17.0.6.256/default/lib/openmpi:1:/opt/intel/ifort-17.0.6.256/system/default/lib/intel64:1:/ccc/products/openmpi-2.0.4/intel--17.0.6.256/default/lib:1:/opt/intel/icc-17.0.6.256/system/default/lib/intel64:2:/ccc/products/gcc-7.3.0/system/default/lib64:3:/ccc/products/nco-4.6.0/intel--17.0.4.196/default/lib:1:/ccc/products/r-3.4.4/intel--17.0.4.196__openmpi--2.0.2/default/lib64/R/lib:1:/ccc/products/geos-3.4.2/intel--17.0.4.196/default/lib:1
++ export LD_LIBRARY_PATH_modshare
++ R_TOOLCHAIN=intel--17.0.4.196__openmpi--2.0.2
++ export R_TOOLCHAIN
++ CCC_CXX=g++
++ export CCC_CXX
++ PROJ_INCDIR=/ccc/products/proj-4.9.1/system/default/include
++ export PROJ_INCDIR
++ TUNING_OPENMPI_TOOLCHAIN=
++ export TUNING_OPENMPI_TOOLCHAIN
++ MXM_TLS=self,shm,ud,rc
++ export MXM_TLS
++ C_GNU_TOOLCHAIN=system
++ export C_GNU_TOOLCHAIN
++ GEOS_CXXFLAGS=-I/ccc/products/geos-3.4.2/intel--17.0.4.196/default/include
++ export GEOS_CXXFLAGS
++ PATH=/ccc/products/geos-3.4.2/intel--17.0.4.196/default/bin:/ccc/products/openmpi-2.0.4/intel--17.0.6.256/default/bin:/ccc/products/gcc-7.3.0/system/default/bin:/ccc/products/cdo-1.7.2rc6/intel--17.0.4.196/default/bin:/ccc/products/icc-17.0.6.256/system/default/17.0.6.256/bin/intel64:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:/opt/ibutils/bin:/ccc/products/ccc_users_env/bin:/ccc/products/ifort-17.0.6.256/system/default/17.0.6.256/bin/intel64:/ccc/products/nco-4.6.0/intel--17.0.4.196/default/bin:/ccc/cont004/home/ineris/gressena/bin:/ccc/products/proj-4.9.1/system/default/bin:/ccc/products/gdal-2.1.0/intel--17.0.4.196/default/bin:/ccc/products/r-3.4.4/intel--17.0.4.196__openmpi--2.0.2/default/bin
++ export PATH
++ CCC_LDFLAGS_modshare=-L/ccc/products/geos-3.4.2/intel--17.0.4.196/default/lib:1::1:-lmkl_intel_lp64:1:-lmkl_core:1:-lpthread:1:-lm:1:-lgeos:1:-L/ccc/products/mkl-17.0.6.256/intel--17.0.6.256__openmpi--2.0.4/default/17.0.6.256/mkl/lib/intel64:1:-lmkl_sequential:1
++ export CCC_LDFLAGS_modshare
++ R_LIBDIR=/ccc/products/r-3.4.4/intel--17.0.4.196__openmpi--2.0.2/default/lib64/R/lib
++ export R_LIBDIR
++ OMPI_F77=ifort
++ export OMPI_F77
++ OPENMPI_EXEDIR=/ccc/products/openmpi-2.0.4/intel--17.0.6.256/default/bin
++ export OPENMPI_EXEDIR
++ LD_LIBRARY_PATH=/ccc/products/geos-3.4.2/intel--17.0.4.196/default/lib:/ccc/products/openmpi-2.0.4/intel--17.0.6.256/default/lib/openmpi:/ccc/products/openmpi-2.0.4/intel--17.0.6.256/default/lib:/ccc/products/gcc-7.3.0/system/default/lib64:/ccc/products/gcc-7.3.0/system/default/lib:/opt/intel/ifort-17.0.6.256/system/default/lib/intel64:/opt/intel/icc-17.0.6.256/system/default/lib/intel64:/ccc/products/mkl-17.0.6.256/intel--17.0.6.256__openmpi--2.0.4/default/17.0.6.256/mkl/lib/intel64:/ccc/products/nco-4.6.0/intel--17.0.4.196/default/lib:/ccc/products/proj-4.9.1/system/default/lib:/ccc/products/gdal-2.1.0/intel--17.0.4.196/default/lib:/ccc/products/r-3.4.4/intel--17.0.4.196__openmpi--2.0.2/default/lib64/R/lib
++ export LD_LIBRARY_PATH
++ CCC_MPICXX=mpicxx
++ export CCC_MPICXX
++ CCC_MPIC=mpicc
++ export CCC_MPIC
++ MXM_LOG_LEVEL=error
++ export MXM_LOG_LEVEL
++ CCC_CXXFLAGS_modshare=-I/ccc/products/icc-17.0.6.256/system/default/17.0.6.256/include:2:-I/ccc/products/mkl-17.0.6.256/intel--17.0.6.256__openmpi--2.0.4/default/17.0.6.256/mkl/include:1:-I/ccc/products/geos-3.4.2/intel--17.0.4.196/default/include:1:-I/ccc/products/mkl-17.0.6.256/intel--17.0.6.256__openmpi--2.0.4/default/17.0.6.256/mkl/include/intel64/lp64:1
++ export CCC_CXXFLAGS_modshare
++ LOADEDMODULES=ccc/1.0:datadir/ineris:datadir/own:dfldatadir/own:flavor/buildtarget/x86_64:feature/mkl/single_node:flavor/buildcompiler/intel/17:feature/mkl/lp64:feature/mkl/sequential:mkl/17.0.6.256:licsrv/intel:c/intel/17.0.6.256:c++/intel/17.0.6.256:fortran/intel/17.0.6.256:intel/17.0.6.256:nco/4.6.0:cdo/1.7.2rc6:c/gnu/7.3.0:c++/gnu/7.3.0:fortran/gnu/7.3.0:gnu/7.3.0:flavor/buildmpi/openmpi/2.0:feature/openmpi/mpi_compiler/intel:flavor/openmpi/standard:feature/openmpi/net/auto:.tuning/openmpi/2.0.4:mpi/openmpi/2.0.4:geos/3.4.2:proj/4.9.1:gdal/2.1.0:r/3.4.4
++ export LOADEDMODULES
++ MPI_LIBDIR=/ccc/products/openmpi-2.0.4/intel--17.0.6.256/default/lib
++ export MPI_LIBDIR
++ GDAL_ROOT=/ccc/products/gdal-2.1.0/intel--17.0.4.196/default
++ export GDAL_ROOT
++ FORTRAN_GNU_TOOLCHAIN=system
++ export FORTRAN_GNU_TOOLCHAIN
++ GDAL_LIBDIR=/ccc/products/gdal-2.1.0/intel--17.0.4.196/default/lib
++ export GDAL_LIBDIR
++ MPI_INCDIR=/ccc/products/openmpi-2.0.4/intel--17.0.6.256/default/include
++ export MPI_INCDIR
++ GEOS_TOOLCHAIN=intel--17.0.4.196
++ export GEOS_TOOLCHAIN
++ GEOS_EXEDIR=/ccc/products/geos-3.4.2/intel--17.0.4.196/default/bin
++ export GEOS_EXEDIR
++ PROJ_ROOT=/ccc/products/proj-4.9.1/system/default
++ export PROJ_ROOT
++ CCC_CXXFLAGS='-I/ccc/products/mkl-17.0.6.256/intel--17.0.6.256__openmpi--2.0.4/default/17.0.6.256/mkl/include -I/ccc/products/mkl-17.0.6.256/intel--17.0.6.256__openmpi--2.0.4/default/17.0.6.256/mkl/include/intel64/lp64 -I/ccc/products/icc-17.0.6.256/system/default/17.0.6.256/include -I/ccc/products/geos-3.4.2/intel--17.0.4.196/default/include'
++ export CCC_CXXFLAGS
++ GDAL_INCDIR=/ccc/products/gdal-2.1.0/intel--17.0.4.196/default/include
++ export GDAL_INCDIR
++ GEOS_ROOT=/ccc/products/geos-3.4.2/intel--17.0.4.196/default
++ export GEOS_ROOT
++ MAN_PATH_modshare=/ccc/products/proj-4.9.1/system/default/share/man:1
++ export MAN_PATH_modshare
++ PROJ_TOOLCHAIN=system
++ export PROJ_TOOLCHAIN
++ CCC_MPICC=mpicc
++ export CCC_MPICC
++ LOADEDMODULES_modshare=r/3.4.4:1:c++/intel/17.0.6.256:1:mkl/17.0.6.256:1:feature/openmpi/net/auto:1:flavor/openmpi/standard:1:c/gnu/7.3.0:1:c/intel/17.0.6.256:1:fortran/intel/17.0.6.256:1:proj/4.9.1:1:datadir/own:1:feature/mkl/lp64:1:feature/openmpi/mpi_compiler/intel:1:nco/4.6.0:1:flavor/buildtarget/x86_64:1:flavor/buildcompiler/intel/17:1:licsrv/intel:1:gnu/7.3.0:1:ccc/1.0:1:intel/17.0.6.256:1:feature/mkl/sequential:1:geos/3.4.2:1:fortran/gnu/7.3.0:1:datadir/ineris:1:dfldatadir/own:1:.tuning/openmpi/2.0.4:1:flavor/buildmpi/openmpi/2.0:1:feature/mkl/single_node:1:gdal/2.1.0:1:mpi/openmpi/2.0.4:1:cdo/1.7.2rc6:1:c++/gnu/7.3.0:1
++ export LOADEDMODULES_modshare
++ OMPI_FC=ifort
++ export OMPI_FC
++ GEOS_CFLAGS=-I/ccc/products/geos-3.4.2/intel--17.0.4.196/default/include
++ export GEOS_CFLAGS
++ _LMPREREQ_='datadir/ineris&ccc:datadir/own&ccc&datadir/ineris:dfldatadir/own&datadir/own&ccc&datadir/ineris:mkl/17.0.6.256&feature/mkl/lp64|feature/mkl/ilp64&feature/mkl/sequential|feature/mkl/multi-threaded&feature/mkl/single_node|feature/mkl/mpi:c/intel/17.0.6.256&licsrv/intel:c++/intel/17.0.6.256&licsrv/intel:fortran/intel/17.0.6.256&licsrv/intel:intel/17.0.6.256&flavor/buildcompiler/intel/17|flavor/buildcompiler&mkl/17.0.6.256|mkl&c/intel/17.0.6.256&c++/intel/17.0.6.256&fortran/intel/17.0.6.256:nco/4.6.0&intel:cdo/1.7.2rc6&intel:c/gnu/7.3.0&flavor/buildcompiler/gcc/7|flavor/buildcompiler:c++/gnu/7.3.0&flavor/buildcompiler/gcc/7|flavor/buildcompiler:fortran/gnu/7.3.0&flavor/buildcompiler/gcc/7|flavor/buildcompiler:gnu/7.3.0&c/gnu/7.3.0&c++/gnu/7.3.0&fortran/gnu/7.3.0:.tuning/openmpi/2.0.4&feature/openmpi/net&c/intel&c++/intel&fortran/intel:mpi/openmpi/2.0.4&feature/mkl/single_node|feature/mkl/mpi&flavor/buildmpi/openmpi/2.0&feature/openmpi/mpi_compiler/intel|feature/openmpi/mpi_compiler&flavor/openmpi&.tuning/openmpi/2.0.4&intel|gnu|pgi:geos/3.4.2&intel:gdal/2.1.0&proj:r/3.4.4&mpi&mkl&geos&gdal'
++ export _LMPREREQ_
++ C_GNU_ROOT=/ccc/products/gcc-7.3.0/system/default
++ export C_GNU_ROOT
++ GDAL_TOOLCHAIN=intel--17.0.4.196
++ export GDAL_TOOLCHAIN
++ MAN_PATH=/ccc/products/proj-4.9.1/system/default/share/man
++ export MAN_PATH
++ FLAVOR_OPENMPI_TOOLCHAIN=
++ export FLAVOR_OPENMPI_TOOLCHAIN
++ FEATURE_OPENMPI_NET_TOOLCHAIN=
++ export FEATURE_OPENMPI_NET_TOOLCHAIN
++ LIBRARY_PATH=/ccc/products/gcc-7.3.0/system/default/lib64
++ export LIBRARY_PATH
++ INFOPATH=/ccc/products/gcc-7.3.0/system/default/share/info
++ export INFOPATH
++ CCC_FC=gfortran
++ export CCC_FC
++ _LMFILES_=/usr/share/modulefiles/environment/ccc/1.0:/usr/share/modulefiles/environment/datadir/ineris:/usr/share/modulefiles/environment/datadir/own:/usr/share/modulefiles/environment/dfldatadir/own:/usr/share/modulefiles/configurations/flavor/buildtarget/x86_64:/usr/share/modulefiles/configurations/feature/mkl/single_node:/usr/share/modulefiles/configurations/flavor/buildcompiler/intel/17:/usr/share/modulefiles/configurations/feature/mkl/lp64:/usr/share/modulefiles/configurations/feature/mkl/sequential:/usr/share/modulefiles/libraries/mkl/17.0.6.256:/usr/share/modulefiles/environment/licsrv/intel:/usr/share/modulefiles/compilers/c/intel/17.0.6.256:/usr/share/modulefiles/compilers/c++/intel/17.0.6.256:/usr/share/modulefiles/compilers/fortran/intel/17.0.6.256:/usr/share/modulefiles/compilers/intel/17.0.6.256:/usr/share/modulefiles/tools/nco/4.6.0:/usr/share/modulefiles/tools/cdo/1.7.2rc6:/usr/share/modulefiles/compilers/c/gnu/7.3.0:/usr/share/modulefiles/compilers/c++/gnu/7.3.0:/usr/share/modulefiles/compilers/fortran/gnu/7.3.0:/usr/share/modulefiles/compilers/gnu/7.3.0:/usr/share/modulefiles/configurations/flavor/buildmpi/openmpi/2.0:/usr/share/modulefiles/configurations/feature/openmpi/mpi_compiler/intel:/usr/share/modulefiles/configurations/flavor/openmpi/standard:/usr/share/modulefiles/configurations/feature/openmpi/net/auto:/usr/share/modulefiles/configurations/.tuning/openmpi/2.0.4:/usr/share/modulefiles/parallel/mpi/openmpi/2.0.4:/usr/share/modulefiles/libraries/geos/3.4.2:/usr/share/modulefiles/libraries/proj/4.9.1:/usr/share/modulefiles/libraries/gdal/2.1.0:/usr/share/modulefiles/tools/r/3.4.4
++ export _LMFILES_
++ LIBRARY_PATH_modshare=/ccc/products/gcc-7.3.0/system/default/lib64:3
++ export LIBRARY_PATH_modshare
++ INFOPATH_modshare=/ccc/products/gcc-7.3.0/system/default/share/info:3
++ export INFOPATH_modshare
++ PROJ_EXEDIR=/ccc/products/proj-4.9.1/system/default/bin
++ export PROJ_EXEDIR
++ OMPI_CXX=icpc
++ export OMPI_CXX
++ _LMSTICKY__modshare=r/3.4.4:1:nco/4.6.0:1:mpi/openmpi/2.0.4:1:cdo/1.7.2rc6:1:flavor/buildtarget/x86_64:1:dfldatadir/own:1:gnu/7.3.0:1:ccc/1.0:1
++ export _LMSTICKY__modshare
++ OMPI_MCA_mca_base_param_files=/opt/ccc_etc/openmpi-2.0.4/openmpi-mca-params.conf
++ export OMPI_MCA_mca_base_param_files
++ _LMCONFLICT_='dfldatadir/own&dfldatadir:flavor/buildtarget/x86_64&flavor/buildtarget:feature/mkl/single_node&feature/mkl/mpi:flavor/buildcompiler/intel/17&flavor/buildcompiler:mkl/17.0.6.256&mkl:c/intel/17.0.6.256&c/intel:c++/intel/17.0.6.256&c++/intel:fortran/intel/17.0.6.256&fortran/intel:intel/17.0.6.256&intel:nco/4.6.0&nco:cdo/1.7.2rc6&cdo:c/gnu/7.3.0&c/gnu:c++/gnu/7.3.0&c++/gnu:fortran/gnu/7.3.0&fortran/gnu:gnu/7.3.0&gnu:flavor/buildmpi/openmpi/2.0&flavor/buildmpi:feature/openmpi/mpi_compiler/intel&feature/openmpi/mpi_compiler&mpi/intelmpi&mpi/mvapich2:feature/openmpi/net/auto&mpi/intelmpi&mpi/mvapich2&feature/openmpi/net:.tuning/openmpi/2.0.4&.tuning/openmpi:mpi/openmpi/2.0.4&mpi/openmpi&mpi/intelmpi:geos/3.4.2&geos:proj/4.9.1&proj:gdal/2.1.0&gdal:r/3.4.4&r'
++ export _LMCONFLICT_
++ MXM_LOG_FILE=/dev/null
++ export MXM_LOG_FILE
++ R_EXEDIR=/ccc/products/r-3.4.4/intel--17.0.4.196__openmpi--2.0.2/default/bin
++ export R_EXEDIR
++ ROMIO_HINTS=/opt/ccc_etc/openmpi-2.0.4/romio-hints
++ export ROMIO_HINTS
++ CCC_C=gcc
++ export CCC_C
++ _LMCONFLICT__modshare='flavor/buildcompiler/intel/17&flavor/buildcompiler:1:fortran/intel/17.0.6.256&fortran/intel:1:mkl/17.0.6.256&mkl:1:c++/intel/17.0.6.256&c++/intel:1:intel/17.0.6.256&intel:1:.tuning/openmpi/2.0.4&.tuning/openmpi:1:nco/4.6.0&nco:1:flavor/buildtarget/x86_64&flavor/buildtarget:1:mpi/openmpi/2.0.4&mpi/openmpi&mpi/intelmpi:1:c/intel/17.0.6.256&c/intel:1:feature/mkl/single_node&feature/mkl/mpi:1:c++/gnu/7.3.0&c++/gnu:1:feature/openmpi/mpi_compiler/intel&feature/openmpi/mpi_compiler&mpi/intelmpi&mpi/mvapich2:1:gnu/7.3.0&gnu:1:dfldatadir/own&dfldatadir:1:r/3.4.4&r:1:feature/openmpi/net/auto&mpi/intelmpi&mpi/mvapich2&feature/openmpi/net:1:flavor/buildmpi/openmpi/2.0&flavor/buildmpi:1:cdo/1.7.2rc6&cdo:1:c/gnu/7.3.0&c/gnu:1:gdal/2.1.0&gdal:1:proj/4.9.1&proj:1:fortran/gnu/7.3.0&fortran/gnu:1:geos/3.4.2&geos:1'
++ export _LMCONFLICT__modshare
++ _LMFILES__modshare=/usr/share/modulefiles/libraries/geos/3.4.2:1:/usr/share/modulefiles/compilers/c/gnu/7.3.0:1:/usr/share/modulefiles/tools/cdo/1.7.2rc6:1:/usr/share/modulefiles/environment/datadir/ineris:1:/usr/share/modulefiles/environment/dfldatadir/own:1:/usr/share/modulefiles/configurations/feature/mkl/lp64:1:/usr/share/modulefiles/configurations/feature/openmpi/mpi_compiler/intel:1:/usr/share/modulefiles/tools/r/3.4.4:1:/usr/share/modulefiles/configurations/flavor/buildtarget/x86_64:1:/usr/share/modulefiles/configurations/flavor/buildcompiler/intel/17:1:/usr/share/modulefiles/libraries/gdal/2.1.0:1:/usr/share/modulefiles/configurations/feature/mkl/sequential:1:/usr/share/modulefiles/compilers/gnu/7.3.0:1:/usr/share/modulefiles/libraries/mkl/17.0.6.256:1:/usr/share/modulefiles/configurations/.tuning/openmpi/2.0.4:1:/usr/share/modulefiles/tools/nco/4.6.0:1:/usr/share/modulefiles/compilers/c++/intel/17.0.6.256:1:/usr/share/modulefiles/configurations/flavor/buildmpi/openmpi/2.0:1:/usr/share/modulefiles/compilers/fortran/gnu/7.3.0:1:/usr/share/modulefiles/configurations/feature/mkl/single_node:1:/usr/share/modulefiles/compilers/c/intel/17.0.6.256:1:/usr/share/modulefiles/libraries/proj/4.9.1:1:/usr/share/modulefiles/compilers/fortran/intel/17.0.6.256:1:/usr/share/modulefiles/environment/datadir/own:1:/usr/share/modulefiles/parallel/mpi/openmpi/2.0.4:1:/usr/share/modulefiles/environment/licsrv/intel:1:/usr/share/modulefiles/configurations/feature/openmpi/net/auto:1:/usr/share/modulefiles/configurations/flavor/openmpi/standard:1:/usr/share/modulefiles/environment/ccc/1.0:1:/usr/share/modulefiles/compilers/c++/gnu/7.3.0:1:/usr/share/modulefiles/compilers/intel/17.0.6.256:1
++ export _LMFILES__modshare
++ CXX_GNU_TOOLCHAIN=system
++ export CXX_GNU_TOOLCHAIN
++ CCC_LDFLAGS='-L/ccc/products/mkl-17.0.6.256/intel--17.0.6.256__openmpi--2.0.4/default/17.0.6.256/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm  -L/ccc/products/geos-3.4.2/intel--17.0.4.196/default/lib -lgeos'
++ export CCC_LDFLAGS
++ R_ROOT=/ccc/products/r-3.4.4/intel--17.0.4.196__openmpi--2.0.2/default
++ export R_ROOT
++ FLAVOR_BUILDMPI_OPENMPI_TOOLCHAIN=
++ export FLAVOR_BUILDMPI_OPENMPI_TOOLCHAIN
++ MPI_TOOLCHAIN=intel--17.0.6.256
++ export MPI_TOOLCHAIN
++ MPI_EXEDIR=/ccc/products/openmpi-2.0.4/intel--17.0.6.256/default/bin
++ export MPI_EXEDIR
++ PATH_modshare=/ccc/products/geos-3.4.2/intel--17.0.4.196/default/bin:1:/usr/bin:1:/ccc/products/ifort-17.0.6.256/system/default/17.0.6.256/bin/intel64:1:/ccc/products/gcc-7.3.0/system/default/bin:3:/usr/local/bin:1:/ccc/products/proj-4.9.1/system/default/bin:1:/ccc/products/gdal-2.1.0/intel--17.0.4.196/default/bin:1:/opt/ibutils/bin:1:/ccc/products/cdo-1.7.2rc6/intel--17.0.4.196/default/bin:1:/ccc/products/r-3.4.4/intel--17.0.4.196__openmpi--2.0.2/default/bin:1:/ccc/products/openmpi-2.0.4/intel--17.0.6.256/default/bin:1:/ccc/cont004/home/ineris/gressena/bin:1:/sbin:1:/ccc/products/nco-4.6.0/intel--17.0.4.196/default/bin:1:/ccc/products/icc-17.0.6.256/system/default/17.0.6.256/bin/intel64:2:/ccc/products/ccc_users_env/bin:1:/usr/sbin:1:/usr/local/sbin:1
++ export PATH_modshare
+ R CMD BATCH --no-save --no-restore mapping_urban_scale_EDK.r /ccc/work/cont004/ineris/gressena/microcapteurs_carto/MAPPING_SENSOR_DIR/03_Mapping/mapping_urban_scale_EDK.Rlog
+ exit 0
