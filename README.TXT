* Instructions for installation
** Development platform
Since this software depends on using the POSIX threads (pthreads), it
can only be built and run under a POSIX system.  This package is
developed on Ubuntu Linux 10.04 x86_64 and is tested on several
different versions of Ubuntu and RedHat Enterprise Linux 6.0.

Since Mac OS X also supports pthreads and MPI, this software should be
able to run on a cluster running such an operating system.  The
developers are aware of a port of pthreads on Microsoft Windows, but
we don't have a cluster computer running on Windows so we are not able
to test such a port.

** Dependency
 1. A working MPI environment, such as MPICH2 or Open MPI with both
    the binary files (mpiboot, mpiexec, etc) and necessary compilation
    support (header files, libraries, and compiler wrappers).
 2. A recent python interpreter (version >= 2.4).  This software is
    not tested under Python 3.0 though.
 3. Python C development files (header files).
 4. The pypar python library, which is used to execute MPI commands
    within python.  Note that this package depends on numpy, a Python
    library which contains many numerical computation routines.
 5. A standard C++ toolchain (e.g., g++ and GNU make) with pthread and
    fPIC support.
 6. The SWIG interface compiler to connect the C++ code with python.

** Install dependencies on a Ubuntu system (10.04 x86_64)
 - MPI :: Two implementations (MPICH2 and Open MPI) can be installed
          from the default repository.
          1. MPICH2.  Use "sudo apt-get install mpich2 libmpich2-dev"
          to install the binaries and headers/libraries/compiler
          wrappers.
          2. Open MPI.  Use "sudo apt-get install libopenmpi-dev
             openmpi-bin" to install the necessary packages.
 - Python :: Ubuntu comes with python as part of the base system.
             Tested version: 2.6.5.  You also need to install
             python-dev package for Python C development support.
 - C++ Compilation tools :: Use this command: "sudo apt-get install
      build-essential" to install a minimum GNU tool chain.
 - SWIG :: Use "sudo apt-get install swig" to install SWIG.  Tested
           version: 1.3.40.
 - numpy :: This Python library is required by pypar.  Install it by:
            "sudo apt-get install python-numpy".
 - pypar :: This python library is not in the Ubuntu repository yet.
   1. Download the source from pypar's homepage
  http://code.google.com/p/pypar/
   2. pypar works with Open MPI out of box.  It expects to call MPI
      library from a file called "libmpi.so.0" in your link path. If
      your MPI library has a different name, such as "libmpi.so" which
      is installed by MPICH2, you have to make a soft link so pypar
      can find this library instead

      sudo ln -s /usr/lib/libmpi.so /usr/lib/libmpi.so.0

      Alternatively, you can edit the pypar source code so it calls
      the right MPI library in your system.
   3. Change to pypar source directory and run the following two
      commands to install it:

      python compile_pypar_locally.py
      sudo python setup.py install

      This will install pypar into this directory:

      /usr/local/lib/python2.6/dist-packages/pypar

   4. make sure environment variable PYTHONPATH includes
      /usr/local/lib/python2.6/dist-packages/pypar.  Alternatively,
      you can create a file "local.pth" in
      /usr/local/lib/python2.6/dist-packages/, which contains just one
      word: pypar.  In this way you don't have to define PYTHONPATH
      manually. because /usr/local/lib/python2.6/dist-packages/ is
      included in the default python search path and any file with
      suffix ".pth" will be parsed and their content added to
      PYTHONPATH.  Tested pypar version: 2.1.4_94
** Compiling and installing this package (corrcperm)
  1. Change to the source directory.
  2. Optional: Use "swig -c++ -python -threads corrcperm.i" to generate file
     "corrcperm_wrap.cxx" if it is not already generated and included
     in the source package.
  3. Run "sudo python setup.py install" to build and install it to
     your system's default python library path, which usually is 
    /usr/local/lib/python2.6/dist-packages/

** Installing this package
   1. Use this command to install python library "corrcperm":

      sudo python setup.py install

      This will install corrcperm into this directory:

      /usr/local/lib/python2.6/dist-packages

   2. [Optional:] Copy the main program "mpi_nstat.py" to system
      search path, such as /usr/local/bin:

      sudo cp mpi_nstat.py /usr/local/bin

      This step is optional because you can run this script from
      anywhere. However, if you would like to run it on a Beowulf type
      of cluster without network shared filesystem, you have to make
      sure that all the dependencies (and corrcperm library) are
      installed on each one of the node, and "mpi_nstat.py" can be
      found in the system search path on each one of them.

* Instructions for running this program:
** On an Open MPI system

To run this program on one node, use the command:

mpiexec -np [number of mpi processes] ./mpi_nstat.py [arguments]

mpi_nstat.py takes a number of arguments, which *must* be terminated with the argument 'last'.  This is something of a hack, because mpiexec passes extra arguments from MPI to the python program, and we need a way to tell which arguments were specified by the user from the ones for MPI.

As an example, suppose this node has 8 computing cores.  To run this program with simulated (random) data on all processors for 20 permutations, use the following command: 

mpiexec -np 9 ./mpi_nstat.py -t --numthreads 1 --permutations 20 last  

After it finishes running, several log files (01.log, 02.log, etc)
will be created, together with a file "pvals", which is a Python
Pickle dump of a vector of p-values.  To load this file in Python:

from cPickle import load
pvalues = load(open('pvals'))

Remark 1: We choose 9 instead of 8 because the master MPI process only
takes care of sending and receiving computing tasks, so it would be
wasteful to assign one node for this process.

Remark 2: Please consult Open MPI's document for instructions about
how to run a job on multiple nodes (file sharing, SSH login, hostfile,
etc).

** On an MPICH2 system
On an MPICH2 system, first the MPI nodes must be started up using
command mpdboot.  For example to run on all 8 nodes, write a valid
hostfile and run:

mpdboot -n 8

After this we can use the same command to run this program, for example

mpiexec -np 33 ./mpi_nstat.py -t --numthreads 1 --permutations 20 last  

Finally, we have to manually turn off mpd.  On the master mode where mpdboot was started, type

mpdallexit  

to stop all mpd processes.

Remark: Again, please consult MPICH2's documentation for a detailed
instruction on running a job on multiple nodes.

** Description of the command line arguments
  -t :: specify whether you want to generate and run the program on simulated data (defaults to false)
  --genes :: number of genes (defaults to 7000, has no effect if -t is not specified)
  --columns :: number of slides in one condition (defaults to 80, has no effect if -t is not specified)
  --groups :: number of groups to use for the N-statistic (defaults to 8)
  --permutations :: number of permutations to run (defaults to 10)
  --kernel :: kernel to use, 1 is an identity kernel, 2 is an n^2 kernel (defaults to 2)
  --numthreads :: number of pthreads to use (defaults to 2)
  --seed :: starting seed for the random number generator (defaults to 12345)
  --file1 :: file that contains the pickle data for the first condition (defaults to 'Hyperdip_n.pck')
  --file2 :: file that contains the pickle data for the second condition (defaults to 'Tel_n.pck')
  --outfile :: filename to write the results to (defaults to pvals)
  --logprefix :: all log files will start with this prefix

** Description of the two included datasets
'Hyperdip_n.pck' and 'Tel_n.pck' are microarray expression data sampled from patients with two subtypes of childhood leukemia.  They are collected and made available by St. Jude Children's Research Hospital.  

To reduce the size of this package, we only included the first 500 genes. More information about this data can be found from:

[1] Rui Hu, Xing Qiu, Galina Glazko, Lev Klebanov, Andrei Yakovlev:
Detecting intergene correlation changes in microarray analysis: a new
approach to gene selection. BMC Bioinformatics 2009, 10:20.

[2] Yeoh EJ, Ross ME, Shurtleff SA, Williams WK, Patel D, Mahfouz R,
Behm FG, Raimondi SC, Relling MV, Patel A, Cheng C, Campana D, Wilkins
D, Zhou X, Li J, Liu H, Pui CH, Evans WE, Naeve C, Wong L, Downing JR:
Classification, subtype discovery, and prediction of outcome in
pediatric acute lymphoblastic leukemia by gene expression
profiling. Cancer Cell 2002, 1(2):133–143.
