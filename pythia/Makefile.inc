# PYTHIA configuration file.
# Generated on Mon Nov 9 20:00:00 EET 2020 with the user supplied options:
# --with-root=/home/heimarry/alice/sw/slc7_x86-64/ROOT/v6-20-02-alice7-1

# Install directory prefixes.
PREFIX_BIN=/home/heimarry/Applications/pythia8303/bin
PREFIX_INCLUDE=/home/heimarry/Applications/pythia8303/include
PREFIX_LIB=/home/heimarry/Applications/pythia8303/lib
PREFIX_SHARE=/home/heimarry/Applications/pythia8303/share/Pythia8

# Compilation flags (see ./configure --help for further documentation).
CXX=g++
CXX_COMMON=-O2 -std=c++11 -pedantic -W -Wall -Wshadow -fPIC
CXX_SHARED=-shared
CXX_SONAME=-Wl,-soname,
LIB_SUFFIX=.so
OBJ_COMMON=

EVTGEN_USE=false
EVTGEN_CONFIG=
EVTGEN_BIN=
EVTGEN_INCLUDE=
EVTGEN_LIB=

FASTJET3_USE=false
FASTJET3_CONFIG=
FASTJET3_BIN=
FASTJET3_INCLUDE=
FASTJET3_LIB=

HEPMC2_USE=false
HEPMC2_CONFIG=
HEPMC2_BIN=
HEPMC2_INCLUDE=
HEPMC2_LIB=

HEPMC3_USE=false
HEPMC3_CONFIG=
HEPMC3_BIN=
HEPMC3_INCLUDE=
HEPMC3_LIB=

LHAPDF5_USE=false
LHAPDF5_CONFIG=
LHAPDF5_BIN=
LHAPDF5_INCLUDE=
LHAPDF5_LIB=

LHAPDF6_USE=false
LHAPDF6_CONFIG=
LHAPDF6_BIN=
LHAPDF6_INCLUDE=
LHAPDF6_LIB=

POWHEG_USE=false
POWHEG_CONFIG=
POWHEG_BIN=
POWHEG_INCLUDE=
POWHEG_LIB=

RIVET_USE=false
RIVET_CONFIG=
RIVET_BIN=
RIVET_INCLUDE=
RIVET_LIB=

ROOT_USE=true
ROOT_BIN=/home/heimarry/alice/sw/slc7_x86-64/ROOT/v6-20-02-alice7-1/bin/
ROOT_CONFIG=$(ROOT_BIN)root-config
ROOT_INCLUDE=-I/home/heimarry/alice/sw/slc7_x86-64/ROOT/v6-20-02-alice7-1/include
ROOT_LIB=-L/home/heimarry/alice/sw/slc7_x86-64/ROOT/v6-20-02-alice7-1/lib -Wl,-rpath,/home/heimarry/alice/sw/slc7_x86-64/ROOT/v6-20-02-alice7-1/lib -lCore

YODA_USE=false
YODA_CONFIG=root-config
YODA_BIN=
YODA_INCLUDE=
YODA_LIB=

GZIP_USE=false
GZIP_CONFIG=root-config
GZIP_BIN=
GZIP_INCLUDE=
GZIP_LIB=

PYTHON_USE=false
PYTHON_CONFIG=root-config
PYTHON_BIN=
PYTHON_INCLUDE=
PYTHON_LIB=

MG5MES_USE=false
MG5MES_CONFIG=root-config
MG5MES_BIN=
MG5MES_INCLUDE=
MG5MES_LIB=

OPENMP_USE=false
OPENMP_CONFIG=root-config
OPENMP_BIN=
OPENMP_INCLUDE=
OPENMP_LIB=
