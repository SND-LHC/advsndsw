# Create a library called "libShipGen" which includes the source files given in
# the array .
# The extension is already found.  Any number of sources could be listed here.

execute_process(COMMAND genie-config --libs OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE GENIE_LIBS)
execute_process(COMMAND genie-config --topsrcdir OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE GENIE_SRCDIR)

set(INCLUDE_DIRECTORIES
${CMAKE_SOURCE_DIR}/shipgen
${CMAKE_SOURCE_DIR}/generators
${CMAKE_SOURCE_DIR}/shipdata
${CMAKE_SOURCE_DIR}/veto
)

set(SYSTEM_INCLUDE_DIRECTORIES
${SYSTEM_INCLUDE_DIRECTORIES}
${PYTHIA8_INCLUDE_DIR}
${PYTHIA8_INCLUDE_DIR}/..
${HEPMC_INCLUDE_DIR}
${EVTGEN_INCLUDE_DIR}
${EVTGEN_INCLUDE_DIR}/..
${GENIE_SRCDIR}
)

include_directories(${INCLUDE_DIRECTORIES} ${VMC_INCLUDE_DIRS} ${FAIRROOT_INCLUDE_DIR} ${FairLogger_INCDIR} ${FMT_INCLUDE_DIR})
include_directories(SYSTEM ${SYSTEM_INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES
${ROOT_LIBRARY_DIR}
${FAIRROOT_LIBRARY_DIR}
${PYTHIA8_LIBRARY_DIR}
${EVTGEN_LIBRARY_DIR}
${SIMPATH}/lib
$ENV{APFEL_ROOT}/lib
$ENV{LOG4_ROOT}/lib
$ENV{LHAPDF_ROOT}/lib
$ENV{LIBXML2_ROOT}/lib
)

link_directories(${LINK_DIRECTORIES})

set(SRCS
Pythia8Generator.cxx
Pythia6Generator.cxx
tPythia6Generator.cxx
HNLPythia8Generator.cxx
DPPythia8Generator.cxx
GenieGenerator.cxx
NuageGenerator.cxx
NtupleGenerator.cxx
NtupleGenerator_FLUKA.cxx
NeutronGenerator_FLUKA.cxx
MuonBackGenerator.cxx
CosmicsGenerator.cxx
MuDISGenerator.cxx
FixedTargetGenerator.cxx
ALPACAGenerator.cxx
)

set(LINKDEF GenLinkDef.h)
set(LIBRARY_NAME ShipGen)

set(DEPENDENCIES Base ShipData EvtGen EvtGenExternal)

GENERATE_LIBRARY()

add_executable(sndlhcGSimpleNtpConverter sndlhcGSimpleNtpConverter.cxx)

# g++ sndlhcGSimpleNtpConverter.cxx -o sndlhcGSimpleNtpConverter `genie-config --libs` `root-config --libs` -lGeom -lEGPythia6 -lMathMore -L${LIBXML2_ROOT}/lib -lxml2 -L${APFEL_ROOT}/lib -lAPFEL -L${LOG4CPP_ROOT}/lib -llog4cpp -L${LHAPDF_ROOT}/lib -lLHAPDF

target_link_libraries(sndlhcGSimpleNtpConverter ${GENIE_LIBS} ${ROOT_LIBRARIES} -lGeom -lEGPythia6 -lMathMore -lxml2 -lAPFEL -llog4cpp -lLHAPDF)