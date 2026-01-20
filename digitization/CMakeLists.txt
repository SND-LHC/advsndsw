set(INCLUDE_DIRECTORIES
${BASE_INCLUDE_DIRECTORIES}
${CMAKE_SOURCE_DIR}/shipdata
${CMAKE_SOURCE_DIR}/shipLHC
${CMAKE_SOURCE_DIR}/digitization
)

include_directories( ${INCLUDE_DIRECTORIES}  ${VMC_INCLUDE_DIRS}  ${FAIRROOT_INCLUDE_DIR} ${FairLogger_INCDIR} ${FMT_INCLUDE_DIR}
)
include_directories(SYSTEM ${SYSTEM_INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES
${ROOT_LIBRARY_DIR}
${FAIRROOT_LIBRARY_DIR}
)
 
link_directories( ${LINK_DIRECTORIES})

set(SRCS
digitize.cxx
)

Set(HEADERS)
Set(LINKDEF digitizationLinkDef.h)
Set(LIBRARY_NAME digitization)
Set(DEPENDENCIES Base ShipData shipLHC GeoBase ParBase Geom Core)

GENERATE_LIBRARY()

add_executable(digitize_mc digitize_mc.cxx)
target_link_libraries(digitize_mc ${ROOT_LIBRARIES} digitization veto FairTools -lROOTTPython)
