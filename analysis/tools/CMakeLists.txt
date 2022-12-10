set(SND_TOOLS_SRCS
sndSciFiTools.cxx
)

add_library(sndTools SHARED ${SND_TOOLS_SRCS})
target_link_libraries(sndTools ${ROOT_LIBRARIES})
