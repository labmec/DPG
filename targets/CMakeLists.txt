add_executable(runDPG DPG-main.cpp)
target_link_libraries(runDPG PRIVATE DPG_lib)

add_executable(cooks-membrane-H1 cooks-membrane-H1.cpp)
target_link_libraries(cooks-membrane-H1 PRIVATE DPG_lib)
