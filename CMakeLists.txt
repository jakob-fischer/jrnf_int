project(odeint_rnet)

add_executable(odeint_rnet main.cpp)
set(CMAKE_CXX_FLAGS "-O3")

include_directories(/home/mit/csb/fischer/libs/include)
link_directories(/home/mit/csb/fischer/libs/lib)