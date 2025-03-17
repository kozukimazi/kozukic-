cmake_minimum_required(VERSION 3.10)
project(JordanWigner)

# Include Eigen headers
include_directories(/home/mazi/Desktop/eigen0/eigen-3.4.0)

# Add executable
add_executable(jw jordanwigner.cc jw.cc)
