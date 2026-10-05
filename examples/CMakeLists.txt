add_executable(basic_usage basic_usage.cpp)
target_link_libraries(basic_usage PRIVATE geometra::geometra)
geometra_enable_warnings(basic_usage)
