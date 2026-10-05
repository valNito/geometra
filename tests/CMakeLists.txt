# Cada archivo test_*.cpp es un ejecutable independiente registrado en CTest.
# El framework de tests es tests/geotest.hpp (sin dependencias externas).
function(geometra_add_test name)
    add_executable(${name} ${name}.cpp)
    target_link_libraries(${name} PRIVATE geometra::geometra)
    target_include_directories(${name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
    geometra_enable_warnings(${name})
    add_test(NAME ${name} COMMAND ${name})
endfunction()

geometra_add_test(test_point)
geometra_add_test(test_predicates)
geometra_add_test(test_convex_hull)
