# Cada archivo test_*.cpp es un ejecutable independiente registrado en CTest.
# El framework de tests es tests/geotest.hpp (sin dependencias externas).
function(geometra_add_test name)
    add_executable(${name} ${name}.cpp)
    target_link_libraries(${name} PRIVATE geometra::geometra)
    target_include_directories(${name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
    geometra_enable_warnings(${name})
    add_test(NAME ${name} COMMAND ${name})
    # Un predicado inconsistente no hace fallar a Delaunay: lo deja en un ciclo infinito.
    # Con este tope el cuelgue se reporta como fallo en vez de esperar los 1500 s de CTest.
    set_tests_properties(${name} PROPERTIES TIMEOUT 300)
endfunction()

geometra_add_test(test_point)
geometra_add_test(test_predicates)
geometra_add_test(test_convex_hull)
geometra_add_test(test_delaunay)
