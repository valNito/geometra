# Warnings de desarrollo para los targets propios de Geometra.
# Son PRIVATE: no se propagan a quien consuma la librería, y solo se activan
# cuando Geometra es el proyecto raíz.
function(geometra_enable_warnings target)
    if(NOT GEOMETRA_IS_TOP_LEVEL)
        return()
    endif()

    target_compile_options(${target} PRIVATE
        $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:
            -Wall;-Wextra;-Wpedantic;-Wconversion;-Wsign-conversion;-Wshadow;
            -Wold-style-cast;-Wcast-align;-Wdouble-promotion;-Wnull-dereference;
            -Wnon-virtual-dtor;-Woverloaded-virtual>
        $<$<CXX_COMPILER_ID:MSVC>:/W4;/permissive->)

    if(GEOMETRA_WARNINGS_AS_ERRORS)
        target_compile_options(${target} PRIVATE
            $<IF:$<CXX_COMPILER_ID:MSVC>,/WX,-Werror>)
    endif()
endfunction()
