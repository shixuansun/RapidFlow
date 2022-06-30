set(UTILITY_CORE_SRC
        relation/global_edge_view.h
        relation/local_edge_view.h
        graphoperations.cpp commandparser.cpp commandparser.h
        computesetintersection.h computesetintersection.cpp QFilter.h
        global_variables.h global_variables.cpp
        mem_usage.cpp mem_usage.h)

add_library(utility_core
        ${UTILITY_CORE_SRC})
