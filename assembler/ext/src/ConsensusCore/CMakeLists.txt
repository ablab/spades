project(ConsensusCore CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${EXT_DIR}/include/ConsensusCore)

add_library(ConsensusCore STATIC 
            AffineAlignment.cpp
            Coverage.cpp
            Mutation.cpp
            PairwiseAlignment.cpp
            Sequence.cpp
            Utils.cpp
            Version.cpp
            Matrix/DenseMatrix.cpp
            Matrix/SparseMatrix.cpp
            Poa/PoaConfig.cpp
            Poa/PoaConsensus.cpp
            Poa/PoaGraph.cpp)
