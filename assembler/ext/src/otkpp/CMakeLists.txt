project(otkpp CXX C)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${EXT_DIR}/include)

SET(OTKPP_INCLUDE_DIRS "${EXT_DIR}/include/otkpp/constraints"
                       "${EXT_DIR}/include/otkpp/objfunc"
                       "${EXT_DIR}/include/otkpp/interpolation"
                       "${EXT_DIR}/include/otkpp/lib"
                       "${EXT_DIR}/include/otkpp/linalg"
                       "${EXT_DIR}/include/otkpp/solvers"
                       "${EXT_DIR}/include/otkpp/solvers/lbfgsb"
                       "${EXT_DIR}/include/otkpp/solvers/native"
                       "${EXT_DIR}/include/otkpp/solvers/native/linmin"
                       "${EXT_DIR}/include/otkpp/solvers/native/trsreg"
                       "${EXT_DIR}/include/otkpp/stopcrit")
include_directories(${OTKPP_INCLUDE_DIRS})

remove_definitions(-Wall -Wextra -Wconversion -Wno-long-long -Wwrite-strings -pedantic)

set(OTKPP_SRCS "constraints/BoundConstraints.cpp"
               "objfunc/EvalCountable.cpp"
               "objfunc/FDiffGradientEvaluator.cpp"
               "objfunc/FDiffHessianEvaluator.cpp"
               "objfunc/Function.cpp"
               "objfunc/FunctionEvaluator.cpp"
               "objfunc/GradientEvaluator.cpp"
               "objfunc/HessianEvaluator.cpp"
               "interpolation/QuadInterp.cpp"
               "linalg/BFGSUpdater.cpp"
               "linalg/LBFGSUpdater.cpp"
               "linalg/InvLBFGSUpdater.cpp"
               "linalg/SR1Updater.cpp"
               "solvers/GradientSolverBase.cpp"
               "solvers/Solver.cpp"
               "solvers/native/linmin/Fletcher.cpp"
               "solvers/native/linmin/LineMinimizer.cpp"
               "solvers/native/linmin/MoreThuente.cpp"
               "solvers/native/linmin/PolyFit.cpp"
               "solvers/native/trsreg/DoglegSolver.cpp"
               "solvers/native/trsreg/SteihaugSolver.cpp"
               "solvers/native/trsreg/TrustRegionSolver.cpp"
               "solvers/native/ConjGradMT.cpp"
               "solvers/native/DoglegBFGS.cpp"
               "solvers/native/DSQA.cpp"
               "solvers/native/HookeJeeves.cpp"
               "solvers/native/LinminBFGS.cpp"
               "solvers/native/LRWWSimplex.cpp"
               "solvers/native/MNewton.cpp"
               "solvers/native/NativeGradientSolver.cpp"
               "solvers/native/NativeSolver.cpp"
               "solvers/native/PARTAN.cpp"
               "solvers/native/SteihaugSR1.cpp"
               "stopcrit/CompoundStoppingCriterion.cpp"
               "stopcrit/FDistToMinTest.cpp"
               "stopcrit/GradNormTest.cpp"
               "stopcrit/MaxNumIterTest.cpp"
               "stopcrit/StoppingCriterion.cpp"
               "stopcrit/XDistToMinTest.cpp")

add_library(otkpp STATIC lib/OTK.cpp ${OTKPP_SRCS})
