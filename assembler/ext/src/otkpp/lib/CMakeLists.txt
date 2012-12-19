
SET(INST_HEADERS "Cloneable.h"
                 "OTK.h")

SET(OTKPP_SRCS "../constraints/BoundConstraints.cpp"
               "../objfunc/EvalCountable.cpp"
               "../objfunc/FDiffGradientEvaluator.cpp"
               "../objfunc/FDiffHessianEvaluator.cpp"
               "../objfunc/Function.cpp"
               "../objfunc/FunctionEvaluator.cpp"
               "../objfunc/GradientEvaluator.cpp"
               "../objfunc/HessianEvaluator.cpp"
               "../interpolation/QuadInterp.cpp"
               "../linalg/BFGSUpdater.cpp"
               "../linalg/LBFGSUpdater.cpp"
               "../linalg/InvLBFGSUpdater.cpp"
               "../linalg/SR1Updater.cpp"
               "../solvers/GradientSolverBase.cpp"
               "../solvers/Solver.cpp"
               "../solvers/native/linmin/Fletcher.cpp"
               "../solvers/native/linmin/LineMinimizer.cpp"
               "../solvers/native/linmin/MoreThuente.cpp"
               "../solvers/native/linmin/PolyFit.cpp"
               "../solvers/native/trsreg/DoglegSolver.cpp"
               "../solvers/native/trsreg/SteihaugSolver.cpp"
               "../solvers/native/trsreg/TrustRegionSolver.cpp"
               "../solvers/native/ConjGradMT.cpp"
               "../solvers/native/DoglegBFGS.cpp"
               "../solvers/native/DSQA.cpp"
               "../solvers/native/HookeJeeves.cpp"
               "../solvers/native/LinminBFGS.cpp"
               "../solvers/native/LRWWSimplex.cpp"
               "../solvers/native/MNewton.cpp"
               "../solvers/native/NativeGradientSolver.cpp"
               "../solvers/native/NativeSolver.cpp"
               "../solvers/native/PARTAN.cpp"
               "../solvers/native/Rosenbrock.cpp"
               "../solvers/native/SteihaugSR1.cpp"
               "../stopcrit/CompoundStoppingCriterion.cpp"
               "../stopcrit/FDistToMinTest.cpp"
               "../stopcrit/GradNormTest.cpp"
               "../stopcrit/MaxNumIterTest.cpp"
               "../stopcrit/StoppingCriterion.cpp"
               "../stopcrit/XDistToMinTest.cpp"
               "../testproblems/MGHTestFunction.cpp")

SET(OTKPP_GSL_SRCS "../solvers/gslsolvers/GSLFSolver.cpp"
                   "../solvers/gslsolvers/GSLFDFSolver.cpp"
                   "../solvers/gslsolvers/GSLUtils.cpp")

IF(WITH_GSL)
  SET(OTKPP_SRCS ${OTKPP_SRCS} ${OTKPP_GSL_SRCS})
ENDIF(WITH_GSL)

SET(OTKPP_MATHEVAL_SRCS "../objfunc/SymbolicFunctionEvaluator.cpp"
                        "../objfunc/SymbolicGradientEvaluator.cpp"
                        "../objfunc/SymbolicHessianEvaluator.cpp")

IF(WITH_LIBMATHEVAL)
  SET(OTKPP_SRCS ${OTKPP_SRCS} ${OTKPP_MATHEVAL_SRCS})
ENDIF(WITH_LIBMATHEVAL)

SET(OTKPP_FORTRAN_SRCS "../solvers/lbfgsb/LBFGSB.cpp"
                       #"../solvers/lmbm/func_wrapper.cpp"
                       #"../solvers/lmbm/LMBM.cpp"
                       #"../solvers/newuoa/NEWUOA.cpp"
                       #"../solvers/newuoa/func_wrapper.cpp"
                       "../solvers/lbfgsb/routines.f"
                       #"../solvers/lmbm/lmbm.f"
                       #"../solvers/lmbm/lmsub.f"
                       #"../solvers/lmbm/matcal.f"
                       #"../solvers/lmbm/tlmbm.f"
                       #"../solvers/lmbm/tnsunc.f"
                       #"../solvers/newuoa/bigden.f"
                       #"../solvers/newuoa/biglag.f"
                       #"../solvers/newuoa/newuoa.f"
                       #"../solvers/newuoa/newuob.f"
                       #"../solvers/newuoa/trsapp.f"
                       #"../solvers/newuoa/update.f"
                       )

INCLUDE_DIRECTORIES(${OTKPP_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIR})

ADD_LIBRARY(otkpp SHARED OTK.cpp ${OTKPP_SRCS})
INSTALL(TARGETS otkpp LIBRARY DESTINATION lib)

IF(WITH_FORTRAN)
  ADD_LIBRARY(otkpp_fortran SHARED ${OTKPP_FORTRAN_SRCS})
  INSTALL(TARGETS otkpp_fortran LIBRARY DESTINATION lib)
  SET_TARGET_PROPERTIES(otkpp_fortran PROPERTIES LINKER_LANGUAGE CXX)
ENDIF(WITH_FORTRAN)

INSTALL(FILES ${INST_HEADERS} DESTINATION include/otkpp/lib)
