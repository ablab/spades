//***************************************************************************
//* Copyright (c) 2025 SPAdes team
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

/**
 * Unit tests for the error code system.
 *
 * Tests verify that:
 * 1. ErrorCode enum values are defined correctly (64-68)
 * 2. FATAL_*_ERROR macros exit with the correct error codes
 * 3. CHECK_FATAL_*_ERROR macros only exit when condition is false
 *
 * These tests use Google Test's death test feature (EXPECT_EXIT) to verify
 * that the macros exit with the expected codes.
 */

#include "utils/logger/logger.hpp"
#include "utils/logger/log_writers.hpp"
#include "utils/logger/error_codes.hpp"

#include <gtest/gtest.h>

// Helper to create a minimal logger for tests
void setup_test_logger() {
    using namespace logging;
    if (!__logger()) {
        logger *lg = create_logger("");
        lg->add_writer(std::make_shared<console_writer>());
        attach_logger(lg);
    }
}

//=============================================================================
// ErrorCode Enum Value Tests
//=============================================================================

class ErrorCodeValuesTest : public ::testing::Test {
protected:
    void SetUp() override {
        setup_test_logger();
    }
};

TEST_F(ErrorCodeValuesTest, InvalidInputFormatValue) {
    // InvalidInputFormat should be 64 (first user-end error code)
    EXPECT_EQ(ErrorCodes::InvalidInputFormat, 64);
}

TEST_F(ErrorCodeValuesTest, InputFileNotFoundValue) {
    // InputFileNotFound should be 65
    EXPECT_EQ(ErrorCodes::InputFileNotFound, 65);
}

TEST_F(ErrorCodeValuesTest, IOErrorValue) {
    // IOError should be 66
    EXPECT_EQ(ErrorCodes::IOError, 66);
}

TEST_F(ErrorCodeValuesTest, InvalidParameterValue) {
    // InvalidParameter should be 67
    EXPECT_EQ(ErrorCodes::InvalidParameter, 67);
}

TEST_F(ErrorCodeValuesTest, MemoryLimitExceededValue) {
    // MemoryLimitExceeded should be 68
    EXPECT_EQ(ErrorCodes::MemoryLimitExceeded, 68);
}

TEST_F(ErrorCodeValuesTest, UserEndErrorRange) {
    // All user-end error codes should be in range 64-127
    // This ensures they won't trigger "report bug" messages
    EXPECT_GE(ErrorCodes::InvalidInputFormat, 64);
    EXPECT_LE(ErrorCodes::InvalidInputFormat, 127);

    EXPECT_GE(ErrorCodes::InputFileNotFound, 64);
    EXPECT_LE(ErrorCodes::InputFileNotFound, 127);

    EXPECT_GE(ErrorCodes::IOError, 64);
    EXPECT_LE(ErrorCodes::IOError, 127);

    EXPECT_GE(ErrorCodes::InvalidParameter, 64);
    EXPECT_LE(ErrorCodes::InvalidParameter, 127);

    EXPECT_GE(ErrorCodes::MemoryLimitExceeded, 64);
    EXPECT_LE(ErrorCodes::MemoryLimitExceeded, 127);
}

//=============================================================================
// FATAL_*_ERROR Macro Death Tests
//=============================================================================

class FatalErrorMacrosDeathTest : public ::testing::Test {
protected:
    void SetUp() override {
        setup_test_logger();
    }
};

TEST_F(FatalErrorMacrosDeathTest, FatalErrorExitsWithMinusOne) {
    // FATAL_ERROR should exit with -1 (which becomes 255 as unsigned)
    EXPECT_EXIT(
        { FATAL_ERROR("Test fatal error"); },
        ::testing::ExitedWithCode(255),  // -1 as unsigned 8-bit
        ""
    );
}

TEST_F(FatalErrorMacrosDeathTest, FatalFormatErrorExitsWithCode64) {
    // FATAL_FORMAT_ERROR should exit with InvalidInputFormat (64)
    EXPECT_EXIT(
        { FATAL_FORMAT_ERROR("Test format error"); },
        ::testing::ExitedWithCode(64),
        ""
    );
}

TEST_F(FatalErrorMacrosDeathTest, FatalFileNotFoundErrorExitsWithCode65) {
    // FATAL_FILE_NOT_FOUND_ERROR should exit with InputFileNotFound (65)
    EXPECT_EXIT(
        { FATAL_FILE_NOT_FOUND_ERROR("Test file not found"); },
        ::testing::ExitedWithCode(65),
        ""
    );
}

TEST_F(FatalErrorMacrosDeathTest, FatalIOErrorExitsWithCode66) {
    // FATAL_IO_ERROR should exit with IOError (66)
    EXPECT_EXIT(
        { FATAL_IO_ERROR("Test IO error"); },
        ::testing::ExitedWithCode(66),
        ""
    );
}

TEST_F(FatalErrorMacrosDeathTest, FatalParamErrorExitsWithCode67) {
    // FATAL_PARAM_ERROR should exit with InvalidParameter (67)
    EXPECT_EXIT(
        { FATAL_PARAM_ERROR("Test param error"); },
        ::testing::ExitedWithCode(67),
        ""
    );
}

TEST_F(FatalErrorMacrosDeathTest, FatalMemExceededErrorExitsWithCode68) {
    // FATAL_MEM_EXCEEDED_ERROR should exit with MemoryLimitExceeded (68)
    EXPECT_EXIT(
        { FATAL_MEM_EXCEEDED_ERROR("Test memory exceeded"); },
        ::testing::ExitedWithCode(68),
        ""
    );
}

TEST_F(FatalErrorMacrosDeathTest, FatalErrorCodeExitsWithCustomCode) {
    // FATAL_ERROR_CODE should exit with the specified code
    EXPECT_EXIT(
        { FATAL_ERROR_CODE("Test custom error", 42); },
        ::testing::ExitedWithCode(42),
        ""
    );
}

//=============================================================================
// CHECK_FATAL_*_ERROR Macro Tests
//=============================================================================

class CheckFatalErrorMacrosTest : public ::testing::Test {
protected:
    void SetUp() override {
        setup_test_logger();
    }
};

TEST_F(CheckFatalErrorMacrosTest, CheckFatalErrorPassesOnTrue) {
    // CHECK_FATAL_ERROR should NOT exit when condition is true
    CHECK_FATAL_ERROR(true, "This should not trigger");
    SUCCEED();  // If we reach here, the test passes
}

TEST_F(CheckFatalErrorMacrosTest, CheckFatalFormatErrorPassesOnTrue) {
    // CHECK_FATAL_FORMAT_ERROR should NOT exit when condition is true
    CHECK_FATAL_FORMAT_ERROR(true, "This should not trigger");
    SUCCEED();
}

TEST_F(CheckFatalErrorMacrosTest, CheckFatalFileNotFoundPassesOnTrue) {
    // CHECK_FATAL_FILE_NOT_FOUND_ERROR should NOT exit when condition is true
    CHECK_FATAL_FILE_NOT_FOUND_ERROR(true, "This should not trigger");
    SUCCEED();
}

TEST_F(CheckFatalErrorMacrosTest, CheckFatalIOErrorPassesOnTrue) {
    // CHECK_FATAL_IO_ERROR should NOT exit when condition is true
    CHECK_FATAL_IO_ERROR(true, "This should not trigger");
    SUCCEED();
}

TEST_F(CheckFatalErrorMacrosTest, CheckFatalParamErrorPassesOnTrue) {
    // CHECK_FATAL_PARAM_ERROR should NOT exit when condition is true
    CHECK_FATAL_PARAM_ERROR(true, "This should not trigger");
    SUCCEED();
}

TEST_F(CheckFatalErrorMacrosTest, CheckFatalMemExceededPassesOnTrue) {
    // CHECK_FATAL_MEM_EXCEEDED_ERROR should NOT exit when condition is true
    CHECK_FATAL_MEM_EXCEEDED_ERROR(true, "This should not trigger");
    SUCCEED();
}

//=============================================================================
// CHECK_FATAL_*_ERROR Death Tests (when condition is false)
//=============================================================================

class CheckFatalErrorMacrosDeathTest : public ::testing::Test {
protected:
    void SetUp() override {
        setup_test_logger();
    }
};

TEST_F(CheckFatalErrorMacrosDeathTest, CheckFatalErrorExitsOnFalse) {
    // CHECK_FATAL_ERROR should exit with -1 when condition is false
    EXPECT_EXIT(
        { CHECK_FATAL_ERROR(false, "Condition failed"); },
        ::testing::ExitedWithCode(255),
        ""
    );
}

TEST_F(CheckFatalErrorMacrosDeathTest, CheckFatalFormatErrorExitsOnFalse) {
    // CHECK_FATAL_FORMAT_ERROR should exit with 64 when condition is false
    EXPECT_EXIT(
        { CHECK_FATAL_FORMAT_ERROR(false, "Format check failed"); },
        ::testing::ExitedWithCode(64),
        ""
    );
}

TEST_F(CheckFatalErrorMacrosDeathTest, CheckFatalFileNotFoundExitsOnFalse) {
    // CHECK_FATAL_FILE_NOT_FOUND_ERROR should exit with 65 when condition is false
    EXPECT_EXIT(
        { CHECK_FATAL_FILE_NOT_FOUND_ERROR(false, "File check failed"); },
        ::testing::ExitedWithCode(65),
        ""
    );
}

TEST_F(CheckFatalErrorMacrosDeathTest, CheckFatalIOErrorExitsOnFalse) {
    // CHECK_FATAL_IO_ERROR should exit with 66 when condition is false
    EXPECT_EXIT(
        { CHECK_FATAL_IO_ERROR(false, "IO check failed"); },
        ::testing::ExitedWithCode(66),
        ""
    );
}

TEST_F(CheckFatalErrorMacrosDeathTest, CheckFatalParamErrorExitsOnFalse) {
    // CHECK_FATAL_PARAM_ERROR should exit with 67 when condition is false
    EXPECT_EXIT(
        { CHECK_FATAL_PARAM_ERROR(false, "Param check failed"); },
        ::testing::ExitedWithCode(67),
        ""
    );
}

TEST_F(CheckFatalErrorMacrosDeathTest, CheckFatalMemExceededExitsOnFalse) {
    // CHECK_FATAL_MEM_EXCEEDED_ERROR should exit with 68 when condition is false
    EXPECT_EXIT(
        { CHECK_FATAL_MEM_EXCEEDED_ERROR(false, "Memory check failed"); },
        ::testing::ExitedWithCode(68),
        ""
    );
}

TEST_F(CheckFatalErrorMacrosDeathTest, CheckFatalErrorCodeExitsOnFalse) {
    // CHECK_FATAL_ERROR_CODE should exit with custom code when condition is false
    EXPECT_EXIT(
        { CHECK_FATAL_ERROR_CODE(false, "Custom check failed", 99); },
        ::testing::ExitedWithCode(99),
        ""
    );
}

//=============================================================================
// Expression Evaluation Tests
//=============================================================================

class CheckMacroExpressionTest : public ::testing::Test {
protected:
    void SetUp() override {
        setup_test_logger();
    }
};

TEST_F(CheckMacroExpressionTest, CheckMacroEvaluatesExpression) {
    // Verify that CHECK macros properly evaluate expressions
    int x = 5;
    CHECK_FATAL_ERROR(x > 0, "x should be positive");
    CHECK_FATAL_FORMAT_ERROR(x == 5, "x should be 5");
    CHECK_FATAL_PARAM_ERROR(x < 10, "x should be less than 10");
    SUCCEED();
}

TEST_F(CheckMacroExpressionTest, CheckMacroWithComplexExpression) {
    // Verify CHECK macros work with more complex expressions
    std::string filename = "test.txt";
    CHECK_FATAL_FILE_NOT_FOUND_ERROR(!filename.empty(), "filename should not be empty");

    int size = 100;
    int limit = 1024;
    CHECK_FATAL_MEM_EXCEEDED_ERROR(size < limit, "size should be within limit");

    SUCCEED();
}
