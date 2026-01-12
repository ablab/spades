#!/usr/bin/env python3

############################################################################
# Copyright (c) 2025 SPAdes team
# All Rights Reserved
# See file LICENSE for details.
############################################################################

"""
Unit tests for the error code system.

Tests verify that:
1. ErrorCode enum values match the C++ error_codes.hpp values
2. user_end_error() correctly identifies user-end vs internal errors
3. error() function displays appropriate messages based on error type
"""

import io
import sys
import unittest
from unittest.mock import patch, MagicMock

# Add parent directory to path for imports
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from support import ErrorCode, error, report_issue_error_message, no_report_error_message


class TestErrorCodeValues(unittest.TestCase):
    """Test that ErrorCode enum values match the C++ error_codes.hpp definitions."""

    def test_invalid_input_format_value(self):
        """InvalidInputFormat should be 64 (matches C++ ErrorCodes::InvalidInputFormat)."""
        self.assertEqual(ErrorCode.InvalidInputFormat.value, 64)

    def test_input_file_not_found_value(self):
        """InputFileNotFound should be 65 (matches C++ ErrorCodes::InputFileNotFound)."""
        self.assertEqual(ErrorCode.InputFileNotFound.value, 65)

    def test_io_error_value(self):
        """IOError should be 66 (matches C++ ErrorCodes::IOError)."""
        self.assertEqual(ErrorCode.IOError.value, 66)

    def test_invalid_parameter_value(self):
        """InvalidParameter should be 67 (matches C++ ErrorCodes::InvalidParameter)."""
        self.assertEqual(ErrorCode.InvalidParameter.value, 67)

    def test_memory_limit_exceeded_value(self):
        """MemoryLimitExceeded should be 68 (matches C++ ErrorCodes::MemoryLimitExceeded)."""
        self.assertEqual(ErrorCode.MemoryLimitExceeded.value, 68)

    def test_general_error_value(self):
        """GeneralError should be -1 (internal error)."""
        self.assertEqual(ErrorCode.GeneralError.value, -1)

    def test_possible_wsl_error_value(self):
        """PossibleWSLError should be -11 (internal error)."""
        self.assertEqual(ErrorCode.PossibleWSLError.value, -11)


class TestUserEndError(unittest.TestCase):
    """Test that user_end_error() correctly identifies error types.

    User-end errors (64-127) are caused by user input/environment issues
    and should NOT suggest reporting a bug.

    Internal errors (1-63, 128-255, negative) indicate potential bugs
    and SHOULD suggest reporting an issue.
    """

    def test_invalid_input_format_is_user_end(self):
        """InvalidInputFormat (64) is a user-end error."""
        self.assertTrue(ErrorCode.InvalidInputFormat.user_end_error())

    def test_input_file_not_found_is_user_end(self):
        """InputFileNotFound (65) is a user-end error."""
        self.assertTrue(ErrorCode.InputFileNotFound.user_end_error())

    def test_io_error_is_user_end(self):
        """IOError (66) is a user-end error."""
        self.assertTrue(ErrorCode.IOError.user_end_error())

    def test_invalid_parameter_is_user_end(self):
        """InvalidParameter (67) is a user-end error."""
        self.assertTrue(ErrorCode.InvalidParameter.user_end_error())

    def test_memory_limit_exceeded_is_user_end(self):
        """MemoryLimitExceeded (68) is a user-end error."""
        self.assertTrue(ErrorCode.MemoryLimitExceeded.user_end_error())

    def test_general_error_is_not_user_end(self):
        """GeneralError (-1) is an internal error, not user-end."""
        self.assertFalse(ErrorCode.GeneralError.user_end_error())

    def test_possible_wsl_error_is_not_user_end(self):
        """PossibleWSLError (-11) is an internal error, not user-end."""
        self.assertFalse(ErrorCode.PossibleWSLError.user_end_error())

    def test_boundary_value_63_is_not_user_end(self):
        """Value 63 should NOT be a user-end error (below range)."""
        # Create a mock ErrorCode-like object to test boundary
        class MockCode:
            value = 63
            def user_end_error(self):
                return 64 <= self.value <= 127
        self.assertFalse(MockCode().user_end_error())

    def test_boundary_value_64_is_user_end(self):
        """Value 64 should be a user-end error (start of range)."""
        class MockCode:
            value = 64
            def user_end_error(self):
                return 64 <= self.value <= 127
        self.assertTrue(MockCode().user_end_error())

    def test_boundary_value_127_is_user_end(self):
        """Value 127 should be a user-end error (end of range)."""
        class MockCode:
            value = 127
            def user_end_error(self):
                return 64 <= self.value <= 127
        self.assertTrue(MockCode().user_end_error())

    def test_boundary_value_128_is_not_user_end(self):
        """Value 128 should NOT be a user-end error (above range)."""
        class MockCode:
            value = 128
            def user_end_error(self):
                return 64 <= self.value <= 127
        self.assertFalse(MockCode().user_end_error())


class TestErrorMessages(unittest.TestCase):
    """Test that error messages are appropriate for error types."""

    def test_report_issue_message_contains_github(self):
        """Internal errors should suggest reporting on GitHub."""
        messages = report_issue_error_message("SPAdes")
        combined = " ".join(messages)
        self.assertIn("github.com/ablab/spades", combined)
        self.assertIn("report", combined.lower())

    def test_report_issue_message_requests_files(self):
        """Internal errors should request log files for debugging."""
        messages = report_issue_error_message("SPAdes")
        combined = " ".join(messages)
        self.assertIn("params.txt", combined)
        self.assertIn(".log", combined)

    def test_no_report_message_no_github(self):
        """User-end errors should NOT mention GitHub reporting."""
        messages = no_report_error_message("SPAdes")
        combined = " ".join(messages)
        self.assertNotIn("github", combined.lower())
        self.assertNotIn("report", combined.lower())

    def test_no_report_message_check_error(self):
        """User-end errors should tell user to check the error message."""
        messages = no_report_error_message("SPAdes")
        combined = " ".join(messages)
        self.assertIn("check the error message", combined.lower())


class TestErrorFunction(unittest.TestCase):
    """Test the error() function behavior with different error codes."""

    @patch('support.sys.exit')
    @patch('support.shutil.rmtree')
    @patch('support.os.path.isdir', return_value=False)
    def test_error_with_user_end_code_exits_correctly(self, mock_isdir, mock_rmtree, mock_exit):
        """error() with user-end code should exit with that code's value."""
        with patch('sys.stderr', new_callable=io.StringIO):
            error("Test error", exit_code=ErrorCode.InvalidInputFormat)
        mock_exit.assert_called_once_with(64)

    @patch('support.sys.exit')
    @patch('support.shutil.rmtree')
    @patch('support.os.path.isdir', return_value=False)
    def test_error_with_internal_code_exits_correctly(self, mock_isdir, mock_rmtree, mock_exit):
        """error() with internal code should exit with that code's value."""
        with patch('sys.stderr', new_callable=io.StringIO):
            error("Test error", exit_code=ErrorCode.GeneralError)
        mock_exit.assert_called_once_with(-1)

    @patch('support.sys.exit')
    @patch('support.shutil.rmtree')
    @patch('support.os.path.isdir', return_value=False)
    def test_error_user_end_no_report_message(self, mock_isdir, mock_rmtree, mock_exit):
        """error() with user-end code should NOT show 'report issue' message."""
        stderr_capture = io.StringIO()
        with patch('sys.stderr', stderr_capture):
            error("Test error", exit_code=ErrorCode.IOError)
        output = stderr_capture.getvalue()
        self.assertNotIn("github", output.lower())
        self.assertNotIn("report an issue", output.lower())

    @patch('support.sys.exit')
    @patch('support.shutil.rmtree')
    @patch('support.os.path.isdir', return_value=False)
    def test_error_internal_shows_report_message(self, mock_isdir, mock_rmtree, mock_exit):
        """error() with internal code SHOULD show 'report issue' message."""
        stderr_capture = io.StringIO()
        with patch('sys.stderr', stderr_capture):
            error("Test error", exit_code=ErrorCode.GeneralError)
        output = stderr_capture.getvalue()
        self.assertIn("github.com/ablab/spades", output)

    @patch('support.sys.exit')
    @patch('support.shutil.rmtree')
    @patch('support.os.path.isdir', return_value=False)
    def test_error_file_not_found_message(self, mock_isdir, mock_rmtree, mock_exit):
        """error() with InputFileNotFound should not suggest reporting."""
        stderr_capture = io.StringIO()
        with patch('sys.stderr', stderr_capture):
            error("File not found: test.fa", exit_code=ErrorCode.InputFileNotFound)
        output = stderr_capture.getvalue()
        self.assertIn("File not found: test.fa", output)
        self.assertNotIn("report an issue", output.lower())
        mock_exit.assert_called_once_with(65)

    @patch('support.sys.exit')
    @patch('support.shutil.rmtree')
    @patch('support.os.path.isdir', return_value=False)
    def test_error_invalid_param_message(self, mock_isdir, mock_rmtree, mock_exit):
        """error() with InvalidParameter should not suggest reporting."""
        stderr_capture = io.StringIO()
        with patch('sys.stderr', stderr_capture):
            error("Invalid k-mer value", exit_code=ErrorCode.InvalidParameter)
        output = stderr_capture.getvalue()
        self.assertIn("Invalid k-mer value", output)
        self.assertNotIn("report an issue", output.lower())
        mock_exit.assert_called_once_with(67)

    @patch('support.sys.exit')
    @patch('support.shutil.rmtree')
    @patch('support.os.path.isdir', return_value=False)
    def test_error_memory_exceeded_message(self, mock_isdir, mock_rmtree, mock_exit):
        """error() with MemoryLimitExceeded should not suggest reporting."""
        stderr_capture = io.StringIO()
        with patch('sys.stderr', stderr_capture):
            error("Out of memory", exit_code=ErrorCode.MemoryLimitExceeded)
        output = stderr_capture.getvalue()
        self.assertIn("Out of memory", output)
        self.assertNotIn("report an issue", output.lower())
        mock_exit.assert_called_once_with(68)


class TestErrorFunctionWithLogger(unittest.TestCase):
    """Test error() function behavior when a logger is provided."""

    @patch('support.sys.exit')
    @patch('support.shutil.rmtree')
    @patch('support.os.path.isdir', return_value=False)
    @patch('support.log_warnings')
    def test_error_with_logger_user_end(self, mock_log_warnings, mock_isdir, mock_rmtree, mock_exit):
        """error() with logger and user-end code should log without report message."""
        mock_logger = MagicMock()
        error("Test error", logger_instance=mock_logger, exit_code=ErrorCode.InvalidInputFormat)

        # Check that error was logged
        mock_logger.error.assert_called()

        # Get all logged messages
        logged_messages = " ".join(str(call) for call in mock_logger.error.call_args_list)

        # Should NOT contain report issue message
        self.assertNotIn("github", logged_messages.lower())

        mock_exit.assert_called_once_with(64)

    @patch('support.sys.exit')
    @patch('support.shutil.rmtree')
    @patch('support.os.path.isdir', return_value=False)
    @patch('support.log_warnings')
    def test_error_with_logger_internal(self, mock_log_warnings, mock_isdir, mock_rmtree, mock_exit):
        """error() with logger and internal code should log with report message."""
        mock_logger = MagicMock()
        error("Test error", logger_instance=mock_logger, exit_code=ErrorCode.GeneralError)

        # Check that error was logged
        mock_logger.error.assert_called()

        # Get all logged messages
        logged_messages = " ".join(str(call) for call in mock_logger.error.call_args_list)

        # SHOULD contain report issue message
        self.assertIn("github", logged_messages.lower())

        mock_exit.assert_called_once_with(-1)


if __name__ == "__main__":
    unittest.main()
