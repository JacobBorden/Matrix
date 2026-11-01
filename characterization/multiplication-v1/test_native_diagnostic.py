"""Compiler-wording variations must not turn unrelated build failures into controls."""
from pathlib import Path
import unittest
from native_diagnostic import is_constructor_rejection


HEADER = Path(__file__).resolve().parent / 'legacy/neuro/src/math/matrix.h'


class ConstructorDiagnostic(unittest.TestCase):
    def diagnostic(self, message, line=422, path=HEADER):
        return f'{path}:{line}:26: error: {message}\n'

    def test_explicit_gcc_constructor_message(self):
        self.assertTrue(is_constructor_rejection(1, self.diagnostic(
            'template-id not allowed for constructor in C++20 [-Wtemplate-id-cdtor]'), HEADER))

    def test_parser_wording_variation(self):
        # Representative parser wording, not a claim that a particular CI emitted it.
        self.assertTrue(is_constructor_rejection(1, self.diagnostic(
            "expected unqualified-id before ')' token"), HEADER))
        self.assertTrue(is_constructor_rejection(1, self.diagnostic(
            "expected ')' before 'row_count'", line=440), HEADER))

    def test_success_is_not_a_negative_control(self):
        self.assertFalse(is_constructor_rejection(0, self.diagnostic('constructor rejected'), HEADER))

    def test_other_file_is_rejected(self):
        self.assertFalse(is_constructor_rejection(1, self.diagnostic(
            'constructor rejected', path=Path('unrelated.cpp')), HEADER))

    def test_other_header_location_is_rejected(self):
        self.assertFalse(is_constructor_rejection(1, self.diagnostic(
            "expected ';' before 'x'", line=10), HEADER))

    def test_missing_header_is_rejected(self):
        self.assertFalse(is_constructor_rejection(1,
            f'{HEADER}:10:1: fatal error: missing.h: No such file or directory\n', HEADER))

    def test_killed_compiler_is_rejected(self):
        self.assertFalse(is_constructor_rejection(1, 'g++: fatal error: Killed signal terminated program cc1plus\n', HEADER))

    def test_unrelated_constructor_line_error_is_rejected(self):
        self.assertFalse(is_constructor_rejection(1, self.diagnostic('division by zero'), HEADER))

    def test_mixed_failure_is_rejected(self):
        output = self.diagnostic('constructor rejected') + self.diagnostic(
            "expected ';' before 'x'", path=Path('unrelated.cpp'))
        self.assertFalse(is_constructor_rejection(1, output, HEADER))


if __name__ == '__main__':
    unittest.main(verbosity=2)
