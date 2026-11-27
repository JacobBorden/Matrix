"""Qualify the pinned constructor syntax failure without depending on one GCC message."""
import re
from pathlib import Path


def is_constructor_rejection(exit_code, output, header):
    if exit_code == 0 or re.search(r'fatal error:|internal compiler error:|Killed signal', output):
        return False
    header = Path(header).resolve()
    constructor_lines = {
        number for number, line in enumerate(header.read_text().splitlines(), 1)
        if re.match(r'^\s*(?:explicit\s+)?Matrix<T>\(', line)
    }
    errors = re.findall(r'^(.+?):(\d+):\d+: error: (.+)$', output, re.MULTILINE)
    # The probe requests one compiler error, so later cascades cannot obscure its cause.
    if len(errors) != 1:
        return False
    path, line, message = errors[0]
    return (Path(path).resolve() == header and int(line) in constructor_lines
            and bool(re.search(r'template-id|constructor|expected .*before', message)))
