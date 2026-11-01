#!/usr/bin/env python3
"""Reproducible, separate legacy builds; verify source pins before compiling."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
from native_diagnostic import is_constructor_rejection

ROOT = Path(__file__).resolve().parent
parser = argparse.ArgumentParser()
parser.add_argument('--output', type=Path, required=True)
args = parser.parse_args()
out = args.output.resolve()
out.mkdir(parents=True, exist_ok=True)
records = json.loads((ROOT/'sources.json').read_text())
for record in records:
    owner = 'neuro' if record['repository'] == 'Neuro-Net-DLL' else 'bitmap'
    path = ROOT/'legacy'/owner/record['path']
    if hashlib.sha256(path.read_bytes()).hexdigest() != record['sha256']:
        sys.exit('Source hash mismatch: '+str(path))
results = []
def run(label, command, runtime=False):
    proc = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=90,
                          env={**os.environ, 'OMP_NUM_THREADS': '2'})
    log = out/(label+'.log'); log.write_text(proc.stdout)
    entry = {'label':label,'command':command,'exit_code':proc.returncode,'log':str(log),'environment':{'OMP_NUM_THREADS':'2'}}
    if runtime:
        m = re.search(r'RESULT passed=(\d+) failed=(\d+)',proc.stdout)
        entry.update(passed=int(m[1]) if m else 0, failed=int(m[2]) if m else None)
    results.append(entry)
    print(label, 'exit='+str(proc.returncode), 'passed='+str(entry.get('passed','n/a')))
    return proc.returncode
run('compiler',['g++','--version'])
run('platform',['uname','-a'])
# Strict C++20 rejects Neuro's template-id constructor syntax; GCC otherwise accepts an extension.
# Keep this native compile probe explicit; do not repair or substitute the source.
headers = {'neuro':'legacy/neuro/src/math/matrix.h','bitmap':'legacy/bitmap/src/matrix/matrix.h'}
base = ['g++','-Wall','-Wextra','-Wno-unknown-pragmas','-ffp-contract=off']
for owner,header in headers.items():
    define = '-DLEGACY_HEADER="'+str(ROOT/header)+'"'
    native = out/(owner+'-cxx20')
    native_flags = ['-fmax-errors=1'] if owner == 'neuro' else []
    run(owner+'-native-cxx20',base+['-std=c++20','-pedantic-errors','-O0']+native_flags+[define,str(ROOT/'multiplication.cpp'),'-o',str(native)])
# Successful baseline builds: pinned Neuro C++17; pinned Bitmap C++20.
for owner,standard in [('neuro','c++17'),('bitmap','c++20')]:
    define = '-DLEGACY_HEADER="'+str(ROOT/headers[owner])+'"'
    for mode,flags in [('debug',['-O0','-g']),('release',['-O2']),('sanitized',['-O1','-g','-fsanitize=address,undefined','-fno-omit-frame-pointer'])]:
        exe = out/(owner+'-'+mode)
        if run(owner+'-'+mode+'-build',base+['-std='+standard]+flags+[define,str(ROOT/'multiplication.cpp'),'-o',str(exe)]) == 0:
            run(owner+'-'+mode+'-tests',[str(exe)],True)
    if owner == 'neuro':
        exe = out/'neuro-openmp'
        if run('neuro-openmp-build',base+['-std=c++17','-O2','-fopenmp',define,str(ROOT/'multiplication.cpp'),'-o',str(exe)]) == 0:
            run('neuro-openmp-tests',[str(exe)],True)
summary = {'sources':records,'results':results,'test_groups_per_execution':16,'expected_native_failure':'neuro-native-cxx20'}
(out/'results.json').write_text(json.dumps(summary,indent=2)+'\n')
expected_failure = [r for r in results if r['label']=='neuro-native-cxx20']
bad = [r for r in results if r['label'] != 'neuro-native-cxx20' and (r['exit_code'] != 0 or ('passed' in r and (r['passed'] != 16 or r['failed'] != 0)))]
if not expected_failure or not is_constructor_rejection(expected_failure[0]['exit_code'], Path(expected_failure[0]['log']).read_text(), ROOT/headers['neuro']):
    sys.exit('Native compatibility outcome changed; review expectations')
sys.exit(1 if bad else 0)
