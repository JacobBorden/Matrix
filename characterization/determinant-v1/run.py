#!/usr/bin/env python3
"""Characterize immutable legacy determinants in separate executables."""
import argparse,hashlib,json,os,re,subprocess,sys
from pathlib import Path
root=Path(__file__).resolve().parent
legacy=root.parent/'multiplication-v1'
p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);args=p.parse_args()
out=args.output.resolve();out.mkdir(parents=True,exist_ok=True)
sources=json.loads((legacy/'sources.json').read_text())
for s in sources:
    owner='neuro' if s['repository']=='Neuro-Net-DLL' else 'bitmap'
    file=legacy/'legacy'/owner/s['path']
    if hashlib.sha256(file.read_bytes()).hexdigest()!=s['sha256']:sys.exit('Source hash mismatch: '+str(file))
results=[]
def run(label,cmd,tests=False):
    proc=subprocess.run(cmd,cwd=root.parents[1],capture_output=True,text=True,timeout=90,env={**os.environ,'OMP_NUM_THREADS':'2'})
    log=out/(label+'.log');log.write_text(proc.stdout+proc.stderr)
    entry=dict(label=label,command=cmd,exit_code=proc.returncode,log=str(log))
    if tests:
        m=re.search(r'RESULT passed=(\d+) failed=(\d+)',proc.stdout)
        entry.update(passed=int(m[1]) if m else 0,failed=int(m[2]) if m else None)
    results.append(entry);print(label,entry['exit_code'],entry.get('passed',''))
    return proc.returncode
run('compiler',['g++','--version']);run('platform',['uname','-a'])
for owner,std,header in [('neuro','c++17','src/math/matrix.h'),('bitmap','c++20','src/matrix/matrix.h')]:
    for mode,flags in [('debug',['-O0','-g']),('release',['-O2','-DNDEBUG']),('sanitized',['-O1','-g','-fsanitize=address,undefined','-fno-omit-frame-pointer'])]:
        exe=out/(owner+'-'+mode)
        cmd=['g++','-std='+std,'-Wall','-Wextra','-Wno-unknown-pragmas','-ffp-contract=off']+flags
        if owner=='neuro':cmd+=['-DNEURO=1']
        cmd+=['-DLEGACY_HEADER="'+str(legacy/'legacy'/owner/header)+'"',str(root/'determinant.cpp'),'-o',str(exe)]
        if run(owner+'-'+mode+'-build',cmd)==0:run(owner+'-'+mode+'-tests',[str(exe)],True)
summary=dict(sources=sources,results=results,expected_groups_per_execution=44,expected_executions=6)
(out/'results.json').write_text(json.dumps(summary,indent=2)+'\n')
bad=any(r['exit_code']!=0 or ('passed' in r and (r['passed']!=44 or r['failed']!=0)) for r in results)
if sum('passed' in r for r in results)!=6:bad=True
sys.exit(1 if bad else 0)
