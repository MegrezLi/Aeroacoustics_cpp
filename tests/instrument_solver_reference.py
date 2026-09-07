"""Add read-only diagnostic output to a pinned OpenFAST development checkout."""
from pathlib import Path
import argparse

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('source',type=Path)
args=p.parse_args()
path=args.source/'modules/openfast-library/src/FAST_Solver.f90'
text=path.read_text()
if 'subroutine TraceTurbine(' not in text:
    backup=path.with_suffix('.f90.trace-original')
    if backup.exists():raise RuntimeError('Existing trace backup; inspect before changing the source')
    backup.write_text(text)
    start=text.index('subroutine FAST_SolverStep(n_t_global')
    end=text.index('\ncontains',start)
    text=text[:end]+'\n   call TraceTurbine(t_global_next, Turbine)\n'+text[end:]
    index=text.lower().rindex('end module')
    text=text[:index]+Path(__file__).with_name('fortran').joinpath('trace_full_solver.inc').read_text()+'\n'+text[index:]
    path.write_text(text)
print(path)
