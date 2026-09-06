"""Replace seven OpenFAST acoustic kernels with ISO_C_BINDING calls to C++.

Generate the reviewable patch first; --apply modifies only the pinned source.
The original OpenFAST outer time loop, states, geometry and file I/O are retained.
"""
import argparse
import difflib
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
COMMON = {1:'ALPSTAR',2:'C',3:'U',4:'THETA',5:'PHI',6:'L',7:'R',
          8:'StallVal',9:'d99Var2',10:'dstarVar1',11:'dstarVar2'}
SPECS = [
    ('LBLVS',1,COMMON,['SPLLAM']),
    ('TBLTE',2,COMMON,['SPLP','SPLS','SPLALPH']),
    ('TIPNOIS',3,{1:'ALPHTIP',2:'C',3:'U',4:'THETA',5:'PHI',7:'R',21:'ALPRAT2'},['SPLTIP']),
    ('InflowNoise',4,{1:'AlphaNoise',2:'Chord',3:'U',4:'THETA',5:'PHI',6:'d',7:'RObs',12:'TINoise'},['SPLti']),
    ('BLUNT',5,{**COMMON,13:'H',14:'PSI'},['SPLBLUNT']),
    ('Simple_Guidati',6,{2:'Chord',3:'U',15:'thick_1p',16:'thick_10p'},['SPLti']),
    ('TBLTE_TNO',7,{3:'U',4:'THETA',5:'PHI',6:'D',7:'R',22:'Cfall(1)',23:'Cfall(2)',
                   24:'d99all(1)',9:'d99all(2)',25:'EdgeVelAll(1)',26:'EdgeVelAll(2)'},['SPLP','SPLS'])]

INTERFACE = '''
   ! Aeroacoustics_cpp: no Python or Fortran numerical kernel at runtime.
   interface
      function cpp_acoustic_kernel(op, n, freq, x, flags, output) &
         bind(C, name="aeroacoustics_kernel") result(status)
         import c_int, c_double
         integer(c_int), value :: op, n
         real(c_double), intent(in) :: freq(*), x(*)
         integer(c_int), intent(in) :: flags(*)
         real(c_double), intent(out) :: output(*)
         integer(c_int) :: status
      end function
      function cpp_acoustic_calls() bind(C, name="aeroacoustics_call_count") result(n)
         import c_long_long
         integer(c_long_long) :: n
      end function
   end interface
'''

def transform(original):
    source=original.replace('module AeroAcoustics\n',
        'module AeroAcoustics\n   use iso_c_binding, only: c_int, c_double, c_long_long\n',1)
    source=source.replace('   contains',INTERFACE+'\n   contains',1)
    for name,op,fields,outputs in SPECS:
        pattern=rf'^SUBROUTINE {name}\(.*?^END SUBROUTINE {name}\b[^\n]*'
        match=re.search(pattern,source,re.M|re.S|re.I)
        if not match:raise ValueError(f'Missing routine {name}')
        old=match.group()
        lines=[old.splitlines()[0]]+[line for line in old.splitlines()
            if not line.lstrip().startswith('!') and re.search(r'\bINTENT\s*\(',line,re.I)]
        lines+=['    real(c_double) :: args(26), cpp_freq(size(p%FreqList))',
                '    real(c_double) :: result(3*size(p%FreqList))',
                '    integer(c_int) :: flags(3), status, nfreq',
                '    args=0.0_c_double',
                '    args(17)=real(p%SpdSound,c_double)',
                '    args(18)=real(p%KinVisc,c_double)',
                '    args(19)=real(p%AirDens,c_double)',
                '    args(20)=real(p%LTurb,c_double)',
                '    flags=[int(p%X_BLMethod,c_int),int(p%ITRIP,c_int),merge(1_c_int,0_c_int,p%ROUND)]']
        lines += [f'    args({i})=real({value},c_double)' for i,value in fields.items()]
        lines += ['    nfreq=int(size(p%FreqList),c_int)',
                  '    cpp_freq=real(p%FreqList,c_double)',
                  f'    status=cpp_acoustic_kernel({op}_c_int,nfreq,cpp_freq,args,flags,result)',
                  f'    if (status /= 0) error stop "Aeroacoustics_cpp kernel failed: {name}"']
        lines += [f'    {value}=real(result({i}*nfreq+1:{i+1}*nfreq),ReKi)' for i,value in enumerate(outputs)]
        lines += [f'END SUBROUTINE {name}']
        source=source[:match.start()]+'\n'.join(lines)+source[match.end():]
    end=re.search(r'^subroutine AA_End\(.*?^END SUBROUTINE AA_End',source,re.M|re.S|re.I)
    body=end.group().replace('    ErrStat = ErrID_None',
        '    write(*,\'(A,I0)\') "C++ acoustic kernel calls: ",cpp_acoustic_calls()\n    ErrStat = ErrID_None',1)
    return source[:end.start()]+body+source[end.end():]

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,default=ROOT/'_work/openfast')
    parser.add_argument('--apply',action='store_true')
    parser.add_argument('--restore',action='store_true',help='Restore the verified original two source files')
    args=parser.parse_args()
    relative=Path('modules/aerodyn/src/AeroAcoustics.f90')
    original=(ROOT/'reference/fortran/AeroAcoustics.f90').read_text(encoding='utf-8')
    updated=transform(original)
    cmake_relative=Path('modules/aerodyn/CMakeLists.txt')
    cmake_path=args.source/cmake_relative
    cmake=cmake_path.read_text(encoding='utf-8')
    cmake_backup=cmake_path.with_suffix('.txt.original')
    if cmake_backup.exists():cmake=cmake_backup.read_text(encoding='utf-8')
    marker='target_link_libraries(aeroacousticslib basicaerolib nwtclibs)'
    assert marker in cmake
    cmake_updated=cmake.replace(marker,'''if(NOT AEROACOUSTICS_CPP_SOURCE)
  message(FATAL_ERROR "Set AEROACOUSTICS_CPP_SOURCE to the Aeroacoustics_cpp checkout")
endif()
add_subdirectory("${AEROACOUSTICS_CPP_SOURCE}" "${CMAKE_BINARY_DIR}/acoustics-cpp")
install(TARGETS aeroacoustics EXPORT OpenFASTLibraries ARCHIVE DESTINATION lib)
install(DIRECTORY "${AEROACOUSTICS_CPP_SOURCE}/include/" DESTINATION include/aeroacoustics)
target_link_libraries(aeroacousticslib basicaerolib nwtclibs aeroacoustics)''',1)
    patch=[]
    for rel,before,after in [(relative,original,updated),(cmake_relative,cmake,cmake_updated)]:
        path=args.source/rel
        actual=path.read_text(encoding='utf-8')
        if actual not in (before,after):raise RuntimeError(f'Refusing to overwrite different source: {path}')
        patch.extend(difflib.unified_diff(before.splitlines(True),after.splitlines(True),
                     fromfile='a/'+rel.as_posix(),tofile='b/'+rel.as_posix()))
        if args.restore:
            path.write_text(before,encoding='utf-8')
        elif args.apply and actual!=after:
            path.with_suffix(path.suffix+'.original').write_text(before,encoding='utf-8')
            path.write_text(after,encoding='utf-8')
    dest=ROOT/'integration/openfast-cpp.patch'
    dest.parent.mkdir(exist_ok=True)
    dest.write_text(''.join(patch),encoding='utf-8')
    print(f'{"Restored originals" if args.restore else "Applied" if args.apply else "Generated"}: {dest}')

if __name__=='__main__':main()
