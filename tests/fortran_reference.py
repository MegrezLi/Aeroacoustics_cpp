"""Untouched upstream numerical bodies plus minimal standalone Fortran types shim."""
from pathlib import Path
import ctypes as ct
import os
import re
import subprocess
ROOT=Path(__file__).resolve().parents[1]
def extract(source,name):
    pat=rf'(?im)^\s*(?:real\(ReKi\)\s+function|subroutine)\s+{name}\(.*?^\s*end (?:function|subroutine)\s+{name}\b[^\n]*'
    return re.search(pat,source,re.S).group(0)

NAMES=['Log10AA','LBLVS','TBLTE','TIPNOIS','InflowNoise','BLUNT','G5COMP',
       'AMIN','AMAX','BMIN','BMAX','A0COMP','THICK','DIRECTH_TE','DIRECTH_LE',
       'DIRECTL','Simple_Guidati','TBLTE_TNO']

def build(fc,precision):
    folder=ROOT/'_work/fortran-reference'/precision
    folder.mkdir(parents=True,exist_ok=True)
    main=(ROOT/'reference/fortran/AeroAcoustics.f90').read_text(encoding='utf-8')
    kind=8 if precision=='double' else 4
    prefix='d' if kind==8 else ''
    library=f'''module NWTC_Library
implicit none
integer, parameter :: ReKi={kind},DbKi=8,R8Ki=8,IntKi=4
real(ReKi),parameter :: pi=acos(-1._ReKi),TwoPi=2*pi
end module
module NWTC_SLATEC
use NWTC_Library
implicit none
contains
subroutine slatec_qk61(f,a,b,answer,abserr,resabs,resasc)
real(ReKi),external :: f
real(ReKi) :: a,b,answer,abserr,resabs,resasc
call {prefix}qk61(f,a,b,answer,abserr,resabs,resasc)
end subroutine
end module
'''
    library+=(ROOT/'reference/fortran/AeroAcoustics_TNO.f90').read_text(encoding='utf-8')
    library+='''
module acoustic_reference
use NWTC_Library
implicit none
real(ReKi),parameter :: AA_EPSILON=1.e-16
integer,parameter :: X_BLMethod_Tables=2,ITRIP_None=0,ITRIP_Heavy=1,ITRIP_Light=2
type AA_ParameterType
real(ReKi) :: SpdSound,KinVisc,AirDens,Lturb
real(ReKi),allocatable :: FreqList(:)
integer :: X_BLMethod,ITRIP
logical :: ROUND
end type
contains
'''
    library+='\n'.join(extract(main,n) for n in NAMES)
    library+='''
subroutine evaluate(x,flags,freq,n,y)
integer :: n,flags(4)
real(ReKi) :: x(22),freq(n),y(n,10)
type(AA_ParameterType) :: p
p%SpdSound=x(17);p%KinVisc=x(18);p%AirDens=x(19);p%Lturb=x(20)
p%X_BLMethod=flags(1);p%ITRIP=flags(2);p%ROUND=flags(3)==1
allocate(p%FreqList(n));p%FreqList=freq
call LBLVS(x(1),x(2),x(3),x(4),x(5),x(6),x(7),p,x(9),x(10),x(11),y(:,1),x(8))
call TBLTE(x(1),x(2),x(3),x(4),x(5),x(6),x(7),p,x(9),x(10),x(11),x(8),y(:,2),y(:,3),y(:,4))
call TIPNOIS(x(1),x(21),x(2),x(3),x(4),x(5),x(7),p,y(:,5))
call InflowNoise(x(1)*pi/180,x(2),x(3),x(4),x(5),x(6),x(7),x(12),p,y(:,6))
call BLUNT(x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(13),x(14),p,x(9),x(10),x(11),y(:,7),x(8))
call Simple_Guidati(x(3),x(2),x(16),x(15),p,y(:,8))
y(:,9:10)=0
if(flags(4)==1) call TBLTE_TNO(x(3),x(4),x(5),x(6),x(7), &
  (/x(22),0.001984380_ReKi/),(/0.01105860_ReKi,x(9)/), &
  (/1.0_ReKi,1.0_ReKi/),p,y(:,9),y(:,10))
end subroutine
subroutine curves(x,y)
real(ReKi) :: x(5),y(9)
type(AA_ParameterType) :: p
y(1)=G5COMP(x(1),x(2));y(2)=AMIN(x(2));y(3)=AMAX(x(2))
y(4)=BMIN(x(2));y(5)=BMAX(x(2));y(6)=A0COMP(x(3))
y(7)=DIRECTH_TE(x(4),x(5),90._ReKi)
y(8)=DIRECTH_LE(x(4),x(5),90._ReKi)
y(9)=DIRECTL(x(4),x(5),90._ReKi)
end subroutine
end module
'''
    (folder/'reference_shim.f90').write_text(library,encoding='utf-8')
    # SLATEC asks only for machine epsilon and minimum normal (indices 4 and 1).
    mname='d1mach' if kind==8 else 'r1mach'
    machine=f'''function {mname}(i) result(v)
integer :: i
real(kind={kind}) :: v
select case(i)
case(1)
v=tiny(v)
case(4)
v=epsilon(v)
case default
stop 'unexpected machine constant'
end select
end function
'''
    (folder/'machine.f90').write_text(machine)
    dll=folder/('reference.dll' if os.name=='nt' else 'reference.so')
    cmd=[fc,'-shared','-O0','-fcheck=all','-frecursive','-ffree-line-length-none']
    if kind==8: cmd+=['-fdefault-real-8','-fdefault-double-8']
    if os.name!='nt':cmd+=['-fPIC']
    cmd+=['reference_shim.f90','machine.f90',str(ROOT/f'reference/fortran/{prefix}qk61.f'),'-o',str(dll)]
    subprocess.run(cmd,cwd=folder,check=True,capture_output=True,text=True)
    return ct.CDLL(str(dll)),cmd