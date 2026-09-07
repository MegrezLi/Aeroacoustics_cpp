"""Development-only differential check against the pinned original Fortran modules.

This compiles a reference test harness, not a dependency of the C++ executable.
It exercises all 30 official polars under prescribed oscillating AoA inputs.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import numpy as np

HARNESS = r'''
program reference_ua
use NWTC_Library
use AirfoilInfo
use AirfoilInfo_Types
use UnsteadyAero
use UnsteadyAero_Types
implicit none
type(AFI_InitInputType) :: ai
type(AFI_ParameterType) :: af(1)
type(UA_InitInputType) :: init
type(UA_InitOutputType) :: initout
type(UA_InputType) :: u(1)
type(UA_ParameterType) :: p
type(UA_ContinuousStateType) :: x
type(UA_DiscreteStateType) :: xd
type(UA_OtherStateType) :: other
type(UA_OutputType) :: y
type(UA_MiscVarType) :: m
integer :: i,n,err,afid,idx(1,1),outunit,bladeunit
real(ReKi) :: span(30),curve,sweep,ca,twist,chord(30),t,alpha,speed
real(DbKi) :: times(1),dt
character(2048) :: root,output,filename,line,msg
character(2) :: number
call get_command_argument(1,root)
call get_command_argument(2,output)
call NWTC_Init(EchoLibVer=.false.)
open(newunit=bladeunit,file=trim(root)//'/RotorSE_FAST_IEA_landBased_RWT_AeroDyn_blade.dat',status='old')
do i=1,6
  read(bladeunit,'(A)') line
end do
do i=1,30
  read(bladeunit,*) span(i),curve,sweep,ca,twist,chord(i),afid
end do
close(bladeunit)
open(newunit=outunit,file=trim(output),status='replace')
write(outunit,'(A)') 'foil,step,alpha_rad,speed,cl,cd,cm'
dt=0.00625_DbKi
do i=1,30
  write(number,'(I2.2)') i-1
  ai%FileName=trim(root)//'/Airfoils/RotorSE_FAST_IEA_landBased_RWT_AeroDyn_Polar_'//number//'.dat'
  ai%AFTabMod=1
  ai%InCol_Alfa=1
  ai%InCol_Cl=2
  ai%InCol_Cd=3
  ai%InCol_Cm=4
  ai%UAMod=3
  call AFI_Init(ai,af(1),err,msg)
  call check()
  init%dt=dt
  init%numBlades=1
  init%nNodesPerBlade=1
  init%UAMod=3
  init%a_s=335
  init%Flookup=.true.
  init%ShedEffect=.true.
  init%OutRootName='ua-reference'
  allocate(init%c(1,1),init%UAOff_innerNode(1),init%UAOff_outerNode(1))
  init%c=chord(i)
  init%UAOff_innerNode=0
  init%UAOff_outerNode=2
  idx=1
  call UA_Init(init,u(1),p,x,xd,other,y,m,dt,af,idx,initout,err,msg)
  call check()
  do n=0,1199
    t=n*dt
    alpha=(8+6*sin(2*Pi*.7_ReKi*t))*D2R
    speed=45+5*sin(.8_ReKi*t)
    u(1)%alpha=alpha
    u(1)%U=speed
    u(1)%Re=speed*chord(i)/1.81206e-5_ReKi
    u(1)%UserProp=0
    u(1)%omega=0
    times=t
    call UA_CalcOutput(1,1,real(t,DbKi),u(1),p,x,xd,other,af(1),y,m,err,msg)
    call check()
    write(outunit,'(I0,A,I0,5(A,ES25.16E3))') i-1,',',n,',',alpha,',',speed,',',y%Cl,',',y%Cd,',',y%Cm
    call UA_UpdateStates(1,1,real(t,DbKi),n,u,times,p,x,xd,other,af(1),m,err,msg)
    call check()
  end do
  call UA_DestroyInitInput(init,err,msg)
  call UA_DestroyParam(p,err,msg)
  call UA_DestroyContState(x,err,msg)
  call UA_DestroyDiscState(xd,err,msg)
  call UA_DestroyOtherState(other,err,msg)
  call UA_DestroyMisc(m,err,msg)
  call UA_DestroyOutput(y,err,msg)
  call UA_DestroyInitOutput(initout,err,msg)
  call AFI_DestroyParam(af(1),err,msg)
end do
close(outunit)
contains
subroutine check()
if(err>=ErrID_Fatal) then
  print *,trim(msg)
  stop 1
end if
end subroutine
end program
'''

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--openfast-build',type=Path,required=True)
    p.add_argument('--case',type=Path,required=True)
    p.add_argument('--cpp-csv',type=Path,required=True)
    p.add_argument('--work',type=Path,required=True)
    p.add_argument('--compiler',default='gfortran')
    p.add_argument('--mkl-lib',type=Path)
    p.add_argument('--link-library',action='append',default=[])
    args=p.parse_args()
    args.work.mkdir(parents=True,exist_ok=True)
    source=args.work/'reference_ua.f90';source.write_text(HARNESS)
    executable=args.work/'reference_ua.exe'
    b=args.openfast_build.resolve()
    libraries=[b/'modules/aerodyn/libaerodynlib.a',b/'modules/aerodyn/libbasicaerolib.a',b/'modules/nwtc-library/libnwtclibs.a']
    subprocess.run([args.compiler,'-ffree-line-length-none','-fdefault-real-8','-fdefault-double-8',f'-I{b / "ftnmods"}',str(source),*[str(x) for x in libraries],*([str(args.mkl_lib)] if args.mkl_lib else args.link_library or ['-llapack','-lblas']),'-lstdc++','-o',str(executable)],check=True)
    reference=args.work/'reference.csv'
    subprocess.run([str(executable.resolve()),str(args.case.resolve()),str(reference.resolve())],check=True)
    cpp=np.genfromtxt(args.cpp_csv,delimiter=',',names=True)
    baseline=np.genfromtxt(reference,delimiter=',',names=True)
    report={'cases':len(cpp),'reference':'original pinned OpenFAST AirfoilInfo and UnsteadyAero modules','max_absolute_error':{}}
    for name in ['foil','step','alpha_rad','speed','cl','cd','cm']:
        error=np.abs(cpp[name]-baseline[name]);worst=int(np.argmax(error))
        report['max_absolute_error'][name]={'value':float(error[worst]),'foil':int(cpp['foil'][worst]),'step':int(cpp['step'][worst])}
    (args.work/'comparison.json').write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2))
    if max(report['max_absolute_error'][x]['value'] for x in ['cl','cd','cm'])>1e-9:
        raise RuntimeError('Unsteady airfoil comparison failed')

if __name__=='__main__':main()
