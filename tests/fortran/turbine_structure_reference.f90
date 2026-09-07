! Development reference only. Links the unmodified pinned ElastoDyn module.
! This executable is never called by the C++ turbine solver.
program structure_reference
use NWTC_Library
use ElastoDyn
use ElastoDyn_Types
use ElastoDyn_Parameters, only: DOF_BF, DOF_BE, DOF_GeAz
implicit none
type(ED_InitInputType) :: init
type(ED_InitOutputType) :: initout
type(ED_InputType) :: u
type(ED_ParameterType) :: p
type(ED_ContinuousStateType) :: x, dx
type(ED_DiscreteStateType) :: xd
type(ED_ConstraintStateType) :: z
type(ED_OtherStateType) :: other
type(ED_OutputType) :: y
type(ED_MiscVarType) :: m
character(2048) :: filename,output,msg
integer :: err,n,b,j,k,outunit,dofs(3)
real(DbKi) :: dt,t,initial_azimuth
real(R8Ki) :: q(3),qd(3),v(12)
call get_command_argument(1,filename)
call get_command_argument(2,output)
init%InputFile=filename
init%RootName=trim(output)//'.ed'
init%CompElast=.true.
init%Gravity=9.80665_ReKi
dt=.00625_DbKi
call ED_Init(init,u,p,x,xd,z,other,y,m,dt,initout,err,msg)
call check()
initial_azimuth=x%QT(DOF_GeAz)
open(newunit=outunit,file=trim(output),status='replace')
write(outunit,'(A)') 'sample,blade,node,ax1,ax2,ax3,x,y,z,vx,vy,vz,nx,ny,nz'
do n=0,99
  t=.037_DbKi*n
  x%QT(DOF_GeAz)=initial_azimuth+p%RotSpeed*t
  do b=1,3
    q=(/.8_R8Ki*sin(t+b),.1_R8Ki*cos(2*t+b),.4_R8Ki*sin(.7_R8Ki*t-b)/)
    qd=(/.8_R8Ki*cos(t+b),-.2_R8Ki*sin(2*t+b),.28_R8Ki*cos(.7_R8Ki*t-b)/)
    dofs=(/DOF_BF(b,1),DOF_BF(b,2),DOF_BE(b,1)/)
    x%QT(dofs)=q
    x%QDT(dofs)=qd
    do j=1,p%BldNodes
      u%BladePtLoads(b)%Force(:,j)=p%DRNodes(j)*(/real(1000+10*j,R8Ki),200*sin(t+j),20*cos(t+j)/)
      u%BladePtLoads(b)%Moment(:,j)=p%DRNodes(j)*(/10*sin(t),5*cos(t),2.0_R8Ki/)
    end do
  end do
  call ED_CalcContStateDeriv(t,u,p,x,xd,z,other,m,dx,err,msg)
  call check()
  do b=1,3
    dofs=(/DOF_BF(b,1),DOF_BF(b,2),DOF_BE(b,1)/)
    do j=0,p%TipNode
      v(1:3)=dx%QDT(dofs)
      v(4:6)=(/m%RtHS%rS(1,b,j),-m%RtHS%rS(3,b,j),m%RtHS%rS(2,b,j)/)
      v(7:9)=(/m%RtHS%LinVelES(1,j,b),-m%RtHS%LinVelES(3,j,b),m%RtHS%LinVelES(2,j,b)/)
      v(10:12)=(/m%CoordSys%n3(b,j,1),-m%CoordSys%n3(b,j,3),m%CoordSys%n3(b,j,2)/)
      write(outunit,'(I0,A,I0,A,I0,12(A,ES25.16E3))') n,',',b-1,',',j,(',',v(k),k=1,12)
    end do
  end do
end do
close(outunit)
contains
subroutine check()
if(err>=ErrID_Fatal)then
  print *,trim(msg)
  stop 1
end if
end subroutine
end program
