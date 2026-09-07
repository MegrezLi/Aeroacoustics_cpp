program bem_reference
use NWTC_Library
use AirfoilInfo
use BEMT_Types
use BEMTUnCoupled
implicit none
type(AFI_InitInputType) :: init
type(AFI_ParameterType),allocatable :: af(:)
type(BEMT_ParameterType) :: p
type(BEMT_InputType) :: u
character(2048) :: filename,output,msg
integer :: fi,fo,err,na,nt,n,i,ad,td,ti
logical :: valid
real(ReKi) :: phi,r,a,ap,k,kp,f
call NWTC_Init(EchoLibVer=.false.)
call get_command_argument(1,filename)
call get_command_argument(2,output)
open(newunit=fi,file=trim(filename),status='old')
open(newunit=fo,file=trim(output),status='replace')
read(fi,*) na,nt
allocate(af(na))
init%AFTabMod=1
init%InCol_Alfa=1
init%InCol_Cl=2
init%InCol_Cd=3
init%InCol_Cm=4
init%UAMod=3
do i=1,na
  read(fi,'(A)') filename
  init%FileName=trim(filename)
  call AFI_Init(init,af(i),err,msg)
  if(err>=ErrID_Fatal)then
    print *,trim(msg)
    stop 1
  endif
enddo
p%numBlades=3
p%BEM_Mod=1
p%kinVisc=1.81206e-5_ReKi
p%useTipLoss=.true.
p%useHubLoss=.true.
allocate(p%chord(1,1),p%hubLossConst(1,1),p%tipLossConst(1,1),p%FixedInductions(1,1))
p%FixedInductions=.false.
allocate(u%rLocal(1,1),u%Vx(1,1),u%Vy(1,1),u%Vz(1,1),u%theta(1,1),u%cantAngle(1,1),u%toeAngle(1,1),u%UserProp(1,1))
u%Vz=0
u%cantAngle=0
u%toeAngle=0
u%UserProp=0
write(fo,'(A)') 'sample,residual,a,ap,k,kp,F,valid'
do n=0,nt-1
  read(fi,*) i,ad,td,ti,u%rLocal(1,1),p%chord(1,1),u%theta(1,1),u%Vx(1,1),u%Vy(1,1),p%hubLossConst(1,1),p%tipLossConst(1,1),phi
  p%useAIDrag=ad==1
  p%useTIDrag=td==1
  p%useTanInd=ti==1
  r=BEMTU_InductionWithResidual(p,u,1,1,phi,af(i),valid,err,msg,a,ap,k,kp,f)
  if(err>=ErrID_Fatal)then
    print *,trim(msg)
    stop 1
  endif
  write(fo,'(I0,6(A,ES25.16E3),A,I0)') n,',',r,',',a,',',ap,',',k,',',kp,',',f,',',merge(1,0,valid)
enddo
close(fi)
close(fo)
end program
