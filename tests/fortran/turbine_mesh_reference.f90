! Original NWTC mesh mapping, used only as a development reference.
program mesh_reference
use NWTC_Library
implicit none
type(MeshType) :: sm,dm,dl,sl,sd
type(MeshMapType) :: motionmap,loadmap
integer :: err,fi,fo,ns,nd,nt,i,j,k,n
character(2048) :: input,output,msg
real(R8Ki) :: pos(3),ori(3,3),v(18)
call get_command_argument(1,input)
call get_command_argument(2,output)
call NWTC_Init()
open(newunit=fi,file=trim(input),status='old')
open(newunit=fo,file=trim(output),status='replace')
read(fi,*) ns,nd,nt
call MeshCreate(sm,COMPONENT_OUTPUT,ns,err,msg,TranslationDisp=.true.,Orientation=.true.,TranslationVel=.true.,RotationVel=.true.)
call check()
call MeshCreate(dm,COMPONENT_INPUT,nd,err,msg,TranslationDisp=.true.,Orientation=.true.,TranslationVel=.true.,RotationVel=.true.)
call check()
call MeshCreate(sl,COMPONENT_INPUT,ns-2,err,msg,Force=.true.,Moment=.true.)
call check()
do j=1,ns
  read(fi,*) pos,((ori(i,k),k=1,3),i=1,3)
  call MeshPositionNode(sm,j,pos,err,msg,ori)
  call check()
  if(j>1.and.j<ns) then
    call MeshPositionNode(sl,j-1,pos,err,msg,ori)
    call check()
    call MeshConstructElement(sl,ELEMENT_POINT,err,msg,j-1)
    call check()
  endif
enddo
do j=1,nd
  read(fi,*) pos,((ori(i,k),k=1,3),i=1,3)
  call MeshPositionNode(dm,j,pos,err,msg,ori)
  call check()
enddo
do j=1,ns-1
  call MeshConstructElement(sm,ELEMENT_LINE2,err,msg,j,j+1)
  call check()
enddo
do j=1,nd-1
  call MeshConstructElement(dm,ELEMENT_LINE2,err,msg,j,j+1)
  call check()
enddo
call MeshCommit(sm,err,msg)
call check()
call MeshCommit(dm,err,msg)
call check()
call MeshCommit(sl,err,msg)
call check()
call MeshCopy(dm,dl,MESH_SIBLING,err,msg,IOS=COMPONENT_OUTPUT,Force=.true.,Moment=.true.)
call check()
call MeshCopy(sl,sd,MESH_SIBLING,err,msg,IOS=COMPONENT_OUTPUT,TranslationDisp=.true.)
call check()
call MeshMapCreate(sm,dm,motionmap,err,msg)
call check()
call MeshMapCreate(dl,sl,loadmap,err,msg)
call check()
write(fo,'(A)',advance='no') 'kind,sample,node'
do i=0,17
  write(fo,'(A,I0)',advance='no') ',v',i
enddo
write(fo,*)
do n=0,nt-1
  do j=1,ns
    read(fi,*) pos,((sm%Orientation(i,k,j),k=1,3),i=1,3),sm%TranslationVel(:,j),sm%RotationVel(:,j)
    sm%TranslationDisp(:,j)=pos-sm%Position(:,j)
    if(j>1.and.j<ns) sd%TranslationDisp(:,j-1)=sm%TranslationDisp(:,j)
  enddo
  call Transfer_Line2_to_Line2(sm,dm,motionmap,err,msg)
  call check()
  do j=1,nd
    read(fi,*) dl%Force(:,j),dl%Moment(:,j)
    v(1:3)=dm%Position(:,j)+dm%TranslationDisp(:,j)
    v(4:6)=dm%TranslationVel(:,j)
    do i=1,3
      v(7+(i-1)*3:9+(i-1)*3)=dm%Orientation(i,:,j)
    enddo
    v(16:18)=dm%RotationVel(:,j)
    write(fo,'(I0,A,I0,A,I0,18(A,ES25.16E3))') 0,',',n,',',j-1,(',',v(k),k=1,18)
  enddo
  call Transfer_Line2_to_Point(dl,sl,loadmap,err,msg,dm,sd)
  call check()
  do j=1,ns-2
    v=0
    v(1:3)=sl%Force(:,j)
    v(4:6)=sl%Moment(:,j)
    write(fo,'(I0,A,I0,A,I0,18(A,ES25.16E3))') 1,',',n,',',j-1,(',',v(k),k=1,18)
  enddo
  sm%RemapFlag=.false.
  dm%RemapFlag=.false.
  dl%RemapFlag=.false.
  sl%RemapFlag=.false.
enddo
close(fi)
close(fo)
contains
subroutine check()
if(err>=ErrID_Fatal)then
  print *,trim(msg)
  stop 1
endif
end subroutine
end program
