! Diagnostic output only; no change to the reference solver equations.
subroutine TraceTurbine(t, Turbine)
   use ElastoDyn_Parameters, only: DOF_BF, DOF_BE
   real(DbKi), intent(in) :: t
   type(FAST_TurbineType), intent(in) :: Turbine
   integer, save :: unit = 0
   integer :: b,j,k,dofs(3)
   real(R8Ki) :: v(23)
   if(unit==0) then
      open(newunit=unit,file='solver-reference.csv',status='replace')
      write(unit,'(A)') 'time,blade,node,q1,q2,q3,qd1,qd2,qd3,phi,alpha,speed,a,ap,cl,cd,cm,x,y,z,fx,fy,fz,mx,my,mz'
   endif
   do b=1,3
      dofs=(/DOF_BF(b,1),DOF_BF(b,2),DOF_BE(b,1)/)
      v(1:3)=Turbine%ED%x(1,STATE_PRED)%QT(dofs)
      v(4:6)=Turbine%ED%x(1,STATE_PRED)%QDT(dofs)
      do j=1,Turbine%AD%p%rotors(1)%NumBlNds
         associate(y=>Turbine%AD%m%rotors(1)%BEMT_y, &
                   motion=>Turbine%AD%Input(INPUT_CURR)%rotors(1)%BladeMotion(b), &
                   loads=>Turbine%AD%y%rotors(1)%BladeLoad(b))
            v(7:14)=(/y%phi(j,b),y%AOA(j,b),y%Vrel(j,b),y%axInduction(j,b),y%tanInduction(j,b),y%Cl(j,b),y%Cd(j,b),y%Cm(j,b)/)
            v(15:17)=motion%Position(:,j)+motion%TranslationDisp(:,j)
            v(18:20)=loads%Force(:,j)
            v(21:23)=loads%Moment(:,j)
         end associate
         write(unit,'(ES25.16E3,A,I0,A,I0,23(A,ES25.16E3))') t,',',b-1,',',j-1,(',',v(k),k=1,23)
      enddo
   enddo
   flush(unit)
end subroutine
