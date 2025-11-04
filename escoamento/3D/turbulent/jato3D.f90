program jato3D
implicit none
integer,parameter				:: ni=271,nj=119,nk=119,itmax=10000000
integer						:: i,j,k,it,nmj,nmk,imax,jmax,kmax,npj,npk,ninf,nsup,nextinf,nextsup,nmj2,nparede,ninj
real(8)						:: dt,dx,dy,dz,L,D,H,re,ret,dx2,dy2,dz2,sch,di
real(8)						:: uin,vin,pin,Tin,win
real(8),dimension(ni,nj,nk)			:: u,v,p,ua,va,pa,ro,roa,qmx,qmy,qmxa,qmya,T,Ta,w,wa,qmz,qmza,mist,mista,mief,Tn,d1,mi0,dife
real(8)						:: duudx,duvdy,duwdz,dvudx,dvvdy,dvwdz,dwudx,dwvdy,dwwdz,d2udx2,d2udy2,d2udz2
real(8)						:: d2vdx2,d2vdy2,d2vdz2,d2wdx2,d2wdy2,d2wdz2,d2vdydx,d2wdzdx,d2udxdy,d2wdzdy,d2udxdz,d2vdydz
real(8)						:: convx,convy,convz,difx,dify,difz,d2d1dx2,d2d1dy2,dpx,dpy,tpv,vtotal,dpdx,dpdy,dpdz
real(8)						:: dqmxzdx,dqmyzdy,d2zdx2,d2zdy2,coef,dqmzzdz,d2zdz2,dqmxdx,dqmydy,dqmzdz
real(8)						:: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,convecz,difuz,mit
real(8)						:: d2ddx2,d2ddy2,d2ddz2,d2qmxudx2,d2qmyvdy2,d2qmzwdz2,d2qmxvdxy,d2qmxwdxz,d2qmywdyz,dpz,dp
real(8),dimension(ni)				:: x
real(8),dimension(nj)				:: y
real(8),dimension(nk)				:: z
real(8)						:: TIME,sbr,dmax,soma,beta,dD,dparede,fat,Diaminte,Diamext,fat1,Linj,betax1,betax2
real(8)						:: roair,miair,mifuel,rofuel,yo22,yf1,mwf,mwo2,zst,Q,cp,alfa,pi,grav,Lc,U00,micomb
real(8)						:: S11,S12,S13,S21,S22,S23,S31,S32,S33,nfrob,delta,Cs,c,b,csit,nit,mitxa
print *, 'Jato Turbulento 3D'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			Constantes				!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	

	re   = 68000.d0
	Cs   = 0.2d0
	dt   = 1.d-3
	L    = 100.d0
	D    = 30.d0
	di   = 1.d-5
	sbr  = 0.97d0
	roair= 1.1614d0
	miair= 1.846d-5
	mifuel=8.1362d-6   !kg/m s
	rofuel=1.96	   !kg/m³
	sch  = (mifuel)/(di*rofuel)
	Lc   =.0052d0
	grav = 9.81d0
	U00  = 56.7d0


	print *, 'Reynolds:',re
	print *, 'Schimidt:',sch
	print *, 'Smagorinsky',Cs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	             Geração da malha				!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
	ninj  =30.d0
	Linj  =10.d0
	betax1=0.7d0
	betax2=1.3d0
!Direção x

	do i=1,ninj
	  fat1 =(1.d0*(ninj-i))/((ninj-1)*1.d0)
	  x(i) =Linj-((dexp(betax1*fat1)-1.d0)/(dexp(betax1)-1.d0))*Linj
!	  print *, x(i)
	end do

 	do i=ninj+1,ni
  	  fat1 =(1.d0*(i-ninj))/(ni-ninj)	
	  x(i) =Linj+(L-Linj)*((dexp(betax2*fat1)-1.d0)/(dexp(betax2)-1.d0))
!	  print *, x(i)
 	end do
!----------------------------------------------------------------
!Direção y
	beta    =2.4d0
	npj     =5
	npk     =5
	nparede =6
	Diaminte=1.d0
	Diamext =1.66d0*Diaminte

	nmj    = ((nj-1)/2 + 1)
	ninf   = nmj - ((npj-1)/2.d0)
	nsup   = nmj + ((npj-1)/2.d0)
	nextinf= nmj - ((npj-1)/2.d0) - (nparede/2.d0)
	nextsup= nmj + ((npj-1)/2.d0) + (nparede/2.d0)
	dD=Diaminte/(((npj-1)*1.d0)+2.d0)
	Dparede=((Diamext-Diaminte)/2.d0)/((nparede/2.d0)-1.d0)
	y=0.d0
	do j=1,nextinf
	  fat  = (1.d0*(nextinf-j)/(nextinf-1.d0))
	  y(j) = -(Diamext/2.d0+(D-Diamext)/2.d0*(((dexp(beta*fat) - 1.d0)/(dexp(beta) -1.d0))))
	end do

	do j=nextinf+1,ninf-1
	  y(j)=y(j-1)+dparede
	end do

	do j=ninf,nsup+1
	  y(j)=y(j-1)+dD
	end do

	do j=nsup+2,nextsup-1
	  y(j)=y(j-1)+dparede
	end do

	do j=nextsup,nj
	  fat  =(1.d0*(j-(nextsup)))/((nj-(nextsup))*1.d0)
	  y(j) =Diamext/2.d0+(D-Diamext)/2.d0*(((dexp(beta*fat) - 1.d0)/(dexp(beta) -1.d0)))
	end do
!----------------------------------------------------------------
!Direção z
	nmk    = ((nk-1)/2 + 1)
	ninf   = nmk - ((npk-1)/2.d0)
	nsup   = nmk + ((npk-1)/2.d0)
	nextinf= nmk - ((npk-1)/2.d0) - (nparede/2.d0)
	nextsup= nmk + ((npk-1)/2.d0) + (nparede/2.d0)
	dD=Diaminte/(((npk-1)*1.d0)+2.d0)
	Dparede=((Diamext-Diaminte)/2.d0)/((nparede/2.d0)-1.d0)
	z=0.d0
	do k=1,nextinf
	  fat  = (1.d0*(nextinf-k)/(nextinf-1.d0))
	  z(k) = -(Diamext/2.d0+(D-Diamext)/2.d0*(((dexp(beta*fat) - 1.d0)/(dexp(beta) -1.d0))))
	end do

	do k=nextinf+1,ninf-1
	  z(k)=z(k-1)+dparede
	end do

	do k=ninf,nsup+1
	  z(k)=z(k-1)+dD
	end do

	do k=nsup+2,nextsup-1
	  z(k)=z(k-1)+dparede
	end do

	do k=nextsup,nk
	  fat  =(1.d0*(k-(nextsup)))/((nk-(nextsup))*1.d0)
	  z(k) =Diamext/2.d0+(D-Diamext)/2.d0*(((dexp(beta*fat) - 1.d0)/(dexp(beta) -1.d0)))
	end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		   Condições iniciais				!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	

	uin = 1.d0
	vin = 0.d0
	win = 0.d0
	pin = 1.d0

	do k=1,nk
	  do j=1,nj
	    do i=1,ni
	    u(i,j,k)    =0.3*(1.d0-(DSQRT(y(j)**2.d0+z(k)**2.d0)/(D/2.d0))**2.d0)*1.5d0
	    v(i,j,k)    =vin
	    w(i,j,k)    =win
	    mist(i,j,k) =1.d-10
	    mista(i,j,k)=mist(i,j,k)
	    end do
	  end do
	end do

	do k=1,nk
	  do j=1,nj
	    do i=1,ni
	    ro(i,j,k)   =1.d0
	    roa(i,j,k)  =ro(i,j,k)
	    p(i,j,k)    =pin
	    ua(i,j,k)   =u(i,j,k)
	    va(i,j,k)   =v(i,j,k)
	    wa(i,j,k)   =w(i,j,k)
	    pa(i,j,k)   =pin
	    qmx(i,j,k)  =ro(i,j,k)*u(i,j,k)
	    qmy(i,j,k)  =ro(i,j,k)*v(i,j,k)
	    qmz(i,j,k)  =ro(i,j,k)*w(i,j,k)
	    qmxa(i,j,k) =qmx(i,j,k)
	    qmya(i,j,k) =qmy(i,j,k)
	    qmza(i,j,k) =qmz(i,j,k)
	    end do
	  end do
	end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	      Equação Quantidade de Movimento			!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
	
do it=1,itmax
  
	do k=2,nk-1
	  do j=2,nj-1
	    do i=2,ni-1

	      dx      =.5d0*(x(i+1)-x(i-1)) 
	      dy      =.5d0*(y(j+1)-y(j-1)) 
	      dz      =.5d0*(z(k+1)-z(k-1))
	      dx2     = dx*dx
	      dy2     = dy*dy
	      dz2     = dz*dz

	      duudx   = ((qmxa(i+1,j,k)*qmxa(i+1,j,k)/roa(i+1,j,k))-(qmxa(i-1,j,k)*qmxa(i-1,j,k)/roa(i-1,j,k)))/(2.d0*dx)
	      duvdy   = ((qmxa(i,j+1,k)*qmya(i,j+1,k)/roa(i,j+1,k))-(qmxa(i,j-1,k)*qmya(i,j-1,k)/roa(i,j-1,k)))/(2.d0*dy)
	      duwdz   = ((qmxa(i,j,k+1)*qmza(i,j,k+1)/roa(i,j,k+1))-(qmxa(i,j,k-1)*qmza(i,j,k-1)/roa(i,j,k-1)))/(2.d0*dz)

	      dvudx   = ((qmya(i+1,j,k)*qmxa(i+1,j,k)/roa(i+1,j,k))-(qmya(i-1,j,k)*qmxa(i-1,j,k)/roa(i-1,j,k)))/(2.d0*dx)
	      dvvdy   = ((qmya(i,j+1,k)*qmya(i,j+1,k)/roa(i,j+1,k))-(qmya(i,j-1,k)*qmya(i,j-1,k)/roa(i,j-1,k)))/(2.d0*dy)
	      dvwdz   = ((qmya(i,j,k+1)*qmza(i,j,k+1)/roa(i,j,k+1))-(qmya(i,j,k-1)*qmza(i,j,k-1)/roa(i,j,k-1)))/(2.d0*dz)

	      dwudx   = ((qmza(i+1,j,k)*qmxa(i+1,j,k)/roa(i+1,j,k))-(qmza(i-1,j,k)*qmxa(i-1,j,k)/roa(i-1,j,k)))/(2.d0*dx)
	      dwvdy   = ((qmza(i,j+1,k)*qmya(i,j+1,k)/roa(i,j+1,k))-(qmza(i,j-1,k)*qmya(i,j-1,k)/roa(i,j-1,k)))/(2.d0*dy)     
	      dwwdz   = ((qmza(i,j,k+1)*qmza(i,j,k+1)/roa(i,j,k+1))-(qmza(i,j,k-1)*qmza(i,j,k-1)/roa(i,j,k-1)))/(2.d0*dz)

	      d2udx2  = ((qmxa(i+1,j,k)/roa(i+1,j,k))-2.d0*(qmxa(i,j,k)/roa(i,j,k))+(qmxa(i-1,j,k)/roa(i-1,j,k)))/(dx2)
	      d2udy2  = ((qmxa(i,j+1,k)/roa(i,j+1,k))-2.d0*(qmxa(i,j,k)/roa(i,j,k))+(qmxa(i,j-1,k)/roa(i,j-1,k)))/(dy2)
	      d2udz2  = ((qmxa(i,j,k+1)/roa(i,j,k+1))-2.d0*(qmxa(i,j,k)/roa(i,j,k))+(qmxa(i,j,k-1)/roa(i,j,k-1)))/(dz2)

	      d2vdx2  = ((qmya(i+1,j,k)/roa(i+1,j,k))-2.d0*(qmya(i,j,k)/roa(i,j,k))+(qmya(i-1,j,k)/roa(i-1,j,k)))/(dx2)
	      d2vdy2  = ((qmya(i,j+1,k)/roa(i,j+1,k))-2.d0*(qmya(i,j,k)/roa(i,j,k))+(qmya(i,j-1,k)/roa(i,j-1,k)))/(dy2)
	      d2vdz2  = ((qmya(i,j,k+1)/roa(i,j,k+1))-2.d0*(qmya(i,j,k)/roa(i,j,k))+(qmya(i,j,k-1)/roa(i,j,k-1)))/(dz2)

	      d2wdx2  = ((qmza(i+1,j,k)/roa(i+1,j,k))-2.d0*(qmza(i,j,k)/roa(i,j,k))+(qmza(i-1,j,k)/roa(i-1,j,k)))/(dx2)
	      d2wdy2  = ((qmza(i,j+1,k)/roa(i,j+1,k))-2.d0*(qmza(i,j,k)/roa(i,j,k))+(qmza(i,j-1,k)/roa(i,j-1,k)))/(dy2)
	      d2wdz2  = ((qmza(i,j,k+1)/roa(i,j,k+1))-2.d0*(qmza(i,j,k)/roa(i,j,k))+(qmza(i,j,k-1)/roa(i,j,k-1)))/(dz2)

	      d2vdydx = ((qmya(i+1,j+1,k)/roa(i+1,j+1,k))-(qmya(i+1,j-1,k)/roa(i+1,j-1,k))-(qmya(i-1,j+1,k)/roa(i-1,j+1,k))&
				  &+(qmya(i-1,j-1,k)/roa(i-1,j-1,k)))*(1.d0/(4.d0*dx*dy))
	      d2wdzdx = ((qmza(i+1,j,k+1)/roa(i+1,j,k+1))-(qmza(i+1,j,k-1)/roa(i+1,j,k-1))-(qmza(i-1,j,k+1)/roa(i-1,j,k+1))&
				  &+(qmza(i-1,j,k-1)/roa(i-1,j,k-1)))*(1.d0/(4.d0*dx*dz))

	      d2udxdy = ((qmxa(i+1,j+1,k)/roa(i+1,j+1,k))-(qmxa(i+1,j-1,k)/roa(i+1,j-1,k))-(qmxa(i-1,j+1,k)/roa(i-1,j+1,k))&
				  &+(qmxa(i-1,j-1,k)/roa(i-1,j-1,k)))*(1.d0/(4.d0*dy*dx))
	      d2wdzdy = ((qmza(i,j+1,k+1)/roa(i,j+1,k+1))-(qmza(i,j+1,k-1)/roa(i,j+1,k-1))-(qmza(i,j-1,k+1)/roa(i,j-1,k+1))&
				  &+(qmza(i,j-1,k-1)/roa(i,j-1,k-1)))*(1.d0/(4.d0*dy*dz))

	      d2udxdz = ((qmxa(i+1,j,k+1)/roa(i+1,j,k+1))-(qmxa(i+1,j,k-1)/roa(i+1,j,k-1))-(qmxa(i-1,j,k+1)/roa(i-1,j,k+1))&
				  &+(qmxa(i-1,j,k-1)/roa(i-1,j,k-1)))*(1.d0/(4.d0*dz*dx))
	      d2vdydz = ((qmya(i,j+1,k+1)/roa(i,j+1,k+1))-(qmya(i,j+1,k-1)/roa(i,j+1,k-1))-(qmya(i,j-1,k+1)/roa(i,j-1,k+1))&
				  &+(qmya(i,j-1,k-1)/roa(i,j-1,k-1)))*(1.d0/(4.d0*dz*dy))

	      dqmxdx    = (qmxa(i+1,j,k)-qmxa(i-1,j,k))/(2.d0*dx)
	      dqmydy    = (qmya(i,j+1,k)-qmya(i,j-1,k))/(2.d0*dy)
	      dqmzdz    = (qmza(i,j,k+1)-qmza(i,j,k-1))/(2.d0*dz)

	      d1(i,j,k) = dqmxdx+dqmydy+dqmzdz

	      dpdx    = (pa(i+1,j,k) - pa(i-1,j,k))/(2.d0*dx)
	      dpdy    = (pa(i,j+1,k) - pa(i,j-1,k))/(2.d0*dy)
	      dpdz    = (pa(i,j,k+1) - pa(i,j,k-1))/(2.d0*dz)

	      convx   = duudx + duvdy + duwdz
	      convy   = dvudx + dvvdy + dvwdz
	      convz   = dwudx + dwvdy + dwwdz

	      dudx    =(qmxa(i+1,j,k)/roa(i+1,j,k)-qmxa(i-1,j,k)/roa(i-1,j,k))/(2.d0*dx) !
	      dudy    =(qmxa(i,j+1,k)/roa(i,j+1,k)-qmxa(i,j-1,k)/roa(i,j-1,k))/(2.d0*dy) !  para j ok
	      dudz    =(qmxa(i,j,k+1)/roa(i,j,k+1)-qmxa(i,j,k-1)/roa(i,j,k-1))/(2.d0*dz) !  para k ok	

	      dvdx    =(qmya(i+1,j,k)/roa(i+1,j,k)-qmya(i-1,j,k)/roa(i-1,j,k))/(2.d0*dx)
	      dvdy    =(qmya(i,j+1,k)/roa(i,j+1,k)-qmya(i,j-1,k)/roa(i,j-1,k))/(2.d0*dy)
	      dvdz    =(qmya(i,j,k+1)/roa(i,j,k+1)-qmya(i,j,k-1)/roa(i,j,k-1))/(2.d0*dz)

	      dwdx    =(qmza(i+1,j,k)/roa(i+1,j,k)-qmza(i-1,j,k)/roa(i-1,j,k))/(2.d0*dx)
	      dwdy    =(qmza(i,j+1,k)/roa(i,j+1,k)-qmza(i,j-1,k)/roa(i,j-1,k))/(2.d0*dy)
	      dwdz    =(qmza(i,j,k+1)/roa(i,j,k+1)-qmza(i,j,k-1)/roa(i,j,k-1))/(2.d0*dz)
	      delta   = (dx*dy*dz)**(1.d0/3.d0)

		 S11		=DABS(2.d0*dudx)
		 S12		=DABS(dudy+dvdx)
		 S13		=DABS(dudz+dwdx)
		 S21		=DABS(dvdx+dudy)
		 S22		=DABS(2.d0*dvdy)
		 S23		=DABS(dvdz+dwdy)
		 S31		=DABS(dwdx+dudz)
		 S32		=DABS(dwdy+dvdy)
		 S33		=DABS(2.d0*dwdz)	

		 nfrob		=DSQRT(0.5d0*(S11*(S11+S12+S13+S21+S22+S23+S31+S32+S33)&
					    &+S12*(S11+S12+S13+S21+S22+S23+S31+S32+S33)&
					    &+S13*(S11+S12+S13+S21+S22+S23+S31+S32+S33)&
					    &+S21*(S11+S12+S13+S21+S22+S23+S31+S32+S33)&
					    &+S22*(S11+S12+S13+S21+S22+S23+S31+S32+S33)&
					    &+S23*(S11+S12+S13+S21+S22+S23+S31+S32+S33)&
					    &+S31*(S11+S12+S13+S21+S22+S23+S31+S32+S33)&
					    &+S32*(S11+S12+S13+S21+S22+S23+S31+S32+S33)&
					    &+S33*(S11+S12+S13+S21+S22+S23+S31+S32+S33)))


	      mit     = ro(i,j,k)*((Cs*delta)**2.d0)*nfrob
	      micomb  = mifuel*mist(i,j,k)+miair*(1.d0-mist(i,j,k))
              mief(i,j,k)=(mit+micomb)/mifuel

	      difx    = -dpdx + (mief(i,j,k)/re)*(d2udx2 + d2udy2 + d2udz2 + (1.d0/3.d0)*(d2udx2 + d2vdydx + d2wdzdx))+&
					&(-grav*Lc/(U00**2.d0))*ro(i,j,k)
	      dify    = -dpdy + (mief(i,j,k)/re)*(d2vdx2 + d2vdy2 + d2vdz2 + (1.d0/3.d0)*(d2udxdy + d2vdy2 + d2wdzdy))
	      difz    = -dpdz + (mief(i,j,k)/re)*(d2wdx2 + d2wdy2 + d2wdz2 + (1.d0/3.d0)*(d2udxdz + d2vdydz + d2wdz2))

	      qmx(i,j,k)=qmxa(i,j,k)+dt*(-convx+difx)
	      qmy(i,j,k)=qmya(i,j,k)+dt*(-convy+dify)
	      qmz(i,j,k)=qmza(i,j,k)+dt*(-convz+difz)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		      Fração de mistura				!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	      dqmxzdx= (qmxa(i+1,j,k)*mista(i+1,j,k) - qmxa(i-1,j,k)*mista(i-1,j,k))/(2.d0*dx)
	      dqmyzdy= (qmya(i,j+1,k)*mista(i,j+1,k) - qmya(i,j-1,k)*mista(i,j-1,k))/(2.d0*dy)
	      dqmzzdz= (qmza(i,j,k+1)*mista(i,j,k+1) - qmza(i,j,k-1)*mista(i,j,k-1))/(2.d0*dz)

	      d2zdx2 = (mista(i+1,j,k) - 2.d0*(mista(i,j,k)) + mista(i-1,j,k))/dx2
	      d2zdy2 = (mista(i,j+1,k) - 2.d0*(mista(i,j,k)) + mista(i,j-1,k))/dy2
	      d2zdz2 = (mista(i,j,k+1) - 2.d0*(mista(i,j,k)) + mista(i,j,k-1))/dz2

	      difuz  = (d2zdx2+d2zdy2+d2zdz2)
	      convecz= (dqmxzdx+dqmyzdy+dqmzzdz) 

	      dife(i,j,k)=mief(i,j,k)/(roair*sch)

	      mist(i,j,k) =(roa(i,j,k)*mista(i,j,k)+dt*(-convecz+difuz*((roa(i,j,k)*dife(i,j,k))/(re*sch))))/roa(i,j,k)
	    end do
	  end do
	end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		   Condições de Contorno			!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	

!Entrada e Saída Externas
  do k=1,nk
      do j=1,nj
        qmx(1,j,k)	=.1736d0*uin*ro(i,j,k)!uin*ro(i,j,k)!*ro(i,j,k)*(1.d0-(y(j)/(D/2.d0))**2.d0)*1.5d0!1.d0!
        qmx(ni,j,k)	=qmx(ni-1,j,k)
        qmy(1,j,k)	=0.d0
        qmy(ni,j,k)	=qmy(ni-1,j,k)
	qmz(1,j,k)	=0.d0
	qmz(ni,j,k)	=qmz(ni-1,j,k)
	mist(1,j,k)     =0.d0
	mist(ni,j,k)    =mist(ni-1,j,k)
      enddo
   enddo
!-----------------------------------------------------------------
! DEFININDO O JATO - Primeiro
!-----------------------------------------------------------------

do i=1,ninj
! PAREDE INFERIOR DO INJETOR 
      do k=nextinf,nextsup	
!	Parede Inferior do Jato
     	   do j=nextinf,ninf-1
     		qmx(i,j,k)	=0.d0
     		qmy(i,j,k)	=0.d0
     		qmz(i,j,k)	=0.d0
		mist(i,j,k)     =0.d0
!		mief(i,j,k)     =0.d0
    	   enddo
!	Pardede Superior do jato
	   do j=nsup+2,nextsup
     		qmx(i,j,k)	=0.d0
     		qmy(i,j,k)	=0.d0
     		qmz(i,j,k)	=0.d0
		mist(i,j,k)     =0.d0
!		mief(i,j,k)     =0.d0
	   enddo
      enddo
!	Pardede Frontal do jato
      do j=nextinf,nextsup
     	   do k=nextinf,ninf-1
     		qmx(i,j,k)	=0.d0
     		qmy(i,j,k)	=0.d0
     		qmz(i,j,k)	=0.d0
		mist(i,j,k)     =0.d0
!		mief(i,j,k)     =0.d0
	   enddo
!	Parede Fundo do jato
	   do k=nsup+2,nextsup
     		qmx(i,j,k)	=0.d0
     		qmy(i,j,k)	=0.d0
     		qmz(i,j,k)	=0.d0
		mist(i,j,k)     =0.d0
!		mief(i,j,k)     =0.d0
	   enddo
      end do

! JATO DE COMBUSTÍVEL 
      do k=ninf-1,nsup+1
	do j=ninf-1,nsup+1
     		qmx(i,j,k)	=uin*ro(i,j,k)*1.2345d0*(1.d0- DSQRT(   y(j)**2.d0  +  z(k)**2.d0  )  /0.71d0 ) **(1.d0/6.8d0)
     		qmy(i,j,k)	=0.d0
     		qmz(i,j,k)	=0.d0
		mist(i,j,k)     =1.d0
	end do
      enddo
end do
!Paredes do Duto
  do k=1,nk
    do i=1,ni
      	qmx(i,1,k)	=0.d0!.75d0*qmx(i,2,k)+.25d0*qmx(i,3,k)!0.d0!
        qmx(i,nj,k)	=0.d0!.75d0*qmx(i,nj-1,k)+.25d0*qmx(i,nj-2,k)!0.d0!
      	qmy(i,1,k)	=0.d0
        qmy(i,nj,k)	=0.d0
      	qmz(i,1,k)	=0.d0
        qmz(i,nj,k)	=0.d0
	mist(i,1,k)     =.75d0*mist(i,2,k)+.25d0*mist(i,3,k)
	mist(i,nj,k)    =.75d0*mist(i,nj-1,k)+.25d0*mist(i,nj-2,k)
    enddo
  enddo
  do j=1,nj
    do i=1,ni
      	qmx(i,j,1)	=0.d0!.75d0*qmx(i,j,2)+.25d0*qmx(i,j,3)!0.d0!
        qmx(i,j,nk)	=0.d0!.75d0*qmx(i,j,nk-1)+.25d0*qmx(i,j,nk-2)!0.d0!
      	qmy(i,j,1)	=0.d0
        qmy(i,j,nk)	=0.d0
      	qmz(i,j,1)	=0.d0
        qmz(i,j,nk)	=0.d0
	mist(i,j,1)     =.75d0*mist(i,j,2)+.25d0*mist(i,j,3)
	mist(i,j,nk)    =.75d0*mist(i,j,k-1)+.25d0*mist(i,j,k-2)
    enddo
  enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	         Equação da Pressão				!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	do k=2,nk-1
	  do j=2,nj-1
	    do i=2,ni-1

		 dx		=.5d0*(x(i+1)-x(i-1)) 
		 dy		=.5d0*(y(j+1)-y(j-1)) 
		 dz		=.5d0*(z(k+1)-z(k-1)) 
		 dx2		=dx*dx
		 dy2		=dy*dy
		 dz2		=dz*dz

		 d2ddx2		=(d1(i+1,j,k)/ro(i+1,j,k)-2.d0*d1(i,j,k)/ro(i,j,k)+d1(i-1,j,k)/ro(i-1,j,k))/(dx**2.d0)
		 d2ddy2		=(d1(i,j+1,k)/ro(i,j+1,k)-2.d0*d1(i,j,k)/ro(i,j,k)+d1(i,j-1,k)/ro(i,j-1,k))/(dy**2.d0)
		 d2ddz2		=(d1(i,j,k+1)/ro(i,j,k+1)-2.d0*d1(i,j,k)/ro(i,j,k)+d1(i,j,k-1)/ro(i,j,k-1))/(dz**2.d0)

		 d2qmxudx2	=((qmx(i+1,j,k)*qmx(i+1,j,k)/ro(i+1,j,k))-2.d0*(qmx(i,j,k)*qmx(i,j,k)/ro(i,j,k))+(qmx(i-1,j,k)*&
				 &qmx(i-1,j,k)/ro(i-1,j,k)))/(dx2)
		 d2qmyvdy2	=((qmy(i,j+1,k)*qmy(i,j+1,k)/ro(i,j+1,k))-2.d0*(qmy(i,j,k)*qmy(i,j,k)/ro(i,j,k))+(qmy(i,j-1,k)*&
				 &qmy(i,j-1,k)/ro(i,j-1,k)))/(dy2)
		 d2qmzwdz2	=((qmz(i,j,k+1)*qmz(i,j,k+1)/ro(i,j,k+1))-2.d0*(qmz(i,j,k)*qmz(i,j,k)/ro(i,j,k))+(qmz(i,j,k-1)*&
				 &qmz(i,j,k-1)/ro(i,j,k-1)))/(dz2)

		 d2qmxvdxy	=((qmx(i+1,j+1,k)*qmy(i+1,j+1,k)/ro(i+1,j+1,k))-(qmx(i+1,j-1,k)*qmy(i+1,j-1,k)/ro(i+1,j-1,k))&
				&-(qmx(i-1,j+1,k)*qmy(i-1,j+1,k)/ro(i-1,j+1,k))+(qmx(i-1,j-1,k)*qmy(i-1,j-1,k)/ro(i-1,j-1,k)))/(4.d0*dx*dy)

		 d2qmxwdxz	=((qmx(i+1,j,k+1)*qmz(i+1,j,k+1)/ro(i+1,j,k+1))-(qmx(i+1,j,k-1)*qmz(i+1,j,k-1)/ro(i+1,j,k-1))&
				&-(qmx(i-1,j,k+1)*qmz(i-1,j,k+1)/ro(i-1,j,k+1))+(qmx(i-1,j,k-1)*qmz(i-1,j,k-1)/ro(i-1,j,k-1)))/(4.d0*dx*dz)

		 d2qmywdyz	=((qmy(i,j+1,k+1)*qmz(i,j+1,k+1)/ro(i,j+1,k+1))-(qmy(i,j+1,k-1)*qmz(i,j+1,k-1)/ro(i,j+1,k-1))&
				&-(qmy(i,j-1,k+1)*qmz(i,j-1,k+1)/ro(i,j-1,k+1))+(qmy(i,j-1,k-1)*qmz(i,j-1,k-1)/ro(i,j-1,k-1)))/(4.d0*dy*dz)

		 dpx		= (pa(i+1,j,k)+pa(i-1,j,k))/(dx2)
		 dpy		= (pa(i,j+1,k)+pa(i,j-1,k))/(dy2)
		 dpz		= (pa(i,j,k+1)+pa(i,j,k-1))/(dz2)

		 coef		=((dx2*dy2*dz2))/(2.d0*(dy2*dz2+dx2*dz2+dx2*dy2))
		 dp		=dpx+dpy+dpz
		 tpv		=d2qmxudx2+d2qmyvdy2+d2qmzwdz2+2.d0*(d2qmxvdxy+d2qmxwdxz+d2qmywdyz)	 

		 p(i,j,k)       =coef*(dp+tpv-(4.d0/(3.d0*re))*(d2ddx2+d2ddy2+d2ddz2)-(1.d0/dt)*d1(i,j,k))
		 p(i,j,k)       =sbr*pa(i,j,k)+(1.d0-sbr)*p(i,j,k)

	      end do
	    end do
	  end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		  Cond Contorno pressão				!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	

!Paredes do Duto
  do k=1,nk
    do i=1,ni
      	p(i,1,k)	=.75d0*p(i,2,k)+.25d0*p(i,3,k)
        p(i,nj,k)	=.75d0*p(i,nj-1,k)+.25d0*p(i,nj-2,k)
    enddo
  enddo
  do j=1,nj
    do i=1,ni
      	p(i,j,1)	=.75d0*p(i,j,2)+.25d0*p(i,j,3)
        p(i,j,nk)	=.75d0*p(i,j,nk-1)+.25d0*p(i,j,nk-2)
    enddo
  enddo

!Entrada e Saída 
  do k=1,nk
      do j=1,nj
        p(1,j,k)	=p(2,j,k)!
        p(ni,j,k)	=pin!
      enddo
  enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	         Renovação das Variáveis			!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	

	do k=1,nk
	  do j=1,nj
	    do i=1,ni
	      mista(i,j,k)= mist(i,j,k) 
              ro(i,j,k)   = 1.d0*mist(i,j,k)+(roair/rofuel)*(1.d0-mist(i,j,k))
	      roa(i,j,k)  = ro(i,j,k)
              pa(i,j,k)	  = p(i,j,k) 
	      qmxa(i,j,k) = qmx(i,j,k)
	      qmya(i,j,k) = qmy(i,j,k)
	      qmza(i,j,k) = qmz(i,j,k)
	      ua(i,j,k)   = qmxa(i,j,k)/ro(i,j,k)
	      va(i,j,k)   = qmya(i,j,k)/ro(i,j,k)
	      wa(i,j,k)   = qmza(i,j,k)/ro(i,j,k)
	    end do
	  end do
	end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	           Impressão dos dados				!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
!o valor dmax indica o erro de acordo com a conservação da massa d1

           if(mod(it,100).eq.0)then
             dmax=0.d0
	     jmax=0
	     imax=0
	     kmax=0
	     do k=1,nk
	       do j=1,nj
	         do i=1,ni
		   if(d1(i,j,k).gt.dmax)then
		     dmax=d1(i,j,k)
		     imax=i
		     jmax=j
	             kmax=k
	           end if
	         end do
	       end do	
	     end do
	         print *, it,dmax,imax,jmax,kmax
!	           if(dmax.le.1.d-20 .or.dmax.ge.1.d3)stop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		   Gravação de dados				!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	

           if(mod(it,100).eq.0)then
	     open(11,file='mesh.dat',status='replace')
	     open(12,file='u.dat',status='replace')
	     open(13,file='v.dat',status='replace')
	     open(14,file='iso.dat',status='replace')
	     open(15,file='p.dat',status='replace')
	     open(16,file='mixture.dat',status='replace')
	     open(17,file='mixcenter.dat',status='replace')
	     open(18,file='velcenter.dat',status='replace')
	     open(21,file='w.dat',status='replace')
	     open(23,file='viscosidade.dat',status='replace')


	     write(11,*)ni,nj,nk
	     write(12,*)ni,nj,1
	     write(13,*)ni,nj,1	  
	     write(14,*)ni,nj,1	  
	     write(15,*)ni,nj,1	  
	     write(16,*)ni,nj,1
	     write(17,*)ni,1
	     write(18,*)ni,1	  
	     write(21,*)ni,nj,1
	     write(23,*)ni,nj,1

	       do j=1,nj
	         do i=1,ni
	           write(11,'(3e20.10)')x(i),y(j),z(nmk)
		 end do
	       end do	           

	       do j=1,nj
	         do i=1,ni
	           vtotal=((ua(i,j,nmk)**2.d0)+(va(i,j,nmk)**2.d0)+(wa(i,j,nmk)**2.d0))**(1.d0/2.d0)
	             write(12,'(3e20.11)')x(i),y(j),ua(i,j,nmk)
	             write(13,'(3e20.11)')x(i),y(j),va(i,j,nmk)
	             write(14,'(3e20.11)')x(i),y(j),vtotal
                     write(15,'(3e20.11)')x(i),y(j),pa(i,j,nmk)
	             write(21,'(3e20.11)')x(i),y(j),wa(i,j,nmk)
		     write(16,'(3e20.11)')x(i),y(j),mista(i,j,nmk)
		     write(23,'(3e20.11)')x(i),y(j),mief(i,j,nmk)
	         end do
	       end do
	       
	       do i=1,ni
	  	 write(18,'(2e20.11)')x(i),u(i,nmj,nmk)
	         write(17,'(2e20.11)')x(i),mist(i,nmj,nmk)
	       end do
	      
	           close(11)
	           close(12)
		   close(13)
		   close(14)
		   close(15)
		   close(16)
		   close(17)
		   close(18)
	           close(21)
	           close(23)
           end if
      end if

end do !término processo iterativo

call CPU_TIME(TIME)
print *, 'Execution Time:', TIME/60.d0, 'minutes'

end program jato3D

