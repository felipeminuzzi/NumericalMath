program jato3D
implicit none
integer,parameter			:: ni=161,nj=51,itmax=3000000
integer						:: i,j,k,it,nmj,nmk,imax,jmax,kmax,npj,npk,ninf,nsup,nextinf,nextsup,nmj2,nparede,ninj
real(8)						:: dt,dt2,dx,dy,dz,L,D,H,re,dx2,dy2,dz2,sch,di,k1x,k1y,k1z,k2x,k2y,k2z,k3x,k3y,k3z,k4x,k4y,k4z
real(8)						:: uin,vin,pin,Tin,win,k1mis,k2mis,k3mis,k4mis
real(8),dimension(ni,nj,nk)			:: u,v,p,ua,va,pa,ro,roa,qmx,qmy,qmxa,qmya,T,Ta,w,wa,qmz,qmza,mist,mista,mief,Tn
real(8)						:: duudx,duvdy,duwdz,dvudx,dvvdy,dvwdz,dwudx,dwvdy,dwwdz,d2udx2,d2udy2,d2udz2,dwdz
real(8)						:: d2vdx2,d2vdy2,d2vdz2,d2wdx2,d2wdy2,d2wdz2,d2vdydx,d2wdzdx,d2udxdy,d2wdzdy,d2udxdz,d2vdydz
real(8)						:: convx,convy,convz,difx,dify,difz,d2d1dx2,d2d1dy2,dpx,dpy,tpv,vtotal,dpdx,dpdy,dpdz
real(8)						:: dqmxzdx,dqmyzdy,d2zdx2,d2zdy2,dudx,dvdy,coef,dqmzzdz,d2zdz2,d2zdy2k4,d2zdz2k4,d2zdz2k2
real(8),dimension(ni,nj,nk)			:: d1,mi0,dife
real(8)						:: d2ddx2,d2ddy2,d2ddz2,d2qmxudx2,d2qmyvdy2,d2qmzwdz2,d2qmxvdxy,d2qmxwdxz,d2qmywdyz,dpz,dp
real(8),dimension(ni)				:: x
real(8),dimension(nj)				:: y
real(8)						:: dqmxzdxk2,dqmyzdyk2,dqmzzdzk2,d2zdx2k2,d2zdy2k2
real(8)						:: dqmxzdxk3,dqmyzdyk3,dqmzzdzk3,d2zdx2k3,d2zdy2k3,d2zdz2k3,dqmxzdxk4,dqmyzdyk4,dqmzzdzk4,d2zdx2k4
real(8)						:: duudxk2,duvdyk2,duwdzk2,dvudxk2,dvvdyk2,dvwdzk2,dwudxk2,dwvdyk2,dwwdzk2,d2udx2k2
real(8)						:: d2udy2k2,d2udz2k2,d2vdx2k2,d2vdy2k2,d2vdz2k2,d2wdx2k2,d2wdy2k2,d2wdz2k2
real(8)						:: duudxk3,duvdyk3,duwdzk3,dvudxk3,dvvdyk3,dvwdzk3,dwudxk3,dwvdyk3,dwwdzk3,d2udx2k3
real(8)						:: d2udy2k3,d2udz2k3,d2vdx2k3,d2vdy2k3,d2vdz2k3,d2wdx2k3,d2wdy2k3,d2wdz2k3
real(8)						:: duudxk4,duvdyk4,duwdzk4,dvudxk4,dvvdyk4,dvwdzk4,dwudxk4,dwvdyk4,dwwdzk4,d2udx2k4
real(8)						:: d2udy2k4,d2udz2k4,d2vdx2k4,d2vdy2k4,d2vdz2k4,d2wdx2k4,d2wdy2k4,d2wdz2k4
real(8)						:: TIME,sbr,dmax,soma,beta,dD,dparede,fat,Diaminte,Diamext,fat1,Linj,betax1
real(8)						:: roh2o,mih2o,yo22,yf1,mwf,mwo2,zst,Q,cp,alfa

print *, 'Jato Laminar 2D'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			Constantes				!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	

	re   = 900.d0
	dt   = 1.d-5
	dt2  = dt/2.d0
	L    = 5.d0
	H    = 1.d0
	sbr  = 0.99d0
	di   = 1.d0!2.27d-5
	roh2o= 1.0 !g/cm3
	mih2o= 0.01 !g/cm s

	print *, 'Reynolds:',re

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	             Geração da malha				!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
	dx = L/ni
	dy = H/nj

	x(1) = 0
	y(1) = 0
!Direção x

	do i=2,ni
	  x(i) = x(i-1) + dx
	end do

!----------------------------------------------------------------
!Direção y

	do j=2,nj
		y(j) = x(j-1) + dy
	end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		   Condições iniciais				!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	

	uin = 1.d0
	vin = 0.d0
	win = 0.d0
	pin = 1.d0
	Tin = 300.d0

	do k=1,nk
	  do j=1,nj
	    do i=1,ni
	    ro(i,j,k)   =1.d0
	    roa(i,j,k)  =ro(i,j,k)
	    u(i,j,k)    =uin!*ro(i,j)*(1.d0-(y(j)/(D/2.d0))**2.d0)*1.5d0
	    v(i,j,k)    =vin
	    w(i,j,k)    =win
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
	    mist(i,j,k) =1.d-10
	    mista(i,j,k)=mist(i,j,k)
	    T(i,j,k)    =Tin
	    Ta(i,j,k)   =Tin
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
!--------------------------------derivadas para k1----------------------------------------------------------------
	      duudx   = ((ua(i+1,j,k)*ua(i+1,j,k)*roa(i+1,j,k))-(ua(i-1,j,k)*ua(i-1,j,k)*roa(i-1,j,k)))/(2.d0*dx)
	      duvdy   = ((ua(i,j+1,k)*va(i,j+1,k)*roa(i,j+1,k))-(ua(i,j-1,k)*va(i,j-1,k)*roa(i,j-1,k)))/(2.d0*dy)
	      duwdz   = ((ua(i,j,k+1)*wa(i,j,k+1)*roa(i,j,k+1))-(ua(i,j,k-1)*wa(i,j,k-1)*roa(i,j,k-1)))/(2.d0*dz)

	      dvudx   = ((va(i+1,j,k)*ua(i+1,j,k)*roa(i+1,j,k))-(va(i-1,j,k)*ua(i-1,j,k)*roa(i-1,j,k)))/(2.d0*dx)
	      dvvdy   = ((va(i,j+1,k)*va(i,j+1,k)*roa(i,j+1,k))-(va(i,j-1,k)*va(i,j-1,k)*roa(i,j-1,k)))/(2.d0*dy)
	      dvwdz   = ((va(i,j,k+1)*wa(i,j,k+1)*roa(i,j,k+1))-(va(i,j,k-1)*wa(i,j,k-1)*roa(i,j,k-1)))/(2.d0*dz)

	      dwudx   = ((wa(i+1,j,k)*ua(i+1,j,k)*roa(i+1,j,k))-(wa(i-1,j,k)*ua(i-1,j,k)*roa(i-1,j,k)))/(2.d0*dx)
	      dwvdy   = ((wa(i,j+1,k)*va(i,j+1,k)*roa(i,j+1,k))-(wa(i,j-1,k)*va(i,j-1,k)*roa(i,j-1,k)))/(2.d0*dy)     
	      dwwdz   = ((wa(i,j,k+1)*wa(i,j,k+1)*roa(i,j,k+1))-(wa(i,j,k-1)*wa(i,j,k-1)*roa(i,j,k-1)))/(2.d0*dz)

	      d2udx2  = (ua(i+1,j,k)-2.d0*(ua(i,j,k))+ua(i-1,j,k))/(dx2)
	      d2udy2  = (ua(i,j+1,k)-2.d0*(ua(i,j,k))+ua(i,j-1,k))/(dy2)
	      d2udz2  = (ua(i,j,k+1)-2.d0*(ua(i,j,k))+ua(i,j,k-1))/(dz2)

	      d2vdx2  = (va(i+1,j,k)-2.d0*(va(i,j,k))+va(i-1,j,k))/(dx2)
	      d2vdy2  = (va(i,j+1,k)-2.d0*(va(i,j,k))+va(i,j-1,k))/(dy2)
	      d2vdz2  = (va(i,j,k+1)-2.d0*(va(i,j,k))+va(i,j,k-1))/(dz2)

	      d2wdx2  = (wa(i+1,j,k)-2.d0*(wa(i,j,k))+wa(i-1,j,k))/(dx2)
	      d2wdy2  = (wa(i,j+1,k)-2.d0*(wa(i,j,k))+wa(i,j-1,k))/(dy2)
	      d2wdz2  = (wa(i,j,k+1)-2.d0*(wa(i,j,k))+wa(i,j,k-1))/(dz2)

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

	      dudx    = (roa(i+1,j,k)*ua(i+1,j,k)-roa(i-1,j,k)*ua(i-1,j,k))/(2.d0*dx)
	      dvdy    = (roa(i,j+1,k)*va(i,j+1,k)-roa(i,j-1,k)*va(i,j-1,k))/(2.d0*dy)
	      dwdz    = (roa(i,j,k+1)*wa(i,j,k+1)-roa(i,j,k-1)*wa(i,j,k-1))/(2.d0*dz)

	      d1(i,j,k) = dudx + dvdy + dwdz

	      dpdx    = (pa(i+1,j,k) - pa(i-1,j,k))/(2.d0*dx)
	      dpdy    = (pa(i,j+1,k) - pa(i,j-1,k))/(2.d0*dy)
	      dpdz    = (pa(i,j,k+1) - pa(i,j,k-1))/(2.d0*dz)

              mi0(i,j,k)   =mief(i,j,k)
              mief(i,j,k)  =mi0(i,j,k)*(sbr)+((T(i,j,k)/300.d0)**1.d0)*(1.d0-sbr)

	      k1x     = -duudx-duvdy-duwdz-dpdx + (mief(i,j,k)/re)*(d2udx2 + d2udy2 + d2udz2 + &
			&(1.d0/3.d0)*(d2udx2 + d2vdydx + d2wdzdx))
	      k1y     = -dvudx-dvvdy-dvwdz-dpdy + (mief(i,j,k)/re)*(d2vdx2 + d2vdy2 + d2vdz2 + &
			&(1.d0/3.d0)*(d2udxdy + d2vdy2 + d2wdzdy))
	      k1z     = -dwudx-dwvdy-dwwdz-dpdz + (mief(i,j,k)/re)*(d2wdx2 + d2wdy2 + d2wdz2 + &
			&(1.d0/3.d0)*(d2udxdz + d2vdydz + d2wdz2))
!--------------------------------derivadas para k2----------------------------------------------------------------
	      duudxk2  = (((ua(i+1,j,k)+(dt2*k1x))*(ua(i+1,j,k)+(dt2*k1x))*roa(i+1,j,k))-((ua(i-1,j,k)+(dt2*k1x))&
	      			&*(ua(i-1,j,k)+(dt2*k1x))*roa(i-1,j,k)))/(2.d0*dx)
	      duvdyk2  = (((ua(i,j+1,k)+(dt2*k1x))*va(i,j+1,k)*roa(i,j+1,k))-((ua(i,j-1,k)+(dt2*k1x))*va(i,j-1,k)*roa(i,j-1,k)))/(2.d0*dy)
	      duwdzk2  = (((ua(i,j,k+1)+(dt2*k1x))*wa(i,j,k+1)*roa(i,j,k+1))-((ua(i,j,k-1)+(dt2*k1x))*wa(i,j,k-1)*roa(i,j,k-1)))/(2.d0*dz)

	      dvudxk2  = (((va(i+1,j,k)+(dt2*k1y))*ua(i+1,j,k)*roa(i+1,j,k))-((va(i-1,j,k)+(dt2*k1y))*ua(i-1,j,k)*roa(i-1,j,k)))/(2.d0*dx)
	      dvvdyk2  = (((va(i,j+1,k)+(dt2*k1y))*(va(i,j+1,k)+(dt2*k1y))*roa(i,j+1,k))-((va(i,j-1,k)+(dt2*k1y))&
	      			&*(va(i,j-1,k)+(dt2*k1y))*roa(i,j-1,k)))/(2.d0*dy)
	      dvwdzk2  = (((va(i,j,k+1)+(dt2*k1y))*wa(i,j,k+1)*roa(i,j,k+1))-((va(i,j,k-1)+(dt2*k1y))*wa(i,j,k-1)*roa(i,j,k-1)))/(2.d0*dz)

	      dwudxk2   = (((wa(i+1,j,k)+(dt2*k1z))*ua(i+1,j,k)*roa(i+1,j,k))-((wa(i-1,j,k)+(dt2*k1z))*ua(i-1,j,k)*roa(i-1,j,k)))/(2.d0*dx)
	      dwvdyk2   = (((wa(i,j+1,k)+(dt2*k1z))*va(i,j+1,k)*roa(i,j+1,k))-((wa(i,j-1,k)+(dt2*k1z))*va(i,j-1,k)*roa(i,j-1,k)))/(2.d0*dy)     
	      dwwdzk2   = (((wa(i,j,k+1)+(dt2*k1z))*(wa(i,j,k+1)+(dt2*k1z))*roa(i,j,k+1))-((wa(i,j,k-1)+(dt2*k1z))&
	      			&*(wa(i,j,k-1)+(dt2*k1z))*roa(i,j,k-1)))/(2.d0*dz)

	      d2udx2k2  = ((ua(i+1,j,k)+(dt2*k1x))-2.d0*(ua(i,j,k)+(dt2*k1x))+(ua(i-1,j,k)+(dt2*k1x)))/(dx2)
	      d2udy2k2  = ((ua(i,j+1,k)+(dt2*k1x))-2.d0*(ua(i,j,k)+(dt2*k1x))+(ua(i,j-1,k)+(dt2*k1x)))/(dy2)
	      d2udz2k2  = ((ua(i,j,k+1)+(dt2*k1x))-2.d0*(ua(i,j,k)+(dt2*k1x))+(ua(i,j,k-1)+(dt2*k1x)))/(dz2)

	      d2vdx2k2  = (va(i+1,j,k)-2.d0*(va(i,j,k))+va(i-1,j,k))/(dx2)
	      d2vdy2k2  = (va(i,j+1,k)-2.d0*(va(i,j,k))+va(i,j-1,k))/(dy2)
	      d2vdz2k2  = (va(i,j,k+1)-2.d0*(va(i,j,k))+va(i,j,k-1))/(dz2)

	      d2wdx2k2  = ((wa(i+1,j,k)+(dt2*k1z))-2.d0*(wa(i,j,k)+(dt2*k1z))+(wa(i-1,j,k)+(dt2*k1z)))/(dx2)
	      d2wdy2k2  = ((wa(i,j+1,k)+(dt2*k1z))-2.d0*(wa(i,j,k)+(dt2*k1z))+(wa(i,j-1,k)+(dt2*k1z)))/(dy2)
	      d2wdz2k2  = ((wa(i,j,k+1)+(dt2*k1z))-2.d0*(wa(i,j,k)+(dt2*k1z))+(wa(i,j,k-1)+(dt2*k1z)))/(dz2)

	      k2x  = -duudxk2-duvdyk2-duwdzk2-dpdx + (mief(i,j,k)/re)*(d2udx2k2 + d2udy2k2 + d2udz2k2 + &
			&(1.d0/3.d0)*(d2udx2k2 + d2vdydx + d2wdzdx))
	      k2y  = -dvudxk2-dvvdyk2-dvwdzk2-dpdy + (mief(i,j,k)/re)*(d2vdx2k2 + d2vdy2k2 + d2vdz2k2 + &
			&(1.d0/3.d0)*(d2udxdy + d2vdy2k2 + d2wdzdy))
	      k2z  = -dwudxk2-dwvdyk2-dwwdzk2-dpdz + (mief(i,j,k)/re)*(d2wdx2k2 + d2wdy2k2 + d2wdz2k2 + &
			&(1.d0/3.d0)*(d2udxdz + d2vdydz + d2wdz2k2))
!--------------------------------derivadas para k3----------------------------------------------------------------
	      duudxk3  = (((ua(i+1,j,k)+(dt2*k2x))*(ua(i+1,j,k)+(dt2*k2x))*roa(i+1,j,k))-((ua(i-1,j,k)+(dt2*k2x))&
	      			&*(ua(i-1,j,k)+(dt2*k2x))*roa(i-1,j,k)))/(2.d0*dx)
	      duvdyk3  = (((ua(i,j+1,k)+(dt2*k2x))*va(i,j+1,k)*roa(i,j+1,k))-((ua(i,j-1,k)+(dt2*k1x))*va(i,j-1,k)*roa(i,j-1,k)))/(2.d0*dy)
	      duwdzk3  = (((ua(i,j,k+1)+(dt2*k2x))*wa(i,j,k+1)*roa(i,j,k+1))-((ua(i,j,k-1)+(dt2*k2x))*wa(i,j,k-1)*roa(i,j,k-1)))/(2.d0*dz)

	      dvudxk3  = (((va(i+1,j,k)+(dt2*k2y))*ua(i+1,j,k)*roa(i+1,j,k))-((va(i-1,j,k)+(dt2*k2y))*ua(i-1,j,k)*roa(i-1,j,k)))/(2.d0*dx)
	      dvvdyk3  = (((va(i,j+1,k)+(dt2*k2y))*(va(i,j+1,k)+(dt2*k2y))*roa(i,j+1,k))-((va(i,j-1,k)+(dt2*k2y))&
	      			&*(va(i,j-1,k)+(dt2*k2y))*roa(i,j-1,k)))/(2.d0*dy)
	      dvwdzk3  = (((va(i,j,k+1)+(dt2*k2y))*wa(i,j,k+1)*roa(i,j,k+1))-((va(i,j,k-1)+(dt2*k2y))*wa(i,j,k-1)*roa(i,j,k-1)))/(2.d0*dz)

	      dwudxk3   = (((wa(i+1,j,k)+(dt2*k2z))*ua(i+1,j,k)*roa(i+1,j,k))-((wa(i-1,j,k)+(dt2*k2z))*ua(i-1,j,k)*roa(i-1,j,k)))/(2.d0*dx)
	      dwvdyk3   = (((wa(i,j+1,k)+(dt2*k2z))*va(i,j+1,k)*roa(i,j+1,k))-((wa(i,j-1,k)+(dt2*k2z))*va(i,j-1,k)*roa(i,j-1,k)))/(2.d0*dy)     
	      dwwdzk3   = (((wa(i,j,k+1)+(dt2*k2z))*(wa(i,j,k+1)+(dt2*k2z))*roa(i,j,k+1))-((wa(i,j,k-1)+(dt2*k2z))&
	      			&*(wa(i,j,k-1)+(dt2*k2z))*roa(i,j,k-1)))/(2.d0*dz)

	      d2udx2k3  = ((ua(i+1,j,k)+(dt2*k2x))-2.d0*(ua(i,j,k)+(dt2*k2x))+(ua(i-1,j,k)+(dt2*k2x)))/(dx2)
	      d2udy2k3  = ((ua(i,j+1,k)+(dt2*k2x))-2.d0*(ua(i,j,k)+(dt2*k2x))+(ua(i,j-1,k)+(dt2*k2x)))/(dy2)
	      d2udz2k3  = ((ua(i,j,k+1)+(dt2*k2x))-2.d0*(ua(i,j,k)+(dt2*k2x))+(ua(i,j,k-1)+(dt2*k2x)))/(dz2)

	      d2vdx2k3  = (va(i+1,j,k)-2.d0*(va(i,j,k))+va(i-1,j,k))/(dx2)
	      d2vdy2k3  = (va(i,j+1,k)-2.d0*(va(i,j,k))+va(i,j-1,k))/(dy2)
	      d2vdz2k3  = (va(i,j,k+1)-2.d0*(va(i,j,k))+va(i,j,k-1))/(dz2)

	      d2wdx2k3  = ((wa(i+1,j,k)+(dt2*k2z))-2.d0*(wa(i,j,k)+(dt2*k2z))+(wa(i-1,j,k)+(dt2*k2z)))/(dx2)
	      d2wdy2k3  = ((wa(i,j+1,k)+(dt2*k2z))-2.d0*(wa(i,j,k)+(dt2*k2z))+(wa(i,j-1,k)+(dt2*k2z)))/(dy2)
	      d2wdz2k3  = ((wa(i,j,k+1)+(dt2*k2z))-2.d0*(wa(i,j,k)+(dt2*k2z))+(wa(i,j,k-1)+(dt2*k2z)))/(dz2)

	      k3x  = -duudxk3-duvdyk3-duwdzk3-dpdx + (mief(i,j,k)/re)*(d2udx2k3 + d2udy2k3 + d2udz2k3 + &
			&(1.d0/3.d0)*(d2udx2k3 + d2vdydx + d2wdzdx))
	      k3y  = -dvudxk3-dvvdyk3-dvwdzk3-dpdy + (mief(i,j,k)/re)*(d2vdx2k3 + d2vdy2k3 + d2vdz2k3 + &
			&(1.d0/3.d0)*(d2udxdy + d2vdy2k3 + d2wdzdy))
	      k3z  = -dwudxk3-dwvdyk3-dwwdzk3-dpdz + (mief(i,j,k)/re)*(d2wdx2k3 + d2wdy2k3 + d2wdz2k3 + &
			&(1.d0/3.d0)*(d2udxdz + d2vdydz + d2wdz2k3))
!--------------------------------derivadas para k4----------------------------------------------------------------
	      duudxk4  = (((ua(i+1,j,k)+(dt*k3x))*(ua(i+1,j,k)+(dt*k3x))*roa(i+1,j,k))-((ua(i-1,j,k)+(dt*k3x))&
	      			&*(ua(i-1,j,k)+(dt*k3x))*roa(i-1,j,k)))/(2.d0*dx)
	      duvdyk4  = (((ua(i,j+1,k)+(dt*k3x))*va(i,j+1,k)*roa(i,j+1,k))-((ua(i,j-1,k)+(dt2*k1x))*va(i,j-1,k)*roa(i,j-1,k)))/(2.d0*dy)
	      duwdzk4  = (((ua(i,j,k+1)+(dt*k3x))*wa(i,j,k+1)*roa(i,j,k+1))-((ua(i,j,k-1)+(dt*k3x))*wa(i,j,k-1)*roa(i,j,k-1)))/(2.d0*dz)

	      dvudxk4  = (((va(i+1,j,k)+(dt*k3y))*ua(i+1,j,k)*roa(i+1,j,k))-((va(i-1,j,k)+(dt*k3y))*ua(i-1,j,k)*roa(i-1,j,k)))/(2.d0*dx)
	      dvvdyk4  = (((va(i,j+1,k)+(dt*k3y))*(va(i,j+1,k)+(dt*k3y))*roa(i,j+1,k))-((va(i,j-1,k)+(dt*k3y))&
	      			&*(va(i,j-1,k)+(dt*k3y))*roa(i,j-1,k)))/(2.d0*dy)
	      dvwdzk4  = (((va(i,j,k+1)+(dt*k3y))*wa(i,j,k+1)*roa(i,j,k+1))-((va(i,j,k-1)+(dt*k3y))*wa(i,j,k-1)*roa(i,j,k-1)))/(2.d0*dz)

	      dwudxk4   = (((wa(i+1,j,k)+(dt*k3z))*ua(i+1,j,k)*roa(i+1,j,k))-((wa(i-1,j,k)+(dt*k3z))*ua(i-1,j,k)*roa(i-1,j,k)))/(2.d0*dx)
	      dwvdyk4   = (((wa(i,j+1,k)+(dt*k3z))*va(i,j+1,k)*roa(i,j+1,k))-((wa(i,j-1,k)+(dt*k3z))*va(i,j-1,k)*roa(i,j-1,k)))/(2.d0*dy)     
	      dwwdzk4   = (((wa(i,j,k+1)+(dt*k3z))*(wa(i,j,k+1)+(dt*k3z))*roa(i,j,k+1))-((wa(i,j,k-1)+(dt*k3z))&
	      			&*(wa(i,j,k-1)+(dt*k3z))*roa(i,j,k-1)))/(2.d0*dz)

	      d2udx2k4  = ((ua(i+1,j,k)+(dt*k3x))-2.d0*(ua(i,j,k)+(dt*k3x))+(ua(i-1,j,k)+(dt*k3x)))/(dx2)
	      d2udy2k4  = ((ua(i,j+1,k)+(dt*k3x))-2.d0*(ua(i,j,k)+(dt*k3x))+(ua(i,j-1,k)+(dt*k3x)))/(dy2)
	      d2udz2k4  = ((ua(i,j,k+1)+(dt*k3x))-2.d0*(ua(i,j,k)+(dt*k3x))+(ua(i,j,k-1)+(dt*k3x)))/(dz2)

	      d2vdx2k4  = (va(i+1,j,k)-2.d0*(va(i,j,k))+va(i-1,j,k))/(dx2)
	      d2vdy2k4  = (va(i,j+1,k)-2.d0*(va(i,j,k))+va(i,j-1,k))/(dy2)
	      d2vdz2k4  = (va(i,j,k+1)-2.d0*(va(i,j,k))+va(i,j,k-1))/(dz2)

	      d2wdx2k4  = ((wa(i+1,j,k)+(dt*k3z))-2.d0*(wa(i,j,k)+(dt*k3z))+(wa(i-1,j,k)+(dt*k3z)))/(dx2)
	      d2wdy2k4  = ((wa(i,j+1,k)+(dt*k3z))-2.d0*(wa(i,j,k)+(dt*k3z))+(wa(i,j-1,k)+(dt*k3z)))/(dy2)
	      d2wdz2k4  = ((wa(i,j,k+1)+(dt*k3z))-2.d0*(wa(i,j,k)+(dt*k3z))+(wa(i,j,k-1)+(dt*k3z)))/(dz2)

	      k4x  = -duudxk4-duvdyk4-duwdzk4-dpdx + (mief(i,j,k)/re)*(d2udx2k4 + d2udy2k4 + d2udz2k4 + &
			&(1.d0/3.d0)*(d2udx2k4 + d2vdydx + d2wdzdx))
	      k4y  = -dvudxk4-dvvdyk4-dvwdzk4-dpdy + (mief(i,j,k)/re)*(d2vdx2k4 + d2vdy2k4 + d2vdz2k4 + &
			&(1.d0/3.d0)*(d2udxdy + d2vdy2k4 + d2wdzdy))
	      k4z  = -dwudxk4-dwvdyk4-dwwdzk4-dpdz + (mief(i,j,k)/re)*(d2wdx2k4 + d2wdy2k4 + d2wdz2k4 + &
			&(1.d0/3.d0)*(d2udxdz + d2vdydz + d2wdz2k4))

!--------------------------------velocidades----------------------------------------------------------------------

	      u(i,j,k)    = ua(i,j,k) + (dt/6.d0)*(k1x+2.d0*k2x+2.d0*k3x+k4x)
	      v(i,j,k)    = va(i,j,k) + (dt/6.d0)*(k1y+2.d0*k2y+2.d0*k3y+k4y)
	      w(i,j,k)    = wa(i,j,k) + (dt/6.d0)*(k1z+2.d0*k2z+2.d0*k3z+k4z)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		      Fração de mistura				!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	      dife(i,j,k) =mief(i,j,k)
!--------------------------------derivadas para k1----------------------------------------------------------------
	      dqmxzdx= (roa(i+1,j,k)*ua(i+1,j,k)*mista(i+1,j,k) - roa(i-1,j,k)*ua(i-1,j,k)*mista(i-1,j,k))/(2.d0*dx)
	      dqmyzdy= (roa(i,j+1,k)*va(i,j+1,k)*mista(i,j+1,k) - roa(i,j-1,k)*va(i,j-1,k)*mista(i,j-1,k))/(2.d0*dy)
	      dqmzzdz= (roa(i,j,k+1)*wa(i,j,k+1)*mista(i,j,k+1) - roa(i,j,k-1)*wa(i,j,k-1)*mista(i,j,k-1))/(2.d0*dz)

	      d2zdx2 = (mista(i+1,j,k) - 2.d0*(mista(i,j,k)) + mista(i-1,j,k))/(dx2)
	      d2zdy2 = (mista(i,j+1,k) - 2.d0*(mista(i,j,k)) + mista(i,j-1,k))/(dy2)
	      d2zdz2 = (mista(i,j,k+1) - 2.d0*(mista(i,j,k)) + mista(i,j,k-1))/(dz2)

	      k1mis  = -dqmxzdx-dqmyzdy-dqmzzdz + (dife(i,j,k)/(re*sch))*(d2zdx2+d2zdy2+d2zdz2)
!--------------------------------derivadas para k2----------------------------------------------------------------
	      dqmxzdxk2= (roa(i+1,j,k)*ua(i+1,j,k)*(mista(i+1,j,k)+(dt2*k1mis)) - roa(i-1,j,k)*ua(i-1,j,k)*(mista(i-1,j,k)+&
			&(dt2*k1mis)))/(2.d0*dx)
	      dqmyzdyk2= (roa(i,j+1,k)*va(i,j+1,k)*(mista(i,j+1,k)+(dt2*k1mis)) - roa(i,j-1,k)*va(i,j-1,k)*(mista(i,j-1,k)+&
			&(dt2*k1mis)))/(2.d0*dy)
	      dqmzzdzk2= (roa(i,j,k+1)*wa(i,j,k+1)*(mista(i,j,k+1)+(dt2*k1mis)) - roa(i,j,k-1)*wa(i,j,k-1)*(mista(i,j,k-1)+&
			&(dt2*k1mis)))/(2.d0*dz)

	      d2zdx2k2 = ((mista(i+1,j,k)+(dt2*k1mis)) - 2.d0*(mista(i,j,k)+(dt2*k1mis)) + (mista(i-1,j,k)+(dt2*k1mis)))/(dx2)
	      d2zdy2k2 = ((mista(i,j+1,k)+(dt2*k1mis)) - 2.d0*(mista(i,j,k)+(dt2*k1mis)) + (mista(i,j-1,k)+(dt2*k1mis)))/(dy2)
	      d2zdz2k2 = ((mista(i,j,k+1)+(dt2*k1mis)) - 2.d0*(mista(i,j,k)+(dt2*k1mis)) + (mista(i,j,k-1)+(dt2*k1mis)))/(dz2)

	      k2mis  = -dqmxzdxk2-dqmyzdyk2-dqmzzdzk2 + (dife(i,j,k)/(re*sch))*(d2zdx2k2+d2zdy2k2+d2zdz2k2)
!--------------------------------derivadas para k3----------------------------------------------------------------
	      dqmxzdxk3= (roa(i+1,j,k)*ua(i+1,j,k)*(mista(i+1,j,k)+(dt2*k2mis)) - roa(i-1,j,k)*ua(i-1,j,k)*(mista(i-1,j,k)+&
			&(dt2*k2mis)))/(2.d0*dx)
	      dqmyzdyk3= (roa(i,j+1,k)*va(i,j+1,k)*(mista(i,j+1,k)+(dt2*k2mis)) - roa(i,j-1,k)*va(i,j-1,k)*(mista(i,j-1,k)+&
			&(dt2*k2mis)))/(2.d0*dy)
	      dqmzzdzk3= (roa(i,j,k+1)*wa(i,j,k+1)*(mista(i,j,k+1)+(dt2*k2mis)) - roa(i,j,k-1)*wa(i,j,k-1)*(mista(i,j,k-1)+&
			&(dt2*k2mis)))/(2.d0*dz)

	      d2zdx2k3 = ((mista(i+1,j,k)+(dt2*k2mis)) - 2.d0*(mista(i,j,k)+(dt2*k2mis)) + (mista(i-1,j,k)+(dt2*k2mis)))/(dx2)
	      d2zdy2k3 = ((mista(i,j+1,k)+(dt2*k2mis)) - 2.d0*(mista(i,j,k)+(dt2*k2mis)) + (mista(i,j-1,k)+(dt2*k2mis)))/(dy2)
	      d2zdz2k3 = ((mista(i,j,k+1)+(dt2*k2mis)) - 2.d0*(mista(i,j,k)+(dt2*k2mis)) + (mista(i,j,k-1)+(dt2*k2mis)))/(dz2)

	      k3mis  = -dqmxzdxk3-dqmyzdyk3-dqmzzdzk3 + (dife(i,j,k)/(re*sch))*(d2zdx2k3+d2zdy2k3+d2zdz2k3)
!--------------------------------derivadas para k4----------------------------------------------------------------
	      dqmxzdxk4= (roa(i+1,j,k)*ua(i+1,j,k)*(mista(i+1,j,k)+(dt*k3mis)) - roa(i-1,j,k)*ua(i-1,j,k)*(mista(i-1,j,k)+&
			&(dt*k3mis)))/(2.d0*dx)
	      dqmyzdyk4= (roa(i,j+1,k)*va(i,j+1,k)*(mista(i,j+1,k)+(dt*k3mis)) - roa(i,j-1,k)*va(i,j-1,k)*(mista(i,j-1,k)+&
			&(dt*k3mis)))/(2.d0*dy)
	      dqmzzdzk4= (roa(i,j,k+1)*wa(i,j,k+1)*(mista(i,j,k+1)+(dt*k3mis)) - roa(i,j,k-1)*wa(i,j,k-1)*(mista(i,j,k-1)+&
			&(dt*k3mis)))/(2.d0*dz)

	      d2zdx2k4 = ((mista(i+1,j,k)+(dt*k3mis)) - 2.d0*(mista(i,j,k)+(dt*k3mis)) + (mista(i-1,j,k)+(dt*k3mis)))/(dx2)
	      d2zdy2k4 = ((mista(i,j+1,k)+(dt*k3mis)) - 2.d0*(mista(i,j,k)+(dt*k3mis)) + (mista(i,j-1,k)+(dt*k3mis)))/(dy2)
	      d2zdz2k4 = ((mista(i,j,k+1)+(dt*k3mis)) - 2.d0*(mista(i,j,k)+(dt*k3mis)) + (mista(i,j,k-1)+(dt*k3mis)))/(dz2)

	      k4mis  = -dqmxzdxk4-dqmyzdyk4-dqmzzdzk4 + (dife(i,j,k)/(re*sch))*(d2zdx2k4+d2zdy2k4+d2zdz2k4)
!--------------------------------fração de mistura----------------------------------------------------------------
	      mist(i,j,k) =mista(i,j,k) + (dt/6.d0)*(k1mis+2.d0*k2mis+2.d0*k3mis+k4mis)
	    end do
	  end do
	end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		   Condições de Contorno			!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	

!Entrada e Saída Externas
  do k=1,nk
      do j=1,nj
        u(1,j,k)	=.01d0!uin*ro(i,j,k)!*ro(i,j,k)*(1.d0-(y(j)/(D/2.d0))**2.d0)*1.5d0!1.d0!
        u(ni,j,k)	=0.d0!qmx(ni-1,j,k)
        v(1,j,k)	=0.d0
        v(ni,j,k)	=0.d0!qmy(ni-1,j,k)
	w(1,j,k)	=0.d0
	w(ni,j,k)	=0.d0!qmz(ni-1,j,k)
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
     		u(i,j,k)	=0.d0
     		v(i,j,k)	=0.d0
     		w(i,j,k)	=0.d0
		mist(i,j,k)     =0.d0
	        ro(i,j,k)       =1.d0
    	   enddo
!	Pardede Superior do jato
	   do j=nsup+2,nextsup
     		u(i,j,k)	=0.d0
     		v(i,j,k)	=0.d0
     		w(i,j,k)	=0.d0
		mist(i,j,k)     =0.d0
	        ro(i,j,k)       =1.d0
	   enddo
      enddo
!	Pardede Frontal do jato
      do j=nextinf,nextsup
     	   do k=nextinf,ninf-1
     		u(i,j,k)	=0.d0
     		v(i,j,k)	=0.d0
     		w(i,j,k)	=0.d0
		mist(i,j,k)     =0.d0
	        ro(i,j,k)       =1.d0
	   enddo
!	Parede Fundo do jato
	   do k=nsup+2,nextsup
     		u(i,j,k)	=0.d0
     		v(i,j,k)	=0.d0
     		w(i,j,k)	=0.d0
		mist(i,j,k)     =0.d0
	        ro(i,j,k)       =1.d0
	   enddo
      end do

! JATO DE COMBUSTÍVEL 
      do k=ninf-1,nsup+1
	do j=ninf-1,nsup+1
     		u(i,j,k)	=1.d0
     		v(i,j,k)	=0.d0
     		w(i,j,k)	=0.d0
		mist(i,j,k)     =1.d0
	        ro(i,j,k)       =1.d0
	end do
      enddo
end do
!Paredes do Duto
  do k=1,nk
    do i=1,ni
      	u(i,1,k)	=0.d0!.75d0*qmx(i,2,k)+.25d0*qmx(i,3,k)!0.d0!
        u(i,nj,k)	=0.d0!.75d0*qmx(i,nj-1,k)+.25d0*qmx(i,nj-2,k)!0.d0!
      	v(i,1,k)	=0.d0
        v(i,nj,k)	=0.d0
      	w(i,1,k)	=0.d0
        w(i,nj,k)	=0.d0
	mist(i,1,k)     =.75d0*mist(i,2,k)+.25d0*mist(i,3,k)
	mist(i,nj,k)    =.75d0*mist(i,nj-1,k)+.25d0*mist(i,nj-2,k)
    enddo
  enddo
  do j=1,nj
    do i=1,ni
      	u(i,j,1)	=0.d0!.75d0*qmx(i,j,2)+.25d0*qmx(i,j,3)!0.d0!
        u(i,j,nk)	=0.d0!.75d0*qmx(i,j,nk-1)+.25d0*qmx(i,j,nk-2)!0.d0!
      	v(i,j,1)	=0.d0
        v(i,j,nk)	=0.d0
      	w(i,j,1)	=0.d0
        w(i,j,nk)	=0.d0
	mist(i,j,1)     =.75d0*mist(i,j,2)+.25d0*mist(i,j,3)
	mist(i,j,nk)    =.75d0*mist(i,j,k-1)+.25d0*mist(i,j,k-2)
    enddo
  enddo
	
	do k=1,nk
	  do j=1,nj
	    do i=1,ni

	      if(mist(i,j,k).le.zst)then
	        T(i,j,k) = Tin + mist(i,j,k)*Tin + ((Q*yf1*mist(i,j,k))/(cp*1.d0*mwf))
	      else
	        T(i,j,k) = Tin + mist(i,j,k)*Tin + ((Q*yo22*(1.d0-mist(i,j,k)))/(cp*5.d0*mwo2))
	      end if

	    end do
	  end do
	end do


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
   	    if (mist(i,j,k)<=0.d0)then
                   mist(i,j,k)=1.d-20
     	    else    
                   if(mist(i,j,k)>=1.0000000001d0)then
                        mist(i,j,k)=1.d0
                   end if
	    end if
	      mista(i,j,k)= mist(i,j,k) 
	    if(mist(i,j,k).le.zst)then
	        T(i,j,k) = Tin + mist(i,j,k)*Tin + (Q*yf1*mist(i,j,k))/(cp*1.d0*mwf)	      
	    else
                T(i,j,k) = Tin + mist(i,j,k)*Tin +  (Q*yo22*(1.d0-mist(i,j,k)))/(cp*5.d0*mwo2)
	    end if

              pa(i,j,k)	  = p(i,j,k) 
	      qmxa(i,j,k) = qmx(i,j,k)
	      qmya(i,j,k) = qmy(i,j,k)
	      qmza(i,j,k) = qmz(i,j,k)
	      ua(i,j,k)   = u(i,j,k)
	      va(i,j,k)   = v(i,j,k)
	      wa(i,j,k)   = w(i,j,k)
	      Ta(i,j,k)   = T(i,j,k)
	      T(i,j,k)    = Ta(i,j,k)-(5.6697d-8/(cp*34.d0))*((Ta(i,j,k)-400.d0)**4.d0)*.021d0 ! Perda de Calor por Radiação

	      Tn(i,j,k)   = (T(i,j,k)-298.d0)/(1900.d0)
              ro(i,j,k)   = 1.d0/(Tn(i,j,k))
	      ro(i,j,k)   = 1.d0/(((alfa*Tn(i,j,k))/(1.d0 - alfa)) +1.d0)
	      roa(i,j,k)  = ro(i,j,k)
	    end do
	  end do
	end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	           Impressão dos dados				!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
!o valor dmax indica o erro de acordo com a conservação da massa d1

           if(mod(it,100).eq.0)then
             dmax=0.d0
	     jmax=0.d0
	     imax=0.d0
	     kmax=0.d0
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

           if(mod(it,1000).eq.0)then
	     open(11,file='mesh.dat',status='replace')
	     open(12,file='u.dat',status='replace')
	     open(13,file='v.dat',status='replace')
	     open(14,file='iso.dat',status='replace')
	     open(15,file='p.dat',status='replace')
	     open(16,file='mixture.dat',status='replace')
	     open(17,file='mixcenter.dat',status='replace')
	     open(18,file='velcenter.dat',status='replace')
	     open(19,file='temperature.dat',status='replace')
	     open(20,file='tempcenter.dat',status='replace')
	     open(21,file='w.dat',status='replace')
	     open(22,file='density.dat',status='replace')
	     open(23,file='viscosidade.dat',status='replace')


	     write(11,*)ni,nj,nk
	     write(12,*)ni-1,nj,1
	     write(13,*)ni-1,nj,1	  
	     write(14,*)ni-1,nj,1	  
	     write(15,*)ni-1,nj,1	  
	     write(16,*)ni-1,nj,1
	     write(17,*)ni,1
	     write(18,*)ni,1	  
	     write(19,*)ni-1,nj,1
	     write(20,*)ni,1
	     write(21,*)ni-1,nj,1
	     write(22,*)ni-1,nj,1
	     write(23,*)ni-1,nj,1

	       do j=1,nj
	         do i=1,ni
	           write(11,'(3e20.10)')x(i),y(j),z(nmk)
		 end do
	       end do	           

	       do j=1,nj
	         do i=1,ni-1
	           vtotal=((ua(i,j,nmk)**2.d0)+(va(i,j,nmk)**2.d0)+(wa(i,j,nmk)**2.d0))**(1.d0/3.d0)
	             write(12,'(3e20.11)')x(i),y(j),ua(i,j,nmk)
	             write(13,'(3e20.11)')x(i),y(j),va(i,j,nmk)
	             write(14,'(3e20.11)')x(i),y(j),vtotal
                     write(15,'(3e20.11)')x(i),y(j),pa(i,j,nmk)
	             write(21,'(3e20.11)')x(i),y(j),wa(i,j,nmk)
		     write(16,'(3e20.11)')x(i),y(j),mista(i,j,nmk)
		     write(19,'(3e20.11)')x(i),y(j),T(i,j,nmk)
	             write(22,'(3e20.11)')x(i),y(j),ro(i,j,nmk)
		     write(23,'(3e20.11)')x(i),y(j),mief(i,j,nmk)
	         end do
	       end do
	       
	       do i=1,ni
	  	 write(18,'(2e20.11)')x(i),u(i,nmj,nmk)
	         write(17,'(2e20.11)')x(i),mist(i,nmj,nmk)
	         write(20,'(2e20.11)')x(i),T(i,nmj,nmk)
	       end do
	      
	           close(11)
	           close(12)
		   close(13)
		   close(14)
		   close(15)
		   close(16)
		   close(17)
		   close(18)
	           close(19)
	           close(20)
	           close(21)
	           close(22)
	           close(23)
           end if
      end if

end do !término processo iterativo

call CPU_TIME(TIME)
print *, 'Execution Time:', TIME/60.d0, 'minutes'

end program jato3D

