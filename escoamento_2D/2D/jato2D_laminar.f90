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
	pin = 1.d0
	Tin = 300.d0

	  do j=1,nj
	    do i=1,ni
	    ro(i.j)   =1.d0
	    roa(i.j)  =ro(i.j)
	    u(i.j)    =uin!*ro(i,j)*(1.d0-(y(j)/(D/2.d0))**2.d0)*1.5d0
	    v(i.j)    =vin
	    p(i.j)    =pin
	    ua(i.j)   =u(i.j)
	    va(i.j)   =v(i.j)
	    pa(i.j)   =pin
	    qmx(i.j)  =ro(i.j)*u(i.j)
	    qmy(i.j)  =ro(i.j)*v(i.j)
	    qmxa(i.j) =qmx(i.j)
	    qmya(i.j) =qmy(i.j)
	    end do
	  end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	      Equação Quantidade de Movimento			!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
	
do it=1,itmax
  
	  do j=2,nj-1
	    do i=2,ni-1

	      dx      =.5d0*(x(i+1)-x(i-1)) 
	      dy      =.5d0*(y(j+1)-y(j-1)) 

	      dx2     = dx*dx
	      dy2     = dy*dy

	      duudx   = ((ua(i+1,j)*ua(i+1,j)*roa(i+1,j))-(ua(i-1,j)*ua(i-1,j)*roa(i-1,j)))/(2.d0*dx)
	      duvdy   = ((ua(i,j+1)*va(i,j+1)*roa(i,j+1))-(ua(i,j-1)*va(i,j-1)*roa(i,j-1)))/(2.d0*dy)

	      dvudx   = ((va(i+1,j)*ua(i+1,j)*roa(i+1,j))-(va(i-1,j)*ua(i-1,j)*roa(i-1,j)))/(2.d0*dx)
	      dvvdy   = ((va(i,j+1)*va(i,j+1)*roa(i,j+1))-(va(i,j-1)*va(i,j-1)*roa(i,j-1)))/(2.d0*dy)

	      d2udx2  = (ua(i+1,j)-2.d0*(ua(i,j))+ua(i-1,j))/(dx2)
	      d2udy2  = (ua(i,j+1)-2.d0*(ua(i,j))+ua(i,j-1))/(dy2)

	      d2vdx2  = (va(i+1,j)-2.d0*(va(i,j))+va(i-1,j))/(dx2)
	      d2vdy2  = (va(i,j+1)-2.d0*(va(i,j))+va(i,j-1))/(dy2)

	      dudx    = (roa(i+1,j)*ua(i+1,j)-roa(i-1,j)*ua(i-1,j))/(2.d0*dx)
	      dvdy    = (roa(i,j+1)*va(i,j+1)-roa(i,j-1)*va(i,j-1))/(2.d0*dy)

	      d1(i,j) = dudx + dvdy

	      dpdx    = (pa(i+1,j) - pa(i-1,j))/(2.d0*dx)
	      dpdy    = (pa(i,j+1) - pa(i,j-1))/(2.d0*dy)

	      k1x     = -duudx-duvdy- dpdx + (1/re)*(d2udx2 + d2udy2)
	      k1y     = -dvudx-dvvdy- dpdy + (1/re)*(d2vdx2 + d2vdy2)

!--------------------------------velocidades----------------------------------------------------------------------
	      u(i,j)    = ua(i,j) + dt*k1x
	      v(i,j)    = va(i,j) + dt*k1y
	    end do
	  end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		   Condições de Contorno			!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	

!Entrada e Saída Externas
  do k=1,nk
      do j=1,nj
        u(1,j)	=.01d0!uin*ro(i,j)!*ro(i,j)*(1.d0-(y(j)/(D/2.d0))**2.d0)*1.5d0!1.d0!
        u(ni,j)	=0.d0!qmx(ni-1,j)
        v(1,j)	=0.d0
        v(ni,j)	=0.d0!qmy(ni-1,j)
      enddo
   enddo

do i=1,ninj
! PAREDE INFERIOR DO INJETOR 
      do k=nextinf,nextsup	
!	Parede Inferior do Jato
     	   do j=nextinf,ninf-1
     		u(i,j)	=0.d0
     		v(i,j)	=0.d0
	        ro(i,j)       =1.d0
    	   enddo
!	Pardede Superior do jato
	   do j=nsup+2,nextsup
     		u(i,j)	=0.d0
     		v(i,j)	=0.d0
	        ro(i,j)       =1.d0
	   enddo
      enddo
!	Pardede Frontal do jato
      do j=nextinf,nextsup
     	   do k=nextinf,ninf-1
     		u(i,j)	=0.d0
     		v(i,j)	=0.d0
	        ro(i,j)       =1.d0
	   enddo
!	Parede Fundo do jato
	   do k=nsup+2,nextsup
     		u(i,j)	=0.d0
     		v(i,j)	=0.d0
	        ro(i,j)       =1.d0
	   enddo
      end do

! JATO DE COMBUSTÍVEL 
      do k=ninf-1,nsup+1
	do j=ninf-1,nsup+1
     		u(i,j)	=1.d0
     		v(i,j)	=0.d0
	        ro(i,j)       =1.d0
	end do
      enddo
end do
!Paredes do Duto
  do k=1,nk
    do i=1,ni
      	u(i,1)	=0.d0!.75d0*qmx(i,2)+.25d0*qmx(i,3)!0.d0!
        u(i,nj)	=0.d0!.75d0*qmx(i,nj-1)+.25d0*qmx(i,nj-2)!0.d0!
      	v(i,1)	=0.d0
        v(i,nj)	=0.d0
    enddo
  enddo
  do j=1,nj
    do i=1,ni
      	u(i,j,1)	=0.d0!.75d0*qmx(i,j,2)+.25d0*qmx(i,j,3)!0.d0!
        u(i,j,nk)	=0.d0!.75d0*qmx(i,j,nk-1)+.25d0*qmx(i,j,nk-2)!0.d0!
      	v(i,j,1)	=0.d0
        v(i,j,nk)	=0.d0
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

		 dx2		=dx*dx
		 dy2		=dy*dy


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

