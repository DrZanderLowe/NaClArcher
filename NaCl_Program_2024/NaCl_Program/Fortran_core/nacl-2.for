!      dhllnparameters are found in a file called 'param.dat'
c
             block data
            implicit double precision (a-h,o-z)
            double precision mr
	    common /refs/tr,pr,mr
            common /refs2/gextr,scorr,sions,sh2otr,gh2otr,sexmr
            data tr,pr,mr/298.15d0,.10d0,6.d0/
            data sh2otr/3.8806d0/
            end

      implicit double precision (a-h,o-z)
      common /dielc/d,dddp,d2ddp2,dddt,d2ddt2,rho,drhodp,drhodt
     1 			,d2rdp2,d2rdt2
      common /cpcalc/cpsave
      common /params/b0(24),b1(24),c0(24),c1(24),xv(10),xc(11)
      common /refs2/gextr,scorr,sions,sh2otr,gh2otr,sexmr
      common /refs/tr,pr,mr
      common /betao/bflag
      double precision num,nux,nu,ms,is,m,iss,mr, t, p
      character*3 answer
      data iac,iak/1,1/


	    open(21, file='param.dat')
! 
            do 1 i = 1,24
	      read(21,*)n,b0(i)
 1          continue
            do 2 i = 1,24
	      read(21,*)n,b1(i)
 2          continue
            do 3 i = 1,24
	      read(21,*)n,c0(i)
 3          continue
            do 13 i = 1,24
	      read(21,*)n,c1(i)
 13          continue
            do 4 i = 1,10
	      read(21,*)n,xv(i)
 4          continue
            do 5 i = 1,11
	      read(21,*)n,xc(i)
 5          continue

          sions=xc(11)

          lsav=4
             
!109      write(*,100)
!100      format('  what are T[K], p[MPa] and m[mol/kg] ? ')
          read *, t,p,m
!101      format(3d8.3)
					!print *, t, p, m
          !t=298.15D0
          !p=1D-1
          !m=1D0
          !print *, t, p, m
          !bflag=0.0
! check for negative inpSut values.
           if (m.lt.0.0d0)then
              print *
              print *,' a negative molality was entered.',
     1      '  please don''t do that.'
             stop 
           endif
           if (p.lt.0.0d0)then
              print *
              print *,' a negative pressure was entered.',
     1      '  please don''t do that.'
              stop 
           elseif (abs(p)<0.000001)then
              print *
              print *,' a zero pressure was entered.',
     1      '  please don''t do that. (pressure should be absolute)'
             stop 
           endif
           if (t.lt.0.0d0)then
              print *
              print *,' a negative temperature was entered.',
     1      '  please don''t do that.'
              stop 
           elseif (abs(t)<0.000001)then
              print *
              print *,' a zero temperature was entered.',
     1      '  please don''t do that. (temp. should be in k)'
             stop 
           endif
!  check the variables for the desired calculation. issue a warning
!  if outside of region of fit.

          iwarn=0
          if(p.gt.120.d0)iwarn=iwarn+3
          if(t.lt. 270.d0)then
              if(p.gt. 1.d0)iwarn=iwarn+1
              if((t .lt. 260.d0).and.(p.gt.0.5d0))iwarn=iwarn+1
          endif
          if((m.gt.7.d0) .and. (t.lt. 373.d0))iwarn=iwarn+1
          if((m.gt.8.5d0) .and. (t.lt.460.d0))iwarn=iwarn+1
          if((m.gt.10.d0) .and. (t.lt.560.d0))iwarn=iwarn+1
          if(m.gt.12.d0) iwarn=iwarn+1 
          if(t.gt. 370.d0)then
              if((m.gt.7.d0) .and. (p.gt.20.d0)) iwarn=iwarn+1
              if((m.gt.7.d0) .and. (p.gt.60.d0)) iwarn=iwarn+1
          endif 

          call dhlln(t,p,ap,ah,av,ac,ak,iac,iak)
            alphaw=-drhodt/rho
            betaw=drhodp/rho
            cpw=cpsave
            rhow=rho
            deltaw=alphaw*alphaw/rhow/cpw

!          bflag=1.0
          call gamcal(t,p,m,gam,phi)
          aw = exp(-phi*2.d0*m*18.0152d0/1000.d0)
!          bflag=0.0

					
          call lphi(t,p,m,hext)

					  
          call vphi(t,p,m,vphic,av)
             vs=(vphic*m+1.d3/rhow)/(1.d3+m*58.443d0)
             rhos=1.d0/vs
  		
           call kphi(t,p,m,comphi,ak)
              comsol=(m*comphi+1.d3*betaw/rhow)/(1.d3+m*58.443d0)
              comsol=comsol*rhos
  
          call cpcalc3(t,p,m,cpphi,ac)
               cpsol=(cpphi*m+1.d3*cpw)/(1.d3+m*58.443d0)

          call ephi(t,p,m,ephic,av)
          	 
             esol=(ephic*m+alphaw*1.d3/rhow)/(1.d3+m*58.443d0)
             esol=esol*rhos
!    calculate difference in the two compressibilities
            deltsl=t*esol*esol/rhos/cpsol
             betais=comsol-deltsl
!             print *,' isentropi! compressibility = ',betais
             corris=deltsl*(1.d0+alphaw/esol)*ephic/esol-deltaw*cpphi/
     1               (cpsol*rhos)+(deltaw-(deltsl*alphaw/esol))*vphic

!              print *,' kphi - kphi(s) = ',corris
!              print *,' kphi(s) = ',comphi-corris


!    do first pass through to get reference values.
             call wsteam(tr,rho,pr,si,dpr,ddt,u,h,s,cv,cp,a,4)
             scorr=sh2otr-s
             gh2otr=h-tr*sh2otr      
             call gamcal(tr,pr,mr,gamma,phic)
             gextr=2.d0*mr*8.3144d-3*tr*(log(gamma)+1.d0-phic)
             call lphi(tr,pr,mr,hexmr)
             sexmr=(hexmr-gextr/mr)/tr

      
                    call delg20(t,p,delg)   
                    call s20(t,p,st)

                    delh=delg+t*st-298.15d0*sions
                    !print *, delg, delh, st, t ,sions
                    dels=st-sions  
              

! conditions, warnings and activities to the screen
       write(6,702)t,p,m
 702      format(1x'variables: '5x,f8.3,' K',5x,f8.3,' MPa',
     1              5x,f9.5,' mol/kg')
	     if (iwarn .gt. 0)write(6,703)iwarn
 703      format(1x,'!!!  warning condition ',i2,' exists.'
     1   '  the larger the warning value the further away from data'
     1   ' the calculation is.')  
	     write(6,800)gam
 800      format(/1x,'solute activity coefficient',10x,f6.4)
	     write(6,801)phi
 801      format(1x,'solvent osmoti! coefficient',10x,f6.4)
	     write(6,802)aw
 802      format(1x,'solvent activity',21x,f8.6)

! apparent molar property results to the screen.
	     write(6,900)
 900      format(/27x,'apparent molar properties')
	     write(6,901)vphic
 901      format(1x,'volume [cm**3/mol]',12x,f9.2)
	     write(6,902)ephic
 902      format(1x,'expansivity [cm**3/mol*K]',7x,1p,e12.4)
	     write(6,903)comphi
 903      format(1x,'compressibility [cm**3/mol*K]',3x,1p,e12.4)
	     write(6,904)hext
 904      format(1x,'relative enthalpy [kJ/mol]',7x,1p,e12.5)
	     write(6,905)cpphi
 905      format(1x,'cp [J/K*mol]',20x,1p,e12.4)


! property results to the screen.
	     write(6,1001)
 1001     format(/27x,'solution property',10x,'solvent property')
	     write(6,1002)rhos,rhow
 1002     format(1x,'density [g/cm**-3]',12x,f8.6,19x,f8.6)
	     write(6,1003)esol,alphaw
 1003     format(1x,'expansivity [1/K]',11x,e12.4,15x,e12.4)
	     write(6,1004)comsol,betaw
 1004     format(1x,'compressibility [1/MPa]',5x,e12.4,15x,e12.4)
	     write(6,1005)cpsol,cpw
 1005     format(1x,'cp [J/K*g]',18x,1p,e12.4,15x,e12.4)
 
! standard-state values to the screen.
	     write(6,1101)
 1101     format(/10x,'g0 - g0(tr, pr)',10x,'h0 - h0(tr, pr)',10x,
     1    's0 - s0(tr, pr)')
	     write(6,1102)
 1102     format(10x,'[kJ/mol-1]'16x,'[kJ/mol]',16x,'[J/K*mol]')
	     write(6,1103)delg,delh,dels*1.d3
 1103	  format(10x,e13.7,11x,e13.7,11x,e13.7)

             write(6,1111)d
 1111     format(/1x,'d',14x,e12.6)
             write(6,1112)dddp
 1112     format(1x,'dddp',11x,e12.6)
             write(6,1113)d2ddp2
 1113     format(1x,'d2ddp2',9x,e12.6)
             write(6,1114)dddt
 1114     format(1x,'dddt',11x,e12.6)
             write(6,1115)d2ddt2
 1115     format(1x,'d2ddt2',9x,e12.6)

             write(6,1120)rho
 1120     format(/1x,'rho',12x,e12.6)
             write(6,1121)drhodp
 1121     format(1x,'drhodp',9x,e12.6)
             write(6,1122)drhodt
 1122     format(1x,'drhodt',9x,e12.6)
             write(6,1123)d2rdp2
 1123     format(1x,'d2rdp2',9x,e12.6)
             write(6,1124)d2rdt2
 1124     format(1x,'d2rdt2',9x,e12.6)

             write(6,1130)ap
 1130     format(/1x,'ap',13x,e12.6)
             write(6,1131)ah
 1131     format(1x,'ah',13x,e12.6)
             write(6,1132)av
 1132     format(1x,'av',13x,e12.6)
             write(6,1133)ac
 1133     format(1x,'ac',13x,e12.6)
             write(6,1134)ak
 1134     format(1x,'ak',13x,e12.6,/10x)

!1200     write(*,102)
!102      format(' another calculation (y) ? ')
!         read(5,103)answer
!103      format(a)
!         if(answer .eq. 'y')goto 109
          stop 
          end        


      subroutine cpcalc3(t,p,m,cpphi,ac)
      implicit double precision (a-h,o-z)
      common /dielc/d,dddp,d2ddp2,dddt,d2ddt2,rho,drhodp,drhodt
     1 ,d2rdp2,d2rdt2
      common /params/b0(24),b1(24),c0(24),c1(24),xv(10),xc(11)
      common /cpcalc/cpsave

      double precision num,nux,nu,ms,is,m,iss,lot,lot2,lot3,lot4
     1 ,lota  
     
     	double precision, dimension (:,:), allocatable :: xb
      double precision, dimension (:), allocatable :: xva, xca
      allocate (xb(4,24), xca(11), xva(10))

       data r,num,nux,nu,zm,zx/8.3144d-3,1.d0,1.d0,2.d0,1.d0,-1.d0/
       data ps,ts,ms/.10d0,298.15d0,6.d0/
       data alpha,alpha2,bb/2.0d0,2.5d0,1.2d0/
          ajt=ac*r
         cp=cpsave

      iss=0.5d0*ms*(num*zm*zm + nux*zx*zx)
      is=(m*(num*zm*zm+nux*zx*zx))/2.d0
      sqis=sqrt(is)
      ai=alpha*sqis
      ai2=alpha2*sqis
       sqiss=sqrt(iss)
      ais=alpha*sqiss
      a2is=alpha2*sqiss 
      if(m.le.0.0d0)then
         ai=0.5d0
         ai2=0.5d0
      endif
      cor=(1.d0-(1.d0+ai)*exp(-ai))/(ai*ai)
      cor2=(6.d0-(6.d0+ai2*(6.d0+ai2*(3.d0+ai2)))*exp(-ai2))/
     1      ((ai2*ai2)*(ai2*ai2))
      corr=(1.d0-(1.d0+ais)*exp(-ais))/(ais*ais)
      corr2=(6.d0-(6.d0+a2is*(6.d0+a2is*(3.d0+a2is)))*exp(-a2is))/
     1      ((a2is*a2is)*(a2is*a2is))

      xca(1)=1.d0
      xca(2)=300.d0/t/t
      xca(3)=t/300.d0
      xca(4)=t*t/300.d0/300.d0
      xca(6)=t*t*t/3.d2/3.d2/3.d2
      xca(7)=1.d2/t
      xca(8)=0.d0
      xca(9)=0.d0
      xca(10)=0.d0
      xca(11)=0.d0


      do 5 i=1,10
      	xva(i)=0.d0
 5		continue
 
      dp=p/100.d0
      dt=t/300.d0
      xva(4)=-t*6.d0*dt*(p-ps)/9.d4/1.d3
      xva(7)=-t*2.d0*(p*p-ps*ps)/9.d4/2.d0/1.d3/100.d0
      xva(10)=0.d0
      

      xb(1,1)=-2.d0*num*nux*r*t*t
      xb(2,1)=xb(1,1)*(2.d0*m*cor-2.d0*ms*corr)
      xb(3,1)=xb(1,1)*(m*m-ms*ms)*num*zm
      xb(4,1)=xb(1,1)*num*zm*4.d0*(m*m*cor2-ms*ms*corr2)
      xb(1,1)=xb(1,1)*(m-ms)

        tt=t/500.d0
        lot=t-225.d0
        lota=t-200.d0
        hit=650.d0-t
        hit2=hit*hit
        hit4=hit2*hit2
        lot2=lot*lot
        lot3=lot2*lot
         lot4=lot3*lot
! next set of lines give dbeta/dt.
      do 1 j = 1,4
         xb(j,2)=xb(j,1)/1.d+3
         xb(j,3)=xb(j,1)*tt*2.d0/500.d0
	 			 xb(j,4)=-xb(j,1)/lota/lota
         xb(j,5)=-3.d0*xb(j,1)*1.d+4/lota/lota/lota/lota
         xb(j,6)=-2.d0*xb(j,1)*1.d2/lota/lota/lota
         xb(j,7)=-2.d0*xb(j,1)*2.d2/t/t/t
         xb(j,8)=3.d0*xb(j,1)*tt*tt/500.d0
         xb(j,9)=0.5d0*xb(j,1)/(hit**1.5d0)
         xb(j,10)=0.0d0
         xb(j,11)=-xb(j,1)*p*1.d-6*2.d2/lot2
         xb(j,12)=3.d0*xb(j,1)*p*1.d2/hit4
         xb(j,13)=xb(j,1)*1.d-5*p/500.d0
         xb(j,14)=xb(j,1)*p*2.d2*1.d-6/hit2
         xb(j,15)=0.0d0
         xb(j,16)=-xb(j,1)*p*p*1.d-8*2.d2/lot2
         xb(j,17)=3.d0*xb(j,1)*p*p/hit4
         xb(j,18)=xb(j,1)*p*p*1.d-7/500.d0
         xb(j,19)=2.d0*xb(j,1)*p*p*tt*1.d-7/500.d0
         xb(j,20)=-2.d0*xb(j,1)*1.d-6*4.d4*p/lot/lot/lot
         xb(j,21)=2.d0*xb(j,1)*1.d-5*p*t/500.d0/5.d2
         xb(j,22)=-xb(j,1)*p*p*p*1.d-10*2.d2/lot2
         xb(j,23)=xb(j,1)*p*p*p*3.d0*1.d-2/hit4
         xb(j,24)=3.d0*xb(j,1)*2.d2/hit4
         xb(j,24)=3.d0*xb(j,1)*2.d2/hit4
  1   continue

!  next loop multiplies dbeta/dt by 2/t.
      do 2 j = 1, 4
        do jj = 2, 24
     		xb(j,jj)=xb(j,jj)*2.d0/t
     		end do
  2		continue

! next loop adds d2beta/dt2.
           hit5=hit4*hit
           hit3=hit2*hit
      do 3 j = 1, 4
         xb(j,2)=xb(j,2)+0.d0 
         xb(j,3)=xb(j,3)+xb(j,1)*2.d0/500.d0/500.d0
         xb(j,4)=xb(j,4)+xb(j,1)*2.d0/lota/lota/lota
         xb(j,5)=xb(j,5)+12.d0*xb(j,1)*1.d+4/lota/lota/lota/lota/lota
         xb(j,6)=xb(j,6)+6.d0*xb(j,1)*1.d2/lota/lota/lota/lota
         xb(j,7)=xb(j,7)+xb(j,1)*6.d0*2.d2/t/t/t/t
         xb(j,8)=xb(j,8)+6.d0*xb(j,1)*tt/500.d0/500.d0
         xb(j,9)=xb(j,9)+0.75d0*xb(j,1)/(hit**2.5d0)
         xb(j,10)=xb(j,10)+0.0d0
         xb(j,11)=xb(j,11)+2.d0*xb(j,1)*200.d-06*p/lot3
         xb(j,12)=xb(j,12)+12.d0*xb(j,1)*1.d2*p/hit5
         xb(j,13)=xb(j,13)+0.0d0
         xb(j,14)=xb(j,14)+2.d0*xb(j,1)*200.d-6*p/hit3
         xb(j,15)=0.0d0
         xb(j,16)=xb(j,16)+2.0d0*xb(j,1)*p*p*200.d-8/lot3
         xb(j,17)=xb(j,17)+12.d0*xb(j,1)*p*p/hit5
         xb(j,18)=xb(j,18)+0.0d0
         xb(j,19)=xb(j,19)+2.d0*xb(j,1)*p*p*1.d-7/500.d0/500.d0
         xb(j,20)=xb(j,20)+6.d0*xb(j,1)*1.d-6*4.d4*p/lot/lot/lot/lot
         xb(j,21)=xb(j,21)+2.d0*xb(j,1)*1.d-5*p/500.d0/5.d2
         xb(j,22)=xb(j,22)+2.d0*xb(j,1)*p*p*p*200.d-10/lot3
         xb(j,23)=xb(j,23)+xb(j,1)*12.d0*p*p*p*1.d-2/hit5
         xb(j,24)=xb(j,24)+12.d0*xb(j,1)*2.d2/hit5
         xb(j,1)=0.d0
 3    continue  


         beta0v=0.0d0
         beta1v=0.0d0
         cphiv=0.0d0
         cphiv1=0.0d0
         dcp0dp=0.0d0
         cpadd=0.0d0
         do 20 i=1, 24
           beta0v=beta0v+xb(1,i)*b0(i)
 20      continue
         do 30 i=1, 24
           beta1v=beta1v+xb(2,i)*b1(i)
 30      continue
         do 40 i=1, 24
           cphiv=cphiv+xb(3,i)*c0(i)
 40      continue
         do 140 i=1, 24
           cphiv1=cphiv1+xb(4,i)*c1(i)
 140      continue
         do 50 i=1,10
           dcp0dp=dcp0dp+xv(i)*xva(i)
 50      continue
         do 60 i=1,7
           cpadd=cpadd+xc(i)*xca(i)
 60      continue
      dhl=abs(zm*zx)*(num+nux)*ajt*(log((1.d0+bb*sqis)/
     1     (1.d0+bb*sqiss)))/(2.d0*bb)
      dhl=dhl+beta0v+beta1v+cphiv+cphiv1+dcp0dp+cpadd

       cpphi=-1.d0*cp/ms+dhl
       cpphi=cpphi*1.d3
			 
			 deallocate (xb, xca, xva)
       return
      end



      subroutine ephi(t,p,m,ephic,av)
      implicit double precision (a-h,o-z)
      common /dielc/d,dddp,d2ddp2,dddt,d2ddt2,rho,drhodp,drhodt
     1 ,d2rdp2,d2rdt2
      common /params/b0(24),b1(24),c0(24),c1(24),xv(10),xc(11)
      common /cpcalc/cpsave
      double precision num,nux,nu,ms,is,m,iss,lot
      
      double precision, dimension (:,:), allocatable :: xb
      double precision, dimension (:), allocatable :: xva, xca
      allocate (xb(4,24), xca(11), xva(10))


       data r,num,nux,nu,zm,zx/8.3144d-3,1.d0,1.d0,2.d0,1.d0,-1.d0/
       data ps,ts,ms/.1d0,298.15d0,6.d0/
       data alpha,alpha2,bb/2.d0,2.5d0,1.2d0/

!  make ae by means of numerical differentiation.
          call dhlln(t-0.01d0,p,ap,ah,av,ac,ak,iac,iak)
          drdpl=drhodp
          avl=av
          call dhlln(t+0.01d0,p,ap,ah,av,ac,ak,iac,iak)
          drdpdt=(drhodp-drdpl)/0.02d0
          avt=(av-avl)/0.02d0
          call dhlln(t,p,ap,ah,av,ac,ak,iac,iak)
 

      iss=0.5d0*ms*(num*zm*zm + nux*zx*zx)
      is=(m*(num*zm*zm+nux*zx*zx))/2.d0
      sqis=sqrt(is)
      ai=alpha*sqis
      ai2=alpha2*sqis
       sqiss=sqrt(iss)
      ais=alpha*sqiss
      a2is=alpha2*sqiss 
      if(abs(ai)<0.000001)then
         ai=0.5
         ai2=0.5
      endif
      cor=(1.d0-(1.d0+ai)*exp(-ai))/(ai*ai)
      corr=(1.d0-(1.d0+ais)*exp(-ais))/(ais*ais)

      dp=p/100.d0
      dt=t/300.d0
      xva(1)=0.0d0
      xva(2)=1.d-2/300.d0
      xva(3)=2.d0*dt/1.d3/300.d0
      xva(4)=3.d0*dt*dt/1.d3/3.d2
      xva(5)=0.0d0
      xva(6)=dp/3.d2/1.d1
      xva(7)=2.d0*dp*dt/1.d3/3.d2
      xva(8)=3.d0*dt*dt*dp/1.d4/3.d2
      xva(9)=dp*dp/1.d3/3.d2
      xva(10)=0.0d0
      xva(3)=(p+10.d0)**1.5d0/3.d2/1.d7
        do 1000 ix=1,10
        	xva(ix)=xva(ix)*xv(ix)
 1000		continue
c*****

c*****
      xb(1,1)=2.d0*num*nux*r*t
      xb(2,1)=xb(1,1)*(2.d0*m*cor-2.d0*ms*corr)
      xb(3,1)=xb(1,1)*(m*m-ms*ms)
      xb(4,1)=xb(1,1)*(m*m*(6.d0+(-ai2*ai2*ai2-3.d0*ai2*ai2
     1         -6.d0*ai2-6.d0)*exp(-ai2))/(ai2**4.d0)-
     2             ms*ms*(6.d0+(-a2is*a2is*a2is-3.d0*a2is*a2is
     3         -6.d0*a2is-6.d0)*exp(-a2is))/(a2is**4.d0))
      xb(1,1)=xb(1,1)*(m-ms)

        tt=t/500.d0
        lot=t-225.d0
        hit=650.d0-t
 


       do 1 j = 1,4
         xb(j,2)=0.d0
         xb(j,3)=0.d0
         xb(j,4)=0.d0
         xb(j,5)=0.d0
	 xb(j,6)=0.d0
         xb(j,7)=0.d0
         xb(j,8)=0.d0
         xb(j,9)=0.d0
         xb(j,10)=xb(j,1)*1.d-5
         xb(j,11)=xb(j,1)*2.d2*1.d-6/lot
         xb(j,12)=xb(j,1)*1.d2/hit/hit/hit
         xb(j,13)=xb(j,1)*tt*1.d-5
         xb(j,14)=xb(j,1)*2.d2*1.d-6/hit
         xb(j,15)=2.d0*xb(j,1)*p*1.d-7
         xb(j,16)=2.d0*xb(j,1)*p*2.d2*1.d-8/lot
         xb(j,17)=2.d0*xb(j,1)*p/hit/hit/hit
         xb(j,18)=2.d0*xb(j,1)*p*tt*1.d-7
         xb(j,19)=2.d0*xb(j,1)*p*tt*tt*1.d-7
         xb(j,20)=xb(j,1)*1.d-6*4.d4/lot/lot
         xb(j,21)=xb(j,1)*1.d-5*t*t/500.d0/5.d2
         xb(j,22)=3.d0*xb(j,1)*p*p*1.d-10*200.d0/lot
         xb(j,23)=3.d0*xb(j,1)*p*p*1.d-2/hit/hit/hit
         xb(j,24)=0.0d0
  1   continue
          do 2 j=1,4
            do jj=2,24
               xb(j,jj)=xb(j,jj)/t
            end do
 2     continue
       do 3 j = 1,4
         xb(j,10)=xb(j,10)+0.0d0
         xb(j,11)=xb(j,11)-xb(j,1)*2.d2*1.d-6/lot/lot
         xb(j,12)=xb(j,12)+3.d0*xb(j,1)*1.d2/hit/hit/hit/hit
         xb(j,13)=xb(j,13)+xb(j,1)*1.d-5/5.d2
         xb(j,14)=xb(j,14)+xb(j,1)*2.d2*1.d-6/hit/hit
         xb(j,15)=xb(j,15)+0.0d0
         xb(j,16)=xb(j,16)-2.d0*xb(j,1)*p*2.d2*1.d-8/lot/lot
         xb(j,17)=xb(j,17)+6.d0*xb(j,1)*p/hit/hit/hit/hit
         xb(j,18)=xb(j,18)+2.d0*xb(j,1)*p*1.d-7/5.d2
         xb(j,19)=xb(j,19)+4.d0*xb(j,1)*p*tt*1.d-7/5.d2
         xb(j,20)=xb(j,20)-2.d0*xb(j,1)*1.d-6*4.d4/lot/lot/lot
         xb(j,21)=xb(j,21)+2.d0*xb(j,1)*1.d-5*t/500.d0/5.d2
         xb(j,22)=xb(j,22)-3.d0*xb(j,1)*p*p*1.d-10*200.d0/lot/lot
         xb(j,23)=xb(j,23)+9.d0*xb(j,1)*p*p*1.d-2/hit/hit/hit/hit
         xb(j,24)=xb(j,24)+0.0d0
!         xb(j,1)=0.d0
  3   continue



         beta0v=0.d0
         beta1v=0.d0
         cphiv1=0.d0
         cphiv2=0.d0
         vms=0.d0


         do 10 ix=10,24
          beta0v=beta0v+xb(1,ix)*b0(ix)
          beta1v=beta1v+xb(2,ix)*b1(ix)
          cphiv1=cphiv1+xb(3,ix)*c0(ix)
          cphiv2=cphiv2+xb(4,ix)*c1(ix)
 10      continue
          do 202 ix=1,10
            vms=vms+xva(ix)
 202       continue

      dhl=abs(zm*zx)*(num+nux)*avt*(log(1.d0+bb*sqis)
     1   -log(1.d0+bb*sqiss))/(2.d0*bb)/1.d3

       dhl=dhl+beta0v+beta1v+cphiv1+cphiv2+vms
       ephic=(1.0d0*drhodt/(ms*rho*rho)+dhl)*1.d3
			 deallocate (xb, xca, xva)
			 
        return
        end

      subroutine kphi(t,p,m,comphi,ak)
      implicit double precision (a-h,o-z)
      common /dielc/d,dddp,d2ddp2,dddt,d2ddt2,rho,drhodp,drhodt
     1 ,d2rdp2,d2rdt2
      common /params/b0(24),b1(24),c0(24),c1(24),xv(10),xc(11)
      common /cpcalc/cpsave
      double precision num,nux,nu,ms,is,m,iss,lot
      
      
      double precision, dimension (:,:), allocatable :: xb
      double precision, dimension (:), allocatable :: xva, xca
      allocate (xb(4,24), xca(11), xva(10))


       data r,num,nux,nu,zm,zx/8.3144d-3,1.d0,1.d0,2.d0,1.d0,-1.d0/
       data ps,ts,ms/.1d0,298.15d0,6.d0/
       data alpha,alpha2,bb/2.d0,2.5d0,1.2d0/

          avt=ak

      iss=0.5d0*ms*(num*zm*zm + nux*zx*zx)
      is=(m*(num*zm*zm+nux*zx*zx))/2.d0
      sqis=sqrt(is)
      ai=alpha*sqis
      ai2=alpha2*sqis
      sqiss=sqrt(iss)
      ais=alpha*sqiss
      a2is=alpha2*sqiss 
      if(ai.le.0.d0)then
         ai=0.5
         ai2=0.5
      endif 
      cor=(1.d0-(1.d0+ai)*exp(-ai))/(ai*ai)
      corr=(1.d0-(1.d0+ais)*exp(-ais))/(ais*ais)


      dp=p/1.d2
      dt=t/3.d2
      xva(1)=0.0d0
      xva(2)=0.0d0
      xva(3)=0.0d0
      xva(4)=0.0d0
      xva(5)=-1.d0/1.d3/1.d2
      xva(6)=-dt/1.d2/1.d1
      xva(7)=-dt*dt/1.d3/1.d2
      xva(8)=-dt*dt*dt/1.d4/1.d2
      xva(9)=-2.d0*dp*dt/1.d3/1.d2
      xva(10)=-1.5d0*sqrt(p+10.d0)/1.d6
      xva(3)=-1.5d0*sqrt(p+10.d0)*dt/1.d7

c*****
      xb(1,1)=2.d0*num*nux*r*t
      xb(2,1)=xb(1,1)*(2.d0*m*cor-2.d0*ms*corr)
      xb(3,1)=xb(1,1)*(m*m-ms*ms)
      xb(4,1)=xb(1,1)*(m*m*(6.d0+(-ai2*ai2*ai2-3.d0*ai2*ai2
     1         -6.d0*ai2-6.d0)*exp(-ai2))/(ai2**4.d0)-
     2         ms*ms*(6.d0+(-a2is*a2is*a2is-3.d0*a2is*a2is
     3         -6.d0*a2is-6.d0)*exp(-a2is))/(a2is**4.d0))
      xb(1,1)=xb(1,1)*(m-ms)

        tt=t/500.d0
        lot=t-225.d0
        hit=650.d0-t
 
       do 1 j = 1,4
         xb(j,2)=0.d0
         xb(j,3)=0.d0
         xb(j,4)=0.d0
         xb(j,5)=0.d0
	 			 xb(j,6)=0.d0
         xb(j,7)=0.d0
         xb(j,8)=0.d0
         xb(j,9)=0.d0
         xb(j,10)=0.d0
         xb(j,11)=0.d0
         xb(j,12)=0.d0
         xb(j,13)=0.d0
         xb(j,14)=0.d0
         xb(j,15)=-2.d0*xb(j,1)*1.d-7
         xb(j,16)=-2.d0*xb(j,1)*2.d2*1.d-8/lot
         xb(j,17)=-2.d0*xb(j,1)/hit/hit/hit
         xb(j,18)=-2.d0*xb(j,1)*tt*1.d-7
         xb(j,19)=-2.d0*xb(j,1)*tt*tt*1.d-7
         xb(j,20)=0.0d0
         xb(j,21)=0.0d0
         xb(j,22)=-6.d0*xb(j,1)*p*1.d-10*200.d0/lot
         xb(j,23)=-6.d0*xb(j,1)*1.d-2*p/hit/hit/hit
         xb(j,24)=0.0d0
         xb(j,1)=0.d0
  1   continue

         beta0v=0.d0
         beta1v=0.d0
         cphiv1=0.d0
         cphiv2=0.d0
         vms=0.d0


         do 10 ix=15,24
          beta0v=beta0v+xb(1,ix)*b0(ix)
          beta1v=beta1v+xb(2,ix)*b1(ix)
          cphiv1=cphiv1+xb(3,ix)*c0(ix)
          cphiv2=cphiv2+xb(4,ix)*c1(ix)
 10      continue
          do 202 ix=1,10
            vms=vms+xva(ix)*xv(ix)
 202       continue

      dhl=abs(zm*zx)*(num+nux)*avt*(log(1.d0+bb*sqis)
     1 -log(1.d0+bb*sqiss))/(2.d0*bb)/1.d3
       dhl=-dhl+beta0v+beta1v+cphiv1+cphiv2+vms

      yy=-1.d0*drhodp/(ms*rho*rho)+dhl
       comphi=yy*1000.d0
			deallocate (xb, xca, xva)
        return
         end



      subroutine vphi(t,p,m,vphic,avt)
      implicit double precision (a-h,o-z)
      common /dielc/d,dddp,d2ddp2,dddt,d2ddt2,rho,drhodp,drhodt
     1 ,d2rdp2,d2rdt2
      common /params/b0(24),b1(24),c0(24),c1(24),xv(10),xc(11)
      common /cpcalc/cpsave
      double precision num,nux,nu,ms,is,m,iss,lot
      
      double precision, dimension (:,:), allocatable :: xb
      double precision, dimension (:), allocatable :: xva, xca
      allocate (xb(4,24), xca(11), xva(10))
      
      data r,num,nux,nu,zm,zx/8.3144d-3,1.d0,1.d0,2.d0,1.d0,-1.d0/
      data ps,ts,ms/.10d0,298.15d0,6.d0/
      data alpha,alpha2,bb/2.d0,2.5d0,1.2d0/

      iss=0.5d0*ms*(num*zm*zm + nux*zx*zx)
      is=(m*(num*zm*zm+nux*zx*zx))/2.d0
      sqis=sqrt(is)
      ai=alpha*sqis
      ai2=alpha2*sqis
       sqiss=sqrt(iss)
      ais=alpha*sqiss
      a2is=alpha2*sqiss 
      if(ai.le.0.d0)then
         ai=0.5
         ai2=0.5
      endif
      cor=(1.d0-(1.d0+ai)*exp(-ai))/(ai*ai)
      corr=(1.d0-(1.d0+ais)*exp(-ais))/(ais*ais)

      dp=p/100.d0
      dt=t/300.d0

      xva(1)=1.d-1
      xva(2)=dt/1.d2
      xva(3)=dt*dt/1.d3
      xva(4)=dt*dt*dt/1.d3
      xva(5)=dp/1.d3
      xva(6)=dp*dt/1.d1
      xva(7)=dp*dt*dt/1.d3
      xva(8)=dt*dt*dt*dp/1.d4
      xva(9)=dp*dp*dt/1.d3
      xva(10)=(p+10.d0)**1.5d0/1.d6
      xva(3)=(p+10.d0)**1.5d0*dt/1.d7

c*****
      xb(1,1)=2.d0*num*nux*r*t
      xb(2,1)=xb(1,1)*(2.d0*m*cor-2.d0*ms*corr)
      xb(3,1)=xb(1,1)*(m*m-ms*ms)
      xb(4,1)=xb(1,1)*(m*m*(6.d0+(-ai2*ai2*ai2-3.d0*ai2*ai2
     1         -6.d0*ai2-6.d0)*exp(-ai2))/(ai2**4.d0)-
     2             ms*ms*(6.d0+(-a2is*a2is*a2is-3.d0*a2is*a2is
     3         -6.d0*a2is-6.d0)*exp(-a2is))/(a2is**4.d0))
      xb(1,1)=xb(1,1)*(m-ms)


      tt=t/500.d0
      lot=t-225.d0
      hit=650.d0-t
 
      do 1 j = 1,4
         xb(j,2)=0.d0
         xb(j,3)=0.d0
         xb(j,4)=0.d0
         xb(j,5)=0.d0
	 			 xb(j,6)=0.d0
         xb(j,7)=0.d0
         xb(j,8)=0.0d0
         xb(j,9)=0.0d0
         xb(j,10)=xb(j,1)*1.d-5
         xb(j,11)=xb(j,1)*2.d2*1.d-6/lot
         xb(j,12)=xb(j,1)*1.d2/hit/hit/hit
         xb(j,13)=xb(j,1)*tt*1.d-5
         xb(j,14)=xb(j,1)*2.d2*1.d-6/hit
         xb(j,15)=2.d0*xb(j,1)*p*1.d-7
         xb(j,16)=2.d0*xb(j,1)*p*2.d2*1.d-8/lot
         xb(j,17)=2.d0*xb(j,1)*p/hit/hit/hit
         xb(j,18)=2.d0*xb(j,1)*p*tt*1.d-7
         xb(j,19)=2.d0*xb(j,1)*p*tt*tt*1.d-7
         xb(j,20)=xb(j,1)*1.d-6*4.d4/lot/lot
         xb(j,21)=xb(j,1)*1.d-5*t*t/500.d0/5.d2
         xb(j,22)=3.d0*xb(j,1)*p*p*1.d-10*200.d0/lot
         xb(j,23)=3.d0*xb(j,1)*p*p*1.d-2/hit/hit/hit
         xb(j,24)=0.0d0
         xb(j,1)=0.d0
  1   continue

         beta0v=0.d0
         beta1v=0.d0
         cphiv1=0.d0
         cphiv2=0.d0
         vms=0.d0


         do 10 ix=10,24
          beta0v=beta0v+xb(1,ix)*b0(ix)
          beta1v=beta1v+xb(2,ix)*b1(ix)
          cphiv1=cphiv1+xb(3,ix)*c0(ix)
          cphiv2=cphiv2+xb(4,ix)*c1(ix)
 10      continue
          do 202 ix=1,10
            vms=vms+xv(ix)*xva(ix)
 202       continue

      dhl=abs(zm*zx)*(num+nux)*avt*(log(1.d0+bb*sqis)
     1   -log(1.d0+bb*sqiss))/(2.d0*bb)/1.d3
      dhl=dhl+beta0v+beta1v+cphiv1+cphiv2+vms
      yy=-1.d0/(ms*rho)+dhl
      vphic=yy*1000.d0
      
      deallocate (xb, xca, xva)  
      return
      end


      subroutine gamcal(t,p,m,gam,phi)
      implicit double precision (a-h,o-z)
      common /params/ b0(24),b1(24),c0(24),c1(24),xv(10),xc(11)
      common /dielc/d,dddp,d2ddp2,dddt,d2ddt2,rho,drhodp,drhodt
     1 ,d2rdp2,d2rdt2
      common /cpcalc/cpsave
      common /betao/bflag
      double precision num,nux,nu,ms,is,m,iss,lot,lota
      
      double precision, dimension (:,:), allocatable :: xb
      allocate (xb(4,24))


       data r,num,nux,nu,zm,zx/8.3144d-3,1.d0,1.d0,2.d0,1.d0,-1.d0/
       data ps,ts,ms/.1d0,298.15d0,6.d0/
       data iac,iak/0,0/
       data alpha,bb/2.0d0,1.2d0/

         if(m.le.0.0d0)then
               gam=1.d0
               phi=1.0d0
               return
         endif

c**************************************************************************
           call dhlln(t,p,ap,ah,av,ac,ak,iac,iak)
           aphit=ap
           is=0.5d0*m*(num*zm*zm + nux*zx*zx)
           ai2=2.5d0*sqrt(m)           
           ais=2.0d0*sqrt(m)
           ai2trm=6.d0-(6.d0+ai2*(6.d0+ai2*(3.d0+ai2*(1.d0-0.5d0*ai2)))
     1         )*exp(-ai2)
              ai2trm=ai2trm/(ai2*ai2*ai2*ai2)
        xb(1,1)=m*4.d0*num*nux/nu
        xb(2,1)=xb(1,1)*
     1     (1.d0-(1.d0+ais-ais*ais*0.5d0)*exp(-ais))/(ais*ais)
        xb(3,1)=xb(1,1)*1.5d0*num*zm*m
        xb(4,1)=xb(1,1)*2.d0*num*zm*m*ai2trm

        t2=t*t/500.d0/500.d0
        lot=t-225.d0
        lota=t-200.d0
        hit=650.d0-t
        hit3=hit*hit*hit 
       do 1 j = 1,4
         xb(j,2)=xb(j,1)*t/1.d+3
         xb(j,3)=xb(j,1)*t2
	 xb(j,4)=xb(j,1)/lota
         xb(j,5)=xb(j,1)*1.d+4/lota/lota/lota
         xb(j,6)=xb(j,1)*1.d2/lota/lota
         xb(j,7)=xb(j,1)*2.d2/t/t
         xb(j,8)=xb(j,1)*t2*t/500.d0
         xb(j,9)=xb(j,1)/sqrt(hit)
         xb(j,9)=xb(j,1)/sqrt(hit)
         xb(j,10)=xb(j,1)*1.d-5*p
         xb(j,11)=xb(j,1)*1.d-6*p*200.d0/lot
         xb(j,12)=xb(j,1)*p*1.d2/hit3
         xb(j,13)=xb(j,1)*1.d-5*p*t/500.d0
         xb(j,14)=xb(j,1)*1.d-6*p*200.d0/hit
         xb(j,15)=xb(j,1)*1.d-7*p*p
         xb(j,16)=xb(j,1)*1.d-8*p*p*200.d0/lot
         xb(j,17)=xb(j,1)*p*p/hit3
         xb(j,18)=xb(j,1)*1.d-7*p*p*t/500.d0
         xb(j,19)=xb(j,1)*p*p*t2*1.d-7
         xb(j,20)=xb(j,1)*1.d-6*4.d4*p/lot/lot
         xb(j,21)=xb(j,1)*1.d-5*p*t*t/500.d0/5.d2
         xb(j,22)=xb(j,1)*1.d-10*p*p*p*2.d2/lot
         xb(j,23)=xb(j,1)*p*p*p*1.d-2/hit3
         xb(j,24)=xb(j,1)*2.d2/hit/hit/hit
  1   continue

         beta0v=0.d0
         beta1v=0.d0
         cphiv=0.d0
         cphiv1=0.d0
         do 2 i=1, 24
           beta0v=beta0v+xb(1,i)*b0(i)
 2       continue
         do 3 i=1, 24
           beta1v=beta1v+xb(2,i)*b1(i)
 3       continue
         do 4 i=1, 24
           cphiv=cphiv+xb(3,i)*c0(i)
 4       continue
         do 5 i=1, 24
           cphiv1=cphiv1+xb(4,i)*c1(i)
 5       continue


      dhl=-abs(zm*zx)*aphit*((sqrt(is)/(1.d0+bb*sqrt(is))) +
     1 2.d0*log(1.d0+bb*sqrt(is))/bb)
      dhl=dhl+beta0v+beta1v+cphiv+cphiv1

      gam=exp(dhl)


!          if((bflag.eq.1.0).and.(m.ne.0.0d0))then
!            b0out=beta0v/xb(1,1)
!            b1out=beta1v/xb(2,1)
!            c0out=cphiv/xb(3,1)
!            c1out=cphiv1/xb(4,1)
!            print *,aphit,b0out,b1out
!             print *,c0out,c1out
!            pause
!          endif

        xb(1,1)=m*2.d0*num*nux/nu
        xb(2,1)=xb(1,1)*exp(-ais)
        xb(3,1)=xb(1,1)*num*zm*m*2.d0
        xb(4,1)=xb(3,1)*exp(-ai2)


       do 10 j = 1,4
         xb(j,2)=xb(j,1)*t/1.d+3
         xb(j,3)=xb(j,1)*t2
	 			 xb(j,4)=xb(j,1)/lota
         xb(j,5)=xb(j,1)*1.d+4/lota/lota/lota
         xb(j,6)=xb(j,1)*1.d2/lota/lota
         xb(j,7)=xb(j,1)*2.d2/t/t

         xb(j,8)=xb(j,1)*t2*t/500.d0
         xb(j,9)=xb(j,1)/sqrt(hit)
         xb(j,10)=xb(j,1)*1.d-5*p
         xb(j,11)=xb(j,1)*1.d-6*p*200.d0/lot
         xb(j,12)=xb(j,1)*p*1.d2/hit3
         xb(j,13)=xb(j,1)*1.d-5*p*t/500.d0
         xb(j,14)=xb(j,1)*1.d-6*p*200.d0/hit
         xb(j,15)=xb(j,1)*1.d-7*p*p
         xb(j,16)=xb(j,1)*1.d-8*p*p*200.d0/lot
         xb(j,17)=xb(j,1)*p*p/hit3
         xb(j,18)=xb(j,1)*1.d-7*p*p*t/500.d0
         xb(j,19)=xb(j,1)*p*p*t2*1.d-7
         xb(j,20)=xb(j,1)*1.d-6*4.d4*p/lot/lot
         xb(j,21)=xb(j,1)*1.d-5*p*t*t/500.d0/5.d2
         xb(j,22)=xb(j,1)*1.d-10*p*p*p*2.d2/lot
         xb(j,23)=xb(j,1)*p*p*p*1.d-2/hit3
         xb(j,24)=xb(j,1)*2.d2/hit/hit/hit
 10     continue

         beta0v=0.d0
         beta1v=0.d0
         cphiv=0.d0
         cphiv1=0.d0
         do 20 i=1, 24
           beta0v=beta0v+xb(1,i)*b0(i)
 20       continue
         do 30 i=1, 24
           beta1v=beta1v+xb(2,i)*b1(i)
 30       continue
         do 40 i=1, 24
           cphiv=cphiv+xb(3,i)*c0(i)
 40       continue
         do 50 i=1, 24
           cphiv1=cphiv1+xb(4,i)*c1(i)
 50       continue

      dhl=-abs(zm*zx)*aphit*sqrt(is)/(1.d0+bb*sqrt(is))
      dhl=dhl + beta0v + beta1v + cphiv + cphiv1
      phi =  1.d0 + dhl
      deallocate (xb)
      return
      end

c*****************************************************************************
          subroutine lphi(t,p,m,hext)
          implicit double precision (a-h,o-z)
          double precision m,num,nux,nu,mr,lot,lota,lot2,lot3
          common /refs/tr,pr,mr
          common /params/b0(24),b1(24),c0(24),c1(24),xv(10),xc(11)
          common /dielc/d,dddp,d2ddp2,dddt,d2ddt2,rho,drhodp,
     1                    drhodt,d2rdp2,d2rdt2

          double precision, dimension (:,:), allocatable :: xb
      		allocate (xb(4,24))
      			
          data r/.0083144d0/
          data num,nux,nu,zm,zx/1.d0,1.d0,2.d0,1.d0,-1.d0/
          data iac,iak/0,0/
          if(m.le.0.0d0)then
                 hext=0.0d0
                 return
          endif

          call dhlln(t,p,ap,ah,av,ac,ak,iac,iak)
           
          ah=ah*r*t
          hext=ah*nu*dlog(1.d0+1.2d0*sqrt(m))/(2.d0*1.2d0)

          ais=2.0d0*sqrt(m)
          ai2=2.5d0*sqrt(m)           
          ai2trm=(6.d0-(6.d0+ai2*(6.d0+ai2*(3.d0+ai2)))*exp(-ai2)
     1       )/((ai2*ai2)*(ai2*ai2))           

      xb(1,1)=-r*t*t*m*2.d0*num*nux
      xb(2,1)=xb(1,1)*2.d0*(1.d0-(1.d0+ais)*exp(-ais))/(ais*ais)
      xb(3,1)=xb(1,1)*num*zm*m
      xb(4,1)=xb(3,1)*4.d0*ai2trm

        tt=t/500.d0
        lot=t-225.d0
        lota=t-200.d0
!        lot=t-227.d0
        hit=650.d0-t
        hit2=hit*hit
        hit4=hit2*hit2
        lot2=lot*lot
        lot3=lot2*lot
      do 1 j = 1,4
         xb(j,2)=xb(j,1)/1.d+3
         xb(j,3)=xb(j,1)*tt*2.d0/500.d0
	       xb(j,4)=-xb(j,1)/lota/lota
!         xb(j,5)=-xb(j,1)*4.d+4/t/t
!         xb(j,6)=-2.d0*xb(j,1)*12.d+6/t/t/t
         xb(j,5)=-3.d0*xb(j,1)*1.d+4/lota/lota/lota/lota
         xb(j,6)=-2.d0*xb(j,1)*1.d+2/lota/lota/lota
         xb(j,7)=-2.d0*xb(j,1)*2.d2/t/t/t
         xb(j,8)=3.d0*xb(j,1)*tt*tt/500.d0
         xb(j,9)=0.5d0*xb(j,1)/(hit**1.5d0)

         xb(j,10)=0.0d0
         xb(j,11)=-xb(j,1)*p*1.d-6*2.d2/lot2
         xb(j,12)=3.d0*xb(j,1)*p*1.d2/hit4
         xb(j,13)=xb(j,1)*1.d-5*p/500.d0
         xb(j,14)=xb(j,1)*p*2.d2*1.d-6/hit2
         xb(j,15)=0.0d0
         xb(j,16)=-xb(j,1)*p*p*1.d-8*2.d2/lot2
         xb(j,17)=3.d0*xb(j,1)*p*p/hit4
         xb(j,18)=xb(j,1)*p*p*1.d-7/500.d0
         xb(j,19)=2.d0*xb(j,1)*p*p*tt*1.d-7/500.d0
!         xb(j,20)=xb(j,1)*0.0d0
         xb(j,20)=-2.d0*xb(j,1)*1.d-6*4.d4*p/lot/lot/lot
         xb(j,21)=2.d0*xb(j,1)*1.d-5*p*t/500.d0/5.d2
!         xb(j,21)=xb(j,1)*p*p*p*1.d-10/5.d2
         xb(j,22)=-xb(j,1)*p*p*p*1.d-10*2.d2/lot2
         xb(j,23)=xb(j,1)*p*p*p*3.d0*1.d-2/hit4
         xb(j,24)=3.d0*xb(j,1)*2.d2/hit4
         xb(j,1)=0.0d0
  1   continue


         beta0=0.d0
         beta1=0.d0
         cphi=0.d0
         cphi1=0.d0
         do 2 i=1, 24
           beta0=beta0+xb(1,i)*b0(i)
 2       continue
         do 3 i=1, 24
           beta1=beta1+xb(2,i)*b1(i)
 3       continue
         do 4 i=1, 24
           cphi=cphi+xb(3,i)*c0(i)
 4       continue
         do 5 i=1, 24
           cphi1=cphi1+xb(4,i)*c1(i)
 5       continue

            
         hext=hext+beta0+beta1+cphi+cphi1
         deallocate (xb)
         return
         end


c*****************************************************************************
            subroutine s20(t,p,s)
            implicit double precision (a-h,o-z)
            double precision mr,m
            common /refs/tr,pr,mr
            common /refs2/gextr,scorr,sions,sh2otr,gh2otr,sexmr
            common /params/b0(24),b1(24),c0(24),c1(24),xv(10),xc(11)
            common /dielc/d,dddp,d2ddp2,dddt,d2ddt2,rho,drhodp,
     1			  drhodt,d2rdp2,d2rdt2
		 				common /debye/aphi,ah,av
            common /dihyd/dgh2ot

! calculate difference in entropy of water (kg).
             call wsteam(t,rho,p,si,dpr,ddt,u,h,s,cv,cp,a,4)
                  sh2ot=s+scorr
                  sh2otl=(sh2otr-sh2ot)/mr
                  

! calculate differences in sex.

                call gamcal(t,p,mr,gamma,phic)        
               gextmr=2.d0*mr*t*.0083144d0*(log(gamma)+1.d0-phic)
               call lphi(t,p,mr,hexmr)
               sext=(hexmr-gextmr/mr)/t
               
               sextl=sexmr-sext
             
!  calculate the g20(t,p) - g20(tr,pr)

!            evaluate integrals
               ctot=0.0d0

!  first int (1/t**2) int cp,pr,mr, dt dt


                tinv=1.d0/t-1.d0/tr

               ctot=ctot+xc(1)*dlog(t/tr)
               ctot=ctot+xc(3)*(t-tr)/3.d2
               ctot=ctot+0.5d0*xc(4)*(t*t-tr*tr)/3.d2/3.d2
               ctot=ctot-xc(7)*tinv*1.d2

               delp=p-pr
               delp2=0.5d0*(p*p-pr*pr)
               delp3=(p*p*p-pr*pr*pr)/3.d0
               delpp=((p+10.d0)**2.5d0-(pr+10.d0)**2.5d0)/2.5d0
               dt=t/300.d0
               vtot=0.0d0
            	 vtot=vtot+xv(2)*delp/3.d2/1.d2
            	 vtot=vtot+3.d0*xv(4)*delp*dt*dt/1.d3/3.d2
            	 vtot=vtot+xv(6)*delp2/1.d2/1.d1/3.d2
            	 vtot=vtot+2.d0*xv(7)*dt*delp2/1.d3/1.d2/3.d2
            	 vtot=vtot+xv(9)*delp3/1.d3/1.d4/3.d2
            	 vtot=vtot+xv(3)*delpp/1.d7/3.d2

               s=sions+sh2otl+sextl-vtot+ctot
               return
               end



c*****************************************************************************
            subroutine delg20(t,p,delg)
            implicit double precision (a-h,o-z)
            double precision mr,m
            common /refs/tr,pr,mr
            common /refs2/gextr,scorr,sions,sh2otr,gh2otr,sexmr
            common /params/b0(24),b1(24),c0(24),c1(24),xv(10),xc(11)
            common /dielc/d,dddp,d2ddp2,dddt,d2ddt2,rho,drhodp,
     1			  drhodt,d2rdp2,d2rdt2
		 				common /debye/aphi,ah,av
            common /dihyd/dgh2ot

            double precision, dimension (:), allocatable :: xca
      			allocate (xca(7))
      			
            call wsteam(t,rho,p,si,dpr,ddt,u,h,s,cv,cp,a,4)
            gh2ot=h-t*(s+scorr)
            dgh2ot=(gh2otr-gh2ot)/mr

             
!  calculate the g20(t,p) - g20(tr,pr)

!           evaluate integrals
            tinv=1.d0/t-1.d0/tr
		
            ts=tr
            xca(1)=t-ts
            xca(2)=-3.d2*(1.d0/t-1.d0/ts)
            xca(3)=0.5d0*(t*t-ts*ts)/300.d0
            xca(4)=(t*t*t-ts*ts*ts)/27.d4
            xca(5)=t/(680.d0-t)-ts/(680.d0-ts)-log((680.d0-ts)
     1                /(680.d0-ts))
            xca(6)=(t*t*t*t-ts*ts*ts*ts)/27.d6/4.d0
            xca(7)=log(t/ts)

            xca(1)=xca(1)-t*log(t/ts)
            xca(2)=xca(2)+t*0.5d0*3.d2*(1.d0/(t*t)-1.d0/(ts*ts))
            xca(3)=xca(3)-t*(t-ts)/3.d2
            xca(4)=xca(4)-t*0.5d0*(t*t-ts*ts)/9.d4
            xca(6)=xca(6)-t*(t*t*t-ts*ts*ts)/27.d6/3.d0
            xca(7)=xca(7)+t*(1.d0/t-1.d0/ts)
            xca(7)=xca(7)*1.d2
            ctot=0.0d0
            do 1010 i=1,7
            		ctot=ctot+xc(i)*xca(i)
 1010       continue


            delp=p-pr
            delp2=0.5d0*(p*p-pr*pr)
            delp3=(p*p*p-pr*pr*pr)/3.d0
            delpp=((p+10.d0)**2.5d0-(pr+10.d0)**2.5d0)/2.5d0
            dt=t/300.d0
            vtot=0.0d0
            vtot=vtot+xv(1)*delp*1.d-1
            vtot=vtot+xv(2)*delp*dt/1.d2
            vtot=vtot+xv(3)*delp*dt*dt/1.d3
            vtot=vtot+xv(4)*delp*dt*dt*dt/1.d3
            vtot=vtot+xv(5)*delp2/1.d3/1.d2
            vtot=vtot+xv(6)*dt*delp2/1.d2/1.d1
            vtot=vtot+xv(7)*dt*dt*delp2/1.d3/1.d2
            vtot=vtot+xv(8)*dt*dt*dt*delp2/1.d2/1.d4
            vtot=vtot+xv(9)*dt*delp3/1.d3/1.d4
            vtot=vtot+xv(10)*delpp/1.d6
            vtot=vtot+xv(3)*delpp*dt/1.d7
						
            call gamcal(t,p,mr,gam,phic)
            gext=2.d0*mr*t*.0083144d0*(log(gam)+1.d0-phic)
            dgexmr=(gextr-gext)/mr



            delg=dgh2ot+dgexmr+vtot+ctot

            delg=delg-(t-tr)*(sions+sh2otr/mr+sexmr)

            deallocate (xca)
            return
            end
