!****   subroutine to convert old dhll calls into new calls.

         subroutine dhll(tin,pin,in,aphi,ah,ac,av,iac)
         implicit double precision (a-h,o-z)
         t=tin
         p=pin
!         iac=1
          iak=0
!         iak=1
         call dhlln(t,p,aphi,ah,av,ac,ak,iac,iak)
!         print *,' returning from dhll  aj = ',ac
         return
         end

!****   subroutine to convert hgk calls into hill calls.
!***    not all variables ar passed.  be careful !!!!!!!
 
        subroutine hgk(t,p,i,rho,cp,dpd,dpdt1,dvwdt,k,d2t,dtr,d2r)
        implicit double precision (a-h,o-z)
!        print *, 
        lsav=4
        call wsteam(t,rho,p,si,dpdr,dpdt,u,h,s,cv,cp,a,lsav) 
          drhodt=-dpdt/dpdr
          dvwdt=-drhodt/rho/rho
          dpd=dpdr
!        write(5,100)t,p,rho
 100    format(1x,' t = ',f6.2,' p = ',e12.6,' rho = ',f10.8)
        return
        end
 
!**********************************************************************
! 	subroutine to calculate debye-huckel slopes.
!     
!                 draft version 1.01, 10/6/89. 
!                 the user of these subroutines assumes all liability 
!                 for the accuracy of the program and the results 
!                 of any calculations performed with this code.
!
!       usage:
!         call dhll(t,p,aphi,ah,av,ac,ak,iac,iak)
!
!       arguments:
!          t       temperature in kelvin
!          p       pressure in mpa
!          aphi    limiting law slope for osmotic coefficient (unitless)
!          ah      limiting law slope for apparent molar enthalpy divided
!                     by rt (unitless)
!          av      limiting law slope for apparent molar volume. (units
!                     are cubic centimeters per mole)
!          ac      limiting law slope for apparent molar heat capacity
!                     divided by r. (unitless)
!          ak      limiting law slope for apparent molar compressibility.
!                     (units are cubic centimeters per mpa per mole)
!          iac     switch for computation of ac (1 = yes); integer
!          iak     switch for computation of ak (1 = yes); integer
! 
!       **    note that, in above, there exists an inverse square root of 
!       **    ionic strength unit, which was not explicitly stated.
!      
!          all arguments other than iac and iak are double precision.
!
!          the two integer switches control calculation of second derivative
!          functions.  these calculations consume more computation time.
!
!          the subroutine contains a common which is named dielc.
!          this common contains the dielectric constant and the 
!          first and second derivatives of the dielectric constant with
!          respect to t and p: d, dddt, d2ddt2, dddp, and d2ddp2 (but
!          not in this order).  the common also contains the density, rho,
!          and its derivatives: drhodt, d2rdt2, drhodp, and d2rdp2 (but
!          not in this order).
!
!          donald g. archer 
!          electrolyte data center
!          national institute of standards and technology
!          gaithersburg, md 20899
!          
!
!**********************************************************************

      subroutine dhlln(t,p,aphi,ah,av,ac,ak,iac,iak)
      implicit double precision (a-h,o-z)
!****
!****    common dielc contains the values of the derivatives of the 
!****    dielectric constant and water density after the call
!****    to subroutine diel.
!****
      common /dielc/ d,dddp,d2ddp2,dddt,d2ddt2,rho,drhodp,drhodt,
     1 d2rdp2,d2rdt2

        double precision k,mu,mu2

        data wm,pi,avno,alpha,k,mu/18.0153,3.14159,6.0221367e+23,
     1  1.444e-24,1.380658e-23,1.84e-18/

        data e,eps0,r/1.60217733e-19,8.8541878e-12,8.314510/

!   get the dielectric constant, density, and their derivatives.

           call diel(t,p,iac,iak)

           aphi=sqrt(2.*pi*rho*avno/1000.)
     1               *((e*e*100./(4.*pi*d*eps0*k*t))**1.5)/3.

!   calculate a few more derivatives.

           alphav=-drhodt/rho
           dalpat=alphav*alphav-d2rdt2/rho
           rd2vdt=2.*drhodt*drhodt/(rho*rho)-d2rdt2/rho
           betav=drhodp/rho
           dbetap=-betav*betav+d2rdp2/rho
           dlnddt=dddt/d
           dlnddp=dddp/d
           d2lndt=-dddt*dddt/(d*d)+d2ddt2/d
           d2lndp=-dddp*dddp/(d*d)+d2ddp2/d

!   calculate some required terms.

           braket=1./t+dlnddt+alphav/3.
           brak2=-1./(t*t)+d2lndt+dalpat/3.
           vbrak=dlnddp-betav/3.

!   calculate remainder of debye-huckel parameters.

           ah=-6.*aphi*t*braket
           av=6.*aphi*r*t*vbrak

           if(iac .eq. 1)then
            ac = 3.*aphi*t*(1./t+2.*dlnddt+5.*t
     1              *dlnddt*dlnddt+2.*t*dlnddt*alphav+2.*alphav/
     2              3.+t*alphav*alphav-2.*t*d2ddt2/d-2.*t
     3              *rd2vdt/3.)
            else
            ac=0.0
           endif
           if(iak .eq. 1)then 
	    ak=6.*r*t*aphi*(-1.5*vbrak*vbrak+d2lndp
     1             -dbetap/3.)
           else
            ak=0.0
           endif

!*****    that's all folks.

        return
        end

      subroutine diel(t,p,iac,iak)
      implicit double precision (a-h,o-z)
      common /dielc/ d,dddp,d2ddp2,dddt,d2ddt2,rho,drhodp,drhodt,
     1 d2rdp2,d2rdt2
       common/cpcalc/cpsave
        dimension x(12)
        double precision k,mu,mu2

!****  assign constants.

        data wm,pi,avno,alpha,k,mu/18.0153,3.1415927,6.0221367e+23,
     1 1.444e-24,1.380658e-16,1.84e-18/

!****  assign parameters.

        data x/0.0,-3.548184,-.4044525e-1,0.0000000,
     1 -1246.311,-.6928953,.2633077e+6,103.6180,75.32165,
     2 -23.23778,-204.4473,0.0/

            tk3=3.*k*t
            mu2=mu*mu
            tr=215.
            lsav=4

!*****   get rho and first deriviatives for t and p.
!*****

             call wsteam(t,rho,p,si,dpdrho,dpdt,u,h,s,cv,cpsave,
     1            a,lsav)  
             drhodp=1./dpdrho
             drhodt=-dpdt*drhodp

!*****   if flagged, get second derivatives of rho with respect to t and p.
!*****
            delp=0.01
            delt=0.001

             if(iak .eq. 1)then
             
                 call wsteam(t,rho1,p-delp,si,dpdr1,dpdt,u,h,s,cv,cp,
     1            a,lsav)  
                 call wsteam(t,rho2,p+delp,si,dpdr2,dpdt,u,h,s,cv,cp,
     1            a,lsav)  
           
                  drdp1=1./dpdr1
                  drdp2=1./dpdr2
                  d2rdp2=(drdp2-drdp1)/(2.*delp)
             else
                  d2rdp2=0.0
             endif

             if(iac .eq. 1)then

                 call wsteam(t-delt,rho1,p,si,dpdr1,dpdt1,u,h,s,cv,cp,
     1            a,lsav)  
                 call wsteam(t+delt,rho2,p,si,dpdr2,dpdt2,u,h,s,cv,cp,
     1            a,lsav)
  
                 drdp1=1./dpdr1
                 drdp2=1./dpdr2
                 drdt1 = -dpdt1*drdp1
                 drdt2 = -dpdt2*drdp2
                 d2rdt2=(drdt2-drdt1)/(2.*delt)
             else
                 d2rdt2=0.0
             endif

!*****   calculate g.
!*****

	 exptrm = exp(x(5)/t+x(6)*p/t+x(7)/(t*t)+
     1		  x(11)*p/(t*t)+x(12)*p)
	 y  =  x(1)+x(2)/((t-tr)**.25)+x(3)*p/t+x(8)/sqrt(t)+
     1         x(9)/((t-tr)**1.00)+x(10)/((t-tr)**0.5)+exptrm
	 g  =  y*rho+1.

!*****  calculate first and second temperature derivatives of y.
!*****

        dpwrdt  = (-x(5)/t-x(6)*p/t-2.*(x(7)+x(11)*p)/(t*t))/t

        dydt    = -x(8)/(2.*(t**1.5))-x(9)/((t-tr)**2.)-x(3)*p
     1              /(t*t)-x(10)/(2.*((t-tr)**1.5))+exptrm*dpwrdt
     2              -x(2)/(4.*((t-tr)**1.25))

        if(iac .eq. 1)then 

          d2prdt  = (2.*x(5)/t+2.*x(6)*p/t+6.*(x(7)+x(11)*p)/(t*t))
     1              /(t*t)
        d2ydt2  = 3.*x(8)/(4.*(t**2.5))+2.*x(9)/((t-tr)**3.)
     1		    +1.5*x(10)/(2.*((t-tr)**2.5))+exptrm*dpwrdt*dpwrdt
     2              +exptrm*d2prdt
     3              +2.*x(3)*p/(t*t*t)+1.25*x(2)/(4.*((t-tr)**2.25))
        else
           d2prdt=0.0
           d2ydt2=0.0
        endif

!*****  calculate first and second pressure derivatives of y.
!*****
        dpwrdp = x(6)/t+x(12)+x(11)/(t*t)
        dydp   = x(3)/t+exptrm*dpwrdp

        if (iak .eq. 1)then
          d2ydp2 = exptrm*dpwrdp*dpwrdp
        else
          d2ydp2=0.0
        endif
        
!*****  calculate first and second temperature derivatives of g.
!*****
        dgdt   = dydt*rho+y*drhodt
        d2gdt2 = d2ydt2*rho+2.*dydt*drhodt+y*d2rdt2

!*****  calculate first and second pressure derivatives of g.
!*****
        dgdp   = dydp*rho+y*drhodp
        d2gdp2 = d2ydp2*rho+2.*dydp*drhodp+y*d2rdp2

!*****   calculate a.
!*****
           a = g*mu*mu/tk3
           a = (a+alpha)*4.*pi*avno*rho/(3.*wm)

!*****   calculate first and second temperature derivatives of a.
!*****

         dadt   = (drhodt*(alpha/mu2+g/tk3)+rho*(dgdt-g/t)/tk3)
     1            *4.*pi*avno*mu2/(3.*wm)

       if (iac .eq. 1)then  
         d2adt2 = (d2rdt2*(alpha/mu2+g/tk3)+2.*drhodt*(dgdt-g/t)
     1            /tk3+rho*(2.*g/t-2.*dgdt+t*d2gdt2)
     2            /(tk3*t))*4.*pi*avno*mu2/(3.*wm)             
       endif

!*****   calculate first and second pressure derivatives of a.
!*****
	 dadp	= (drhodp*(alpha/mu2+g/tk3)+rho*(dgdp/tk3))*4.*
     1            pi*avno*mu2/(3.*wm)
      if(iak .eq. 1)then
	 d2adp2 = (d2rdp2*(alpha/mu2+g/tk3)+drhodp*2.*dgdp/tk3+
     1		  rho*d2gdp2/tk3)*4.*pi*avno*mu2/(3.*wm)
      endif

!*****    calculate d from a.
!***** 
          apoly = 9.*a*a+2.*a+1.
          d     = 0.25*(1.+9.*a+3.*sqrt(apoly))

!*****    calculate first and second derivatives of d w.r.t. a.
!*****

          ddda   = (9. + 1.5*(18.*a+2.)/sqrt(apoly))/4.
       if(iac .eq. 1)then
          d2dda2 = (-0.75*(18.*a+2.)*(18.*a+2.)/(apoly**1.5)+27./
     1             sqrt(apoly))/4.
       endif

!*****   calculate first and second temperature derivatives of d.
!*****
          dddt   = ddda*dadt

          if(iac .eq. 1) d2ddt2 = ddda*d2adt2 + dadt*dadt*d2dda2

!*****   calculate first and second pressure derivatives of d.
!*****
          dddp   = ddda*dadp

          if(iak.eq.1) d2ddp2 = ddda*d2adp2 + dadp*dadp*d2dda2

          if(iak .ne. 1)d2ddp2=0.0
          if(iac .ne. 1)d2ddt2=0.0

!****  finished.  

	  return
           end
      

!********************************************************************
!      remainder of this file is the equation of state of hill.
!      this code compiles and executes properly under ms fortran and 
!      lahey f77l compilers.  errors in some calculations have been
!      found when executing under lahey's em32 compiler (v. 1.01).
!                                                    d.g.a. 4/26/89 
!********************************************************************



      subroutine wsteam(temk,rho,press,si,dpdr,dpdt,u,h,s,cv,cp,
     1      a,lsav)
!
!
!********************************************************************
!
!         variables in call statement
!
!********************************************************************
!         temk         temperature                     k
!         rho          density                         kj/dm3
!         press        pressure                        mpa
!         si           helmholtz free energy           kj/kg.k
!         dpdr         pressure-density derivative     mpa.dm3/kg
!         dpdt         pressure-temperature deriv.     mpa/k
!         u            internal energy                 kj/kg
!         h            enthalpy                        kj/kg
!         s            entropy                         kj/kg.k
!         cv           isochoric specific heat         kj/kg.k
!         cp           isobaric specific heat          kj/kg.k
!         a            isentropic sonic speed          m/s
!
!           additional values passed in common /carga/
!    
!           call routine swvir for virial coeff values
!
!********************************************************************
!
!         options
!
!********************************************************************
!
!         l = 1        no operation
!                           (code value used by automatic initialization
!                            procedure on first call to wsteam)
!
!         l = 2        given rho and temk calculate all variables
!         l = 3        given press and temk ( vapour side )
!                            calculate all variables
!         l = 4        given press and temk ( liquid side )
!                            calculate all variables
!         l = 5        given temk calculate saturation values
!
!                 [iterated sat values passed in labelled common]
!                 /satliq/  and  /satvap/
!
!
!         l = 6        given temk and press, calculate psat, then
!                             calculate all derivatives on liq or vap
!                             side, depending on input press value
!
!
!         notes:  error messages printed on fortran unit 6
!                approx sat values could be obtained using fpswag
!                correlation for ps with l=3(vap) or l=4 (liq)
!
!                for l=3 with temk where press .gt. psat wsteam may
!                calculate metastable states, (and
!                 similarly for l=4)
!           release 2.0   (june 1989)
!                comments and suggestions for improvement of these
!                routines should be directed to
!                         p.g. hill,
!                         dept. mech. engineering,
!                         university of british columbia
!                         vancouver, b.c., canada
!
!                this program is written in fortran/77
!
!                the approximate range of the input data used
!                to determine the coefficients is:
!
!                         273.15 <  temk  <  1273.15  k
!                          0.0   < press  <  1000.    mpa
!                          0.0   <  rho   <  1.       kj/dm3
!
!              (also included is a limited amount of data
!               at 293 k, up to approx 40,000 mpa)
!
!********************************************************************
!
!         normalizing constants
!
!********************************************************************
!
!         gas constant            rval = 0.46151 kj/kg.k
!         critical temperature      tc = 647.067 k
!         critical density        rhoc = 0.322778 kj/dm3
!
!*******************************************************************
!
!         dimensionless variables
!
!********************************************************************
!
!         r = rho/rhoc                dr = r - 1
!         t = -tc/temk                dt = t + 1
!         p = press/(rhoc*rval*temk)
!        si = si/(rval*temk)
!         u = u/(rval*tc)
!         h = h/(rval*tc)
!         s = s/rval
!        cv = cv/rval
!        cp = cp/rval
!         a = a/sqrt(rval*tc)
!
!********************************************************************
!
!        the dimensionless equation of state
!
!********************************************************************
!
!        si = sif + f*( sin - sif )
!        sif = si0 + dlog(r) + w1 + e*w2 + g*w3 + h*w4
!                in which
!        si0 = ideal gas function defined by subroutine si0f
!        sin =revised and extended scaling function of levelt-sengers
!             kamgar-parsi, balfour, and sengers (sept 1981) with
!             adjustment of the two arbitrary constants
!             evaluated in subroutine sinfun
!        w1 = sum over i and j of c1(i,j)*r**(i-2)*t**(j-1)*(1-e)
!         except when i=2 c1(i,j)*(ln(r)(1-e-r*r)+r*r/2)*t**(j-1)
!        w2 = sum over i and j of c2(i,j)*r**(i)*t**(j-1)*e
!        w3 = sum over i and j of c3(i,j)*r**(i+1)*t**(j+1)*g
!        w4 = sum over i and j 0f c4(i,j)*r**(i-2)*1-e)*t**(j-1)*h
!         except when i=2 c4(i,j)*(ln(r)(1-e-r*r)+r*r/2)*t**(j-1)
!        w1, w2, w3,,w4 and their derivatives evaluated
!        by subroutine wfuncs
!        e is defined by subroutine efunc
!        f is defined by subroutine ffunc
!        g is defined by subroutine gfunc
!        h is defined by subroutine hfunc
!
!*******************************************************************
!
!        derivatives of the dimensionless helmholtz function
!
!********************************************************************
!
!         notation     sir = first  partial derivative of si wrt r
!                     sirr = second partial derivative of si wrt r
!                     sirt = second partial derivative of si wrt r and t
!
!         pressure           p = r*r*sir
!                           pr = 2*r*sir + r*r*sirr
!                           pt = r*r*sirt
!         internal energy    u = -sit
!         enthalpy           h = -sit - r*sir/t
!         entropy            s = -si + t*sit
!         isochoric sp ht   cv = -t*t*sitt
!         isobaric sp ht    cp = cv + r*r*pr*( t*pt - p )**2
!         isentropic speed   a = sqrt( -cp/( cv*t ) * pr )
!           of sound
!
!********************************************************************
!
       implicit double precision(a-h,o-z)
!      deleted unused arrays d(3),dd(3)
       double precision w(4,6),c0(8)
       double precision qv(20)
       common/cdel/cdelta
       common/ccread/cread
       common/cerr/cerrtk,cerrho,cerrp,cerrl,cerrls
       common/cqv/qv
       common/efgval/eval,dr0,dt0,delta,aval,bval,cval,dval,hval
       common/norms/tc,rhoc,rval
       common/iter/iterqk
       common/osg/thet48,posg,rfosg,rgosg,hfosg,hgosg
       common/sat/psat,rf,rg,hf,hg,sf,sg
       common/czero/c0
       common/cinit/c000,c001
       common/satliq/ctemkf,crhof,cpf,csif,cdpdrf,cdpdtf,cuf,chf,
     *               csf,ccvf,ccpf,caf
       common/satvap/ctemkg,crhog,cpg,csig,cdpdrg,cdpdtg,cug,chg,
     *               csg,ccvg,ccpg,cag
!
       common/carga/dhdpt,dtdph,dtdps,expis,x5,x6,x7,x8,x9,x10
       common/cargf/dhdptf,dtdphf,dtdpsf,expisf,xf5,xf6,xf7,xf8,xf9,xf10
       common/cargg/dhdptg,dtdphg,dtdpsg,expisg,xg5,xg6,xg7,xg8,xg9,xg10
       common/cscnt/icscnt
       common/cwcnt/icwcnt
!
       data inival/0/
!
!      cread set to 1.0d0 if coefficients read from fortran unit 10
!
!       Coef132 is not used
!      open(unit=10,file='coef132.dat')
       psp=0.0d0
          l=lsav
       if(l.eq.1)return
       pc=22.046d0
       if(l.ne.2) rho= 0.0d0
       if((l.eq.2).or.(l.eq.5)) press= 0.0d0
       cerrp= 0.0d0
       cerrtk=temk
       cerrho=rho
       if((l.eq.3).or.(l.eq.4).or.(l.eq.6)) cerrp=press
       cerrl= l
!
!
       if(inival.eq.0) iterqk=1
       if(inival.eq.0) then
       call winit(q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11,q12,1)
            inival=1
            endif
!
!       check value of call arguments
!
       if (temk.le.0.0d0) then
                          write(6,9901) temk
                          return
                          endif
9901   format('0','temk .le. 0.0 ',e25.15)
!
       if((l.eq.2).and.(rho.lt.0.0d0)) then
                                        write(6,9902)rho
                                       return
                                       endif
9902   format('0','rho .lt. 0.0 ',e25.15)
!
!
!        compare input press argument with wsteam spinodal
!                pressure to guard against unsuccessful or
!                incorrect iteration in wrho routine to
!                supersaturated metastable states
!
!        (note that for l=4, where the supersaturated pressure can be
!                           negative, wrho can be called
!                           with  press < 0.0)
!
       if((l.ne.3).and.(l.ne.4).and.(l.ne.6)) go to 500
       if((l.eq.6).and.(press.gt.0.0d0))go to 500
       if((l.eq.6).and.(press.le.0.0d0)) go to 599
       if((l.eq.3).and.(press.le.0.0d0)) go to 599
       if((l.eq.3).and.(temk.ge.tc)) go to 500
       if((l.eq.4).and.(press.le.0.d0).and.(temk.ge.tc))go to 599
       if((l.eq.4).and.(temk.ge.tc)) go to 500
       psw= fpswag(temk)
       if((l.eq.3).and.(temk.lt.tc).and.(press.le.psw))go to 500
       if(l.eq.3)then
       sping= wsping(temk)
       call wderiv(temk,sping,psp,q4,q5,q6,q7,q8,q9,q10,q11,q12,2)
            endif
       if((l.eq.3).and.(temk.lt.tc).and.(press.gt.psp*1.0001))then
                     write(6,9599) temk,press,psp
9599      format(' ','wsteam: input press gt p(spinodal)',
     &              ' for l=3',3f15.8)
                     write(6,9597)
9597      format(' is l value correct?')
                     rho=-10.
                     return
                     endif
       if(l.eq.3) go to 500
!
       if((l.eq.4).and.(temk.lt.tc).and.(press.ge.psw))go to 500
       if(l.eq.4)then
       spinf= wspinf(temk)
       if((temk.ge.311.15).and.(temk.le.423.15))spinf=spinf*1.001
       call wderiv(temk,spinf,psp,q4,q5,q6,q7,q8,q9,q10,q11,q12,2)
                   endif
!      
       if((l.eq.4).and.(temk.lt.tc).and.(press.lt.psp*.999))then
                     write(6,9598) temk,press,psp
9598      format(' ','wsteam: input press lt p(spinodal)',
     &              ' for l=4',3f15.8)
                     write(6,9597)
                     rho=-10.
                     return
                     endif
       if(l.eq.4) go to 500
!
599    continue
       write(6,9903)temk,press,l
9903   format(' ','wsteam: press value out of range; ',2e15.6,i3)
       rho= -10.0
       return
!
500    continue
!
!
          iterqk = 0
!      set iterqk = 1 to reduce iteration time for wrho(p,t) 
          if( l .ge. 3 .and. l .le. 5 ) iterqk = 1
          tsave = 0.d0
!
       if(l.eq.2) then 
      call wderiv(temk,rho,press,si,dpdr,dpdt,u,h,s,cv,cp,a,l)
           return
           endif
!
       if((l.eq.3).or.(l.eq.4)) then
      call wrho(temk,rho,press,si,dpdr,dpdt,u,h,s,cv,cp,a,l)
            return
            endif
!
       if(l.eq.5) then
      call wsat(temk,rho,press,si,dpdr,dpdt,u,h,s,cv,cp,a,l)
                 return
                 endif
!
       if((l.eq.6).and.(temk.gt.tc)) then
      call wrho(temk,rho,press,si,dpdr,dpdt,u,h,s,cv,cp,a,4)
                  return
                  endif
!
!
       if(l.eq.6) 
     *  call wsat(temk,rho,psc,si,dpdr,dpdt,u,h,s,cv,cp,a,5)
       if((l.eq.6).and.(press.ge.psc))
     *  call wrho(temk,rho,press,si,dpdr,dpdt,u,h,s,cv,cp,a,4)
!
       if((l.eq.6).and.(press.lt.psc))
     *  call wrho(temk,rho,press,si,dpdr,dpdt,u,h,s,cv,cp,a,3)
       if(l.eq.6) return
!
       write(6,999) l
999    format(' ','error - l value action undefined: ',i3)
       return
       end
       
       
       subroutine winit(q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11,q12,l)
          implicit double precision(a-h,o-z)
!         deleted unused arrays d(3),dd(3)
          double precision w(4,6),c0(8)
          double precision qv(20)
       common/ccread/cread
        common/cerr/cerrtk,cerrho,cerrp,cerrl,cerrls
          common/cqv/qv
          common/efgval/eval,dr0,dt0,delta,aval,bval,cval,dval,hval
          common/norms/tc,rhoc,rval
          common/iter/iterqk
          common/osg/thet48,posg,rfosg,rgosg,hfosg,hgosg
          common/sat/psat,rf,rg,hf,hg,sf,sg
        common/czero/c0
       common/cinit/c000,c001
       common/satliq/ctemkf,crhof,cpf,csif,cdpdrf,cdpdtf,cuf,chf,
     *               csf,ccvf,ccpf,caf
       common/satvap/ctemkg,crhog,cpg,csig,cdpdrg,cdpdtg,cug,chg,
     *               csg,ccvg,ccpg,cag
!
       common/carga/dhdpt,dtdph,dtdps,expis,x5,x6,x7,x8,x9,x10
       common/cargf/dhdptf,dtdphf,dtdpsf,expisf,xf5,xf6,xf7,xf8,xf9,xf10
       common/cargg/dhdptg,dtdphg,dtdpsg,expisg,xg5,xg6,xg7,xg8,xg9,xg10
!
!********************************************************************
!
!        initialize
!
!********************************************************************
!
!      constants for triple liq ref state
!
!       c0(1)= 0.707510440802d1
!       c0(2)= -.834258572635d1
!       c000= 0.168157893281d-3
!       c001=  -.390307421339d-5
!
!
 1        rval = 0.46151d0
          tc = 647.067d0
          rhoc = 0.322778d0
!      constant for e function
          eval = 1.d0
!      constants for f function
          dr0 = 0.23d0
          dt0 = 0.05d0
          delta = 1.028667d0
!      constants for g function
          aval = 80.d0
          bval = 1.d0
          cval = 130.d0
          dval = 12.d0
          hval = 4.0d0
!
!      constants for ideal gas function
!       Coef132 is not used
!      if(abs(cread-1.d0)<1.0E-26)
!    *read(10,1010) c0(1),c0(2), (c0(i), i = 3,8 )
!1010  format(1x,2d18.12,6d14.8)
!      write(5,1010)c0
!
!
!      coefficients for w1, w2, w3, and w4 functions
          call wfuncs( q,q1,w,1)
!       Coef132 is not used
!      additions to arbitrary constants of scaling function
!      if(abs(cread-1.d0)<1.0E-26)  read(10,99) c000, c001
!99       format(1x,10d16.10)
!
!
!
!************************  c000, c001, c0(1), c0(2) were all entered
!                           into block data statement as they appeared
!                           in hill's steam table.  thus, next line of 
!                           code is commented so that the correction is
!                           made to the coefficients.
!
!       if(cread.ne.1.0d0) return
!************************************************************************
!
!      this section of code calculates
!      adjustments to triple point liquid reference state
!       (required when values of coefficients are changed)
        temk = 273.16d0
       call wsat(temk,rho,press,si,dpdr,dpdt,u,h,s,cv,cp,a,l)
 100        uref = hf - press/rf
        sref = sf
        c0(1) = c0(1) - uref/rval/tc
        c0(2) = c0(2) + sref/rval
        c000 = c000 + sref/rval - uref/rval/tc
        c001 = c001 + uref/rval/tc
!        write(6,101)
!        write(6,102) c0(1),c0(2)
!        write(6,103) c000,c001
 101        format(10x,'arbitrary constants after adjustment'
     1 'to triple point liquid state')
 102        format(15x,'c0(1) = ',d18.12,5x,'c0(2) = ',d18.12)
 103      format(15x,'c000  = ',d18.12,5x,'c001  = ',d18.12)
          return
       end
       
       subroutine wderiv(temk,rho,pcalc,si,dpdr,dpdt,u,h,s,cv,cp,
     *                   a,l)
          implicit double precision(a-h,o-z)
!         deleted unused arrays d(3),dd(3)
          double precision w(4,6),c0(8)
          double precision qv(20)
       common/ccread/cread
        common/cerr/cerrtk,cerrho,cerrp,cerrl,cerrls
          common/cqv/qv
       common/fxxx/fxval
          common/efgval/eval,dr0,dt0,delta,aval,bval,cval,dval,hval
          common/norms/tc,rhoc,rval
          common/iter/iterqk
          common/osg/thet48,posg,rfosg,rgosg,hfosg,hgosg
          common/sat/psat,rf,rg,hf,hg,sf,sg
        common/czero/c0
       common/cinit/c000,c001
       common/satliq/ctemkf,crhof,cpf,csif,cdpdrf,cdpdtf,cuf,chf,
     *               csf,ccvf,ccpf,caf
       common/satvap/ctemkg,crhog,cpg,csig,cdpdrg,cdpdtg,cug,chg,
     *               csg,ccvg,ccpg,cag
!
       common/carga/dhdpt,dtdph,dtdps,expis,x5,x6,x7,x8,x9,x10
       common/cargf/dhdptf,dtdphf,dtdpsf,expisf,xf5,xf6,xf7,xf8,xf9,xf10
       common/cargg/dhdptg,dtdphg,dtdpsg,expisg,xg5,xg6,xg7,xg8,xg9,xg10
!
!********************************************************************
!
!        given rho and temperature derive all properties (l=2)
!
!********************************************************************
!
 20     r = rho/rhoc
       cerrho=rho
        t = -tc/temk
        dt = t + 1.d0
!
       if(abs(rho)<1.0E-26)then
            pcalc=0.0d0
            dpdt=0.0d0
            dpdr=1.0d0*rval*temk
            call si0f(t,si0,si0t,si0tt,2)
            cv=-t*t*si0tt
            cp= cv+1.d0
            u= -si0t
            h= (u-1.d0/t)
! 
            cv= cv*rval
            cp=cp*rval
            u= u *rval*tc
            h= h * rval*tc
            a= dsqrt(cp/cv*rval*temk)
            return
        endif

        call wfuncs( t,r,w,2 )
        w1   = w(1,1)
        w1r  = w(1,2)
        w1rr = w(1,3)
        w1t  = w(1,4)
        w1tt = w(1,5)
        w1rt = w(1,6)

        w2   = w(2,1)
        w2r  = w(2,2)
        w2rr = w(2,3)
        w2t  = w(2,4)
        w2tt = w(2,5)
        w2rt = w(2,6)

        w3   = w(3,1)
        w3r  = w(3,2)
        w3rr = w(3,3)
        w3t  = w(3,4)
        w3tt = w(3,5)
        w3rt = w(3,6)

        w4   = w(4,1)
        w4r  = w(4,2)
        w4rr = w(4,3)
        w4t  = w(4,4)
        w4tt = w(4,5)
        w4rt = w(4,6)


        call efunc( r,e,er,er2 )
        call ffunc( r,t,f,fr,frr,ft,ftt,frt )
        call gfunc( r,t,g,gr,grr,gt,gtt,grt )
        call hfunc( r,t,h,hr,hrr,ht,htt,hrt )
        call si0f( t,si0,si0t,si0tt,2 )

!
!********************************************************************
!
!       far field functions
!
!********************************************************************
!
         rlog = -170.d0
         if( r .gt. 1.5d-37 ) rlog = dlog(r)
        sif = si0 + rlog + w1 + e*w2 + g*w3 + h*w4
        pf = r + r*r*( w1r + er*w2 + e*w2r + gr*w3 +g*w3r 
     1                 +hr*w4 + h*w4r )
        pfr = 1.d0 + 2.d0*r*( w1r + er*w2 + e*w2r + gr*w3 + g*w3r 
     1        +hr*w4 + h*w4r )
     2        +r*r*( w1rr + er2*w2 +2.d0*er*w2r + e*w2rr
     3                    + grr*w3 +2.d0*gr*w3r + g*w3rr 
     4                    + hrr*w4 +2.d0*hr*w4r + h*w4rr )
!
        pft = r*r*( w1rt + er*w2t + e*w2rt + grt*w3 + gr*w3t
     1                                     + gt*w3r + g*w3rt 
     2                                     + ht*w4r + h*w4rt )
        uf = -si0t - w1t - e*w2t - gt*w3 - g*w3t -ht*w4 -h*w4t
        cvf = -t*t*( si0tt + w1tt + e*w2tt + gtt*w3 + 2.d0*gt*w3t
     1                                              + g*w3tt 
     2              +htt*w4 + 2.d0*ht*w4t + h*w4tt )
!
       if(abs(f)>0.1E-26) then 
         call sinfun(temk,rho,sin,pn,un,pnr,pnt,cvn,2)
          go to 381
          endif
!
        si = sif
        p = pf
        pr = pfr
        pt = pft
        u = uf
        cv = cvf
       pcalc= p*rhoc*rval*temk
       si= si*rval*temk
       dpdr= pr*rval*temk
       if(iterqk.eq.1)return
        go to 389
!
!
!********************************************************************
!
!       combined functions
!
!*******************************************************************
!
 381    continue
        sin = sin + c000 + c001*dt
        un = un - c001
        si = sif + f*( sin - sif)

! p is pbar = p/(rho(crit)*r*t) where r is the gas constant.
!  in hill's code r is the reduced density.

        p = pf + f*( pn -pf ) + r*r*fr*( sin - sif )
        term1 = pfr
        term2 = f*( pnr - pfr )
        term3 = 2.d0*fr*( pn - pf )
        term4 = ( 2.d0*r*fr + r*r*frr )*( sin -sif )

!   pr is (dpbar/drho) at constant t

        pr = term1 + term2 + term3 + term4
       pcalc= p*rhoc*rval*temk
       si= si*rval*temk
       dpdr= pr*rval*temk
       if(iterqk.eq.1) return
!
!
        term5 = pft
        term6 = f*( pnt -pft ) + ft*( pn - pf )
        term7 = -r*r*fr*( un - uf )
        term8 = r*r*frt*( sin - sif )

!  pt is (dpbar/dt) at constant rho

        pt = term5 + term6 + term7 + term8
        u = uf + f*( un -uf ) - ft*( sin - sif )
        term9 = cvf
        term10 = f*( cvn - cvf )
        term11 = 2.d0*t*t*ft*( un -uf )
        term12 = -t*t*ftt*( sin -sif )
        cv = term9 + term10 + term11 + term12
 389    continue
        h = u - p/r/t
      cp = 1000.d10
        a = 0.d0
      delr = r - 1.d0
      if(dabs(dt).le.1.d-7.and.dabs(delr).le.6.d-4)go to 3891
!      if(((dabs(pr).le.1.d-20).or.(dabs(pr).gt.1.d20)) .or.
!     &  (dabs(t*pt-p) .gt. 1.d15)) then
!       print *,' dabs(pr) d-20 <> d+20 ',dabs(pr)
!       print *,' dabs(t*pt-p) < d+15 ',dabs(t*pt-p)

      if(((dabs(pr).le.1.d-30).or.(dabs(pr).gt.1.d30)) .or.
     &  (dabs(t*pt-p) .gt. 1.d20)) then
                   write(6,9999)temk,rho,pcalc,l
                   return
                   endif
        cp = cv + ( t*pt - p )**2/( r*r*pr )
         a = ( cp/cv*(-pr/t) )
       if(a.le.0.0d0) then
!          print *,' a > 0.0 ',a
                      a=-10.d0
!                      write(6,9999) temk,rho,pcalc,l
! next line commented 12/5/89.  a is dimensionless speed of sound.
!  in the remainder of this subroutine's code further calculations of
!   a are now made conditional on sign of a.
!                      return
                      endif
9999   format('0','error wderiv: cannot calculate a',3e20.10,i3/
     *    ' tk or rho out of range or in unstable supersat region')
       if(a .ge. 0.d0)a= dsqrt(a)
 3891   continue
        s = -si/rval/temk - u*t
!
!
!        restore from dimensionless variables
!
        dpdt = ( pt*rhoc*rval*tc + pcalc )/temk
        u = u*rval*tc
        h = h*rval*tc
        s = s*rval
        cv = cv*rval
        cp = cp*rval
        if(a .ge. 0.d0)a = a*dsqrt( 1000.d0*rval*tc )
!
!
!       calculation additional quantities
!               dhdpt:  joule thomson coeff (isothermal)
!               dtdph:  joule thomson coeff (isenthalpic)
!               dtdps:  pressure temperature coeff (isentropic)
!               expis:  isentropic exponent (rho/p*dpdrs)
!
!               qv(20): array of values used by plot programs
!
!      dpdr= 0.0 at critical point)(
!
       if(abs(dpdr)<1.0E-26) return
       dtdph= 1.0d0/rho-temk/rho**2*dpdt/dpdr
       dtdph= -dtdph/cp
       dtdps=  temk/cp/rho**2*dpdt/dpdr
       dhdpt= dtdph*(-cp)
       expis= rho/pcalc*dpdr*cp/cv
!
       qv(1)= dtdph
       qv(2)= dtdps
       qv(6)= dhdpt
        return
       end
!
!*******************************************************************
!
!       given pressure and temperature find rho etc(l = 3 or 4)
!
!********************************************************************
!
       subroutine wrho(temk,rho,press,si,dpdr,dpdt,u,h,s,cv,cp,
     *                   a,l)
       implicit double precision(a-h,o-z)
!      deleted unused arrays d(3),dd(3)
       double precision w(4,6),c0(8)
       double precision qv(20)
       common/cqv/qv
       common/efgval/eval,dr0,dt0,delta,aval,bval,cval,dval,hval
       common/norms/tc,rhoc,rval
       common/iter/iterqk
       common/osg/thet48,posg,rfosg,rgosg,hfosg,hgosg
       common/sat/psat,rf,rg,hf,hg,sf,sg
       common/czero/c0
       common/cinit/c000,c001
       common/satliq/ctemkf,crhof,cpf,csif,cdpdrf,cdpdtf,cuf,chf,
     *               csf,ccvf,ccpf,caf
       common/satvap/ctemkg,crhog,cpg,csig,cdpdrg,cdpdtg,cug,chg,
     *               csg,ccvg,ccpg,cag
!
       common/carga/dhdpt,dtdph,dtdps,expis,x5,x6,x7,x8,x9,x10
       common/cargf/dhdptf,dtdphf,dtdpsf,expisf,xf5,xf6,xf7,xf8,xf9,xf10
       common/cargg/dhdptg,dtdphg,dtdpsg,expisg,xg5,xg6,xg7,xg8,xg9,xg10
       common/cwcnt/icwcnt
 600    continue
       rhogss=press/temk/.4
       if(rhogss.gt. 1.0) rhogss=1.0
       if((l.eq.4).and.(temk.lt.tc))rhogss=frfwag(temk)
!
!      iteration to rho from press input arg
!
       spinf= wspinf(temk)
       if((temk.ge.311.15).and.(temk.le.321.15))
     &     spinf=spinf*.001 +spinf
       sping= wsping(temk)
!
       do 500 j=1,40
       call wderiv(temk,rhogss,pcalc,si,dpdr,q6,q7,q8,q9,q10,q11,q12,l)
       if(abs(press)>1.0E-26) then
       		if(dabs(1.0d0-pcalc/press).lt.0.1d-8)then
       			go to 510
       			endif
          endif
       if((l.eq.4).and.(dabs(pcalc-press).lt.0.1d-8))go to 510
!
!      check for dpdr negative
!
       if((temk.ge.tc).and.(dpdr.lt.0.0d0))go to 599
       if(temk.ge.tc) go to 520
       if((l.eq.3).and.(rhogss.lt.sping).and.(dpdr.ge.0.0d0))
     &    go to 520
       if((l.eq.4).and.(rhogss.gt.spinf).and.(dpdr.gt.0.0d0))
     &    go to 520
!
       if((l.eq.3).and.(rhogss.gt.sping))rhogss=sping
       if((l.eq.4).and.(rhogss.le.spinf))rhogss=spinf
       if(l.eq.3) rhogss=rhogss*.98
       if(l.eq.4) rhogss=rhogss*1.02
       go to 500
520    continue
!
       dpdx= dpdr*1.1
       dpdlim= .1
      if(((press.ge.22.046d0).and.(press.le.23.))
     *  .and.((temk.ge.647.067d0).and.(temk.lt.678.15d0)))
     *     dpdlim= .1d-8
       if(dpdx.lt.dpdlim) dpdx=dpdlim
       x=(press-pcalc)/dpdx
       if(dabs(x) .gt. .1d0)x= x*.1/dabs(x)
       rhogss=rhogss+x
       if(rhogss.le..1d-12)rhogss=.1d-12
500    continue
!
599    continue
       write(6,9999)l,temk,press,pcalc,rhogss,dpdr
       rho=-10.
9999   format('wrho fail convg:l,temk,press,pcalc,rgss: ',
     *   'dpdr: '/i2,f7.2,4e20.10)
       write(6,9996)
9996   format(' ','is l set correctly?')
       return
!
!
510    continue
       rho=rhogss
       icwcnt=j
       if(l.eq.5) return
       if(l.eq.1)return
!
!
!      (if l=1, or l=5, only require calc of rho,si,dpdr from press)
!
       iterqk=0
       call wderiv(temk,rho,pcalc,si,dpdr,dpdt,u,h,s,cv,cp,a,l)
       return
       end
       
       
      subroutine wsat(temk,rho,press,si,dpdr,dpdt,u,h,s,cv,cp,
     *                   a,l)
      implicit double precision(a-h,o-z)
!     deleted unused arrays d(3),dd(3)
      double precision w(4,6),c0(8)
      double precision qv(20)
      common/cdel/cdelta
      common/ccread/cread
      common/cerr/cerrtk,cerrho,cerrp,cerrl,cerrls
      common/cqv/qv
      common/efgval/eval,dr0,dt0,delta,aval,bval,cval,dval,hval
      common/norms/tc,rhoc,rval
      common/iter/iterqk
      common/osg/thet48,posg,rfosg,rgosg,hfosg,hgosg
      common/sat/psat,rf,rg,hf,hg,sf,sg
      common/czero/c0
      common/cinit/c000,c001
      common/satliq/ctemkf,crhof,cpf,csif,cdpdrf,cdpdtf,cuf,chf,
     *               csf,ccvf,ccpf,caf
      common/satvap/ctemkg,crhog,cpg,csig,cdpdrg,cdpdtg,cug,chg,
     *               csg,ccvg,ccpg,cag
!
      common/carga/dhdpt,dtdph,dtdps,expis,x5,x6,x7,x8,x9,x10
      common/cargf/dhdptf,dtdphf,dtdpsf,expisf,xf5,xf6,xf7,xf8,xf9,xf10
      common/cargg/dhdptg,dtdphg,dtdpsg,expisg,xg5,xg6,xg7,xg8,xg9,xg10
      common/cscnt/icscnt
!
!********************************************************************
!
!       given temperature find saturation states(l=5)
!
!********************************************************************
!
      if(temk.gt.tc) then
       	write(6,705)temk,tc
705     format('0','wsat: input temk ',f10.3,' above tc ',f10.3,
     *             ' for l=5')
      	return
      endif
!
!
!
      time = -tc/temk
      beta= .325d0
      dt = time + 1.d0
      cf = 0.07857d0
      cg = 0.07805d0
      tcdiff = tc -temk
      if(abs(tcdiff)<1.0E-26) then
          rf=rhoc
          rg=rhoc
          go to 760
      endif
      if(tcdiff.le.0.0001d0) then
          rf = rhoc + cf*tcdiff**beta
          rg = rhoc - cg*tcdiff**beta
      		go to 760
      endif
!
!      iteration to saturation (gdiff= 0.0)
!
       rf= frfwag(temk)
       rg= frgwag(temk)
       !print *, rf ,rg

       do 750 jsat=1,50
       			rho=rf
       			call wderiv(temk,rho,pf,sif,dpdrf,q6,q7,q8,q9,q10,q11,q12,2)

       			rho=rg
       			call wderiv(temk,rho,pg,sig,dpdrg,q6,q7,q8,q9,q10,q11,q12,2)
						!print *, pf ,pg
       			dp=pg-pf
       			dg= (sig+pg/rg) -(sif+pf/rf)

!							Zmniejszono tolerancję z e-12 do e-11
       			if(dabs(dg) .lt. 0.1d-11*rval*temk) go to 760

       			dp= -dp
       			dg= -dg
       			aa= (1.d0/rg) - (1.d0/rf)
       			drg= (dg-dp/rf)/(aa*dpdrg)
       			drf= (dg-dp/rg)/(aa*dpdrf)
       			!print *, rf, rg, drg, drf
       			rg= rg+drg
       			rf= rf+drf
       			!print *, rf, rg
       			
750    continue

       write(6,999) temk,pf,pg,dp,dg
999    format('wsat non-convg: tk,pf,pg,dp,dg; ',5e15.8)

760    continue
			 !print *, "Po zbieglo", temk,pf,pg,dp,dg

       icscnt= jsat
       iterqk = 0
       rho = rg
       call wderiv(temk,rho,press,si,dpdr,dpdt,u,h,s,cv,cp,a,l)
       psat = press
       rg = rho
       hg = h
       sg = s
       ag = a
       ctemkg=temk
       crhog=rho
       cpg=press
       csig=si
       cdpdrg=dpdr
       cdpdtg=dpdt
       cug=u
       chg=h
       csg=s
       ccvg=cv
       ccpg=cp
       cag=a
       dhdptg=dhdpt
       dtdphg=dtdph
       dtdpsg=dtdps
       expisg=expis
!
       rho = rf
       call wderiv(temk,rho,pcalcf,si,dpdr,dpdt,u,h,s,cv,cp,a,l)
       rf = rho
       hf = h
       sf = s
       af = a
       ctemkf=temk
       crhof=rho
       cpf=press
       csif=si
       cdpdrf= dpdr
       cdpdtf=dpdt
       cuf=u
       chf=h
       csf=s
       ccvf=cv
       ccpf=cp
       caf=a
       dhdptf=dhdpt
       dtdphf=dtdph
       dtdpsf=dtdps
       expisf=expis
        return
        end
!
!********************************************************************
!********************************************************************
!
        subroutine wfuncs( t,r,w,m )
!
!********************************************************************
!
!       options
!
!       m = 1       initialize
!       m = 2       obtain w1, w2, w3, and all their derivatives
!                                     ( iterqk = 0 )
!                       or w1, w2, w3, and w1r, w2r, w3r
!                                     ( iterqk = 1 )
!
!********************************************************************
!
        implicit double precision(a-h,o-z)
!
        double precision w(4,6),c1(10,10),c2(10,12),c3(10,10),c4(10,10),
     1         r1(10),r1r(10),r1rr(10),r2(10),r2r(10),r2rr(10),
     2         t1(12),t1t(12),t1tt(12),t2(12),t2t(12),t2tt(12),
     3         r3(10),r3r(10),r3rr(10),r4(10),r4r(10),r4rr(10)
       common/ccread/cread
       common/cerr/cerrtk,cerrho,cerrp,cerrl,cerrls
       common /cvcc/c1,c2,c3,c4
       common/norms/tc,rhoc,rval
       common/iter/iterqk
       common/efgval/eval,dr0,dt0,delta,aval,bval,cval,dval,hval

! 		 double precision, dimension (:), allocatable :: r1,r1r,r1rr,r2,
!    *   r2r,r2rr,t1,t1t,t1tt,t2,t2t,t2tt,r3,r3r,r3rr,r4,r4r,r4rr
! 			allocate (r1(10),r1r(10),r1rr(10),r2(10),r2r(10),r2rr(10),
!    *         t1(12),t1t(12),t1tt(12),t2(12),t2t(12),t2tt(12),
!    *         r3(10),r3r(10),r3rr(10),r4(10),r4r(10),r4rr(10))
     

     
     
       data tsave/0.0d0/
			 nr1 = 7
       nt1 = 7
       nr2 = 7
       nt2 = 12
       nr3 = 5
       nt3 = 5
       nr4 = 5
       nt4 = 10
       nrmax = 8
       ntmax = 12
       go to (1,200),m
 1     if( abs(cread-1.d0)>1.0E-26)return
!
!       n.b. nrmax must be .ge. nr1, nr2, and nr3
!            ntmax must be .ge. nt1, nt2, and nt3
!            otherwise rx(i) and tx(j) will be wrong

!
!*****************************************************************
!
!      read non-dimensional coefficients
!
!********************************************************************
      write(6,101)
      do 2 j = 1,nt1
!       Coef132 is not used
!     		if(cread.eq.1.d0)read(10,100) (c1(i,j),i=1,nr1)
      		do 120 i=1,nr1
      				if( abs(c1(i,j))< 1.0E-26 ) go to 120
          		write(*,105) i, j, c1(i,j)
 120  		continue
 2    continue
      write(*,102)
      do 3 j = 1,nt2
!       Coef132 is not used
!     		if(cread.eq.1.d0)read(10,100) (c2(i,j),i=1,nr2)
      		do 130 i=1,nr2
      				if( abs(c2(i,j))< 1.0E-26 ) go to 130
        			write(*,105) i, j, c2(i,j)
 130  		continue
 3    continue
      write(*,103)
      do 4 j = 1,nt3
!       Coef132 is not used
!     		if(cread.eq.1.d0)read(10,100) (c3(i,j),i=1,nr3)
      		do 140 i=1,nr3
      				if( abs(c3(i,j))< 1.0E-26 ) go to 140
          		write(*,105) i, j, c3(i,j)
 140  		continue
 4    continue
      write(*,103)
      do 5 j = 1,nt4
!       Coef132 is not used
!     		if(cread.eq. 1.d0)read(10,100) (c4(i,j),i=1,nr4)
      		do 150 i=1,nr4
      				if( abs(c4(i,j))< 1.0E-26 ) go to 150
      				write(*,105) i, j, c4(i,j)
 150  		continue
 5    continue
 100    format(1x,10d16.10)
 101    format(1x,'   i   j        c1(i,j)  ')
 102    format(1x,'   i   j        c2(i,j)  ')
 103    format(1x,'   i   j        c3(i,j)  ')
 105    format(1x,i4,i4,3x,d16.10)
!
        return
!
!*******************************************************************
!
!       evaluate w functions
!
!********************************************************************
 200  call efunc( r,e,er,er2 )
      eterm = 1.d0 - e
      if( r .lt. 0.1d-4 ) eterm = eval*r*r
      do 34 i=1,nrmax
       zi=i
      if(i .eq. 2) go to 31
      r1(i) = eterm*r**(i-2)
      r1r(i) = eterm*(zi-2.d0)*r**(i-3) - er*r**(i-2)
      r1rr(i) = eterm*((zi-2.d0)*(zi-3.d0))*r**(i-4)
     1         -2.d0*er*(zi-2.d0)*r**(i-3)
     2         -er2*r**(i-2)
      go to 32
 31   continue
!
       if(r.le.0.0d0)write(6,999)cerrtk,cerrho
999    format('0','wfuncs r .le. 0.0; tk,rho: ',2e25.15)
      dlogr= dlog(r)
!
      r1(i) = (1.d0 - e)*dlogr -eval*r*r*dlogr + eval*r*r*0.5d0
      r1r(i) = (1.d0 - e)/r - er*dlogr - 2.d0*eval*r*dlogr
      r1rr(i) = - (1.d0 - e)/r/r - 2.d0*er/r - er2*dlogr
     1          - 2.d0*eval*dlogr - 2.d0*eval
 32   continue
!
      r2(i) = r**(i)
      r2r(i) = zi*r**(i-1)
      r2rr(i) = (zi*(zi-1.d0))*r**(i-2)
!
      r3(i) = r**(i+1)
      r3r(i) = (zi+1.d0)*r**(i)
      r3rr(i) = ((zi+1.d0)*zi)*r**(i-1)
      r4(i) = r1(i)
      r4r(i) = r1r(i)
      r4rr(i)= r1rr(i)
 34   continue
!
!     if(abs(t-tsave)<1.0E-26) go to 250

      do 37 j = 1,ntmax
      		t1(j) = t**(j-1)
      		t1t(j) = dfloat(j-1)*t**(j-2)
        	t1tt(j) = dfloat((j-1)*(j-2))*t**(j-3)
      		t2(j) = t**(j+1)
      		t2t(j) = dfloat(j+1)*t**(j)
      		t2tt(j) = dfloat((j+1)*j)*t**(j-1)
 37   continue
 			
        tsave = t
 250    continue

        do 301 i = 1,4
        		do 300 j = 1,6
        				w(i,j) = 0.d0
 300    		continue
 301    continue


        do 321 i = 1,nr1
        		do 320 j = 1,nt1
        				if( abs(c1(i,j))< 1.0E-26 ) go to 319
        				w(1,1) = w(1,1) + c1(i,j)*r1(i)  *t1(j)
        				w(1,2) = w(1,2) + c1(i,j)*r1r(i) *t1(j)
        				w(1,3) = w(1,3) + c1(i,j)*r1rr(i)*t1(j)
        				if( iterqk .eq. 1 ) go to 319
        				w(1,4) = w(1,4) + c1(i,j)*r1(i)  *t1t(j)
        				w(1,5) = w(1,5) + c1(i,j)*r1(i)  *t1tt(j)
        				w(1,6) = w(1,6) + c1(i,j)*r1r(i) *t1t(j)
 319    		continue
 320    		continue
 321    continue
!
        do 331 i = 1,nr2
        		do 330 j = 1,nt2
        				if( abs(c2(i,j))< 1.0E-26 ) go to 329
        				w(2,1) = w(2,1) + c2(i,j)*r2(i)  *t1(j)
        				w(2,2) = w(2,2) + c2(i,j)*r2r(i) *t1(j)
        				w(2,3) = w(2,3) + c2(i,j)*r2rr(i)*t1(j)
        				if( iterqk .eq. 1 ) go to 329
        				w(2,4) = w(2,4) + c2(i,j)*r2(i)  *t1t(j)
        				w(2,5) = w(2,5) + c2(i,j)*r2(i)  *t1tt(j)
        				w(2,6) = w(2,6) + c2(i,j)*r2r(i) *t1t(j)
 329    		continue
 330    		continue
 331    continue

        do 341 i = 1,nr3
        		do 340 j = 1,nt3
        				if( abs(c3(i,j))< 1.0E-26 ) go to 339
        				w(3,1) = w(3,1) + c3(i,j)*r3(i)  *t2(j)
        				w(3,2) = w(3,2) + c3(i,j)*r3r(i) *t2(j)
        				w(3,3) = w(3,3) + c3(i,j)*r3rr(i)*t2(j)
        				if( iterqk .eq. 1 ) go to 339
        				w(3,4) = w(3,4) + c3(i,j)*r3(i)  *t2t(j)
        				w(3,5) = w(3,5) + c3(i,j)*r3(i)  *t2tt(j)
        				w(3,6) = w(3,6) + c3(i,j)*r3r(i) *t2t(j)
 339    		continue
 340    		continue
 341    continue
!
        do 351 i = 1,nr4
        do 350 j = 1,nt4
        if( abs(c4(i,j))< 1.0E-26  ) go to 349
        w(4,1) = w(4,1) + c4(i,j)*r4(i)  *t1(j)
        w(4,2) = w(4,2) + c4(i,j)*r4r(i) *t1(j)
        w(4,3) = w(4,3) + c4(i,j)*r4rr(i)*t1(j)
        if( iterqk .eq. 1 ) go to 349
        w(4,4) = w(4,4) + c4(i,j)*r4(i)  *t1t(j)
        w(4,5) = w(4,5) + c4(i,j)*r4(i)  *t1tt(j)
        w(4,6) = w(4,6) + c4(i,j)*r4r(i) *t1t(j)
 349    continue
 350    continue
 351    continue
  
! 			deallocate (r1,r1r,r1rr,r2,r2r,r2rr,t1,t1t,t1tt,t2,t2t,
!    1         t2tt,r3,r3r,r3rr,r4,r4r,r4rr) 
        return
        end
!
!********************************************************************
!********************************************************************
!
        subroutine efunc( r,e,er,er2 )
!
!*******************************************************************
!
        implicit double precision(a-h,o-z)
        common/efgval/eval,dr0,dt0,delta,aval,bval,cval,dval,hval
!
        arge = eval*r*r
        if( arge .gt. 170.d0 ) arge = 170.d0
        e = ddexp( -arge )
        er = -2.d0*eval*r*e
        er2 = eval*( -2.d0 + 4.d0*eval*r*r )*e
!
        return
        end
canames
!********************************************************************
!
        subroutine ffunc( r,t,f,fr,frr,ft,ftt,frt )
!
!*******************************************************************
!
        implicit double precision(a-h,o-z)
       common/fxxx/fxval
        common/efgval/eval,dr0,dt0,delta,aval,bval,cval,dval,hval
       dr = r - 1.d0
       dt = t + 1.d0
       s = dsqrt( (dr/dr0)**2 + (dt/dt0)**2)
      if( s .lt. 0.003d0 ) s = 0.003d0
!
       power = 4.d0
       pow1  = (power - 1.d0)/power
       pow2  = (power - 2.d0)/power
       pow3  = power/(power - 1.d0)
       pow4  = power*(power - 1.d0)
       args  = (s/delta)**power
       if(args .gt. 23.5d0) args = 23.5d0
       caps  = ddexp( args ) - 1.d0
       argc  = 1.d0/caps
       if( argc .gt. 35.d0 ) argc = 35.d0
       exps  = ddexp( -argc )
       fs1   = dlog( caps + 1.d0 )
       fs2   = (1.d0 - caps - caps*caps )/caps/caps
       fs3   = (1.d0 + caps )/caps/caps
!
!
       f     = 1.d0 - exps
       if( f .lt. 0.1d-9 ) f = 0.d0
       dfds  = -exps/delta*fs3*power*fs1**pow1
       if( dabs(dfds) .lt. 0.1d-12 ) dfds = 0.d0
       d2fds2= -exps/delta/delta*fs3*pow4*fs1**pow2
     1         *(1.d0 + fs2*pow3*fs1)
       if( dabs(d2fds2) .lt. 0.1d-12 ) d2fds2 = 0.d0
!
       dsdr = dr/s/dr0/dr0
       dsdt = dt/s/dt0/dt0
       dtfac = (dt/dt0)**2
       drfac = (dr/dr0)**2
       d2sdt2 = (1.d0/s - dtfac/(s)**3)/dt0/dt0
       d2sdr2 = (1.d0/s -drfac/(s)**3)/dr0/dr0
       d2sdtr = -dt*dr/(s*dt0*dr0)**2/s
!
       ft = dfds*dsdt
       fr = dfds*dsdr
       ftt = d2fds2*(dsdt)**2 + dfds*d2sdt2
       frr = d2fds2*(dsdr)**2 + dfds*d2sdr2
       frt = dfds*d2sdtr +d2fds2*dsdt*dsdr
       if( dabs(ft) .lt. 0.1d-12 ) ft = 0.d0
       if( dabs(fr) .lt. 0.1d-12 ) fr = 0.d0
       if( dabs(ftt).lt. 0.1d-12 ) ftt = 0.d0
       if( dabs(frr).lt. 0.1d-12 ) frr = 0.d0
       if( dabs(frt).lt. 0.1d-12 ) frt = 0.d0
       fxval=f
        return
        end
!
!********************************************************************
!
        subroutine gfunc( r,t,g,gr,grr,gt,gtt,grt )
!
!********************************************************************
!
        implicit double precision(a-h,o-z)
        common/efgval/eval,dr0,dt0,delta,aval,bval,cval,dval,hval
!
        dr = r - 1.d0
        dt = t + 1.d0
        dt2=dt*dt
        dr2=dr*dr
        a2=aval*aval
        b2=bval*bval
        s=dsqrt(a2*dt2+b2*dr2)
      if( s .lt. 0.0001d0 ) s = 0.0001d0
        scu=s*s*s
        top1 = a2*dt 
        top2 = b2*dr 
        dsdt=top1/s
        dsdr=top2/s
        d2sdt=a2/s - top1**2/scu
        d2sdr=b2/s - top2**2/scu
        dsdrt= -top1*top2/scu
!
!
!
        gfac = -aval*dt - bval*dr -cval*dt2 - dval*dr2
        if(dabs(gfac) .gt. 100.d0) gfac=-100.d0
        g = ddexp(gfac)
!
!
!
        t1 = -aval - 2.d0*cval*dt
        t2 = -bval - 2.d0*dval*dr
        gt=g*t1
        gr=g*t2
!
!
!
         gtt = gt*t1 - 2.d0*cval*g
         grr = gr*t2 - 2.d0*dval*g
         grt = gr*t1
      if( dabs(g)  .le. 0.1d-12 ) g = 0.d0
      if( dabs(gr) .le. 0.1d-12 ) gr = 0.d0
      if( dabs(grr).le. 0.1d-12 ) grr= 0.d0
      if( dabs(gt) .le. 0.1d-12 ) gt = 0.d0
      if( dabs(gtt).le. 0.1d-12 ) gtt= 0.d0
      if( dabs(grt).le. 0.1d-12 ) grt= 0.d0
        return
        end
!
       double precision function ddexp(x)
       implicit double precision (a-h,o-z)
!
!      threshold to set dexp function to zero
!        for large negative arguments  (-35)
!
         double precision x
       if(x.le.-35.) ddexp= 0.0d0
       if(x.le.-35.) return
       ddexp= dexp(x)
         return
         end
!
!
!*********************************************************************
!
         subroutine hfunc( r,t,h,hr,hrr,ht,htt,hrt )
!
!*********************************************************************
!
          implicit double precision(a-h,o-z)
        common/efgval/eval,dr0,dt0,delta,aval,bval,cval,dval,hval
      argh = hval*(t + 3.d0)
      h = ddexp(-argh )
      hr = 0.d0
      hrr = 0.d0
      ht = -hval*h
      htt = hval*hval*h
      hrt = 0.d0
         return
         end
!**********************************************************************
!
        subroutine si0f( t,si0,si0t,si0tt,l )
!
!********************************************************************
!
        implicit double precision(a-h,o-z)
        common/norms/tc,rhoc,rval
      common/czero/c0
        double precision c0(8)
!
!
 100    si0 = 0.d0
        si0t = 0.d0
        si0tt = 0.d0
        do 110 i = 1,6
        power = ( -t )**(2-i)
        si0 = si0 + c0(i)*power
        si0t = si0t - c0(i)*dfloat(2-i)*power/(-t)
        si0tt = si0tt + c0(i)*dfloat(2-i)*dfloat(1-i)*power/t/t
 110    continue
        si0 = si0 + ( c0(7)*t + c0(8) )*dlog(-t)
        si0t = si0t + c0(7)*( dlog(-t) + 1.d0 ) +c0(8)/t
        si0tt = si0tt + c0(7)/t - c0(8)/t/t
        return
        end
!
!****************************************************************
!****************************************************************
!****************************************************************
!****************************************************************
!****************************************************************
!****************************************************************
!
!     subroutines for the revised and extended scaling solution
!         to calculate values near the critical point
!
!         ref:  j.m.h. levelt sengers, b. kamgar-parsi,
!               f.w. balfour, j.v. sengers
!               "thermodynamic properties of steam in the 
!          critical region", j.phys.chem.ref.data vol 12,p 1(1983)
!
!*****************************************************************
!
      subroutine sinfun(t,dgmpcc,si,p,u,dpdr,dpdt,cv,l)
!
!
!      given the temperature t(k) and density d(gm/cc), this routine
!     calculates pressure p(mpa), energy u and enthalpy h(kj/kg),
!     entropy, specific heats cp and cv(kj/kg-k), velocity of sound
!     cs(m/s), and compressibility comp(1/mpa).
!
      implicit double precision (a-h,o-z)
      character(len=6) anames,qnames
      common/coefs/a(20),q(20),anames(20),qnames(20)
      common/crits/tc,rhoc,pc,pcon,ucon,scon,rval
      common/secder/d2pdm2
      common/ccmeta/cunst
      dimension s(2),xk(2),sd(2)
      equivalence (pw1,a(5)),(pw2,a(4)),(pw3,a(2)),
     1   (am0,a(13)),(am1,a(14)),(am2,a(15)),(am3,a(16)),
     2   (p00,q(11)),(p20,q(12)),(p40,q(13)),
     3   (p01,q(18)),(p21,q(19)),(p41,q(20)),
     4   (aa,a(10)),(xk0,a(7)),(xk1,a(12)),(pw11,q(9)),
     5   (alpha,q(10)),(alhi,q(15)),(besq,a(9))
!
        data inival/0/
!
       if(l.eq.0)return
      if  (inival .ne. 0)  go to 50
      call congen
      q(3)   = a(8) * q(7)
      tc     = q(4) + a(8)
      rhoc   = q(5) + a(3)
      pc     = q(6) + q(3)
      pcon   = pc / tc
      ucon   = 1.d3 * pc
      scon   = ucon / tc
      dpcon  = pcon / rhoc / rhoc
      rhocc  = 0.322778d0
      rval   = 0.46151d0
       inival=1
!
!
 50   continue
      d  = dgmpcc * 1000.d0
      xk(1)=xk0
      xk(2)=xk1
      tee=(t-tc)/tc
      tw=-tc/t
      dtw=1.d0+tw
      rho=d/rhoc
!
      delr = rho - 1.d0
      if(dabs(tee).le.1.d-7.and.dabs(delr).le.0.0006d0)go to 55
      dtstin= 1.0d0 - (tc/t)
      call consim(rho,tee,amu,th1,r1,rho1,s,rhodi,unstab,dtstin)
!      write(9,9099) th1,r1
9099   format(' ','theta,r: ',2e20.8)
      deltmu = amu
      tt1=th1*th1
      tt2=tt1*tt1
      call ss(r1,th1,s,sd)
      pw0=1.d0+dtw*(pw1+dtw*(pw2+dtw*pw3))
      pwmu=amu*rhodi
!
!       pwmu = dmu * (1 + p11*dt)
!
      p0th  = p00+p20*tt1+p40*tt2
      p1th  = p01+p21*tt1+p41*tt2
      dpw0  = xk0*p0th*r1**(2.d0-alpha)
      dpw1  = xk1*p1th*r1**(2.d0-alhi)
      dpw   = aa*(dpw0+dpw1)
      pw    = pw0+pwmu+dpw
      p     = pw * pcon * t
      dp0dt = pw1+dtw*(2.d0*pw2+3.d0*pw3*dtw)
      dm0dt = am1+dtw*(2.d0*am2+3.d0*am3*dtw)
      uw    = dp0dt- rho * dm0dt + pw11 * amu+s(1)+s(2)
      hw    = pw - tw * uw
      amw   = amu + am0 + dtw*(am1+dtw*(am2+dtw*am3))
      sw    = hw - rho * amw
      u     = uw * ucon / d
      scond = scon / d
      entrop= sw * scond
      si    = u - t * entrop
      pdimen= p
      p     = p / rhocc / rval / t
      u     = u / rval / tc
      si    = si / rval / t
      if  (l .eq. 1)  return
!
!
      d2p0dt= 2.d0*pw2+6.d0*pw3*dtw
      d2m0dt= 2.d0*am2+6.d0*am3*dtw
      call aux(r1,th1,d2pdt2,d2pdmt,d2pdm2,aa,xk,sd,cvcoex)
      dpdtcd=dp0dt+pw11*(amu-rho/d2pdm2)+s(1)+s(2)-d2pdmt*rho/d2pdm2
      dpwdtw= pw-tw*dpdtcd
      cvitw2= d2p0dt-rho*d2m0dt+d2pdt2-(pw11+d2pdmt)**2/d2pdm2
      dpdd  = dpcon * d * t / d2pdm2
      dpdr  = 1000.d0 * dpdd
!     write(6,101) t,d,pcon,dpcon,d2pdm2,dpdd,dpdr
      dpdt  = pcon * dpwdtw
      cvw   = cvitw2 * tw * tw
      cv    = cvw * scond
      dpdr  = dpdr / rval / t * unstab
!     write(6,101) t,d,dpcon,d2pdm2,dpdd,dpdr
c101  format(1x,'t,d,pcon,dpcon,d2pdm2,dpdd,dpdr ',7d14.6)
      dpdt  = (dpdt * t - pdimen) / rhocc / rval / tc
      cv    = cv / rval
!
!
      return
 55   continue
      p = 22.0460d0/rval/tc/rhocc
      u = 2014.82d0 - 468.03d0*delr
      u = u/rval/tc
      ent = 4.40498d0 - 0.8231d0*delr
      si = u - tc*ent/rval/tc
      dpdr = 0.d0
      dpdt = 0.d0
      cv = 1000.d0
      return
      end
       block data clsblk
      character(len=6) anames,qnames
      double precision a,q
       common/coefs/a(20),q(20),anames(20),qnames(20)
!
!     this subroutine supplies the parameters used in the equation
!     of state.
!
      data a/-.017762d0,5.2380d0,0.d0,-25.4915d0,6.8445d0,.325d0,1.4403
     1d0,0.d0,1.3757d0,23.6666d0,4.82d0,.2942d0,-11.2317d0,
     2-22.6407d0,-17.7876d0,-4.9332d0,0.d0,0.d0,0.d0,0.d0/
      data q/-.006d0,-.003d0,0.d0,647.067d0,322.778d0,22.0460d0
     1,.267d0,-1.6d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0
     2,0.d0,0.d0,0.d0,0.d0/
      data anames/'c' ,'cot3' ,'delroc' ,'cot2' ,'dpdt' ,'beta',
     1            'k0' ,'deltc' ,'b*b' ,'a' ,'delta' ,'k1' ,
     2            'muc' ,'mu1' ,'mu2' ,'mu3' ,'s00' ,'s20' ,
     3            's01' ,'s21' /
      data qnames/2*'blank'     ,'delpc' ,'tc'    ,'rhoc'  ,'pc'    ,
     1              'dpcdtc','slopdi','p11'   ,'alpha' ,'p00'   ,
     2              'p20'   ,'p40'   ,'deltai','alphai','betai',
     3              'gammai','p01'   ,'p21'   ,'p41'   /
      end
      subroutine congen
!
!   this subroutine calculates all quantities not dependent on
!  r or theta
!
      implicit double precision (a-h,o-z)
			character(len=6) anames,qnames
      common /coefs/ a(20),q(20),anames(20),qnames(20)
!
!
      equivalence (deli,q(14)),(drhoc,a(3)),(beta,a(6)),(delta,a(11))
     1,(cotee2,a(4)),(dpdt,a(5)),(deltc,a(8)),(besq,a(9))
     2,(cc,a(1)),(aa,a(10)),(tck,q(4)),(rhoc,q(5)),(pc,q(6)),(xk0,a(7))
     3,(alhi,q(15)),(beti,q(16)),(gami,q(17)),(p01,q(18)),(p21,q(19))
     4,(alpha,q(10)),(p00,q(11)),(p20,q(12)),(p40,q(13)),(p41,q(20))
     5,(s00,a(17)),(s20,a(18)),(s01,a(19)),(s21,a(20))
!
!
      data delpp,delpm,delpph,delpmh/1.d-6,-1.d-6,1.000001d0,.999999d0/
      alpha = 2.d0 - a(6)*(a(11) + 1.d0   )
      if  (dabs(alpha) .lt. 1.d-6)   alpha=1.d-6
      gamma = beta*(delta - 1.d0)
      deli = 0.50d0
      if (deli .lt. -beta)   deli=-beta
      alhi = alpha - deli
      beti = beta + deli
      gami = gamma - deli
      if(alhi.le.delpph.and.alhi.ge.delpmh) alhi=1.000001d0
      if (alhi.le.delpp.and.alhi.ge.delpm) alhi = .000001d0
      delpc = q(7)*deltc
      er2 = 2.d0*beta*delta - 1.d0
      p00 = (beta*(delta-3.d0)-besq*alpha*gamma)/(2.d0*
     1besq*besq*(2.d0-alpha)*(1.d0-alpha)*alpha)
      p20=-(beta*(delta-3.d0)-besq*alpha*er2)
     1 / (2.d0 * besq * (1.d0 - alpha) * (alpha))
      p40=(er2-2.d0)/2.d0/alpha
      s00 = (2.d0-alpha)*p00
      s20 = -beta*(delta-3.d0)/2.d0/besq/alpha
      da=q(1)
      db=q(2)
      ra = da/(1.d0 - besq)
      rb = db/(1.d0 - besq)
      sw0 = s00+s20
      dr0=sw0*(ra**(1.d0-alpha)-rb**(1.d0-alpha))
      p01 = (beta*(delta-3.d0)-3.d0*deli-a(9)*alhi*gami)
     1 /(2.d0*a(9)*a(9)*(2.d0-alhi)*(1.d0-alhi)*alhi)
      p21 = -((beta*(delta-3.d0)-3.d0*deli-a(9)*alhi*er2)/
     1 (2.d0*a(9)*(1.d0-alhi)*alhi))
      p41 = (.5d0*er2 -1.d0)/alhi
      s01 = (2.d0-alhi)*p01
      s21 = -(beta*delta-3.d0*beti)/2.d0/besq/alhi
      ptw = p01 + p21 + p41
      sw1 = s01+s21
      dr1=sw1*(ra**(1.d0-alhi)-rb**(1.d0-alhi))
      if(xk0.le.0.d0) xk0=0.0001d0
      d1=cc*aa*xk0*dr0
      d2=cc*aa*a(12)*dr1
      q(9) = (q(8)*(1.d0/(1.d0-db)-1.d0/
     1   (1.d0 - da)) + d1 + d2) / (db - da)
      return
      end
      subroutine ss(r,th,s,sd)
      implicit double precision (a-h,o-z)
      character(len=6) anames,qnames
      dimension s(2),sd(2)
      common /coefs/ a(20),q(20),anames(20),qnames(20)
      equivalence (alpha,q(10)),(beta,a(6)),(besq,a(9)),(delta,a(11))
     1,(deli,q(14)),(alhi,q(15)),(beti,q(16)),(gami,q(17)),(p00,q(11))
     2,(p01,q(18)),(s00,a(17)),(s20,a(18)),(s01,a(19)),(s21,a(20))
      tt = th*th
      s(1) = s00 + s20*tt
      sd(1) = 2.d0*s20*th
      s(2) = s01 + s21*tt
      sd(2) = 2.d0*s21*th
      s(1)=s(1)*a(10)*a(7)*r**(1.d0-alpha)
      s(2)=s(2)*a(10)*a(12)*r**(1.d0-alhi)
      return
      end
      subroutine aux(r1,th1,d2pdt2,d2pdmt,d2pdm2,aa,xk,sd,cvcoex)
!
!     this routine calculates second derivatives of the scaled
!     equation of state
!
      implicit double precision (a-h,o-z)
      character(len=6) anames,qnames
      common /coefs/ a(20),q(20),anames(20),qnames(20)
      dimension xk(2),s(2), sd(2),w(2),y(2),z(2),coex(2)
      equivalence (cc,a(1)),(beta,a(6)),(besq,a(9)),(delta,a(11))
     1,(alpha,q(10)),(s00,a(17)),(s20,a(18)),(s01,a(19)),(s21,a(20))
      deli = 0.d0
      s(1)=s00+s20*th1*th1
      s(2)=s01+s21*th1*th1
      sd(1) = 2.d0*th1*s20
      sd(2) = 2.d0*th1*s21
      ww = 0.d0
      yy = 0.d0
      zz = 0.d0
      gamma = beta*(delta -1.d0)
      tt1 = th1*th1
      ter = 2.d0*beta*delta - 1.d0
      g = (1.d0+(besq*ter-3.d0)*tt1 - besq*(ter-2.d0)*tt1*tt1)
      cvcoex = 0.d0
      do 30 i = 1,2
      alhi = alpha - deli
      beti = beta + deli
      gami = gamma - deli
      w(i)=(1.d0-alhi)*(1.d0-3.d0*tt1)*s(i)-
     1  beta * delta * (1.d0 - tt1) * th1 * sd(i)
      w(i) = (w(i) *(r1**(-alhi)))/g
      w(i) = w(i)*xk(i)

      ww = ww + w(i)
      y(i) = beti*(1.d0-3.d0*tt1)*th1 - beta*delta*(1.d0-tt1)*th1
      y(i) = (y(i) * (r1**(beti - 1.d0)))*xk(i)/g
      yy = yy + y(i)
      z(i) = 1.d0 - besq*(1.d0 - (2.d0*beti))*tt1
      z(i) = (z(i) *(r1**(-gami)))*xk(i)/g
      zz = zz + z(i)
      a1 = (beta*(delta-3.d0)-3.d0*deli-besq*alhi*gami)
     1/(2.d0*besq*besq*(2.d0-alhi)*(1.d0-alhi)*alhi)
      a2 = 1.d0+((beta*(delta-3.d0)-3.d0*deli-besq*alhi*ter)/
     1(2.d0*besq*(1.d0-alhi)*alhi))
      a2 = - a2
      a4 = 1.d0+((ter-2.d0)/(2.d0*alhi))
      f1 = a1 + a2 + a4
      coex(i) = ((2.d0 - alhi)*(1.d0 - alhi)*(r1**(-alhi))*f1*xk(i))
      cvcoex = cvcoex + coex(i)
      deli = 0.5d0
 30   continue
      d2pdt2 = aa*ww
      d2pdmt = yy + aa*cc*ww
      d2pdm2 = zz/aa + 2.d0*cc*yy + (cc**2)*aa*ww
      return
      end
       subroutine consim(rm,t,amu,theta,r,rho1s,s1,rhodi,unstab,dtstin)
!
!      this routine converts calls in old conver form to
!                                     new conver form
!
       implicit double precision (a-h,o-z)
       character(len=6) anames,qnames
       common/cerr/cerrtk,cerrho,cerrp,cerrl,cerls
       common/coefs/a(20),q(20),anames(20),qnames(20)
       common/fxxx/fxval
       common/ccr/cr,ctheta
!
       dimension s1(2),sd(2)
!
       aa= a(10)
       beta=a(6)
       delta=a(11)
       p11=q(9)
       tstar= t+1.0d0
!
       rhodi=  1.d0+p11*dtstin
       delrm= rm - rhodi
!
       call conver(dtstin,delrm,r,theta,unstab)
!
!      write(9,9900)r,theta
9900   format(2f10.5)
       tt= theta*theta
       if(r.gt. 0.0d0)
     & amu=  aa*r**(beta*delta)*theta*(1.0d0-tt)
!
!
!      write(6,200) t,rm,r,theta
200    format(' t,rm,r,theta  ',4e20.8)
       cr=r
       ctheta=theta
!
       return
       end
      subroutine conver(t,m,r,theta,unstab) 
!     this routine converts the reduced temperature and 
!     density to revised and extended parametric variables.

      implicit double precision (a-h,o-z)
      logical ltest 
      double precision m,k0,k1 
      parameter (itrmx2=100,smlno=1.0d-08) 
       common/cerr/cerrtk,cerrho,cerrp,cerrl,cerls
      common /coeffs/ a(20),q(10) 
      common /params/ alpha,beta,gamma,delta,weg,bb,p00,p02, 
     >                p04,p10,p12,p14,s00,s20,s01,s21 
       data inival/0/
!
      s0(th)=s00+s20*(th**2) 
      s1(th)=s01+s21*(th**2) 
!
       unstab=1.0d0
       if(inival.eq.0)then
                      call const
                      inival=1
!                     write(6,9898)s00,s20,s01,s21
9898   format(' ','s00,s20,s01,s21 ',4e20.8)
                      endif
      aa=a(1) 
      k0=a(2) 
      k1=a(3) 
      c =a(4) 
      if(k0.le.smlno)then 
         print *,  'invalid value for k0 --conver aborting' 
         return  
      endif 
      if(k1.lt.0.0d0) k1=-k1 
***   get initial estimates for r & theta *** 
      call rtheta(t,m,k0,r0,th0,unstab) 
       if(abs(unstab+1.0)<1.0E-8)go to 400
*** 
***   now get actual r & theta *** 
*** 
      do 200 k=1,itrmx2 
!      write(6,9875)t,aa,c,r0,beta,delta,th0 
9875   format(7e15.7)
         if(r0.le.0.0d0) go to 200
         p=t+aa*c*(r0**(beta*delta))*th0*(1.d0-th0**2) 
         qq=(m-aa*c*(r0**(1.d0-alpha))* 
     >      (k0*s0(th0)+k1*s1(th0)*r0**weg))/ 
     >                     (1.d0+k1*(r0**weg)/k0) 
         call rtheta(p,qq,k0,r1,th1,unstab) 
         if(abs(unstab+1.0)<1.0E-8) go to 400
         if(dabs(r1-r0)/(1.d0+r0).lt.1.d-12.and. 
     >      dabs(th1-th0)/(1.d0+dabs(th0)).lt.1.d-12)then 
       r=(r1+r0)/2.d0 
             theta=(th1+th0)/2.d0 
             return
         else 
            th0=th1 
      r0=r1 
         endif 
200   continue 
      print *,'can''t get routine to converge for r & theta' 
      print *,'when t=',t,'   m=',m,'   rinit=',r0,'   thetainit=',th0 
      write(6,150)cerrtk,cerrho,cerrp,cerrl
150   format(' ','orig call args: tk,rho,p,l: ',3f15.8,f3.0)
      write(6,151)
151   format(' ','(in unstable supersat region)')
400      continue
         r=dabs(r0)
         theta=th0
300   return  
      end 
      subroutine rtheta(t,m,k0,r,theta,unstab) 
       implicit double precision(a-h,o-z)
       common/cerr/cerrtk,cerrho,cerrp,cerrl,cerls
      double precision m,k0 
      parameter (eps=1.0d-12) 
      parameter (beta=0.325d0,bb=1.3757d0) 
****  this subroutine converts t and m into the parametric 
***   coordinates r & theta.  the particular model used here 
***   is the linear model of schofield; the reference used is 
***   m.r. moldover, j.res.nat.bur.std. 83, 329 (july-august 1978). 
***   note that this routine may be used iteratively to compute 
***   models which deviate slightly from the symmetric linear case. 
***   t is the reduced temperature 
***   m is the order parameter 
***   k0 is the non-universal constant appearing in the linear model 
***   r is the parametric variable describing the distance from 
*     the critical point 
***   theta is the parametric variable specifying the 
*     thermodynamic path 
***   "*" is a statement number in the calling program to which the 
*     rtheta routine should return to if an error occurrs. 
      dabst=dabs(t) 
      dabsm=dabs(m) 
      b=dsqrt(bb) 
      if(k0.lt.0.0d0)k0=-k0 
      if(bb.le.1.0d0)then 
         print *,'b is < or = 1.0; rtheta routine aborting' 
         return 
      endif 
***   if the system is on the critical isochore 
      if(dabsm.lt.eps)then 
         if(dabst.lt.eps)then 
            theta=0.0d0 
            r=0.0d0 
            return 
         else if (t.le.-eps) then 
            theta=dsign(1.0d0,m) 
            r=t/(1.d0-bb) 
            return 
         else if (t.ge.eps) then 
            theta=0.0d0 
            r=t 
            return 
         endif 
!     if, however, we have a non-zero order parameter 
      else 
!        critical isotherm 
         if(dabst.lt.eps)then 
            theta=dsign(1.d0,m)/b 
            r=(m/(k0*theta))**(1.0d0/beta) 
!           write(6,9897)m,r,theta
9897        format(' m,r,theta: ',4e20.8)
            return 
         else if (t.le.-eps)then 
!           check for 2-phase region 
             zold=(-t*(k0/dabsm)**(1.d0/beta))
             zold=  dsqrt(zold+1.0d0)*1.00234/1.07
             if(zold.gt.b) then
!					       if((-t*(k0/dabsm)**(1.d0/beta)).gt.(bb-1.0d0))then 
               theta=dsign(1.0d0,m) 
               r=t/(1.d0-bb) 
               unstab=-1.0
               return 
            endif 
!           get initial estimate by linear interpolation 
!           z=b*theta 
            zinit=1.d0-(1.d0-b)*t*(k0/dabsm)**(1.d0/beta)/(1.d0-bb) 
         else if (t.ge.eps)then 
!           get initial estimate by rational interpolation 
!	      z=b*theta 
         zinit=(1.d0+t*(k0/(b*dabsm))**(1.d0/beta))**(-beta) 
      endif 
      c=-m*b/(k0*dabst**beta) 
      z=dsign(zinit,m) 
!     use newton-raphson method to get exact z 
!      write(6,9898)z,c,t,m
9898   format(3x,'z,c,t,m: ',6e20.8)
      do 10 i=1,20 
         	dz=(1.d0-z**2)*(z+c*dabs(1.d0-z**2)**beta)/ 
     &			((2.d0*beta-1.d0)*z**2+1.d0) 
      		z=z-dz 
      		if(dabs(dz)/(1.d0+dabs(z)).lt.eps)then 
      				r=t/(1.0d0-z**2) 
      				theta=z/b 
    					if(dabs(theta).gt.(1.0d0-eps))theta=dsign(1.d0,theta) 
     			return 
     			endif 
10    		continue 
      endif 
      print *,' rtheta failed to converge' 
      print *,' t=',t,'   m=',m,' dz=',z 
      write(6,150)cerrtk,cerrho,cerrp,cerrl
150   format(' ','orig call args: tk,rho,p,l: ',3f15.8,f3.0)
      write(6,151)
151   format(' ','(in unstable supersat region)')
       unstab=-1.0d0
       r= dabs(r)
      return 
      end 
      
      subroutine const 
!     this routine calculates the various constants needed 
!     for the new revised conver. 
       implicit double precision (a-h,o-z)
      common /params/ alpha,beta,gamma,delta,weg,bb,p00 
     &		,p02,p04,p10,p12,p14,s00,s20,s01,s21 
!
!
       alpha= 0.1085d0
       beta= 0.325d0
       gamma=1.2415d0
       delta= 4.82d0
       weg= .5d0
       bb= 1.3757d0
!
      p00=(beta*(delta-3.d0)-bb*alpha*gamma)/ 
     &		(2.d0*(bb**2)*(2.d0-alpha)*(1.d0-alpha)*alpha) 
      p20=-(beta*(delta-3.d0)-bb*alpha*(2.d0*beta*delta-1.d0))/ 
     &		(2.d0*bb*(1.d0-alpha)*alpha) 
      p40=(2.d0*beta*delta-3.d0)/(2.d0*alpha) 
      p01=(beta*(delta-3.d0)-3.d0*weg-bb*(alpha-weg)*(gamma-weg))/ 
     &		(2.d0*(bb**2)*(2.d0-alpha+weg)*(1.d0-alpha+weg)*(alpha-weg)) 
      p21=-(beta*(delta-3.d0)-3.d0*weg-bb*(alpha-weg)* 
     &		(2.d0*beta*delta-1.d0))/ 
     &		(2.d0*bb*(1.d0-alpha+weg)*(alpha-weg)) 
      p41=(2.d0*beta*delta-3.d0)/(2.d0*(alpha-weg)) 
      s00=(2.d0-alpha)*p00 
      s20=-beta*(delta-3.d0)/(2.d0*bb*alpha) 
      s01=(2.d0-alpha+weg)*p01 
      s21=-(beta*delta-3.d0*(beta+weg))/ 
     &		(2.d0*bb*(alpha-weg)) 
      return 
      end 
      
      block data bcoeff 
      implicit double precision (a-h,o-z)
      common /coeffs/ a(20),q(10) 
!   	this program initializes the parameters needed in the 
!   	revised and extended parametric model. 
!
!      note that the constants here are for the new,revised
!           conver, and the order is different from that in
!           blockdata clsblk
!     (actually, only a(1) thru a(4) are referred to in new conver)
!
*     a(1)=a 
*     a(2)=k0 
*     a(3)=k1 
*     a(4)=c 
*     a(5)=p1 
*     a(6)=p2 
*     a(7)=p3 
*     a(8)=p11 
*     a(9)=delta tc 
*     a(10)=delta rhoc 
*     a(11)=delta pc 
*     a(12)=mu0 
*     a(13)=mu1 
*     a(14)=mu2 
*     a(15)=mu3 
*     a(16) - a(20) are unused 
*     q(1)=tc 
*     q(2)=rhoc 
*     q(3)=pc 
*     q(4)=dp/dt evaluated at the critical point in actual units 
*     q(5)=delta pc calculated from shift in tc and dp/dt at !.p. 
*     q(6)=slope of the diameter in the classical region 
*     q(7) - q(10) are unused 
       data a/23.6666d0,1.4403d0,.2942d0,-.017762d0,6.8445d0,
     *        -25.4915d0,5.238d0,0.4918d0,0.0d0,0.0d0,
     *         0.0d0,-11.2326d0,-22.6547d0,-17.8876d0,-4.9332d0,
     *         5*0.0d0/
!
       data q/647.067d0,322.778d0,22.046d0,.267d0,.0d0,
     *      -1.6d0,4*0.d0/
      end 
        subroutine swvir(tkd,bz,cz)
!
!       subroutine to calculate 
!             bz:  second virial coefficient
!             cz:  third virial  coefficient
!
!
       implicit double precision(a-h,o-z)
       double precision c1(10,10),c2(10,12),c3(10,10)
       double precision c4(10,10)
       common/cvcc/c1,c2,c3,c4
       common/efgval/eval,dr0,dt0,delta,aval,abval,bval,cval,dval
       common/norms/tc,rhoc,rval
       nr1 = 7
       nt1 = 7
       nr2 = 7
       nt2 = 12
       nr3 = 5
       nt3 = 5
       nr4 = 5
       nt4 = 10
       nrmax = 8
       ntmax = 12
!
!
        tbar= -tc/tkd
        dtbar= tbar + 1.0d0
!
!       summation where i=1 for b, i=2 for c;
!       coefficients c1 correspond to k=1
!                    c2 ~ k=2
!                    c3 ~ k=3
!
!
       bz= 0.0
       cz= 0.0
!
       do 200 j=1,nt1
       ba=  eval * c1(1,j) * tbar**(j-1)
       bz= bz+ba
       tbj= tbar**(j-1)
200    continue
!
       do 220 j=1,nt2
       bb= c2(1,j) * tbar**(j-1)
       tbj= tbar**(j-1)
       bz= bz +bb
       cz= cz + c2(2,j) * tbar**(j-1)
220    continue
!
       eaval= (-aval*dtbar-cval*dtbar**2)
       dexlim= -35.
       darg= eaval
       if(darg.ge.dexlim) expval=dexp(darg)
       if(darg.lt.dexlim) expval=0.0d0
       eaval= expval
       ecval= dexp(-dval)
!
!      do 230 j=1,nt3
!      bc= c3(1,j)*ecval*eaval*tbar**(j-1)*tbar**2
!      cz= cz + c3(2,j)*ecval*eaval*tbar**(j-1)*tbar**2
!      bz= bz+bc
!      tbj= tbar**(j-1)
c230    continue
!
       h= dexp(-4.0d0*(tbar+3.0d0))
!
       do 240 j=1,nt4
       bd= (eval*h*c4(1,j))*tbar**(j-1)
       cd= (eval*h*c4(2,j))*tbar**(j-1)
       bz= bz+bd
       cz=cz+cd
240    continue
!
       bz= bz / rhoc/1000.
       cz= cz /(rhoc*1000.)**2
       return
       end
!
!   ********************
!
!    saul,wagner (1985) h2o saturation functions
!
!     ref:  wagner, w, sengers,j
!           draft release on saturation properties
!           technical report bn 1051
!
!
!                       [fpwag,fdpwag,rfwag,rgwag,fhfwag,fhgwag,fawag]
!
!   ********************************************************
!
       double precision function fpswag(tk)
       implicit double precision (a-h,o-y)
       tc= 647.14d0
       if(tk.ge.tc) then
                    fpswag=22.064d0
                    write(6,999)tk
                    return
999                 format(' fpswag: tk ge tc ',f20.10)
                    endif
       tau=  1.0d0 - tk/tc
       w= -7.858230d0 * tau
       w= w + 1.839910d0 * tau**1.5
       w= w + (-11.78110d0 *tau**3)
       w= w + (22.6705d0 * tau**3.5)
       w= w + (-15.9393d0 * tau**4)
       w= w + (1.77516d0 * tau**7.5)
       dexlim= -35.
       darg= tc/tk*w
       if(darg.ge.dexlim) w=dexp(darg)
       if(darg.lt.dexlim) w=0.0d0
       fpswag= 22.064d0 * w
       return
       end
       
       
       double precision function fdpsw(tk)
       implicit double precision (a-h,o-y)
       ps= fpswag(tk)
       tc= 647.14
       if(tk.ge.tc) then
           fdpsw= 22.064d0
           write(6,999) tk
           return
999        format(' fdpsw: tk gt tc ',f20.10)
       endif
       tau=  1.0 - tk/tc
       w= -7.858230 
       w= w + 1.839910 * tau**0.5*1.5
       w= w + (-11.78110 *tau**2*3.0)
       w= w + (22.6705 * tau**2.5*3.5)
       w= w + (-15.9393 * tau**3*4.0)
       w= w + (1.77516 * tau**6.5*7.5)
       w= dlog(ps/22.064d0) + w
       fdpsw= w *(-ps/tk)
       return
       end
       
       double precision function frfwag(tk)
       implicit double precision (a-h,o-y)
       tkc= 647.14
       rc= .322
       if(tk.ge.tkc) then
       		frfwag=-10
       		print *, "Too high temp"
			 else
       		tau= 1.0d0 - tk/tkc
       		w= 1.99206 *tau**(1.0/3.0)
       		w= w + 1.10123* tau**(2.0/3.0)
       		w= w + (-0.512506 * tau**(5./3.))
       		w= w + (-1.75263 * tau**(16./3.))
       		w= w + (-45.4485 * tau**(43./3.))
       		w= w + (-675615. * tau**(110./3.))
       		frfwag= (w+1.0)*rc
       end if
       return
       end
       
       
       double precision function frgwag(tk)
       implicit double precision (a-h,o-y)
       tc= 647.14
       if(tk.ge.tc) then
       		frgwag=-10
       		print *, "Too high temp"
       else
       		tau= 1.0 - tk/tc
       		w=  -2.029570 * tau**(2.0/6.0)
       		w=  w + (-2.68781*tau**(4./6.))
       		w= w +  (-5.381070*tau**(8./6.))
       		w= w +  (-17.31510*tau**(18./6.))
       		w= w +  (-44.6384 *tau**(37./6.))
       		w= w +  (-64.3486 *tau**(71./6.))
       		dexlim= -35.
       		darg= w
       		if(darg.ge.dexlim) w=dexp(darg)
       		if(darg.lt.dexlim) w=0.0d0
       		frgwag= .322 * w
       end if
       return
       end
       double precision function fbzlnr(tk)
       implicit double precision (a-h,o-z)
!
!      ref:  lefvre,e, nightingale,m, rose, j
!            j.mech.eng.sci. vol 17,no 5 (1975)
!
       b= 1.0d0/(1.0d0+tk/10000.d0)*.0015d0
       darg= -1500.d0/tk
       dexlim= -35.
       if(darg.ge.dexlim) expval=dexp(darg)
       if(darg.lt.dexlim) expval=0.0d0
       c= (1.d0-expval)**(5.d0/2.d0)
       c= c*dsqrt(tk/1500.d0)*dexp(1500.d0/tk)*(-.000942d0)
       d= 1500.d0/tk * (-.0004882d0)
       fbzlnr= b + c +d
       return
       end
!
!  ********************************************
!
!      spinodal rho (liq and vap) correlations
!         of wsteam tk and rho values where dpdr=0.0 
!         for use in iteration to rho from pressure 
!         value arguments in the supersaturated
!         region, 
!
!         in addition these functions are used to 
!         check input pressure values (l=3 or l=4)
!         to give an error message when the rho
!         for a given pressure lies in the unstable
!         supersaturated region 
!
!  ***********************************************
       double precision function wspinf(tk)
       implicit double precision (a-h,o-z)
       tkc= 647.067d0
       rc=  .322778d0
       wspinf=0.0d0
      if(tk.ge.tkc)return
!
       tau= 1.0d0 - tk/tkc
       w= -.5948463865d0 *tau**(1.0d0/3.0d0)
       w= w + .5080947833d2* tau**(2.0d0/3.0d0)
       w= w + (-.3797642336d3 * tau**(3.d0/3.0d0))
       w= w + (.1249347735d4 * tau**(4.d0/3.0d0))
       w= w + (-.1860121245d4 * tau**(5.d0/3.0d0))
       w= w + (.1031487945d4 * tau**(6.d0/3.0d0))
       w= w + (-.1574103153d4 * tau**(17.d0/3.0d0))
       w= w + (.9435352115d4  * tau**(21.d0/3.0d0))
       w= w + (-.8168784707d4 * tau**(22.d0/3.0d0))
       wspinf= (w+1.0d0)*rc
       return
       end
       double precision function wsping(tk)
       implicit double precision (a-h,o-z)
       tc= 647.067d0
       rc=  .322778d0
       wsping=0.0d0
       if(tk.ge.tc) return
       tau= 1.0 - tk/tc
       w=  -.3346401461d1 * tau**(1.0d0/3.0d0)
       w=  w + (.1268114311d2*tau**(2.d0/3.0d0))
       w= w +  (-.2316813580d2*tau**(3.d0/3.0d0))
       w= w +  (-.6567012350d2*tau**(4.d0/3.0d0))
       w= w +  (0.2166943752d3 *tau**(5.d0/3.0d0))
       w= w +  (-.1663688422d3 *tau**(6.d0/3.0d0))
       w= w +  (.1409844499d3 *tau**(16.d0/3.0d0))
       w= w +  (-.1501900170d4 *tau**(26.d0/3.0d0))
       w= w +  (.1722928004d4   *tau**(29.d0/3.0d0))
       dexlim= -35.
       darg= w
       if(darg.ge.dexlim) w=dexp(darg)
       if(darg.lt.dexlim) w=0.0d0
       wsping= rc * w
       return
       end
       
       
       block data chiblk
       implicit double precision(a-h,o-z)
!      deleted unused array w(4,6)
       double precision c1(10,10),c2(10,12),c3(10,10),c4(10,10)
       double precision c0(8)
       common/czero/c0
       common/ccread/cread
       common/cinit/c000,c001
       common /cvcc/c1,c2,c3,c4
       data cread/0.0d0/
       data c0/.707501277532d1,-.834240573672d1,-.36460138d0
     *         ,-.36897043d-1,0.30338153d-2,
     *         0.39010933d-3,0.11359287d0,.24131785d1/
       data c000/-.3474470158d-4/
       data c001/-.3040231678d-4/
        data c1/
     &  3*0.d0,0.3384249125d+0,-.7153393406d-1,0.d0,.549368081d-3
     &  ,3*0d0,
     &  2*0.d0,.4933218501d-1,2*0.d0,-.2328491212d-1,.2402095181d-2
     &  ,3*0.d0,
     &  .7529422956d+0,0.d+0,-.2280260070d+1,7*0.d+0,
     &  0.d+0,.1142004144d+1,-.2619059624d+1,0.d+0,.4395237702d+0,
     &  -.3161046646d-1,.6814467692d-3,3*0.d+0,
     &  -.3924227294d+0,0.d+0,-.2738770648d+0,7*0.d+0,
     &  3*0.d+0,-.1943443857d-01,.3048860434d-2,5*0.d+0,
     &  2*0.d+0,.3946510403d-2,7*0.d+0
     &  ,30*0.d0/
       data c2/
     &  .22436101314d0,.1193250201d0,2*0.d0,.6582959348d-1,
     &   5*0.d0,
     &  4*0.d0,.1651430628d0,5*0.d0,
     &  -.2178969357d1,.2674090542d0,.8647490995d0,7*0.d0,
     &  -.1530432257d0,0.d0,.2059881454d1,2*0.d0,
     &  -.4888628703d0,.1375328753d0,3*0.d0,
     &  4*0.d0,-.9015180666d0,-.1444258609d0,.1558046279d0,3*0.d0
     &  ,3*0.d0,-.2740652563d1,0.d0,.4983771706d0,4*0.d0,
     &  3*0.d0,-.3261978564d1,.1609338784d1,5*0.d0,
     &  .3484674963d-1,-.1537646434d1,2*0.d0,.2316225257d0,5*0.d0,
     &  0.d0,-.1419249232d1,.7969984635d0,7*0.d0,
     &  4*0.d0,.7510544627d-2,5*0.d0,
     &  10*0.d0,
     &	.5364384732d-3,9*0.d0/
        data c3/
     & .6109381296d0,0.d0,-.1906644459d-1,0.d0,.7976092188d-2,
     & 5*0.d0,
     & .1934466766d1,9*0.d0,
     & .1921820547d1,0.d0,-.4410105919d-1,7*0.d0,
     & .6130354419d0,-.2855258689d0,2*0.d0,.2526137080d-1,5*0.d0,
     & 0.d0,-.2374074642d+0,0.d+0,.3855866402d-1,.8041672150d-2,
     & 55*0.d0/
       data c4/
     & 2*0.d0,-.1635439033d2,7*0.d0,
     & -.5025818675d2,9*0.d0,
     & 10*0.d0,
     & 0.d0,.1649003040d0,8*0.d0,
     & -.8499893502d0,9*0.d0,
     & 30*0.d0,
     & 0.8314382544d-2,.8781327858d-3,8*0.d0,
     & 0.d0,.1537391213d-2,-.9016873786d-3,0.d0,.3326628664d-3,
     & 5*0.d0/
       end

