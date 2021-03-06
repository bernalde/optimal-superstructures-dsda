* Code to replicate Figure 6 in Optimal design of superstructures for placing 
* units and streams with multiple and ordered available locations. 
* Part I: A new mathematical framework

*------------------------SETS---------------------------------------------------
set N "Superstructure units" /1*5/;
set I "Chemical components" /A,B/;
alias(N,N1);

*-----------------------PARAMETERS AND SCALARS----------------------------------
***Kinetic parameters
scalar m1 "Partial order of reaction 1 []" /1/;
scalar m2 "Partial order of reaction 2 []" /1/;
scalar k "kientic constant [L/(mol*s)]"/2/;
***Feed parameters
scalar QF0 "Inlet volumetric flow [L/s]" /1/;
parameter C0(i) "Initial concentration of reagents [mol/L]"
/
A 0.99
B 0.01
/
;
parameter F0(i) "Inlet molar flow [mol/s]";
F0(i)=C0(i)*QF0;

*-----------------------BINARY TERMS--------------------------------------------
***Independent binary terms
parameter yf(n) "Location of the last reactor (feed reactor) in the sueprstructure";
parameter yr(n) "Location of the recycle flow";
***Dependent binary terms
parameter yp(n) "Unit operation. 1: Reaction. 0: Simple input-output";
*-----------------------REAL VARIABLES------------------------------------------
***Network variables
positive variable Q(n) "Outlet flow rate of the superstricture unit [L/s]";
positive variable F(i,n) "Molar flow [mol/s]";
variable  rate(i,n) "Reaction rate [mol/(L*s)]";
positive variable V(n) "Reactor volume [L]";
***Splitter variables
positive variable QR "Recycle flow rate  [L/s]";
positive variable QP "Product flow rate  [L/s]";
positive variable R(i) "Recycle molar flow [mol/s]";
positive variable P(i) "Product molar flow [mol/s]";

*-----------------------SUPERSTRUCTURE CONSTRAINTS------------------------------

***Kinetic constraints
equation net_rate(i,n) "Reactor Network: Reaction rates";
net_rate(i,n)..(rate('A',n)*((Q(n))**m1)*((Q(n))**m2)+k*((F('A',n))**m1)*((F('B',n))**m2))$(ord(i) eq 1)+(rate('B',n)+rate('A',n))$(ord(i) eq 2)=e=0;

***Network constraints
equation net_comp(i,n) "Reactor Network: Component molar balance";
equation net_cont(n) "Reactor Network: continuity equation";
net_comp(i,n)$(ord(n) ne card(n))..F(i,n+1)+yr(n)*R(i)-F(i,n)+yp(n)*rate(i,n)*V(n)=e=0;
net_cont(n)$(ord(n) ne card(n))..Q(n+1)+yr(n)*QR-Q(n)=e=0;


***Feed unit constraints
equation feed_comp(i,n) "Feed unit: Component molar balance";
equation feed_cont(n) "Feed unit: continuity equation";
feed_comp(i,n)$(ord(n) eq card(n))..F0(i)+yr(n)*R(i)-F(i,n)+yp(n)*rate(i,n)*V(n)=e=0;
feed_cont(n)$(ord(n) eq card(n))..QF0+yr(n)*QR-Q(n)=e=0;

***Splitter constraints
equation spl_comp(i) "Splitter: Component molar balance";
equation spl_cont "Splitter: continuity equation";
equation spl_ad(i) "Splitter: additional splitter constraints";
spl_comp(i)..F(i,'1')-P(i)-R(i)=e=0;
spl_cont..Q('1')-QP-QR=e=0;
spl_ad(i)..P(i)*Q('1')-F(i,'1')*QP=e=0;

*---------------------VOLUMEN CONSTRAINT----------------------------------------
equation eqvol(n);
eqvol(n)$(ord(n) ne 1)..V(n)=e=V(n-1);

*---------------------PRODUCT QUALITY CONSTRAINT--------------------------------
equation qual;
qual..QP*0.95=e=P('B');
*----------------------OBJECTIVE FUNCTION---------------------------------------

variables zobj;

equation Fobj;
Fobj..zobj=e=sum(n,V(n)*yp(n));
*----------------------MODELS---------------------------------------------------
***CSTR model
model SUP_CSTR /all/;

***External variables model (for feasible initialization only)
integer variables variablex1,variablex2;
equations feas;
feas..variablex2-variablex1+1=l=0;
variable z_ext;
equation ob_ext;
ob_ext..z_ext=e=1;
model EXT_INIT /all-SUP_CSTR/;
*----------------------VARIABLE INITIALIZATION AND BOUNDS-----------------------
***CSTR model: bounds
Q.up(n)=10;
QR.up=10;
QP.up=10;
R.up(i)=10;
P.up(i)=10;
F.up(i,n)=10;
rate.lo(i,n)=-10;
rate.up(i,n)=10;
V.up(n)=10;
***External variables model: bounds (for feasible initialization only)
variablex1.lo=1;
variablex1.up=card(N);
variablex2.lo=0;
***CSTR model: Initialization
Q.l(n)=0;
F.l(i,n)=0;
rate.l(i,n)=0;
V.l(n)=0;
QR.l=0;
QP.l=0;
R.l(i)=0;
P.l(i)=0;
***External variables model: initialization (for feasible initialization only)
variablex1.l=0;
variablex2.l=0;
*---------------------SUBSOLVERS------------------------------------------------
***NLP solver for A3
option nlp=msnlp;
$onecho > baron.opt
EpsR 0.0001
$offecho
SUP_CSTR.OptFile = 1;

***MIP solver for external variables model (for feasible initialization only)
option mip=cplex;
*--------------------FEASIBLE INITIALIZATION------------------------------------
***Feasible initialization of external variables
solve EXT_INIT using mip minimizing z_ext;
***Feasible initialization of CSTR model variables
yf(n)=0;
yf(n)$(ord(n) eq floor(variablex1.l))=1-mod(variablex1.l,1);
yf(n)$(ord(n) eq 1+floor(variablex1.l))=mod(variablex1.l,1);
yr(n)=0;
yr(n)$(ord(n) eq floor(variablex1.l-variablex2.l))=1-mod(variablex1.l-variablex2.l,1);
yr(n)$(ord(n) eq 1+floor(variablex1.l-variablex2.l))=mod(variablex1.l-variablex2.l,1);
yp(n)=1-(sum(n1$(ord(n1) <= ord(n)),yf(n1))-yf(n));
solve SUP_CSTR using nlp minimizing zobj;
execute_unload "in";

*--------------------SOLUTION WITH FEASIBLE INITIALIZATION----------------------
***ITER SETS
set iter "iteraciones" /i1*i100/;
set iner1 "iteraciones internas 1" /i1*i10/;
set iner2 "iteraciones internas 2" /i1*i100/;
set extvar "variables externas" /x1*x2/;
set extineq "restricciones desigualdad externas" /g1*g4/;

***INFINITY-NEIGHBORHOOD
set neigh /d1*d4/;
table directions(neigh,extvar)
    x1       x2
d1  -1        0
d2  1        0
d3  0        -1
d4  0        1
;

**DEFINED PARAMETERS
scalar hstep "paso" /1/;
parameter hiner1(iner1) "paso iteraiones internas 1";
loop(iner1,
if(ord(iner1) eq 1,
hiner1(iner1)=0;
else
hiner1(iner1)=hiner1(iner1-1)+(hstep/((card(iner1))-1));
);
);
parameter hiner2(iner2) "paso iteraciones internas 2";
loop(iner2,
if(ord(iner2) eq 1,
hiner2(iner2)=hstep;
else
hiner2(iner2)=hiner2(iner2-1)+hstep;
);
);
parameter extvarhstep(extvar,neigh) "parametro usado en la evaluacion del gradiente";
extvarhstep(extvar,neigh)=hstep*directions(neigh,extvar);

**PARAMETROS DEFINIDOS EN PROCESO ITERATIVO
scalar stop1 "criterio de parada principal";
scalar stop2 "criterio de parada para evaluacion del gradiente por infactibilidad";
scalar stop3 "criterio de parada cuando ya se tiene el gradiente";
scalar stop4 "criterio de parada cuando ya se tiene el gradiente por infactibilidad";
scalar count "contador";
parameter xvalue(iter,extvar) "valor de x en cada iteracion";
parameter xvalueselect(iter,extvar) "valor de x seleccionado en la siguiente iteracion";
parameter gval(extineq) "valor de g en x";
parameter dvs(iter,neigh) "valores d";
parameter dmin(iter) "valor d minimo en la iteracion actual";
parameter fvalue(iter) "valor de fobj en cada iteracion";
parameter fplushvalue(iter,neigh) "fobj para calcular valores d";
parameter fvalueiner(iter,iner2) "valor de fobj en iteraciones internas";
parameter selectd(iter,neigh) "direccion seleccioanda para optimizar";
parameter mstatS2(iter)"model status de la seccion S2";
parameter mstatS4(iter,neigh,iner1) "model status de la seccion S4";
parameter mstatS5(iter,neigh,iner1) "model status de la seccion S5" ;
parameter mstatS7(iter,iner2,iner1) "model status de la seccion S7";
parameter mstatS8(iter,iner2,iner1) "model status de la seccion S8";
parameter objval "Objective function value";
parameter xvalinit(extvar) "valor de inicalizacion para cada etapa de boil up";
parameter CPUtime;
scalar CPUtimeactual;

CPUtimeactual=0;
xvalue(iter,extvar)=0;
xvalueselect(iter,extvar)=0;
gval(extineq)=0;
dvs(iter,neigh)=0;
dmin(iter)=0;
fvalue(iter)=0;
fplushvalue(iter,neigh)=0;
fvalueiner(iter,iner2)=0;
selectd(iter,neigh)=0;
mstatS2(iter)=0;
mstatS4(iter,neigh,iner1)=0;
mstatS5(iter,neigh,iner1)=0;
mstatS7(iter,iner2,iner1)=0;
mstatS8(iter,iner2,iner1)=0;



xvalinit("x1")=variablex1.l;
xvalinit("x2")=variablex2.l;
xvalue('i1',extvar)=xvalinit(extvar);


stop1=0;
loop(iter$((stop1 eq 0) and (ord(iter) ne card(iter))),

***S1_evaluacion de inicializacion factible
         if(ord(iter) eq 1,
         gval("g1")=1-xvalue(iter,"x1");
         gval("g2")=xvalue(iter,"x1")-card(n);
         gval("g3")=1-(xvalue(iter,"x1")-xvalue(iter,"x2"));
         gval("g4")=(xvalue(iter,"x1")-xvalue(iter,"x2"))-xvalue(iter,"x1");

                 if( sum(extineq$(gval(extineq) gt 0 ),gval(extineq)) gt 0,
                 stop1=1;
                 );

         );
***S2_Calculo de fvalue(iter)
         if(ord(iter) eq 1,
         execute_loadpoint "in";
         else
         execute_loadpoint "DLR";
         );

yf(n)=0;
yf(n)$(ord(n) eq floor(xvalue(iter,"x1")))=1-mod(xvalue(iter,"x1"),1);
yf(n)$(ord(n) eq 1+floor(xvalue(iter,"x1")))=mod(xvalue(iter,"x1"),1);
yr(n)=0;
yr(n)$(ord(n) eq floor(xvalue(iter,"x1")-xvalue(iter,"x2")))=1-mod(xvalue(iter,"x1")-xvalue(iter,"x2"),1);
yr(n)$(ord(n) eq 1+floor(xvalue(iter,"x1")-xvalue(iter,"x2")))=mod(xvalue(iter,"x1")-xvalue(iter,"x2"),1);
yp(n)=1-(sum(n1$(ord(n1) <= ord(n)),yf(n1))-yf(n));

solve SUP_CSTR using nlp minimizing zobj;
mstatS2(iter)=SUP_CSTR.modelstat;
         if((mstatS2(iter) eq 3 or mstatS2(iter) eq 4 or mstatS2(iter) eq 5 or mstatS2(iter) eq 6 or mstatS2(iter) eq 11 or mstatS2(iter) eq 12 or mstatS2(iter) eq 13 or mstatS2(iter) eq 14 or mstatS2(iter) eq 18 or mstatS2(iter) eq 19) and ord(iter) eq 1,
         stop1=1;
         elseif (mstatS2(iter) eq 3 or mstatS2(iter) eq 4 or mstatS2(iter) eq 5 or mstatS2(iter) eq 6 or mstatS2(iter) eq 11 or mstatS2(iter) eq 12 or mstatS2(iter) eq 13 or mstatS2(iter) eq 14 or mstatS2(iter) eq 18 or mstatS2(iter) eq 19),
         fvalue(iter)=zobj.l;
         else
         execute_unload "DLR";
         fvalue(iter)=zobj.l;
         );


         loop(neigh,
***S3_Calculo de factibilidad de variables externas para los diferentes valors d

                 gval("g1")=1-(xvalue(iter,"x1")+extvarhstep("x1",neigh));
                 gval("g2")=(xvalue(iter,"x1")+extvarhstep("x1",neigh))-card(n);
                 gval("g3")=1-((xvalue(iter,"x1")+extvarhstep("x1",neigh))-(xvalue(iter,"x2")+extvarhstep("x2",neigh)));
                 gval("g4")=((xvalue(iter,"x1")+extvarhstep("x1",neigh))-(xvalue(iter,"x2")+extvarhstep("x2",neigh)))-(xvalue(iter,"x1")+extvarhstep("x1",neigh));

                          if( sum(extineq$(gval(extineq) gt 0 ),gval(extineq)) gt 0 ,
                          dvs(iter,neigh)=1;
                          else
                          dvs(iter,neigh)=-1;
                          );

                          if(dvs(iter,neigh) le 0,
                          stop2=0;

                                 loop(iner1$(stop2 eq 0),

                                         if(ord(iner1) eq 1,
***S4_Calculo de valores d cuando no hay problemas de factibilidad
                                         execute_loadpoint "DLR";

                                         yf(n)=0;
                                         yf(n)$(ord(n) eq floor(xvalue(iter,"x1")+extvarhstep("x1",neigh)))=1-mod(xvalue(iter,"x1")+extvarhstep("x1",neigh),1);
                                         yf(n)$(ord(n) eq 1+floor(xvalue(iter,"x1")+extvarhstep("x1",neigh)))=mod(xvalue(iter,"x1")+extvarhstep("x1",neigh),1);
                                         yr(n)=0;
                                         yr(n)$(ord(n) eq floor((xvalue(iter,"x1")+extvarhstep("x1",neigh))-(xvalue(iter,"x2")+extvarhstep("x2",neigh))))=1-mod((xvalue(iter,"x1")+extvarhstep("x1",neigh))-(xvalue(iter,"x2")+extvarhstep("x2",neigh)),1);
                                         yr(n)$(ord(n) eq 1+floor((xvalue(iter,"x1")+extvarhstep("x1",neigh))-(xvalue(iter,"x2")+extvarhstep("x2",neigh))))=mod((xvalue(iter,"x1")+extvarhstep("x1",neigh))-(xvalue(iter,"x2")+extvarhstep("x2",neigh)),1);
                                         yp(n)=1-(sum(n1$(ord(n1) <= ord(n)),yf(n1))-yf(n));

                                         solve SUP_CSTR using nlp minimizing zobj;
                                         mstatS4(iter,neigh,iner1)=SUP_CSTR.modelstat;

                                                 if(mstatS4(iter,neigh,iner1) eq 3 or mstatS4(iter,neigh,iner1) eq 4 or mstatS4(iter,neigh,iner1) eq 5 or mstatS4(iter,neigh,iner1) eq 6 or mstatS4(iter,neigh,iner1) eq 11 or mstatS4(iter,neigh,iner1) eq 12 or mstatS4(iter,neigh,iner1) eq 13 or mstatS4(iter,neigh,iner1) eq 14 or mstatS4(iter,neigh,iner1) eq 18 or mstatS4(iter,neigh,iner1) eq 19,
                                                 stop2=0;
                                                 else
                                                 stop2=1;
                                                 fplushvalue(iter,neigh)=zobj.l;
                                                 dvs(iter,neigh)=(fplushvalue(iter,neigh)-fvalue(iter))/(hstep);
                                                 );

                                         else
***S5_Calculo de valores d cuando hay problemas de factibilidad

                                                 if(ord(iner1) eq 2,
                                                 execute_loadpoint "DLR";
                                                 else
                                                 execute_loadpoint "DLR1";
                                                 );
                                         yf(n)=0;
                                         yf(n)$(ord(n) eq floor(xvalue(iter,"x1")+((extvarhstep("x1",neigh))/(hstep))*(hiner1(iner1))))=1-mod(xvalue(iter,"x1")+((extvarhstep("x1",neigh))/(hstep))*(hiner1(iner1)),1);
                                         yf(n)$(ord(n) eq 1+floor(xvalue(iter,"x1")+((extvarhstep("x1",neigh))/(hstep))*(hiner1(iner1))))=mod(xvalue(iter,"x1")+((extvarhstep("x1",neigh))/(hstep))*(hiner1(iner1)),1);
                                         yr(n)=0;
                                         yr(n)$(ord(n) eq floor(xvalue(iter,"x1")+((extvarhstep("x1",neigh))/(hstep))*(hiner1(iner1))-(xvalue(iter,"x2")+((extvarhstep("x2",neigh))/(hstep))*(hiner1(iner1)))))=1-mod(xvalue(iter,"x1")+((extvarhstep("x1",neigh))/(hstep))*(hiner1(iner1))-(xvalue(iter,"x2")+((extvarhstep("x2",neigh))/(hstep))*(hiner1(iner1))),1);
                                         yr(n)$(ord(n) eq 1+floor(xvalue(iter,"x1")+((extvarhstep("x1",neigh))/(hstep))*(hiner1(iner1))-(xvalue(iter,"x2")+((extvarhstep("x2",neigh))/(hstep))*(hiner1(iner1)))))=mod(xvalue(iter,"x1")+((extvarhstep("x1",neigh))/(hstep))*(hiner1(iner1))-(xvalue(iter,"x2")+((extvarhstep("x2",neigh))/(hstep))*(hiner1(iner1))),1);
                                         yp(n)=1-(sum(n1$(ord(n1) <= ord(n)),yf(n1))-yf(n));

                                         solve SUP_CSTR using nlp minimizing zobj;
                                         mstatS5(iter,neigh,iner1)=SUP_CSTR.modelstat;
                                         execute_unload "DLR1";

                                                 if(ord(iner1) eq card(iner1),
                                                         if(mstatS5(iter,neigh,iner1) eq 3 or mstatS5(iter,neigh,iner1) eq 4 or mstatS5(iter,neigh,iner1) eq 5 or mstatS5(iter,neigh,iner1) eq 6 or mstatS5(iter,neigh,iner1) eq 11 or mstatS5(iter,neigh,iner1) eq 12 or mstatS5(iter,neigh,iner1) eq 13 or mstatS5(iter,neigh,iner1) eq 14 or mstatS5(iter,neigh,iner1) eq 18 or mstatS5(iter,neigh,iner1) eq 19 ,
                                                         dvs(iter,neigh)=1;
                                                         else
                                                         fplushvalue(iter,neigh)=zobj.l;
                                                         dvs(iter,neigh)=(fplushvalue(iter,neigh)-fvalue(iter))/(hstep);
                                                         );
                                                 );

                                         );

                                 );

                         );

         );

***S6_verificacion del criterio principal de parada y seleccion de la direccion de minimizacion
*Avoid numerical errors
dvs(iter,neigh)=round(dvs(iter,neigh),8);
dmin(iter)=smin((neigh),dvs(iter,neigh));
         if(dmin(iter)  ge 0 ,
         stop1=1;
         );

count=0;
         loop(neigh,
                         if(dvs(iter,neigh) eq dmin(iter) and count eq 0,
                         selectd(iter,neigh)=1;
                         count=count+1;
                         );
         );

stop3=0;
         loop(iner2$(stop1=0 and stop3=0),
         stop4=0;
         gval("g1")=1-(xvalue(iter,"x1")+(sum(neigh,selectd(iter,neigh)*directions(neigh,"x1")))*hiner2(iner2));
         gval("g2")=(xvalue(iter,"x1")+(sum(neigh,selectd(iter,neigh)*directions(neigh,"x1")))*hiner2(iner2))-card(n);
         gval("g3")=1-((xvalue(iter,"x1")+(sum(neigh,selectd(iter,neigh)*directions(neigh,"x1")))*hiner2(iner2))-(xvalue(iter,"x2")+(sum(neigh,selectd(iter,neigh)*directions(neigh,"x2")))*hiner2(iner2)));
         gval("g4")=((xvalue(iter,"x1")+(sum(neigh,selectd(iter,neigh)*directions(neigh,"x1")))*hiner2(iner2))-(xvalue(iter,"x2")+(sum(neigh,selectd(iter,neigh)*directions(neigh,"x2")))*hiner2(iner2)))-(xvalue(iter,"x1")+(sum(neigh,selectd(iter,neigh)*directions(neigh,"x1")))*hiner2(iner2));

                 if( sum(extineq$(gval(extineq) gt 0 ),gval(extineq)) gt 0  ,
                 stop4=1;
                 stop3=1;
                 );

                 loop(iner1$(stop4 eq 0),
                         if(ord(iner1) eq 1 and ord(iner2) eq 1 ,
                         execute_loadpoint "DLR";
                         elseif ord(iner1) eq 2 and ord(iner2) eq 1,
                         execute_loadpoint "DLR";
                         else
                         execute_loadpoint "DLR2";
                         );


                           if(ord(iner1) eq 1 ,
***S7_minimizacion en la direccion seleccionada cuando no hay problemas de convergencia
                           yf(n)=0;
                           yf(n)$(ord(n) eq floor(xvalue(iter,"x1")+(sum(neigh,selectd(iter,neigh)*directions(neigh,"x1")))*hiner2(iner2)))=1-mod(xvalue(iter,"x1")+(sum(neigh,selectd(iter,neigh)*directions(neigh,"x1")))*hiner2(iner2),1);
                           yf(n)$(ord(n) eq 1+floor(xvalue(iter,"x1")+(sum(neigh,selectd(iter,neigh)*directions(neigh,"x1")))*hiner2(iner2)))=mod(xvalue(iter,"x1")+(sum(neigh,selectd(iter,neigh)*directions(neigh,"x1")))*hiner2(iner2),1);
                           yr(n)=0;
                           yr(n)$(ord(n) eq floor((xvalue(iter,"x1")+(sum(neigh,selectd(iter,neigh)*directions(neigh,"x1")))*hiner2(iner2))-(xvalue(iter,"x2")+(sum(neigh,selectd(iter,neigh)*directions(neigh,"x2")))*hiner2(iner2))))=1-mod((xvalue(iter,"x1")+(sum(neigh,selectd(iter,neigh)*directions(neigh,"x1")))*hiner2(iner2))-(xvalue(iter,"x2")+(sum(neigh,selectd(iter,neigh)*directions(neigh,"x2")))*hiner2(iner2)),1);
                           yr(n)$(ord(n) eq 1+floor((xvalue(iter,"x1")+(sum(neigh,selectd(iter,neigh)*directions(neigh,"x1")))*hiner2(iner2))-(xvalue(iter,"x2")+(sum(neigh,selectd(iter,neigh)*directions(neigh,"x2")))*hiner2(iner2))))=mod((xvalue(iter,"x1")+(sum(neigh,selectd(iter,neigh)*directions(neigh,"x1")))*hiner2(iner2))-(xvalue(iter,"x2")+(sum(neigh,selectd(iter,neigh)*directions(neigh,"x2")))*hiner2(iner2)),1);
                           yp(n)=1-(sum(n1$(ord(n1) <= ord(n)),yf(n1))-yf(n));

                           solve SUP_CSTR using nlp minimizing zobj;
                           mstatS7(iter,iner2,iner1)=SUP_CSTR.modelstat;
                           execute_unload "DLR2";

                                 if( mstatS7(iter,iner2,iner1) eq 3 or mstatS7(iter,iner2,iner1) eq 4 or mstatS7(iter,iner2,iner1) eq 5 or mstatS7(iter,iner2,iner1) eq 6 or mstatS7(iter,iner2,iner1) eq 11 or mstatS7(iter,iner2,iner1) eq 12 or mstatS7(iter,iner2,iner1) eq 13 or mstatS7(iter,iner2,iner1) eq 14 or mstatS7(iter,iner2,iner1) eq 18 or mstatS7(iter,iner2,iner1) eq 19  ,
                                 stop4=0;
                                 else
                                 stop4=1;
                                 fvalueiner(iter,iner2)=zobj.l;
                                 );
***S8_minimizacion en la direccion seleccionada cuando hay problemas de convergencia
                           else
                           yf(n)=0;
                           yf(n)$(ord(n) eq floor(xvalue(iter,"x1")+(sum(neigh,selectd(iter,neigh)*directions(neigh,"x1")))*(hiner2(iner2)+hiner1(iner1)-hstep)))=1-mod(xvalue(iter,"x1")+(sum(neigh,selectd(iter,neigh)*directions(neigh,"x1")))*(hiner2(iner2)+hiner1(iner1)-hstep),1);
                           yf(n)$(ord(n) eq 1+floor(xvalue(iter,"x1")+(sum(neigh,selectd(iter,neigh)*directions(neigh,"x1")))*(hiner2(iner2)+hiner1(iner1)-hstep)))=mod(xvalue(iter,"x1")+(sum(neigh,selectd(iter,neigh)*directions(neigh,"x1")))*(hiner2(iner2)+hiner1(iner1)-hstep),1);
                           yr(n)=0;
                           yr(n)$(ord(n) eq floor(xvalue(iter,"x1")+(sum(neigh,selectd(iter,neigh)*directions(neigh,"x1")))*(hiner2(iner2)+hiner1(iner1)-hstep)-(xvalue(iter,"x2")+(sum(neigh,selectd(iter,neigh)*directions(neigh,"x2")))*(hiner2(iner2)+hiner1(iner1)-hstep))))=1-mod(xvalue(iter,"x1")+(sum(neigh,selectd(iter,neigh)*directions(neigh,"x1")))*(hiner2(iner2)+hiner1(iner1)-hstep)-(xvalue(iter,"x2")+(sum(neigh,selectd(iter,neigh)*directions(neigh,"x2")))*(hiner2(iner2)+hiner1(iner1)-hstep)),1);
                           yr(n)$(ord(n) eq 1+floor(xvalue(iter,"x1")+(sum(neigh,selectd(iter,neigh)*directions(neigh,"x1")))*(hiner2(iner2)+hiner1(iner1)-hstep)-(xvalue(iter,"x2")+(sum(neigh,selectd(iter,neigh)*directions(neigh,"x2")))*(hiner2(iner2)+hiner1(iner1)-hstep))))=mod(xvalue(iter,"x1")+(sum(neigh,selectd(iter,neigh)*directions(neigh,"x1")))*(hiner2(iner2)+hiner1(iner1)-hstep)-(xvalue(iter,"x2")+(sum(neigh,selectd(iter,neigh)*directions(neigh,"x2")))*(hiner2(iner2)+hiner1(iner1)-hstep)),1);
                           yp(n)=1-(sum(n1$(ord(n1) <= ord(n)),yf(n1))-yf(n));

                           solve SUP_CSTR using nlp minimizing zobj;
                           mstatS8(iter,iner2,iner1)=SUP_CSTR.modelstat;
                           execute_unload "DLR2";

                                 if(ord(iner1) eq card(iner1),
                                 stop4=1;

                                         if(mstatS8(iter,iner2,iner1) eq 3 or mstatS8(iter,iner2,iner1) eq 4 or mstatS8(iter,iner2,iner1) eq 5 or mstatS8(iter,iner2,iner1) eq 6 or mstatS8(iter,iner2,iner1) eq 11 or mstatS8(iter,iner2,iner1) eq 12 or mstatS8(iter,iner2,iner1) eq 13 or mstatS8(iter,iner2,iner1) eq 14 or mstatS8(iter,iner2,iner1) eq 18 or mstatS8(iter,iner2,iner1) eq 19  ,
                                         fvalueiner(iter,iner2)=1E+10;
                                         else
                                         fvalueiner(iter,iner2)=zobj.l;
                                         );

                                 );

                           );
                 );

                 if(ord(iner2) ne 1 and  fvalueiner(iter,iner2) ge fvalueiner(iter,iner2-1),
                 stop3=1;
                 );
                 if(stop3=0 and stop1=0,
                 objval=zobj.l;
                 execute_unload "DLR";
                 execute_unload "EVR2";
                         loop(extvar,
                          xvalueselect(iter,extvar)=xvalue(iter,extvar)+(sum(neigh,selectd(iter,neigh)*directions(neigh,extvar)))*hiner2(iner2);
                         );
                 );


         );




         if(stop1 =0,
         xvalue(iter+1,extvar)= xvalueselect(iter,extvar);
         else
         xvalue(iter+1,extvar)=xvalue(iter,extvar);
         );



);
CPUtime=timeElapsed-CPUtimeactual;
CPUtimeactual=timeElapsed;


