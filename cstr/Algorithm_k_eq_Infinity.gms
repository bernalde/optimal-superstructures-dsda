* Code to replicate Figure 4a in Optimal design of superstructures for placing 
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
feas..variablex2-variablex1=l=0;
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
variablex2.lo=1;
variablex2.up=card(N);
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
***NLP sub-solver for the algorithm
option optcr = 0;
option threads = 1;
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
yr(n)$(ord(n) eq floor(variablex2.l))=1-mod(variablex2.l,1);
yr(n)$(ord(n) eq 1+floor(variablex2.l))=mod(variablex2.l,1);
yp(n)=1-(sum(n1$(ord(n1) <= ord(n)),yf(n1))-yf(n));
solve SUP_CSTR using nlp minimizing zobj;
execute_unload "in";

*--------------------SOLUTION WITH FEASIBLE INITIALIZATION----------------------
***Iterations Sets
*iter: Select the maximum value of iterations allowed. If the polyhedron is known
*to be bounded, this should be as large as possible.
set iter "Iterations" /i1*i100/;
*iner1: This corresponds to the set It=In.
set iner1 "Inner iterations 1" /i1*i10/;
*iner2: Select the maximum value of inner iterations in the line search. If iner2
*=/i1/, the line search procedure wont be executed.
set iner2 "Inner iterations 2" /i1*i100/;
set extvar "External variables" /x1*x2/;
set extineq "Convex polyhedron inequality equations" /g1*g4/;

***Neighborhood
*In this case, the infinity neighborhood is used.
set neigh /d1*d8/;
table directions(neigh,extvar)
    x1       x2
d1 -1        -1
d2 -1        0
d3 -1        1
d4 0        -1
d5 0        1
d6 1        -1
d7 1        1
d8 1        0
;

***Steps defintion
*hstep: So far, the infinity-norm and 2-norms have been compared with an hstep of
*1. Considering other values of hstep is left as future work
scalar hstep "step" /1/;
parameter hiner1(iner1) "step: inner iterations 1";
loop(iner1,
if(ord(iner1) eq 1,
hiner1(iner1)=0;
else
hiner1(iner1)=hiner1(iner1-1)+(hstep/((card(iner1))-1));
);
);
parameter hiner2(iner2) "step: inner iterations 2";
loop(iner2,
if(ord(iner2) eq 1,
hiner2(iner2)=hstep;
else
hiner2(iner2)=hiner2(iner2-1)+hstep;
);
);
parameter extvarhstep(extvar,neigh) "Neighborhood considering hstep";
extvarhstep(extvar,neigh)=hstep*directions(neigh,extvar);

***Iterative process parameters and scalars
*stop1: stopping criterion in S4 for local optimality
scalar stop1 "Main stopping criterion";
*stop2: Stopping criterion that decides if the objective function value
*of a neighbor was properly calculated (e.g, if the nlp solver finds an infeasible
*solution, the additional convergence procedure is executed).
scalar stop2 "Stoppong criterion 2";
*stop3: Stopping criterion for the line search
scalar stop3 "Stoppong criterion 3";
*stop4: Stopping criterion that decides if the objective function value
*in the line search was properly calculated.
scalar stop4 "Stoppong criterion 4";
scalar count;
parameter xvalue(iter,extvar) "value of x at each iteration";
parameter xvalueselect(iter,extvar) "selected value of x for the next iteration";
parameter gval(extineq) "value of inequality constraints at x";
parameter dvs(iter,neigh) "f(neigh)-f(x) at the iteration";
parameter dmin(iter) "Minimum value of dvs";
parameter fvalue(iter) "Objective function value";
parameter fplushvalue(iter,neigh) "Objective function to calculate dvs";
parameter fvalueiner(iter,iner2) "Objective function: inner iterations";
parameter selectd(iter,neigh) "ds: Selected direction for line search";
parameter mstatS2(iter)"model status: Feasibility of initialization point";
parameter mstatS4(iter,neigh,iner1) "model status: Feasibility of neighbors";
parameter mstatS5(iter,neigh,iner1) "model status: Feasibility of neighnors with convergence procedure" ;
parameter mstatS7(iter,iner2,iner1) "model status: Feasibility in line search";
parameter mstatS8(iter,iner2,iner1) "model status: Feasibility in line search with convergence procedure";
parameter objval "Objective function value";
parameter xvalinit(extvar) "Initialization of external variables";
parameter CPUtime;
scalar CPUtimeactual;
parameter Starttime;
set step /1*6/
parameter timesneigh(step,iter,neigh) "time takes to solve each step, indexed by iteration and neighbor";
parameter timesiner(step,iter,iner2) "time takes to solve each step, indexed by iteration and inner line search";

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


***S1. Initialization of external variabkes
xvalinit("x1")=variablex1.l;
xvalinit("x2")=variablex2.l;
xvalue('i1',extvar)=xvalinit(extvar);


stop1=0;
loop(iter$((stop1 eq 0) and (ord(iter) ne card(iter))),

***S2. Feasibility of external variables
         if(ord(iter) eq 1,
         gval("g1")=1-xvalue(iter,"x1");
         gval("g2")=xvalue(iter,"x1")-card(n);
         gval("g3")=1-xvalue(iter,"x2");
         gval("g4")=xvalue(iter,"x2")-xvalue(iter,"x1");

                 if( sum(extineq$(gval(extineq) gt 0 ),gval(extineq)) gt 0,
                 stop1=1;
                 );

         );
***S1. Initialization if continuous variables
         if(ord(iter) eq 1,
         execute_loadpoint "in";
         else
         execute_loadpoint "DLR";
         );
***S2. Feasibility of continuous variables
yf(n)=0;
yf(n)$(ord(n) eq floor(xvalue(iter,"x1")))=1-mod(xvalue(iter,"x1"),1);
yf(n)$(ord(n) eq 1+floor(xvalue(iter,"x1")))=mod(xvalue(iter,"x1"),1);
yr(n)=0;
yr(n)$(ord(n) eq floor(xvalue(iter,"x2")))=1-mod(xvalue(iter,"x2"),1);
yr(n)$(ord(n) eq 1+floor(xvalue(iter,"x2")))=mod(xvalue(iter,"x2"),1);
yp(n)=1-(sum(n1$(ord(n1) <= ord(n)),yf(n1))-yf(n));

Starttime = timeElapsed;
solve SUP_CSTR using nlp minimizing zobj;
timesneigh('2','i1','d1') = timeElapsed - Starttime;
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
***S3. k-Neighborhood

                 gval("g1")=1-(xvalue(iter,"x1")+extvarhstep("x1",neigh));
                 gval("g2")=(xvalue(iter,"x1")+extvarhstep("x1",neigh))-card(n);
                 gval("g3")=1-(xvalue(iter,"x2")+extvarhstep("x2",neigh));
                 gval("g4")=(xvalue(iter,"x2")+extvarhstep("x2",neigh))-(xvalue(iter,"x1")+extvarhstep("x1",neigh));

                          if( sum(extineq$(gval(extineq) gt 0 ),gval(extineq)) gt 0 ,
                          dvs(iter,neigh)=1;
                          else
                          dvs(iter,neigh)=-1;
                          );

                          if(dvs(iter,neigh) le 0,
                          stop2=0;

                                 loop(iner1$(stop2 eq 0),

                                         if(ord(iner1) eq 1,
***S4. Local optimality: Computation of the objective function value for the neighbors
                                         execute_loadpoint "DLR";

                                         yf(n)=0;
                                         yf(n)$(ord(n) eq floor(xvalue(iter,"x1")+extvarhstep("x1",neigh)))=1-mod(xvalue(iter,"x1")+extvarhstep("x1",neigh),1);
                                         yf(n)$(ord(n) eq 1+floor(xvalue(iter,"x1")+extvarhstep("x1",neigh)))=mod(xvalue(iter,"x1")+extvarhstep("x1",neigh),1);
                                         yr(n)=0;
                                         yr(n)$(ord(n) eq floor(xvalue(iter,"x2")+extvarhstep("x2",neigh)))=1-mod(xvalue(iter,"x2")+extvarhstep("x2",neigh),1);
                                         yr(n)$(ord(n) eq 1+floor(xvalue(iter,"x2")+extvarhstep("x2",neigh)))=mod(xvalue(iter,"x2")+extvarhstep("x2",neigh),1);
                                         yp(n)=1-(sum(n1$(ord(n1) <= ord(n)),yf(n1))-yf(n));

                                         Starttime = timeElapsed;
                                         solve SUP_CSTR using nlp minimizing zobj;
                                         timesneigh('4',iter,neigh) = timeElapsed - Starttime;
                                         mstatS4(iter,neigh,iner1)=SUP_CSTR.modelstat;

                                                 if(mstatS4(iter,neigh,iner1) eq 3 or mstatS4(iter,neigh,iner1) eq 4 or mstatS4(iter,neigh,iner1) eq 5 or mstatS4(iter,neigh,iner1) eq 6 or mstatS4(iter,neigh,iner1) eq 11 or mstatS4(iter,neigh,iner1) eq 12 or mstatS4(iter,neigh,iner1) eq 13 or mstatS4(iter,neigh,iner1) eq 14 or mstatS4(iter,neigh,iner1) eq 18 or mstatS4(iter,neigh,iner1) eq 19,
                                                 stop2=0;
                                                 else
                                                 stop2=1;
                                                 fplushvalue(iter,neigh)=zobj.l;
                                                 dvs(iter,neigh)=(fplushvalue(iter,neigh)-fvalue(iter))/(hstep);
                                                 );

                                         else
***S4. Local optimality: Computation of the objective function value for the neighbors with the convergence procedure

                                                 if(ord(iner1) eq 2,
                                                 execute_loadpoint "DLR";
                                                 else
                                                 execute_loadpoint "DLR1";
                                                 );
                                         yf(n)=0;
                                         yf(n)$(ord(n) eq floor(xvalue(iter,"x1")+((extvarhstep("x1",neigh))/(hstep))*(hiner1(iner1))))=1-mod(xvalue(iter,"x1")+((extvarhstep("x1",neigh))/(hstep))*(hiner1(iner1)),1);
                                         yf(n)$(ord(n) eq 1+floor(xvalue(iter,"x1")+((extvarhstep("x1",neigh))/(hstep))*(hiner1(iner1))))=mod(xvalue(iter,"x1")+((extvarhstep("x1",neigh))/(hstep))*(hiner1(iner1)),1);
                                         yr(n)=0;
                                         yr(n)$(ord(n) eq floor(xvalue(iter,"x2")+((extvarhstep("x2",neigh))/(hstep))*(hiner1(iner1))))=1-mod(xvalue(iter,"x2")+((extvarhstep("x2",neigh))/(hstep))*(hiner1(iner1)),1);
                                         yr(n)$(ord(n) eq 1+floor(xvalue(iter,"x2")+((extvarhstep("x2",neigh))/(hstep))*(hiner1(iner1))))=mod(xvalue(iter,"x2")+((extvarhstep("x2",neigh))/(hstep))*(hiner1(iner1)),1);
                                         yp(n)=1-(sum(n1$(ord(n1) <= ord(n)),yf(n1))-yf(n));

                                         Starttime = timeElapsed;
                                         solve SUP_CSTR using nlp minimizing zobj;
                                         timesneigh('4',iter,neigh) = timeElapsed - Starttime;
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

***S4. Local optimality: Stopping criterion verification
dvs(iter,neigh)=round(dvs(iter,neigh),8);
dmin(iter)=smin((neigh),dvs(iter,neigh));
         if(dmin(iter)  ge 0 ,
         stop1=1;
         );
***S5. Steepest descent: Selection of ds
count=0;
         loop(neigh,
                         if(dvs(iter,neigh) eq dmin(iter) and count eq 0,
                         selectd(iter,neigh)=1;
                         count=count+1;
                         );
         );

stop3=0;
         loop(iner2$(stop1=0 and stop3=0),
***S6. Line search: Feasibility of external variables
         stop4=0;
         gval("g1")=1-(xvalue(iter,"x1")+(sum(neigh,selectd(iter,neigh)*directions(neigh,"x1")))*hiner2(iner2));
         gval("g2")=(xvalue(iter,"x1")+(sum(neigh,selectd(iter,neigh)*directions(neigh,"x1")))*hiner2(iner2))-card(n);
         gval("g3")=1-(xvalue(iter,"x2")+(sum(neigh,selectd(iter,neigh)*directions(neigh,"x2")))*hiner2(iner2));
         gval("g4")=(xvalue(iter,"x2")+(sum(neigh,selectd(iter,neigh)*directions(neigh,"x2")))*hiner2(iner2))-(xvalue(iter,"x1")+(sum(neigh,selectd(iter,neigh)*directions(neigh,"x1")))*hiner2(iner2));

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
***S6. Line search: Computation of the objective function value
                           yf(n)=0;
                           yf(n)$(ord(n) eq floor(xvalue(iter,"x1")+(sum(neigh,selectd(iter,neigh)*directions(neigh,"x1")))*hiner2(iner2)))=1-mod(xvalue(iter,"x1")+(sum(neigh,selectd(iter,neigh)*directions(neigh,"x1")))*hiner2(iner2),1);
                           yf(n)$(ord(n) eq 1+floor(xvalue(iter,"x1")+(sum(neigh,selectd(iter,neigh)*directions(neigh,"x1")))*hiner2(iner2)))=mod(xvalue(iter,"x1")+(sum(neigh,selectd(iter,neigh)*directions(neigh,"x1")))*hiner2(iner2),1);
                           yr(n)=0;
                           yr(n)$(ord(n) eq floor(xvalue(iter,"x2")+(sum(neigh,selectd(iter,neigh)*directions(neigh,"x2")))*hiner2(iner2)))=1-mod(xvalue(iter,"x2")+(sum(neigh,selectd(iter,neigh)*directions(neigh,"x2")))*hiner2(iner2),1);
                           yr(n)$(ord(n) eq 1+floor(xvalue(iter,"x2")+(sum(neigh,selectd(iter,neigh)*directions(neigh,"x2")))*hiner2(iner2)))=mod(xvalue(iter,"x2")+(sum(neigh,selectd(iter,neigh)*directions(neigh,"x2")))*hiner2(iner2),1);
                           yp(n)=1-(sum(n1$(ord(n1) <= ord(n)),yf(n1))-yf(n));

                           Starttime = timeElapsed;
                           solve SUP_CSTR using nlp minimizing zobj;
                           timesiner('6',iter,iner2) = timeElapsed - Starttime;
                           mstatS7(iter,iner2,iner1)=SUP_CSTR.modelstat;
                           execute_unload "DLR2";

                                 if( mstatS7(iter,iner2,iner1) eq 3 or mstatS7(iter,iner2,iner1) eq 4 or mstatS7(iter,iner2,iner1) eq 5 or mstatS7(iter,iner2,iner1) eq 6 or mstatS7(iter,iner2,iner1) eq 11 or mstatS7(iter,iner2,iner1) eq 12 or mstatS7(iter,iner2,iner1) eq 13 or mstatS7(iter,iner2,iner1) eq 14 or mstatS7(iter,iner2,iner1) eq 18 or mstatS7(iter,iner2,iner1) eq 19  ,
                                 stop4=0;
                                 else
                                 stop4=1;
                                 fvalueiner(iter,iner2)=zobj.l;
                                 );
***S6. Line search: Computation of the objective function value  with the convergence procedure
                           else
                           yf(n)=0;
                           yf(n)$(ord(n) eq floor(xvalue(iter,"x1")+(sum(neigh,selectd(iter,neigh)*directions(neigh,"x1")))*(hiner2(iner2)+hiner1(iner1)-hstep)))=1-mod(xvalue(iter,"x1")+(sum(neigh,selectd(iter,neigh)*directions(neigh,"x1")))*(hiner2(iner2)+hiner1(iner1)-hstep),1);
                           yf(n)$(ord(n) eq 1+floor(xvalue(iter,"x1")+(sum(neigh,selectd(iter,neigh)*directions(neigh,"x1")))*(hiner2(iner2)+hiner1(iner1)-hstep)))=mod(xvalue(iter,"x1")+(sum(neigh,selectd(iter,neigh)*directions(neigh,"x1")))*(hiner2(iner2)+hiner1(iner1)-hstep),1);
                           yr(n)=0;
                           yr(n)$(ord(n) eq floor(xvalue(iter,"x2")+(sum(neigh,selectd(iter,neigh)*directions(neigh,"x2")))*(hiner2(iner2)+hiner1(iner1)-hstep)))=1-mod(xvalue(iter,"x2")+(sum(neigh,selectd(iter,neigh)*directions(neigh,"x2")))*(hiner2(iner2)+hiner1(iner1)-hstep),1);
                           yr(n)$(ord(n) eq 1+floor(xvalue(iter,"x2")+(sum(neigh,selectd(iter,neigh)*directions(neigh,"x2")))*(hiner2(iner2)+hiner1(iner1)-hstep)))=mod(xvalue(iter,"x2")+(sum(neigh,selectd(iter,neigh)*directions(neigh,"x2")))*(hiner2(iner2)+hiner1(iner1)-hstep),1);
                           yp(n)=1-(sum(n1$(ord(n1) <= ord(n)),yf(n1))-yf(n));

                           Starttime = timeElapsed;
                           solve SUP_CSTR using nlp minimizing zobj;
                           timesiner('6',iter,iner2) = timeElapsed - Starttime;
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
***S6. Line search: Stopping criterion
                 if(ord(iner2) ne 1 and  fvalueiner(iter,iner2) ge fvalueiner(iter,iner2-1),
                 stop3=1;
                 );
                 if(stop3=0 and stop1=0,
                 objval=zobj.l;
                 execute_unload "DLR";
                 execute_unload "Alg_k_eq_infinity";
                         loop(extvar,
                          xvalueselect(iter,extvar)=xvalue(iter,extvar)+(sum(neigh,selectd(iter,neigh)*directions(neigh,extvar)))*hiner2(iner2);
                         );
                 );


         );



***Initial value of x for the next iteration
         if(stop1 =0,
         xvalue(iter+1,extvar)= xvalueselect(iter,extvar);
         else
         CPUtimeactual=timeElapsed;
         xvalue(iter+1,extvar)=xvalue(iter,extvar);
         yf(n)=0;
         yf(n)$(ord(n) eq floor(xvalue(iter,"x1")))=1-mod(xvalue(iter,"x1"),1);
         yf(n)$(ord(n) eq 1+floor(xvalue(iter,"x1")))=mod(xvalue(iter,"x1"),1);
         yr(n)=0;
         yr(n)$(ord(n) eq floor(xvalue(iter,"x2")))=1-mod(xvalue(iter,"x2"),1);
         yr(n)$(ord(n) eq 1+floor(xvalue(iter,"x2")))=mod(xvalue(iter,"x2"),1);
         yp(n)=1-(sum(n1$(ord(n1) <= ord(n)),yf(n1))-yf(n));

         solve SUP_CSTR using nlp minimizing zobj;
         );



);
CPUtime=timeElapsed-CPUtimeactual;
execute_unload "Alg_k_eq_infinity_INFO" CPUtimeactual, CPUtime, dvs, fvalue, fplushvalue, fvalueiner, selectd, xvalue, F, rate, V, P, zobj, timesneigh, timesiner;


