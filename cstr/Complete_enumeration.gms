* Code to replicate Figure 3 in Optimal design of superstructures for placing 
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
yr(n)$(ord(n) eq floor(variablex2.l))=1-mod(variablex2.l,1);
yr(n)$(ord(n) eq 1+floor(variablex2.l))=mod(variablex2.l,1);
yp(n)=1-(sum(n1$(ord(n1) <= ord(n)),yf(n1))-yf(n));
solve SUP_CSTR using nlp minimizing zobj;
execute_unload "in";

*--------------------SOLUTION FOR EVERY CONFIGURATION OF EXTERNAL VARIABLES-----
set case /1*15/;
parameter par_x1(case)
/
1        1
2        2
3        3
4        4
5        5
6        2
7        3
8        4
9        5
10        3
11        4
12        5
13        4
14        5
15        5
/
;
parameter par_x2(case)
/
1        1
2        1
3        1
4        1
5        1
6        2
7        2
8        2
9        2
10        3
11        3
12        3
13        4
14        4
15        5
/
;
parameter objval(case);
parameter mstat(case);
parameter QRval(case);


loop(case,
yf(n)=0;
yf(n)$(ord(n) eq par_x1(case))=1;
yr(n)=0;
yr(n)$(ord(n) eq par_x2(case))=1;
yp(n)=1-(sum(n1$(ord(n1) <= ord(n)),yf(n1))-yf(n));



execute_loadpoint "in";
solve SUP_CSTR using nlp minimizing zobj;
mstat(case)=SUP_CSTR.modelstat;
objval(case)=zobj.l;
QRval(case)=QR.l;

);
execute_unload "Enumeration"

