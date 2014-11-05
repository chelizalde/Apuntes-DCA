%----------------------------
%    Control optimo
%---------------------------

clear 
clc 
syms s t

Id=eye(3);
vcero=zeros(3,1);
%--------------------------
% Representacion de estado
%--------------------------

A= [0 1 0;
	0 0 1;
	0 0 0]
b=[0; 0; 1]
c=[-1;0;1]
xo=[1;0;0]

%-------------------------
% Funcion de transferencia variable de estado
%-------------------------

FTv=simplify(c.'*inv(s*Id-A)*b)

%--------------------------
% Polos --> valores propios
%--------------------------

vpA = eig(A)

%-------------------------
% Ceros --> raices del determinante de la Matriz Sistema
%-------------------------

MS=[(s*Id-A) b; -c.' 0]
pMS=factor(det(MS))

%-------------------------
% Matriz de contrlabilidad
%-------------------------

MC=[b A*b A^2*b]
dMC = det(MC)

%-------------------------
% Matriz de obsevabilidad
%-------------------------

MO=[c.'; c.'*A;c.'*A^2]
dMO=det(MO)

%-------------------------
% Proposicion de asignacion de polos FT=((s-1)(s+1))
%-------------------------

fp= [-1;-3;-3]

% Polos lazo cerrado

Afp=A+b*fp.'
vpAfp=eig(Afp)

% Ceros lazo cerrado

MSLcp=[(s*Id-Afp) b; -c.' 0]
pMSLcp=factor(det(MSLcp))

% Forma de Jordan del Sistema en lazo cerrado 

[Tp,Jfp]=jordan(Afp)
inv(Tp)*Afp*Tp
bJp=inv(Tp)*b
cJp=Tp.'*c

% Matriz sistem del sistema en lazo cerrado en su forma de jordan

MSLcJp=[(s*Id-Jfp) bJp; -cJp.' 0]
pretty(MSLcJp)

% Funcion de transferencia en lazo cerrado

FTLcp=factor(cJp.'*inv(s*Id-Jfp)*bJp)

%-------------------------
% Proposicion Control optimo : J = \int_{0}^{\infty}(y^2 + \rho u^2)
%-------------------------

Q=c*c.'
vpQ=eig(Q)
r=1
%r=10
%r=1/10
%-------------------------
% Ecuacion de ricati : A.'*P+ P*A - P*b*r^(-1)*b.' + Q =0
%-------------------------

[P,L,G]=care(A,b,Q,r,vcero,Id)
vpP=eig(P)

%Retroalimentacion de estado control optimo

fo=G.'

%Polos lazo cerrado

Afo=A-b*fo.'
pcLco=factor(det(s*Id-Afo)) 
vpAfo=eig(Afo)

%Ceros lazo cerrado

MSLco = [(s*Id-Afo) b; -c.' 0]
pMSLco = factor(det(MSLco))

% Forma de jordan del sistema en lazo cerrado

[To,Jfo]=jordan(Afo)

inv(To)*Afo*To
bJo = inv(To)*b
cJo= To.'*c 

%-------------------------
% Matriz sistema del sistema en lazo cerrado en su forma de Jordan
%-------------------------


MSLcJo = vpa([(s*Id-Jfo) bJo; -cJo.' 0],2)
pretty(MSLcJo)

%Funcion de transferencia en lazo cerrado

FTLco = vpa(factor(cJo.'*inv(s*Id-Jfo)*bJo),2)


%-------------------
%Comparacion polos en lazo cerrado

vpAfp
r
vpAfo
%----------------------
%Realizar simulacion mdl (simulink)



