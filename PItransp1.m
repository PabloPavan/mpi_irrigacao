% TRANSPORTE DE ÁGUA NO SOLO EM COORDENDAS CILÍNDRICAS 

%--------DEFINIÇÃO DA MALHA ----------
clear; clc;
format long
nn=100;  %iterações temporais, varia com k
L=16;   % no. de pontos em r, varia com i
M=L;   % no. de pontos em z, varia com j
R=0.15;   % raio do tubo,m
h=0.33;   % altura do solo,m
dr=R/(L-1);
dz=h/(M-1);
dt=0.3; % horas

%-----------curva de retenção de água -----------
a=0.467; %0.674;%0.526;%;  % coeficiente do Psi na equação do potencial */
m=0.934; %0.768;%0.519;%;  % expoente na equação do potencial */
n=2/(1-m); % expoente na equação do potencial, Burdines 
Tr=0.07;  %teor de umidade residual, dimensional
Ts=0.76; %teor de umidade de saturação, dimensional
tempo=nn*dt
Ko=0.000004; %condutividade hidráulica de saturação

%-----------DADOS EXPERIMENTAIS -----------
%---t = 0   CONDIÇOES INICIAIS -------

T0=0.65;   % teor de umidade vol do solo em t=0 , dimensional
Tr=0.01;
Ts=0.76;
T00=(T0-Tr)/(Ts-Tr); %teor de umidade inicial adimensional

for j=1:M
   for i=1:L
       T(i,j)=T00;  %condição inicial para toda a malha
   end
end 

%---t = 0 h -------
To(1,2)=T00;To(1,5)=T00;To(1,8)=T00;
To(3,2)=T00;To(3,5)=T00;To(3,8)=T00;
To(6,2)=T00;To(6,5)=T00;To(6,8)=T00;
To(11,2)=T00;To(11,5)=T00;To(11,8)=T00;

%---t = 10 h -------
T1(1,2)=0.75*T00;T1(1,5)=0.75*T00;T1(1,8)=0.75*T00;
T1(3,2)=0.78*T00;T1(3,5)=0.8*T00;T1(3,8)=0.8*T00;
T1(6,2)=0.7*T00;T1(6,5)=0.8*T00;T1(6,8)=0.9*T00;
T1(11,2)=0.98*T00;T1(11,5)=0.85*T00;T1(11,8)=0.95*T00;

%---t = 20 h -------
T2(1,2)=0.5*T00;T2(1,5)=0.5*T00;T2(1,8)=0.5*T00;
T2(3,2)=0.65*T00;T2(3,5)=0.7*T00;T2(3,8)=0.75*T00;
T2(6,2)=0.45*T00;T2(6,5)=0.7*T00;T2(6,8)=0.8*T00;
T2(11,2)=0.97*T00;T2(11,5)=0.8*T00;T2(11,8)=0.95*T00;

%---t = 30 h -------
T3(1,2)=0.3*T00;T3(1,5)=0.3*T00;T3(1,8)=0.35*T00;
T3(3,2)=0.7*T00;T3(3,5)=0.6*T00;T3(3,8)=0.75*T00;
T3(6,2)=0.6*T00;T3(6,5)=0.6*T00;T3(6,8)=0.75*T00;
T3(11,2)=0.96*T00;T3(11,5)=0.75*T00;T3(11,8)=0.9*T00;



%============= CÁLCULO DOS PARÂMETROS DA EVAPORAÇÃO ==============

        %---AJUSTE DE CURVAS EM Z=0  EVAPORAÇÃO-----
n=4; % número de pontos
tt=[0;10;20;30];
TEv=[To(1,2);T1(1,2);T2(1,2);T3(1,2)];
for i=1:n
    A(i,1)=tt(i);
    A(i,2)=1;
end
x=inv(A'*A)*A'*TEv;


%------------PROBLEMA INVERSO E DIRETO---------------------
Menor=100000000;
nc=59;  %número de q testados
qo=0.035; %fluxo de transpiração mínimo (m3/hora)
qM=0.09; % %fluxo de transpiração máximo
dq=(qM-qo)/nc;

for ii=1:nc 	% contador dos chutes do PI
   q=qo+(ii-1)*dq;   
   for j=1:M
   	for i=1:L
       	T(i,j)=T00;  %condição inicial para toda a malha
   	end
	end 
   
for j=1:L
	for i=1:M
  		Sf(i,j)=0;  %fonte zero para toda a malha
	end
end 

   for kk=1:nn   % TEMPORAL
      te(kk)=(kk-1)*dt;
      qq=q*exp(-0.05*(te(kk)-12)^2);
      for i=2:L-1
      	Tsup=x(1)*te(kk)+x(2); %teta em z=0 pela evaporaçao
      	for i1=1:L
      		T(1,i1)=Tsup; %condição de fronteira em z=0
      	end  % do i1
         for j=2:M-1   
            
            
           
   		%------------FONTE = TRANSPIRAÇÃO ----------------
			q1=qq;
			q2=qq;
			q3=qq;
			
			Sf(2,2)=qq*(T(2,2)-Tr);
			Sf(3,2)=qq*(T(3,2)-Tr);
			Sf(4,2)=qq*(T(4,2)-Tr);Sf(4,3)=q1*(T(4,3)-Tr);
			Sf(5,2)=qq*(T(5,2)-Tr);Sf(5,3)=q1*(T(5,3)-Tr);Sf(5,4)=q2*(T(5,4)-Tr);
			Sf(6,2)=qq*(T(6,2)-Tr);Sf(6,3)=q1*(T(6,3)-Tr);Sf(6,4)=q2*(T(6,4)-Tr);Sf(6,5)=q3*(T(6,5)-Tr);
			Sf(7,2)=qq*(T(7,2)-Tr);Sf(7,3)=q1*(T(7,3)-Tr);Sf(7,4)=q2*(T(7,4)-Tr);
			Sf(8,2)=qq*(T(8,2)-Tr);Sf(8,3)=q1*(T(8,3)-Tr);
			Sf(9,2)=qq*(T(9,2)-Tr);
			Sf(10,2)=qq*(T(10,2)-Tr);

            
         	 r=(j-1)*dr+0.5*dr;
                             
             KE=KA(Ko,T(i,j+1));
             KW=KA(Ko,T(i,j-1));
             KP=KA(Ko,T(i,j));
             KN=KA(Ko,T(i+1,j));
			 	 KS=KA(Ko,T(i-1,j));
 
            AE=dt/dr*((KE-KW)/(4*dr)+KP/(2*r)+KP/dr);
        		AW=dt/dr*(-(KE-KW)/(4*dr)-KP/(2*r)+KP/dr);
        		AN=dt/(4*dz^2)*(KN-KS+4*KP);
        		AS=dt/(4*dz^2)*(-KN+KS+4*KP);
        		AP0=-2*dt*KP*(1/dr^2+1/dz^2);
               
            PME=POT(T(i,j+1));
            PMW=POT(T(i,j-1));
            PMP=POT(T(i,j));
            PMN=POT(T(i+1,j));
            PMS=POT(T(i-1,j));
            
            TN(i,j)=(AE*PME+AW*PMW+AN*PMN+AS*PMS+AP0*PMP)+(KN-KS)/(2*dz)+T(i,j)-Sf(i,j);
         
        end % do for do j
		  TN(i,1)=TN(i,2); % CONDIÇÃO DE FRONTEIRA EM r=R, isolamento*/
        TN(i,L)=TN(i,L-1); % CONDIÇÃO DE FRONTEIRA EM r=R, isolamento*/

     end % do for do i
     
     for j1=1:L % condição de fronteira em z=h, isolamento
         TN(M,j1)=TN(M-1,j1);
     end
     
     for i1=2:L    %envelhecimento da TN
         for j1=1:M
         	T(i1,j1)=TN(i1,j1);
         end
     end
     
     %-------PROBLEMA INVERSO-----------
     
     if ( abs(te(kk)-10) <= dt)
         SQE10=(T1(1,2)-T(1,2))^2+(T1(1,5)-T(1,5))^2+(T1(1,8)-T(1,8))^2+(T1(3,2)-T(3,2))^2+(T1(3,5)-T(3,5))^2+(T1(3,8)-T(3,8))^2+(T1(6,2)-T(6,2))^2+(T1(6,5)-T(6,5))^2+(T1(6,8)-T(6,8))^2+(T1(11,2)-T(11,2))^2+(T1(11,5)-T(11,5))^2+(T1(11,8)-T(11,8))^2;    
     end % do if do 10h
     if ( abs(te(kk)-20) <= dt)
       	SQE20=(T2(1,2)-T(1,2))^2+(T2(1,5)-T(1,5))^2+(T2(1,8)-T(1,8))^2+(T2(3,2)-T(3,2))^2+(T2(3,5)-T(3,5))^2+(T2(3,8)-T(3,8))^2+(T2(6,2)-T(6,2))^2+(T2(6,5)-T(6,5))^2+(T2(6,8)-T(6,8))^2+(T2(11,2)-T(11,2))^2+(T2(11,5)-T(11,5))^2+(T2(11,8)-T(11,8))^2;  
     end % do if do 20h
	  if ( abs(te(kk)-29) <= dt)
     		SQE30=(T3(1,2)-T(1,2))^2+(T3(1,5)-T(1,5))^2+(T3(1,8)-T(1,8))^2+(T3(3,2)-T(3,2))^2+(T3(3,5)-T(3,5))^2+(T3(3,8)-T(3,8))^2+(T3(6,2)-T(6,2))^2+(T3(6,5)-T(6,5))^2+(T3(6,8)-T(6,8))^2+(T3(11,2)-T(11,2))^2+(T3(11,5)-T(11,5))^2+(T3(11,8)-T(11,8))^2;
     end
             	
	end  % do for do kk

	Erro(ii)=SQE10+SQE20+SQE30;
   
   if(Erro(ii) <= Menor)       
		qot=q;
   	iot=ii;
   	Menor=Erro(ii);
	else
   	Menor=Menor;
   end

	

end % do ii, Problema inverso

%---escolha do menor Erro---

qot
iot



S1=(T1(1,2)+T1(1,5)+T1(1,8)+T1(3,2)+T1(3,5)+T1(3,8)+T1(6,2)+T1(6,5)+T1(6,8)+T1(11,2)+T1(11,5)+T1(11,8))/12;    
S2=(T2(1,2)+T2(1,5)+T2(1,8)+T2(3,2)+T2(3,5)+T2(3,8)+T2(6,2)+T2(6,5)+T2(6,8)+T2(11,2)+T2(11,5)+T2(11,8))/12;    
S3=(T3(1,2)+T3(1,5)+T3(1,8)+T3(3,2)+T3(3,5)+T3(3,8)+T3(6,2)+T3(6,5)+T3(6,8)+T3(11,2)+T3(11,5)+T3(11,8))/12;    

Med=(S1+S2+S3)/3;
SQT1=(T1(1,2)-Med)^2+(T1(1,5)-Med)^2+(T1(1,8)-Med)^2+(T1(3,2)-Med)^2+(T1(3,5)-Med)^2+(T1(3,8)-Med)^2+(T1(6,2)-Med)^2+(T1(6,5)-Med)^2+(T1(6,8)-Med)^2+(T1(11,2)-Med)^2+(T1(11,5)-Med)^2+(T1(11,8)-Med)^2;    
SQT2=(T2(1,2)-Med)^2+(T2(1,5)-Med)^2+(T2(1,8)-Med)^2+(T2(3,2)-Med)^2+(T2(3,5)-Med)^2+(T2(3,8)-Med)^2+(T2(6,2)-Med)^2+(T2(6,5)-Med)^2+(T2(6,8)-Med)^2+(T2(11,2)-Med)^2+(T2(11,5)-Med)^2+(T2(11,8)-Med)^2;  
SQT3=(T3(1,2)-Med)^2+(T3(1,5)-Med)^2+(T3(1,8)-Med)^2+(T3(3,2)-Med)^2+(T3(3,5)-Med)^2+(T3(3,8)-Med)^2+(T3(6,2)-Med)^2+(T3(6,5)-Med)^2+(T3(6,8)-Med)^2+(T3(11,2)-Med)^2+(T3(11,5)-Med)^2+(T3(11,8)-Med)^2;

SQT=(SQT1+SQT2+SQT3)^2;
R2=1-Erro(iot)/SQT

%=====SOLUÇÃO ÓTIMA======

q=qot; 
   %------------FONTE = TRANSPIRAÇÃO ----------------
	q1=q;
	q2=q;
	q3=q;
	for j=1:L
   	for i=1:M
         Sf(i,j)=0;  %condição inicial para toda a malha
   	end
	end 
	Sf(2,2)=q;
	Sf(3,2)=q;
	Sf(4,2)=q;Sf(4,3)=q1;
	Sf(5,2)=q;Sf(5,3)=q1;Sf(5,4)=q2;
	Sf(6,2)=q;Sf(6,3)=q1;Sf(6,4)=q2;Sf(6,5)=q3;
	Sf(7,2)=q;Sf(7,3)=q1;Sf(7,4)=q2;
	Sf(8,2)=q;Sf(8,3)=q1;
	Sf(9,2)=q;
	Sf(10,2)=q;
   
   for j=1:M
   	for i=1:L
       	T(i,j)=T00;  %condição inicial para toda a malha
   	end
	end 

   for kk=1:nn   % TEMPORAL
      te(kk)=(kk-1)*dt;
      qq=q*exp(-0.05*(te(kk)-12)^2);

   	for i=2:L-1
      	Tsup=x(1)*te(kk)+x(2); %teta em z=0 pela evaporaçao
      	for i1=1:L
      		T(1,i1)=Tsup; %condição de fronteira em z=0
      	end  % do i1
         for j=2:M-1    
            
            
            %------------FONTE = TRANSPIRAÇÃO ----------------
			q1=qq;
			q2=qq;
			q3=qq;
			
			Sf(2,2)=qq*(T(2,2)-Tr);
			Sf(3,2)=qq*(T(3,2)-Tr);
			Sf(4,2)=qq*(T(4,2)-Tr);Sf(4,3)=q1*(T(4,3)-Tr);
			Sf(5,2)=qq*(T(5,2)-Tr);Sf(5,3)=q1*(T(5,3)-Tr);Sf(5,4)=q2*(T(5,4)-Tr);
			Sf(6,2)=qq*(T(6,2)-Tr);Sf(6,3)=q1*(T(6,3)-Tr);Sf(6,4)=q2*(T(6,4)-Tr);Sf(6,5)=q3*(T(6,5)-Tr);
			Sf(7,2)=qq*(T(7,2)-Tr);Sf(7,3)=q1*(T(7,3)-Tr);Sf(7,4)=q2*(T(7,4)-Tr);
			Sf(8,2)=qq*(T(8,2)-Tr);Sf(8,3)=q1*(T(8,3)-Tr);
			Sf(9,2)=qq*(T(9,2)-Tr);
			Sf(10,2)=qq*(T(10,2)-Tr);

         	 r=(j-1)*dr+0.5*dr;
                             
             KE=KA(Ko,T(i,j+1));
             KW=KA(Ko,T(i,j-1));
             KP=KA(Ko,T(i,j));
             KN=KA(Ko,T(i+1,j));
			 	 KS=KA(Ko,T(i-1,j));
 
            AE=dt/dr*((KE-KW)/(4*dr)+KP/(2*r)+KP/dr);
        		AW=dt/dr*(-(KE-KW)/(4*dr)-KP/(2*r)+KP/dr);
        		AN=dt/(4*dz^2)*(KN-KS+4*KP);
        		AS=dt/(4*dz^2)*(-KN+KS+4*KP);
        		AP0=-2*dt*KP*(1/dr^2+1/dz^2);
               
            PME=POT(T(i,j+1));
            PMW=POT(T(i,j-1));
            PMP=POT(T(i,j));
            PMN=POT(T(i+1,j));
            PMS=POT(T(i-1,j));
            
            TN(i,j)=(AE*PME+AW*PMW+AN*PMN+AS*PMS+AP0*PMP)+(KN-KS)/(2*dz)+T(i,j)-Sf(i,j);
         
        end % do for do j
		  TN(i,1)=TN(i,2); % CONDIÇÃO DE FRONTEIRA EM r=R, isolamento*/
        TN(i,L)=TN(i,L-1); % CONDIÇÃO DE FRONTEIRA EM r=R, isolamento*/

     end % do for do i
     
     for j1=1:L % condição de fronteira em z=h, isolamento
         TN(M,j1)=TN(M-1,j1);
     end
     
     for i1=2:L    %envelhecimento da TN
         for j1=1:M
         	T(i1,j1)=TN(i1,j1);
         end
     end
                  
     %--------DADOS PARA GRÁFICO DAS CURVAS ------------- 
    
	n1=fix(0.08/dz);
	n2=fix(0.15/dz);
	n3=fix(0.25/dz);
	m3=fix(0.02/dr);

	C1(kk)=T(1,m3)   ;%T(1,6);   
	C2(kk)=T(n1,m3);    
	C3(kk)=T(n2,m3); 
	C4(kk)=T(n3,m3); 

	
	end  % do for do kk

%-------------GRAFICOS--------------

figure(1);
r=linspace(0,R,L);
y=linspace(0,h,M);
[X,Y]=meshgrid(r,y);
surf(X,Y,T);
xlabel('raio(m)');
ylabel('altura(m)');
zlabel('teor de água (ad)');

figure(2);
contour(X,Y,TN,20);
xlabel('raio(m)');
ylabel('altura(m)');

%--Dados experimentais
T12=[To(1,2),T1(1,2),T2(1,2),T3(1,2)];
T32=[To(3,2),T1(3,2),T2(3,2),T3(3,2)];
T62=[To(6,2),T1(6,2),T2(6,2),T3(6,2)];
T112=[To(11,2),T1(11,2),T2(11,2),T3(11,2)];

figure(3);  %curva teta X tempo
tp=linspace(0,tempo,nn);
plot(tt,T12,'or',tt,T32,'xg',tt,T62,'*',tt,T112,'+k',tp,C1,'r',tp,C2,'g',tp,C3,'b',tp,C4,'k');
xlabel('tempo (horas)');
ylabel('teor de água (ad)');
legend('0 m exp','0.08m exp','0.15m exp','0.25m exp','0 m cal','0.08m cal','0.15m cal','0.25m cal');






