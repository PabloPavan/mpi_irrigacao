% FUN��O PARA CALCULAR O POTENCIAL MATRICIAL COM BASE NO TETA
% FUN��O CHAMADA PELO PROGRAMA SOLOCILCOMP

function y=POT(X)
	%y � o potencial
	%X � o teta

a=0.674;%0.467; ;%0.526;%;  % coeficiente do Psi na equa��o do potencial */
m=0.768;%0.934; %%0.519;%;  % expoente na equa��o do potencial */
n=2/(1-m); % expoente na equa��o do potencial, Burdines */

if (X <=0)
   y=0.35;
else

	if (X <=1)
   	y=-10^(1/a*(X^(-1/m)-1)^(1/n)); % pot matricial 
	else
   	y=0;
	end
end



%function y=POT(X,Y1,Y2,Y3)
	%y � o potencial
	%X � o teta
    %Y1 � o a
    %Y2 � o m
    %Y3 � o n

	
%   y=-1/(10^((-1/Y1*(X^(-1/Y2)-1)^(1/Y3))*0.434294481)); % pot matricial 
  % y=-1/(10^((-1/Y1*(X^(-1/Y2)-1)^(1/Y3)))); % pot matricial	
%    if (abs(y) >= 30000)
%        y=-30000;
%    end