% FUNÇÃO PARA CALCULAR A CONDITIVIDADE HIDRÁULICA COM BASE NO TETA
% FUNÇÃO CHAMADA PELO PROGRAMA SOLOCILCOMP

function y=KA(X1,X2)
	%y é a condutividade térmica
	%X1 É O Ko
	%X2 É O TETA
a=0.467;%0.674;% %%0.526;%;  % coeficiente do Psi na equação do potencial */
m=0.934;%0.768;%; %%0.519;%;  % expoente na equação do potencial */
n=2/(1-m); % expoente na equação do potencial, Burdines */

if (X2 <=0)
   y=0;
else

	if (X2 <=1)
   	y=X1*X2^0.5*(1-(1-X2^(1/m))^m)^2; %Bordine%
	else
   	y=X1;
	end
end




%function y=KA(X1,X2,X3)
	%y é a condutividade térmica
	%X1 É O Kz OU Kr
	%X2 É O TETA
    %X3 é O m



%if (X2 <=0)
%   y=0;
%else

%	if (X2 <=1)
%   	y=X1*X2^0.5*(1-(1-X2^(1/X3))^X3)^2; %Bordine

%	else
%   	y=X1;
%	end
%end
