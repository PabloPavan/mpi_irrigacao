% FUN��O PARA CALCULAR A CONDITIVIDADE HIDR�ULICA COM BASE NO TETA
% FUN��O CHAMADA PELO PROGRAMA SOLOCILCOMP

function y=KA(X1,X2)
	%y � a condutividade t�rmica
	%X1 � O Ko
	%X2 � O TETA
a=0.467;%0.674;% %%0.526;%;  % coeficiente do Psi na equa��o do potencial */
m=0.934;%0.768;%; %%0.519;%;  % expoente na equa��o do potencial */
n=2/(1-m); % expoente na equa��o do potencial, Burdines */

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
	%y � a condutividade t�rmica
	%X1 � O Kz OU Kr
	%X2 � O TETA
    %X3 � O m



%if (X2 <=0)
%   y=0;
%else

%	if (X2 <=1)
%   	y=X1*X2^0.5*(1-(1-X2^(1/X3))^X3)^2; %Bordine

%	else
%   	y=X1;
%	end
%end
