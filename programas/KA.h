//FUN�O PARA CALCULAR A CONDITIVIDADE HIDR�LICA COM BASE NO TETA
//FUN�O CHAMADA PELO PROGRAMA SOLOCILCOMP

//FUN�O PARA CALCULAR A CONDITIVIDADE HIDR�LICA COM BASE NO TETA
//FUN�O CHAMADA PELO PROGRAMA SOLOCILCOMP

double KA(double X1,double X2)
{
double eleva_KA (double b,double e);
    //y �a condutividade t�mica
	//X1 �O Kz OU Kr
	//X2 �O TETA
double a=0.467, //0.674;0.526; coeficiente do Psi na equa�o do potencial 
      m=0.934, //0.768;0.519; expoente na equa�o do potencial 
      n=2/(1-m), // expoente na equa�o do potencial, Burdines 
      y=0;

if (X2 <=0)
   y=0;
else

	if (X2 <=1){
      y=X1*(X2*X2)*eleva_KA (X2,m);//y=X1*(X2*X2)2*(1-(1-X2^(1/m))^m)^2; //Bordine
      }
	else
   	y=X1;
return (y);
}

double eleva_KA (double b,double e) {  // (1-(1-X2^(1/m))^m)^2;

     double C=0; 
        C=pow(b,1/e); 
        C=1-C;
        C=pow(C,e);
        C=1-C;
        C=pow(C,2);
   
  return (C);       
   }


double POT (double X)
{
double eleva_POT (double Y, double M, double N);
   //y �o potencial
   //X �o teta

double a=0.467, //0.674;0.526; coeficiente do Psi na equa�o do potencial */
      m=0.934, //0.768;0.519; expoente na equa�o do potencial */
      n=2/(1-m), // expoente na equa�o do potencial, Burdines */
      y=0,
      p=0;

if (X <=0)
   y=0.35;
else

	if (X <=1) {
     p=eleva_POT (X,m,n);     
   	 y=-1/a*p*0.098038;   //y=-1/a*(X^(-1/m)-1)^(1/n)*0.098038; pot matricial 
     }
	else
   	y=0;
return (y);
}

double eleva_POT (double Y, double M, double N) {  // (X^(-1/m)-1)^(1/n);

double C=0; 
        
        
        C=pow(Y,-1/M); 
        C=C-1;
        C=pow(C,1/N);
    
  return (C);       
   }

