//#include <cstdlib>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

using namespace std;

// TRANSPORTE DE ÁGUA NO SOLO EM COORDENDAS CILÍNDRICAS 

// --------variaveis da função evaporação---------
int tt[4][1];
double TEv[4][1];
double A[4][2];
double A2[2][2];
double mat1[2][4];
double mat2[2][2];
double mat3[2][4];
double x[2][1];

//--------DEFINIÇÃO DA MALHA ----------
int nn=100;  //iterações temporais, varia com k
int L=16;   // no. de pontos em r, varia com i
int M=16;   // no. de pontos em z, varia com j
double R=0.15;   // raio do tubo,m
double h=0.33;   //altura do solo,m
double dr=R/(L-1);
double dz=h/(M-1);
double dt=0.3; // horas

//-----------curva de retenção de água -----------
double a=0.467; //%0.674;%0.526;%;  % coeficiente do Psi na equação do potencial */
double m=0.934; //%0.768;%0.519;%;  % expoente na equação do potencial */
double n=2/(1-m); //expoente na equação do potencial, Burdines 
double Tr=0.07;  //teor de umidade residual, dimensional
double Ts=0.76; //teor de umidade de saturação, dimensional
double tempo=nn*dt;
double Ko=0.000004; //condutividade hidráulica de saturação

//----------DADOS EXPERIMENTAIS -----------
//---t = 0   CONDIÇOES INICIAIS -------

double T0=0.65;   //teor de umidade vol do solo em t=0 , dimensional
//double Tr=0.01; 
//double Ts=0.76; 
double T00=(T0-Tr)/(Ts-Tr); //teor de umidade inicial adimensional
double T[16][16];
double To[16][16];
double T1[16][16];
double T2[16][16];
double T3[16][16];
double Sf[16][16];
double TN[16][16];
double te [100][1];
double Erro[59][1];

int nc = 0,
    iot = 0;

double Menor = 0,
       Tsup = 0,
       SQE10 = 0,
       SQE20 = 0,
       SQE30 = 0,
       S1 = 0, 
       S2 = 0, 
       S3 = 0, 
       Med = 0, 
       SQT = 0, 
       SQT1 = 0, 
       SQT2 = 0, 
       SQT3 = 0, 
       R2 = 0,
       pow1 = 0, 
       pow2 = 0, 
       pow3 = 0,
       qo = 0,
       qM = 0,
       dq = 0,
       q = 0,
       qq = 0,
       qot = 0,
       q1 = 0,
       q2 = 0,
       q3 = 0,
       r = 0,
       KE = 0,
       KW = 0,
       KP = 0,
       KN = 0,
       KS = 0,
       AE = 0,
       AW = 0,
       AN = 0,
       AS = 0,
       APO = 0,
       PME = 0,
       PMW = 0,
       PMP = 0,
       PMN = 0,
       PMS = 0,
       ero = 0;

// FUNÇÃO PARA CALCULAR A CONDITIVIDADE HIDRÁULICA COM BASE NO TETA

double KA(double x1, double x2){    // X1 É O Ko    // X2 É O TETA
       double y=0;
       double a=0.467; // coeficiente do Psi na equação do potencial 
       double m=0.934; // expoente na equação do potencial 
       double n=2/(1-m); // expoente na equação do potencial, Burdines 
       if(x2 <= 0)
            y=0;
       else{
           if (x2 <=1){
   	           //Bordines
               pow1 = pow(1-x2,(1/m));          
               pow2 = pow(1 - pow1,m);
               pow3 = pow(0.5 * pow2, 2);
               y = pow(x1*x2, pow3);
           }
           else
   	           y=x1;
       }
       return(y);      
}

// FUNÇÃO PARA CALCULAR O POTENCIAL MATRICIAL COM BASE NO TETA

double POT(double x){  // X é o teta
       double y=0;
       double a=0.674; //coeficiente do Psi na equação do potencial 
       double m=0.768; // expoente na equação do potencial 
       double n=2/(1-m); // expoente na equação do potencial, Burdines 
       if (x <=0)
           y=0.35;
       else{      
           if (x <=1){
               pow1 = (pow(x,(-1/m))-1);
        	   pow2 = pow(((1/a)*pow1),1/n);
               y=pow(-10,pow2); // pot matricial 
           }
           else
   	           y=0;
       }
       return(y);
}

// FUNÇÃO PARA CÁLCULO DOS PARÂMETROS DA EVAPORAÇÃO---AJUSTE DE CURVAS EM Z=0  EVAPORAÇÃO-----

void evaporacao(){
    double valor=0;
    double linha1[4];
    double linha2[4];
    double laux[4];
    double soma=0;
    int nump=4; // número de pontos
    tt[1][1]=0;
    tt[2][1]=10;
    tt[3][1]=20;
    tt[4][1]=30;
    TEv[1][1]=To[1][2];
    TEv[2][1]=T1[1][2];
    TEv[3][1]=T2[1][2];
    TEv[4][1]=T3[1][2];
    for (int i=1; i <=nump; i++){
        A[i][1]=tt[i][1];
        A[i][2]=1;
    }
    // x=inv(A'*A)*A'*TEv A SEGUIR ESTA FORMULA
    
    for (int y0 = 1; y0 <=4; y0++){         //
         for (int z = 1; z <=2; z++){    // tranposta de A'
              mat1[z][y0] = A[y0][z];      //
         }
    }
    for (int y1 = 1; y1 <=4; y1++){                   //    
         soma = soma + (mat1[1][y1] * A[y1][1]);      //    
    }                                                 //
    A2[1][1] = soma;                                  //
    soma = 0;                                         //
    for (int y2 = 1; y2 <=4; y2++){                   //
         soma = soma + (mat1[1][y2] * A[y2][2]);      //
    }                                                 //
    A2[1][2] = soma;                                  // multiplicação de A'*A
    soma = 0;                                         //  A2 = (mat1 * A) 
    for (int y3 = 1; y3 <=4; y3++){                   //
         soma = soma + (mat1[2][y3] * A[y3][1]);      //
    }                                                 //
    A2[2][1] = soma;                                  //     
    soma = 0;                                         //
    for (int y4 = 1; y4 <=4; y4++){                   //
         soma = soma + (mat1[2][y4] * A[y4][2]);      //
    }                                                 //
    A2[2][2] = soma;  
    soma = 0;                                //
    // A SEGUIR MATRIZ INVERSA DE A2 E GRAVA EM MAT2
    linha1[0] = A2[1][1];
    linha1[1] = A2[1][2];
    linha1[2] = 1;
    linha1[3] = 0;          //definindo a malha
    linha2[0] = A2[2][1];
    linha2[1] = A2[2][2];
    linha2[2] = 0;
    linha2[3] = 1;
    //1º regra - o valor de linha1[0] deverá ser 1 - capturar o valor e dividir linha1 por ele
    valor = linha1[0];
    linha1[0] = linha1[0] / valor;
    linha1[1] = linha1[1] / valor;
    linha1[2] = linha1[2] / valor;
    linha1[3] = linha1[3] / valor;
    //2º regra - zerar o valor de linha2[0] - capturar o valor e multiplicar por (-) ele
    //toda a linha1 e somar com a linha2;
    valor = linha2[0];
    laux[0] = linha1[0] * (-valor);
    laux[1] = linha1[1] * (-valor);
    laux[2] = linha1[2] * (-valor);
    laux[3] = linha1[3] * (-valor);
    linha2[0] = laux[0] + linha2[0];
    linha2[1] = laux[1] + linha2[1];
    linha2[2] = laux[2] + linha2[2];
    linha2[3] = laux[3] + linha2[3];
    //3º regra - o valor de linha2[1] deverá ser 1 - capturar o valor e dividir linha2 por ele
    valor = linha2[1];
    linha2[0] = linha2[0] / valor;
    linha2[1] = linha2[1] / valor;
    linha2[2] = linha2[2] / valor;
    linha2[3] = linha2[3] / valor;
    //4º regra - zerar o valor de linha1[1] - capturar o valor e multiplicar por (-) ele
    //toda a linha2 e somar com a linha1;
    valor = linha1[1];
    laux[0] = linha2[0] * (-valor);
    laux[1] = linha2[1] * (-valor);
    laux[2] = linha2[2] * (-valor);
    laux[3] = linha2[3] * (-valor);
    linha1[0] = laux[0] + linha1[0];
    linha1[1] = laux[1] + linha1[1];
    linha1[2] = laux[2] + linha1[2];
    linha1[3] = laux[3] + linha1[3];
    //GRAVANDO EM MAT2
    mat2[1][1] = linha1[2];
    mat2[1][2] = linha1[3];
    mat2[2][1] = linha2[2];
    mat2[2][2] = linha2[3];
    
    //inv(A'*A)*A' ou seja mat2 * mat1
    for (int y5 = 1; y5 <=2; y5++){                   //    
         soma = soma + (mat2[1][y5] * mat1[y5][1]);   //    
    }                                                 //
    mat3[1][1] = soma;                                // 
    soma = 0;                                         //
    for (int y6 = 1; y6 <=2; y6++){                   //    
         soma = soma + (mat2[1][y6] * mat1[y6][2]);   //    
    }                                                 //
    mat3[1][2] = soma;                                // 
    soma = 0;                                         //
    for (int y7 = 1; y7 <=2; y7++){                   //    
         soma = soma + (mat2[1][y7] * mat1[y7][3]);   //    
    }                                                 //
    mat3[1][3] = soma;                                // 
    soma = 0;                                         //    
    for (int y8 = 1; y8 <=2; y8++){                   //    multiplicação
         soma = soma + (mat2[1][y8] * mat1[y8][4]);   //    mat3 =inv(A'*A)*A' 
    }                                                 //
    mat3[1][4] = soma;                                // 
    soma = 0;                                         //    
    for (int y9 = 1; y9 <=2; y9++){                   //    
         soma = soma + (mat2[2][y9] * mat1[y9][1]);   //    
    }                                                 //
    mat3[2][1] = soma;                                // 
    soma = 0;                                         //    
    for (int y10 = 1; y10 <=2; y10++){                //    
         soma = soma + (mat2[2][y10] * mat1[y10][2]); //    
    }                                                 //
    mat3[2][2] = soma;                                // 
    soma = 0;                                         //    
    for (int y11 = 1; y11 <=2; y11++){                //    
         soma = soma + (mat2[2][y11] * mat1[y11][3]); //    
    }                                                 //
    mat3[2][3] = soma;                                // 
    soma = 0;                                         //    
    for (int y12 = 1; y12 <=2; y12++){                //    
         soma = soma + (mat2[2][y12] * mat1[y12][4]); //    
    }                                                 //
    mat3[2][4] = soma;                                // 
    soma = 0;                                         //    

    //FINALIZANDO x = mat3 * TEv;
    for (int y13 = 1; y13 <=4; y13++){                //    
         soma = soma + (mat3[1][y13] * TEv[y13][1]);  //    
    }                                                 //
    x[1][1] = soma;                                   // 
    soma = 0;                                         //  x=inv(A'*A)*A'*TEv 
    for (int y14 = 1; y14 <=4; y14++){                //    
         soma = soma + (mat3[2][y14] * TEv[y14][1]);  //    
    }                                                 //
    x[2][1] = soma;                                   //

}

int main()
{
    for (int j=1; j<=M; j++){
       for (int i=1; i<=L; i++){
            T[i][j]=T00;  //condição inicial para toda a malha
       }
    } 

    //---t = 0 h -------
    To[1][2]=T00; To[1][5]=T00; To[1][8]=T00;
    To[3][2]=T00; To[3][5]=T00; To[3][8]=T00;
    To[6][2]=T00; To[6][5]=T00; To[6][8]=T00;
    To[11][2]=T00; To[11][5]=T00; To[11][8]=T00;

    //---t = 10 h -------
    T1[1][2]=0.75*T00; T1[1][5]=0.75*T00; T1[1][8]=0.75*T00;
    T1[3][2]=0.78*T00; T1[3][5]=0.8*T00; T1[3][8]=0.8*T00;
    T1[6][2]=0.7*T00; T1[6][5]=0.8*T00; T1[6][8]=0.9*T00;
    T1[11][2]=0.98*T00; T1[11][5]=0.85*T00; T1[11][8]=0.95*T00;

    //---t = 20 h -------
    T2[1][2]=0.5*T00; T2[1][5]=0.5*T00; T2[1][8]=0.5*T00;
    T2[3][2]=0.65*T00; T2[3][5]=0.7*T00; T2[3][8]=0.75*T00;
    T2[6][2]=0.45*T00; T2[6][5]=0.7*T00; T2[6][8]=0.8*T00;
    T2[11][2]=0.97*T00; T2[11][5]=0.8*T00; T2[11][8]=0.95*T00;

    //---t = 30 h -------
    T3[1][2]=0.3*T00; T3[1][5]=0.3*T00; T3[1][8]=0.35*T00;
    T3[3][2]=0.7*T00; T3[3][5]=0.6*T00; T3[3][8]=0.75*T00;
    T3[6][2]=0.6*T00; T3[6][5]=0.6*T00; T3[6][8]=0.75*T00;
    T3[11][2]=0.96*T00; T3[11][5]=0.75*T00; T3[11][8]=0.9*T00;

    evaporacao();
      
    //------------PROBLEMA INVERSO E DIRETO---------------------
    Menor=100000000;
    nc=59;  //número de q testados
    qo=0.035; //fluxo de transpiração mínimo (m3/hora)
    qM=0.09; //fluxo de transpiração máximo
    dq=(qM-qo)/nc;        




    for (int ii=1; ii<=nc; ii++){ 	// contador dos chutes do PI



       q=qo+(ii-1)*dq;   
       for (int j=1; j<=M; j++){
          	for (int i=1; i<=L; i++){
                 T[i][j]=T00;  //condição inicial para toda a malha 
          	}
      	}
	      
        for (int j=1; j<=L; j++){
	         for (int i=1; i<=M; i++) {
  			      Sf[i][j]=0;  //fonte zero para toda a malha 
	         }
        }
        
        for (int kk=1; kk<=nn; kk++){   //TEMPORAL
             te[kk][1]=(kk-1)*dt;
             qq=q*exp(pow((-0.05*(te[kk][1]-12)),2)); 
             for (int i=2; i<=(L-1); i++){
  	              Tsup=x[1][1]*te[kk][1]+x[2][1]; //teta em z=0 pela evaporaçao 
      	          for (int i1=1; i1<=L; i1++){
      		           T[1][i1]=Tsup; //condição de fronteira em z=0
      	      	  }    	      	   	    	      	        	  
 	         	  for (int j=2; j<=(M-1); j++){   
                        
              		   //------------FONTE = TRANSPIRAÇÃO ----------------
			           q1=qq;
			           q2=qq;
			           q3=qq;
			
   			           Sf[2][2]= qq*(T[2][2]-Tr); Sf[3][2]= qq*(T[3][2]-Tr);
   			           Sf[4][2]= qq*(T[4][2]-Tr); Sf[4][3]=q1*(T[4][3]-Tr);
   			           Sf[5][2]= qq*(T[5][2]-Tr); Sf[5][3]=q1*(T[5][3]-Tr); Sf[5][4]=q2*(T[5][4]-Tr);
   			           Sf[6][2]= qq*(T[6][2]-Tr); Sf[6][3]=q1*(T[6][3]-Tr); Sf[6][4]=q2*(T[6][4]-Tr); Sf[6][5]=q3*(T[6][5]-Tr);
   			           Sf[7][2]= qq*(T[7][2]-Tr); Sf[7][3]=q1*(T[7][3]-Tr); Sf[7][4]=q2*(T[7][4]-Tr);
   			           Sf[8][2]= qq*(T[8][2]-Tr); Sf[8][3]=q1*(T[8][3]-Tr);
   			           Sf[9][2]= qq*(T[9][2]-Tr); Sf[10][2]= qq*(T[10][2]-Tr);

                   	   r=(j-1)*dr + 0.5 * dr;
                             
                       KE= KA(Ko,T[i][j+1]);
                       KW= KA(Ko,T[i][j-1]);
                       KP= KA(Ko,T[i][j]);
                       KN= KA(Ko,T[i+1][j]);
			 	       KS= KA(Ko,T[i-1][j]);
 
                       AE= dt/dr*((KE-KW)/(4*dr)+KP/(2*r)+KP/dr);
                       AW= dt/dr*(-(KE-KW)/(4*dr)-KP/(2*r)+KP/dr);
        		       AN= dt/pow(4*dz,2)*(KN-KS+4*KP);  
         		       AS= dt/pow(4*dz,2)*(-KN+KS+4*KP); 
        		       APO= -2*dt*KP*(pow(1/dr,2)+pow(1/dz,2));
               
                       PME=POT(T[i][j+1]);
                       PMW=POT(T[i][j-1]);
                       PMP=POT(T[i][j]);
                       PMN=POT(T[i+1][j]);
                       PMS=POT(T[i-1][j]);
                       TN[i][j]=(AE*PME+AW*PMW+AN*PMN+AS*PMS+APO*PMP)+(KN-KS) / (2*dz)+T[i][j]-Sf[i][j];
                         
                  } //end do for do j

		          TN[i][1]= TN[i][2]; // CONDIÇÃO DE FRONTEIRA EM r=R, isolamento
                  TN[i][L]= TN[i][L-1]; // CONDIÇÃO DE FRONTEIRA EM r=R, isolamento

             }//end do for do i
     
             for (int j1=1; j1<=L; j1++){ // condição de fronteira em z=h, isolamento
                  TN[M][j1]= TN[M-1][j1];
             }
             for (int i1=2; i1<=L; i1++){    //envelhecimento da TN
                  for (int j1=1; j1<=M; j1++){
         	           T[i1][j1]= TN[i1][j1];
                  }
             }

             //-------PROBLEMA INVERSO-----------

             if (fabs(te[kk][1]-10) <= dt){ //10 horas     
                 SQE10=pow((T1[1][2]-T[1][2]),2)+pow((T1[1][5]-T[1][5]),2)+pow((T1[1][8]-T[1][8]),2)+pow((T1[3][2]-T[3][2]),2)+pow((T1[3][5]-T[3][5]),2)+pow((T1[3][8]-T[3][8]),2)+pow((T1[6][2]-T[6][2]),2)+pow((T1[6][5]-T[6][5]),2)+pow((T1[6][8]-T[6][8]),2)+pow((T1[11][2]-T[11][2]),2)+pow((T1[11][5]-T[11][5]),2)+pow((T1[11][8]-T[11][8]),2);    
             }
             if (fabs(te[kk][1]-20) <= dt){ //20 horas
       	         SQE20=pow((T2[1][2]-T[1][2]),2)+pow((T2[1][5]-T[1][5]),2)+pow((T2[1][8]-T[1][8]),2)+pow((T2[3][2]-T[3][2]),2)+pow((T2[3][5]-T[3][5]),2)+pow((T2[3][8]-T[3][8]),2)+pow((T2[6][2]-T[6][2]),2)+pow((T2[6][5]-T[6][5]),2)+pow((T2[6][8]-T[6][8]),2)+pow((T2[11][2]-T[11][2]),2)+pow((T2[11][5]-T[11][5]),2)+pow((T2[11][8]-T[11][8]),2);  
             }
             if (fabs(te[kk][1]-29) <= dt){ //30 horas
     		     SQE30=pow((T3[1][2]-T[1][2]),2)+pow((T3[1][5]-T[1][5]),2)+pow((T3[1][8]-T[1][8]),2)+pow((T3[3][2]-T[3][2]),2)+pow((T3[3][5]-T[3][5]),2)+pow((T3[3][8]-T[3][8]),2)+pow((T3[6][2]-T[6][2]),2)+pow((T3[6][5]-T[6][5]),2)+pow((T3[6][8]-T[6][8]),2)+pow((T3[11][2]-T[11][2]),2)+pow((T3[11][5]-T[11][5]),2)+pow((T3[11][8]-T[11][8]),2);
             }
        } //end do for do kk

 cout << ii;

	Erro[ii][1]=SQE10+SQE20+SQE30;

        if(Erro[ii][1] <= Menor){       
	       qot=q;
   	       iot=ii;
   	       Menor= Erro[ii][1];
	    }
        else{
   	       Menor= Menor;
        }
 cout << ii;

    } //end do ii, Problema inverso

    //---escolha do menor Erro---

    //printf(" %f",qot);
    //printf(" %f",iot);
         
    S1=(T1[1][2]+T1[1][5]+T1[1][8]+T1[3][2]+T1[3][5]+T1[3][8]+T1[6][2]+T1[6][5]+T1[6][8]+T1[11][2]+T1[11][5]+T1[11][8])/12;    
    S2=(T2[1][2]+T2[1][5]+T2[1][8]+T2[3][2]+T2[3][5]+T2[3][8]+T2[6][2]+T2[6][5]+T2[6][8]+T2[11][2]+T2[11][5]+T2[11][8])/12;    
    S3=(T3[1][2]+T3[1][5]+T3[1][8]+T3[3][2]+T3[3][5]+T3[3][8]+T3[6][2]+T3[6][5]+T3[6][8]+T3[11][2]+T3[11][5]+T3[11][8])/12;    
    Med=(S1+S2+S3)/3;
    SQT1=pow((T1[1][2]-Med),2)+pow((T1[1][5]-Med),2)+pow((T1[1][8]-Med),2)+pow((T1[3][2]-Med),2)+pow((T1[3][5]-Med),2)+pow((T1[3][8]-Med),2)+pow((T1[6][2]-Med),2)+pow((T1[6][5]-Med),2)+pow((T1[6][8]-Med),2)+pow((T1[11][2]-Med),2)+pow((T1[11][5]-Med),2)+pow((T1[11][8]-Med),2);    
    SQT2=pow((T2[1][2]-Med),2)+pow((T2[1][5]-Med),2)+pow((T2[1][8]-Med),2)+pow((T2[3][2]-Med),2)+pow((T2[3][5]-Med),2)+pow((T2[3][8]-Med),2)+pow((T2[6][2]-Med),2)+pow((T2[6][5]-Med),2)+pow((T2[6][8]-Med),2)+pow((T2[11][2]-Med),2)+pow((T2[11][5]-Med),2)+pow((T2[11][8]-Med),2);  
    SQT3=pow((T3[1][2]-Med),2)+pow((T3[1][5]-Med),2)+pow((T3[1][8]-Med),2)+pow((T3[3][2]-Med),2)+pow((T3[3][5]-Med),2)+pow((T3[3][8]-Med),2)+pow((T3[6][2]-Med),2)+pow((T3[6][5]-Med),2)+pow((T3[6][8]-Med),2)+pow((T3[11][2]-Med),2)+pow((T3[11][5]-Med),2)+pow((T3[11][8]-Med),2);
    SQT=pow((SQT1+SQT2+SQT3),2);

    cout << iot;

    ero = Erro[iot][1];
    R2=(1-ero)/SQT;

    //=====SOLUÇÃO ÓTIMA======

    q=qot; 
   
    //------------FONTE = TRANSPIRAÇÃO ----------------
	q1=q;
	q2=q;
	q3=q;
	for (int j=1; j<=L; j++){
	   	for (int i=1; i<=M; i++){
             Sf[i][j]=0;  //condição inicial para toda a malha
        }
	} 
	Sf[2][2]=q;
	Sf[3][2]=q;
	Sf[4][2]=q;Sf[4][3]=q1;
	Sf[5][2]=q;Sf[5][3]=q1;Sf[5][4]=q2;
	Sf[6][2]=q;Sf[6][3]=q1;Sf[6][4]=q2;Sf[6][5]=q3;
	Sf[7][2]=q;Sf[7][3]=q1;Sf[7][4]=q2;
	Sf[8][2]=q;Sf[8][3]=q1;
	Sf[9][2]=q;
	Sf[10][2]=q;
   
    for (int j=1; j<=M; j++){
   	    for(int i=1; i<=L; i++){
       	     T[i][j]=T00;  //condição inicial para toda a malha
   	    }
	}

   for (int kk=1; kk<=nn; kk++){   // TEMPORAL
        te[kk][1]=(kk-1)*dt;
        qq=q*exp(pow((-0.05*(te[kk][1]-12)),2));

   	    for (int i=2; i<=(L-1); i++){
    	     Tsup=x[1][1]*te[kk][1]+x[2][1]; //teta em z=0 pela evaporaçao 
 	         for (int i1=1; i1<=L; i1++){
      		      T[1][i1]=Tsup; //condição de fronteira em z=0
      	     }
             for (int j=2; j<=(M-1); j++){    
                      
                  //------------FONTE = TRANSPIRAÇÃO ----------------
			      q1=qq;
			      q2=qq;
			      q3=qq;
			
   			      Sf[2][2]=qq*(T[2][2]-Tr);
			      Sf[3][2]=qq*(T[3][2]-Tr);
			      Sf[4][2]=qq*(T[4][2]-Tr);Sf[4][3]=q1*(T[4][3]-Tr);
			      Sf[5][2]=qq*(T[5][2]-Tr);Sf[5][3]=q1*(T[5][3]-Tr);Sf[5][4]=q2*(T[5][4]-Tr);
			      Sf[6][2]=qq*(T[6][2]-Tr);Sf[6][3]=q1*(T[6][3]-Tr);Sf[6][4]=q2*(T[6][4]-Tr);Sf[6][5]=q3*(T[6][5]-Tr);
			      Sf[7][2]=qq*(T[7][2]-Tr);Sf[7][3]=q1*(T[7][3]-Tr);Sf[7][4]=q2*(T[7][4]-Tr);
			      Sf[8][2]=qq*(T[8][2]-Tr);Sf[8][3]=q1*(T[8][3]-Tr);
			      Sf[9][2]=qq*(T[9][2]-Tr);
			      Sf[10][2]=qq*(T[10][2]-Tr);

         	      r=(j-1)*dr+0.5*dr;
                                               
                  KE=KA(Ko,T[i][j+1]);
                  KW=KA(Ko,T[i][j-1]);
                  KP=KA(Ko,T[i][j]);
                  KN=KA(Ko,T[i+1][j]);
			 	  KS=KA(Ko,T[i-1][j]);
 
                  AE=dt/dr*((KE-KW)/(4*dr)+KP/(2*r)+KP/dr);
        		  AW=dt/dr*(-(KE-KW)/(4*dr)-KP/(2*r)+KP/dr);
        		  AN=dt/pow(4*dz,2)*(KN-KS+4*KP);
        		  AS=dt/pow(4*dz,2)*(-KN+KS+4*KP);
        		  APO=-2*dt*KP*(pow(1/dr,2)+pow(1/dz,2));
                                           
                  PME=POT(T[i][j+1]);
                  PMW=POT(T[i][j-1]);
                  PMP=POT(T[i][j]);
                  PMN=POT(T[i+1][j]);
                  PMS=POT(T[i-1][j]);
            
                  TN[i][j]=(AE*PME+AW*PMW+AN*PMN+AS*PMS+APO*PMP)+(KN-KS)/(2*dz)+T[i][j]-Sf[i][j];
         
              } // end do for do j
		  
              TN[i][1]=TN[i][2]; //CONDIÇÃO DE FRONTEIRA EM r=R, isolamento*/
              TN[i][L]=TN[i][L-1]; // CONDIÇÃO DE FRONTEIRA EM r=R, isolamento*/

           } //end do for do i
     
           for (int j1=1; j1<=L; j1++){ //condição de fronteira em z=h, isolamento
                TN[M][j1]=TN[M-1][j1];
           }
           for (int i1=2; i1<=L; i1++){    //envelhecimento da TN
                for (int j1=1; j1<=M; j1++){
         	         T[i1][j1]=TN[i1][j1];
                }
           }
	}  // end do for do kk
    //system("PAUSE");
  
}
