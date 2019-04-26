#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "KA.h"
#include "libMatematicas.h"

// TRANSPORTE DE �UA NO SOLO EM COORDENDAS CIL�DRICAS - M�ODO EXPL�ITO


FILE *arq ;

//--------DEFINI�O DA MALHA ----------


 int L=102,    // no. de pontos em r, varia com j
     M=0,     // no. de pontos em z, varia com i
     q=1;

double n=0, 
      n1=0,
      n2=0,
      n3=0,
      m3=0,
      m5=4,    //constante para teste de converg�cia (depende do L)
      ME=0,
      R=0.2,   // raio do tubo,m
      h=0.4,   // altura do solo,m
      S=0,
      PE=0,
      Pi=0,//inverso do cosseno *
      AE=0,
      AW=0,
      AN=0,
      AS=0,
      AP=0,
      PN=0,
      PW=0,
      PP=0,
      PS=0,
      ri=0.06,
      dr=0,
      dz=0,
      tempo=0,
      dt=0,
      Dr=0.15,
      K0=0.00003, // condutividade hid saturado
      KE=0,
      KN=0,
      Kr=0,
      KS=0,
      Kz=0,
      Kf=0,
      Kfz=0,
      a=0.467, // 0.674; 0.526;   coeficiente do Psi na equa�o do potencial 
      m=0.934, // 0.768; 0.519;  expoente na equa�o do potencial
      PPE=0,
      PPN=0,
      PPP=0,
      PPS=0,
      MP=0,
      MS=0,
      KP=0,
      KPr=0,
      KPz=0,
      GE=0,
      GN=0,
      GP=0,
      GS=0,
      x=0;
    





//-----------CONDI�ES INICIAIS -----------

double  Tr=0.07,  // teor de umidade residual, dimensional  
       Ts=0.76, // teor de umidade de satura�o, dimensional
       T0=0.1069,  // teor de umidade vol do solo em t=0 , dimensional
       T00=0, // teor de umidade inicial adimensional
       T1=0,
       Tf=0,
       te=0;
     
//------------DADOS PARA O C�CULO DA IRRIGA�O ---------------------

double ti=2100,   // seg, tempo de irriga�o 
    GW=0,
    PPW=0,
    MN=0,
    MW=0,
    MWi=0,
    MW0=0,
    KW=0,
    TT=0,
    Vt=0.008,   // volume total de �ua de irriga�o, m3 
    VB=0, 
    Vs=0,
    Vw=0,
    Vwe=0,
    VwE=0,
    Vwn=0,
    VwN=0,
    Vws=0,
    VwS=0,
    VwsE=0,
    VwsN=0,
    VwsW=0;

//------------VARI�EIS PARA LA�S DE REPETI�O ---------------------    

int nn=100,  // itera�es temporais, varia com k
    i=0,
    i1=0,
    j=0,
    kk=0,
    jj=0;

//------------DECLARA�O DOS TIPOS DA FUN�O
 double KA ( double X1, double X2);
 double POT (double X);

main () {

   Kr=K0; //------------C�CULO DE TETA ---------------------
   Kz=K0; 
   Kf=K0;   //0.0001; em r
   Kfz=Kf;
   T00=(T0-Tr)/(Ts-Tr);
   M=L;
   dr=R/(L-1);
   dz=h/(M-1);
   dt=2*Dr*(dr*dr)/K0; // segundos
   tempo=dt*nn;
   n=2/(1-m); // expoente na equa�o do potencial, Burdines */
   Vw=Vt*dt/ti; // volume (m3) de �ua que entra no solo por intervalo dt
   Pi=3.1416;



double T[L][L], TT1[L],TN[L-1][L-1];

for (i=0;i<L;i++) {      //condi�o inicial para toda a malha
    TT1[i]=0;
   for (j=0;j<L;j++) {
      T[i][j]=T00;
     
   }
 }

for (i=0;i<L-1;i++) {      //condi�o inicial para toda a malha
   for (j=0;j<L-1;j++) {
      TN[i][j]=0;
     
   }
 } 


tempo1();

for (kk=1;kk<=nn;kk++) {  // 1o. for
   te=(kk-1)*dt;

           // -------CONDI�O DE FRONTEIRA EM Z=0 e Z=H------------------
        // --1o. VC--  
 Vs=Pi*(dr*dr)*dz; // no 1o. anel
   if (te <= ti){ // tempo de irriga�o    
      Vwe=Vt*dt*(dr*dr)/((ri*ri)*ti);
      }
   else {
      Vwe=0;
    }
   
   VwsW=0;
   VwsE=KA(Kf,T[0][0])*(T[0][1]-T[0][0])*2*Pi*dz*dt; // volume de �ua que sai/entra em E  
   VwsN=KA(Kfz,T[0][0])*(T[1][0]-T[0][0])*Pi*(dr*dr)*dt/dz;// volume de �ua que sai/entra em S
   VB=Vwe+VwsW+VwsE+VwsN; // balan� de �ua em cada VC
	T1=T[0][0]*(Ts-Tr)+Tr; // teta dimensional
   T[0][0]=(T1+VB/Vs-Tr)/(Ts-Tr); // teta adimensional	
   TT1[0]=T[0][0];

      //-------outros VC(1,jj) 
 for (jj=2;jj<=(L-1);jj++) {
      Vs=Pi*(2*jj-1)*(dr*dr)*dz;  // volume do anel de solo 
 
    if (te <= ti) { // tempo de irriga�o    
      	if ((jj-1)*dr <=ri) {  // zona de irriga�o     
         	Vwe=Vt*dt*(2*jj-1)*(dr*dr)/((ri*ri)*ti); // volume de �ua que entra, irriga�o, m3
          }
      	else {
         	Vwe=0;
             }
          }
      else {
        Vwe=0;
      }
     
      VwsW=KA(Kf,T[0][jj-1])*(T[0][jj-2]-T[0][jj-1])*2*Pi*jj*dz*dt; // volume de �ua que sai/entra em W
      VwsE=KA(Kf,T[0][jj-1])*(T[0][jj]-T[0][jj-1])*2*Pi*(jj+1)*dz*dt; // volume de �ua que sai/entra em E
      VwsN=KA(Kfz,T[0][jj-1])*(T[1][jj-1]-T[0][jj-1])*Pi*(2*jj-1)*(dr*dr)*dt/dz; // volume de �ua que sai/entra em S
      VB=Vwe+VwsW+VwsE+VwsN; // balan� de �ua em cada VC
         
   Tf=T[0][jj-1]*(Ts-Tr)+Tr; // teta dimensional

  if ((jj-1)*dr >ri) {
      T[0][jj-1]=(T[0][jj-2]+T[0][jj]+T[1][jj-1])/3; //T[1][jj]
      }
   else {
      T[0][jj-1]=(Tf+VB/Vs-Tr)/(Ts-Tr); //  teta adimensional
      }    
      
   TT1[jj-1]=T[0][jj-1]; 
    }// do for com jj    

   T[0][L-1]=T[0][L-2];//  maquiagem fronteira r=R
   TT1[L-1]=TT1[L-2]; 


// -----------SOLVE--------------- 
   for (i=2;i<=(M-1);i++) {   // 2o. for, dire�o Z 
        
           // ****C�culo de teta em e r=0***
          VwE=KA(Kf,T[i-1][1])*(T[i-1][1]-T[i-1][0])*2*Pi*dz*dt;  //volume de �ua que sai/entra em E
          VwS=KA(Kf,T[i-2][0])*(T[i-2][0])-T[i-1][0]*Pi*(dr*dr)*dt/dz; //volume de �ua que sai/entra em S
      	  VwN=KA(Kf,T[i-1][0])*(T[i][0]-T[i-1][0])*Pi*(dr*dr)*dt/dz; //volume de �ua que sai/entra em S
     	  VB=VwE+VwN+VwS; //balan� de �ua em cada VC
         
        Tf=T[i-1][0]*(Ts-Tr)+Tr; //teta dimensional
        

  	for (j=2;j<=(L-1);j++) {        // 3o. for, dire�o r
   	   
            // ----C�culo dos Potenciais ---------------------------
            // ----Potenciais p/ W
            GW=-9.8*dz*(i-1)*0.001*(T[i-1][j-2]*(Ts-Tr)+Tr); // Pot grav, em KPa

             if(T[i-2][j-2]>=1) { // Caso sat     
          		TT=T[i-2][j-2];  // in�io do c�culo do Pot Press�
          		q=1;
          
              while (TT >= 1) {
                 if(q >= (i-1)) {
                    TT=0.5;
                    } 
                 else {
                 	  TT=T[i-2-q][j-2];                    
                 	  q=q+1;
                      }
                    }  
                         	  
          	       PPW=q*dz*9.8*0.001; //Potencial de press� para W
                   MW=0; // Pot mat p/ sol sat
                   KW=K0;  // difusividade solo sat        
                  } 
    			else { //  Caso n� sat
       			  MW=-POT(T[i-1][j-2]);
                  KW=KA(Kr,T[i-1][j-2]);
                  PPW=0; // Pot press� p/ solo n� sat
                    }
 

				// ----Potenciais p/ P
           GP=-9.8*dz*(i-1)*0.001*(T[i-1][j-1]*(Ts-Tr)+Tr); // Pot grav, em KPa

            if(T[i-2][j-1]>=1) {  //  Caso sat
          		TT=T[i-2][j-1];  // in�io do c�culo do Pot Press�
          		q=1;
          		while (TT >= 1) {
                 if(q >= (i-1))  {
                    TT=0.5;
                    }
                 else {
                    TT=T[i-2-q][j-1];                    
                    q=q+1;
                    }
                   }
                          		    
                PPP=q*dz*9.8*0.001; // Potencial de press� para P
          		MP=0; //  Pot mat p/ sol sat
          		KP=K0;  // difusividade solo sat   
                    }               
    			else { //  Caso n� sat
    		       MP=-POT(T[i-1][j-1]);
                   KPr=KA(Kr,T[i-1][j-1]);
                   KPz=KA(Kz,T[i-1][j-1]);
                   PPP=0; // Pot press� p/ solo n� sat
               	}
             
				// ----Potenciais p/ E
           GE=-9.8*dz*(i-1)*0.001*(T[i-1][j]*(Ts-Tr)+Tr); // Pot grav, em KPa
            
          
            if(T[i-2][j]>=1) { //  Caso sat
          		TT=T[i-2][j];  // in�io do c�culo do Pot Press�
          		q=1;
          		while (TT >= 1){ 
                 if(q >= (i-1)) {
                    TT=0.5;
                    }
                 else {
                   TT=T[i-2-q][j];                    
                   q=q+1;
                   }
                  }
                PPE=q*dz*9.8*0.001;// Potencial de press� para P
          		ME=0; // Pot mat p/ sol sat
          		KE=K0;  // difusividade solo sat  
                  }             
    			else { // Caso n� sat
    			ME=-POT(T[i-1][j]);
                KE=KA(Kr,T[i-1][j]);
                PPE=0; // Pot press� p/ solo n� sat
       			
       		  }
                          
            // ----Potenciais p/ N
             GN=-9.8*dz*(i-1)*0.001*(T[i][j-1]*(Ts-Tr)+Tr); // Pot grav, em KPa
             
          		
            if (i==2)    {
              if(T[i-2][j-1]>=1) { //  Caso sat
          		PPN=dz*9.8*0.001; // Potencial de press� para P
          		MN=0; // Pot mat p/ sol sat
          		KN=K0;  // difusividade solo sat
                   }                  
    		  else { // Caso n� sat
                MN=-POT(T[i][j-1]);
                KN=KA(Kz,T[i][j-1]);
                PPN=0; // Pot press� p/ solo n� sat
                   }
                  }                                                   
              else {  // do 1o. if  
               
                if(T[i-3][j-1]>=1) { //  Caso sat
                  TT=T[i-2][j-1];  // in�io do c�culo do Pot Press�
          	      q=1;
                   
                   while (TT >= 1) {
                     TT=T[i-2-q][j-1];                    
                      q=q+1;
                  
                     if(q >= (i-1)) {
                         TT=0.5;
                          }
                        }
                      PPN=q*dz*9.8*0.001; // Potencial de press� para P
          	          MN=0; // Pot mat p/ sol sat
          		      KN=K0;  // difusividade solo sat
                       }    
                                                       
                  
          		else { // Caso n� sat
       			 MN=-POT(T[i-2][j-1]);
                 KN=KA(Kz,T[i-2][j-1]);
                 PPN=0; // Pot press� p/ solo n� sat
                 }       
   			   } // do 1o. if      	
        
			// ----Potenciais p/ S
         GS=-9.8*dz*(i-1)*0.001*(T[i][j-1]*(Ts-Tr)+Tr); //Pot grav, em KPa

         if(T[i-1][j-1]>=1) {  // Caso sat
          		TT=T[i-1][j-1];  // in�io do c�culo do Pot Press�
          		q=1;
          		while (TT >= 1) {
                 TT=T[i-q-1][j-1];                    
                 q=q+1;
                 if(q >= (i-1)) {
                    TT=0.5;
                      }
     	            }
                   PPS=q*dz*9.8*0.001; // Potencial de press� para P
                   MS=0; //  Pot mat p/ sol sat
                   KS=K0;  // difusividade solo sat
                   }   		              
    			else { //  Caso n� sat
       		    MS=-POT(T[i-2][j-1]);
                    KS=KA(Kz,T[i-2][j-1]);
                    PPS=0; // Pot press� p/ solo n� sat
       	 	  }


    			PW=MW+GW+PPW;
    			PP=MP+GP+PPP;
    			PE=ME+GE+PPE;
    			PN=MN+GN+PPN;
    			PS=MS+GS+PPS;
  
    			x=(j-1)*dr;  //x �o r
    			AE=KE*dt*(dr+2*x)/(2*x*(dr*dr));
    			AW=KW*dt*(2*x-dr)/(2*x*(dr*dr));
    			AN=KN*dt/(dz*dz);
    			AS=KS*dt/(dz*dz);
	 		AP=-2*dt*(KPr/(dr*dr)+KPz/(dz*dr));

   			TN[i-1][j-1]=(AE*PE+AW*PW+AP*PP+AN*PN+AS*PS)+T[i-1][j-1];

      		} // fim do 3o. for, j


	
}  //  fim do 2o. for, i 
	//novo 
	for (i=2;i<=(M-1);i++) {
       		T[i-1][L-1]=TN[i-1][L-2]; //  CONDI�O DE FRONTEIRA EM r=R, isolamento
		T[i-1][0]=(T[i-1][0]+T[i-2][0]+T[i][0]+T[i-1][1])/4;   //(Tf+VB/Vs-Tr)/(Ts-Tr);  teta adimensional	
	}
// int eu,ele;
// for (eu=0;eu<L;eu++){
// for (ele=0;ele<L;ele++){
// printf ("T %.11f\n",T[eu][ele]);
// }
// }

   for (i1=2;i1<=(M-1);i1++) {    // envelhecimento da TN
         for (j=2;j<=(L-1);j++) {
         	T[i1-1][j-1]=TN[i1-1][j-1];
          }
       }  
          
   for (j=2;j<=(L-1);j++){
		T[M-1][j-1]=TN[M-2][j-1]; //  maquiagem em Z=H
    }
   
      T[M-1][L-1]=TN[M-2][L-2];
      T[M-1][0]=T[M-2][0];
  }

tempo2();
tempoFinal("", "arquivo", MSGLOG);
getchar ();

 arq=fopen ("arquivo.txt", "at");   

for (i=0;i<L;i++) {   //Escrevendo no arquivo
   for (j=0;j<L;j++) {
    if (j==L-1){
                 fprintf (arq,"%.11f \n",T[i][j]);
               }
    else fprintf (arq,"%.11f \n",T[i][j]);
    }
    
   }     
  
 
fclose(arq); 
}

