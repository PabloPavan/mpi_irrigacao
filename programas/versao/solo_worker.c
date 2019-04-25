#include <stdio.h>
#include <math.h>
#include "pvm3.h"
#include "libMatematicas.h"
#include "KA.h"

#define NUMPROCESSOS 32

char buffer[100];

int mytid;       /* my task id */
int tids[NUMPROCESSOS];    /* task ids   */
int idpai;
int n, nworkers;
int tam, numlins, linnum=0, extra, Li, Lf ;

char arqlog[20]="solo_worker";

//Dados da malha

double *T, *TN;

int L=0,    // no. de pontos em r, varia com j
    M=0,     // no. de pontos em z, varia com i
    q=1,
    i=0,
    j=0,
    kk=0,
    nn=0,
    ini,
    ctr;


double PPE=0, PPN=0, PPP=0, PPS=0, PPW=0,
       AE=0, AW=0, AN=0, AS=0, AP=0,
       PW=0, PP=0, PS=0, PN=0, PE=0,
       dr=0, dz=0, dt=0,TT=0, x=0, Dr=0.15,
       Pi=3.1416,
       R=0.2, h=0.4, // raio do tubo,m
       VwE=0, VwN=0, VwS=0,
       MS=0, MP=0, ME=0, MN=0, MW=0,
       KE=0, KN=0, KW=0, KP=0, KS=0, K0=0.00003, // condutividade hid saturado
       KPr=0, KPz=0, Kfz=0,
       Kr=0, Kz=0, Kf=0,
       GE=0, GN=0, GP=0, GS=0, GW=0,
       VB=0,
       Tf=0,
       Tr=0.07,  // teor de umidade residual, dimensional
       Ts=0.76, // teor de umidade de satura�o, dimensional
       T0=0.1069;  // teor de umidade vol do solo em t=0 , dimensional */


void aloca() {
  int tam1=0;
    tam1=tam-1;

  // Aloca espaco para um vetor coluna
  if((T = (double *) malloc(tam*tam*sizeof(double)))==NULL) {
         puts("Nao existe memoria suficiente para este sistema.");
         exit(1);
  }

  if((TN = (double *) malloc (tam1*tam1*sizeof(double)))==NULL) {
         puts("Nao existe memoria suficiente para este sistema.");
         exit(1); 
  }
}


void recebe_parametrosIniciais () {
    mytid = pvm_mytid();

    /* Receive data from master */
        pvm_recv( -1, 1 );	// Recebe=1
	pvm_upkint(&nworkers, 1, 1);
	pvm_upkint(tids, nworkers, 1);
	pvm_upkint(&tam, 1, 1);
        pvm_upkint(&nn, 1 , 1);
        L=tam;
}

void recebe_dados() {
    //pvm_initsend(PvmDataDefault);
    pvm_recv(idpai, 2); 			// Recebe=2
    pvm_upkint(&Li,1 , 1);
    pvm_upkint(&Lf,1 , 1);
  
}

void recebe_TAtual() {
    //pvm_initsend(PvmDataDefault);
        pvm_recv(idpai , 3); 			// Recebe=3    - mudar
//     pvm_upkdouble(T, (tam*tam), 1); 
        pvm_upkint  (&ini,1 ,1 );
	pvm_upkint  (&ctr,1 ,1 );
 	pvm_upkdouble (&T[ini], ctr ,1);

}

void recebe_TNAtual() {
   // pvm_initsend(PvmDataDefault);
    pvm_recv(idpai, 4);		// Recebe=4
    pvm_upkdouble(TN, (tam-1)*(tam-1), 1); 

}



void envia_TNMestre(int ini) {

    pvm_initsend(PvmDataDefault);
    pvm_pkint (&ini,1,1);
    pvm_pkdouble(&TN[ini],L-2,1);
//     pvm_pkdouble(TN,(L-1)*(L-1),1);
    pvm_send(idpai,5);	// Envia=5
}

void envia_TMestre (int ind){
    pvm_initsend(PvmDataDefault);
    pvm_pkint (&ind,1,1);
    pvm_pkdouble(&T[ind],1,1);
//     pvm_pkdouble(T,tam*tam,1);
    pvm_send(idpai,6);	// Envia=6
}




void iniciaMatriz (){


for (i=0;i<(L-1)*(L-1);i++) {      //condi�o inicial para toda a malha
      TN[i]=0;
  } 
}


void calcula () {
int cont,cont1;
  // -----------SOLVE--------------- 


for (kk=1;kk<=nn;kk++) {

 recebe_TAtual ();

     for (i=Li;i<=(Lf-1);i++) {
 
//  recebe_TAtual ();

    	for (j=2;j<=(L-1);j++) {        // 3o. for, dire�o r

            // ----C�culo dos Potenciais ---------------------------
            // ----Potenciais p/ W
            GW=-9.8*dz*(i-1)*0.001*(T[L*(i-1)+(j-2)]*(Ts-Tr)+Tr); // Pot grav, em KPa
            
             if(T[L*(i-2)+(j-2)]>=1) { // Caso sat     
          		TT=T[L*(i-2)+(j-2)];  // in�io do c�culo do Pot Press�
          		q=1;
          
              while (TT >= 1) {
                 if(q >= (i-1)) {
                    TT=0.5;
                    } 
                 else {
                 	  TT=T[L*(i-2-q)+(j-2)];                    
                 	  q=q+1;
                      }
                    }  
                         	  
          	   PPW=q*dz*9.8*0.001; //Potencial de press� para W
                   MW=0; // Pot mat p/ sol sat
                   KW=K0;  // difusividade solo sat        
                  } 
    			else { //  Caso n� sat
       			  MW=-POT(T[L*(i-1)+(j-2)]);
                  KW=KA(Kr,T[L*(i-1)+(j-2)]);
                  PPW=0; // Pot press� p/ solo n� sat
                    }

 
				// ----Potenciais p/ P
           GP=-9.8*dz*(i-1)*0.001*(T[L*(i-1)+(j-1)]*(Ts-Tr)+Tr); // Pot grav, em KPa
           
             if(T[L*(i-2)+(j-1)]>=1) {  //  Caso sat
          		TT=T[L*(i-2)+(j-1)];  // in�io do c�culo do Pot Press�
          		q=1;
          		while (TT >= 1) {
                 if(q >= (i-1))  {
                    TT=0.5;
                    }
                 else {
                    TT=T[L*(i-2-q)+(j-1)];                    
                    q=q+1;
                    }
                   }
                          		    
                PPP=q*dz*9.8*0.001; // Potencial de press� para P
          		MP=0; //  Pot mat p/ sol sat
          		KP=K0;  // difusividade solo sat   
                    }               
    			else { //  Caso n� sat
    		       MP=-POT(T[L*(i-1)+(j-1)]);
                   KPr=KA(Kr,T[L*(i-1)+(j-1)]);
                   KPz=KA(Kz,T[L*(i-1)+(j-1)]);
                   PPP=0; // Pot press� p/ solo n� sat
               	}
             
				// ----Potenciais p/ E
           GE=-9.8*dz*(i-1)*0.001*(T[L*(i-1)+j]*(Ts-Tr)+Tr); // Pot grav, em KPa
            
          
            if(T[L*(i-2)+j]>=1) { //  Caso sat
          		TT=T[L*(i-2)+j];  // in�io do c�culo do Pot Press�
          		q=1;
          		while (TT >= 1){ 
                 if(q >= (i-1)) {
                    TT=0.5;
                    }
                 else {
                   TT=T[L*(i-2-q)+j];                    
                   q=q+1;
                   }
                  }
                PPE=q*dz*9.8*0.001;// Potencial de press� para P
          		ME=0; // Pot mat p/ sol sat
          		KE=K0;  // difusividade solo sat  
                  }             
    			else { // Caso n� sat
    			ME=-POT(T[L*(i-1)+j]);
                KE=KA(Kr,T[L*(i-1)+j]);
                PPE=0; // Pot press� p/ solo n� sat
       			
       		  }
                          
            // ----Potenciais p/ N
            GN=-9.8*dz*(i-1)*0.001*(T[L*(i)+(j-1)]*(Ts-Tr)+Tr); // Pot grav, em KPa
             
          		
            if (i==2)    {
              if(T[L*(i-2)+(j-1)]>=1) { //  Caso sat
          		PPN=dz*9.8*0.001; // Potencial de press� para P
          		MN=0; // Pot mat p/ sol sat
          		KN=K0;  // difusividade solo sat
                   }                  
    		  else { // Caso n� sat
                MN=-POT(T[L*(i)+(j-1)]);
                KN=KA(Kz,T[L*(i)+(j-1)]);
                PPN=0; // Pot press� p/ solo n� sat
                   }
                  }                                                   
              else {  // do 1o. if  
               
                if(T[L*(i-3)+(j-1)]>=1) { //  Caso sat
                  TT=T[L*(i-2)+(j-1)];  // in�io do c�culo do Pot Press�
          	      q=1;
                   
                   while (TT >= 1) {
                     TT=T[L*(i-2-q)+(j-1)];                    
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
       		 MN=-POT(T[L*(i-2)+(j-1)]);
                 KN=KA(Kz,T[L*(i-2)+(j-1)]);
                 PPN=0; // Pot press� p/ solo n� sat
                 }       
   			   } // do 1o. if      	
        
			// ----Potenciais p/ S
         GS=-9.8*dz*(i-1)*0.001*(T[L*(i)+(j-1)]*(Ts-Tr)+Tr); //Pot grav, em KPa

         if(T[L*(i-1)+(j-1)]>=1) {  // Caso sat
          		TT=T[L*(i-1)+(j-1)];  // in�io do c�culo do Pot Press�
          		q=1;
          		while (TT >= 1) {
                 TT=T[L*(i-q-1)+(j-1)];
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
       		    MS=-POT(T[L*(i-2)+(j-1)]);
                KS=KA(Kz,T[L*(i-2)+(j-1)]);
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

                TN[(L-1)*(i-1)+(j-1)] = (AE*PE+AW*PW+AP*PP+AN*PN+AS*PS)+T[L*(i-1)+(j-1)];
      		} // fim do 3o. for, j


          cont=((L-1)*(i-1))+1;
         envia_TNMestre (cont);
       //  T[L*(i-1)+(L-1)]=TN[(L-1)*(i-1)+(L-2)]; //  CONDI�O DE FRONTEIRA EM r=R, isolamento

 //       cont1=(L*(i-1))+(L-1); //Controlador para determinar tamanho da MATRIZ
//printf ("cont %d\n",cont1);
 //       envia_TMestre (cont1);
    }

  }
}


main (int argc, char *argv[]) {

   idpai = pvm_parent ();
   recebe_parametrosIniciais ();
   aloca ();
   iniciaMatriz();
   M=L;
   Kr=K0; //------------C�CULO DE TETA ---------------------
   Kz=K0; 
   Kf=K0;   //0.0001; em r
   Kfz=Kf;
   dr=R/(L-1);
   dz=h/(M-1);
   dt=2*Dr*(dr*dr)/K0; // segundos
   Pi=3.1416;
   recebe_dados ();
   calcula ();



pvm_exit ();
  
}


