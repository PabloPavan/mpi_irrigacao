#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "libMatematicas.h"
#include "pvm3.h"
#include "KA.h"

#define MAXHOSTS 32
#define NUMPROCESSOS 32


FILE *arq ;

int mytid, nproc, narch,procAtivos;
int info;
int tids[MAXHOSTS];
int nworkers=0;
struct pvmhostinfo hosts[MAXHOSTS];
char buffer[100];

int tam, Li, Lf;


double *T,*TN,*TT1;


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
    jj=0,
    ini=0,
    ctr=0,
    ind=0;

//------------DECLARA�O DOS TIPOS DA FUN�O
 double KA ( double X1, double X2);
 double POT (double X);

void aloca(){



if((T=(double *) malloc((L*L)*sizeof(double)))==NULL){
   puts("Nao existe memoria suficiente para este sistema.");
   exit(1);
	}

if((TN=(double *) malloc((L-1)*(L-1)*sizeof(double)))==NULL){
   puts("Nao existe memoria suficiente para este sistema.");
   exit(1);
	}

if((TT1=(double *) malloc(L*L*sizeof(double)))==NULL){
   puts("Nao existe memoria suficiente para este sistema.");
   exit(1);
	}
}

void iniciaMatriz (){


for (i=0;i<(L*L);i++) {      //condi�o inicial para toda a malha
      T[i]=T00;
   }
 
}

void dispara_processos(int numProcessos) {
    int resultado = 0;
    mytid = pvm_mytid();

    info = pvm_config(&nproc, &narch, (struct pvmhostinfo**)&hosts);
    // nworkers = nproc > 1 ? nproc-1 : 1;

    nworkers = numProcessos;
    resultado = pvm_spawn("solo_worker", (char**)0, PvmTaskDefault, "", nworkers, tids);
    if (resultado != nworkers) {
        printf("\nResultado de pvm_spawn=%d", resultado);
        exit(1);
 }
}

void envia_parametrosIniciais() {

int numlins, extra, linnum=0, k=2;  // distribui linhas aos clientes
    Li = k;
    numlins = (L-2) / (nworkers);
    extra = (L-2) % (nworkers);

    printf ("numero de linhas  %d  - extra %d \n", numlins, extra);

    pvm_initsend(PvmDataDefault);

    // Envia numero de trabalhadores, identificacao deles e o tam da matriz
    pvm_pkint(&nworkers, 1, 1);
    pvm_pkint(tids, nworkers, 1);
    pvm_pkint(&tam, 1, 1);
    pvm_pkint(&nn, 1 , 1);
    pvm_mcast(tids, nworkers, 1);	// Envio=1
 
    // Envia numero de linhas e a linha inicial; envia matriz T
    for(i=0; i < nworkers; i++) {
	pvm_initsend(PvmDataDefault);



	if (k == 2) {  // I inicial
		Lf = Li + numlins;
		if (extra > 0) {
			Lf ++;
			extra --; 
		}
		k = Lf;
       	}
       	else{
        	Li = Lf;
          	Lf = Lf + numlins;
		if (extra > 0) {
			Lf ++;
			extra --; 
		}
       	}



       pvm_pkint(&Li,1 , 1);
       pvm_pkint(&Lf,1 , 1);

       pvm_send(tids[i],2);		// Envio=2
       printf ("enviou para nodo %d - li = %d  - lf = %d - e = %d \n", tids[i], Li, Lf, extra);
    }

}

void envia_TAtual () {
    pvm_initsend(PvmDataDefault);
    pvm_pkdouble (T ,(tam*tam) ,1);
    pvm_mcast (tids,nworkers,3);	// Envio=3
}

void envia_TNAtual () {
    pvm_initsend(PvmDataDefault);
    pvm_pkdouble (TN ,(tam-1)*(tam-1) ,1);
    pvm_mcast (tids,nworkers,4);	// Envio=4
}



void recebe_TN(){
   int h,P,U;

    pvm_recv(-1,5);      // Recebe=5
    pvm_upkint (&ini,1,1);
    pvm_upkdouble(&TN[ini],L-2,1);
//     pvm_upkdouble(TN,(L-1)*(L-1),1);

}

void recebe_TWorker(){
    pvm_recv(-1,6);    // Recebe=5
    pvm_upkint (&ind,1,1);
    pvm_upkdouble(&T[ind],1,1);

//      pvm_upkdouble(T,tam*tam,1);
}


void calcula_mestre () {

//======================================================
//=====================================================
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
   VwsE=KA(Kf,T[0])*(T[1]-T[0])*2*Pi*dz*dt; // volume de �ua que sai/entra em E  
   VwsN=KA(Kfz,T[0])*(T[L]-T[0])*Pi*(dr*dr)*dt/dz;// volume de �ua que sai/entra em S
   VB=Vwe+VwsW+VwsE+VwsN; // balan� de �ua em cada VC
	T1=T[0]*(Ts-Tr)+Tr; // teta dimensional
   T[0]=(T1+VB/Vs-Tr)/(Ts-Tr); // teta adimensional	
   TT1[0]=T[0];

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
     
      VwsW=KA(Kf,T[jj-1])*(T[jj-2]-T[jj-1])*2*Pi*jj*dz*dt; // volume de �ua que sai/entra em W
      VwsE=KA(Kf,T[jj-1])*(T[jj]-T[jj-1])*2*Pi*(jj+1)*dz*dt; // volume de �ua que sai/entra em E
      VwsN=KA(Kfz,T[jj-1])*(T[L+(jj-1)]-T[jj-1])*Pi*(2*jj-1)*(dr*dr)*dt/dz; // volume de �ua que sai/entra em S
      VB=Vwe+VwsW+VwsE+VwsN; // balan� de �ua em cada VC
         
   Tf=T[jj-1]*(Ts-Tr)+Tr; // teta dimensional

  if ((jj-1)*dr >ri) {
      T[jj-1]=(T[jj-2]+T[jj]+T[L+(jj-1)])/3; //T[1][jj]
      }
   else {
      T[jj-1]=(Tf+VB/Vs-Tr)/(Ts-Tr); //  teta adimensional
      }    
      
   TT1[jj-1]=T[jj-1]; 
    }// do for com jj    

   T[L-1]=T[L-2];//  maquiagem fronteira r=R
   TT1[L-1]=TT1[L-2]; 

 envia_TAtual ();

//=================================================
//================================================
for (i=2;i<=(L-1);i++) { // 2o. for

         // ****C�culo de teta em e r=0***
          VwE=KA(Kf,T[L*(i-1)+1])*(T[L*(i-1)+1]-T[L*(i-1)])*2*Pi*dz*dt;  //volume de �ua que sai/entra em E
          VwS=KA(Kf,T[L*(i-2)])*(T[L*(i-2)])-T[L*(i-1)]*Pi*(dr*dr)*dt/dz; //volume de �ua que sai/entra em S
      	  VwN=KA(Kf,T[L*(i-1)])*(T[L*i]-T[L*(i-1)])*Pi*(dr*dr)*dt/dz; //volume de �ua que sai/entra em S
     	  VB=VwE+VwN+VwS; //balan� de �ua em cada VC
         
        Tf=T[L*(i-1)]*(Ts-Tr)+Tr; //teta dimensional
    //    T[L*(i-1)]=(T[L*(i-1)]+T[L*(i-2)]+T[L*i]+T[L*(i-1)+1])/4;   //(Tf+VB/Vs-Tr)/(Ts-Tr);  teta adimensional */


recebe_TN ();	

} 

//for (i=1;i<=nworkers;i++)



for (i=2;i<=(M-1);i++) {
		T[L*(i-1)+(L-1)]=TN[(L-1)*(i-1)+(L-2)];
		T[L*(i-1)]=(T[L*(i-1)]+T[L*(i-2)]+T[L*i]+T[L*(i-1)+1])/4;
	}


  for (i1=2;i1<=(M-1);i1++) {    // envelhecimento da TN
         for (j=2;j<=(L-1);j++) {
         	T[M*(i1-1)+(j-1)]=TN[(M-1)*(i1-1)+(j-1)];
          }
       }  
  
  for (j=2;j<=(L-1);j++){
		T[L*(M-1)+(j-1)]=TN[(L-1)*(M-2)+(j-1)]; //  maquiagem em Z=H
    }
   
      T[L*(M-1)+(L-1)]=TN[(L-1)*(M-2)+(L-2)];
      T[L*(M-1)]=T[L*(M-2)];
  }  //fim do 3o. 
}




main (int argc, char *argv[] ) {

   Kr = K0; //------------C�CULO DE TETA ---------------------
   Kz = K0; 
   Kf = K0;   //0.0001; em r
   Kfz = Kf;
   T00 = (T0-Tr)/(Ts-Tr);
   M = L;
   dr = R/(L-1);
   dz = h/(M-1);
   dt = 2*Dr*(dr*dr)/K0; // segundos
   tempo = dt*nn;
   n = 2/(1-m); // expoente na equa�o do potencial, Burdines */
   Vw = Vt*dt/ti; // volume (m3) de �ua que entra no solo por intervalo dt
   Pi = 3.1416;
   tam=L;


        if(argc<3) {
	   puts("digite: solo_main <arquivo de saida> numProcessos");
           exit(1);
        }

   aloca (); //Aloca memoria para o sistema
   iniciaMatriz (); //Define valores de saída da Matriz
   tempo1();
   dispara_processos (atoi(argv[2]));
   envia_parametrosIniciais (); //Envia dados para os Escravos
   calcula_mestre ();
   tempo2();
   tempoFinal("", argv[0], MSGLOG);

if ((arq = fopen(argv[1], "at")) == NULL) {
  	puts("Erro ao criar arquivo de saida");
  	exit(1);
  	}
   getchar();

arq=fopen (argv[1], "at");

for (i=0;i<((L)*(L));i++) {   //Escrevendo no arquivo
   fprintf (arq,"%.11f \n",T[i]);
  
  
  }
}







