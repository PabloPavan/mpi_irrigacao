#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <omp.h>

#define  L 448   // no. de pontos em r, varia com i
#define  M L  // no. de pontos em z, varia com j
	
	// TRANSPORTE DE AGUA NO SOLO EM COORDENDAS CILDRICAS 
	
	// --------variaveis da funo evaporao---------
	int tt[4][1];
	double TEv[4][1];
	double A[4][2];
	double A2[2][2];
	double mataux[2][4];
	double mat1[2][4];
	double mat2[2][2];
	double mat3[2][4];
	double x[2][1];
	
	//--------DEFINICAO DA MALHA ----------
	int nn=100;  //iteraes temporais, varia com k
	int nc = 3072;  //nmero de q testados
	double R=0.015;   // raio do tubo,m
	double h=0.033;   //altura do solo,m
	double dr;
	double dz;
	double dt=0.3; // horas
	double T0 = 0.65;  
	double Tr=0.01;
	double Ts=0.76; //teor de umidade de saturao, dimensional
	double T00;	

	
	//-------x----curva de reteno de ua -----------
	double a=0.467; //%0.674;%0.526;%;  % coeficiente do Psi na equao do potencial */
	double m=0.934; //%0.768;%0.519;%;  % expoente na equao do potencial */
	double n;	
	double tempo;	
	double Ko=0.000004; //condutividade hidrlica de saturao
	double T[L][M];
	double T1[L][M];
	double T2[L][M];
	double T3[L][M];
	double Sf[L][M];
	double TN[L][M];
	double te[100][1];
	double erro[60];
	double vqot[60];
	double Menor = 100000000;
	double qo = 0.035; //fluxo de transpirao mimo (m3/hora)
	double qM = 0.09; //fluxo de transpirao mimo
	double dq;
	double startwtime = 0.0, endwtime;
		
	int kk = 0,
	i = 0,
	j = 0,
	iot = 0,
	tamanho = 0;
	
	double ero = 0,
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
	menor_erro = 0;


//prototipação das funcoes
	double KA(double x1, double x2);
	double POT(double x);
	double evaporacao();
	int worker();

int main(int argc, char *argv[])
{
	double start = omp_get_wtime();
			
	dr=R/(L-1);
	dz=h/(M-1);
	T00=(T0-Tr)/(Ts-Tr); //teor de umidade inicial adimensional
	tempo=nn*dt;
	dq=(qM-qo)/59; 
	n=2/(1-m); //expoente na equao do potencial, Burdines 
	int iv,j,i,i1,j1;
	int nProc, proc;
	int posi = 0;
	MPI_Status status;
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&nProc);
//	printf(" \n\n Processos criados %d", nProc); 	

	MPI_Comm_rank(MPI_COMM_WORLD,&proc);
//	printf(" \n\n Eu suo o processo %d", proc); 	



	if (proc == 0) /*mestre*/
	{




		evaporacao();
		double evap1 = x[1][1];
		double evap2 = x[2][1];
		int pos_ini = 1;
		int pos_fim = 0;
		tamanho = nc / (nProc-1); //tamanho do intervalo de chutes (chutes / processos)
		double vetor_erro[tamanho];
		double vetor_qot[tamanho];
		int patual,ratual,zera,posm = 0, posi = 0;
		
		
		
		//Envia parametros (evaporacao, posicoes inicial e final)
		for (patual = 1; patual <= (nProc-1); patual++){
//			printf(" \n\n enviou para comp %d", patual); 
//			startwtime = MPI_Wtime();
//			printf("\n start clock time = %f\n",startwtime);	
			pos_fim = tamanho + pos_fim;
			MPI_Send(&tamanho, 1, MPI_INT, patual, 0, MPI_COMM_WORLD);
			MPI_Send(&evap1, 1, MPI_DOUBLE, patual, 0, MPI_COMM_WORLD);
			MPI_Send(&evap2, 1, MPI_DOUBLE, patual, 0, MPI_COMM_WORLD);
			MPI_Send(&pos_ini, 1, MPI_INT, patual, 0, MPI_COMM_WORLD);
			MPI_Send(&pos_fim, 1, MPI_INT, patual, 0, MPI_COMM_WORLD);
			pos_ini = pos_fim + 1;
			
	
		}
		
		pos_ini = 0;
		pos_fim = 0;
		for (ratual = 1; ratual <=(nProc-1); ratual++){
            MPI_Recv(vetor_erro, tamanho, MPI_DOUBLE, ratual, 10, MPI_COMM_WORLD, &status);
           	MPI_Recv(vetor_qot, tamanho, MPI_DOUBLE, ratual, 10, MPI_COMM_WORLD, &status);
//			printf(" \n\n recebeu de %d", ratual);
//			endwtime = MPI_Wtime();
//			printf("\nwall clock time = %f\n",endwtime);
			pos_fim = (tamanho) + pos_fim;
			for (posm = pos_ini; posm < pos_fim; posm++){
				erro[posm] = vetor_erro[posi];
				vqot[posm] = vetor_qot[posi];
//				printf(" \n\n iot= %d", posm);
//				printf(" \n\n qot= %f", erro[1][posm]);
				posi++;
			}
			posi = 0;
			pos_ini = pos_fim;
		}


		//tempo2();
		//tempoFinal("mili segundos", argv[0], MSGLOG);


	} else {		
		worker();
	}




//////////////////////////////////////////////////////////	
	if (proc == 0) /*mestre*/
	{
		qot=0;
		for (iv = 0; iv <= 59; iv++){
			menor_erro = erro[iv];      
			if(menor_erro < Menor){       
			qot = vqot[iv];
			Menor = menor_erro;
			ero = menor_erro;
			}
			else{
				Menor = Menor;
			}
		}

		printf(" \n\n Erro: %f \n\n", qot);

		//---escolha do menor Erro---
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
		
		S1=(T1[1][2]+T1[1][5]+T1[1][8]+T1[3][2]+T1[3][5]+T1[3][8]+T1[6][2]+T1[6][5]+T1[6][8]+T1[11][2]
		+T1[11][5]+T1[11][8])/12;
    		S2=(T2[1][2]+T2[1][5]+T2[1][8]+T2[3][2]+T2[3][5]+T2[3][8]+T2[6][2]+T2[6][5]+T2[6][8]+T2[11][2]
		+T2[11][5]+T2[11][8])/12;
    		S3=(T3[1][2]+T3[1][5]+T3[1][8]+T3[3][2]+T3[3][5]+T3[3][8]+T3[6][2]+T3[6][5]+T3[6][8]+T3[11][2]
		+T3[11][5]+T3[11][8])/12;
		Med=(S1+S2+S3)/3;
	
		SQT1=(T1[1][2]-Med)*(T1[1][2]-Med)+(T1[1][5]-Med)*(T1[1][5]-Med)+(T1[1][8]-Med)*
		(T1[1][8]-Med)+(T1[3][2]-Med)*(T1[3][2]-Med)+(T1[3][5]-Med)*(T1[3][5]-Med)+(T1[3][8]-Med)
		*(T1[3][8]-Med)+(T1[6][2]-Med)*(T1[6][2]-Med)+(T1[6][5]-Med)*(T1[6][5]-Med)+(T1[6][8]-Med)
		*(T1[6][8]-Med)+(T1[11][2]-Med)*(T1[11][2]-Med)+(T1[11][5]-Med)*(T1[11][5]-Med)+
		(T1[11][8]-Med)*(T1[11][8]-Med);
	
		SQT2=(T2[1][2]-Med)*(T2[1][2]-Med)+(T2[1][5]-Med)*(T2[1][5]-Med)+(T2[1][8]-Med)
		*(T2[1][8]-Med)+(T2[3][2]-Med)*	(T2[3][2]-Med)+(T2[3][5]-Med)*(T2[3][5]-Med)+(T2[3][8]-Med)*
		(T2[3][8]-Med)+(T2[6][2]-Med)*(T2[6][2]-Med)+(T2[6][5]-Med)*(T2[6][5]-Med)+(T2[6][8]-Med)
		*(T2[6][8]-Med)+(T2[11][2]-Med)*(T2[11][2]-Med)+(T2[11][5]-Med)*(T2[11][5]-Med)+
		(T2[11][8]-Med)*(T2[11][8]-Med);
	
		SQT3=(T3[1][2]-Med)*(T3[1][2]-Med)+(T3[1][5]-Med)*(T3[1][5]-Med)+(T3[1][8]-Med)
		*(T3[1][8]-Med)+(T3[3][2]-Med)*(T3[3][2]-Med)+(T3[3][5]-Med)*(T3[3][5]-Med)+(T3[3][8]-Med)
		*(T3[3][8]-Med)+(T3[6][2]-Med)*(T3[6][2]-Med)+(T3[6][5]-Med)*(T3[6][5]-Med)+(T3[6][8]-Med)
		*(T3[6][8]-Med)+(T3[11][2]-Med)*(T3[11][2]-Med)+(T3[11][5]-Med)*(T3[11][5]-Med)+
		(T3[11][8]-Med)*(T3[11][8]-Med);
	
		SQT=(SQT1+SQT2+SQT3)*(SQT1+SQT2+SQT3);
		R2=1-(ero/SQT);
		//printf ("\n %ld", R2); 		//cout << "\n" << R2;
		
		//=====SOLUO IMA======
		q=qot; 
		//------------FONTE = TRANSPIRAO ----------------
		q1=q;
		q2=q;
		q3=q;
		for (j=1; j<=L; j++){
			for (i=1; i<=M; i++){
				Sf[i][j]=0;  //condio inicial para toda a malha
			}
		} 
		Sf[2][2]=q;Sf[3][2]=q;Sf[4][2]=q;Sf[4][3]=q;Sf[5][2]=q;Sf[5][3]=q;Sf[5][4]=q;Sf[6][2]=q;
		Sf[6][3]=q;Sf[6][4]=q;Sf[6][5]=q;Sf[7][2]=q;Sf[7][3]=q;Sf[7][4]=q;Sf[8][2]=q;
		Sf[8][3]=q;Sf[9][2]=q;Sf[10][2]=q;
	
		for (j=1; j<=M; j++){
			for(i=1; i<=L; i++){
				T[i][j]=T00;  //condio inicial para toda a malha
			}
		}
	
		for (kk=1; kk<=nn; kk++){   // TEMPORAL      
			te[kk][1]=(kk-1)*dt;
			qq=q*exp(-0.05*((te[kk][1]-12)*(te[kk][1]-12))); 
			for (i=2; i<=(L-1); i++){
			Tsup=x[1][1]*te[kk][1]+x[2][1]; //teta em z=0 pela evaporao 
				for (i1=1; i1<=L; i1++){
					T[1][i1]=Tsup; //condio de fronteira em z=0
				}
				for (j=2; j<=(M-1); j++){    
					//------------FONTE = TRANSPIRAO ----------------
					q1=qq;
					q2=qq;
					q3=qq;
				Sf[2][2]=qq*(T[2][2]-Tr);Sf[3][2]=qq*(T[3][2]-Tr);Sf[4][2]=qq*(T[4][2]-Tr);
				Sf[4][3]=q1*(T[4][3]-Tr);Sf[5][2]=qq*(T[5][2]-Tr);Sf[5][3]=q1*(T[5][3]-Tr);
				Sf[5][4]=q2*(T[5][4]-Tr);Sf[6][2]=qq*(T[6][2]-Tr);Sf[6][3]=q1*(T[6][3]-Tr);
				Sf[6][4]=q2*(T[6][4]-Tr);Sf[6][5]=q3*(T[6][5]-Tr);Sf[7][2]=qq*(T[7][2]-Tr);
				Sf[7][3]=q1*(T[7][3]-Tr);Sf[7][4]=q2*(T[7][4]-Tr);Sf[8][2]=qq*(T[8][2]-Tr);
				Sf[8][3]=q1*(T[8][3]-Tr);Sf[9][2]=qq*(T[9][2]-Tr);Sf[10][2]=qq*(T[10][2]-Tr);
					r=(j-1)*dr+0.5*dr;                                               
					KE=KA(Ko,T[i][j+1]);
					KW=KA(Ko,T[i][j-1]);
					KP=KA(Ko,T[i][j]);
					KN=KA(Ko,T[i+1][j]);
					KS=KA(Ko,T[i-1][j]);
					AE=dt/dr*((KE-KW)/(4*dr)+KP/(2*r)+KP/dr);
					AW=dt/dr*(-(KE-KW)/(4*dr)-KP/(2*r)+KP/dr);
					AN= dt/(4*(dz*dz))*(KN-KS+4*KP);  
					AS= dt/(4*(dz*dz))*(-KN+KS+4*KP); 
					APO= -(2*dt*KP*(((1/dr)*(1/dr))+((1/dz)*(1/dz))));  
					PME=POT(T[i][j+1]);
					PMW=POT(T[i][j-1]);
					PMP=POT(T[i][j]);
					PMN=POT(T[i+1][j]);
					PMS=POT(T[i-1][j]);
					TN[i][j]=(AE*PME+AW*PMW+AN*PMN+AS*PMS+APO*PMP)+(KN-KS)
					/(2*dz)+T[i][j]-Sf[i][j];
				} // end do for do j
			
				TN[i][1]=TN[i][2]; //CONDIO DE FRONTEIRA EM r=R, isolamento
				TN[i][L]=TN[i][L-1]; // CONDIO DE FRONTEIRA EM r=R, isolamento
			} //end do for do i

			for (j1=1; j1<=L; j1++){ //condio de fronteira em z=h, isolamento
				TN[M][j1]=TN[M-1][j1];
			}
			for (i1=2; i1<=L; i1++){    //envelhecimento da TN
				for (j1=1; j1<=M; j1++){
					T[i1][j1]=TN[i1][j1];
				}
			}
		}  // end do for do kk



	
	} // end do master
	
	MPI_Finalize();

	double end = omp_get_wtime();
	if(proc == 0){
		printf("nc=%d\n", nc);
		printf("%.5lf\n", end - start);
	}
} //end main















//***************************************************** KA *******************************************************
// FUNO PARA CALCULAR A CONDITIVIDADE HIDRLICA COM BASE NO TETA

double KA(double x1, double x2){    
	double y=0;
	double a=0.467; // coeficiente do Psi na equao do potencial 
	double m=0.934; // expoente na equao do potencial 
	double n=2/(1-m); // expoente na equao do potencial, Burdines 
	if(x2 <= 0)
		y=0;
	else{
		if (x2 <=1){
		//Bordines
		pow1 = 1 - pow(x2,(1/m));          
		pow2 = 1 - pow(pow1,m);
		pow3 = pow2*pow2;
		y = x1 * pow(x2,0.5) * pow3;
		}
		else
			y=x1;
	}
	return(y);
}

//**************************************************** POT ********************************************************
// FUNO PARA CALCULAR O POTENCIAL MATRICIAL COM BASE NO TETA

double POT(double x){
	double y=0;
	double a=0.674; //coeficiente do Psi na equao do potencial 
	double m=0.768; // expoente na equao do potencial 
	double n=2/(1-m); // expoente na equao do potencial, Burdines 
	if (x <= 0)
		y= 0.35;
	else{      
		if (x <=1){
		pow1 = pow(x,(-1/m))-1;
			pow2 = pow((pow1),1/n);
		y=-pow(10,(1/a*pow2)); // pot matricial 
		}
		else
			y=0;
	}
	return(y);
}

//************************************************** evaporação ****************************************************
// FUNO PARA CCULO DOS PARETROS DA EVAPORAO---AJUSTE DE CURVAS EM Z=0  EVAPORAO-----

double evaporacao(){
	double valor=0;
	double linha1[4];
	double linha2[4];
	double laux[4];
	double soma=0, soma1=0, soma2=0;
	int nump=4; // nmero de pontos
	int i,z,y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13;	
	tt[1][1]=0;
	tt[2][1]=10;
	tt[3][1]=20;
	tt[4][1]=30;
	TEv[1][1]=T00;
	TEv[2][1]=0.75*T00;
	TEv[3][1]=0.5*T00;
	TEv[4][1]=0.3*T00;
	
	for (i=1; i <=nump; i++){
		A[i][1]=tt[i][1];
		A[i][2]=1;
	}
	// x=inv(A'*A)*A'*TEv 
	for (y0 = 1; y0 <=4; y0++){         
		for (z = 1; z <=2; z++){                    // tranposta de A'
			mat1[z][y0] = A[y0][z];
			mataux[z][y0] = A[y0][z];
		}
	}
	for (y1 = 1; y1 <=4; y1++){                       //
		soma = soma + (mat1[1][y1] * A[y1][1]);   //
	}                                                 //
	A2[1][1] = soma;                                  //
	soma = 0;                                         //
	for (y2 = 1; y2 <=4; y2++){                       //
		soma = soma + (mat1[1][y2] * A[y2][2]);   //
	}                                                 //
	A2[1][2] = soma;                                  // multiplicao de A'*A
	soma = 0;                                         //  A2 = (mat1 * A) 
	for (y3 = 1; y3 <=4; y3++){                       //
		soma = soma + (mat1[2][y3] * A[y3][1]);   //
	}                                                 //
	A2[2][1] = soma;                                  //
	soma = 0;                                         //
	for (y4 = 1; y4 <=4; y4++){                       //
		soma = soma + (mat1[2][y4] * A[y4][2]);   //
	}                                                 //
	A2[2][2] = soma;                                  //
	soma = 0;                                         //
								
	// A SEGUIR MATRIZ INVERSA DE A2 E GRAVA EM MAT2
	linha1[0] = A2[1][1];
	linha1[1] = A2[1][2];
	linha1[2] = 1;
	linha1[3] = 0;                           //definindo a malha
	linha2[0] = A2[2][1];
	linha2[1] = A2[2][2];
	linha2[2] = 0;
	linha2[3] = 1;
	//1 regra - o valor de linha1[0] deverser 1 - capturar o valor e dividir linha1 por ele
	valor = linha1[0];
	linha1[0] = linha1[0] / valor;
	linha1[1] = linha1[1] / valor;
	linha1[2] = linha1[2] / valor;
	linha1[3] = linha1[3] / valor;
	//2 regra - zerar o valor de linha2[0] - capturar o valor e multiplicar por (-) ele
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
	//3 regra - o valor de linha2[1] deverser 1 - capturar o valor e dividir linha2 por ele
	valor = linha2[1];
	linha2[0] = linha2[0] / valor;
	linha2[1] = linha2[1] / valor;
	linha2[2] = linha2[2] / valor;
	linha2[3] = linha2[3] / valor;
	//4 regra - zerar o valor de linha1[1] - capturar o valor e multiplicar por (-) ele
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
	for (y5 = 1; y5 <=2; y5++){                              
		soma = soma + (mataux[y5][1] * mat2[y5][1]);     
	}                                                    
	mat3[1][1] = soma;                                    
	soma = 0;                                            
	for (y6 = 1; y6 <=2; y6++){                              
		soma = soma + (mataux[y6][2] * mat2[y6][1]);     
	}                                                    
	mat3[1][2] = soma;                                    
	soma = 0;                                            
	for (y7 = 1; y7 <=2; y7++){                              
		soma = soma + (mataux[y7][3] * mat2[y7][1]);     
	}                                                    
	mat3[1][3] = soma;                                    
	soma = 0;                                                
	for (y8 = 1; y8 <=2; y8++){                                   //    multiplicao
		soma = soma + (mataux[y8][4] * mat2[y8][1]);          //    mat3 =inv(A'*A)*A' 
	}                                                    
	mat3[1][4] = soma;                                    
	soma = 0;                                                
	for (y9 = 1; y9 <=2; y9++){                              
		soma = soma + (mataux[y9][1] * mat2[y9][2]);     
	}
	mat3[2][1] = soma;                                    
	soma = 0;                                                
	for (y10 = 1; y10 <=2; y10++){                          
		soma = soma + (mataux[y10][2]*mat2[y10][2]); 
	}                      
	mat3[2][2] = soma;                                   
	soma = 0;                                         
	for (y11 = 1; y11 <=2; y11++){                
		soma = soma + (mataux[y11][3]*mat2[y11][2]);  
	}                                                 
	mat3[2][3] = soma;                                 
	soma = 0;                                         
	for (y12 = 1; y12 <=2; y12++){                
		soma = soma + (mataux[y12][4]*mat2[y12][2]); 
	}                                                 
	mat3[2][4] = soma;                                
	soma = 0;                                         
		
	//FINALIZANDO x = mat3 * TEv;
	for (y13 = 1; y13 <=4; y13++){                    
		soma1 = soma1 + (mat3[1][y13]*TEv[y13][1]);      
		soma2 = soma2 + (mat3[1+1][y13]*TEv[y13][1]);    
	}                                                 
	x[1][1] = soma1;                                   
	x[2][1] = soma2;           
}

//****************************************************** Worker ***************************************************

int worker(){
	double To[L][M];
	double evap1 = 0;
	double evap2 = 0;
	int pos_ini = 0;
	int pos_fim = 0;
	int tamanho = 0;
	int posi = 0;
	int j0, j1, j2, kk, i, i1, jj, ii;	

	for (j0=1; j0<=M; j0++){
		for (i=1; i<=L; i++){
			T[i][j0]=T00;  //condio inicial para toda a malha
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
	
	//------------PROBLEMA INVERSO E DIRETO---------------------
	
	MPI_Status status;
	MPI_Recv(&tamanho, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
	MPI_Recv(&evap1, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
	MPI_Recv(&evap2, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
	MPI_Recv(&pos_ini, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
	MPI_Recv(&pos_fim, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
//	printf(" \n\n tamanho %d", tamanho); 
//	printf(" \n\n evap1 %f", evap1); 
//	printf(" \n\n evap2 %f", evap2); 
//	printf(" \n\n pos_ini %d", pos_ini); 
//	printf(" \n\n pos_fim %d", pos_fim); 
	double vetor_erro[tamanho];
	double vetor_qot[tamanho];

	for (ii=pos_ini; ii <= pos_fim; ii++){ 	// contador dos chutes do PI
		q=qo+(ii-1)*dq;
		for (j1=1; j1<=M; j1++){
			for (i=1; i<=L; i++){
			T[i][j1]=T00;  //condio inicial para toda a malha 
		}        
		}
		for (j2=1; j2<=L; j2++){
			for (i=1; i<=M; i++) {
				Sf[i][j2]=0;  //fonte zero para toda a malha 
			}
		}
		for (kk=1; kk<=nn; kk++){   //TEMPORAL
		te[kk][1]=(kk-1)*dt;
		qq=q*exp(-0.05*((te[kk][1]-12)*(te[kk][1]-12))); 
		for (i=2; i<=(L-1); i++){
			Tsup=evap1*te[kk][1]+evap2; //teta em z=0 pela evaporao 
			for (i1=1; i1<=L; i1++){
				T[1][i1]=Tsup; //condio de fronteira em z=0                          
			}
			for (j=2; j<=(M-1); j++){                       
			//------------FONTE = TRANSPIRAO ----------------
				q1=qq;
				q2=qq;
				q3=qq;
				Sf[2][2]= qq*(T[2][2]-Tr); Sf[3][2]= qq*(T[3][2]-Tr);		
				Sf[4][2]= qq*(T[4][2]-Tr); Sf[4][3]=q1*(T[4][3]-Tr);
				Sf[5][2]= qq*(T[5][2]-Tr); Sf[5][3]=q1*(T[5][3]-Tr); Sf[5][4]=q2*(T[5][4]-Tr);
				Sf[6][2]= qq*(T[6][2]-Tr); Sf[6][3]=q1*(T[6][3]-Tr);
				Sf[6][4]=q2*(T[6][4]-Tr);Sf[6][5]=q3*(T[6][5]-Tr);
				Sf[7][2]= qq*(T[7][2]-Tr); Sf[7][3]=q1*(T[7][3]-Tr); Sf[7][4]=q2*(T[7][4]-Tr);
				Sf[8][2]= qq*(T[8][2]-Tr); Sf[8][3]=q1*(T[8][3]-Tr);
				Sf[9][2]= qq*(T[9][2]-Tr); Sf[10][2]= qq*(T[10][2]-Tr);
				r=((j-1)*dr)+(0.5 * dr);   		           
				KE= KA(Ko,T[i][j+1]);
				KW= KA(Ko,T[i][j-1]);
				KP= KA(Ko,T[i][j]);
				KN= KA(Ko,T[i+1][j]);
				KS= KA(Ko,T[i-1][j]);
				AE= dt/dr*((KE-KW)/(4*dr)+KP/(2*r)+KP/dr);
				AW= dt/dr*(-(KE-KW)/(4*dr)-KP/(2*r)+KP/dr);
				AN= dt/(4*(dz*dz))*(KN-KS+4*KP);  
				AS= dt/(4*(dz*dz))*(-KN+KS+4*KP); 
				APO= -(2*dt*KP*(((1/dr)*(1/dr))+((1/dz)*(1/dz))));            
				PME=POT(T[i][j+1]);
				PMW=POT(T[i][j-1]);
				PMP=POT(T[i][j]);
				PMN=POT(T[i+1][j]);
				PMS=POT(T[i-1][j]);
				TN[i][j]=(AE*PME+AW*PMW+AN*PMN+AS*PMS+APO*PMP)+(KN-KS)/(2*dz)+T[i][j]-Sf[i][j]; 
			} //end do for do j
			
			TN[i][1]= TN[i][2]; // CONDIO DE FRONTEIRA EM r=R, isolamento
			TN[i][L]= TN[i][L-1]; // CONDIO DE FRONTEIRA EM r=R, isolamento        
			
		}//end do for do i
		for (j1=1; j1<=L; j1++){ // condio de fronteira em z=h, isolamento
			TN[M][j1]= TN[M-1][j1];    
		}
		for (i1=2; i1<=L; i1++){    //envelhecimento da TN
			for (j1=1; j1<=M; j1++){
				T[i1][j1]= TN[i1][j1];   
			}
		}
		//-------PROBLEMA INVERSO-----------                             
		if (fabs(te[kk][1]-10) <= dt){ //10 horas     
			SQE10=(T1[1][2]-T[1][2])*(T1[1][2]-T[1][2])+
			(T1[1][5]-T[1][5])*(T1[1][5]-T[1][5])+(T1[1][8]-T[1][8])*(T1[1][8]-T[1][8])+
			(T1[3][2]-T[3][2])*(T1[3][2]-T[3][2])+(T1[3][5]-T[3][5])*(T1[3][5]-T[3][5])+
			(T1[3][8]-T[3][8])*(T1[3][8]-T[3][8])+(T1[6][2]-T[6][2])*(T1[6][2]-T[6][2])+
			(T1[6][5]-T[6][5])*(T1[6][5]-T[6][5])+(T1[6][8]-T[6][8])*(T1[6][8]-T[6][8])+
			(T1[11][2]-T[11][2])*(T1[11][2]-T[11][2])+(T1[11][5]-T[11][5])*(T1[11][5]-T[11][5])+
			(T1[11][8]-T[11][8])*(T1[11][8]-T[11][8]);    
		}
		if (fabs(te[kk][1]-20) <= dt){ //20 horas
			SQE20=(T2[1][2]-T[1][2])*(T2[1][2]-T[1][2])+
			(T2[1][5]-T[1][5])*(T2[1][5]-T[1][5])+(T2[1][8]-T[1][8])*(T2[1][8]-T[1][8])+
			(T2[3][2]-T[3][2])*(T2[3][2]-T[3][2])+(T2[3][5]-T[3][5])*(T2[3][5]-T[3][5])+
			(T2[3][8]-T[3][8])*(T2[3][8]-T[3][8])+(T2[6][2]-T[6][2])*(T2[6][2]-T[6][2])+
			(T2[6][5]-T[6][5])*(T2[6][5]-T[6][5])+(T2[6][8]-T[6][8])*(T2[6][8]-T[6][8])+
			(T2[11][2]-T[11][2])*(T2[11][2]-T[11][2])+(T2[11][5]-T[11][5])*(T2[11][5]-T[11][5])+
			(T2[11][8]-T[11][8])*(T2[11][8]-T[11][8]);    
		}
		if (fabs(te[kk][1]-29) <= dt){ //30 horas
			SQE30=(T3[1][2]-T[1][2])*(T3[1][2]-T[1][2])+
			(T3[1][5]-T[1][5])*(T3[1][5]-T[1][5])+(T3[1][8]-T[1][8])*(T3[1][8]-T[1][8])+
			(T3[3][2]-T[3][2])*(T3[3][2]-T[3][2])+(T3[3][5]-T[3][5])*(T3[3][5]-T[3][5])+
			(T3[3][8]-T[3][8])*(T3[3][8]-T[3][8])+(T3[6][2]-T[6][2])*(T3[6][2]-T[6][2])+
			(T3[6][5]-T[6][5])*(T3[6][5]-T[6][5])+(T3[6][8]-T[6][8])*(T3[6][8]-T[6][8])+
			(T3[11][2]-T[11][2])*(T3[11][2]-T[11][2])+(T3[11][5]-T[11][5])*(T3[11][5]-T[11][5])+
			(T3[11][8]-T[11][8])*(T3[11][8]-T[11][8]);  
		}
		} //end do for do kk
		vetor_erro[posi] = SQE10+SQE20+SQE30;
		vetor_qot[posi] = q;
		posi++;
	} //end do ii, Problema inverso 
	//retornar vetor
	MPI_Send(vetor_erro, tamanho, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
	MPI_Send(vetor_qot, tamanho, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
}
