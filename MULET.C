//**************************************************
//**  Auteur : SCHUTZ Antony
//**  Section : Licence de physique
//**  Session : 2000 / 2001
//**  But : Les étoiles multiples
//**  Prg : mulet.c
//**  Responsable : Mr Aristidi
//**************************************************
//**************************************************
//**  DECLARATION DES LIBRAIRIES
//**************************************************
#include <stdio.h>
#include <math.h>
#include "/usr/include/cpgplot.h"
//** G En annee lumiere au cube par masse solaire
//** et par million d'annee au carre
#define G  0.15656
#define MAXINT (pow(2.,31.)-1)
//** perspective
#define PI 3.141592654
//** saisie max et min de la souris
#define MAX(a,b) ( (a)>(b) ? (a) : (b)) 
#define MIN(a,b) ( (a)<(b) ? (a) : (b)) 
//**************************************************
//**  DECLARATIONS DES VARIABLES GLOBALES
//**************************************************
//**  Les caracteristiques des Etoiles
struct CaracteristiqueDesCorps
{   double Masse,Position[3],Vitesse[3]; };             
//**  Les Vecteurs intermediaires    
struct Intermediaire
{   double Axe[3]; };
//**  Les parametres de Runge Kutta
struct ParametreDeRungeKutta
{    double X[4],Y[4],Z[4];};
//** Instance de structure
typedef struct CaracteristiqueDesCorps Corps;
typedef struct Intermediaire Vecteur;
typedef struct ParametreDeRungeKutta RK;
//**************************************************
//**  PROTOTYPES DES FONCTIONS
//**************************************************
void main ( void ) ;
void Process ( int choix, int NombreDeCorps, Corps * Etoile,Vecteur * VectPost, Vecteur * VectForce, RK * KX, RK * KV );
//** calcul
//*********
double CarreDe ( double );
double CubeDe ( double );
void acceleration ( int EnCour, Corps * Etoile, int NombreDeCorps, Vecteur * VectPost, Vecteur * VectForce );
void RungeKutta ( int NombreDeCorps, double PasDeTemp,Vecteur * VectPost, Vecteur * VectForce, RK * KX, RK * KV ,Corps * Etoile);
double EnergieCinetique ( int NombreDeCorps , Corps * Etoile );
double Potentiel ( int NombreDeCorps , Corps * Etoile , double * Corps1 );
double EnergiePotentiel ( int NombreDeCorps , Corps * Etoile );
void InitAlea ( int NombreDeCorps, int alea, Corps * Etoile );
float InitSystConnu ( int choix, Corps * Etoile );
double SommeDesMasse ( int NombreDeCorps , Corps * Etoile );
double DistanceMax ( int choix, int NombreDeCorps , Corps * Etoile );
double CalculduTempTho (  int NombreDeCorps , Corps * Etoile );
//** affichage graphique
//**********************
void TraceEnergie ( int NombreDeCorps,int NbExecut, double PasDeTemp, Vecteur * VectPost, Vecteur * VectForce, RK * KX, RK * KV ,Corps * Etoile );
void Systeme ( int NbExecut, int NombreDeCorps, int Tcolor, double PasDeTemp,int Type,float max,Vecteur * VectPost, Vecteur * VectForce, RK * KX, RK * KV ,Corps * Etoile );
void TraceRefresh ( int Tcolor, int NombreDeCorps, float * x,float * y, float * z );
void TraceCercle( int Tcolor, int NombreDeCorps,  float * x,float * y, float * z, float * m ) ;
void AffichInfo ( int IndiceDuCorps,int NombreDeCorps,double Tho,double PasDeTemps, Corps * Etoile, float xmin, float ymax );
void FondSyst ( float max );
void FondEner ( void ) ;
void FondOpt ( float max );
void TraceOpt (  int NbExecut, int NombreDeCorps, int Tcolor, double PasDeTemp,float max,Vecteur * VectPost, Vecteur * VectForce, RK * KX, RK * KV ,Corps * Etoile ) ;
//** affichage et initalisation
//*****************************
void Presentation ( void );
void Menu ( void );
int  SaisieNCorps ( void );
void InitCorps ( int NombreDeCorps , Corps * Etoile);
void ConfirmInit ( int NombreDeCorps , Corps * Etoile);
double InitPasDeTemp ( void );
float InitDim ( void );
int TypeInit ( void );
int TypeAffich ( void );
int TypeCorps ( void ) ;
int Nbiteration ( void );
void MenuChoixSyst( Corps * Etoile,Vecteur * VectPost, Vecteur * VectForce, RK * KX, RK * KV  );
//**************************************************
//**  PROGRAMME PRINCIPALE
//**************************************************
void main ( void )
{
  int Choix,Choix1,Nombre,NombreDeCorps;
  //** Pointeur sur structure
  Corps * Etoile;
  Vecteur * VectPost, * VectForce; 
  RK * KX, *KV;
  //  debut du programme
  Presentation();
  for(;;)
    {
      Menu();
      scanf("%d",&Choix);
      switch(Choix)
	{
	case(1) : 
	  NombreDeCorps=2;
	  Process(Choix,NombreDeCorps,Etoile,VectPost,VectForce,KX,KV);
	  break;	  	  
	case(2) :   
	  NombreDeCorps=SaisieNCorps();
	  Process(Choix,NombreDeCorps,Etoile,VectPost,VectForce,KX,KV);
	  break;	 	  
	case(3):
	  NombreDeCorps=SaisieNCorps();
	  Process(Choix,NombreDeCorps,Etoile,VectPost,VectForce,KX,KV);
	  break;	  
	case(4):
	  NombreDeCorps=SaisieNCorps();
	  Process(Choix,NombreDeCorps,Etoile,VectPost,VectForce,KX,KV);
	  break;
	case(5):
	  MenuChoixSyst(Etoile,VectPost,VectForce,KX,KV);
	  break;
	case(6):
	  exit(0);
	  break;  
	default : 
	  printf("**\n** hein hein hein essaye encore \n**");
	  break;  
	}
    }
}
//**  procedure generale
//**********************
void Process ( int choix, int NombreDeCorps, Corps * Etoile,Vecteur * VectPost, Vecteur * VectForce, RK * KX, RK * KV  )
{
  int Tcolor,Type,NbExecut,verif,Init;
  double PasDeTemp;
  float dim=200;
  Etoile=(Corps *)calloc(NombreDeCorps,sizeof(Corps));
  KX=(RK *)calloc(NombreDeCorps,sizeof(RK));
  KV=(RK *)calloc(NombreDeCorps,sizeof(RK));
  VectPost=(Vecteur *)calloc(1,sizeof(Vecteur));
  VectForce=(Vecteur *)calloc(1,sizeof(Vecteur));
  Init=TypeInit();
  PasDeTemp=InitPasDeTemp();
  NbExecut=Nbiteration();
  if(Init==1)
    InitAlea(NombreDeCorps,(int)rand(),Etoile);
else
  InitCorps (NombreDeCorps ,Etoile);
  switch(choix)
    {
    case(1):
      Tcolor=TypeAffich();
      Type=TypeCorps();
      dim=InitDim();
      Systeme(NbExecut,NombreDeCorps,Tcolor,PasDeTemp,Type,dim,VectPost,VectForce,KX,KV ,Etoile);
      break;
    case(2):
      Tcolor=TypeAffich();
      Type=TypeCorps();
      dim=InitDim();
      Systeme(NbExecut,NombreDeCorps,Tcolor,PasDeTemp,Type,dim,VectPost,VectForce,KX,KV ,Etoile);
      break;
    case(3):
      TraceEnergie(NombreDeCorps,NbExecut,PasDeTemp,VectPost,VectForce,KX,KV,Etoile);
      break;
    case(4):
      dim=InitDim();
      Tcolor=TypeAffich();
      TraceOpt (NbExecut,NombreDeCorps,Tcolor,PasDeTemp,dim,VectPost,VectForce,KX,KV,Etoile);  
      break;
    }
  free(Etoile);
  free(KX);
  free(KV);
  free(VectPost);
  free(VectForce);
}
//**************************************************
//**  FONCTIONS CALCUL
//**************************************************
//** calcul du carre d un nombre
//******************************
double CarreDe ( double x )
{return(x*x);}
//**  calcul du cube d un nombre
//******************************
double CubeDe ( double x )
{return(x*x*x);}
//**  Calcul de la force subit par un corp
//****************************************
void acceleration ( int EnCour, Corps * Etoile, int NombreDeCorps, Vecteur * VectPost, Vecteur * VectForce )
{
  int indice;
  double Distance, Tampon;
  for(indice=0;indice<3;indice++)  //  initialisation
    VectPost[0].Axe[indice]=0;            //  du tableau
  for(indice=0;indice<NombreDeCorps;indice++)
    {
      if(indice!=EnCour)  //  si ce n est pas elle meme
	{
	  Distance=sqrt(CarreDe(Etoile[indice].Position[0]-VectForce[0].Axe[0])+CarreDe(Etoile[indice].Position[1]-VectForce[0].Axe[1])+CarreDe(Etoile[indice].Position[2]-VectForce[0].Axe[2]));  
	  // Tampon de calcul
	  Tampon= ( G * Etoile[EnCour].Masse * Etoile[indice].Masse ) / CubeDe(Distance);  
	  // modification de la position
	  VectPost[0].Axe[0]+=(Etoile[indice].Position[0]-VectForce[0].Axe[0]) * Tampon ;
	  VectPost[0].Axe[1]+=(Etoile[indice].Position[1]-VectForce[0].Axe[1]) * Tampon ;
	  VectPost[0].Axe[2]+=(Etoile[indice].Position[2]-VectForce[0].Axe[2]) * Tampon ;	  
	}
    }
}
//**  Algorithme de Runge Kutta
//*****************************
void RungeKutta ( int NombreDeCorps, double PasDeTemp,Vecteur * VectPost, Vecteur * VectForce, RK * KX, RK * KV ,Corps * Etoile)
{
  int indice,i;
  for(indice=0;indice<NombreDeCorps;indice++)
    {
      //**//
      KX[indice].X[0]=Etoile[indice].Vitesse[0]*PasDeTemp;
      KX[indice].Y[0]=Etoile[indice].Vitesse[1]*PasDeTemp;
      KX[indice].Z[0]=Etoile[indice].Vitesse[2]*PasDeTemp;
      for(i=0;i<3;i++)
	VectForce[0].Axe[i]=Etoile[indice].Position[i];
      acceleration(indice,Etoile,NombreDeCorps,VectPost,VectForce);
      KV[indice].X[0]=VectPost[0].Axe[0]*PasDeTemp;
      KV[indice].Y[0]=VectPost[0].Axe[1]*PasDeTemp;
      KV[indice].Z[0]=VectPost[0].Axe[2]*PasDeTemp;
      //**//
      KX[indice].X[1]=(KX[indice].X[0]+(KV[indice].X[0]/2))*PasDeTemp;
      KX[indice].Y[1]=(KX[indice].Y[0]+(KV[indice].Y[0]/2))*PasDeTemp;
      KX[indice].Z[1]=(KX[indice].Z[0]+(KV[indice].Z[0]/2))*PasDeTemp;
      VectForce[0].Axe[0]=((Etoile[indice].Position[0])+(KX[indice].X[0]/2));
      VectForce[0].Axe[1]=((Etoile[indice].Position[1])+(KX[indice].Y[0]/2));
      VectForce[0].Axe[2]=((Etoile[indice].Position[2])+(KX[indice].Z[0]/2));
      acceleration(indice,Etoile,NombreDeCorps,VectPost,VectForce);
      KV[indice].X[1]=VectPost[0].Axe[0]*PasDeTemp;
      KV[indice].Y[1]=VectPost[0].Axe[1]*PasDeTemp;
      KV[indice].Z[1]=VectPost[0].Axe[2]*PasDeTemp;
      //**//
      KX[indice].X[2]=(KX[indice].X[0]+(KV[indice].X[1]/2))*PasDeTemp;
      KX[indice].Y[2]=(KX[indice].Y[0]+(KV[indice].Y[1]/2))*PasDeTemp;
      KX[indice].Z[2]=(KX[indice].Z[0]+(KV[indice].Z[1]/2))*PasDeTemp;
      VectForce[0].Axe[0]=((Etoile[indice].Position[0])+(KX[indice].X[1]/2));
      VectForce[0].Axe[1]=((Etoile[indice].Position[1])+(KX[indice].Y[1]/2));
      VectForce[0].Axe[2]=((Etoile[indice].Position[2])+(KX[indice].Z[1]/2));
      acceleration(indice,Etoile,NombreDeCorps,VectPost,VectForce);
      KV[indice].X[2]=VectPost[0].Axe[0]*PasDeTemp;
      KV[indice].Y[2]=VectPost[0].Axe[1]*PasDeTemp;
      KV[indice].Z[2]=VectPost[0].Axe[2]*PasDeTemp;
      //**//
      KX[indice].X[3]=(KX[indice].X[0]+(KV[indice].X[2]/2))*PasDeTemp;
      KX[indice].Y[3]=(KX[indice].Y[0]+(KV[indice].Y[2]/2))*PasDeTemp;
      KX[indice].Z[3]=(KX[indice].Z[0]+(KV[indice].Z[2]/2))*PasDeTemp;
      VectForce[0].Axe[0]=((Etoile[indice].Position[0])+(KX[indice].X[2]));
      VectForce[0].Axe[1]=((Etoile[indice].Position[1])+(KX[indice].Y[2]));
      VectForce[0].Axe[2]=((Etoile[indice].Position[2])+(KX[indice].Z[2]));
      acceleration(indice,Etoile,NombreDeCorps,VectPost,VectForce);
      KV[indice].X[3]=VectPost[0].Axe[0]*PasDeTemp;
      KV[indice].Y[3]=VectPost[0].Axe[1]*PasDeTemp;
      KV[indice].Z[3]=VectPost[0].Axe[2]*PasDeTemp;
    }
  for(indice=0;indice<NombreDeCorps;indice++)
    {
      // modification de la position      
      Etoile[indice].Position[0]+=(KX[indice].X[0]+KX[indice].X[3])/6+(KX[indice].X[1]+KX[indice].X[2])/3;      
      Etoile[indice].Position[1]+=(KX[indice].Y[0]+KX[indice].Y[3])/6+(KX[indice].Y[1]+KX[indice].Y[2])/3;      
      Etoile[indice].Position[2]+=(KX[indice].Z[0]+KX[indice].Z[3])/6+(KX[indice].Z[1]+KX[indice].Z[2])/3;
      // modification de la vitesse      
      Etoile[indice].Vitesse[0]+=(KV[indice].X[0]+KV[indice].X[3])/6+(KV[indice].X[1]+KV[indice].X[2])/3;      
      Etoile[indice].Vitesse[1]+=(KV[indice].Y[0]+KV[indice].Y[3])/6+(KV[indice].Y[1]+KV[indice].Y[2])/3;      
      Etoile[indice].Vitesse[2]+=(KV[indice].Z[0]+KV[indice].Z[3])/6+(KV[indice].Z[1]+KV[indice].Z[2])/3;
    }
}
//**  calcul de l energie cinetique
//*********************************
double EnergieCinetique ( int NombreDeCorps , Corps * Etoile )
{
  int indice;
  double Ec=0;
  for(indice=0;indice<NombreDeCorps;indice++)
    {           
      Ec+=(double)0.5*(Etoile[indice].Masse)*((CarreDe(Etoile[indice].Vitesse[0]))+(CarreDe(Etoile[indice].Vitesse[1]))+(CarreDe(Etoile[indice].Vitesse[2])));
    }
  return(Ec);
}
//**  calcul du potentielle en un point
//*************************************
double Potentiel ( int NombreDeCorps , Corps * Etoile , double * Corps1 )
{
  int indice;
  double Potentiel;
  double * Corps2;
  Corps2=(double*)calloc(3,sizeof(double));
  for(indice=0;indice<NombreDeCorps;indice++)
    {
      Corps2[0]=Etoile[indice].Position[0];
      Corps2[1]=Etoile[indice].Position[1];
      Corps2[2]=Etoile[indice].Position[2];
      if(((Corps1[0]-Corps2[0])!=0)&&((Corps1[1]-Corps2[1])!=0)&&((Corps1[2]-Corps2[2])!=0))
	Potentiel+=((Etoile[indice].Masse)/sqrt(CarreDe(Corps1[0]-Corps2[0])+CarreDe(Corps1[1]-Corps2[1])+CarreDe(Corps1[2]-Corps2[2])));
    }
  return(Potentiel*G);
}
//**  calcul de l energie potentielle
//***********************************
double EnergiePotentiel ( int NombreDeCorps , Corps * Etoile )
{
  int indice;
  double Ep;
  double * Corps1;
  Corps1=(double*)calloc(3,sizeof(double));
  for(indice=0;indice<NombreDeCorps;indice++)
    {
      Corps1[0]=Etoile[indice].Position[0];
      Corps1[1]=Etoile[indice].Position[1];
      Corps1[2]=Etoile[indice].Position[2];      
      Ep+= Potentiel(NombreDeCorps, Etoile, Corps1 ) * Etoile[indice].Masse ;      
    }
  return(Ep);
}
//**  initialisation aleatoire des coordonnees
//********************************************
void InitAlea( int NombreDeCorps, int alea, Corps * Etoile )
{
  int indice;
  float a,b,c;
  srandom(getpid()+alea);
  for (indice=0;indice<NombreDeCorps;indice++)
    { 
      a=(float)rand()/MAXINT*200; 
      b=(float)rand()/MAXINT*3.14;
      c=(float)rand()/MAXINT*6.28;
      Etoile[indice].Position[0] = a * sin (b) * cos (c);
      Etoile[indice].Position[1] = a * sin (b) * sin (c);
      Etoile[indice].Position[2] = a * cos (b);
      Etoile[indice].Vitesse[0]=(float)rand()/MAXINT*5-2.5; 
      Etoile[indice].Vitesse[1]=(float)rand()/MAXINT*5-2.5; 
      Etoile[indice].Vitesse[2]=(float)rand()/MAXINT*5-2.5;
      Etoile[indice].Masse=(float)rand()/MAXINT*27.7+0.3;
    }                                    
}
//** calcul de tho ou temps de revolution pour deux corps
//*******************************************************
double CalculduTempTho ( int NombreDeCorps , Corps * Etoile )
{
  double SDMasse;
  double xmax,ymax,zmax,rmax;
  double Tho;
  xmax= DistanceMax(0,NombreDeCorps,Etoile);
  ymax= DistanceMax(1,NombreDeCorps,Etoile);
  zmax= DistanceMax(2,NombreDeCorps,Etoile);
  rmax=sqrt(CarreDe(xmax)+CarreDe(ymax)+CarreDe(zmax));
  printf("\nxmax : %lg",xmax);
  printf("\nymax : %lg",ymax);
  printf("\nzmax : %lg",zmax);
  printf("\nrmax : %lg",rmax);	
  SDMasse=SommeDesMasse (NombreDeCorps,Etoile);
  printf("\ndsm : %lg",SDMasse);
  printf("\ng : %lg",G);
  Tho=sqrt( (4*CarreDe(PI)*CubeDe(rmax) ) / (G*SDMasse) );
  printf("\ntho : %lg",Tho);
  return(Tho);
}
//**  Calcul la somme des masses
//******************************
double SommeDesMasse ( int NombreDeCorps , Corps * Etoile )
{
  int indice;
  double MasseTotale=0;
  for(indice=0;indice<NombreDeCorps;indice++)
    {
      MasseTotale+=Etoile[indice].Masse; 
      printf("\n%d masse : %lg",indice,Etoile[indice].Masse);
    }
  return(MasseTotale);
}
//** calcul la distance maximum entre deux corps
//**********************************************
double DistanceMax ( int choix, int NombreDeCorps , Corps * Etoile )
{
  int indice;
  double DMax,DMin;
      switch(choix)
	{
	case(0):DMax=Etoile[0].Position[0];DMin=DMax;break;
	case(1):DMax=Etoile[0].Position[1];DMin=DMax;break;
	case(2):DMax=Etoile[0].Position[2];DMin=DMax;break;
	}
  for(indice=1;indice<NombreDeCorps;indice++)
    {
      switch(choix)
	{
	case(0):
	  DMax=MAX(DMax,Etoile[indice].Position[0]);
	  DMin=MIN(DMin,Etoile[indice].Position[0]);
	  break;
	case(1):
	  DMax=MAX(DMax,Etoile[indice].Position[1]);
	  DMin=MIN(DMin,Etoile[indice].Position[1]);
	  break;
	case(2):
	  DMax=MAX(DMax,Etoile[indice].Position[2]);
	  DMin=MIN(DMin,Etoile[indice].Position[2]);
	  break;
	}
    }
  DMax=DMax-DMin;
  return(DMax);
}
//**************************************************
//**  FONCTIONS D'AFFICHAGE GRAPHIQUES
//**************************************************
//** tracer de la conservation de l'energie
//*****************************************
void TraceEnergie (int NombreDeCorps,int NbExecut,double PasDeTemp,Vecteur * VectPost, Vecteur * VectForce, RK * KX, RK * KV , Corps * Etoile)
{
  int indice,ind,style,choix;
  float Ec,Ep,abs[1],ord[1],Et,ttemp=0;
  float * x, * y, * z, * x1, * y1, * x2, * y2, * x3, * y3;
  char Info0[100],Info1[100],Info2[100],temp[50],thot[50];
  double Thot;
  x=(float *)calloc(NombreDeCorps,sizeof(float));
  y=(float *)calloc(NombreDeCorps,sizeof(float));
  z=(float *)calloc(NombreDeCorps,sizeof(float));
  x1=(float *)calloc(1000,sizeof(float));
  y1=(float *)calloc(1000,sizeof(float));
  x2=(float *)calloc(1000,sizeof(float));
  y2=(float *)calloc(1000,sizeof(float));
  x3=(float *)calloc(1000,sizeof(float));
  y3=(float *)calloc(1000,sizeof(float));
  style=4;
  printf("**  projection x y ( 1 ) ou perspective ( 0 ) : ");
  scanf("%d",&choix);
  if(choix==0)
    {
      for(indice=0;indice<1000;indice++)
	{
	  x1[indice]=0.5*((indice)*cos(PI/6));
	  y1[indice]=0.000005*((indice)*sin(PI/6)-1000);
	  x2[indice]=0.5*(indice+(1000)*cos(PI/6));
	  y2[indice]=0.000005*((1000)*sin(PI/6)-1000);
	  x3[indice]=0.5*((1000)*cos(PI/6));
	  y3[indice]=0.000005*(indice+(1000)*sin(PI/6)-1000);
	}
    }
  cpgopen("/xw");
  cpgbeg(0,"/xw",2,2);
  cpgscf(2.0);
  FondEner();
  Et=(float)EnergieCinetique(NombreDeCorps,Etoile)+(float)EnergiePotentiel(NombreDeCorps,Etoile);
  Thot=CalculduTempTho(NombreDeCorps,Etoile);
  for(indice=0;indice<NbExecut;indice++)
    {
      Ec=(float)EnergieCinetique(NombreDeCorps,Etoile);
      Ep=(float)EnergiePotentiel(NombreDeCorps,Etoile);
      ttemp=indice*PasDeTemp;
      for(ind=0;ind<NombreDeCorps;ind++)
	{
	  ord[0]=(Ec+Ep-Et)/Et;
	  abs[0]=(indice+1)*(1000/NbExecut);
	  sleep(0.5);
	  //  efface les dernieres info avant d'afficher les nouvelles
	  if(indice!=0)
	    {
	      cpgpanl(1,2); 
	      cpgsci(0);
	      cpgscf(2);
	      cpgsch(1.5);
	      cpgtext(50,0.003,Info0);
	      cpgtext(50,0.0015,Info1);
	      cpgtext(50,-0,Info2);
	      cpgtext(50,-0.0015,temp);
	      cpgtext(50,-0.003,thot);
	    }
	  sprintf(Info0,"Energie Cinetique : %lg",Ec);
	  sprintf(Info1,"Energie potentiel : %lg",Ep);
	  sprintf(Info2,"Conservation de l'energie : %7.5lg ",Ec+Ep);
	  sprintf(temp,"temp d'execution: %lg ",ttemp);
	  sprintf(thot,"temp Thot: %lg ",Thot);
	  ord[0]=((Ec+Ep-Et)/Et)/200;// renormalisation meilleur affichage
	  cpgpanl(1,1);
	  cpgsci(2);
	  cpgscf(1);
	  cpgsch(1.0);
	  if(ord[0]<-0.005) { ord[0]=-0.005; cpgsci(8); }
	  if(ord[0]> 0.005) { ord[0]= 0.005; cpgsci(8); }
	  cpgpt(1,abs,ord,1);
	  
	  cpgpanl(1,2);
	  cpgsci(1);
	  cpgscf(2);
	  cpgsch(1.5);
	  cpgtext(50,0.003,Info0);
	  cpgtext(50,0.0015,Info1);
	  cpgtext(50,0,Info2);
	  cpgtext(50,-0.0015,temp);
	  cpgtext(50,-0.003,thot);
	  
	  cpgpanl(2,1);
	  ord[0]=  ord[0] * 200;
	  if(ord[0]<-0.005) { ord[0]=-0.005; cpgsci(8); }
	  if(ord[0]> 0.005) { ord[0]= 0.005; cpgsci(8); }
	  cpgsci(2);
	  cpgscf(1);
	  cpgsch(1.0);
	  cpgpt(1,abs,ord,1);
	  
	  cpgpanl(2,2);
	  cpgsci(4);
	  cpgsch(0.1);
	  cpgpt(1000,x1,y1,style);
	  cpgpt(1000,x2,y2,style);
	  cpgpt(1000,x3,y3,style);
	  cpgsch(1.0);
	  cpgsci(0);
	  if(indice!=0)
	    cpgpt(NombreDeCorps,x,y,style);
	  cpgsci(7);
	  if(choix==1)
	    {
	      x[ind]=((float)Etoile[ind].Position[0]+1000)*0.5;
	      y[ind]= (float)Etoile[ind].Position[1]* 0.000005;
	    }
	  if(choix==0)
	    {
	      z[ind]= (float)Etoile[ind].Position[2];
	      x[ind]=0.5*(((float)Etoile[ind].Position[0]+1000)+((float)Etoile[ind].Position[2]*cos(PI/6)));
	      y[ind]= 0.000005 *((float)Etoile[ind].Position[1]+((float)Etoile[ind].Position[2]*sin(PI/6)));
	    }
	  cpgpt(NombreDeCorps,x,y,style);
	  cpgask(0);
	  RungeKutta(NombreDeCorps,PasDeTemp,VectPost,VectForce,KX,KV,Etoile);
	}
    }
  cpgask(1);
  cpgend();
  free(x);
  free(y);
  free(z);
  free(x1);
  free(y1);
  free(x2);
  free(y2);
  free(x3);
  free(y3);
}                                
//** mise en forme des données et gestion de l'affichage
//******************************************************
void Systeme ( int NbExecut, int NombreDeCorps, int Tcolor, double PasDeTemp,int Type,float max,Vecteur * VectPost, Vecteur * VectForce, RK * KX, RK * KV ,Corps * Etoile )
{   
  int indice,i,IndiceDuCorps;
  float * x,* y,* z,* m;
  double Thot,temp=0;
  x=(float *)calloc(NombreDeCorps,sizeof(float));
  y=(float *)calloc(NombreDeCorps,sizeof(float));
  z=(float *)calloc(NombreDeCorps,sizeof(float));
  m=(float *)calloc(NombreDeCorps,sizeof(float));
  printf("**  De quel corps desirez vous les informations : ");
  scanf("%d",&IndiceDuCorps);
  Thot=CalculduTempTho(NombreDeCorps,Etoile);
  cpgopen("/xw");                                        // ouverture
  FondSyst(max);  
  switch(Type)
    {
    case(0):
      for(i=0;i<NbExecut;i++)	
	for(indice=0;indice<NombreDeCorps;indice++)
	  {
	    TraceRefresh(Tcolor,NombreDeCorps,x,y,z);
	    x[indice]= (float)Etoile[indice].Position[0] ;
	    y[indice]= (float)Etoile[indice].Position[1] ;
	    z[indice]= (float)Etoile[indice].Position[2] ;
	    TraceRefresh(7,NombreDeCorps,x,y,z);
	    RungeKutta(NombreDeCorps,PasDeTemp,VectPost,VectForce,KX,KV,Etoile);
	    AffichInfo(IndiceDuCorps,NombreDeCorps,Thot,temp,Etoile,-max,max);
	    temp+=PasDeTemp;
	    cpgask(1);
	  }
      break;
    case(1):
      for(i=0;i<NbExecut;i++)
	for(indice=0;indice<NombreDeCorps;indice++)
	  {
	    TraceCercle(Tcolor,NombreDeCorps,x,y,z,m);
	    m[indice]= (float)Etoile[indice].Masse ;
	    x[indice]= (float)Etoile[indice].Position[0] ;
	    y[indice]= (float)Etoile[indice].Position[1] ;
	    z[indice]= (float)Etoile[indice].Position[2] ;
	    TraceCercle(7,NombreDeCorps,x,y,z,m);
	    RungeKutta(NombreDeCorps,PasDeTemp,VectPost,VectForce,KX,KV,Etoile);
	    AffichInfo(IndiceDuCorps,NombreDeCorps,Thot,temp,Etoile,-max,max);
	    temp+=PasDeTemp;
	    cpgask(1);
	  }
      break;
    }
  cpgend();
  free(x);
  free(y);
  free(z);
  free(m);
}
//** trace les points et enleve les anciens points
//************************************************
void TraceRefresh ( int Tcolor, int NombreDeCorps, float * x,float * y, float * z )
{
  int style=3;
  cpgscf(1.0);
  cpgsch(1.0);
  // trace les points et suivant la couleur affiche le point ,la trajectoire ou
  // les effaces ...
  cpgsci(Tcolor);
  cpgpanl(1,1);
  cpgpt(NombreDeCorps,x,y,style);
  cpgpanl(2,1);
  cpgpt(NombreDeCorps,x,z,style);
  cpgpanl(1,2);
  cpgpt(NombreDeCorps,z,y,style);
  if(Tcolor==7)
    sleep(0.1);
}
//** Trace des cercles
//********************
void TraceCercle (int Tcolor, int NombreDeCorps, float * x,float * y, float * z, float * m )
{
  int indice;
  float CentreX, CentreY, CentreZ, Rayon;
  cpgsfs(1);
  cpgsci(Tcolor);
  for(indice=0;indice<NombreDeCorps;indice++)
    {
      CentreX=x[indice];
      CentreY=y[indice];
      CentreZ=z[indice];
      Rayon = m[indice];
      if(Tcolor!=8)
	{
	  cpgpanl(1,1);
	  cpgcirc(CentreX,CentreY,Rayon);
	  cpgpanl(2,1);
	  cpgcirc(CentreX,CentreZ,Rayon);
	  cpgpanl(1,2);
	  cpgcirc(CentreZ,CentreY,Rayon);
	}
      if(Tcolor==8)
	{
	  cpgsci(0);
	  cpgpanl(1,1);
	  cpgcirc(CentreX,CentreY,Rayon);
	  cpgpanl(2,1);
	  cpgcirc(CentreX,CentreZ,Rayon);
	  cpgpanl(1,2);
	  cpgcirc(CentreZ,CentreY,Rayon);
	  cpgsci(8);
	  Rayon=Rayon*0.95;
	  cpgpanl(1,1);
	  cpgcirc(CentreX,CentreY,Rayon);
	  cpgpanl(2,1);
	  cpgcirc(CentreX,CentreZ,Rayon);
	  cpgpanl(1,2);
	  cpgcirc(CentreZ,CentreY,Rayon);
	}
    }
  if(Tcolor==7)
    sleep(0.1);      
}
//**  affichage des coordonnées d'un corps
//****************************************
void AffichInfo ( int IndiceDuCorps,int NombreDeCorps,double Tho,double PasDeTemps, Corps * Etoile, float xmin, float ymax  )
{
  char Masse[50],VitX[50],VitY[50],VitZ[50],DepX[50],DepY[50],DepZ[50],Temp[50],Eci[50],Epo[50],Eto[50],TTho[50];
  double Et,Ec,Ep;
  //  calcul et saisie prealable
  Ec=(float)EnergieCinetique(NombreDeCorps,Etoile);
  Ep=(float)EnergiePotentiel(NombreDeCorps,Etoile);
  Et=(float)EnergieCinetique(NombreDeCorps,Etoile)+(float)EnergiePotentiel(NombreDeCorps,Etoile);
  sprintf(Masse," Masse du corp %d : %03lg ",IndiceDuCorps,Etoile[IndiceDuCorps].Masse);
  sprintf(DepX, " Position sur x du corp %d : %06lg ",IndiceDuCorps,Etoile[IndiceDuCorps].Position[0]);
  sprintf(DepY, " Position sur y du corp %d : %06lg ",IndiceDuCorps,Etoile[IndiceDuCorps].Position[1]);
  sprintf(DepZ, " Position sur z du corp %d : %06lg ",IndiceDuCorps,Etoile[IndiceDuCorps].Position[2]);
  sprintf(VitX," Vitesse sur x du corp %d : %06lg ",IndiceDuCorps,Etoile[IndiceDuCorps].Vitesse[0]);
  sprintf(VitY," Vitesse sur y du corp %d : %06lg ",IndiceDuCorps,Etoile[IndiceDuCorps].Vitesse[1]);
  sprintf(VitZ," Vitesse sur z du corp %d : %06lg ",IndiceDuCorps,Etoile[IndiceDuCorps].Vitesse[2]);
  sprintf(Temp," Temps d'execution en M.Annee  : %04lg ",PasDeTemps);
  sprintf(TTho," Tho : %06lg ",Tho);
  sprintf(Eci," Energie cinetique : %06lg ",Ec);
  sprintf(Epo," Energie potentiel : %06lg ",Ep);
  sprintf(Eto," Energie totale : %06lg ",Et);
  //  affichage des informations
  cpgsci(1);         // | definition de la couleur du texte
  cpgscf(1);         // | definition de la police des textes 
  cpgsch(1.5);       // | definition de la taille de la police
  cpgpanl(2,2);      // | choix de la fenetre
  cpgeras();
  cpgtext(xmin,ymax*0.95,Masse);
  cpgtext(xmin,ymax*0.8,DepX);
  cpgtext(xmin,ymax*0.65,DepY);
  cpgtext(xmin,ymax*0.50,DepZ);
  cpgtext(xmin,ymax*0.35,VitX);
  cpgtext(xmin,ymax*0.20,VitY);
  cpgtext(xmin,ymax*0.05,VitZ);
  cpgtext(xmin,-ymax*0.25,Temp);
  cpgtext(xmin,-ymax*0.40,TTho);
  cpgtext(xmin,-ymax*0.55,Eci);
  cpgtext(xmin,-ymax*0.70,Epo);
  cpgtext(xmin,-ymax*0.85,Eto);
}
//** affichage du fond, partie stable
//***********************************
void FondSyst ( float max )
{
  float min=-max;
  // plusieurs fenetres
  cpgbeg(0,"/xw",2,2); 
  cpgscf(2);         // definition de la police des textes 
  cpgsch(1.8);       // definition de la taille de la police
  // fenetre 1-> haut gauche
  cpgpanl(2,2);
  cpgsci(4);
  cpgenv(min,max,min,max,0,-1); 
  cpgsci(1);
  cpglab("X","Y","   Projection sur X Y  ");
  // fenetre 2-> haut droite
  cpgpanl(1,1);                      
  cpgsci(4);  
  cpgenv(min,max,min,max,0,-1);
  cpgsci(1);
  cpglab("X","Z","   Projection sur X Z  ");
  // fenetre 3-> bas gauche
  cpgpanl(2,1);
  cpgsci(4);
  cpgenv(min,max,min,max,0,-1); 
  cpgsci(1);
  cpglab("Z","Y","   Projection sur Z Y  "); 
  // fenetre 4-> bas droite
  cpgpanl(1,2);
  cpgsci(0);
  cpgenv(min,max,min,max,0,-1);
  cpgsci(1);
  cpglab("","","  Information sur les corps : ");
}
//**  dessine le fond du trace energie
//************************************
void FondEner ( void ) 
{
  //**
  cpgpanl(2,2);
  cpgsci(1);
  cpgsch(1);
  cpgenv(0,1000,-1,1,0,2);
  cpgsci(4);
  cpgsch(1.5);
  cpglab("","Energie totale","Conservation de l'Energie Totale");
  //**
  cpgpanl(2,1);
  cpgsci(8);
  cpgsch(1);
  cpgenv(0,1000,0,1000,0,-1);
  cpgsch(1.5);
  cpglab("","","Information");
  //**
  cpgpanl(1,2);
  cpgsci(8);
  cpgsch(1);
  cpgenv(-1000,1000,-1000,1000,-1,0);
  cpgsch(1.5);
  cpglab("X","Y","Z");
  //**
  cpgpanl(1,1);
  cpgsci(1);
  cpgsch(1);
  cpgenv(0,1000,-0.005,0.005,0,2);
  cpgsci(4);
  cpgsch(1.5);
  cpglab("","Energie totale","Conservation de l'Energie Totale");
}
//** affichage du fond, partie stable
//***********************************
void FondOpt ( float max )
{
  float min=-max;
  // plusieurs fenetres
  cpgbeg(0,"/xw",2,2); 
  cpgscf(2);         // definition de la police des textes 
  cpgsch(1.8);       // definition de la taille de la police
  // fenetre 1-> haut gauche
  cpgpanl(2,2);
  cpgsci(4);
  cpgenv(min,max,min,max,0,0); 
  cpgsci(1);
  cpglab("X","Y","   Projection sur X Y  ");
  // fenetre 2-> haut droite
  cpgpanl(1,1);                      
  cpgsci(4);  
  cpgenv(min,max,min,max,0,0);
  cpgsci(1);
  cpglab("X","Z","   Projection sur X Z  ");
  // fenetre 3-> bas gauche
  cpgpanl(2,1);
  cpgsci(4);
  cpgenv(min,max,min,max,0,0); 
  cpgsci(1);
  cpglab("Z","Y","   Projection sur Z Y  "); 

  // fenetre 3-> bas gauche
  cpgpanl(1,2);
  cpgsci(4);
  cpgenv(min,max,min,max,0,0); 
  cpgsci(1);
  cpglab("X","Y","   Perspective Z fuit"); 
}
//**  tracer du systeme optimise
//******************************
void TraceOpt (  int NbExecut, int NombreDeCorps, int Tcolor, double PasDeTemp,float max,Vecteur * VectPost, Vecteur * VectForce, RK * KX, RK * KV ,Corps * Etoile ) 
{
  int indice,i,a,style=4;
  float * x,* y,* z, * x1, * y1, * x2, * y2, * x3, * y3, * xp, * yp;
  a=max;
  xp=(float *)calloc(NombreDeCorps,sizeof(float));
  yp=(float *)calloc(NombreDeCorps,sizeof(float));
  x=(float *)calloc(NombreDeCorps,sizeof(float));
  y=(float *)calloc(NombreDeCorps,sizeof(float));
  z=(float *)calloc(NombreDeCorps,sizeof(float));
  x1=(float *)calloc(a,sizeof(float));
  y1=(float *)calloc(a,sizeof(float));
  x2=(float *)calloc(a,sizeof(float));
  y2=(float *)calloc(a,sizeof(float));
  x3=(float *)calloc(a,sizeof(float));
  y3=(float *)calloc(a,sizeof(float));
  for(indice=0;indice<a;indice++)
    {
      x1[indice]= indice * cos(PI/6) - a;
      y1[indice]= indice * sin(PI/6) - a;

      x2[indice]= indice + a * ( cos(PI/6) -1 );
      y2[indice]= a * ( sin(PI/6) - 1 );

      x3[indice]= a * ( cos(PI/6) - 1 );
      y3[indice]= indice + a * ( sin(PI/6) -1 );
    }
  cpgopen("/xw");                                      
  FondOpt(max);
  cpgscf(1);
  cpgsch(0.5);  
  cpgpanl(2,2);
  cpgsci(4);
  cpgpt(max,x1,y1,3);
  cpgpt(max,x2,y2,3);
  cpgpt(max,x3,y3,3);
  cpgscf(1);
  cpgsch(1.0);  
  for(i=0;i<NbExecut;i++)	
    for(indice=0;indice<NombreDeCorps;indice++)
      {
	x[indice]=(float)Etoile[indice].Position[0] ;
	y[indice]=(float)Etoile[indice].Position[1] ;
	z[indice]=(float)Etoile[indice].Position[2] ;
	xp[indice]=x[indice]+(float)Etoile[indice].Position[2]*cos(PI/6);
	yp[indice]=y[indice]+(float)Etoile[indice].Position[2]*sin(PI/6);
	// ajoute
	cpgsci(7);
	cpgpanl(1,2);
	cpgpt(NombreDeCorps,z,y,style);
	cpgpanl(1,1);
	cpgpt(NombreDeCorps,x,y,style);
	cpgpanl(2,1);
	cpgpt(NombreDeCorps,x,z,style);
	cpgpanl(2,2);
	cpgpt(NombreDeCorps,xp,yp,style);

	RungeKutta(NombreDeCorps,PasDeTemp,VectPost,VectForce,KX,KV,Etoile);

	// enleve
	cpgsci(Tcolor);
	cpgpanl(1,2);
	cpgpt(NombreDeCorps,z,y,style);
	cpgpanl(1,1);
	cpgpt(NombreDeCorps,x,y,style);
	cpgpanl(2,1);
	cpgpt(NombreDeCorps,x,z,style);
	cpgpanl(2,2);
	cpgpt(NombreDeCorps,xp,yp,style);

	cpgask(0);
      }       
  cpgask(1);
  cpgend();
  free(x);
  free(y);
  free(z);
}
//**************************************************
//**  PROCEDURES ET FONCTIONS 
//**  D'AFFICHAGE ET D'INITIALISATION
//**************************************************
//**  Presentation
//*****************
void Presentation (void)
{
  printf("\n");
  printf("\n*****************************");
  printf("\n**  LES ETOILES MULTIPLES  **");
  printf("\n*****************************");
  printf("\n** Session : 2000 / 2001   **");
  printf("\n** Responsable : M Aristidi *");
  printf("\n*****************************");
}
//**  Menu
//********
void Menu (void)
{
  printf("\n**  MENU");
  printf("\n*****************************");
  printf("\n**  Etude d'un systeme ");
  printf("\n** 1. à deux corps");
  printf("\n** 2. à N corps ");
  printf("\n** 3. et tracer de l'énergie");
  printf("\n** 4. affichage optimise");
  printf("\n** 5. connu");
  printf("\n** 6. sortie du programme");
  printf("\n*****  Entrez votre choix : ");
}
//**  Saisie du nombre de corps
//*****************************
int SaisieNCorps ( void )
{ 
  int NombreDeCorps;
  printf("**  Entrez le nombre de corps que vous desirez : ");
  scanf("%d",&NombreDeCorps);
  return(NombreDeCorps);
}
//**  Initalisation des corps
//***************************
void InitCorps ( int NombreDeCorps , Corps * Etoile)
{
  int indice;
  printf("**  Il y a %d corps dans le systeme",NombreDeCorps);
  for (indice=0;indice<NombreDeCorps;indice++)
    {
      printf("\n**  Corps numero : %d",indice);
      printf("\n**  En masse solaire, en année lumière et million d'année");
      printf("\n**  A%02d.  Masse : ",indice);
      scanf("%lg",&Etoile[indice].Masse);
      printf("**  B%02d.  Position sur X : ",indice);
      scanf("%lg",&Etoile[indice].Position[0]);
      printf("**                 sur Y : ");
      scanf("%lg",&Etoile[indice].Position[1]);
      printf("**                 sur Z : ");
      scanf("%lg",&Etoile[indice].Position[2]);
      printf("**  C%02d.  Vitesse  sur X : ",indice);
      scanf("%lg",&Etoile[indice].Vitesse[0]);
      printf("**                 sur Y : ");
      scanf("%lg",&Etoile[indice].Vitesse[1]);
      printf("**                 sur Z : ");
      scanf("%lg",&Etoile[indice].Vitesse[2]);
      printf("**");
    }
  ConfirmInit(NombreDeCorps,Etoile);
}
//** Confirmation de l'initalisation
//**********************************
void  ConfirmInit ( int NombreDeCorps , Corps * Etoile)
{
  int indice,choix;
  for(indice=0;indice<NombreDeCorps;indice++)
    {
      printf("\n**");
      printf("\n**  Corps : %d, Masse : %lg, ",indice,Etoile[indice].Masse);
      printf("Position sur X : %lg",Etoile[indice].Position[0]);    
      printf(", sur Y : %lg",Etoile[indice].Position[1]);
      printf(", sur Z : %lg",Etoile[indice].Position[2]);
      printf("\n**  Vitesse sur X : %lg",Etoile[indice].Vitesse[0]);
      printf(", sur Y : %lg",Etoile[indice].Vitesse[1]);
      printf(", sur Z : %lg",Etoile[indice].Vitesse[2]);
    }
  printf("\n**  Si vous n'etes pas d'accord tapez '1' : ");
  scanf("%d",&choix);
  if(choix==1)
    { 
      printf("\n**");
      InitCorps(NombreDeCorps,Etoile);
    }
}
//**  Initialisation du pas de temp
//*********************************
double InitPasDeTemp ( void )
{
  double PasDeTemp=0.05;
  int choix;
  printf("**  Le pas de temps est de 0.05 M annee le configurer ( 1 ) : ");
  scanf("%d",&choix);
  if(choix==1)
    {
      printf("**  entrez le pas de temps : ");
      scanf("%lg",&PasDeTemp);
    }
  return(PasDeTemp);
}
//**  Initialisation des dimensions
//*********************************
float InitDim ( void )
{
  float dim=1000;
  int choix;
  printf("**  Les dimensions sont : -%06.02f - %06.02f Annee lumiere, les configurers ( 1 ) : ",dim,dim);
  scanf("%d",&choix);
  if(choix==1)
    {
      printf("**  entrez un extremum : ");
      scanf("%f",&dim);
    }
  return(dim);
}
//**  choix du type d initialisation
//**********************************
int TypeInit ( void )
{
  int choix;
  printf("**  Initialisation aleatoire ( 1 ) : ");
  scanf("%d",&choix);
  if(choix==1)
    return(1);
  return(0);
}
//**  choix du type de trace
//**************************
int TypeAffich ( void )
{
  int choix;
  printf("**  si vous désirez le tracé de la trajectoire tapez 1 : ");
  scanf("%d",&choix);
  if(choix==1)
    return(4);
  return(0);
}
//**  choix du type de corps
//**************************
int TypeCorps ( void )
{
  int choix;
  printf("**  si vous le désirez les corps peuvent etre des cercles\n**  Attention les masses > 10 seront peu exploitable, tapez 1 : ");
  scanf("%d",&choix);
  if(choix==1)
    return(1);
  return(0);
}
//**  definition du nombre d'iteration
//************************************
int Nbiteration ( void )
{
  int NbExecut;
  printf("**  Entrez le nombre d'iteration que vous desirez : ");
  scanf("%d",&NbExecut);
  return(NbExecut);
}
//**  Menu choix d un systeme connu
//*********************************
void MenuChoixSyst( Corps * Etoile,Vecteur * VectPost, Vecteur * VectForce, RK * KX, RK * KV  )
{
  int choix,NbExecut,NombreDeCorps,Tcolor,verif=0;
  float dim;
  double PasDeTemp;
  printf("**  Quels systeme desirez vous observer ");
  printf("\n**  0. deux corps instables");
  printf("\n**  1. deux corps stables");
  printf("\n**  Entrez votre choix : ");
  scanf("%d",&choix);
  Tcolor=TypeAffich();
  switch(choix)
    {
    case(0):
      NombreDeCorps=2;
      PasDeTemp=0.5;
      NbExecut=2000;
      break;
    case(1):
      NombreDeCorps=2;
      PasDeTemp=0.5;
      NbExecut=12000;
      break;
    default : 
      printf("**\n**  Hein hein hein essaye encore\n**");
      verif=1;
      break;
    }
  if(verif!=1)
    {
      Etoile=(Corps *)calloc(NombreDeCorps,sizeof(Corps));
      KX=(RK *)calloc(NombreDeCorps,sizeof(RK));
      KV=(RK *)calloc(NombreDeCorps,sizeof(RK));
      VectPost=(Vecteur *)calloc(1,sizeof(Vecteur));
      VectForce=(Vecteur *)calloc(1,sizeof(Vecteur));
            
      dim=InitSystConnu(choix,Etoile);  

      TraceOpt (NbExecut,NombreDeCorps,Tcolor,PasDeTemp,dim,VectPost,VectForce,KX,KV,Etoile); 
      
      free(Etoile);
      free(KX);
      free(KV);
      free(VectPost);
      free(VectForce);
    }
}
//**  initialisation du systeme connu
//***********************************
float InitSystConnu ( int choix, Corps * Etoile )
{
  switch(choix)
    {

    case(0):
      Etoile[0].Masse=20;
      Etoile[0].Position[0] =-100; Etoile[0].Position[1] =0; Etoile[0].Position[2] =-200;
      Etoile[0].Vitesse[0]=0.5;   Etoile[0].Vitesse[1]=0.2;   Etoile[0].Vitesse[2]=1;

      Etoile[1].Masse=20;
      Etoile[1].Position[0] =100; Etoile[1].Position[1] =0; Etoile[1].Position[2] =200;
      Etoile[1].Vitesse[0]=-0.5;   Etoile[1].Vitesse[1]=-0.2;   Etoile[1].Vitesse[2]=-1; 
      return(300);
      break;

    case(1):
      Etoile[0].Masse=80;
      Etoile[0].Position[0] =-100; Etoile[0].Position[1] =0; Etoile[0].Position[2] =-200;
      Etoile[0].Vitesse[0]=0.1;   Etoile[0].Vitesse[1]=1;   Etoile[0].Vitesse[2]=0.2;

      Etoile[1].Masse=80;
      Etoile[1].Position[0] =100; Etoile[1].Position[1] =0; Etoile[1].Position[2] =200;
      Etoile[1].Vitesse[0]=-0.1;   Etoile[1].Vitesse[1]=-1;   Etoile[1].Vitesse[2]=-0.2; 
      return(300);
      break;

    }
}
