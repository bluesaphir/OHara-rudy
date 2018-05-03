//atr
#include <stdio.h>

#define NTABLE 1000
#define ETMIN -100.
#define ETMAX 100.

typedef double REAL;

struct Table { int n1; float frac; } T1;
struct Tlocal { REAL to; REAL uo; REAL tn; REAL un; };

struct State {
    double v;
    double nai;
    double nass;
    double ki;
    double kss;
    double cai;
    double cass;
    double cansr;
    double cajsr;
    double m;
    double hf;
    double hs;
    double j;
    double hsp;
    double jp;
    double mL;
    double hL;
    double hLp;
    double a;
    double iF;
    double iS;
    double ap;
    double iFp;
    double iSp;
    double d;
    double ff;
    double fs;
    double tff;
    double tfs;
    double fcaf;
    double fcas;
    double jca;
    double nca;
    double ffp;
    double fcafp;
    double xrf;
    double xrs;
    double xs1;
    double xs2;
    double xura;
    double xuri;
    double xk1;
    double Jrelnp;
    double Jrelp;
    double CaMKt;
    double fp;
    double f;
    double dss;
    double fss;
} *CurrentState, *NextState;

int iter;
REAL ee;

// Coefficients  for ventricle 10 23 12+/ 4 17 16 +/5 23 16+/ 11 2 16+/2 19 17+/3 11 17/ 9 30 16
//'SCN5A', 'KCNA4', 'CACNA1C', 'KCNH2', 'KCNQ1', 'KCNJ2', 'SLC8A1', 'ATP1A1' ,'KCNJ11' ,'ATP2B4' ,'ATP2A2', 'CALM1', 'RYR2', 'CAMK2D'
	
//ref

double const c_gnal = 1;
double const c_gto =1; 
double const c_pca = 1;
double const c_gkr = 1;
double const c_gks =1;
double const c_gk1 =1;
double const c_gncx =1;
double const c_pnak =1;
double const c_gkb = 1;
double const c_gkatp = 1;
double const c_gpca = 1;
double const c_jrelinf =1;
double const c_jup = 1;
double const c_CMDN = 1;
double const c_arel = 1;
double const c_CaMKo =1; 
	
//means diff	
	
double const c_m_gnal =	0.736; // SCN5A
double const c_m_gto =	0.287; // ventr SLOW endocard for KCNA4 
double const c_m_pca = 	1.106; // CACNA1C 
double const c_m_gkr = 	0.896; // KCNH2 
double const c_m_gks = 	1.174; // KCNQ1 
double const c_m_gk1 = 	0.151; // KCNJ2 
double const c_m_gncx = 0.983; // SLC8A1 
double const c_m_pnak = 0.607; // ATP1A1 
double const c_m_gkatp =1.129; // KCNJ11 
double const c_m_gpca = 0.521; // ATP2B4 
double const c_m_jup = 	1.695; // ATP2A2
double const c_m_CMDN = 0.857; // CALM1
 
double const c_m_arel =	0.937; // RYR2 
double const c_m_CaMKo =0.733; // (rudy CAMK2B) Calcium/Calmodulin-Dependent Protein Kinase (CaMK)
								// CAMK2D
	
						
// 10 23 12						

/*4 17 16 						
*/
						
/*5 23 16 
 */
						
/*11_2_16
						
 */																
									
/*2_19_17
						 */
											
/*3_11_17
						 */						 	
/*9_30_16

*/	
						 	




//iso coefficients
double const c_INa_dV_act=1.046;
double const c_INa_dV_inact=0.543;
double const c_ICaL_dV_act=0.959;
double const c_ICaL_dV_inact=1.228;
double const c_Knai=1.026;
double const c_GKb_iso=0.617;
double const c_SERCA=0.995;
double const c_GNa_ISO=1.012;
double const c_Tnl=1.155;
double const c_IKs_ISO=0.958;//1.652;
double const c_txs_1_iso=0.459;//1.987;
double const c_a_rel=1.051;
double const c_tau_rel=0.993;

//constants
double const nao=140.0;//extracellular sodium in mM
double const cao=1.8;//extracellular calcium in mM
double const ko=5.4;//extracellular potassium in mM

//buffer paramaters
double const BSRmax=0.047;
double const KmBSR=0.00087;
double const BSLmax=1.124;
double const KmBSL=0.0087;
double const cmdnmax=0.05;
double const kmcmdn=0.00238;
double const trpnmax=0.07;
double const kmtrpn=0.0005;
double const csqnmax=10.0;
double const kmcsqn=0.8;

//CaMK paramaters
double const aCaMK=0.05;
double const bCaMK=0.00068;
double const CaMKo=c_m_CaMKo * c_CaMKo*0.05;
double const KmCaM=0.0015;
double const KmCaMK=0.15;

//physical constants
double const R=8314.0;
double const T=310.0;
double const F=96485.0;

//cell geometry
double const L=0.01;
double const rad=0.0011;
double const vcell=1000*3.14*rad*rad*L;
double const Ageo=2*3.14*rad*rad+2*3.14*rad*L;
double const Acap=2*Ageo;
double const vmyo=0.68*vcell;
double const vmito=0.26*vcell;
double const vsr=0.06*vcell;
double const vnsr=0.0552*vcell;
double const vjsr=0.0048*vcell;
double const vss=0.02*vcell;

//introduce varaibles for reversal potentials, currents, fluxes, and CaMK
double ENa,EK,EKs;
double dss,fss;
double INa,Ito,ICaL,ICaNa,ICaK,IKr,IKur,IKs,IK1,INaCa_i,INaCa_ss,INaCa,INaK,IKb,INab,IpCa,ICab,Ist,IKatp;
double Jrel,Jup,Jtr,Jdiff,JdiffNa,JdiffK,Jleak;
double CaMKa,CaMKb;

//introduce APD, timing, and counting parameters
int APD_flag=0;
double APD;
double t_vdot_max;
double vrest;
double vo = 20;
int step =10; //
int basis = -60; // change this parameter if ypu want different height in the step experiment
//int v1 = -60;
double v2 = -80;
double dt=0.005; //0.005
double t0=0.;
double t=0;
double dto;
double vdot_old;
double vdot=0;
double vdot_max;
int p=1;
int n_stim=0;
int n = 0;
int Count=0;


const double amp = - 80;//-250; //-250;//250;//stimulus amplitude in uA/uF
const double start = 0;//start time of the stimulus, relative to each beat
const double duration = 0.5;//duration of teh stimulus in ms 0.5 //COMMENT both for step and basis
//const double duration = 300 // both for step and basis
const double duration1 = 3000; // basis length,  both for step and basis



double g_gap_junc=5.0;
