#include "reco.h"
#include "helix/Helix.h"
#include EVTCLS_H //R2 distribution
#include "toolbox/Thrust.h"
#include "toolbox/FuncPtr.h"
#include "benergy/BeamEnergy.h"
#include "particle/utility.h" //This include have most functions for fitting data, setting MC data
#include "particle/combination.h"  //This include have functions for combining particles
#include "mdst/findKs.h"
#include "kid/atc_pid.h"
#include "toolbox/Thrust.h"
#include <toolbox/FoxWolfr.h>
 
#if defined(BELLE_NAMESPACE)


namespace Belle {
#endif

  using namespace std;

  void User_reco::hist_def( void )
  {
    //Get Address of BASF histogram manager
    extern BelleTupleManager* BASF_Histogram;
    btm=BASF_Histogram;
    //Make new ntuple (first parameter is key, second is string with all variables separated by space)
    std::string rs=t_par("DsP")+t_par("Ds17")+t_par("DsI")+t_par("DsII")+t_par("phi1")+t_par("phi2")+t_par("pi0")+t_par("K1")+t_par("K2")+t_par("K3")+t_par("K4")+"egam1 egam2 ip_x ip_y ip_z mctr cos0 thrr thrs r2 DsPart count hel";
    std::string rs2="m_pi0 p_pi0 egam";
    //t0=btm->ntuple("0","chi2 asym ppi0");
    t1=btm->ntuple("1",rs);
    //t2=btm->ntuple("2",rs2);
    //t3=btm->ntuple("3", "PID");
    //t4=btm->ntuple("4", "c2_phi");
    //t5=btm->ntuple("5", "c2_Ds");
    //t6=btm->ntuple("6", "c2_2317");
    };


  //===================================================================================

  void User_reco::event ( BelleEvent* evptr, int* status ) {
    int DEBUG=1;
    //FOR DEBUG ONLY
    if (IpProfile::usable() == 0) {
      if (DEBUG) {
	printf("[debugging routine] Ip Profile unusable -> skip current event.\n");
      }
      return;
    }

    if (IpProfile::b_life_smeared_usable() == 0) {
      if (DEBUG) {
	printf("[debugging routine] b_life_smeared unusable -> skip current event.\n");
      }
      return;
    }
    //status(skip histogram writing if 0)
    //hadron manager
    //Evtcls_hadron_info_Manager&  ehimgr = Evtcls_hadron_info_Manager::get_manager();
    //position and error of interaction point
    HepPoint3D ip_position = IpProfile::position(1);
    const HepSymMatrix& runIp_err = IpProfile::position_err(1);
    //printf("IP position acq.\n");
    //Evtcls_hadron_info_Manager::iterator iti = ehimgr.begin();
    //Evtcls_hadronic_flag_Manager::iterator eti = ehadfl.begin();

    Mdst_charged_Manager& mdstchgmgr =
    Mdst_charged_Manager::get_manager();//e,mu,pi,k
    Mdst_gamma_Manager & GamMgr = Mdst_gamma_Manager::get_manager();

    // double r2=0;
    //int ntrk=0;
    //double evis=0;
    //double Pz=0;
    //if (iti!=ehimgr.end()){
    //r2 = iti->R2();
    //ntrk = iti->Ntrk();
    //evis = iti->Evis();
    //Pz   = iti->Pz();
    //}
    getEventInfo( expNo, runNo, evtNo, isMC);//Get event info (exp, run, event, MC)
    //printf("Got evt info\n");
    //   double ecm = BeamEnergy::Ecm();
    //   double elec = BeamEnergy::E_HER();
    //   double posi = BeamEnergy::E_LER();

    // R2 Fox-folfram value, event shape cut
    Evtcls_hadron_info_Manager& HadMgr =   Evtcls_hadron_info_Manager::get_manager();
    Evtcls_hadron_info& Had = *HadMgr.begin();
    double r2 = Had.R2();

    //*************** Make particle lists ********************************

    p_list k_p, k_m, pi_p, pi_m, k_s, pi_0,gamma;

    makeKPi(k_p, k_m, pi_p, pi_m, 1);
    /*
    for(unsigned it=0; it<pi_p.size(); ++it){
      atc_pid sel_kpi(3,1,5,3,2);
      float kid=sel_kpi.prob(& (k_p[it].mdstCharged()));
      t2->column("PID", kid);
      float piid=sel_kpi.prob(& (pi_p[it].mdstCharged()));
      t2->column("PID", piid);
      t2->dumpData();}
    */

    //MakeKPi, last parameter define, if we select good charged particles
    //printf("Kpi asq., N(K_p)=%i, N(K_m)=%i,N(pi_p)=%i,N(pi_m)=%i\n",(int)k_p.size(),(int)k_m.size(),(int)pi_p.size(),(int)pi_m.size());

    //*************** Apply PID cuts **************************************

    double minProbPID = 0.5;  // preselected - probability
    double minProbPi = 0.9;

    //If each plist element is not within atc_pID.prob >= prob,
    //its element is removed from list.
    // Likelihood cuts
    /*
    for(int it=0; it<k_p.size(); ++it){
      atc_pid sel_kpi(3,1,5,3,2);
      float kid=sel_kpi.prob(& (k_p[it].mdstCharged()));
      if(kid>minProbPID){k_p.erase(k_p.begin()+it);--it;}
    }

    for(int it=0; it<k_m.size(); ++it){
      atc_pid sel_kpi(3,1,5,3,2);
      float kid=sel_kpi.prob(& (k_m[it].mdstCharged()));
      if(kid>minProbPID){k_m.erase(k_m.begin()+it);--it;}
    }


    for(int it=0; it<pi_p.size(); ++it){
      atc_pid sel_kpi(3,1,5,3,2);
      float piid=sel_kpi.prob(& (pi_p[it].mdstCharged()));
      if(piid>minProbPID){pi_p.erase(pi_p.begin()+it);--it;}
    }

    for(int it=0; it<pi_m.size(); ++it){
      atc_pid sel_kpi(3,1,5,3,2);
      float piid=sel_kpi.prob(& (pi_m[it].mdstCharged()));
      if(piid>minProbPID){pi_m.erase(pi_m.begin()+it);--it;}
    }
    */
    //First way to implement kaon PID
    //compare to pion
    //withKaonId(k_p, 0.5, 3, 1, 5, 3, 2);//K+ vs bg pi
    //withKaonId(k_m, 0.2, 3, 1, 5, 3, 2);//K- vs bg pi
    //compare to proton
    //withKaonId(k_p, minProbPID, 3, 1, 5, 3, 4);//K+ vs bg p
    //withKaonId(k_m, minProbPID, 3, 1, 5, 3, 4);//K- vs bg p

    //If each plist element is not within atc_pID.prob < prob,
    //its element is removed from plist.

    //First way to implement pion PID
    //compare to kaon
    //withPionId(pi_p, 0.9, 3, 1, 5, 3, 2);//pi+ vs bg k
    //withPionId(pi_m, 0.9, 3, 1, 5, 3, 2);//pi- vs bg k
    //compare to proton
    //     withPionId(pi_p, 1. - minProbPID, 3, 1, 5, 4, 2);//pi+ vs bg p
    //     withPionId(pi_m, 1. - minProbPID, 3, 1, 5, 4, 2);//pi- vs bg p


    //Second way to implement kaon-pion particle identification
    //(using external function)
    //pidKvspi(k_p,minProbPID);
    //pidKvspi(k_m,0.2);
    //pidpivsK(pi_p,minProbPi);
    //pidpivsK(pi_m,minProbPi);


    //Apply Dr and Dz cuts
    //double maxdr=100000;
    //double maxdz=100000;
    double maxdr=0.5;
    double maxdz=3;
    //printf("size1=%f", pi_m.size());
    withDrDz(pi_p,maxdr,maxdz);
    withDrDz(pi_m,maxdr,maxdz);
    withDrDz(k_p,maxdr,maxdz);
    withDrDz(k_m,maxdr,maxdz);
    //printf("size2=%f", pi_m.size());

    //*************** Make basic particles **********************************

    makePi0(pi_0);

    std::vector<double> chi2_vec;
    std::vector<double> asym_vec;
    std::vector<double> ppi0_vec;


    //Apply cuts on created po_0
    p_list new_pi0;
    for (p_list_it pi0_iterator=pi_0.begin(); pi0_iterator!=pi_0.end(); ++pi0_iterator) {

      const Mdst_pi0& pi0_mdst = pi0_iterator->mdstPi0();

      //double pi0_mass = pi0_mdst.mass();
      double pi0_chisq = pi0_mdst.chisq();

      Mdst_gamma& gam0 = pi0_mdst.gamma(0);
      Mdst_gamma& gam1 = pi0_mdst.gamma(1);

      double egam0 = sqrt(gam0.px()*gam0.px()+gam0.py()*gam0.py()+gam0.pz()*gam0.pz());
      HepLorentzVector gam0_vector (gam0.px(), gam0.py(), gam0.pz(), egam0);

      double egam1 = sqrt(gam1.px()*gam1.px()+gam1.py()*gam1.py()+gam1.pz()*gam1.pz());
    HepLorentzVector gam1_vector (gam1.px(), gam1.py(), gam1.pz(), egam1);

    HepLorentzVector pi0_vector = gam0_vector + gam1_vector;

    double pi0_theta0 = pi0_mdst.gamma(0).ecl().theta();
    double pi0_theta1 = pi0_mdst.gamma(1).ecl().theta();

    double thr0 = 100.0;
    double thr1 = 100.0;

    if ((pi0_theta0>0.2096)&&(pi0_theta0<0.5473)) thr0 = 0.1; // endcap(forward)
    if ((pi0_theta0>0.5619)&&(pi0_theta0<2.2466)) thr0 = 0.05; // barrel
    if ((pi0_theta0>2.2951)&&(pi0_theta0<2.7416)) thr0 = 0.1; // endcap(backward)

    if ((pi0_theta1>0.2096)&&(pi0_theta1<0.5473)) thr1 = 0.1;
    if ((pi0_theta1>0.5619)&&(pi0_theta1<2.2466)) thr1 = 0.05;
    if ((pi0_theta1>2.2951)&&(pi0_theta1<2.7416)) thr1 = 0.1;

    double pi0_chi2 = pi0_chisq;
    double pi0_asym = fabs((egam0-egam1)/(egam0+egam1));
    double ppi0 = pi0_vector.vect().mag(); // pi0 momentum

    chi2_vec.push_back(pi0_chi2);
    asym_vec.push_back(pi0_asym);
    ppi0_vec.push_back(ppi0);

    HepLorentzVector pi0_4 (pi0_mdst.px(), pi0_mdst.py(), pi0_mdst.pz(), pi0_mdst.energy());

    // this is the pi0 BEFORE FIT
    Particle pi0_unfitted (pi0_vector, Ptype("PI0"));

    // this is the pi0 AFTER FIT (NOMINAL mass), but fitted value is already in MDST
    Particle pi0_fitted (pi0_4, Ptype("PI0"));

    pi0_fitted.relation().append(pi0_iterator->child(0));
    pi0_fitted.relation().append(pi0_iterator->child(1));

    //if(pi0_chi2<200){new_pi0.push_back(pi0_fitted);}
    //if(pi0_asym>0.6){new_pi0.push_back(pi0_fitted);}
    //if(ppi0>0.1){new_pi0.push_back(pi0_fitted);}
    //new_pi0.push_back(pi0_fitted);

    double minegam=0.1;
    if(pi0_chi2<20 && pi0_asym<0.9 && ppi0>0.15 && egam0>minegam && egam1>minegam){
    //if(pi0_chi2<200 && pi0_asym<1.0 && ppi0>0.15 && egam0>minegam && egam1>minegam){ 
      new_pi0.push_back(pi0_fitted);
    //double pi0_m=(pi0_iterator->child(0).p()+pi0_iterator->child(1).p()).m();
    //t0->column("pi0_invm", pi0_m);
    //t0->column("chi2", pi0_chi2);
    //t0->column("asym",pi0_asym);
    //t0->column("ppi0",ppi0);
    //t0->dumpData();
    }
  }

//   printf("Parameters\n LENGTH1=%i\n LENGTH2=%i\n", (int)new_pi0.size(), (int)pi_0.size());
  pi_0.swap(new_pi0);

  setUserInfo(pi_0,1);
  for(unsigned i=0;i<pi_0.size();++i){
    setGammasError(pi_0[i], ip_position, runIp_err);
    //if(makeKvFit(pi_0[i],0) || pi_0[i].mass()<0.122 || pi_0[i].mass()>0.148  ||  makeKmvFit(pi_0[i])) {pi_0.erase(pi_0.begin()+i);--i;}
    //if(makeKvFit(pi_0[i],0)) {pi_0.erase(pi_0.begin()+i);--i;}
  }

  makeKvFit(pi_0,1);
  //maxCh2(pi_0, 6);
  minP(pi_0, 0.15);
  /*
  for (int it=0; it<pi_0.size();++it){
    t2->column("m_pi0", pi_0[it].mass());
    t2->column("p_pi0", pi_0[it].p().vect().mag());
    UserInfo & info =dynamic_cast<UserInfo&>(pi_0[it].userInfo());
    t2->column("c2_pi0", info.vchisq());
    t2->column("egam", pi_0[it].child(0).p().e());
    t2->column("egam", pi_0[it].child(1).p().e());
    t2->dumpData();
  }

  double minegam=0.1;
  for (int it=0; it<pi_0.size();++it){
    if(pi_0[it].child(0).p().e()<minegam || pi_0[it].child(1).p().e()<minegam){pi_0.erase(pi_0.begin()+it);--it;}
  }
  */
  withMassCut(pi_0, 0.125, 0.145);
  makeKmvFit(pi_0);

  /*
  makeGamma(gamma);
  for(p_list_it it=gamma.begin();it!=gamma.end();++it){
      setGammaError(*it, ip_position, runIp_err);
  }
  setUserInfo(gamma,1);
  */

//Get particles info
  if(isMC){
    setGenHepInfoF(pi_p);
    setGenHepInfoF(pi_m);
    setGenHepInfoF(k_p);
    setGenHepInfoF(k_m);
    setGenHepInfoKs(k_s);
    setGenHepInfoP(pi_0);
    setGenHepInfoG(gamma);
  }

//*************** Make phi particles ****************************

  double m_phi=Ptype(333).mass();
  double phi_w=0.01;
  p_list PHI1, PHI2, PHI3, PHI4;
  combination(PHI1, Ptype(333), k_m, k_p, phi_w*40);
  combination(PHI2, Ptype(333), k_m, k_p, phi_w*40);
  combination(PHI3, Ptype(333), k_m, k_p, phi_w*40);
  combination(PHI4, Ptype(333), k_m, k_p, phi_w*40);
  setUserInfo(PHI1,1);
  setUserInfo(PHI2,1);
  setUserInfo(PHI3,1);
  setUserInfo(PHI4,1);
  makeKvFit(PHI1,1);
  makeKvFit(PHI2,1);
  makeKvFit(PHI3,1);
  makeKvFit(PHI4,1);
  withMassCut(PHI1,m_phi-phi_w,m_phi+phi_w);
  withMassCut(PHI2,m_phi-phi_w,m_phi+phi_w);
  withMassCut(PHI3,m_phi-phi_w,m_phi+phi_w);
  withMassCut(PHI4,m_phi-phi_w,m_phi+phi_w);
  makeKmvFit(PHI1);
  makeKmvFit(PHI2);
  makeKmvFit(PHI3);
  makeKmvFit(PHI4);
if(isMC) {setGenHepInfoT(PHI1);}
if(isMC) {setGenHepInfoT(PHI2);}
if(isMC) {setGenHepInfoT(PHI3);}
if(isMC) {setGenHepInfoT(PHI4);}

/*
for (int it=0; it<PHI1.size();++it){
  UserInfo & info =dynamic_cast<UserInfo&>(PHI1[it].userInfo());
  t4->column("c2_phi", info.vchisq());
  t4->dumpData();
}
*/

//*************** Make Dsp2 particles ****************************

  //double m_Dsp2=Ptype(431).mass();
  double Dsp2_w=0.03;
  p_list DSP2;
  combination(DSP2, Ptype(431), PHI3, pi_p, Dsp2_w*10);
  setUserInfo(DSP2,1);
  makeKvFit(DSP2,1);
  withMassCut(DSP2, 1.9585, 1.9785);
  makeKmvFit(DSP2);
if(isMC) {setGenHepInfoT(DSP2);}

//*************** Make Dsm2 particles ****************************

  //double m_Dsm2=Ptype(-431).mass();
  double Dsm2_w=0.03;
  p_list DSM2;
  combination(DSM2, Ptype(-431), PHI4, pi_m, Dsm2_w*10);
  setUserInfo(DSM2,1);
  makeKvFit(DSM2,1);
  withMassCut(DSM2, 1.9585, 1.9785);
  makeKmvFit(DSM2);
if(isMC) {setGenHepInfoT(DSM2);}

/*
for (int it=0; it<DSM2.size();++it){
  UserInfo & info =dynamic_cast<UserInfo&>(DSM2[it].userInfo());
  t5->column("c2_Ds", info.vchisq());
  t5->dumpData();
}
*/

//Cuts on Dsp(m)2
//minCh2(DSP2, 0.001);
//minCh2(DSM2, 0.001);

//*************** Make Ds(2317)+ particles ****************************

  //double m_DspS=Ptype(10431).mass();
  double DspS_w=0.04;
  p_list DSPS;
  combination(DSPS, Ptype(10431), DSP2, pi_0, DspS_w*10);
  setUserInfo(DSPS,1);
  makeKvFit(DSPS,1);
  //withMassCut(DSPS, m_DspS-DspS_w, m_DspS+DspS_w);
  //makeKmvFit(DSPS);
if(isMC) {setGenHepInfoT(DSPS);}

//*************** Make Ds(2317)- particles ****************************

  //double m_DsmS=Ptype(-10431).mass();
  double DsmS_w=0.04;
  p_list DSMS;
  combination(DSMS, Ptype(-10431), DSM2, pi_0, DsmS_w*10);
  setUserInfo(DSMS,1);
  makeKvFit(DSMS,1);
  //withMassCut(DSMS, m_DsmS-DsmS_w, m_DsmS+DsmS_w);
  //makeKmvFit(DSMS);
if(isMC) {setGenHepInfoT(DSMS);}

/*
for (int it=0; it<DSMS.size();++it){
  UserInfo & info =dynamic_cast<UserInfo&>(DSMS[it].userInfo());
  t6->column("c2_2317", info.vchisq());
  t6->dumpData();
}
*/

// Cuts on Ds(2317)+(-)
/*
 for(unsigned i=0;i<DSPS.size();++i)
   if(DSPS[i].p().vect().mag()<0.3){DSPS.erase(DSPS.begin()+i);--i;}
 for(unsigned i=0;i<DSMS.size();++i)
   if(DSMS[i].p().vect().mag()<0.3){DSMS.erase(DSMS.begin()+i);--i;}
*/
//*************** Make Dsp1 particles ****************************

  //double m_Dsp1=Ptype(431).mass();
  double Dsp1_w=0.03;
  p_list DSP1;
  combination(DSP1, Ptype(431), PHI2, pi_p, Dsp1_w*10);
  setUserInfo(DSP1,1);
  makeKvFit(DSP1,1);
  //withMassCut(DSP1, 1.955, 1.979);
  //makeKmvFit(DSP1);
if(isMC) {setGenHepInfoT(DSP1);}

//*************** Make Dsm1 particles ****************************

  //double m_Dsm1=Ptype(-431).mass();
  double Dsm1_w=0.03;
  p_list DSM1;
  combination(DSM1, Ptype(-431), PHI1, pi_m, Dsm1_w*10);
  setUserInfo(DSM1,1);
  makeKvFit(DSM1,1);
  //withMassCut(DSM1, 1.955, 1.979);
  //makeKmvFit(DSM1);
if(isMC) {setGenHepInfoT(DSM1);}

//Cuts on Dsp(m)1
//minCh2(DSP1, 0.001);
//minCh2(DSM1, 0.001);
//*************** Make mother particle (Dsprt) ****************************

p_list DSPRT;
combination(DSPRT, Ptype("UPS4"), DSPS, DSM1);
setUserInfo(DSPRT, 1);
combination(DSPRT, Ptype("UPS4"), DSMS, DSP1);
setUserInfo(DSPRT, 2);
makeKvFit(DSPRT,1);
//withMassCut(DSPRT,m_Dsm1+m_DspS-0.5, m_Dsm1+m_DspS+0.5);
if(isMC) {setGenHepInfoT(DSPRT);}

*status=1;
if(DSPRT.empty()) *status=0;

//R2 cut
// for(int i=0;i<DSMS.size();++i)
//   if(r2>0.9){DSMS.erase(DSPS.begin()+i);--i;}

/*
 p_list newDSPRT;
 double minEnGam=0;
 double minEnPi0=0;
 double minPPi0=0;
 double minPstVal=0;
 minEnGam=0.1;
 minPPi0=0.2;

 for(int it=0; it<DSPRT.size();it++){
   UserInfo &info_Dsprt=dynamic_cast<UserInfo&>(DSPRT[it].userInfo());
   int Dsprt_cn=info_Dsprt.channel();
   if((minEn(DSPRT[it].child(0).child(1).child(0), minEnGam)==1) || (minEn(DSPRT[it].child(0).child(1).child(1), minEnGam)==1))
    continue;
   if(minP(DSPRT[it].child(0).child(1), minPPi0)==1)
     continue;
   else
     newDSPRT.push_back(DSPRT[it]);
 }
 */

 /*
  int bst=0;
  double del_min=0;
  if(DSPRT.size()!=0) del_min = (DSPRT[0].child(0).momentum().p().mag()-2.317)*(DSPRT[0].child(0).momentum().p().mag()-2.317) + (DSPRT[0].child(1).momentum().p().mag()-1.969)*(DSPRT[0].child(1).momentum().p().mag()-1.969);
  // if(DSPRT.size()!=0) cout << del_min << endl;
  //  if(DSPRT.size()!=0) cout << "NMP=" << DSPRT.size() << endl;
   if(DSPRT.size()>1){
     for (int i=1; i<DSPRT.size(); ++i){
       if ((DSPRT[i].child(0).momentum().p().mag()-2.317)*(DSPRT[i].child(0).momentum().p().mag()-2.317) +\
           (DSPRT[i].child(1).momentum().p().mag()-1.969)*(DSPRT[i].child(1).momentum().p().mag()-1.969)<del_min){
 	del_min=(DSPRT[i].child(0).momentum().p().mag()-2.317)*(DSPRT[i].child(0).momentum().p().mag()-2.317) +\
 	  (DSPRT[i].child(1).momentum().p().mag()-1.969)*(DSPRT[i].child(1).momentum().p().mag()-1.969);
 	bst=i;}
       //      cout << bst << "~" << abs(DSPRT[i].momentum().mass()-4.287) << "  ";
       }}
   // if(DSPRT.size()!=0)  cout << endl;
   //if(DSPRT.size()!=0) cout << "BestDel=" << del_min << endl;

   if(DSPRT.size()!=0) {
     p_list temp;
     temp.push_back(DSPRT[bst]);
     DSPRT.clear();
     DSPRT.push_back(temp[0]);
   }

   //if(DSPRT.size()!=0) cout << abs(DSPRT[0].momentum().mass()-4.287) << endl;
   //if(DSPRT.size()!=0) cout << "newSIZE=" << DSPRT.size() << endl;
   */

 Hep3Vector  CMBoost   = - BeamEnergy::CMBoost();
 for(unsigned i=0;i<DSPRT.size();++i){
   std::vector<Hep3Vector> sVic;
   for (unsigned j=0; j < DSPRT[i].relation().nFinalStateParticles();++j){
     HepLorentzVector MPc(DSPRT[i].relation().finalStateParticle(j).p());
     MPc.boost(CMBoost);
     Hep3Vector MP3(MPc.vect());
     sVic.push_back(MP3);
   }
   Vector3 thrust_signal = thrust(sVic.begin(), sVic.end(), SelfFunc(Vector3()));

   std::vector<Hep3Vector> oVic;
   Vector4 v4_o(0);
   Ptype ptype(321);
   for(std::vector<Mdst_charged>::iterator it=mdstchgmgr.begin();
     it!=mdstchgmgr.end(); ++it){
     //Mdst_charged& chg = *it;
     //Mdst_trk& Trk = chg.trk();
     Particle  temp(*it, ptype);
     HepLorentzVector MP1(temp.momentum().p());
     MP1.boost(CMBoost);
     Hep3Vector temp1(MP1.vect());
     if(it->get_ID() == DSPRT[i].child(0).child(0).child(0).child(0).mdstCharged().get_ID()) continue;
     if(it->get_ID() == DSPRT[i].child(0).child(0).child(0).child(1).mdstCharged().get_ID()) continue;
     if(it->get_ID() == DSPRT[i].child(0).child(1).mdstCharged().get_ID()) continue;
     if(it->get_ID() == DSPRT[i].child(1).child(1).mdstCharged().get_ID()) continue;
     if(it->get_ID() == DSPRT[i].child(1).child(0).child(0).mdstCharged().get_ID()) continue;
     if(it->get_ID() == DSPRT[i].child(1).child(0).child(1).mdstCharged().get_ID()) continue;
     oVic.push_back(temp1);
   }
   ptype = Ptype(-321);
   for(std::vector<Mdst_charged>::iterator it=mdstchgmgr.begin();
       it!=mdstchgmgr.end(); ++it){
     //Mdst_charged& chg = *it;
     //Mdst_trk& Trk = chg.trk();
     Particle  temp(*it, ptype);
     HepLorentzVector MP1(temp.momentum().p());
     MP1.boost(CMBoost);
     Hep3Vector temp1(MP1.vect());
     if(it->get_ID() == DSPRT[i].child(0).child(0).child(0).child(0).mdstCharged().get_ID()) continue;
     if(it->get_ID() == DSPRT[i].child(0).child(0).child(0).child(1).mdstCharged().get_ID()) continue;
     if(it->get_ID() == DSPRT[i].child(0).child(1).mdstCharged().get_ID()) continue;
     if(it->get_ID() == DSPRT[i].child(1).child(1).mdstCharged().get_ID()) continue;
     if(it->get_ID() == DSPRT[i].child(1).child(0).child(0).mdstCharged().get_ID()) continue;
     if(it->get_ID() == DSPRT[i].child(1).child(0).child(1).mdstCharged().get_ID()) continue;
     oVic.push_back(temp1);
   }
   ptype = Ptype(211);
   for(std::vector<Mdst_charged>::iterator it=mdstchgmgr.begin();
       it!=mdstchgmgr.end(); ++it){
     //Mdst_charged& chg = *it;
     //Mdst_trk& Trk = chg.trk();
     Particle  temp(*it, ptype);
     HepLorentzVector MP1(temp.momentum().p());
     MP1.boost(CMBoost);
     Hep3Vector temp1(MP1.vect());
     if(it->get_ID() == DSPRT[i].child(0).child(0).child(0).child(0).mdstCharged().get_ID()) continue;
     if(it->get_ID() == DSPRT[i].child(0).child(0).child(0).child(1).mdstCharged().get_ID()) continue;
     if(it->get_ID() == DSPRT[i].child(0).child(1).mdstCharged().get_ID()) continue;
     if(it->get_ID() == DSPRT[i].child(1).child(1).mdstCharged().get_ID()) continue;
     if(it->get_ID() == DSPRT[i].child(1).child(0).child(0).mdstCharged().get_ID()) continue;
     if(it->get_ID() == DSPRT[i].child(1).child(0).child(1).mdstCharged().get_ID()) continue;
     oVic.push_back(temp1);
   }
   ptype = Ptype(-211);
   for(std::vector<Mdst_charged>::iterator it=mdstchgmgr.begin();
       it!=mdstchgmgr.end(); ++it){
     //Mdst_charged& chg = *it;
     //Mdst_trk& Trk = chg.trk();
     Particle  temp(*it, ptype);
     HepLorentzVector MP1(temp.momentum().p());
     MP1.boost(CMBoost);
     Hep3Vector temp1(MP1.vect());
     if(it->get_ID() == DSPRT[i].child(0).child(0).child(0).child(0).mdstCharged().get_ID()) continue;
     if(it->get_ID() == DSPRT[i].child(0).child(0).child(0).child(1).mdstCharged().get_ID()) continue;
     if(it->get_ID() == DSPRT[i].child(0).child(1).mdstCharged().get_ID()) continue;
     if(it->get_ID() == DSPRT[i].child(1).child(1).mdstCharged().get_ID()) continue;
     if(it->get_ID() == DSPRT[i].child(1).child(0).child(0).mdstCharged().get_ID()) continue;
     if(it->get_ID() == DSPRT[i].child(1).child(0).child(1).mdstCharged().get_ID()) continue;
     oVic.push_back(temp1);
   }

   // Noe the neutral gamma
   for (std::vector<Mdst_gamma>::const_iterator ig = GamMgr.begin();
	ig != GamMgr.end(); ig++)
     { //Mdst_ecl& ecl_gam = ig -> ecl();
       Particle gama(*ig);
       HepLorentzVector MP2(gama.p());
       MP2.boost(CMBoost);
       Hep3Vector temp2(MP2.vect());
       if( ig->get_ID()==DSPRT[i].child(0).child(1).child(0).mdstGamma().get_ID()) continue;
       if( ig->get_ID()==DSPRT[i].child(0).child(1).child(1).mdstGamma().get_ID()) continue;
       oVic.push_back(temp2);
     }
   Vector3 thrust_other = thrust(oVic.begin(), oVic.end(), SelfFunc(Vector3()));
   double costhr = thrust_signal.unit().dot(thrust_other.unit());
   double m_thr_roe = thrust_other.mag();
   double m_thr_sig = thrust_signal.mag();
   t1->column("cos0", costhr);
   t1->column("thrr", m_thr_roe);
   t1->column("thrs", m_thr_sig);
   t1->column("r2", r2);
   //   if (abs(Get.MiD(DSPRT[i].child(0).child(0),1))==10431 || abs(Get.MiD(DSPRT[i].child(1),1))==920){
   //     t1->column("mctr",1);}
   //   else
   //     t1->column("mctr",-1);
      t_fill_ups(DSPRT[i], t1);
 }

btm->clearAllData();
return;

}

//===================================================================================

//*************** Set info for particle (WTF section) ***************

void setUserInfo(Particle &p, unsigned ch)
{
  if(! &p.userInfo()){p.userInfo(UserInfo(ch));}
}

//===================================================================================

//*************** Set info for list (WTF section) *****************
void setUserInfo(p_list &p, unsigned ch)
{ for(p_list_it it=p.begin();it!=p.end();++it) setUserInfo(*it,ch); }

// UserInfo Class

UserInfo::UserInfo()
  : m_vchisq(-1.),
    m_dgf(-1.),
    m_vmass(0.),
    m_channel(0)
  //  m_index(all_infos.size())
{
//all_infos.push_back(this);
}

UserInfo::UserInfo(unsigned ch)
  : m_vchisq(-1.),
    m_dgf(-1),
    m_vmass(0.),
    m_channel(ch)
   // m_index(all_infos.size())
{
//all_infos.push_back(this);
}

UserInfo::UserInfo(const UserInfo &x)
  : m_vchisq(x.m_vchisq),
    m_dgf(x.m_dgf),
    m_vmass(x.m_vmass),
    m_channel(x.m_channel)
  //  m_index(all_infos.size())
{
//all_infos.push_back(this);
}
UserInfo::~UserInfo()
{
}

UserInfo*
UserInfo::clone(void) const
{
  UserInfo *x = new UserInfo(*this);
  //all_infos.push_back(x);
  return x;
}

UserInfo &
UserInfo::operator = (const UserInfo &x)
{
  m_vchisq = x.m_vchisq;
  m_dgf = x.m_dgf;
  m_vmass  = x.m_vmass;
  m_channel= x.m_channel;
  return *this;
}

//===================================================================================

//*************** Make the fits *************************************

//Get ID from hepevt (for MC)
int IDhep(Particle &part) {
    if(!part.genHepevt()) return 0;
    return part.genHepevt().idhep();
}

//===================================================================================

//Do vertex fit
int makeKvFit(Particle &p, int flg) {
      if(p.nChildren()<2) return 1;
//if no userinfo, create one
  if(! &p.userInfo() ) setUserInfo(p);
//save old mass in vmass variable
  dynamic_cast<UserInfo&>(p.userInfo()).vmass(p.mass());
  double v_chi=-1;
  double prob_chi=-1;
//   if( dynamic_cast<UserInfo&>(p.userInfo()).vchisq() != -1. ) return;
  //create fitter
  kvertexfitter kv;
  //int nfit=0;
//add tracks
  for(unsigned i=0; i<p.nChildren(); i++)
      addTrack2fit(kv,p.child(i));
  if(flg && IpProfile::usable()) addBeam2fit(kv);
  //make fit and set chi2, make mother particle from fitter, if fit was successful (err!=0)
  int err=kv.fit();
  if(!err){
    v_chi=kv.chisq();
    prob_chi=kv.dgf();
    dynamic_cast<UserInfo&>(p.userInfo()).vchisq(v_chi);
    dynamic_cast<UserInfo&>(p.userInfo()).prob_vchisq(prob_chi);
    makeMother(kv,p);
    return 0;
  }else return 1;
//   dynamic_cast<UserInfo&>(p.userInfo()).vchisq(v_chi);
}

//===================================================================================

//Make fit of list of particles
void makeKvFit(vector<Particle> &plist, int flg) {
  for(unsigned i=0;i<plist.size();++i) if(makeKvFit(plist[i], flg)) {plist.erase(plist.begin()+i);--i;}
}

//===================================================================================

int makeKmvFit(Particle &p) {
//if no userinfo, create one
      if(p.nChildren()<2) return 1;

  if(! &p.userInfo() ) setUserInfo(p);
  //save old mass in vmass variable
  dynamic_cast<UserInfo&>(p.userInfo()).vmass(p.mass());
  double v_chi=-1;
  double prob_chi=-1;
  //if( dynamic_cast<UserInfo&>(p.userInfo()).vchisq() != -1. ) return;
  //create fitter
  kmassvertexfitter kmv;
  kmv.invariantMass(p.pType().mass());
  //int nfit=0;
  //add tracks
  for(unsigned i=0; i<p.nChildren(); i++)
      addTrack2fit(kmv,p.child(i));
  //make fit and set chi2, make mother particle from fitter, if fit was successful (err!=0)
  int err=kmv.fit();
  if(!err){
    v_chi=kmv.chisq();
    prob_chi=kmv.dgf();
    makeMother(kmv,p);
    dynamic_cast<UserInfo&>(p.userInfo()).vchisq(v_chi);
    dynamic_cast<UserInfo&>(p.userInfo()).prob_vchisq(prob_chi);
    return 0;
  } else return 1;

}

//===================================================================================

//*************** Filling ntuples ***********************************

void makeKmvFit(vector<Particle> &plist) {
  for(unsigned i=0;i<plist.size();++i) if (makeKmvFit(plist[i])) {plist.erase(plist.begin()+i);--i;}
}
    std::string User_reco::t_par(std::string tag){
        const int npars=15;
        const int maxl=4;
        const char vals[npars][maxl]={"m","px","py","pz","p","e","c2","dgf","cn","mc","mv","vx","vy","vz","pst"};
        string res;
        for(int i=0;i<npars;++i){
            res+=vals[i];
            res+="_";
            res+=tag;
            res+=" ";
        }
        return res;
    }
//===================================================================================

    void User_reco::t_fill(Particle &p,BelleTuple *t,std::string tg){
        const int npars=15;
        const int maxl=4;
        const char vals[npars][maxl]={"m","px","py","pz","p","e","c2","dgf","cn","mc","mv","vx","vy","vz","pst"};
        std::string tag="_"+tg;
        t->column(vals[0]+tag,p.momentum().mass());
        t->column(vals[1]+tag,p.momentum().p().px());
        t->column(vals[2]+tag,p.momentum().p().py());
        t->column(vals[3]+tag,p.momentum().p().pz());
        t->column(vals[4]+tag,p.momentum().p().vect().mag());
        t->column(vals[5]+tag,p.momentum().p().e());

   if(& p.userInfo()){
	    UserInfo & info =dynamic_cast<UserInfo&>(p.userInfo());
            t->column(vals[6]+tag,info.vchisq());
            t->column(vals[7]+tag,info.prob_vchisq());
            t->column(vals[8]+tag,info.channel());
            t->column(vals[10]+tag,info.vmass());
        }
   else {
            t->column(vals[6]+tag,-1);
            t->column(vals[7]+tag,-1);
	    t->column(vals[8]+tag,p.pType().lund());
            t->column(vals[10]+tag,p.momentum().mass());
        }
   //t->column(vals[9]+tag,IDhep(p));
        if(& p.momentum().decayVertex()){
            t->column(vals[11]+tag,p.momentum().decayVertex().x());
            t->column(vals[12]+tag,p.momentum().decayVertex().y());
            t->column(vals[13]+tag,p.momentum().decayVertex().z());
        }
        t->column(vals[14]+tag,pStar(p,BeamEnergy::E_HER(),BeamEnergy::E_LER(),BeamEnergy::Cross_angle()).vect().mag());

    }

    //===================================================================================

    void User_reco::t_fill_ups(Particle &p,BelleTuple *t){
        if(IpProfile::usable()){
	  const HepPoint3D & ip= IpProfile::position(1);
	  t->column("ip_x",ip.x());
	  t->column("ip_y",ip.y());
	  t->column("ip_z",ip.z());
	}
	t_fill(p,t,"DsP");
	t_fill(p.child(0),t,"Ds17");
	t_fill(p.child(0).child(0),t,"DsII");
	t_fill(p.child(0).child(0).child(0),t,"phi1");
	t_fill(p.child(0).child(0).child(0).child(0),t,"K1");
	t_fill(p.child(0).child(0).child(0).child(1),t,"K2");
	t_fill(p.child(0).child(1),t,"pi0");
	t->column("egam1",p.child(0).child(1).child(0).p().e());
	t->column("egam2",p.child(0).child(1).child(1).p().e());
	t_fill(p.child(1),t,"DsI");
	t_fill(p.child(1).child(0),t,"phi2");
	t_fill(p.child(1).child(0).child(0),t,"K3");
	t_fill(p.child(1).child(0).child(0),t,"K4");
	// Helicity angle
	HepLorentzVector ds_part=p.child(0).child(0).p();
	HepLorentzVector phi_part=p.child(0).child(0).child(0).p();
	HepLorentzVector k_part=p.child(0).child(0).child(0).child(0).p();
	Hep3Vector bst=-phi_part.boostVector();
	phi_part.boost(bst);
	k_part.boost(bst);
	double hel=(ds_part.vect().x()*k_part.vect().x()+ds_part.vect().y()*k_part.vect().y()+ds_part.vect().z()*k_part.vect().z())/ds_part.vect().mag()/k_part.vect().mag();
	t->column("hel", hel);
        t->dumpData();
    }

//===================================================================================

int pidKvspi(Particle &p, float prob){
    atc_pid sel_kpi(3,1,5,3,2);
    float kid=sel_kpi.prob(& (p.mdstCharged()));
    if(kid<prob) return 1;
    return 0;
}

//===================================================================================

void pidKvspi(p_list &plist, float prob){
    for(unsigned i=0;i<plist.size();++i)
        if(pidKvspi(plist[i],prob)) {plist.erase(plist.begin()+i);--i;}
}

//===================================================================================

int pidpivsK(Particle &p, float prob){
    atc_pid sel_kpi(3,1,5,3,2);
    float kid=sel_kpi.prob(& (p.mdstCharged()));
    if(kid>prob) return 1;
    return 0;
}

//===================================================================================

void pidpivsK(p_list &plist, float prob){
    for(unsigned i=0;i<plist.size();++i)
        if(pidpivsK(plist[i],prob)) {plist.erase(plist.begin()+i);--i;}
}

  HepVector imp_par_ip(const Mdst_charged& charged, int hyp){
  const HepPoint3D pivot(charged.trk().mhyp(hyp).pivot_x(),
                           charged.trk().mhyp(hyp).pivot_y(),
                           charged.trk().mhyp(hyp).pivot_z());
    HepVector  a(5);
    a[0] = charged.trk().mhyp(hyp).helix(0);
    a[1] = charged.trk().mhyp(hyp).helix(1);
    a[2] = charged.trk().mhyp(hyp).helix(2);
    a[3] = charged.trk().mhyp(hyp).helix(3);
    a[4] = charged.trk().mhyp(hyp).helix(4);
    HepSymMatrix Ea(5,0);
    Ea[0][0] = charged.trk().mhyp(hyp).error(0);
    Ea[1][0] = charged.trk().mhyp(hyp).error(1);
    Ea[1][1] = charged.trk().mhyp(hyp).error(2);
    Ea[2][0] = charged.trk().mhyp(hyp).error(3);
    Ea[2][1] = charged.trk().mhyp(hyp).error(4);
    Ea[2][2] = charged.trk().mhyp(hyp).error(5);
    Ea[3][0] = charged.trk().mhyp(hyp).error(6);
    Ea[3][1] = charged.trk().mhyp(hyp).error(7);
    Ea[3][2] = charged.trk().mhyp(hyp).error(8);
    Ea[3][3] = charged.trk().mhyp(hyp).error(9);
    Ea[4][0] = charged.trk().mhyp(hyp).error(10);
    Ea[4][1] = charged.trk().mhyp(hyp).error(11);
    Ea[4][2] = charged.trk().mhyp(hyp).error(12);
    Ea[4][3] = charged.trk().mhyp(hyp).error(13);
    Ea[4][4] = charged.trk().mhyp(hyp).error(14);
    Helix helix(pivot, a, Ea);
    const Hep3Vector&   IP = IpProfile::position();
    if (IP.mag())
      helix.pivot(IP);
    return helix.a();
    }

 //===================================================================================

 int withDrDz(Particle &p, double maxdr, double maxdz){
        int typ=-1;
        int pcode=abs(p.pType().lund());
        if(pcode==11) typ=0;
        if(pcode==13) typ=1;
        if(pcode==211)typ=2;
        if(pcode==321)typ=3;
        if(pcode==2212)typ=4;
        if(typ<0) return 0;
        HepVector I_p=imp_par_ip(p.mdstCharged(),typ);
        double dz=I_p[3];
        double dr=I_p[0];
        if(abs(dz)>maxdz || abs(dr)>maxdr) return 1;
        return 0;
    }

//===================================================================================

void withDrDz(p_list &plist, double maxdr, double maxdz){
        for(unsigned i=0;i<plist.size();++i) if(withDrDz(plist[i], maxdr,maxdz)) {plist.erase(plist.begin()+i);--i;}
    }

//===================================================================================

int minEn(Particle &p, double minE){
//     if(! & p.mdstGamma()) return 0;
//     Hep3Vector v(gam.px(),gam.py(),gam.pz());
    return p.p().e()<minE;
}

//===================================================================================

void minEn(p_list &plist, double minE){
    for(unsigned i=0;i<plist.size();++i)
        if(minEn(plist[i],minE)){plist.erase(plist.begin()+i);--i;}
}

//===================================================================================

int minP(Particle &p, double minp){
  //     if(! & p.mdstGamma()) return 0;
  //     Hep3Vector v(gam.px(),gam.py(),gam.pz());
  return p.p().vect().mag()<minp;
}

//===================================================================================

void minP(p_list &plist, double minp){
  for(unsigned i=0;i<plist.size();++i)
    if(minP(plist[i],minp)){plist.erase(plist.begin()+i);--i;}
}

//===================================================================================
int minCh2(Particle &p, double minch2){
  UserInfo & info =dynamic_cast<UserInfo&>(p.userInfo());
  return info.vchisq()<minch2;
}

//===================================================================================

void minCh2(p_list &plist, double minch2){
  for(unsigned i=0;i<plist.size();++i)
    if(minCh2(plist[i],minch2)){plist.erase(plist.begin()+i);--i;}
}

//===================================================================================

void combination5p(p_list &plist,const Ptype &ptype, Particle &p1, Particle &p2, Particle &p3, Particle &p4, Particle &p5){
    if(checkSame(p1,p2)||checkSame(p1,p3)||checkSame(p1,p4)||checkSame(p1,p5)|| checkSame(p2,p3) ||checkSame(p2,p4)||checkSame(p2,p5)||checkSame(p3,p4)||checkSame(p3,p5)||checkSame(p4,p5)) return;
    HepLorentzVector vect=p1.momentum().p()+p2.momentum().p()+p3.momentum().p()+p4.momentum().p()+p5.momentum().p();
    Particle p(vect,ptype);
    p.relation().append(p1);
    p.relation().append(p2);
    p.relation().append(p3);
    p.relation().append(p4);
    p.relation().append(p5);
    plist.push_back(p);
}

void combination5p(p_list &plist, const Ptype &ptype, p_list &p1, p_list &p2, p_list &p3, p_list &p4, p_list &p5){
    for(unsigned long  i=0;i<p1.size();++i)
        for(unsigned long  j=0;j<p2.size();++j)
            for(unsigned long  k=0;k<p3.size();++k)
                for(unsigned long  l=0;l<p4.size();++l)
                    for(unsigned long m=0;m<p5.size();++m)
                        combination5p(plist,ptype,p1[i],p2[j],p3[k],p4[l],p5[m]);
}

//===================================================================================

void combination5p(p_list &plist, const Ptype &ptype, Particle &p1, Particle &p2, Particle &p3, Particle &p4, Particle &p5, const double &larg, const double &rarg){
    if(checkSame(p1,p2)||checkSame(p1,p3)||checkSame(p1,p4)||checkSame(p1,p5)|| checkSame(p2,p3) ||checkSame(p2,p4)||checkSame(p2,p5)||checkSame(p3,p4)||checkSame(p3,p5)||checkSame(p4,p5)) return;
    HepLorentzVector vect=p1.momentum().p()+p2.momentum().p()+p3.momentum().p()+p4.momentum().p()+p5.momentum().p();
    double mass_dif=vect.mag()-ptype.mass();
    if(mass_dif<-larg || mass_dif>rarg) return;
    Particle p(vect,ptype);
    p.relation().append(p1);
    p.relation().append(p2);
    p.relation().append(p3);
    p.relation().append(p4);
    p.relation().append(p5);
    plist.push_back(p);
}

void combination5p(p_list &plist, const Ptype &ptype, p_list &p1, p_list &p2, p_list &p3, p_list &p4, p_list &p5, const double &larg, const double &rarg){
    for(unsigned long i=0;i<p1.size();++i)
        for(unsigned long  j=0;j<p2.size();++j)
            for(unsigned long  k=0;k<p3.size();++k)
                for(unsigned long  l=0;l<p4.size();++l)
                    for(unsigned long  m=0;m<p5.size();++m)
                        combination5p(plist,ptype,p1[i],p2[j],p3[k],p4[l],p5[m],larg,rarg);
}

//===================================================================================

void combination5p(p_list &plist, const Ptype &ptype, Particle &p1, Particle &p2, Particle &p3, Particle &p4, Particle &p5, const double &marg){
    combination5p(plist,ptype,p1,p2,p3,p4,p5,marg,marg);
}

void combination5p(p_list &plist, const Ptype &ptype, p_list &p1, p_list &p2, p_list &p3, p_list &p4, p_list &p5, const double &marg){
    combination5p(plist,ptype,p1,p2,p3,p4,p5,marg,marg);
}

//===================================================================================

void combination6p(p_list &plist,const Ptype &ptype, Particle &p1, Particle &p2, Particle &p3, Particle &p4, Particle &p5, Particle &p6){
    if(checkSame(p1,p2)||checkSame(p1,p3)||checkSame(p1,p4)||checkSame(p1,p5)|| checkSame(p1,p6)|| checkSame(p2,p3) ||checkSame(p2,p4)||checkSame(p2,p5)||checkSame(p2,p6)||checkSame(p3,p4)||checkSame(p3,p5) ||checkSame(p3,p6) ||checkSame(p4,p5)||checkSame(p4,p6) ||checkSame(p5,p6)) return;
    HepLorentzVector vect=p1.momentum().p()+p2.momentum().p()+p3.momentum().p()+p4.momentum().p()+p5.momentum().p()+p6.momentum().p();
    Particle p(vect,ptype);
    p.relation().append(p1);
    p.relation().append(p2);
    p.relation().append(p3);
    p.relation().append(p4);
    p.relation().append(p5);
    p.relation().append(p6);
    plist.push_back(p);
}

void combination6p(p_list &plist, const Ptype &ptype, p_list &p1, p_list &p2, p_list &p3, p_list &p4, p_list &p5, p_list &p6){
    for(unsigned long  i=0;i<p1.size();++i)
        for(unsigned long  j=0;j<p2.size();++j)
            for(unsigned long  k=0;k<p3.size();++k)
                for(unsigned long  l=0;l<p4.size();++l)
                    for(unsigned long m=0;m<p5.size();++m)
                        for(unsigned long n=0;n<p6.size();++n)
                            combination6p(plist,ptype,p1[i],p2[j],p3[k],p4[l],p5[m],p6[n]);
}

//===================================================================================

void combination6p(p_list &plist, const Ptype &ptype, Particle &p1, Particle &p2, Particle &p3, Particle &p4, Particle &p5, Particle &p6, const double &larg, const double &rarg){
     if(checkSame(p1,p2)||checkSame(p1,p3)||checkSame(p1,p4)||checkSame(p1,p5)|| checkSame(p1,p6)|| checkSame(p2,p3) ||checkSame(p2,p4)||checkSame(p2,p5)||checkSame(p2,p6)||checkSame(p3,p4)||checkSame(p3,p5) ||checkSame(p3,p6) ||checkSame(p4,p5)||checkSame(p4,p6) ||checkSame(p5,p6)) return;
    HepLorentzVector vect=p1.momentum().p()+p2.momentum().p()+p3.momentum().p()+p4.momentum().p()+p5.momentum().p()+p6.momentum().p();
    Particle p(vect,ptype);
    p.relation().append(p1);
    p.relation().append(p2);
    p.relation().append(p3);
    p.relation().append(p4);
    p.relation().append(p5);
    p.relation().append(p6);
    plist.push_back(p);
}

void combination6p(p_list &plist, const Ptype &ptype, p_list &p1, p_list &p2, p_list &p3, p_list &p4, p_list &p5, p_list &p6, const double &larg, const double &rarg){
    for(unsigned long  i=0;i<p1.size();++i)
        for(unsigned long  j=0;j<p2.size();++j)
            for(unsigned long  k=0;k<p3.size();++k)
                for(unsigned long  l=0;l<p4.size();++l)
                    for(unsigned long m=0;m<p5.size();++m)
                        for(unsigned long n=0;n<p6.size();++n)
                            combination6p(plist,ptype,p1[i],p2[j],p3[k],p4[l],p5[m],p6[n]);
}

//===================================================================================

void combination6p(p_list &plist, const Ptype &ptype, Particle &p1, Particle &p2, Particle &p3, Particle &p4, Particle &p5, Particle &p6, const double &marg){
    combination6p(plist,ptype,p1,p2,p3,p4,p5,p6,marg,marg);
}

void combination6p(p_list &plist, const Ptype &ptype, p_list &p1, p_list &p2, p_list &p3, p_list &p4, p_list &p5, p_list &p6, const double &marg){
    combination6p(plist,ptype,p1,p2,p3,p4,p5,p6,marg,marg);
}

//===================================================================================

#if defined(BELLE_NAMESPACE)
}
#endif
