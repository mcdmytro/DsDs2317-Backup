//___________________________________________________________________
//       UTILITY written again to make my code much simpler 
//         & clean...    
//   AUTHOR --- Vishal ... ^_^
//___________________________________________________________________

#include"belle.h"
#include "myutl.h"
#include <strstream>
#include <iomanip>
#include BELLETDF_H
#include MDST_H
#include HEPEVT_H
#include "particle/Ptype.h"
#include <vector>

#if defined(BELLE_NAMESPACE)
namespace Belle{
#endif

  //________________
  //________________         MY VERTEX FIT
  
  unsigned myVfit::VXfit( Particle & p, double & confLevel, int debug) {
        kvertexfitter kv;
    for(unsigned j=0;j<p.relation().nChildren();++j){
      
      Particle child = p.relation().child(j);
      kv.addTrack( child.momentum().p(),
		   child.momentum().x(),
		   child.momentum().dpx(),
		   child.pType().charge(),
		   child.pType().mass() );
    }

    unsigned merr = kv.fit();
    ferr=kv.fit();
    if( merr ) {
      return 0;
    }
    confLevel = kv.cl();
   
    return myVfit::makevMother( kv, p );
  }

  unsigned myVfit::makevMother( kvertexfitter &kv, Particle &mother )
{
  unsigned n = kv.tracks();
  kmakemother kmm;
  for(unsigned i=0;i<n;++i){
    kmm.addTrack(kv.momentum(i),
                 kv.position(i),
                 kv.error(i),
                 mother.relation().child(i).pType().charge());
    kmm.errVertexTrack(kv.errVertexTrack(i));
    for(unsigned j=i+1;j<n;++j){
      kmm.correlation(kv.correlation(i,j));
      HepMatrix tmp1(kv.correlation(i,j));
    }
  }
  
  kmm.vertex(kv.vertex());
  kmm.errVertex(kv.errVertex());
  unsigned err = kmm.make();
  if(err != 0)return 0;
  mother.momentum().momentumPosition(kmm.momentum(),
                                     kmm.position(),
                                     kmm.error());
  mother.momentum().decayVertex(kv.vertex(), kv.errVertex());
  return 1;
}

 
  

  
  int myVfit::f_err() {
    return ferr; }
  
  HepPoint3D myVfit::VX(){
    return kv.vertex();}
  
  HepSymMatrix myVfit::erVX(){
    return kv.errVertex();}
  
  double myVfit::chisq(){
    return kv.chisq();}
  
  int myVfit::dgf(){
    return kv.dgf();}
  
  //____________________________________________________________________
  //                         MY MASS FIT
  //____________________________________________________________________



  

 
  unsigned myMfit::VMfit( Particle & p, double & confLevel, int flag, double massy ) {
    kmassfitter km;
    if(flag==1){
      km.invariantMass(massy);
    } else {
      km.invariantMass(p.pType().mass()); }
  
  /*if( debug > 0 ) cout << "doKmFit: Parent: " << p.pType().name()
		       << ", No. Children: "
		       << p.relation().nChildren() <<endl;*/
  
  for(unsigned j=0;j<p.relation().nChildren();++j){
    
    Particle child = p.relation().child(j);
    km.addTrack( child.momentum().p(),
		 child.momentum().x(),
		 child.momentum().dpx(),
		 child.pType().charge(),
		 child.pType().mass() );
  }
  km.notDecayPoint();
  
  unsigned merr = km.fit();
  if( merr ) {
    return 0;
  }
  confLevel = km.cl();
  
  return myMfit::makemMother( km, p );
}
  
  double myMfit::chisq(){
    return km.chisq();}

  double myMfit::m_err(){
    return merr;}
  double myMfit::mcl(){
    return km.cl();}
  

unsigned  myMfit::makemMother( kmassfitter &kv, Particle &mother ) {
  
  unsigned n = kv.tracks();
  kmakemother kmm;
  
  
  for( unsigned i = 0; i<n; ++i ) {
    kmm.addTrack( kv.momentum( i ),
		  kv.position( i ),
		  kv.error(    i ),
		  mother.relation().child( i ).pType().charge(),
		  KF_AFTER_FIT );
    kmm.errVertexTrack( kv.errVertexTrack( i ));
    
    for( unsigned j=i+1; j<n; ++j ) {
      kmm.correlation( kv.correlation(i,j) );
      HepMatrix tmp1( kv.correlation(i,j) );
    }
    
  }
  
  
  kmm.vertex( kv.vertex() );
  kmm.errVertex( kv.errVertex() );
  unsigned err = kmm.make();
  if(err != 0) return 0;
  
  mother.momentum().momentumPosition(kmm.momentum(),
				     kmm.position(),
				     kmm.error());
  mother.momentum().decayVertex( kv.vertex(), kv.errVertex() );
  return 1;
  
}

  //----------------------------------------------------------
  //                     PARTY
  //__________________________________________________________
  double party::getpT(Particle &pa){
    return sqrt (pa.p().px() * pa.p().px() + pa.p().py()*pa.p().py());
    //        return pa.p().vect().mag() * sin(pa.p().theta());
  }

  double party::getpT(HepLorentzVector &pa){
        return sqrt(pa.px()*pa.px() + pa.py()*pa.py());
  }
  //for genhep
  HepLorentzVector party::mon4(Gen_hepevt & a){
    HepLorentzVector p(a.E(), a.PX(), a.PY(), a.PZ());
    return p;
  }
  //for genhep
  HepVector party::vex(Gen_hepevt & a){
    HepVector p(3);
    p(0)=a.VX();
    p(1)=a.VY();
    p(2)=a.VZ();
    return p;
  }
  //for genhep
  double party::genpTt(Gen_hepevt & a){
    HepLorentzVector p(a.E(), a.PX(), a.PY(), a.PZ());
    //  return p.vect().mag() * sin ( acos( a.VZ()/ (a.VX() *a.VX() + a.VY()*a.VY() + a.VZ()*a.VZ())  )) ;
    return p.vect().mag() * sin( atan (a.VY()/a.VX()));
  }
  
  //for genhep
  double party::genpT(Gen_hepevt & a){
    // HepLorentzVector p(a.E(), a.PX(), a.PY(), a.PZ());
    return sqrt (a.PX()*a.PX() + a.PY()*a.PY());
  }

  // __________________________________________________________________
  //                        MOTHER IDs
  //___________________________________________________________________
  int motherid::oiD(Particle &pa) {
    Gen_hepevt_Manager& hepevt_mag = Gen_hepevt_Manager::get_manager();
    testid =  pa.relation().genHepevt().get_ID();
    if (abs(testid) != 0) {
    
      return pa.relation().genHepevt().idhep(); }
    else return 9876789;
}

   
  int motherid::MiD( Particle &pa, int gen) {
    Gen_hepevt_Manager& hepevt_mag = Gen_hepevt_Manager::get_manager();
    
    
      testid = pa.relation().genHepevt().get_ID();
      if (abs(testid) !=0) {    	ma_id = pa.relation().genHepevt().idhep(); 
       } else { ma_id = 9876789; }
      for ( int p=1; p<= gen; p++) {
	if ( abs(testid) !=0 && abs(ma_id) !=300553 && abs(ma_id) != 10022  ) {
	for(std::vector<Gen_hepevt>::iterator it = hepevt_mag.begin();
	    it != hepevt_mag.end(); ++it){
	  Gen_hepevt& particle = *it;
	  if (particle.get_ID()== testid) {
	    ma_id = particle.mother().idhep();
	    testid =particle.mother().get_ID();
	  }
	}
	}
      }
      return ma_id;
  }



  double motherid::helPi0(Particle &PI0){
    double e_pi0 = PI0.momentum().p().e();
    double p_pi0 = PI0.momentum().p().vect().mag();
    double e_1 = PI0.mdstPi0().gamma(0).ecl().energy();
    double e_2 = PI0.mdstPi0().gamma(1).ecl().energy();
    return (e_pi0/p_pi0)*((e_1-e_2)/(e_1+e_2))    ;
    

  }


  int motherid::pi0trumaa(Particle &PI0){
    Mdst_sim_ecl_Manager & sime = Mdst_sim_ecl_Manager::get_manager();
    int ma1 =0; int ma2 =0;
    int id1 =0; int id2 =0;
    int eclid1 = PI0.mdstPi0().gamma(0).ecl().get_ID();
    int eclid2 = PI0.mdstPi0().gamma(1).ecl().get_ID();
    for(std::vector<Mdst_sim_ecl>::iterator i=sime.begin();
	i!=sime.end(); i++){
      Mdst_sim_ecl& simec = *i;
      if( simec.ecl().get_ID()== eclid1) {id1= simec.hepevt().get_ID();
	if (id1 !=0 && simec.hepevt().mother().get_ID()!=0) {
	  ma1= simec.hepevt().mother().idhep();} }
      if( simec.ecl().get_ID()== eclid2) {id2= simec.hepevt().get_ID();
	if (id2 !=0&& simec.hepevt().mother().get_ID()!=0) {
	  ma2= simec.hepevt().mother().idhep();} }
    }
    
    if( (abs(ma1) ==111  && abs(ma2)!=111 ) || (abs(ma2) ==111 && abs(ma1) !=111))
      {return 1; }
    else if(abs(ma1)==111 && abs(ma2)==111) {return 2; }
    else {return 0; }
  }
  
  
  int motherid::pi0Tinfo( Particle &PI0, int nkid, int noTri) {
    Mdst_ecl_aux_Manager & EcL_AuX = Mdst_ecl_aux_Manager::get_manager();
    int eclid = PI0.mdstPi0().gamma(nkid).ecl().get_ID();
    double nhit  = EcL_AuX[eclid -1].property(6);
    
    float propy1 = EcL_AuX[eclid-1].property(1);
    float propy2 = EcL_AuX[eclid-1].property(2);
    float propy3 = EcL_AuX[eclid-1].property(3);
    
    int pop1 = *(int*)&propy1;
    int pop2 = *(int*)&propy2;
    int pop3 = *(int*)&propy3;
    
    int tcid1 = (pop1 & 0x3ff);
    int tsid1 = (pop1 >> 10) & 0x1;
    int tcid2 = (pop1 >> 11) & 0x3ff;
    int tsid2 = (pop1 >> 21) & 0x1;

    int tdccount1 = pop2 & 0xffff;
    int elinkcls1 = (pop2 >> 16) & 0x7f;
    int elinktc1 = (pop2 >> 23) & 0x7f;
    int hecls1 = (pop2 >> 30) & 0x1;
    int hetc1 = (pop2 >> 31) & 0x1;

    int tdccount2 = pop3 & 0xffff;
    int elinkcls2 = (pop3 >> 16) & 0x7f;
    int elinktc2 = (pop3 >> 23) & 0x7f;
    int hecls2 = (pop3 >> 30) & 0x1;
    int hetc2 = (pop3 >> 31) & 0x1;

    if(noTri==1){
      return tdccount1;}
    else if(noTri ==2){
      return tdccount2;}
    else{
      std::cout << "FATAL ERROR MAA says\t PLEASE PUT CORRECT WHICH TRIGGER you want\n correct input 1 or 2 " << std::endl;
      return 0; }
    // if( tcid1!=0 && tcid1==tcid2 && tsid1==tsid2){

    //return nhit;
  }
 

  //mother of the gamma daughter of pi0                                         
  //Particle Pi0                                                                
  //  nkid is kid 0, 1                                                          
  //gen is the generation of mother                                             
  // 1 mother                                                                   
  // 2 mother                                                                   
  //3 mother                                                                    
  // you can go till infnite mother as permitted or code allow                  
  // make good logic inorder to avoid Segmentation violation                    

  int motherid::pi0Gmaa(Particle &PI0, int nkid, int gen){
    Mdst_sim_ecl_Manager & sime = Mdst_sim_ecl_Manager::get_manager();
    Gen_hepevt_Manager& hepevt_mag = Gen_hepevt_Manager::get_manager();
    int eclid = PI0.mdstPi0().gamma(nkid).ecl().get_ID();
    int maaid=0; int maa =0;
    for(std::vector<Mdst_sim_ecl>::iterator i=sime.begin();
        i!=sime.end(); i++){
      Mdst_sim_ecl& simec = *i;
      if( simec.ecl().get_ID()== eclid) {
        maaid = simec.hepevt().mother().get_ID();
        if (maaid !=0) { maa = simec.hepevt().mother().idhep(); }
      }
    }

    if(maaid ==0) return 9876789;

    if (gen >1) {
      for ( int p=0; p< gen-1; p++) {
        if ( abs(maaid) !=0 && abs(maa) !=300553 && abs(maa) != 10022) {
          for(std::vector<Gen_hepevt>::iterator it = hepevt_mag.begin();
              it != hepevt_mag.end(); ++it){
            Gen_hepevt& particle = *it;
            if (particle.get_ID()== maaid) {
              maa = particle.mother().idhep();
              maaid = particle.mother().get_ID();
              //      for(int ss=0; ss<gen;  ss++){  // for fun                 
	      //                std::cout << "::::::::"; }                      
              //  std::cout << endl;                                            
            }
          }
        } else{ return 9876789;}      }
    }

    return maa;
  }





 
  int motherid::GMiD(int eclid, int gen){
    Mdst_sim_ecl_Manager & sime = Mdst_sim_ecl_Manager::get_manager();
    Gen_hepevt_Manager& hepevt_mag = Gen_hepevt_Manager::get_manager();
    int maaid=0; int maa =0;
    for(std::vector<Mdst_sim_ecl>::iterator i=sime.begin();
	i!=sime.end(); i++){
      Mdst_sim_ecl& simec = *i;
      if( simec.ecl().get_ID()== eclid) {
	maaid = simec.hepevt().mother().get_ID();
	if (maaid !=0) { maa = simec.hepevt().mother().idhep(); }
      }
    }
    
    if(maaid ==0) return 9876789;	
    
    if (gen >1) {
      for ( int p=0; p< gen-1; p++) {
	if ( abs(maaid) !=0 && abs(maa) !=300553 && abs(maa) != 10022) {
	  for(std::vector<Gen_hepevt>::iterator it = hepevt_mag.begin();
	      it != hepevt_mag.end(); ++it){
	    Gen_hepevt& particle = *it;
	    if (particle.get_ID()== maaid) {
	      maa = particle.mother().idhep();
	      maaid = particle.mother().get_ID();
	      //      for(int ss=0; ss<gen;  ss++){  // for fun
	      //		std::cout << "::::::::"; }
	      //  std::cout << endl;
	    }
	  }
	} else{ return 9876789;}      }
    }
    
    return maa;
  }



  int motherid::Gisthep( int eclid){
    Mdst_ecl_aux_Manager &EcL_AuX = Mdst_ecl_aux_Manager::get_manager();
    Mdst_sim_ecl_Manager & sime = Mdst_sim_ecl_Manager::get_manager();
    Gen_hepevt_Manager& hepevt_mag = Gen_hepevt_Manager::get_manager();
    int ownid=0; int isthep =0;
    for(std::vector<Mdst_sim_ecl>::iterator i=sime.begin();
        i!=sime.end(); i++){
      Mdst_sim_ecl& simec = *i;
      if( simec.ecl().get_ID()== eclid) {
        ownid = simec.hepevt().get_ID();
        if (ownid !=0) { isthep = simec.hepevt().isthep(); }
      }
    }
    return isthep;
  }

  int motherid::GTinfo( int eclid, int noTri) {
    Mdst_ecl_aux_Manager & EcL_AuX = Mdst_ecl_aux_Manager::get_manager();

    double nhit  = EcL_AuX[eclid -1].property(6);

    float propy1 = EcL_AuX[eclid-1].property(1);
    float propy2 = EcL_AuX[eclid-1].property(2);
    float propy3 = EcL_AuX[eclid-1].property(3);
    
    int pop1 = *(int*)&propy1;
    int pop2 = *(int*)&propy2;
    int pop3 = *(int*)&propy3;
    
    int tcid1 = (pop1 & 0x3ff);
    int tsid1 = (pop1 >> 10) & 0x1;
    int tcid2 = (pop1 >> 11) & 0x3ff;
    int tsid2 = (pop1 >> 21) & 0x1;

    int tdccount1 = pop2 & 0xffff;
    int elinkcls1 = (pop2 >> 16) & 0x7f;
    int elinktc1 = (pop2 >> 23) & 0x7f;
    int hecls1 = (pop2 >> 30) & 0x1;
    int hetc1 = (pop2 >> 31) & 0x1;

    int tdccount2 = pop3 & 0xffff;
    int elinkcls2 = (pop3 >> 16) & 0x7f;
    int elinktc2 = (pop3 >> 23) & 0x7f;
    int hecls2 = (pop3 >> 30) & 0x1;
    int hetc2 = (pop3 >> 31) & 0x1;

    if(noTri==1){
      return tdccount1;}
    else if(noTri ==2){
      return tdccount2;}
    else{
      std::cout << "FATAL ERROR MAA says\t PLEASE PUT CORRECT WHICH TRIGGER you want\n correct input 1 or 2 " << std::endl;
      return 0; }
    // if( tcid1!=0 && tcid1==tcid2 && tsid1==tsid2){

    //return nhit;
  }










  //mother of the pi+ of KS
  //Particle Ks0
  //  nkid is kid 0, 1
  //gen is the generation of mother
  // 1 mother
  // 2 mother
  // 3 mother
  // you can go till infnite mother as permitted or code allow
  // make good logic inorder to avoid Segmentation violation
  
  int motherid::KsPimaa(Particle &K_S, int nkid, int gen){
    Mdst_sim_trk_Manager & sime = Mdst_sim_trk_Manager::get_manager();
    Gen_hepevt_Manager& hepevt_mag = Gen_hepevt_Manager::get_manager();
    int eclid = K_S.mdstVee2().chgd(nkid).trk().get_ID();
    int maaid=0; int maa =0;
    for(std::vector<Mdst_sim_trk>::iterator i=sime.begin();
        i!=sime.end(); i++){
      Mdst_sim_trk& simec = *i;
      if( simec.trk().get_ID()== eclid) {
        maaid = simec.hepevt().mother().get_ID();
        if (maaid !=0) { maa = simec.hepevt().mother().idhep(); }
      }
    }

    if(maaid ==0) return 9876789;

    if (gen >1) { 
      for ( int p=0; p< gen-1; p++) {
        if ( abs(maaid) !=0 && abs(maa) !=300553 && abs(maa) != 10022) {
          for(std::vector<Gen_hepevt>::iterator it = hepevt_mag.begin();
              it != hepevt_mag.end(); ++it){
            Gen_hepevt& particle = *it;
            if (particle.get_ID()== maaid) {
              maa = particle.mother().idhep();
              maaid = particle.mother().get_ID();
              //      for(int ss=0; ss<gen;  ss++){  // for fun
              //                std::cout << "::::::::"; }
              //  std::cout << endl;
            }
          }
        } else{ return 9876789;}      }
    }

    return maa;
  }

	

  //____________________  FAMILY CHAIN ________________

  int fachain::daus( int pa, int da1, int da2) {
   
    Gen_hepevt_Manager& hepevt_mag = Gen_hepevt_Manager::get_manager();
    
    testok= 0;
    for(std::vector<Gen_hepevt>::iterator it = hepevt_mag.begin();
	it != hepevt_mag.end(); ++it){
      Gen_hepevt& particle = *it; 
     
      if( abs(particle.idhep()) == pa ){
	if(particle.da(1)-particle.da(0)!=1) { testok = 0; }
	else {
	  Gen_hepevt &child1 = *(hepevt_mag.begin() - 1 + 
				 particle.da(0));
	  Gen_hepevt &child2 = *(hepevt_mag.begin() - 1 +
				 particle.da(1));
	  if((abs(child1.idhep())==da1 && abs(child2.idhep())==da2) ||
	     (abs(child1.idhep())==da2 && abs(child2.idhep())==da1)){
	    testok=1; }
	  else  { testok = 0 ;}
	}
      } 
      
    }

    return testok;
  }   


  //_________________________________________________
  //     HELICITY ANGLE CALCULATION
  //_________________________________________________


  double hell::hel(Particle &pa, Particle &da, Particle &un) {
    
    HepLorentzVector Particle(da.momentum().p());
    double gamE = Particle.e();
    HepLorentzVector parent(pa.momentum().p());
    HepLorentzVector friendparent(un.momentum().p());
    
    HepLorentzVector grandparent;
    grandparent = friendparent + parent;
    Hep3Vector boosttoparent = -(parent.boostVector());
    Particle.boost(boosttoparent);
    grandparent.boost(boosttoparent);
    Hep3Vector particle3 = Particle.vect();
    Hep3Vector grandparent3 = grandparent.vect();
    double num = particle3.dot(grandparent3);
    double den = (particle3.mag())*(grandparent3.mag());
    double  hell = num/den;

    return hell;
  }



  double Kall::enrGY(Particle & partu, double angle){
    HepLorentzVector me(partu.momentum().p());
    sum_energy =me.e();
    /*Mdst_charged_Manager & ChgMgr= Mdst_charged_Manager::get_manager();
    Mdst_gamma_Manager & GamMgr = Mdst_gamma_Manager::get_manager();
    std:: vector<Particle> any;
    for(Mdst_charged_Manager::iterator p = ChgMgr.begin();
	p!=ChgMgr.end(); p++){
      Mdst_charged& Chg = *p;
      
      Particle ANY(*p, ptype_kaon)((*p).charge()<0.0) ? ptype_kaon_minus : ptype_kaon_plus, ipProf);
	    kaon.us
      // Mdst_trk& Trk = Chg.trk
      //Mdst_trk_fit& TrkFit0 = Trk.mhyp(0) ;//for electron
      // Mdst_trk_fit& TrkFit1 = Trk.mhyp(1) ;//for muon
      // Mdst_trk_fit& TrkFit2 = Trk.mhyp(2) ;//for pion
      //Mdst_trk_fit& TrkFit3 = Trk.mhyp(3) ;//for kaon
      //Mdst_trk_fit& TrkFit4 = Trk.mhyp(4) ;//for proton
      
      if (checkSame(*partu, *p)) continue;

      HepLorentzVector any((*p).momentum().p());


    }*/

    return sum_energy;
  }




dChain dChain::operator | (dChain rhs)
  {
    std::vector<int> temp;
    temp.push_back(da1);
    temp.push_back(rhs.da1);
    return temp;
  }
  
  dChain dChain::operator > (dChain da) {
    std::vector<int> temp;
    int sum=0;
    temp.reserve(1+ da.size());
    for(int i =0; i< da.size(); i++) {
      temp.push_back(da.value(i));
      if(da.value(i) <0 && da.value(i)!=-99) sum=sum-da.value(i)+1;
      if(da.value(i)==-99) sum = sum +1;
    }
    temp.push_back(-da.size() +sum );
    temp.push_back(da1);
    return temp; }

  dChain dChain::operator * (dChain da1) {
  std::vector<int> temp;
  temp.reserve(da.size()+ da1.size());
  temp.insert(temp.end(), da.begin(), da.end());
  for(int i =0; i< da1.size(); i++) {
    temp.push_back(da1.value(i));
  }
  //  temp.push_back(-99);
  return temp; }

  dChain dChain::operator + (dChain da){
    std::vector<int> temp;
    temp.reserve(1+ da.size());
    for(int i =0; i< da.size(); i++) {
      temp.push_back(da.value(i));
    }
    temp.push_back(da1);
    return temp; }
  


  int BeeList::c2_body(int cda1, int cda2, int da1, int da2 ){
    if ((cda1==da1 && cda2==da2) || (cda1==da2 && cda2==da1)) 
      { return 1;}
    else return 0;
  }
  
  
  int BeeList::c3_body(int cda1, int cda2, int cda3,
		       int da1, int da2, int da3 ){
    if (cda1==da1 && cda2==da2 && cda3==da3 ||
	cda1==da1 && cda3==da2 && cda2==da3 ||
	cda2==da1 && cda1==da2 && cda3==da3 ||
	cda2==da1 && cda3==da2 && cda1==da3 ||
	cda3==da1 && cda1==da2 && cda2==da3 || 
	cda3==da1 && cda2==da2 && cda1==da3) 
      { return 1;}
    else return 0;
  }
  
  
  int BeeList::c4_body(int cda1, int cda2, int cda3, int cda4,
		       int da1, int da2, int da3, int da4 ){
    
    if (cda1==da1 && cda2==da2 && cda3==da3 && cda4==da4  ||
	cda1==da1 && cda2==da2 && cda4==da3 && cda3==da4  ||
	cda1==da1 && cda3==da2 && cda2==da3 && cda4==da4  ||
	cda1==da1 && cda3==da2 && cda4==da3 && cda2==da4  ||
	cda1==da1 && cda4==da2 && cda3==da3 && cda2==da4  ||
	cda1==da1 && cda4==da2 && cda2==da3 && cda3==da4  ||
	
	cda2==da1 && cda1==da2 && cda3==da3 && cda4==da4  ||
	cda2==da1 && cda1==da2 && cda4==da3 && cda3==da4  ||
	cda2==da1 && cda3==da2 && cda1==da3 && cda4==da4  ||
	cda2==da1 && cda3==da2 && cda4==da3 && cda1==da4  ||
	cda2==da1 && cda4==da2 && cda3==da3 && cda1==da4  ||
	cda2==da1 && cda4==da2 && cda1==da3 && cda3==da4  ||
	
	cda3==da1 && cda2==da2 && cda1==da3 && cda4==da4  ||
	cda3==da1 && cda2==da2 && cda4==da3 && cda1==da4  ||
	cda3==da1 && cda1==da2 && cda2==da3 && cda4==da4  ||
	cda3==da1 && cda1==da2 && cda4==da3 && cda2==da4  ||
	cda3==da1 && cda4==da2 && cda1==da3 && cda2==da4  ||
	cda3==da1 && cda4==da2 && cda2==da3 && cda1==da4  ||
	
	cda4==da1 && cda2==da2 && cda3==da3 && cda1==da4  ||
	cda4==da1 && cda2==da2 && cda1==da3 && cda3==da4  ||
	cda4==da1 && cda3==da2 && cda2==da3 && cda1==da4  ||
	cda4==da1 && cda3==da2 && cda1==da3 && cda2==da4  ||
	cda4==da1 && cda1==da2 && cda3==da3 && cda2==da4  ||
	cda4==da1 && cda1==da2 && cda2==da3 && cda3==da4  )
      { return 1;}
    else return 0;
  }
  
  
  
  int BeeList::nofDaug(Gen_hepevt &mother){
    return mother.da(1)-mother.da(0) +1  ;}
  
  int BeeList::sameDaugMC(Gen_hepevt &mother, std::vector<int> daug){
    
    Gen_hepevt_Manager & hepevt_mag = Gen_hepevt_Manager::get_manager();
    int no_dau =nofDaug(mother);
    if (no_dau == 1) return 0;
    else if (no_dau ==2){
      Gen_hepevt & child1  =  *(hepevt_mag.begin() - 1 + 
				mother.da(0));
      Gen_hepevt & child2  =  *(hepevt_mag.begin() - 1 +
				mother.da(1));
      cda1=abs(child1.idhep());
      cda2=abs(child2.idhep());
      
      return c2_body(cda1,cda2,daug[0],daug[1]);
    }
    else if (no_dau==3){
      Gen_hepevt &child1= *(hepevt_mag.begin() -1 +
			    mother.da(0));
      Gen_hepevt &child2= *(hepevt_mag.begin() + mother.da(0));
      Gen_hepevt &child3= *(hepevt_mag.begin() -1 +
			    mother.da(1));
      cda1=abs(child1.idhep());
      cda2=abs(child2.idhep());
      cda3=abs(child3.idhep());
      
      return c3_body(cda1,cda2,cda3,daug[0],daug[1],daug[2]);
    }
    else if (no_dau==4){
      Gen_hepevt &child1= *(hepevt_mag.begin() -1 +
			    mother.da(0));
      Gen_hepevt &child2= *(hepevt_mag.begin() +mother.da(0));
      Gen_hepevt &child3= *(hepevt_mag.begin() -2 +
			    mother.da(1));
      Gen_hepevt &child4= *(hepevt_mag.begin() -1 +
			    mother.da(1));
      cda1=abs(child1.idhep());
      cda2=abs(child2.idhep());
      cda3=abs(child3.idhep());
      cda4=abs(child4.idhep());

      return c4_body(cda1, cda2,cda3,cda4,
			   daug[0],daug[1],daug[2],daug[3]);
    }
    else return 0;

  }


  std::vector<int> BeeList::fillD(int place, dChain pa){
   
    int listdau = -pa.value(place);
    int check = 0; int skip=0;
    std::vector<int>daug;
    daug.resize(listdau);
    for ( int i=0; i < listdau; i++){
      place = place -1;
      if(pa.value(place) > 0) {
	daug[i]=pa.value(place); }
      else {
	check = listdau - (i+1) - pa.value(place);
	
      minuscheck:
	for (int j =0; j< check; j ++){
	  place=place-1;
	  if( pa.value(place) < 0)
	    {check = check - (j+1) - pa.value(place);
	      goto minuscheck;
	    }
	}
	place = place -1;
	daug[i]=pa.value(place);
      }
    }
    return  daug; 
  }
  
  std::vector<int> BeeList::FillB(Particle &B, int dau){
    std::vector<int>daug;
    daug.resize(dau);
    int didU; 
    for (int i =0; i < dau ; i++){
      didU = abs(MAA.oiD(B.relation().child(i))) ;
      daug[i]=didU;
    }
    return daug;
  }



  int BeeList::compare(std::vector<int> A, std::vector<int> B){
    
    if(A.size() != B.size()) {return 0;}
    else{
      for (int i=0; i < A.size(); i++){
	for(int j=0;j<B.size();j++){
	  if(A[i]==B[j]) {B.erase(B.begin()+j);--j; }
	}
	}
      if (B.size() ==0 ) {return 1;}
      else return 0;
    }
  }


  int  BeeList::match(Particle & Bee, dChain pa){
    Gen_hepevt_Manager & hepevt_mag = Gen_hepevt_Manager::get_manager();
    int Beeid =  Bee.relation().genHepevt().get_ID();
    int Entry = pa.size();
    //std::cout <<  "BUGGGGG    "<<Entry << "\t"  <<Beeid << std::endl;
    //if (pa.value(Entry -1) != Beeid) {return 0;}
    
    int listdau = -pa.value(Entry-2);
   
    std::vector<int>daug;
    daug.resize(listdau);
    daug=fillD(Entry-2, pa);
   
    std::vector<int>GBdaug; 
    GBdaug.resize(listdau);
    GBdaug=FillB(Bee,listdau); 
    int check = compare(daug,GBdaug);
    //std::cout << daug.size()   << "\t" << GBdaug.size() << std::endl;
    //for(int i =0; i < daug.size() ; i++) {
    //  std::cout << "LOOPED  " << daug[i] << "\t" << GBdaug[i] << endl; }
    daug.clear(); GBdaug.clear();
    
    if(check!=1) return 0;
    //  std::cout << "PASSED THIS " << std::endl;
    int dauloop = 0;
    dauloop = Entry-4;
    
    int smarter = 0;
    int checker = 0;
    
    int gen_1,gen_2,gen_3,gen_4,gen_5,gen_6;

    if (Entry -4 > 0) {
      for (int pp = 0; pp< Entry-4; pp++){
	//std::cout << "LOOP BUG  " << pp << endl;
	int gen_1,gen_2,gen_3,gen_4,gen_5;
  
	dauloop = dauloop- 1;
	int partid=0;
	if(pa.value(dauloop) < 0){
	  daug.erase(daug.begin(),daug.end());
	  daug.resize(-pa.value(dauloop));
	  daug=fillD(dauloop,pa);
	  partid = pa.value(dauloop+1);
	  smarter =smarter+1; 
	  checker = -pa.value(dauloop) + checker; }
	else {checker = checker + 1;
	  if(checker ==  -pa.value(dauloop) +1 ) smarter = smarter -1 ; }
	
		
	for(int i=0; i<-pa.value(dauloop);i++){
	  GBdaug.erase(GBdaug.begin(),GBdaug.end());
	  GBdaug.resize(-pa.value(dauloop));
	  if(smarter ==1){  
	    if (Bee.relation().child(i).relation().genHepevt().idhep() == partid)
	      { gen_1 = i;  //std::cout << "\n GENERATION I \n";
	      GBdaug=FillB(Bee.relation().child(gen_1),-pa.value(dauloop));}}
	  else if (smarter==2){
	    if(Bee.relation().child(gen_1).relation().child(i).relation().genHepevt().idhep()==partid)
	      {gen_2=i; //std::cout << "\n GENERATION II \n";
		GBdaug=FillB(Bee.relation().child(gen_1).relation().child(gen_2), -pa.value(dauloop));}}
	  else if (smarter==3){//std::cout << "\n GENERATION III \n";
	    if(Bee.relation().child(gen_1).relation().child(gen_2).relation().child(i).relation().genHepevt().idhep()==partid)
	      { gen_3=i; 
		GBdaug=FillB((Bee.relation().child(gen_1).relation().child(gen_2).relation().child(gen_3)), -pa.value(dauloop));} }
	  else if(smarter==4){//std::cout << "\n GENERATION IV \n";
	    if(Bee.relation().child(gen_1).relation().child(gen_2).relation().child(gen_3).relation().child(i).relation().genHepevt().idhep()==partid) 
	      {gen_4 =i;
		GBdaug=FillB((Bee.relation().child(gen_1).relation().child(gen_2).relation().child(gen_3).relation().child(gen_4)), -pa.value(dauloop));}}
	   else if (smarter==5){
	    if(Bee.relation().child(gen_1).relation().child(gen_2).relation().child(gen_3).relation().child(gen_4).relation().child(i).relation().genHepevt().idhep()==partid) 
	      {gen_5 =i; //std::cout << "\n GENERATION V \n";
		GBdaug=FillB((Bee.relation().child(gen_1).relation().child(gen_2).relation().child(gen_3).relation().child(gen_4).relation().child(gen_5)), -pa.value(dauloop));}}
	    else if (smarter==6){
	    if(Bee.relation().child(gen_1).relation().child(gen_2).relation().child(gen_3).relation().child(gen_4).relation().child(gen_5).relation().child(i).relation().genHepevt().idhep()==partid) 
	      {gen_6 =i; //std::cout << "\n GENERATION VI \n";
	      GBdaug=FillB((Bee.relation().child(gen_1).relation().child(gen_2).relation().child(gen_3).relation().child(gen_4).relation().child(gen_5).relation().child(gen_6)), -pa.value(dauloop));}}
	  else {return 0; }
	  
      } // dauloop 
	
	if (smarter != 0)
	  {
	    
	    ///comp(daug,GBdaug);
	    
	    //	  check = compare(daug,GBdaug);    
	    
	    if(check !=1) return 0; }
	
      } //Loop over all 

      
      //as code survives here return 1;
    } //entry-4=0; 
    return 1;
  }
















#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif


  
