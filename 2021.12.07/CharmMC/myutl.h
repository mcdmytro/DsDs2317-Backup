#include "belle.h"
#include <string>
#include <vector>

#include "belleCLHEP/Vector/LorentzVector.h"

#include "tuple/BelleTupleManager.h"

#include "particle/Particle.h"

#include "kfitter/kvertexfitter.h"
#include "kfitter/kmassfitter.h"
#include "kfitter/kmassvertexfitter.h"
#include "kfitter/kmakemother.h"
#if defined(BELLE_NAMESPACE)
namespace Belle{
#endif


  class myVfit 
  {
  public:
    myVfit(){};
    ~myVfit(void){};
    unsigned VXfit(Particle & p, double & confLevel, int debug);
    HepPoint3D VX();
    HepSymMatrix erVX();
    double chisq();
    int dgf();
    int f_err();
  private:
    kvertexfitter kv;
    int ferr;
    unsigned makevMother( kvertexfitter     &kv, Particle &mother );
    kmakemother kmm;
  };
  
  
  class myMfit
  {
  public:
    myMfit(){};
    ~myMfit(void){};
    unsigned VMfit( Particle & p, double & confLevel, int flag,double mass );
    double chisq();
    double m_err();
    double mcl();
    
    
  private:
     kmassfitter km;
     int merr;
     unsigned makemMother( kmassfitter &kv, Particle &mother );
  }; 

  class party
  { public:
    party(){};
    ~party(void){};
    double getpT(HepLorentzVector &a);
    double getpT(Particle & p);
    HepLorentzVector mon4(Gen_hepevt & a);
    HepVector vex(Gen_hepevt & a);
    double genpT(Gen_hepevt &a);
    double genpTt(Gen_hepevt &a);
 
  };


  class motherid
  {
  public:
    motherid(){};
    ~motherid(void){};
    double thetagama(Particle & pa);
    int MiD( Particle & pa , int genration);
    int oiD( Particle & pa);
    int pi0Tinfo (Particle &da, int nkid, int noTri);
    int GTinfo(int eclid, int noTri);
    int Gisthep(int eclid);
    int pi0trumaa(Particle &PI0);
    int pi0Gmaa (Particle & da, int nkid, int gen);

    int GMiD(int eclid, int gen);
    int KsPimaa(Particle &K_S, int nkid, int gen);
    double helPi0(Particle &PI0);
   private:
    int my_id;
    int testid;
    int ma_id;
  };
    

  class fachain
  {
  public:
    fachain(){};
    ~fachain(void){};
    int daus (int ma, int da1, int da2);
    /*int daus (particle & pa; int da1, int da2, int da3);
    int daus (particle & pa; int da1, int da2, int da3, int da4);
    int daus (particle & pa; int da1, int da2, int da3, int da4, int da5);
    */
    private:
    int okay; int testok;
  };

  class hell
  { 
  public:
    hell(){};
    ~hell(void){};
    double hel(Particle & parent, Particle & child, Particle & uncle);
  private:
    HepLorentzVector temp;
  };


  class Kall
  {
  public:
    Kall(){};
    ~Kall(void){};
    double enrGY(Particle & partu, double angle);
  private:
    double sum_energy;
  };






 class dChain
  {
    public:
      dChain(){}
    dChain(std::vector<int> vi){ 
      da.reserve(vi.size());
      da.insert(da.end(), vi.begin(), vi.end());
    }
    
  
    dChain(int x){da1 = x;}
    
    int size(){
      return da.size(); }
    
    int value(int i){
      return da[i];}

    void show(){//std::cout << da.size() << std::endl;
      for(int i=0; i< int(da.size()); i++){
	std::cout << da[i] << std::endl; 
      }}
  


    dChain operator | (dChain rhs);
    dChain operator + (dChain da);
    dChain operator > (dChain da);
    dChain operator * (dChain da1);
  private:
  
    int da1,da2;
    int size1, size2;
    std::vector<int> da;
    std::vector<int> das;

  };
    

  class BeeList{
  public:
    friend class dChain;
    motherid MAA;
    int match(Particle & Bee, dChain pa);
    int nofDaug(Gen_hepevt &mother);
    int sameDaugMC(Gen_hepevt& mother, std::vector<int> daug);
    int sameDaug(Particle & Bee, std::vector<int> daug);
    int matchdau(int place, std::vector<int>daug);
    int getDaug(int x, int y);
    int c2_body(int cda1, int cda2, int da1, int da2);
    int c3_body(int cda1, int cda2, int cda3, int da1, int da2, int da3);
    int c4_body(int cda1, int cda2, int cda3, int cda4, int da1, int da2, int da3, int da4);
    std::vector <int> fillD(int place, dChain pa);
    std::vector <int> FillB(Particle & B, int dau);
    int compare(std::vector<int> A, std::vector<int>B);
  private:
    int cda1, cda2, cda3, cda4, da1, da2, da3, da4;
    int flag;
    std::vector<int> list;
  };





#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
