#include "myutl.h"
#include "belle.h"
#include "event/BelleEvent.h"
#include "tuple/BelleTupleManager.h"
#include "basf/module.h"
#include "basf/module_descr.h"
#include "panther/panther.h"
#include <iostream>
#include MDST_H
#include BELLETDF_H
#include HEPEVT_H
#include "particle/utility.h"
#include "particle/combination.h"
#include "kid/atc_pid.h"
#include "mdst/mdst.h"
#include "mdst/Muid_mdst.h"
#include "eid/eid.h"
#include "ip/IpProfile.h"
#include "particle/Particle.h"
#include "particle/Ptype.h"
#include "particle/ParticleUserInfo.h"

using namespace std;



#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif
typedef std::vector<Particle> p_list;
typedef p_list::iterator p_list_it;

class User_reco : public Module {
public:

  User_reco ( void );
  ~User_reco ( void ){};
  void init ( int* ){}; //called when basf is being initialized
  void term ( void );//called when basf is being terminated
  void disp_stat ( const char* ){};//auxillary function
  void hist_def ( void );//called for histogram definition (descibed in basf class)
  void event ( BelleEvent*, int* ); //called at start of each event
  void begin_run ( BelleEvent*, int* );//called at start of each run
  void end_run   ( BelleEvent*, int* ){};//called at end of each run
  void other ( int*, BelleEvent*, int* ){};//auxillary function


//functions above are mandatory for module
    void fill_tup_d0(Particle &p);
    std::string t_par(std::string tag);
    void t_fill(Particle &p,BelleTuple *t,std::string tg);
    void t_fill_ups(Particle &p,BelleTuple *t);


    //void makePi0();
    //void makeKs(p_list &part);




  /*******************************************
    Define things here that you want to use
    inside and outside of the event function
   *******************************************/
  BelleTuple *t0;
  BelleTuple *t1;
  BelleTuple *t2;
  BelleTuple *t3;
  BelleTuple *t4;
  BelleTuple *t5;
  BelleTuple *t6;
  BelleTuple *t7;
  BelleTuple *t8;
  BelleTuple *t9;
  BelleTuple *t10;


  BelleTupleManager *btm;


private:


   int expNo;
   int runNo;
   int evtNo;
   bool isMC;
   motherid Get;
};

//Function below is MANDATORY, it created User_reco class, and make BASF recognize shared library as module
extern "C" Module_descr *mdcl_User_reco ()
{ /* main */
  User_reco *module = new User_reco;
  Module_descr *dscr = new Module_descr ( "User_reco", module );
  return dscr;
};

/* Constructor  */
User_reco::User_reco ( void )
    : expNo(-1),
    runNo(-1),
    evtNo(-1),
    isMC(0)
{
  /*******************************************
   Define things here that you want to use
   inside and outside of the event function
   ******************************************/



};

//Run begin
void User_reco::begin_run ( BelleEvent*, int* )
{
  eid::init_data();
  IpProfile::begin_run();//Get info about IP
  IpProfile::dump();//Same

}

void User_reco::term(void){
}


void setUserInfo(Particle &p, unsigned ch=0);
void setUserInfo(std::vector<Particle> &p, unsigned ch=0);


//This class can store additional parameters for each particle (non-vertexed mass, channel, dgf, etc)
class UserInfo : public ParticleUserInfo
{
public:
  /// Default constructor
  UserInfo();

  UserInfo(unsigned);

  /// Copy constructor
  UserInfo(const UserInfo &);

  /// Destructor
  virtual ~UserInfo();

  /// constructs self object.
  UserInfo * clone(void) const;

  /// Copy operator
  UserInfo & operator = (const UserInfo &);

public:
  void vchisq(const double &v) { m_vchisq = v; }
  void prob_vchisq(const double &v) { probm_vchisq = v; }
  void dgf(const double &v) {m_dgf = v;}
  const double & vchisq(void) const { return m_vchisq; }
  const double & prob_vchisq(void) const { return probm_vchisq; }
  const double & dgf(void) const { return m_dgf;}
  void vmass(const double &v) { m_vmass = v; }
  const double & vmass(void) const { return m_vmass; }
  //const unsigned & index(void) const {return m_index;}
  void channel(const unsigned &v) { m_channel = v; }
  const unsigned & channel(void) const { return m_channel; }


private:
  double   m_vchisq;
  double   probm_vchisq;
  double   m_dgf;
  double   m_vmass;
  unsigned m_channel;
 // unsigned m_index;
};

int IDhep(Particle &part);//Called for getting MC particle id
int makeKvFit(Particle &p, int flg);//vertex fitter
void makeKvFit(p_list &plist, int flg);
int makeKmvFit(Particle &p);//vertex fitter
void makeKmvFit(p_list &plist);
int pidKvspi(Particle &p, float prob);
int pidpivsK(Particle &p, float prob);
void pidKvspi(p_list &plist, float prob);
void pidpivsK(p_list &plist, float prob);
int withDrDz(Particle &p, double maxdr, double maxdz);
void withDrDz(p_list &plist, double maxdr, double maxdz);
void minEn(p_list &plist, double minE);
int minEn(Particle &p, double minE);
void minP(p_list &plist, double minp);
int minP(Particle &p, double minp);
void minCh2(p_list &plist, double minch2);
int minCh2(Particle &p, double minch2);
void maxCh2(p_list &plist, double maxch2);
int maxCh2(Particle &p, double maxch2);
int minPst(Particle &p, double minPst);
void minPst(p_list &plist, double minPst);
int IDcheck(Particle &p);
void IDcheck (p_list &plist);

void combination5p(p_list &plist, const Ptype &ptype, p_list &p1, p_list &p2, p_list &p3, p_list &p4, p_list &p5);
void combination5p(p_list &plist, const Ptype &ptype, p_list &p1, p_list &p2, p_list &p3, p_list &p4, p_list &p5,const double &marg);
void combination5p(p_list &plist, const Ptype &ptype, p_list &p1, p_list &p2, p_list &p3, p_list &p4, p_list &p5, const double &larg, const double &rarg);
void combination5p(p_list &plist, const Ptype &ptype, Particle &p1, Particle &p2, Particle &p3, Particle &p4, Particle &p5);
void combination5p(p_list &plist, const Ptype &ptype, Particle &p1, Particle &p2, Particle &p3, Particle &p4, Particle &p5,const double &marg);
void combination5p(p_list &plist, const Ptype &ptype, Particle &p1, Particle &p2, Particle &p3, Particle &p4, Particle &p5, const double &larg,const double &rarg);

void combination6p(p_list &plist, const Ptype &ptype, p_list &p1, p_list &p2, p_list &p3, p_list &p4, p_list &p5, p_list &p6);
void combination6p(p_list &plist, const Ptype &ptype, p_list &p1, p_list &p2, p_list &p3, p_list &p4, p_list &p5, p_list &p6,const double &marg);
void combination6p(p_list &plist, const Ptype &ptype, p_list &p1, p_list &p2, p_list &p3, p_list &p4, p_list &p5, p_list &p6, const double &larg, const double &rarg);
void combination6p(p_list &plist, const Ptype &ptype, Particle &p1, Particle &p2, Particle &p3, Particle &p4, Particle &p5, p_list &p6);
void combination6p(p_list &plist, const Ptype &ptype, Particle &p1, Particle &p2, Particle &p3, Particle &p4, Particle &p5, p_list &p6,const double &marg);
void combination6p(p_list &plist, const Ptype &ptype, Particle &p1, Particle &p2, Particle &p3, Particle &p4, Particle &p5, p_list &p6, const double &larg,const double &rarg);
HepVector imp_par_ip(const Mdst_charged& charged, int hyp);
#if defined(BELLE_NAMESPACE)
}
#endif
