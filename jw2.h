#import "niqlow"

/** pseudonyms for 0/1. @name Answers **/ 	enum{no,yes,Answers}
/** pseudonyms for 0/1. @name Gender **/ 	enum{male,female,Ngender}
/** Marital status. @name Status **/ 		enum{single,married,divorced,Nstatus}
/** Sub-periods while adult. @name SubT **/	enum{mingle,nest,SubT}
/** ages while adult. @name Age **/			enum{young,old,Ages}  //names for lifecycle
/** Dimensions . @name Dims **/				enum{Tinit=1,T = Tinit+Ages*SubT,Nx=7,Nh=2}
/** Stages of model solution .**/  			enum{EQEndow,EQMatch,EQCouple,Nstages}

struct AGGTauchen : Tauchen {
	static const decl
		/** &rho;<sub>x</sub>, autocorrelation.**/ rho=0.7,
		/** &sigma;<sub>x</sub>, variance .**/ 	sig =0.4,
			c2 = <0.5;0.0>,
			EE=1.0;  //4.2  4.2
	static decl g, PrevQ, mnQ, means;
	static MyMean();
	static OthMean();
	static Initialize();
	virtual Update();
	}
	
/** A SimpleJump random variable that remains constant between ages and if the match survives. **/
struct SkipJump : SimpleJump {
	const decl mstat;
	SkipJump(mstat);
	Transit(FeasA);
	}
	
/** Marital status determined by both genders (fixed effects) saying yes to the match.**/
struct MaritalStatus : Random {
	const decl props, gen;
	decl bothyes;
	MaritalStatus(gen,props);
	EvalEquilibrium(printit);
	Transit(FeasA);
	}

struct MyProductivity : AGGTauchen {
	MyProductivity();
	Transit(FeasA);
	Update();
	}

struct Match : AGGTauchen {
	static decl Previous, StillSingle;
	const decl status;
	Match(status);
	Transit(FeasA);
	Update();
	}
	
struct AGG : ExPostSmoothing {
	static const decl
			TimeEndowment=1.0,
			delt = .67,
			alph = 0.4, //0.4,
			/**AGG &lt; &theta;<sub>1</sub>,&delta;<sub>1</sub> &gt; . **/
				xtheta=<0.1,0.5>,	
			/**AGG &lt; &theta;<sub>2</sub>,&delta;<sub>2</sub> &gt; . **/
				xdelta=<0.7,0.9>,
				Smthing = 3.0,
			/** Welfare benefit as fraction of average earnings.**/	
			WELFRATE = 0.22,
			/** Divorced father child support rate.**/	
			CSUPP = 0.10;
	static 	decl
				  viter,
				  TaxRate,
				  AVGINCOME,
				  eqdist,
				  eqqty,
				  last,
				  Medx,
				  MINcons;
	/** Variables that are set and used by Utility; their
		values should NOT be used outside utility functions. **/
	static decl lcont, myx, qmtch, wedded, split, ow, q, icont, c;

	static 	decl
	/** Gender Fixed Effect. **/ 					 g,
	/** Own productivity. **/						 myz,
	/** Match productivity. **/						 othz,
	/** Propose/Accept choice variable. **/ 		 ask,
	/** Accept welfare choice variable. **/			 welf, 
	/** AGG call this &gamma;, bliss value. **/		 h,
													 status,
	/** lagged marital status. **/					 lm;

	decl
	/** tracking index of matched person **/				mtchind,
	/** prob. matched person says yes **/					saysyes,
	/** gross pre-tax personal income. **/					grossinc,
															welfrcpt,
	/** child learning. set by female **/				 	Q,
	/** Mother's time teaching children. **/			 	t,
	/** personal incomes; used for equilibrium
		level of welfare rates.**/ 							inc;

	static Setup();
	static Reachable();
	static IsNest();
	static IsChild();
	static Age();
	static Discount();
	static MaritalEquilibrium();
	static EQ();
	
	AGG();
	Utility();
	Smooth(EV);
	FeasibleActions(A);
	ConditionChoice(p);
	Predict(ps,tod);
	MaleU();
	FemaleU();
	OutputValue();
	}
	