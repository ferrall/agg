#import "niqlow"

/** pseudonyms for 0/1. @name Answers **/
enum{no,yes,Answers}
/** pseudonyms for 0/1. @name Gender **/
enum{male,female,Ngender}

/** Marital status. @name Status **/
enum{single,married,divorced,Nstatus}

/** Sub-periods while adult. @name SubT **/
enum{mingle,nest,SubT}

/** ages while adult. @name Age **/
enum{young,old,Ages}  //names for lifecycle

/** Dimensions of the state space. @name Dims **/
enum{Tinit=1,T = Tinit+Ages*SubT,Nx=17,Nh=2}

/** A SimpleJump random variable that remains constant between
ages and if the match survives. **/
struct SkipJump : SimpleJump {
	const decl mstat;
	SkipJump(mstat);
	Transit(FeasA);
	}

/** Marital status determined by both genders (fixed effects) saying yes to the match.**/
struct MaritalStatus : Random {
	enum{cur,last,Neqiter}
	static const decl AdjRate = 0.8;
	const decl props, gen;
	decl bothyes, eqiter, matters,
		/** vector of proposal acceptance probabilities;
		set in `AGG::Predict` and used by `AGG::status`
		**/											Pacc;
	MaritalStatus(gen,props);
	InitializeP(Nstates);
	EvalEquilibrium(printit);
	Adjust(ii,newp);
	Transit(FeasA);
	}

/**A Tauchen variable that changes only between ages, not between sub-periods.**/	
struct AGGTauchen : Tauchen {
	AGGTauchen(L,N,M,mu, sig,rho);
	Transit(FeasA);
	}
	
struct AGG : ExPostSmoothing {
	static const decl
			TimeEndowment=1.0,
			/** &rho;<sub>x</sub>, autocorrelation.**/ rho=0.7,
			/** &sigma;<sub>x</sub>, variance .**/ 	sig =0.4,
			/** &beta;, discount factor. **/		delt=0.67,
//			/** .**/ 								zeta=0.325,
//			eta=0.5,
			c2 = 4.2, EE=4.2,
			alph = 0.4,
			/**AGG &lt; &theta;<sub>1</sub>,&delta;<sub>1</sub> &gt; . **/
				xtheta=<0.1,0.5>,	
			/**AGG &lt; &theta;<sub>2</sub>,&delta;<sub>2</sub> &gt; . **/
				xdelta=<0.7,0.9>,
			/** Welfare benefit as fraction of average earnings.**/	
			WELFRATE = 0.22,
			/** Divorced father child support rate.**/	
			CSUPP = 0.10, TaxRate=.03,
			GuessAvgIncome = 20.0;

	static 	decl
			/** the solution method object. **/ viter,
				  Medx,
				  eqiter,
				  MINcons,
				  eqdist,
				  muz,
				  newmu,
				  avginc,
				  lastolg,
				  lastlap;
	/** Variables that are set and used by Utility; their
		values should NOT be used outside utility functions. **/
	static decl myx, wedded, split, ow, q, lcont, c;

	static 	decl
	/** Gender Fixed Effect. **/ 					 g,
	/** Propose/Accept choice variable. **/ 		 ask,
	/** Accept welfare choice variable. **/			 welf, 
	/** array of gender-specific productivities. **/ z,
	/** AGG call this &gamma;, bliss value. **/		 h,
	/** marital status state. **/					 status;

	decl
	/** child learning. set by female **/				 	Q,
	/** vector of labor supplies for married states. **/ 	lsup,
	/** Mother's time teaching children. **/			 	t,
	/** vector of personal incomes; used for child
		support and equilibrium level of welfare rates.**/ 	inc;

	static Setup();
	static Reachable();
	static IsNest();
	static IsChild();
	static Discount();
	static mnzM();
	static mnzF();
	static MaritalEquilibrium();
	static EQ();
	
	AGG();
	Utility();
	Smooth(EV);
	FeasibleActions(A);
	Predict(ps,tod);
	MaleU();
	FemaleU();
	}
	