#include "jw.h"

main() {
	AGG::Setup();
	AGG::MaritalEquilibrium();
	}
	
/** indicator for t&lt;<code>Tinit</code>. **/
AGG::IsChild() { return curt<Tinit;}
/** indicator that t is either childhood or second sub-period when
	labour supply is decided; otherwise the period is mingling
	and the marriage choices are made. **/
AGG::IsNest() {	return IsChild() || imod(curt-Tinit,SubT);	}

/** returns the discount factor &delta;, called &beta; by AGG;
No discounting between sub-periods. **/
AGG::Discount() { return IsNest() ? delt : 1.0; } 

/** Called when creating a new point in the state space, &theta;;
This initializes the state-specific (non-static) variables used in
the model.  **/
AGG::AGG() {
	t = Q = 0.0;  //initial value does not matter because female problem solved first
	if (IsNest()){
		if (CV(status)==married || CV(status)==divorced ) {
			lsup = constant(0.5,Ngender,1);  //wives initially believe husband works 50%
			inc = zeros(1,Ngender);			
			}
		else	  // when not married women can choose welfare, so income depends on choice and gender.
			inc = zeros(Answers,Ngender); //divorced women initially believe father is a deadbeat
		}
	else {  //mingling 
		inc = 0.0;   
		}
	}

/** A virtual replacement for built-in Bellman::Predict().
@param ps probaiblity of arriving at the current state (computed internally
		by moving forward from initial conditions.)
@param tod Outcomes today
@comment call Bellman::Predict first, then compute state's contribution
to aggregate outcomes.
**/
AGG::Predict(ps,tod){
	Bellman::Predict(ps,tod);
	if (IsNest()) { //only do this for labour supply states
		newmu += ps*sumc(pandv[0][].*Q);  			//average human capital transmission
		avginc+= ps*sumc(pandv[0][].*meanr(inc)); //every state covers 2 people, so average across gender
		}
	}

decl SubToPrint;	
ToScreen(p) {
	return ( (p.t==0 || imod(p.t-Tinit,SubT))==SubToPrint);
	}

/** Static function that is assigned to DP::PostRESolve; This will
be called after each fixed effect (gender) specific problem is fixed.
Only want to something when marital equilibrium is found, and only needs
to be done once (after female) because each state handles two people
**/
AGG::EQ() {
	if (lastlap && CV(g)==female) {
		avginc = newmu = 0.0;
		decl p = new PanelPrediction(0);
		p->Predict(T);
		println("Average personal earnings = ",double(avginc));
		println("Average human capital transfer = ",double(newmu));
		if (lastolg) {
			SubToPrint=nest;
			p->Histogram(status,ToScreen,TRUE);
			p->Histogram(welf,ToScreen,TRUE);
			SubToPrint=mingle;
			p->Histogram(ask,ToScreen,TRUE);
			}
		delete p;
		}
	}

/** Check for Nash behaviour in marriages, the update human capital and average income.**/
AGG::MaritalEquilibrium() {
	decl toler = 1E-4, notdone, olgdone,k ;
	lastolg=FALSE;
	do {  //human capital loop
		lastlap = FALSE;
		do { //married behaviour loop
			eqdist = 0.0;
			viter->Solve(AllFixed,0);
			eqdist += status->EvalEquilibrium(TRUE);
			notdone = !lastlap;  		//done after lastlap is run
			lastlap = (eqdist<toler);	//lastlap when tolerance met
			} while (notdone);
		olgdone = lastolg;
		lastolg = sqr(muz-newmu)<toler;  //change in mean productivity
		muz = newmu;	//update productivity.
		} while (!olgdone);
	}
	
/** Virtual replacement for Bellman::Smooth.
This will be called after value of state has been computed.
The built in smooth will compute optimal choice probabilities.
So this is the time to record chance that match will say yes and to
compute difference between their beliefs and my optimal choice probability.
**/
AGG::Smooth(EV) {
	ExPostSmoothing::Smooth(EV);
	if (IsChild()) return;	
	if (!IsNest() && CV(status)!=divorced) status->Adjust(ind[tracking],pandv[0][yes]);
	}

/** See AGG page 218. **/
AGG::mnzF() { return EE*muz; }
/** See AGG page 218. **/
AGG::mnzM() { return c2+ mnzF(); 	}
	
AGG::Setup() {
	Medx = int( (Nx-1)/2 );
	Initialize(Reachable,0,0);
	MINcons = DBL_EPSILON;			// too close to 0 to take log()
	muz = -1.1;
	avginc = exp(mnzM())*0.8; // initial avginc for welfare
	println("Initial avginc = ",avginc);
	SetClock(NormalAging,T);
	SetDelta(Discount);
	GroupVariables(g = new FixedEffect("g",Ngender));
	Actions(ask = new ActionVariable("IDo",Answers),
			welf = new ActionVariable("welf",Answers));
	status = new MaritalStatus(g,ask);
	h = new SkipJump(status);
	z = new array[Ngender];
		z[male] = new AGGTauchen("zm",Nx,2,mnzM,sig,rho);
		z[female] = new AGGTauchen("zf",Nx,2,mnzF,sig,rho);
	EndogenousStates(status,h,z[male],z[female]);
	CreateSpaces(NoSmoothing,2.0);
	status->InitializeP(SS[tracking].size);
	viter  = new ValueIteration(0);
	PostRESolve=EQ;
	}

/** Virtual replacement function.
@param A matrix of all possible discrete choices. Each row is a combination
of choices; each column is a choice variable. An action variable's column is
stored in its <code>pos</code>.
@return a column of 0s/1s, where 1 indicates the choice combination is
feasible at the current state.

State variables have been set so that their current and actual (CV(),AV()) are
stored.  This is called once for each reachable state by CreateSpaces().

Feasible action sets CANNOT depend on the value of fixed effects (gender).
So formally single men have the option of going on welfare.  The male utility
has to ensure that they do not choose this option.

**/	
AGG::FeasibleActions(A) {
	decl nesting=IsNest();
	if (curt<Tinit  //can't choose anything if a child 
		|| (!nesting&&CV(status)==divorced) //or divorced and mingling
		|| (nesting&&CV(status)==married) //or married and mingling
		)  
		return  A[][ask.pos].==no .&& A[][welf.pos] .==no;
	if (!nesting||CV(status)==married) return A[][welf.pos].==no; //can't choose welfare when mingling or married, 
	return A[][ask.pos].==no;  //can't accept when nesting
	}

/** Return a new &theta; if the current value of state variables
cannot be reached from initial conditions. **/
AGG::Reachable() {
	if ( curt<Tinit && (CV(status)!=single || CV(z[male])!=Medx || CV(z[female])!=Medx || CV(h)!=0 ) )
		return 0;	/* At birth, everyone is at the mean/median productivity, so
						other productivities are unreachable
						household bliss irrelevant */
	if (curt==Tinit && (CV(status)!=single || CV(h)!=0))
		return 0;  // start out single and bliss irrelevant
	return new AGG();
	}
	
AGG::Utility() {
	if (IsChild())	return <0.0>; 			//  no choice/utility if child
	if (!IsNest())	return zeros(Answers,1); //no utility associated with mingle choice
	decl myg = CV(g);
	myx = exp(AV(z[myg]));			//these variables used by both genders
	wedded = CV(status)==married;
	split = CV(status)==divorced;
	ow = aa(welf);	
	c =  wedded ? lsup[!myg]*exp(AV(z[!myg]))-AV(h) : 0.0;
	q =  ((myg==female)||(wedded))*(1-alph)*xtheta[myg];
	decl U = (myg==male) ? MaleU() : FemaleU();
	if (wedded|| (split&&myg==male) ) {
//		eqdist += sqr(lsup[myg]-lcont);
		lsup[myg] = lcont;				//store my labour supply
//		println(myg,curt,wedded,split,CV(h),CV(z[0]),CV(z[1])," ",lcont,lsup,U);
		}
	return U;
	}

AGG::MaleU() {
	lcont = setbounds(
			( myx*(1+ q) - xdelta[male]*c )	/  (myx*(1+q+xdelta[male])),
			0.0,1-MINcons );
	inc[][male] = (1-ow)*(1-TaxRate)*myx*lcont;		 //CF added TaxRate
	
	c += inc[][male]*(1-CSUPP*split) ;
	c = setbounds(c,MINcons,+.Inf);
	decl Qcont = wedded ? ( (1-alph)*log(c) + alph*log(t) ) : 0.0;
	return 		log(c) 										// U
			+	xtheta[male]*Qcont			 			// V
			+	xdelta[male]*log(TimeEndowment-lcont);  //R
	}

AGG::FemaleU() {
	decl tcont, Qcont;	
	decl A = xtheta[female]*alph / (xdelta[female]+xtheta[female]*alph), 
		B = (1+q)*myx;
	lcont = (1-ow)*
			setbounds(
				((1-A)*B - xdelta[female]*c) / (myx*((1-A)*(1+q)+xdelta[female])),
				0.0, 1-MINcons);
	inc[][female] = (1-TaxRate)*myx*lcont;
	c += ow*WELFRATE*avginc + inc[][female] ;
	if (split) c += CSUPP*lsup[male]*exp(AV(z[male]));
	c = setbounds(c,MINcons,+.Inf);
	t = tcont = A*(TimeEndowment-lcont); 
	Q = Qcont = (1-alph)*log(c) + alph*log(t);
	return 	log(c)	// U
			+xtheta[female]*Qcont			 // V
			+xdelta[female]*log(TimeEndowment-lcont-tcont);  //R
	}

SkipJump::SkipJump(mstat) {
	SimpleJump("h",Nh);
	this.mstat = mstat;
	actual = <0.0;2.6>; // See AGG Table 1
	}
	
SkipJump::Transit(FeasA) {
	return (!AGG::IsChild() && AGG::IsNest()&&CV(mstat)!=married)
			? SimpleJump::Transit(FeasA)
	 		:  UnChanged(FeasA);
	}

MaritalStatus::MaritalStatus(gen,props) {
	StateVariable("mstat",Nstatus);
	this.gen = gen;
	this.props = props;
	}

MaritalStatus::InitializeP(Nstates) {
	decl i,g;
	Pacc = new array[Neqiter];
	for(i=0;i<Neqiter;++i) Pacc[i] = ones(Nstates,gen.N);
	matters = ones(Nstates,1);
	eqiter = last;
	}

MaritalStatus::EvalEquilibrium(printit) {
	decl dist,i;
	for(i=0,dist=0.0;i<gen.N;++i) dist+=norm(Pacc[last][][i]-Pacc[cur][][i],2);
	if (printit) println("dist = ",dist,MyMoments(selectifr(Pacc[eqiter],matters)));
	eqiter = !eqiter;
	return dist;
	}

MaritalStatus::Transit(FeasA) {
	decl ii = DP::ind[tracking];
	if (v==divorced) matters[ii] = FALSE;
	if (AGG::IsChild()||AGG::IsNest()) {
		matters[ii] = FALSE;
		return UnChanged(FeasA);  //marital status decided in adult mingle period
		}
	bothyes = Pacc[!eqiter][ii][!CV(gen)] * (FeasA[][props.pos].==yes);
	return { (v==married) ? <divorced,married>
						  : <single,married> , (1-bothyes)~bothyes };
	}

MaritalStatus::Adjust(ii,newp) {
	Pacc[eqiter][ii][CV(gen)] = (1-AdjRate)*Pacc[eqiter][ii][CV(gen)] + AdjRate*newp;
	}
	
AGGTauchen::AGGTauchen(L,N,M,mu, sig,rho) {
	Tauchen(L,N,M,mu,sig,rho);
	}

AGGTauchen::Transit(FeasA) {  //nesting includes childhood
	decl myg = CV(AGG::g);
	return	(AGG::IsNest()) ? Tauchen::Transit(FeasA) : UnChanged(FeasA);
	}
		