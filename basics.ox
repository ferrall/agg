
/** indicator for t&lt;<code>Tinit</code>. **/
AGG::IsChild() { return curt<Tinit;}
AGG::IsNest() {	return IsChild() || imod(curt-Tinit,SubT);	}
AGG::Age() { return IsChild() ? -1 : idiv(curt-Tinit,SubT); }
AGG::Discount() { return IsNest() ? delt : 1.0; } 
AGG::AGG() {
	mtchind = ind[tracking];
	mtchind += OO[tracking][othz.pos]*(CV(myz)-CV(othz));
	mtchind += OO[tracking][myz.pos]*(CV(othz)-CV(myz));
	t = Q = 0.0;  //initial value does not matter because female problem solved first
	if (IsNest()){
		if (CV(status)==married) inc = constant(0.5,Ngender,1);
		else if (CV(status)==divorced && CV(g)==male) 
			inc = 0.0;
		}
	else {  //mingling 
		if (CV(status)!=divorced) saysyes = ones(Ngender,1);
		}
	}

AGG::Setup() {
	Medx = int( (Nx-1)/2 );
	TaxRate =.03;
	Initialize(Reachable,0,0);
	MINcons = DBL_EPSILON;			// too close to 0 to take log()
	SetClock(NormalAging,T);
	SetDelta(Discount);
	g = new FixedEffect("g",Ngender);
	AGGTauchen::Initialize();
	AVGINCOME = exp(AGGTauchen::EE);
	GroupVariables(g);
	myz = new MyProductivity();	
	Actions(ask = new ActionVariable("IDo",Answers),
			welf = new ActionVariable("welf",Answers));
	status = new MaritalStatus(g,ask);
	h = new SkipJump(status);
	othz = new Match(status);
	EndogenousStates(h,status,myz,othz); //,lm
	CreateSpaces(LogitKernel,Smthing);
	viter  = new ValueIteration(0);
	PostRESolve=EQ;
	println("Parameter Values",
		"\n EE=",AGGTauchen::EE,
		"\n c2=",AGGTauchen::c2[male],
		"\n rho=",AGGTauchen::rho,
		"\n sig=",AGGTauchen::sig,
		"\n alph=",alph,
		"\n xtheta=",xtheta,
		"\n xdelta=",xdelta,
		"\n bliss = ",h.actual',
		"\n welfare rate=",WELFRATE,
		"\n Child Support rate=",CSUPP,
		"\n choice smoothing=",Smthing);	
	}

decl SubToPrint;	
ToScreen(p) {	return ( (p.t==0 || imod(p.t-Tinit,SubT))==SubToPrint);	}

SkipJump::SkipJump(mstat) {
	SimpleJump("h",Nh);
	this.mstat = mstat;
	actual = <0.0;2.6>; // See AGG Table 1 2.6
	}
	
SkipJump::Transit(FeasA) {
	if (AGG::IsChild()|| AGG::IsNest() || CV(mstat)!=single)  
		return  UnChanged(FeasA);
	decl nv, p;
	[nv,p] = SimpleJump::Transit(FeasA);
	return {nv,mstat.bothyes.*p+(1-mstat.bothyes).*(v.==nv)};
	}

	
AGG::FeasibleActions(A) {
	decl nesting=IsNest(),civ = CV(status);
	if (curt<Tinit  //can't choose anything if a child 
		|| (!nesting &&civ==divorced) //or divorced and mingling
		|| (nesting&&civ==married) //or married and nesting
		)  
		return  A[][ask.pos].==no .&& A[][welf.pos] .==no;
	if (!nesting||civ==married) return A[][welf.pos].==no; //can't choose welfare when mingling or married, 
	return A[][ask.pos].==no;  //can't accept when nesting
	}

AGG::Reachable() {
	if ( curt<Tinit && (CV(status)!=single || CV(myz)!=Medx || CV(othz)!=Medx || CV(h)!=0 ) )
		return 0;	/* At birth, everyone is at the mean/median productivity household bliss irrelevant */
	if (curt==Tinit && (CV(status)!=single || CV(h)!=0))
		return 0;  // start out single and bliss irrelevant
	if (curt<T-1 && (CV(status)==divorced))
		return 0;  // can't be divorced until last period
	return new AGG();
	}
	
AGG::Utility() {
	if (IsChild())	return <0.0>; 			//  no choice/utility if child
	if (!IsNest())	return zeros(Answers,1); //no utility associated with mingle choice
	decl myg = CV(g);
	qmtch = Settheta(mtchind);
	wedded = CV(status)==married;
	split = CV(status)==divorced;
	myx = (1-TaxRate)*(1-CSUPP*split*(myg==male))*exp(AV(myz));	 //net wage rate
	ow = aa(welf);
	c =  wedded ? qmtch.inc[!myg] - AV(h) : 0.0;
	q =  ((myg==female)||(wedded))*(1-alph)*xtheta[myg];
	decl U = (myg==male) ? (1-ow)*MaleU() + ow*DBL_MIN_EXP : FemaleU();
	grossinc = exp(AV(myz))*lcont;
	if (wedded) inc[myg] = icont;
	else if (split&&myg==male) inc = icont;
	return U;
	}

AGG::MaleU() {
	lcont = setbounds(
			( myx*(1+ q) - xdelta[male]*c )	/  (myx*(1+q+xdelta[male])),
			0.0,1-MINcons );
	icont = (myx*lcont);
	c += icont;
	c = setbounds(c,MINcons,+.Inf);
	welfrcpt = 0.0;
	return 		log(c) 										// U
			+	xtheta[male]*(wedded ? qmtch.Q : 0.0)  		// V
			+	xdelta[male]*log(TimeEndowment-lcont);  	//R
	}

AGG::FemaleU() {
	decl A = xtheta[female]*alph / (xdelta[female]+xtheta[female]*alph), 
		B = (1+q)*myx;
	lcont = (1-ow)*		//can't work if on welfare
			setbounds(
				((1-A)*B - xdelta[female]*c) / (myx*((1-A)*(1+q)+xdelta[female])),
				0.0, 1-MINcons);
	icont = myx*lcont;
	welfrcpt = ow*WELFRATE*AVGINCOME;
	c += welfrcpt  + icont ;   
	if (split) c += CSUPP * qmtch.inc;
	c = setbounds(c,MINcons,+.Inf);
	t = A*(TimeEndowment-lcont); 
	Q = (1-alph)*log(c) + alph*log(t);
	return 	log(c)	// U
			+xtheta[female]*Q			 // V
			+xdelta[female]*log(TimeEndowment-lcont-t);  //R
	}

MaritalStatus::MaritalStatus(gen,props) {
	StateVariable("mstat",Nstatus);
	this.gen = gen;
	this.props = props;
	}

MaritalStatus::Transit(FeasA) {
	if (AGG::IsChild()||AGG::IsNest()) 
		return UnChanged(FeasA);  //marital status decided in adult mingle periods
	bothyes = DP::Settheta(DP::ind[tracking]).saysyes[!CV(gen)] * (FeasA[][props.pos].==yes);
	return { (v==married) ? <divorced,married>
						  : <single,married> , (1-bothyes)~bothyes };
	}

AGGTauchen::Initialize() {	g = AGG::g;	PrevQ = mnQ = 1.0;	}
AGGTauchen::MyMean() {	return c2[CV(g)] + EE*mnQ;}
AGGTauchen::OthMean() {	return c2[!CV(g)]+ EE*mnQ;}
	
MyProductivity::MyProductivity() { Tauchen("myz",Nx,2.0,MyMean,sig,rho);	}
Match::Match(status) {
	Tauchen("othz",Nx,2.0,OthMean,sig,rho);
	this.status = status;
	Previous = StillSingle = constant(1/Nx,2,Nx); 
	}

AGGTauchen::Update() {	Tauchen::Update();	}
	
MyProductivity::Update() {
	AGGTauchen::Update();
	//	mnQ = 0.5*mnQ + 0.5*PrevQ;
	mnQ = PrevQ;
	}

Match::Update() 		 {
	AGGTauchen::Update();
	//	StillSingle = 0.5*StillSingle + 0.5*Previous;
	StillSingle = Previous;
	}	
MyProductivity::Transit(FeasA) {
	if (AGG::IsNest()) return Tauchen::Transit(<0>);
	return UnChanged(FeasA); //no change while mingling
	}

Match::Transit(FeasA) {
	if (!AGG::IsNest())  return UnChanged(FeasA); //no change while mingling
	if (AGG::IsChild()||(CV(status)!=single)) return Tauchen::Transit(<0>);
	//single person gets a new draw from equilibrium 
	return	{vals,StillSingle[!CV(g)][]};	
	}
	