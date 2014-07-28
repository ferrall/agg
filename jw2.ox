#include "jw2.h"
#include "basics.ox"

main() {
	AGG::Setup();
	AGG::MaritalEquilibrium();
	}

AGG::MaritalEquilibrium() {
	decl toler = <1E-1,1E-5,1E-4>, notdone;
	eqdist = notdone = last = zeros(Nstages,1);
	eqqty = new array[Nstages];
	eqqty[EQEndow] = new array[2];
	do {  
		last[EQMatch] = FALSE;
		do {
			last[EQCouple] = FALSE;
			do { //married behaviour loop
				eqdist[]=0.0;
				viter->Solve(AllFixed,0);
				notdone[EQCouple] = !last[EQCouple];
				last[EQCouple] = sqrt(eqdist[EQCouple][0])<toler[EQCouple];
				} while (notdone[EQCouple]);
			notdone[EQMatch] = !last[EQMatch];
			last[EQMatch] = sqrt(eqdist[EQMatch])<toler[EQMatch];  
			} while (notdone[EQMatch]);
		notdone[EQEndow] = !last[EQEndow];
		println("@ ",double(eqdist[EQEndow]));
		last[EQEndow] = sqrt(eqdist[EQEndow])<toler[EQEndow];  
		} while (notdone[EQEndow]);
	println("Equilibrium Singles Matching ","%c",{"v","Male Values","Females expect","Female Values","Males expect"},
					"%cf",{"%3.0f","%13.2f","%14.4f","%13.2f","%14.4f"},
					myz.vals'~exp(myz.actual')~Match::StillSingle[0][]'~exp(othz.actual')~Match::StillSingle[1][]');
	println("Equibrium mnQ ",AGGTauchen::mnQ);
	println("Equilibrium Tax Rate ",TaxRate,", Avg Income ",AVGINCOME);
	viter.Volume = LOUD;
	decl sim = new Panel ( 0 , viter , TRUE);
	sim->Simulate(100,T,0,FALSE);
	sim->Print("simdata.dta");
	}	
	

/** .**/
AGG::Predict(ps,tod){
	decl hold;
	Bellman::Predict(ps,tod);
	if (curt==3) eqqty[EQMatch][CV(othz)] += ps*(CV(status)==single);
	if (!IsChild() && IsNest()) { //only do this for labour supply states
		if (CV(g)==female) eqqty[EQEndow][0] += ps*(pandv[0]'*exp(Q));
		eqqty[EQEndow][1] += ps*sumc(pandv[0].*(grossinc~welfrcpt)); 
		}
	}

/** Static function that is assigned to DP::PostRESolve **/
AGG::EQ() {
	if (!last[EQCouple]) return;
	decl p = new PanelPrediction(0),cg= CV(g),tot;
	eqqty[EQMatch] = zeros(1,othz.N);
	if (cg==female) {
		eqqty[EQEndow][0] = 0.0;
		eqqty[EQEndow][1] = zeros(1,2);
		}
	p->Predict(T-1);
	tot = double(sumr(eqqty[EQMatch]));
	eqqty[EQMatch] /= tot;
	eqdist[EQMatch] += norm(eqqty[EQMatch]-Match::StillSingle[!CV(g)][],2);
	Match::Previous[!CV(g)][] = eqqty[EQMatch];
	if (last[EQMatch]) {
		decl reqrt = eqqty[EQEndow][1][0][1]/eqqty[EQEndow][1][0][0];
		if (CV(g)==female) AGGTauchen::PrevQ = eqqty[EQEndow][0];
		eqdist[EQEndow] += sqr(reqrt-TaxRate)+sqr(eqqty[EQEndow][0]-AGGTauchen::mnQ);
		TaxRate = reqrt;
		AVGINCOME = eqqty[EQEndow][1][0][0];
		if (last[EQEndow]) {
			SubToPrint=nest;
			p->Histogram(status,ToScreen,TRUE);
			p->Histogram(welf,ToScreen,TRUE);
			p->Histogram(0,ToScreen,TRUE);
			SubToPrint=mingle;
			p->Histogram(ask,ToScreen,TRUE);
			p->Histogram(myz,ToScreen,TRUE);
			p->Histogram(othz,ToScreen,TRUE);
			}
		}
	delete p;
	}

AGG::Smooth(EV) {
	ExPostSmoothing::Smooth(EV);
	if (IsChild()) return;
	if (!IsNest() && CV(status)!=divorced) {
		decl h = Settheta(mtchind),myg=CV(g);
		eqdist[EQCouple] += sqr(h.saysyes[myg]-pandv[0][yes]);
		h.saysyes[myg] = pandv[0][yes];
		};
	}
	
AGG::OutputValue() {
	return exp(Q);
	}