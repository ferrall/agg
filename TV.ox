#import "DDP"

main() {
	decl v, TV = new Tauchen("T",9,2,1.0,2.3,0.7);
	TV->Update();
	for (v=0;v<TV.N;++v) {
		TV.v = v;
		println(v," ",AV(TV),TV->Transit(<0>));
		}
	}