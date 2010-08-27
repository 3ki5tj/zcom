// Calculate the exact Partition function, energy, and heat capacity
// for a finite 2d Ising model on square lattice
#include <stdio.h>
#include <math.h>
#include <assert.h>

#define LOG_BIG 50.0
// c = log(exp(a)+exp(b))
static double log_add_(double a, double b){
	double c;
	if(a<b){c=a;a=b;b=c;} // now a >= b
	if((c=a-b)>LOG_BIG) return a;
	else return a+log(1+exp(-c));
}

// c = log(exp(a)+b);
static double log_addn_(double a, double b){
	if(a>LOG_BIG) return a;
	else return a+log(1+b*exp(-a));
}

// c = log(exp(a)-exp(b)), only works for a>b
static double log_min_(double a, double b){
	double c;
	assert(a>=b);
	if((c=a-b)>LOG_BIG) return a;
	else return a+log(1-exp(-c));
}

double Ising2DlnZ(double beta, double *E, double *Cv, int m, int n){
	double lnZ, ex, factor, gamma, gamma0, th,sech;
	double lnZ1, lnZ2, lnZ3, lnZ4;
	double Z21,Z31,Z41;// Z21= Z2/Z1
	double dr1, dr2, dr3, dr4;
	double ddr1, ddr2, ddr3, ddr4;
	double gamder, gamdder, gamder0;
	double exp2b,exp2bi,sinh2b,coth2b;
	double lncosh2b, lncc2b, lncl, lnsl, clderfactor, lncldder;
	double Zder, Zdder;
	//double cosh2b,sl, cc2b, cl, clder;
	int r;
	int sgn4=1;

	lnZ1=lnZ2=lnZ3=lnZ4=0;
	dr1=dr2=dr3=dr4=0;
	ddr1=ddr2=ddr3=ddr4=0;

	exp2b=exp(2*beta);
	exp2bi=1.0/exp2b;
	//cosh2b=0.5*(exp2b+exp2bi);
	lncosh2b=log_add_(2*beta, -2*beta)-log(2);
	coth2b=(1+exp2bi*exp2bi)/(1-exp2bi*exp2bi);
	//cc2b=cosh2b*coth2b; // this is always greater than 1
	lncc2b=lncosh2b+log(coth2b);
	gamma0=2*beta+log(1-exp2bi)-log(1+exp2bi);
	sgn4=((gamma0>=0)?1:-1) ;

	sinh2b=0.5*(exp2b-exp2bi);
	gamder0=2+2/sinh2b;
	// cl'=2*cosh2b*(1-1.0/sinh2b^2);
	clderfactor=2.0*(1-1.0/(sinh2b*sinh2b));
	// cl''=4.0*(cosh2b^2/sinh2b+2/sinh2b^3);
	lncldder=log_addn_(lncc2b, 2.0/(sinh2b*sinh2b*sinh2b))+log(4.0); // = log(cl'')
	//lnclder=log_add_(2*beta, -2*beta)+log(fabs(1-1.0/(sinh2b*sinh2b) ));

	for(r=0; r<=n-1; r++){
		// for odd number
		// cl=cc2b-cos((2*r+1)*M_PI/n);
		lncl = log_addn_(lncc2b, -cos((2*r+1)*M_PI/n));
		// sl=sqrt(cl*cl-1);
		lnsl=lncl+0.5*log(1-exp(-2*lncl));
		// gamma=log(cl+sl);
		gamma=log_add_(lncl, lnsl);
		if(gamma<=0){
			printf("Ising2DlnZ fail at Beta=%g\n", beta);
			*E=*Cv=0;
			return 0;
		}
		factor = 0.5*m*gamma;
		lnZ1+=log_add_(factor, -factor); // log(exp(factor)+exp(-factor)) = log 2 cosh(factor)
		lnZ2+=log_min_(factor, -factor);

		ex=exp(-factor);
		th=(1-ex*ex)/(1+ex*ex);
		// gamma'=cl'/sl;
		// we split cl' because it has a large part, yet, it is not always positive
		gamder=exp(lncosh2b-lnsl)*clderfactor;
		dr1+=0.5*m*gamder*th;
		dr2+=0.5*m*gamder/th;

		//gamma''=cl''/sl - cl' ^2 *cl/sl^3;
		gamdder=exp(lncldder-lnsl)-exp(lncosh2b*2+lncl-3*lnsl)*(clderfactor*clderfactor);
		sech=2.0*gamder/(ex+1.0/ex); // gamma' * sech(0.5*m*gamma)
		ddr1+=0.5*m*(gamdder*th+0.5*m*(sech*sech));
		sech=2.0*gamder/(ex-1.0/ex); // gamma' * sech(0.5*m*gamma)
		ddr2+=0.5*m*(gamdder/th-0.5*m*(sech*sech));


		// for even numbers
		if(r==0){
			gamma=gamma0;
		}else{
			// cl=cc2b-cos(2*r*M_PI/n);
			// sl=sqrt(cl*cl-1);
			// gamma=log(cl+sl);

			lncl = log_addn_(lncc2b, -cos((2*r)*M_PI/n));
			lnsl=lncl+0.5*log(1-exp(-2*lncl));
			gamma=log_add_(lncl, lnsl);
			assert(gamma>0);
		}
		factor = 0.5*m*gamma;
		lnZ3+=log_add_(factor, -factor); // log(exp(factor)+exp(-factor)) = log 2 cosh(factor)

		ex=exp(-factor);
		th=(1-ex*ex)/(1+ex*ex);
		if(factor<0)factor=-factor; // to avoid negative gamma0, at high T
		lnZ4+=log_min_(factor, -factor);

		if(r==0){
			gamder=gamder0;
		}else{
			// cl' = 2*cosh2b*(1-1.0/(sinh2b*sinh2b));
			// gamma' =cl'/sl;
			gamder=exp(lncosh2b-lnsl)*clderfactor;
		}

		dr3+=0.5*m*gamder*th;
		dr4+=0.5*m*gamder/th;

		if(r==0){
			gamdder=-4*coth2b*coth2b*exp(-lncosh2b);
		}else{
			// gam'' =cl''/sl - cl' ^2 *cl/sl^3;
			gamdder=exp(lncldder-lnsl)-exp(lncosh2b*2+lncl-3*lnsl)*(clderfactor*clderfactor);
		}
		sech=2.0*gamder/(ex+1.0/ex);
		ddr3+=0.5*m*(gamdder*th+0.5*m*(sech*sech));
		sech=2.0*gamder/(ex-1.0/ex);
		ddr4+=0.5*m*(gamdder/th-0.5*m*(sech*sech));
	}

	Z21=exp(lnZ2-lnZ1);
	Z31=exp(lnZ3-lnZ1);
	Z41=sgn4*exp(lnZ4-lnZ1);


	lnZ=lnZ1+log(1+Z21+Z31+Z41)+0.5*m*n*log(exp2b-exp2bi)-log(2);

	Zder=(dr1+Z21*dr2+Z31*dr3+Z41*dr4)/(1+Z21+Z31+Z41);
	*E=-m*n*coth2b-Zder;

	ddr1+=dr1*dr1;
	ddr2+=dr2*dr2;
	ddr3+=dr3*dr3;
	ddr4+=dr4*dr4;
	Zdder=(ddr1+Z21*ddr2+Z31*ddr3+Z41*ddr4)/(1+Z21+Z31+Z41);
	*Cv=beta*beta*(-2*m*n/(sinh2b*sinh2b)+Zdder-Zder*Zder);

	return lnZ;
}
/*
int main(void){
	double Beta=1.0/2.3, E, lnZ, Cv;
	lnZ=Ising2DlnZ(Beta, &E, &Cv, 8, 8);
	printf("%16.12f, %16.12f %16.12f\n", lnZ, E, Cv);
	return 0;
}
*/
