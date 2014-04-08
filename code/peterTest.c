#include <stdio.h>
#include "cosmo.h"
#include "BEMRI.h"

int main()
{

// testing Peters' lifetime calculation

G = 6.67e-11;
c = 2.99e8;

// peters value 1
double T = peters_lifetime(1e-5, 10*RSUN, MSUN, MSUN);
printf("e = 1e-5, a = 10*RSUN, m1=m2=MSUN: T = %.5es = %.5e yrs\n",T,T/(365.0*24.0*3600.0));

// peters value 1
T = peters_lifetime(0, 1e8, MSUN, MSUN);
printf("e = 1e-5, a = 10^10cm, m1=m2=MSUN: T = %.5es = %.5e yrs\n",T,T/(365.0*24.0*3600.0));

// My sim values
T = peters_lifetime(0,10*RSUN,10*MSUN,10*MSUN);
printf("BEMRI binary: e = 0, a = 10*RSUN, m1=m2=10*MSUN : T = %.5e\n",T);
T = peters_lifetime(0.147,1.015*10*RSUN,10*MSUN,10*MSUN);
printf("BEMRI binary: e = 0.147, a = 1.1015*10*RSUN, m1=m2=10*MSUN : T = %.5e\n",T);

return 0;
}
