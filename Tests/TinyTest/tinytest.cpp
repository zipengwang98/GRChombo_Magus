#include <iostream>
#include <math.h>

int main(){
        double V;
        std::cin >> V;
        double a=5.0;
        double M=1.0;
        double x=1.0;
        double y=1.0;
        double z=1.0;

        const double t = 0.0;
        const double v2 = V * V;
        const double boost2 = 1.0 / (1 - v2);
        double boost = sqrt(boost2);
        double rp = M + sqrt(M * M - a * a);
        double rm = M - sqrt(M * M - a * a);

        // the isotropic radius (boosted)
        double x_boost = x * boost;
        double r2 = x_boost * x_boost + y * y + z * z;
        double r = sqrt(r2);
        double rBL = r * pow((1 + rp / (4.0 * r)),2);
        double d_rBL = 1.0 - rp * rp / (16. * r * r);
 
        
        //Helper functions:
        double func0 = 42;
        double func1 = 42;
        double func2 = 42;
        double func3 = 42;

        double dx_r = x_boost/r; double dy_r = y/r; double dz_r = z/r;

        double dx_func0 = 42;
        double dy_func0 = 42;
        double dz_func0 = 42;
        double dx_func1 = 42;
        double dy_func1 = 42;
        double dz_func1 = 42;
        double dx_func2 = 42;
        double dy_func2 = 42;
        double dz_func2 = 42;
        double dx_func3 = 42;
        double dy_func3 = 42;
        double dz_func3 = 42;

        double alpha;
        double d1betaU[3][3];
        double d1gammaLL[3][3][3];
        double K[3][3];

        #include "BoostedIsotropicKerrBHFixedBG_Coutput/alpha.txt"
        #include "BoostedIsotropicKerrBHFixedBG_Coutput/d1betaU00.txt"
        #include "BoostedIsotropicKerrBHFixedBG_Coutput/d1betaU01.txt"
        #include "BoostedIsotropicKerrBHFixedBG_Coutput/d1betaU02.txt"
        #include "BoostedIsotropicKerrBHFixedBG_Coutput/d1betaU10.txt"
        #include "BoostedIsotropicKerrBHFixedBG_Coutput/d1betaU11.txt"
        #include "BoostedIsotropicKerrBHFixedBG_Coutput/d1betaU12.txt"
        #include "BoostedIsotropicKerrBHFixedBG_Coutput/d1betaU20.txt"
        #include "BoostedIsotropicKerrBHFixedBG_Coutput/d1betaU21.txt"
        #include "BoostedIsotropicKerrBHFixedBG_Coutput/d1betaU22.txt"
        #include "BoostedIsotropicKerrBHFixedBG_Coutput/d1gammaLL000.txt"
        #include "BoostedIsotropicKerrBHFixedBG_Coutput/d1gammaLL001.txt"
        #include "BoostedIsotropicKerrBHFixedBG_Coutput/d1gammaLL002.txt"
        #include "BoostedIsotropicKerrBHFixedBG_Coutput/d1gammaLL010.txt"
        #include "BoostedIsotropicKerrBHFixedBG_Coutput/d1gammaLL011.txt"
        #include "BoostedIsotropicKerrBHFixedBG_Coutput/d1gammaLL012.txt"
        #include "BoostedIsotropicKerrBHFixedBG_Coutput/d1gammaLL020.txt"
        #include "BoostedIsotropicKerrBHFixedBG_Coutput/d1gammaLL021.txt"
        #include "BoostedIsotropicKerrBHFixedBG_Coutput/d1gammaLL022.txt"
        #include "BoostedIsotropicKerrBHFixedBG_Coutput/K00.txt"
        #include "BoostedIsotropicKerrBHFixedBG_Coutput/K01.txt"
        #include "BoostedIsotropicKerrBHFixedBG_Coutput/K02.txt"
        #include "BoostedIsotropicKerrBHFixedBG_Coutput/K11.txt"
        #include "BoostedIsotropicKerrBHFixedBG_Coutput/K12.txt"
        #include "BoostedIsotropicKerrBHFixedBG_Coutput/K22.txt"

        std::cout << "K:" << K[2][2] << std::endl;
    return 0;
}