/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(BOOSTEDISOTROPICKERR_HPP_)
#error "This file should only be included through BoostedIsotropicKerr.hpp"
#endif

#ifndef BOOSTEDISOTROPICKERR_IMPL_HPP_
#define BOOSTEDISOTROPICKERR_IMPL_HPP_

#include "DimensionDefinitions.hpp"



    template <class data_t>
    data_t BoostedIsotropicKerr::get_func(double M, double a, double V, data_t x, double y, double z, int index) const{
        const double t = 0.0;
        const double v2 = V * V;
        const double boost2 = 1.0 / (1 - v2);
        double boost = sqrt(boost2);
        double rp = M + sqrt(M * M - a * a);
        double rm = M - sqrt(M * M - a * a);

        // the isotropic radius (boosted)
        data_t x_boost = x * boost;
        data_t r2 = x_boost * x_boost + y * y + z * z;
        data_t r = sqrt(r2);
        data_t rBL = r * pow((1 + rp / (4.0 * r)),2);
        data_t d_rBL = 1.0 - rp * rp / (16. * r * r);

        if (index == 0){
            return pow(a,2)*pow(z,2) + pow(r,2)*pow(rBL,2);
        }

        if (index == 1){
            return pow(z,2)/(pow(r,6)*(pow(x_boost,2) + pow(y,2)));
        }

        if (index == 2){
            return pow(r + rp/4.,2)/(pow(pow(r,2),2.5)*(-rm + rBL));
        }

        if (index == 3){
            return (pow(pow(a,2) + pow(rBL,2),2) - (pow(a,2)*(pow(x_boost,2) + pow(y,2))*(pow(a,2) - 2*M*rBL + pow(rBL,2)))/pow(r,2))/(pow(x_boost,2) + pow(y,2));
        }
    }

    template <class data_t>
    data_t BoostedIsotropicKerr::get_d_func(double M, double a, double V, data_t x, double y, double z, int index, int d_index) const{
        const double t = 0.0;
        const double v2 = V * V;
        const double boost2 = 1.0 / (1 - v2);
        double boost = sqrt(boost2);
        double rp = M + sqrt(M * M - a * a);
        double rm = M - sqrt(M * M - a * a);

        // the isotropic radius (boosted)
        data_t x_boost = x * boost;
        data_t r2 = x_boost * x_boost + y * y + z * z;
        data_t r = sqrt(r2);
        data_t rBL = r * pow((1 + rp / (4.0 * r)),2);
        data_t d_rBL = 1.0 - rp * rp / (16. * r * r);

        if (index == 0){
            if (d_index == 1) {
                return 2*x_boost*rBL*(rBL + r*d_rBL);
            }

            if (d_index == 2) {
                return 2*y*rBL*(rBL + r*d_rBL);
            }

            if (d_index == 3) {
                return 2*z*(pow(a,2) + pow(rBL,2) + r*rBL*d_rBL);
            }
        }

        if (index == 1){
            if (d_index == 1) {
                return (-2*x_boost*pow(z,2)*(4*pow(x_boost,2) + 4*pow(y,2) + pow(z,2)))/(pow(r,8)*pow(pow(x_boost,2) + pow(y,2),2));
            }

            if (d_index == 2) {
                return (-2*y*pow(z,2)*(4*pow(x_boost,2) + 4*pow(y,2) + pow(z,2)))/(pow(r,8)*pow(pow(x_boost,2) + pow(y,2),2));
            }

            if (d_index == 3) {
                return (2*z*(pow(x_boost,2) + pow(y,2) - 2*pow(z,2)))/(pow(r,8)*(pow(x_boost,2) + pow(y,2)));
            }
        }

        if (index == 2){
            if (d_index == 1) {
                return -((4*r + rp)*x_boost*(-(rm*(12*r + 5*rp)) + (12*r + 5*rp)*rBL + (r*rp + 4*pow(x_boost,2) + 4*pow(y,2) + 4*pow(z,2))*d_rBL))/(16.*pow(pow(r,2),3.5)*pow(rm - rBL,2));
            }

            if (d_index == 2) {
                return -((4*r + rp)*y*(-(rm*(12*r + 5*rp)) + (12*r + 5*rp)*rBL + (r*rp + 4*pow(x_boost,2) + 4*pow(y,2) + 4*pow(z,2))*d_rBL))/(16.*pow(pow(r,2),3.5)*pow(rm - rBL,2));
            }

            if (d_index == 3) {
                return -((4*r + rp)*z*(-(rm*(12*r + 5*rp)) + (12*r + 5*rp)*rBL + (r*rp + 4*pow(x_boost,2) + 4*pow(y,2) + 4*pow(z,2))*d_rBL))/(16.*pow(pow(r,2),3.5)*pow(rm - rBL,2));
            }
        }

        if (index == 3){
            if (d_index == 1) {
                return (-2*(pow(r,2)*x_boost*(pow(r,2)*pow(pow(a,2) + pow(rBL,2),2) - pow(a,2)*(pow(x_boost,2) + pow(y,2))*(pow(a,2) - 2*M*rBL + pow(rBL,2))) - x_boost*(pow(x_boost,2) + pow(y,2))*(-(pow(a,4)*pow(z,2)) - pow(a,2)*pow(z,2)*pow(rBL,2) + pow(a,2)*M*r*(pow(x_boost,2) + pow(y,2))*d_rBL + 2*pow(pow(r,2),1.5)*pow(rBL,3)*d_rBL + pow(a,2)*rBL*(2*M*pow(z,2) + r*(pow(x_boost,2) + pow(y,2) + 2*pow(z,2))*d_rBL))))/(pow(r,4)*pow(pow(x_boost,2) + pow(y,2),2));
            }

            if (d_index == 2) {
                return (-2*(pow(r,2)*y*(pow(r,2)*pow(pow(a,2) + pow(rBL,2),2) - pow(a,2)*(pow(x_boost,2) + pow(y,2))*(pow(a,2) - 2*M*rBL + pow(rBL,2))) - y*(pow(x_boost,2) + pow(y,2))*(-(pow(a,4)*pow(z,2)) - pow(a,2)*pow(z,2)*pow(rBL,2) + pow(a,2)*M*r*(pow(x_boost,2) + pow(y,2))*d_rBL + 2*pow(pow(r,2),1.5)*pow(rBL,3)*d_rBL + pow(a,2)*rBL*(2*M*pow(z,2) + r*(pow(x_boost,2) + pow(y,2) + 2*pow(z,2))*d_rBL))))/(pow(r,4)*pow(pow(x_boost,2) + pow(y,2),2));
            }

            if (d_index == 3) {
                return 2*z*((pow(a,2)*(pow(a,2) - 2*M*rBL + pow(rBL,2)))/pow(r,4) + (pow(a,2)*(M - rBL)*d_rBL)/pow(pow(r,2),1.5) + (2*rBL*(pow(a,2) + pow(rBL,2))*d_rBL)/(sqrt(pow(r,2))*(pow(x_boost,2) + pow(y,2))));
            }
        }
    }

    template <class data_t>
    data_t BoostedIsotropicKerr::get_alpha(double M, double a, double V, data_t x, double y, double z) const {

        const double t = 0.0;
        const double v2 = V * V;
        const double boost2 = 1.0 / (1 - v2);
        double boost = sqrt(boost2);
        double rp = M + sqrt(M * M - a * a);
        double rm = M - sqrt(M * M - a * a);

        // the isotropic radius (boosted)
        data_t x_boost = x * boost;
        data_t r2 = x_boost * x_boost + y * y + z * z;
        data_t r = sqrt(r2);
        data_t rBL = r * pow((1 + rp / (4.0 * r)),2);
        data_t d_rBL = 1.0 - rp * rp / (16. * r * r);
        data_t rb = r;
 

        //Helper functions:
        data_t func0 = get_func(M, a, V, x, y, z, 0);
        data_t func1 = get_func(M, a, V, x, y, z, 1);
        data_t func2 = get_func(M, a, V, x, y, z, 2);
        data_t func3 = get_func(M, a, V, x, y, z, 3);



        data_t dx_r = x_boost/r; data_t dy_r = y/r; data_t dz_r = z/r;

        data_t dx_func0 = get_d_func(M, a, V, x, y, z, 0, 1);
        data_t dy_func0 = get_d_func(M, a, V, x, y, z, 0, 2);
        data_t dz_func0 = get_d_func(M, a, V, x, y, z, 0, 3);
        data_t dx_func1 = get_d_func(M, a, V, x, y, z, 1, 1);
        data_t dy_func1 = get_d_func(M, a, V, x, y, z, 1, 2);
        data_t dz_func1 = get_d_func(M, a, V, x, y, z, 1, 3);
        data_t dx_func2 = get_d_func(M, a, V, x, y, z, 2, 1);
        data_t dy_func2 = get_d_func(M, a, V, x, y, z, 2, 2);
        data_t dz_func2 = get_d_func(M, a, V, x, y, z, 2, 3);
        data_t dx_func3 = get_d_func(M, a, V, x, y, z, 3, 1);
        data_t dy_func3 = get_d_func(M, a, V, x, y, z, 3, 2);
        data_t dz_func3 = get_d_func(M, a, V, x, y, z, 3, 3);
        
        data_t alpha;
        #include "BoostedIsotropicKerrBHFixedBG_Coutput/alpha.txt"
        return alpha;
    }

    template<class data_t>
    Tensor<1,data_t> BoostedIsotropicKerr::get_beta(double M, double a, double V, data_t x, double y, double z) const {

        const double t = 0.0;
        const double v2 = V * V;
        const double boost2 = 1.0 / (1 - v2);
        double boost = sqrt(boost2);
        double rp = M + sqrt(M * M - a * a);
        double rm = M - sqrt(M * M - a * a);

        // the isotropic radius (boosted)
        data_t x_boost = x * boost;
        data_t r2 = x_boost * x_boost + y * y + z * z;
        data_t r = sqrt(r2);
        data_t rBL = r * pow((1 + rp / (4.0 * r)),2);
        data_t d_rBL = 1.0 - rp * rp / (16. * r * r);
        data_t rb = r;
 
        
        //Helper functions:
        data_t func0 = get_func(M, a, V, x, y, z, 0);
        data_t func1 = get_func(M, a, V, x, y, z, 1);
        data_t func2 = get_func(M, a, V, x, y, z, 2);
        data_t func3 = get_func(M, a, V, x, y, z, 3);

        data_t dx_r = x_boost/r; data_t dy_r = y/r; data_t dz_r = z/r;

        data_t dx_func0 = get_d_func(M, a, V, x, y, z, 0, 1);
        data_t dy_func0 = get_d_func(M, a, V, x, y, z, 0, 2);
        data_t dz_func0 = get_d_func(M, a, V, x, y, z, 0, 3);
        data_t dx_func1 = get_d_func(M, a, V, x, y, z, 1, 1);
        data_t dy_func1 = get_d_func(M, a, V, x, y, z, 1, 2);
        data_t dz_func1 = get_d_func(M, a, V, x, y, z, 1, 3);
        data_t dx_func2 = get_d_func(M, a, V, x, y, z, 2, 1);
        data_t dy_func2 = get_d_func(M, a, V, x, y, z, 2, 2);
        data_t dz_func2 = get_d_func(M, a, V, x, y, z, 2, 3);
        data_t dx_func3 = get_d_func(M, a, V, x, y, z, 3, 1);
        data_t dy_func3 = get_d_func(M, a, V, x, y, z, 3, 2);
        data_t dz_func3 = get_d_func(M, a, V, x, y, z, 3, 3);
        

        Tensor<1,data_t> betaU;
        #include "BoostedIsotropicKerrBHFixedBG_Coutput/beta0.txt"
        #include "BoostedIsotropicKerrBHFixedBG_Coutput/beta1.txt"
        #include "BoostedIsotropicKerrBHFixedBG_Coutput/beta2.txt"
        return betaU;
    }

    template<class data_t>
    Tensor<2,data_t> BoostedIsotropicKerr::get_gamma(double M, double a, double V, data_t x, double y, double z) const {

        const double t = 0.0;
        const double v2 = V * V;
        const double boost2 = 1.0 / (1 - v2);
        double boost = sqrt(boost2);
        double rp = M + sqrt(M * M - a * a);
        double rm = M - sqrt(M * M - a * a);

        // the isotropic radius (boosted)
        data_t x_boost = x * boost;
        data_t r2 = x_boost * x_boost + y * y + z * z;
        data_t r = sqrt(r2);
        data_t rBL = r * pow((1 + rp / (4.0 * r)),2);
        data_t d_rBL = 1.0 - rp * rp / (16. * r * r);
        data_t rb = r;
 
        
        //Helper functions:
        data_t func0 = get_func(M, a, V, x, y, z, 0);
        data_t func1 = get_func(M, a, V, x, y, z, 1);
        data_t func2 = get_func(M, a, V, x, y, z, 2);
        data_t func3 = get_func(M, a, V, x, y, z, 3);

        data_t dx_r = x_boost/r; data_t dy_r = y/r; data_t dz_r = z/r;

        data_t dx_func0 = get_d_func(M, a, V, x, y, z, 0, 1);
        data_t dy_func0 = get_d_func(M, a, V, x, y, z, 0, 2);
        data_t dz_func0 = get_d_func(M, a, V, x, y, z, 0, 3);
        data_t dx_func1 = get_d_func(M, a, V, x, y, z, 1, 1);
        data_t dy_func1 = get_d_func(M, a, V, x, y, z, 1, 2);
        data_t dz_func1 = get_d_func(M, a, V, x, y, z, 1, 3);
        data_t dx_func2 = get_d_func(M, a, V, x, y, z, 2, 1);
        data_t dy_func2 = get_d_func(M, a, V, x, y, z, 2, 2);
        data_t dz_func2 = get_d_func(M, a, V, x, y, z, 2, 3);
        data_t dx_func3 = get_d_func(M, a, V, x, y, z, 3, 1);
        data_t dy_func3 = get_d_func(M, a, V, x, y, z, 3, 2);
        data_t dz_func3 = get_d_func(M, a, V, x, y, z, 3, 3);
        

        Tensor<2,data_t> gammaLL;
        #include "BoostedIsotropicKerrBHFixedBG_Coutput/gamma00.txt"
        #include "BoostedIsotropicKerrBHFixedBG_Coutput/gamma01.txt"
        #include "BoostedIsotropicKerrBHFixedBG_Coutput/gamma02.txt"
        #include "BoostedIsotropicKerrBHFixedBG_Coutput/gamma11.txt"
        #include "BoostedIsotropicKerrBHFixedBG_Coutput/gamma12.txt"
        #include "BoostedIsotropicKerrBHFixedBG_Coutput/gamma22.txt"
        gammaLL[1][0] = gammaLL[0][1];
        gammaLL[2][0] = gammaLL[0][2];
        gammaLL[2][1] = gammaLL[1][2];
        return gammaLL;
    }

    /*template<class data_t>
    Tensor<1,data_t> BoostedIsotropicKerr::get_d1alpha(double M, double a, double V, data_t x, double y, double z) const {

        const double t = 0.0;
        const double v2 = V * V;
        const double boost2 = 1.0 / (1 - v2);
        double boost = sqrt(boost2);
        double rp = M + sqrt(M * M - a * a);
        double rm = M - sqrt(M * M - a * a);

        // the isotropic radius (boosted)
        data_t x_boost = x * boost;
        data_t r2 = x_boost * x_boost + y * y + z * z;
        data_t r = sqrt(r2);
        data_t rBL = r * pow((1 + rp / (4.0 * r)),2);
        data_t d_rBL = 1.0 - rp * rp / (16. * r * r);
        data_t rb = r;
 
        
        //Helper functions:
        data_t func0 = get_func(M, a, V, x, y, z, 0);
        data_t func1 = get_func(M, a, V, x, y, z, 1);
        data_t func2 = get_func(M, a, V, x, y, z, 2);
        data_t func3 = get_func(M, a, V, x, y, z, 3);

        data_t dx_r = x_boost/r; data_t dy_r = y/r; data_t dz_r = z/r;

        data_t dx_func0 = get_d_func(M, a, V, x, y, z, 0, 1);
        data_t dy_func0 = get_d_func(M, a, V, x, y, z, 0, 2);
        data_t dz_func0 = get_d_func(M, a, V, x, y, z, 0, 3);
        data_t dx_func1 = get_d_func(M, a, V, x, y, z, 1, 1);
        data_t dy_func1 = get_d_func(M, a, V, x, y, z, 1, 2);
        data_t dz_func1 = get_d_func(M, a, V, x, y, z, 1, 3);
        data_t dx_func2 = get_d_func(M, a, V, x, y, z, 2, 1);
        data_t dy_func2 = get_d_func(M, a, V, x, y, z, 2, 2);
        data_t dz_func2 = get_d_func(M, a, V, x, y, z, 2, 3);
        data_t dx_func3 = get_d_func(M, a, V, x, y, z, 3, 1);
        data_t dy_func3 = get_d_func(M, a, V, x, y, z, 3, 2);
        data_t dz_func3 = get_d_func(M, a, V, x, y, z, 3, 3);
        

        Tensor<1,data_t> d1alpha;
        #include "BoostedIsotropicKerrBHFixedBG_Coutput/d1alpha0.txt"
        #include "BoostedIsotropicKerrBHFixedBG_Coutput/d1alpha1.txt"
        #include "BoostedIsotropicKerrBHFixedBG_Coutput/d1alpha2.txt"


        return d1alpha;
    }

    template<class data_t>
    Tensor<2,data_t> BoostedIsotropicKerr::get_d1beta(double M, double a, double V, data_t x, double y, double z) const {

        const double t = 0.0;
        const double v2 = V * V;
        const double boost2 = 1.0 / (1 - v2);
        double boost = sqrt(boost2);
        double rp = M + sqrt(M * M - a * a);
        double rm = M - sqrt(M * M - a * a);

        // the isotropic radius (boosted)
        data_t x_boost = x * boost;
        data_t r2 = x_boost * x_boost + y * y + z * z;
        data_t r = sqrt(r2);
        data_t rBL = r * pow((1 + rp / (4.0 * r)),2);
        data_t d_rBL = 1.0 - rp * rp / (16. * r * r);
        data_t rb = r;

 
        
        //Helper functions:
        data_t func0 = get_func(M, a, V, x, y, z, 0);
        data_t func1 = get_func(M, a, V, x, y, z, 1);
        data_t func2 = get_func(M, a, V, x, y, z, 2);
        data_t func3 = get_func(M, a, V, x, y, z, 3);

        data_t dx_r = x_boost/r; data_t dy_r = y/r; data_t dz_r = z/r;

        data_t dx_func0 = get_d_func(M, a, V, x, y, z, 0, 1);
        data_t dy_func0 = get_d_func(M, a, V, x, y, z, 0, 2);
        data_t dz_func0 = get_d_func(M, a, V, x, y, z, 0, 3);
        data_t dx_func1 = get_d_func(M, a, V, x, y, z, 1, 1);
        data_t dy_func1 = get_d_func(M, a, V, x, y, z, 1, 2);
        data_t dz_func1 = get_d_func(M, a, V, x, y, z, 1, 3);
        data_t dx_func2 = get_d_func(M, a, V, x, y, z, 2, 1);
        data_t dy_func2 = get_d_func(M, a, V, x, y, z, 2, 2);
        data_t dz_func2 = get_d_func(M, a, V, x, y, z, 2, 3);
        data_t dx_func3 = get_d_func(M, a, V, x, y, z, 3, 1);
        data_t dy_func3 = get_d_func(M, a, V, x, y, z, 3, 2);
        data_t dz_func3 = get_d_func(M, a, V, x, y, z, 3, 3);
        

        Tensor<2,data_t> d1betaU;
        #include "BoostedIsotropicKerrBHFixedBG_Coutput/d1betaU00.txt"
        #include "BoostedIsotropicKerrBHFixedBG_Coutput/d1betaU01.txt"
        #include "BoostedIsotropicKerrBHFixedBG_Coutput/d1betaU02.txt"
        #include "BoostedIsotropicKerrBHFixedBG_Coutput/d1betaU10.txt"
        #include "BoostedIsotropicKerrBHFixedBG_Coutput/d1betaU11.txt"
        #include "BoostedIsotropicKerrBHFixedBG_Coutput/d1betaU12.txt"
        #include "BoostedIsotropicKerrBHFixedBG_Coutput/d1betaU20.txt"
        #include "BoostedIsotropicKerrBHFixedBG_Coutput/d1betaU21.txt"
        #include "BoostedIsotropicKerrBHFixedBG_Coutput/d1betaU22.txt"



        return d1betaU;
    }

    template<class data_t>
    Tensor<3,data_t> BoostedIsotropicKerr::get_d1gamma(double M, double a, double V, data_t x, double y, double z) const {

        const double t = 0.0;
        const double v2 = V * V;
        const double boost2 = 1.0 / (1 - v2);
        double boost = sqrt(boost2);
        double rp = M + sqrt(M * M - a * a);
        double rm = M - sqrt(M * M - a * a);

        // the isotropic radius (boosted)
        data_t x_boost = x * boost;
        data_t r2 = x_boost * x_boost + y * y + z * z;
        data_t r = sqrt(r2);
        data_t rBL = r * pow((1 + rp / (4.0 * r)),2);
        data_t d_rBL = 1.0 - rp * rp / (16. * r * r);
        data_t rb = r;
 
        
        //Helper functions:
        data_t func0 = get_func(M, a, V, x, y, z, 0);
        data_t func1 = get_func(M, a, V, x, y, z, 1);
        data_t func2 = get_func(M, a, V, x, y, z, 2);
        data_t func3 = get_func(M, a, V, x, y, z, 3);

        data_t dx_r = x_boost/r; data_t dy_r = y/r; data_t dz_r = z/r;

        data_t dx_func0 = get_d_func(M, a, V, x, y, z, 0, 1);
        data_t dy_func0 = get_d_func(M, a, V, x, y, z, 0, 2);
        data_t dz_func0 = get_d_func(M, a, V, x, y, z, 0, 3);
        data_t dx_func1 = get_d_func(M, a, V, x, y, z, 1, 1);
        data_t dy_func1 = get_d_func(M, a, V, x, y, z, 1, 2);
        data_t dz_func1 = get_d_func(M, a, V, x, y, z, 1, 3);
        data_t dx_func2 = get_d_func(M, a, V, x, y, z, 2, 1);
        data_t dy_func2 = get_d_func(M, a, V, x, y, z, 2, 2);
        data_t dz_func2 = get_d_func(M, a, V, x, y, z, 2, 3);
        data_t dx_func3 = get_d_func(M, a, V, x, y, z, 3, 1);
        data_t dy_func3 = get_d_func(M, a, V, x, y, z, 3, 2);
        data_t dz_func3 = get_d_func(M, a, V, x, y, z, 3, 3);
        

        Tensor<3,data_t> d1gammaLL;
        #include "BoostedIsotropicKerrBHFixedBG_Coutput/d1gammaLL000.txt"
        #include "BoostedIsotropicKerrBHFixedBG_Coutput/d1gammaLL001.txt"
        #include "BoostedIsotropicKerrBHFixedBG_Coutput/d1gammaLL002.txt"
        #include "BoostedIsotropicKerrBHFixedBG_Coutput/d1gammaLL010.txt"
        #include "BoostedIsotropicKerrBHFixedBG_Coutput/d1gammaLL011.txt"
        #include "BoostedIsotropicKerrBHFixedBG_Coutput/d1gammaLL012.txt"
        #include "BoostedIsotropicKerrBHFixedBG_Coutput/d1gammaLL020.txt"
        #include "BoostedIsotropicKerrBHFixedBG_Coutput/d1gammaLL021.txt"
        #include "BoostedIsotropicKerrBHFixedBG_Coutput/d1gammaLL022.txt"
        FOR1(i){
            d1gammaLL[1][0][i] = d1gammaLL[0][1][i];
            d1gammaLL[2][0][i] = d1gammaLL[0][2][i];
            d1gammaLL[2][1][i] = d1gammaLL[1][2][i];
        }

        return d1gammaLL;
    }
    */
    template<class data_t>
    const Tensor<2,data_t> BoostedIsotropicKerr::get_K(double M, double a, double V, data_t x, double y, double z) const {

        const double t = 0.0;
        const double v2 = V * V;
        const double boost2 = 1.0 / (1 - v2);
        double boost = sqrt(boost2);
        double rp = M + sqrt(M * M - a * a);
        double rm = M - sqrt(M * M - a * a);

        // the isotropic radius (boosted)
        data_t x_boost = x * boost;
        data_t r2 = x_boost * x_boost + y * y + z * z;
        data_t r = sqrt(r2);
        data_t rBL = r * pow((1 + rp / (4.0 * r)),2);
        data_t d_rBL = 1.0 - rp * rp / (16. * r * r);
        data_t rb = r;
 
        
        //Helper functions:
        data_t func0 = get_func(M, a, V, x, y, z, 0);
        data_t func1 = get_func(M, a, V, x, y, z, 1);
        data_t func2 = get_func(M, a, V, x, y, z, 2);
        data_t func3 = get_func(M, a, V, x, y, z, 3);

        data_t dx_r = x_boost/r; data_t dy_r = y/r; data_t dz_r = z/r;

        data_t dx_func0 = get_d_func(M, a, V, x, y, z, 0, 1);
        data_t dy_func0 = get_d_func(M, a, V, x, y, z, 0, 2);
        data_t dz_func0 = get_d_func(M, a, V, x, y, z, 0, 3);
        data_t dx_func1 = get_d_func(M, a, V, x, y, z, 1, 1);
        data_t dy_func1 = get_d_func(M, a, V, x, y, z, 1, 2);
        data_t dz_func1 = get_d_func(M, a, V, x, y, z, 1, 3);
        data_t dx_func2 = get_d_func(M, a, V, x, y, z, 2, 1);
        data_t dy_func2 = get_d_func(M, a, V, x, y, z, 2, 2);
        data_t dz_func2 = get_d_func(M, a, V, x, y, z, 2, 3);
        data_t dx_func3 = get_d_func(M, a, V, x, y, z, 3, 1);
        data_t dy_func3 = get_d_func(M, a, V, x, y, z, 3, 2);
        data_t dz_func3 = get_d_func(M, a, V, x, y, z, 3, 3);
        


        Tensor<2,data_t> K;
        #include "BoostedIsotropicKerrBHFixedBG_Coutput/K00.txt"
        #include "BoostedIsotropicKerrBHFixedBG_Coutput/K01.txt"
        #include "BoostedIsotropicKerrBHFixedBG_Coutput/K02.txt"
        #include "BoostedIsotropicKerrBHFixedBG_Coutput/K11.txt"
        #include "BoostedIsotropicKerrBHFixedBG_Coutput/K12.txt"
        #include "BoostedIsotropicKerrBHFixedBG_Coutput/K22.txt"
        K[1][0] = K[0][1];
        K[2][0] = K[0][2];
        K[2][1] = K[1][2];
        return K;
    }

    /// This just calculates chi which helps with regridding, debug etc
    /// it is only done once on setup as the BG is fixed
    template <class data_t> void BoostedIsotropicKerr::compute(Cell<data_t> current_cell) const
    {
        // get position and set vars

        Vars<data_t> vars;
        compute_metric_background(vars, current_cell);
        current_cell.store_vars(vars);
    }

    /// Schwarzschild boosted solution as above
    template <class data_t, template <typename> class vars_t>
    void BoostedIsotropicKerr::compute_metric_background(vars_t<data_t> &vars,
                                   const Cell<data_t> &current_cell) const
    {
        // where am i?
        const Coordinates<data_t> coords(current_cell, m_dx, m_params.center);

        // black hole params - mass M and boost v
        // "boost" is the gamma factor for the boost
        const double M = m_params.mass;
        const double V = m_params.velocity;
        const double a = m_params.spin;
        //std::cout << "V: " << V << std::endl;
        // work out where we are on the grid including effect of boost
        // on x direction (length contraction)
        //const double t = 0.0;
        data_t x = coords.x;
        double y = coords.y;
        double z = coords.z;

        //DEBUG!!!
        //data_t xnum = 0.8; double ynum = 0.4; double znum = 0.5;
        // populate ADM vars
        const data_t alpha = get_alpha(M, a, V, x, y, z);
        Tensor<1,data_t> betaU = get_beta(M, a, V, x, y, z);
        Tensor<2,data_t> gamma = get_gamma(M, a, V, x, y, z);
        Tensor<2,data_t> K = get_K(M, a, V, x, y, z);
        
        vars.lapse = alpha;
        FOR1(i){ vars.shift[i] = betaU[i]; }     
        FOR2(i,j){ vars.A[i][j] = K[i][j]; }

        using namespace TensorAlgebra;

        const auto gamma_UU = compute_inverse_sym(gamma);
                
        vars.K = compute_trace(K, gamma_UU );

        data_t det_gamma = compute_determinant_sym(gamma);
        vars.chi = pow(det_gamma, -1.0 / 3.0);

        vars.K = compute_trace(vars.A, gamma_UU);
        make_trace_free(vars.A, gamma, gamma_UU);
        FOR(i, j)
        {
            vars.h[i][j] = gamma[i][j] *  vars.chi;
            vars.A[i][j] *= vars.chi;
        }
        /*
        Tensor<1,data_t> d1alpha = get_d1alpha(M, a, V, x, y, z);
        Tensor<2,data_t> d1betaU = get_d1beta(M, a, V, x, y, z);
        Tensor<3,data_t> d1gamma = get_d1gamma(M, a, V, x, y, z);
    

        // Calculate partial derivative of spatial metric
        FOR3(i, j, k) { vars.d1_gamma[i][j][k] = d1gamma[i][j][k]; }

        // calculate derivs of lapse and shift
        FOR1(i)
        {
            vars.d1_lapse[i] = d1alpha[i];
        }

        FOR2(i, j)
        {
            vars.d1_shift[i][j] = d1betaU[i][j];
        }
        */
        
    }

    // used to decide when to excise - ie when within the horizon of the BH
    // note that this is not templated over data_t
    /*double BoostedIsotropicKerr::excise(const Cell<double> &current_cell) const
    {
        // black hole params - mass M and boost v
        // "boost" is the gamma factor for the boost
        const Coordinates<double> coords(current_cell, m_dx, m_params.center);
        const double M = m_params.mass;
        const double a = m_params.spin;
        const double boost =
            pow(1.0 - m_params.velocity * m_params.velocity, -0.5);

        // work out where we are on the grid including effect of boost
        // on x direction (length contraction)
        const double x_p = coords.x * boost;
        const double y = coords.y;
        const double z = coords.z;

        // the coordinate radius (boosted)
        const double r2 = x_p * x_p + y * y + z * z;

        // compare this to horizon in isotropic coords
        double rp = M + sqrt(M * M - a * a);
        const double r_horizon = rp;

        return sqrt(r2) / r_horizon;
    }*/

#endif /* BOOSTEDISOTROPICKERRBHFIXEDBG_HPP_ */
