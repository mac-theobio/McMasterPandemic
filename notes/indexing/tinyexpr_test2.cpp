// A more realistic example
// Here, we calcuate math expression:
//     (Ia) * (beta0) * (1/N) * (Ca) + (Ip) * (beta0) * (1/N) * (Cp)
// in two ways --- by direct calculation and by parser. Their results are the same.
// So, we conclude that tinyexpr works.
// NOTE: I have tested generic math expressions with log, sin, etc in another simple example.

#include <iostream>
#include <vector>
#include <map>
#include "tinyexpr.h"

int main()
{
    int error;

    // vector of state names (passed in from R side)
    std::vector<std::string> state_paramNames {"S", "E", "Ia", "Ip", "Im", "Is", "H", "H2", \
                                "ICUs", "ICUd", "D", "R", "X", "V", \
                                "beta0", "Ca", "Cp", "Cm", "Cs", \
                                "alpha", "sigma", "gamma_a", "gamma_m", "gamma_s", \
                                "gamma_p", "rho", "delta", "mu", "N", \
                                "E0", "nonhosp_mort", "iso_m", "iso_s", "phi1", \
                                "phi2", "psi1", "psi2", "psi3", "c_prop", \
                                "c_delay_mean", "c_delay_cv", "proc_disp", "zeta"};

    // vector of state and parameter variables (passed in from R side)
    // MUST be double not float
    double state_param [] =  {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, \
                             11, 12, 13, 14, 15, 16, 17, 18, 19, 10, 1, 2, 3, 4, \
                             5, 6, 7, 8, 9, 0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9};

    te_variable vars[14+29];

    // Build a map from name to index
    std::map<std::string, int> name2IndexMap;
    for (int i= 0; i< state_paramNames.size(); i++) {
        name2IndexMap[state_paramNames[i]] = i;
        vars[i] = { state_paramNames[i].c_str(), &state_param[i], 0, 0 };
        std::cout << "adding " << i << ", " << state_paramNames[i].c_str() \
                  << ", " << state_param[i] << std::endl;
    }

    std::string expr("(Ia) * (beta0) * (1/N) * (Ca) + (Ip) * (beta0) * (1/N) * (Cp)");
    std::cout << "Let me calculate " << expr << " ... " << std::endl;

    te_expr *n = te_compile(expr.c_str(), vars, 14+29, &error);
    std::cout << "calculating ... " << std::endl;

    // Check operands
    std::cout << "Ia    = " << state_param[name2IndexMap["Ia"]] << std::endl;
    std::cout << "beta0 = " << state_param[name2IndexMap["beta0"]] << std::endl;
    std::cout << "N     = " << state_param[name2IndexMap["N"]] << std::endl;
    std::cout << "Ca    = " << state_param[name2IndexMap["Ca"]] << std::endl;
    std::cout << "Ip    = " << state_param[name2IndexMap["Ip"]] << std::endl;
    std::cout << "Cp    = " << state_param[name2IndexMap["Cp"]] << std::endl;

    if (n)
        std::cout << "result by parser = " << te_eval(n) << std::endl;
    else
        std::cout << "Something wrong! Error code: " << error << std::endl;

    // The calculate the expression directly
    //float groudtruth = (Ia) * (beta0) * (1/N) * (Ca) + (Ip) * (beta0) * (1/N) * (Cp)
    float groudtruth = state_param[name2IndexMap["Ia"]] * 
                       state_param[name2IndexMap["beta0"]] * 
                       (1/state_param[name2IndexMap["N"]]) * 
                       state_param[name2IndexMap["Ca"]] + 
                       state_param[name2IndexMap["Ip"]] * 
                       state_param[name2IndexMap["beta0"]] * 
                       (1/state_param[name2IndexMap["N"]]) * 
                       state_param[name2IndexMap["Cp"]] ;
    
    std::cout << "result by direct calc = " << groudtruth << std::endl;    
}
