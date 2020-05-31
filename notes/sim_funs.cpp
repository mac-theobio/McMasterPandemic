// TMB version

enum param_pos {
		beta0 = 0,		
		mu = 13,
		N = 14,
		E0 = 15,
		phi1 = 19,
		zeta = 28
}
  
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_MATRIX(M);
  DATA_INTEGER(nt);
  DATA_SCALAR(dt);
  double foi;
  // make initial state

  for (int i=1; i<nt; i++) {

    
    // update_foi
    foi = beta0/params[pN]*(params[pCa]*params[pIa]+params[pCp]*params[pIp]+
			    (1-params[piso_m])*params[pCm]*params[pIm]+
			    (1-params[piso_s])*params[pCs]*params[pIs]);
    foi *= pow((S/N),zeta);

      M[rS,rE] = foi;
    
    // do_step
