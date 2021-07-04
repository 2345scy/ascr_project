#include <TMB.hpp>

template<class Type>
Type Gamma(Type z){
  if(z < 0.5){
    z = 1 - z;
    Type y = sqrt(2 * 3.14159265358979323846) * pow(z + 6.5, (z - 0.5)) * exp((-1) * z - 6.5) * 
      (0.99999999999980993 + 676.5203681218851/(z) + (-1259.1392167224028)/(z + 1) + 771.32342877765313/(z + 2) +
      (-176.61502916214059)/(z + 3) + 12.507343278686905/(z + 4) + (-0.13857109526572012)/(z + 5));
    //+ 9.9843695780195716e-6/(z + 6) + 1.5056327351493116e-7/(z + 7));
    
    return 3.14159265358979323846 / (sin(3.14159265358979323846 * (z - 1)) * y);
  } else {
    Type y = sqrt(2 * 3.14159265358979323846) * pow(z + 6.5, (z - 0.5)) * exp((-1) * z - 6.5) * 
      (0.99999999999980993 + 676.5203681218851/(z) + (-1259.1392167224028)/(z + 1) + 771.32342877765313/(z + 2) +
      (-176.61502916214059)/(z + 3) + 12.507343278686905/(z + 4) + (-0.13857109526572012)/(z + 5)); 
    //+ 9.9843695780195716e-6/(z + 6) + 1.5056327351493116e-7/(z + 7));
    
    return y;
  }
}



template<class Type>
Type det_hn(const Type &dx, const Type &g0, const Type &sigma){
  Type ans = g0 * exp(-1 * pow(dx, 2) / (2 * pow(sigma, 2)));
  return ans;
}

//I could confirm this function is correct, because I tested it by comparing with
//the build-in besselI(). However, when I used this function for nll, the optimization
//does not converge.
template<class Type>
// Bessel function.
Type bessi0 (Type x){
    Type ax;
    Type ans;
    Type y;
    if ((ax = fabs(x)) < 3.75){
      y = x/3.75;
      y *= y;
      ans = 1.0 + y*(3.5156229 + y*(3.0899424 + y*(1.2067492 + y*(0.2659732 + y*(0.360768e-1 + y*0.45813e-2)))));
    } else {
      y = 3.75/ax;
      ans = (exp(ax)/sqrt(ax))*(0.39894228 + y*(0.1328592e-1 + y*(0.225319e-2 + 
        y*(-0.157565e-2 + y*(0.916281e-2 + y*(-0.2057706e-1 + y*(0.2635537e-1 + 
        y*(-0.1647633e-1 + y*0.392377e-2))))))));
    }
    return ans;
}

int lookup_data_full(const int &is_ani, const int &s, const int &a, const int &i, const int &t,
                     const vector<int> &n_a, const vector<int> &n_i,
                     const vector<int> &n_t, const vector<int> &n_i_each_a){
  int ans = 0;
  if(is_ani == 1){
    ans = -1;
  } else {
    if(s > 1){
      for(int s_index = 1; s_index < s; s_index++){
        ans = ans + n_i(s_index - 1) * n_t(s_index - 1);
      }
    }
    ans = ans + (i - 1) * n_t(s - 1) + t - 1;
  }
  
  return ans;
  
}

int lookup_data_dist_theta(const int &s, const int &t, const int &m,
                           const vector<int> &n_t, const vector<int> &n_m){
  int ans = 0;
  if(s > 1){
    for(int s_index = 1; s_index < s; s_index++){
      ans = ans + n_t(s_index - 1) * n_m(s_index - 1);
    }
  }
  ans = ans + (m - 1) * n_t(s - 1) + t - 1 ;
  
  return ans;
}


int lookup_bincapt_uid(const int &s, const int &uid, const int &t, const vector<int> &n_uid_session, 
	const vector<int> &n_t){
	
	int ans = 0;
	if(s > 1){
		for(int s_index = 1; s_index < s; s_index++){
			ans = ans + n_uid_session(s_index - 1) * n_t(s_index - 1);
		}
	}

	ans = ans + (uid - 1) * n_t(s - 1) + t - 1;

	return ans;
}

int lookup_n_det_and_id_uid(const int &s, const int &uid, const vector<int> &n_uid_session){
  	int ans = 0;

	if(s > 1){
		for(int s_index = 1; s_index < s; s_index++){
		ans = ans + n_uid_session(s_index - 1);
		}
	}

	ans = ans + uid - 1;
  
  	return ans;
}

vector<int> look_up_ids_uid(const int &s, const int &uid, const vector<int> &n_i, 
	const vector<int> n_ids_uid, const vector<int> &n_uid_session, const vector<int> &u_id_match){
	int index = 0;

	if(s > 1){
		for(int s_index = 1; s_index < s; s_index++){
			index += n_i(s_index - 1);
		}
	}

	int index_n_id_each_uid = lookup_n_det_and_id_uid(s, 1, n_uid_session);

	if(uid > 1){
		for(int uid_index = 1; uid_index < uid; uid_index++){
			index += n_ids_uid[index_n_id_each_uid];
			index_n_id_each_uid++;
		}
	}

	int n_ids = n_ids_uid[index_n_id_each_uid];

	vector<int> ans(n_ids);

	for(int i = 0; i < n_ids; i++){
		ans(i) = u_id_match(index);
		index++;
	}

	return ans;
}

int lookup_data_uid(const int &s, const int &uid, const int &t,
                     const vector<int> &n_uid, const vector<int> &n_t){
  int ans = 0;

	if(s > 1){
		for(int s_index = 1; s_index < s; s_index++){
		ans = ans + n_uid(s_index - 1) * n_t(s_index - 1);
		}
	}
	ans = ans + (uid - 1) * n_t(s - 1) + t - 1;
  
  
  return ans;
}

vector<int> look_up_traps_uid(const int &s, const int &uid, const int &index_n_det_uid, const int &n_det,
	const vector<int> n_det_uid, const vector<int> &index_traps_uid){
	
	vector<int> ans(n_det);
	int index = 0;
	
	if(index_n_det_uid > 0){
		for(int j = 0; j < index_n_det_uid; j++){
			index += n_det_uid[j];
		}
	}

	for(int k = 0; k < n_det; k++){
		ans(k) = index_traps_uid(index);
		index++;
	}
	
	return ans;
}

template<class Type>
Type objective_function<Type>::operator() (){
  Type nll = Type(0.0);
  //Type Inf = std::numeric_limits<double>::infinity();

  DATA_INTEGER(n_sessions);
  DATA_IVECTOR(n_IDs);
  DATA_IVECTOR(n_IDs_for_datafull); 
  DATA_IVECTOR(n_traps);
  DATA_IVECTOR(n_masks);
  DATA_IVECTOR(n_detection);
  DATA_VECTOR(A);


  DATA_VECTOR(capt_bin);
  DATA_VECTOR(capt_dist);
  DATA_VECTOR(dx);


  DATA_IVECTOR(n_uid_session);
  DATA_IVECTOR(n_ids_each_uid);
  DATA_IVECTOR(u_id_match);
  DATA_VECTOR(capt_bin_uid);
	DATA_IVECTOR(index_traps_uid);
  DATA_IVECTOR(n_detection_uid);

  PARAMETER(g0);
  PARAMETER(sigma);
  PARAMETER(alpha);
  PARAMETER(D);

  DATA_MATRIX(kappa_bound);

  Type g0_og = exp(g0)/(1 + exp(g0));
  Type sigma_og = exp(sigma);
  Type D_og = exp(D);
  Type alpha_og = exp(alpha);


	//if(kappa < kappa_bound(0, 0) || kappa > kappa_bound(1, 0)) nll += Inf;



  Type lambda_theta;
  Type l_i;
  Type D_tem;
  Type fy_log;
  //Type fy;
  Type bin_tem;
  int index_dx;
  int index_data_full;
  int n_i;
  int n_m;
  int n_t;
  Type area_unit;


  int n_uid;
  int id;
  int index_data_uid;
  int index_data_uid_inital;
  int index_n_det_and_id_uid;
  int n_det;

  for(int s = 1; s <= n_sessions; s++){
    n_m = n_masks(s - 1);
    n_t = n_traps(s - 1);
    area_unit = A(s - 1);

    n_uid = n_uid_session(s - 1);


    D_tem = D_og * area_unit;

    lambda_theta = Type(0.0);
    vector<Type> p_dot(n_m);
    matrix<Type> p_k(n_m, n_t);

    index_dx = lookup_data_dist_theta(s, 1, 1, n_traps, n_masks);


    for(int m = 1; m < n_m; m++){
      p_dot(m - 1) = Type(1.0);

      for(int t = 1; t <= n_t; t++){
        p_k(m - 1, t - 1) = det_hn(dx[index_dx], g0_og, sigma_og);
        index_dx++;
        p_dot(m - 1) *= 1 - p_k(m - 1, t - 1);
      }

      p_dot(m - 1) = 1 - p_dot(m - 1);
      lambda_theta += D_tem * p_dot(m - 1);
    }

    nll += lambda_theta;

    
    if(n_uid > 0){
      for(int uid = 1; uid <= n_uid; uid++){
        vector<Type> fw(n_m);
        index_data_uid_inital = lookup_data_uid(s, uid, 1, n_uid_session, n_traps);
        for(int m = 1; m <= n_m; m++){
          index_data_uid = index_data_uid_inital;
          fw(m - 1) = Type(1.0);
          for(int t = 1; t <= n_t; t++){
            bin_tem = capt_bin_uid(index_data_uid);
            fw(m - 1) *= pow(p_k(m - 1, t - 1), bin_tem) * 
              pow((1 - p_k(m - 1, t - 1)), (1 - bin_tem));
            
            index_data_uid++;
          }
        }

        
        index_n_det_and_id_uid = lookup_n_det_and_id_uid(s, uid, n_uid_session);
        n_i = n_ids_each_uid(index_n_det_and_id_uid);
        n_det = n_detection_uid(index_n_det_and_id_uid);
        
        vector<int> ids = look_up_ids_uid(s, uid, n_IDs, n_ids_each_uid, n_uid_session, u_id_match);
        vector<int> traps = look_up_traps_uid(s, uid, index_n_det_and_id_uid, n_det,
	                    n_detection_uid, index_traps_uid);

        for(int i = 0; i < n_i; i++){
          id = ids(i);
          l_i = Type(0.0);

          for(int m = 1; m <= n_m; m++){
            fy_log = Type(0.0);
            //fy = Type(1.0);

            for(int t = 0; t < n_det; t++){
              index_data_full = lookup_data_full(0, s, 0, id, traps(t),
                    0, n_IDs_for_datafull, n_traps, 0);
              index_dx = lookup_data_dist_theta(s, traps(t), m, n_traps, n_masks);

              //fy *=  pow(pow(dx(index_dx) / alpha_og, alpha_og) * Gamma(alpha_og), -1) *  
                //pow(capt_dist(index_data_full), (alpha_og - 1)) * 
                //exp((-1) * alpha_og * capt_dist(index_data_full) / dx(index_dx));
              
              fy_log +=  ((-1) * (alpha_og * (log(dx(index_dx)) - log(alpha_og)) + log(Gamma(alpha_og))) + (alpha_og - 1) * 
                log(capt_dist(index_data_full)) - alpha_og * capt_dist(index_data_full) / dx(index_dx));				
              
              
            }

            l_i += fw(m - 1) * exp(fy_log);
          }
          nll -= log(l_i * D_tem);
        }
        
      }
    }

  }

  return nll;
}
