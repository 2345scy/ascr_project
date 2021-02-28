#include <TMB.hpp>
template<class Type>
vector<Type> extractsubvector(vector<Type> v, int a, int b){
  int sublength = b - a + 1;
  vector<Type> subv(sublength);
  for(int m = 0; m < sublength; m++) subv[m] = v[a + m];
  return subv;
}


vector<int> incheck(vector<int> a, vector<int> b){
  int lengtha = a.size();
  int lengthb = b.size();
  vector<int> ans(lengtha);
  for(int i = 0; i < lengtha; i++){
    ans(i) = 0;
    for(int j = 0; j < lengthb; j++){
      if(a(i) == b(j)){
        ans(i) = 1;
        break;
      }
    }
  }
  return ans;
}

int incheck_scalar(int a, vector<int> b){
  int lengthb = b.size();
  int ans = 0;
  for(int i = 0; i < lengthb; i++){
    if(a == b(i)){
      ans = 1;
      break;
    }
  }
  return ans;
}

template<class Type>
vector<Type> isNA(vector<Type> x){
  int len = x.size();
  vector<Type> ans(len);
  for(int i = 0; i < len; i++){
    if(R_IsNA(asDouble(x(i)))){
      ans(i) = 1;
    } else {
      ans(i) = 0;
    }
  }
  return ans;
}




//here n_i is n_IDs_for_datafull
int lookup_data_full(int is_ani, int s, int a, int i, int t,
                     vector<int> n_a, vector<int> n_i,
                     vector<int> n_t, vector<int> n_i_each_a){
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


int lookup_data_dist_theta(int s, int t, int m,
                           vector<int> n_t, vector<int> n_m){
  int ans = 0;
  if(s > 1){
    for(int s_index = 1; s_index < s; s_index++){
      ans = ans + n_t(s_index - 1) * n_m(s_index - 1);
    }
  }
  ans = ans + (t - 1) * n_m(s - 1) + m - 1 ;
  
  return ans;
}


int lookup_data_mask(int s, int m, vector<int> n_m){
  int ans = 0;
  if(s > 1){
    for(int s_index = 1; s_index < s; s_index++){
      ans = ans + n_m(s_index - 1);
    }
  }
  ans = ans + m - 1;
  
  return ans;
}

//here n_i is n_IDs

int lookup_data_IDmask(int is_ani, int s, int a, int i, int m,
                       vector<int> n_a, vector<int> n_i, 
                       vector<int> n_i_each_a, vector<int> n_m){
  int ans = 0;
  if(is_ani == 1){
    ans = -1;
  } else {
    if(s > 1){
      for(int s_index = 1; s_index < s; s_index++){
        ans = ans + n_i(s_index - 1) * n_m(s_index - 1);
      }
    }
    
    ans = ans + (i - 1) * n_m(s - 1) + m - 1;
  }
  
  return ans;
}

//here n_i is n_IDs as well
int lookup_n_detection(int is_ani, int s, int a, int i, vector<int> n_a,
                       vector<int> n_i, vector<int> n_i_each_a){
  int ans = 0;
  if(is_ani == 1){
    ans = -1;
  } else {
    if(s > 1){
      for(int s_index = 1; s_index < s; s_index++){
        ans = ans + n_i(s_index - 1);
      }
    }
    
    ans = ans + i - 1;
  }
  
  return ans;
}


template<class Type>
Type det_hn(Type dx, Type g0, Type sigma){
  Type ans = 0.0;
  ans = g0 * exp(-1 * pow(dx, 2) / (2 * pow(sigma, 2)));
  return ans;
}

template<class Type>
Type trans(Type x, int link){
  Type ans = 0.0;
  if(link == 1) {
    ans = x;
  }
  if(link == 2){
    ans = exp(x);
  }
  if(link == 3){
    ans = exp(x)/(1 + exp(x));
  }
  
  return ans;
  
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(n_sessions);
  DATA_IVECTOR(n_animals);
  DATA_IVECTOR(n_IDs);
  DATA_IVECTOR(n_IDs_for_datafull); 
  DATA_IVECTOR(n_traps);
  DATA_IVECTOR(n_masks);
  DATA_IVECTOR(n_detection);
  DATA_IVECTOR(n_calls_each_animal);
  
  
  DATA_INTEGER(nrow_data_full);
  DATA_INTEGER(nrow_data_mask);
  DATA_INTEGER(nrow_dx);
  DATA_INTEGER(nrow_id_mask);
  
  
  DATA_VECTOR(A);
  DATA_VECTOR(survey_length);
  DATA_SCALAR(sound_speed);
  
  
  DATA_INTEGER(het_method);
  DATA_INTEGER(n_dir_quadpoints);
  DATA_INTEGER(n_het_source_quadpoints);
  DATA_VECTOR(het_source_nodes);
  DATA_VECTOR(het_source_weights);
  DATA_SCALAR(cutoff);
  
  
  DATA_INTEGER(detfn_index);
  DATA_IVECTOR(param_og);
  DATA_IMATRIX(par_n_col);
  DATA_IVECTOR(par_link);
  
  
  DATA_INTEGER(is_animalID);
  DATA_INTEGER(is_ss_origin);
  DATA_INTEGER(is_ss_het);
  DATA_INTEGER(is_ss_dir);
  DATA_INTEGER(is_bearing);
  DATA_INTEGER(is_toa);
  DATA_INTEGER(is_dist);
  DATA_INTEGER(is_local);
  DATA_INTEGER(is_freqs);
  
  
  DATA_VECTOR(capt_bin);
  DATA_VECTOR(capt_bearing);
  DATA_VECTOR(capt_dist);
  DATA_VECTOR(capt_ss);
  DATA_VECTOR(capt_toa);
  DATA_VECTOR(dx);
  DATA_VECTOR(theta);
  DATA_VECTOR(index_local);
  DATA_VECTOR(toa_ssq);
  
  
  DATA_MATRIX(g0_DX);
  DATA_MATRIX(sigma_DX);
  DATA_MATRIX(lambda0_DX);
  DATA_MATRIX(z_DX);
  DATA_MATRIX(shape_1_DX);
  DATA_MATRIX(shape_2_DX);
  DATA_MATRIX(shape_DX);
  DATA_MATRIX(scale_DX);
  DATA_MATRIX(b0_ss_DX);
  DATA_MATRIX(b1_ss_DX);
  DATA_MATRIX(b2_ss_DX);
  DATA_MATRIX(sigma_ss_DX);
  DATA_MATRIX(kappa_DX);
  DATA_MATRIX(alpha_DX);
  DATA_MATRIX(sigma_toa_DX);
  DATA_MATRIX(sigma_b0_ss_DX);
  DATA_MATRIX(D_DX);
  
  
  DATA_MATRIX(g0_DX_mask);
  DATA_MATRIX(sigma_DX_mask);
  DATA_MATRIX(lambda0_DX_mask);
  DATA_MATRIX(z_DX_mask);
  DATA_MATRIX(shape_1_DX_mask);
  DATA_MATRIX(shape_2_DX_mask);
  DATA_MATRIX(shape_DX_mask);
  DATA_MATRIX(scale_DX_mask);
  DATA_MATRIX(b0_ss_DX_mask);
  DATA_MATRIX(b1_ss_DX_mask);
  DATA_MATRIX(b2_ss_DX_mask);
  DATA_MATRIX(sigma_ss_DX_mask);
  DATA_MATRIX(kappa_DX_mask);
  DATA_MATRIX(alpha_DX_mask);
  DATA_MATRIX(sigma_toa_DX_mask);
  DATA_MATRIX(sigma_b0_ss_DX_mask);
  DATA_MATRIX(D_DX_mask);
  
  
  DATA_MATRIX(g0_bound);
  DATA_MATRIX(sigma_bound);
  DATA_MATRIX(lambda0_bound);
  DATA_MATRIX(z_bound);
  DATA_MATRIX(shape_1_bound);
  DATA_MATRIX(shape_2_bound);
  DATA_MATRIX(shape_bound);
  DATA_MATRIX(scale_bound);
  DATA_MATRIX(b0_ss_bound);
  DATA_MATRIX(b1_ss_bound);
  DATA_MATRIX(b2_ss_bound);
  DATA_MATRIX(sigma_ss_bound);
  DATA_MATRIX(kappa_bound);
  DATA_MATRIX(alpha_bound);
  DATA_MATRIX(sigma_toa_bound);
  DATA_MATRIX(sigma_b0_ss_bound);
  DATA_MATRIX(D_bound);
  
  
  //settle with "D" as it will be used regardless type of model
  PARAMETER_VECTOR(D);
  ADREPORT(D);
  
  vector<Type> D_full = D.head(par_n_col(16, 0));
  vector<Type> D_mask = D.tail(par_n_col(16, 1));
  vector<Type> D_vec_full = D_DX * D_full;
  vector<Type> D_vec_mask(nrow_data_mask);
  
  if(par_n_col(16, 1) == 0){
    D_vec_mask.setZero();
  } else {
    D_vec_mask = D_DX_mask * D_mask;
  }
  

  PARAMETER_VECTOR(kappa);
  vector<Type> kappa_full = kappa.head(par_n_col(12, 0));
  vector<Type> kappa_mask = kappa.tail(par_n_col(12, 1));
  vector<Type> kappa_vec_full = kappa_DX * kappa_full;
  vector<Type> kappa_vec_mask(nrow_data_mask);
  
  if(par_n_col(12, 1) == 0){
    kappa_vec_mask.setZero();
  } else {
    kappa_vec_mask = kappa_DX_mask * kappa_mask;
  }
  
  int check = incheck_scalar(int tem = 13, param_og);
  if(check == 1) ADREPORT(kappa);
  
  
  PARAMETER_VECTOR(alpha);
  vector<Type> alpha_full = alpha.head(par_n_col(13, 0));
  vector<Type> alpha_mask = alpha.tail(par_n_col(13, 1));
  vector<Type> alpha_vec_full = alpha_DX * alpha_full;
  vector<Type> alpha_vec_mask(nrow_data_mask);
  
  if(par_n_col(13, 1) == 0){
    alpha_vec_mask.setZero();
  } else {
    alpha_vec_mask = alpha_DX_mask * alpha_mask;
  }
  
  int check = incheck_scalar(int tem = 14, param_og);
  if(check == 1) ADREPORT(alpha);
  
  
  PARAMETER_VECTOR(sigma_toa);
  vector<Type> sigma_toa_full = sigma_toa.head(par_n_col(14, 0));
  vector<Type> sigma_toa_mask = sigma_toa.tail(par_n_col(14, 1));
  vector<Type> sigma_toa_vec_full = sigma_toa_DX * sigma_toa_full;
  vector<Type> sigma_toa_vec_mask(nrow_data_mask);

  if(par_n_col(14, 1) == 0){
    sigma_toa_vec_mask.setZero();
  } else {
    sigma_toa_vec_mask = sigma_toa_DX_mask * sigma_toa_mask;
  }

  int check = incheck_scalar(int tem = 15, param_og);
  if(check == 1) ADREPORT(sigma_toa);
  
  Type nll = Type(0.0);
  
  //model with detfn == 'hn'
  if(detfn_index == 1){
    PARAMETER_VECTOR(g0);
    ADREPORT(g0);
    
    vector<Type> g0_full = g0.head(par_n_col(0, 0));
    vector<Type> g0_mask = g0.tail(par_n_col(0, 1));
    vector<Type> g0_vec_full = g0_DX * g0_full;
    vector<Type> g0_vec_mask(nrow_data_mask);

    if(par_n_col(0, 1) == 0){
      g0_vec_mask.setZero();
    } else {
      vector<Type> g0_vec_mask = g0_DX_mask * g0_mask;
    }
    
    
    PARAMETER_VECTOR(sigma);
    ADREPORT(sigma);
    
    vector<Type> sigma_full = sigma.head(par_n_col(1, 0));
    vector<Type> sigma_mask = sigma.tail(par_n_col(1, 1));
    vector<Type> sigma_vec_full = sigma_DX * sigma_full;
    vector<Type> sigma_vec_mask(nrow_data_mask);

    if(par_n_col(1, 1) == 0){
      sigma_vec_mask.setZero();
    } else {
      vector<Type> sigma_vec_mask = sigma_DX_mask * sigma_mask;
    }
    
    //------------------------------------------------------------------

    for(int s = 1; s <= n_sessions; s++){
      int n_i = n_IDs(s - 1);
      int n_m = n_masks(s - 1);
      int n_t = n_traps(s - 1);
      Type area_unit = A(s - 1);
      
      
      //firstly calculate lambda(theta), the rate of Poisson distribution
      
      //currently we don't have ID - level extension, so p_dot is a vector
      //with length of n_mask and p_k is a matrix with n_mask * n_trap
      //in the future, if we add ID - level extension, both of them should
      //add one dimension
      
      //and incidentally calculate all of the p_k and p_dot
      Type lambda_theta = Type(0.0);
      //one p_dot for each mask point
      vector<Type> p_dot(n_m);
      //row index of p_k is for mask, colnum index is for trap/detector
      matrix<Type> p_k(n_m, n_t);

      for(int m = 1; m <= n_m; m++){
        int index_data_mask = lookup_data_mask(s, m, n_masks);
        //"D" is special, not related to trap nor ID, so we set id = 1 and t = 1
        int index_data_full_D = lookup_data_full(is_animalID, s, 0, 1, 1, 0, n_IDs_for_datafull, n_traps, 0);
        Type D_tem_full = D_vec_full[index_data_full_D];
        Type D_tem_mask = D_vec_mask[index_data_mask];
        Type D_tem = D_tem_full + D_tem_mask;
        D_tem = trans(D_tem, par_link(16));
        p_dot(m - 1) = Type(1.0);

        for(int t = 1; t <= n_t; t++){
          //the arguments in the function of "int a", "vector<int> n_a/n_i_each_a",
          //are set to 0, and "int i" is set to 1 because currently we do not have
          //"ID-level" parameter extension
          int index_data_full = lookup_data_full(is_animalID, s, 0, 1, t,
                                                 0, n_IDs_for_datafull, 
                                                 n_traps, 0);
          
          int index_data_dist_theta = lookup_data_dist_theta(s, t, m, n_traps, n_masks);
          
          
          Type g0_tem = g0_vec_full(index_data_full) + g0_vec_mask(index_data_mask);
          g0_tem = trans(g0_tem, par_link(0));
          
          Type sigma_tem = sigma_vec_full(index_data_full) + sigma_vec_mask(index_data_mask);
          sigma_tem = trans(sigma_tem, par_link(1));
          
          p_k(m - 1, t - 1) = det_hn(dx(index_data_dist_theta), g0_tem, sigma_tem);
          p_dot(m - 1) *= 1 - p_k(m - 1, t - 1);
          //end for trap t
        }
        
        p_dot(m - 1) = 1 - p_dot(m - 1);
        lambda_theta += D_tem * p_dot(m - 1);
        //end for mask m
      }
      
      //end of lambda_theta calculation
      
      
      //canceled out original likelihood: nll -= dpois(Type(n_i), lambda_theta, true);
      nll += lambda_theta;
      
      
      if(n_i != 0){
        for(int i = 1; i <= n_i; i++){
          //Z_i is the number of detections for ID i
          int Z_i = lookup_n_detection(is_animalID, s, 0, i, 0, n_IDs, 0);
          Type ll_i = Type(0.0);
          for(int m = 1; m <= n_m; m++){
            int index_data_mask = lookup_data_mask(s, m, n_masks);
            //"D" is special, not related to trap nor ID, so we set id = 1 and t = 1
            int index_data_full_D = lookup_data_full(is_animalID, s, 0, 1, 1,
                                                     0, n_IDs_for_datafull, n_traps, 0);
            Type D_tem = D_vec_full(index_data_full_D) + D_vec_mask(index_data_mask);
            D_tem = trans(D_tem, par_link(16));
            //since we sum up original likelihood for each mask, so the inital value
            //should be 1 instead of 0
            Type fx = Type(1.0);
            Type fw = Type(1.0);
            Type fy = Type(1.0);
            
            
            //for fx=f(x|n;theta)
            //canceled out p_dot & lambda from original likelihood: 
            //fx = D_tem * p_dot(m - 1) / lambda_theta;
            fx = D_tem;
            
            //for fw=f(w_i|x,n;theta)
            for(int t = 1; t <= n_t; t++){
              int index_data_full = lookup_data_full(is_animalID, s, 0, i, t,
                                                     0, n_IDs_for_datafull, n_traps, 0);
              fw *= pow(p_k(m - 1, t - 1), capt_bin(index_data_full)) * 
                pow((1 - p_k(m - 1, t - 1)), (1 - capt_bin(index_data_full)));
            }
            //cancelled out p_dot from original likelihood: fw /= p_dot(m - 1);

            //for fy=f(y_i|w_i,x,n;theta)
            
            
            
            
            
            //we sum up likelihood (not log-likelihood) of each mask 
            ll_i += fw * fx * fy;;
            //end for mask m
          }
          
          nll -= log(ll_i);
          std::cout << "nll until ID " << i << " : " << nll << std::endl;
          //end for ID i
        }
        
        //end for if(n_i != 0)
      }
      std::cout << "nll until session " << s << " : " << nll << std::endl;
      //end for session s
    }
    
    //end for model with detfn == "hn"
  }
  
  return nll;
}
