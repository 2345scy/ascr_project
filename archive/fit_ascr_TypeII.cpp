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
  ans = ans + (m - 1) * n_t(s - 1) + t - 1 ;
  
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


int look_up_u(int s, int i, vector<int> n_i){
	int ans = 0;
	if(s > 1){
		for(int s_index = 1; s_index < s; s_index++){
			ans += n_i(s_index - 1);
		}
	}

	ans += i - 1;
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

//detect functions

//hn
template<class Type>
Type det_hn(Type dx, vector<Type> param){
  Type ans = 0.0;
  Type g0 = param(0);
  Type sigma = param(1);
  ans = g0 * exp(-1 * pow(dx, 2) / (2 * pow(sigma, 2)));
  return ans;
}

//hhn
template<class Type>
Type det_hhn(Type dx, vector<Type> param){
	Type ans = 0.0;
	Type lambda0 = param(0);
	Type sigma = param(1);
	Type lambda_d = lambda0 * exp(pow(dx, 2) * (-0.5) / pow(sigma, 2));
	ans = 1 - exp((-1) * lambda_d);
	return ans;
}

//hr
template<class Type>
Type det_hr(Type dx, vector<Type> param){
	Type ans = 0.0;
	Type g0 = param(0);
  	Type sigma = param(1);
	Type z = param(2);
	ans = g0 * (1 - exp((-1) * pow(dx / sigma, (-z))));
	return ans;
}

//lth
template<class Type>
Type det_lth(Type dx, vector<Type> param){
	Type ans = 0.0;
	Type shape_1 = param(0);
	Type shape_2 = param(1);
	Type scale = param(2);
	ans = 0.5 - 0.5 * erf(shape_1 - exp(shape_2 - scale * dx));
	return ans;
}

//th
template<class Type>
Type det_th(Type dx, vector<Type> param){
	Type ans = 0.0;
	Type shape = param(0);
	Type scale = param(1);
	ans = 0.5 - 0.5 * erf(dx / scale - shape);
	return ans;
}

//ss
//in ss, we do something different, here we use
//the same input, but output is back-transformed mu
//instead of the probability be detected
template<class Type>
Type mu_ss_identical(Type dx, vector<Type> param){
	Type b0_ss = param(0);
	Type b1_ss = param(1);
	Type mu = b0_ss - b1_ss * dx;
	return mu;
}

template<class Type>
Type mu_ss_log(Type dx, vector<Type> param){
	Type b0_ss = param(0);
	Type b1_ss = param(1);
	Type mu = exp(b0_ss - b1_ss * dx);
	return mu;
}

template<class Type>
Type mu_ss_spherical(Type dx, vector<Type> param){
	Type b0_ss = param(0);
	Type b1_ss = param(1);
	Type mu = 0.0;
	if(dx > 1){
		mu += b0_ss - 20 * log10(dx) - b1_ss * (dx - 1);
	} else {
		mu += b0_ss;
	}
	return mu;
}

//ss_dir
template<class Type>
Type mu_ss_dir_identical(Type dx, vector<Type> param){
	Type b0_ss = param(0);
	Type b1_ss = param(1);
	Type b2_ss = param(2);
	Type theta = param(3);
	Type mu = b0_ss - (b1_ss - b2_ss * (cos(theta) - 1)) * dx;
	return mu;
}

template<class Type>
Type mu_ss_dir_log(Type dx, vector<Type> param){
	Type b0_ss = param(0);
	Type b1_ss = param(1);
	Type b2_ss = param(2);
	Type theta = param(3);
	Type mu = exp(b0_ss - (b1_ss - b2_ss * (cos(theta) - 1)) * dx);
	return mu;
}

template<class Type>
Type mu_ss_dir_spherical(Type dx, vector<Type> param){
	Type b0_ss = param(0);
	Type b1_ss = param(1);
	Type b2_ss = param(2);
	Type theta = param(3);
	Type mu = 0.0;
	if(dx > 1){
		mu += b0_ss;
	} else {
		mu += b0_ss - 20 * log10(dx) - (b1_ss - b2_ss * (cos(theta) - 1)) * (dx - 1);
	}
	return mu;
}

//ss_het
template<class Type>
Type mu_ss_het(Type dx, vector<Type> param){
	Type b0_ss = param(0);
	Type b1_ss = param(1);
	Type u = param(2);
	Type mu = b0_ss - b1_ss * dx + u;
	return mu;
}

template<class Type>
Type objective_function<Type>::operator() ()
{
	Type nll = Type(0.0);
	Type Inf = std::numeric_limits<double>::infinity();
	
	
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
	DATA_INTEGER(ss_link);
  
  
	DATA_INTEGER(is_animalID);
	DATA_INTEGER(is_ss);
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
	//sigma_b0_ss is not extentable
	//DATA_MATRIX(sigma_b0_ss_DX);
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
	//sigma_b0_ss is not extentable
	//DATA_MATRIX(sigma_b0_ss_DX_mask);
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
  
	//define a pointer variable which point to a detect function
	//based on our input "detfn_index"
	//and also point out the number of parameters for that detfn

	Type (*detfn)(Type, vector<Type>);
	int n_detfn_param;
	if(detfn_index == 1){
		detfn = det_hn;
		n_detfn_param = 2;
	} else if(detfn_index == 2){
		detfn = det_hhn;
		n_detfn_param = 2;
	} else if(detfn_index == 3){
		detfn = det_hr;
		n_detfn_param = 3;
	} else if(detfn_index == 5){
		detfn = det_lth;
		n_detfn_param = 3;
	} else if(detfn_index == 4){
		detfn = det_th;
		n_detfn_param = 2;
	} else if(detfn_index == 6){
		if(ss_link == 1){
			detfn = mu_ss_identical;
		} else if(ss_link == 2){
			detfn = mu_ss_log;
		} else if(ss_link == 4){
			detfn = mu_ss_spherical;
		}
		n_detfn_param = 2;
	} else if(detfn_index == 7){
		if(ss_link == 1){
			detfn = mu_ss_dir_identical;
		} else if(ss_link == 2){
			detfn = mu_ss_dir_log;
		} else if(ss_link == 4){
			detfn = mu_ss_dir_spherical;
		}
		n_detfn_param = 4;
	} else if(detfn_index == 8){
		detfn = mu_ss_het;
		n_detfn_param = 3;
	}
	//define the parameter vector which will be used later
	vector<Type> detfn_param(n_detfn_param);

	//declare all parameters and generate its corresponding
	//data vectors based on data_full and data_mask
	//and only report this parameter if it is used in the model
	//to avoid NaN in std, outside of this .cpp, the unused
	//parameters should be fixed in MakeADFun() by "map=list(...)"


	//in addition, very important point, although all parameters
	//are used here in their original names, such as g0, however,
	//all of them are transformed by their corresponding link function
	//in the likelihood, they need to be backtransformed

	//g0
	PARAMETER_VECTOR(g0);

	for(int i = 0; i < (par_n_col(0, 0) + par_n_col(0, 1)); i++){
		if(g0(i) < g0_bound(0, i) || g0(i) > g0_bound(1, i)) nll += Inf; //infinity is Inf or INFINITY?
	}

	vector<Type> g0_full = g0.head(par_n_col(0, 0));
	vector<Type> g0_mask = g0.tail(par_n_col(0, 1));
	vector<Type> g0_vec_full = g0_DX * g0_full;
	vector<Type> g0_vec_mask(nrow_data_mask);

	if(par_n_col(0, 1) == 0){
		g0_vec_mask.setZero();
	} else {
		g0_vec_mask = g0_DX_mask * g0_mask;
	}
	if(incheck_scalar(1, param_og) == 1) ADREPORT(g0);

	//sigma
	PARAMETER_VECTOR(sigma);

	for(int i = 0; i < (par_n_col(1, 0) + par_n_col(1, 1)); i++){
		if(sigma(i) < sigma_bound(0, i) || sigma(i) > sigma_bound(1, i)) nll += Inf;
	}

	vector<Type> sigma_full = sigma.head(par_n_col(1, 0));
	vector<Type> sigma_mask = sigma.tail(par_n_col(1, 1));
	vector<Type> sigma_vec_full = sigma_DX * sigma_full;
	vector<Type> sigma_vec_mask(nrow_data_mask);

	if(par_n_col(1, 1) == 0){
		sigma_vec_mask.setZero();
	} else {
		sigma_vec_mask = sigma_DX_mask * sigma_mask;
	}
	if(incheck_scalar(2, param_og) == 1) ADREPORT(sigma);
	
	//lambda0
	PARAMETER_VECTOR(lambda0);

	for(int i = 0; i < (par_n_col(2, 0) + par_n_col(2, 1)); i++){
		if(lambda0(i) < lambda0_bound(0, i) || lambda0(i) > lambda0_bound(1, i)) nll += Inf;
	}

	vector<Type> lambda0_full = lambda0.head(par_n_col(2, 0));
	vector<Type> lambda0_mask = lambda0.tail(par_n_col(2, 1));
	vector<Type> lambda0_vec_full = lambda0_DX * lambda0_full;
	vector<Type> lambda0_vec_mask(nrow_data_mask);

	if(par_n_col(2, 1) == 0){
		lambda0_vec_mask.setZero();
	} else {
		lambda0_vec_mask = lambda0_DX_mask * lambda0_mask;
	}
	if(incheck_scalar(3, param_og) == 1) ADREPORT(lambda0);
  
	//z
	PARAMETER_VECTOR(z);

	for(int i = 0; i < (par_n_col(3, 0) + par_n_col(3, 1)); i++){
		if(z(i) < z_bound(0, i) || z(i) > z_bound(1, i)) nll += Inf;
	}


	vector<Type> z_full = z.head(par_n_col(3, 0));
	vector<Type> z_mask = z.tail(par_n_col(3, 1));
	vector<Type> z_vec_full = z_DX * z_full;
	vector<Type> z_vec_mask(nrow_data_mask);

	if(par_n_col(3, 1) == 0){
		z_vec_mask.setZero();
	} else {
		z_vec_mask = z_DX_mask * z_mask;
	}
	if(incheck_scalar(4, param_og) == 1) ADREPORT(z);

	//shape_1
	PARAMETER_VECTOR(shape_1);

	for(int i = 0; i < (par_n_col(4, 0) + par_n_col(4, 1)); i++){
		if(shape_1(i) < shape_1_bound(0, i) || shape_1(i) > shape_1_bound(1, i)) nll += Inf;
	}

	vector<Type> shape_1_full = shape_1.head(par_n_col(4, 0));
	vector<Type> shape_1_mask = shape_1.tail(par_n_col(4, 1));
	vector<Type> shape_1_vec_full = shape_1_DX * shape_1_full;
	vector<Type> shape_1_vec_mask(nrow_data_mask);

	if(par_n_col(4, 1) == 0){
		shape_1_vec_mask.setZero();
	} else {
		shape_1_vec_mask = shape_1_DX_mask * shape_1_mask;
	}
	if(incheck_scalar(5, param_og) == 1) ADREPORT(shape_1);

	//shape_2
	PARAMETER_VECTOR(shape_2);

	for(int i = 0; i < (par_n_col(5, 0) + par_n_col(5, 1)); i++){
		if(shape_2(i) < shape_2_bound(0, i) || shape_2(i) > shape_2_bound(1, i)) nll += Inf;
	}

	vector<Type> shape_2_full = shape_2.head(par_n_col(5, 0));
	vector<Type> shape_2_mask = shape_2.tail(par_n_col(5, 1));
	vector<Type> shape_2_vec_full = shape_2_DX * shape_2_full;
	vector<Type> shape_2_vec_mask(nrow_data_mask);

	if(par_n_col(5, 1) == 0){
		shape_2_vec_mask.setZero();
	} else {
		shape_2_vec_mask = shape_2_DX_mask * shape_2_mask;
	}
	if(incheck_scalar(6, param_og) == 1) ADREPORT(shape_2);

	//shape
	PARAMETER_VECTOR(shape);

	for(int i = 0; i < (par_n_col(6, 0) + par_n_col(6, 1)); i++){
		if(shape(i) < shape_bound(0, i) || shape(i) > shape_bound(1, i)) nll += Inf;
	}

	vector<Type> shape_full = shape.head(par_n_col(6, 0));
	vector<Type> shape_mask = shape.tail(par_n_col(6, 1));
	vector<Type> shape_vec_full = shape_DX * shape_full;
	vector<Type> shape_vec_mask(nrow_data_mask);

	if(par_n_col(6, 1) == 0){
		shape_vec_mask.setZero();
	} else {
		shape_vec_mask = shape_DX_mask * shape_mask;
	}
	if(incheck_scalar(7, param_og) == 1) ADREPORT(shape);

	//scale
	PARAMETER_VECTOR(scale);

	for(int i = 0; i < (par_n_col(7, 0) + par_n_col(7, 1)); i++){
		if(scale(i) < scale_bound(0, i) || scale(i) > scale_bound(1, i)) nll += Inf;
	}

	vector<Type> scale_full = scale.head(par_n_col(7, 0));
	vector<Type> scale_mask = scale.tail(par_n_col(7, 1));
	vector<Type> scale_vec_full = scale_DX * scale_full;
	vector<Type> scale_vec_mask(nrow_data_mask);

	if(par_n_col(7, 1) == 0){
		scale_vec_mask.setZero();
	} else {
		scale_vec_mask = scale_DX_mask * scale_mask;
	}
	if(incheck_scalar(8, param_og) == 1) ADREPORT(scale);

	//b0_ss
	PARAMETER_VECTOR(b0_ss);

	for(int i = 0; i < (par_n_col(8, 0) + par_n_col(8, 1)); i++){
		if(b0_ss(i) < b0_ss_bound(0, i) || b0_ss(i) > b0_ss_bound(1, i)) nll += Inf;
	}

	vector<Type> b0_ss_full = b0_ss.head(par_n_col(8, 0));
	vector<Type> b0_ss_mask = b0_ss.tail(par_n_col(8, 1));
	vector<Type> b0_ss_vec_full = b0_ss_DX * b0_ss_full;
	vector<Type> b0_ss_vec_mask(nrow_data_mask);

	if(par_n_col(8, 1) == 0){
		b0_ss_vec_mask.setZero();
	} else {
		b0_ss_vec_mask = b0_ss_DX_mask * b0_ss_mask;
	}
	if(incheck_scalar(9, param_og) == 1) ADREPORT(b0_ss);

	//b1_ss
	PARAMETER_VECTOR(b1_ss);

	for(int i = 0; i < (par_n_col(9, 0) + par_n_col(9, 1)); i++){
		if(b1_ss(i) < b1_ss_bound(0, i) || b1_ss(i) > b1_ss_bound(1, i)) nll += Inf;
	}

	vector<Type> b1_ss_full = b1_ss.head(par_n_col(9, 0));
	vector<Type> b1_ss_mask = b1_ss.tail(par_n_col(9, 1));
	vector<Type> b1_ss_vec_full = b1_ss_DX * b1_ss_full;
	vector<Type> b1_ss_vec_mask(nrow_data_mask);

	if(par_n_col(9, 1) == 0){
		b1_ss_vec_mask.setZero();
	} else {
		b1_ss_vec_mask = b1_ss_DX_mask * b1_ss_mask;
	}
	if(incheck_scalar(10, param_og) == 1) ADREPORT(b1_ss);

	//b2_ss
	PARAMETER_VECTOR(b2_ss);

	for(int i = 0; i < (par_n_col(10, 0) + par_n_col(10, 1)); i++){
		if(b2_ss(i) < b2_ss_bound(0, i) || b2_ss(i) > b2_ss_bound(1, i)) nll += Inf;
	}

	vector<Type> b2_ss_full = b2_ss.head(par_n_col(10, 0));
	vector<Type> b2_ss_mask = b2_ss.tail(par_n_col(10, 1));
	vector<Type> b2_ss_vec_full = b2_ss_DX * b2_ss_full;
	vector<Type> b2_ss_vec_mask(nrow_data_mask);

	if(par_n_col(10, 1) == 0){
		b2_ss_vec_mask.setZero();
	} else {
		b2_ss_vec_mask = b2_ss_DX_mask * b2_ss_mask;
	}
	if(incheck_scalar(11, param_og) == 1) ADREPORT(b2_ss);


	//sigma_ss
	PARAMETER_VECTOR(sigma_ss);

	for(int i = 0; i < (par_n_col(11, 0) + par_n_col(11, 1)); i++){
		if(sigma_ss(i) < sigma_ss_bound(0, i) || sigma_ss(i) > sigma_ss_bound(1, i)) nll += Inf;
	}

	vector<Type> sigma_ss_full = sigma_ss.head(par_n_col(11, 0));
	vector<Type> sigma_ss_mask = sigma_ss.tail(par_n_col(11, 1));
	vector<Type> sigma_ss_vec_full = sigma_ss_DX * sigma_ss_full;
	vector<Type> sigma_ss_vec_mask(nrow_data_mask);

	if(par_n_col(11, 1) == 0){
		sigma_ss_vec_mask.setZero();
	} else {
		sigma_ss_vec_mask = sigma_ss_DX_mask * sigma_ss_mask;
	}
	if(incheck_scalar(12, param_og) == 1) ADREPORT(sigma_ss);

	//kappa
	PARAMETER_VECTOR(kappa);

	for(int i = 0; i < (par_n_col(12, 0) + par_n_col(12, 1)); i++){
		if(kappa(i) < kappa_bound(0, i) || kappa(i) > kappa_bound(1, i)) nll += Inf;
	}


	vector<Type> kappa_full = kappa.head(par_n_col(12, 0));
	vector<Type> kappa_mask = kappa.tail(par_n_col(12, 1));
	vector<Type> kappa_vec_full = kappa_DX * kappa_full;
	vector<Type> kappa_vec_mask(nrow_data_mask);
  
	if(par_n_col(12, 1) == 0){
		kappa_vec_mask.setZero();
	} else {
		kappa_vec_mask = kappa_DX_mask * kappa_mask;
	}
	if(incheck_scalar(13, param_og) == 1) ADREPORT(kappa);

	//alpha
	PARAMETER_VECTOR(alpha);

	for(int i = 0; i < (par_n_col(13, 0) + par_n_col(13, 1)); i++){
		if(alpha(i) < alpha_bound(0, i) || alpha(i) > alpha_bound(1, i)) nll += Inf;
	}


	vector<Type> alpha_full = alpha.head(par_n_col(13, 0));
	vector<Type> alpha_mask = alpha.tail(par_n_col(13, 1));
	vector<Type> alpha_vec_full = alpha_DX * alpha_full;
	vector<Type> alpha_vec_mask(nrow_data_mask);
  
	if(par_n_col(13, 1) == 0){
		alpha_vec_mask.setZero();
	} else {
		alpha_vec_mask = alpha_DX_mask * alpha_mask;
	}
	if(incheck_scalar(14, param_og) == 1) ADREPORT(alpha);
  
	//sigma_toa
	PARAMETER_VECTOR(sigma_toa);

	for(int i = 0; i < (par_n_col(14, 0) + par_n_col(14, 1)); i++){
		if(sigma_toa(i) < sigma_toa_bound(0, i) || sigma_toa(i) > sigma_toa_bound(1, i)) nll += Inf;
	}

	vector<Type> sigma_toa_full = sigma_toa.head(par_n_col(14, 0));
	vector<Type> sigma_toa_mask = sigma_toa.tail(par_n_col(14, 1));
	vector<Type> sigma_toa_vec_full = sigma_toa_DX * sigma_toa_full;
	vector<Type> sigma_toa_vec_mask(nrow_data_mask);
  
	if(par_n_col(14, 1) == 0){
		sigma_toa_vec_mask.setZero();
	} else {
		sigma_toa_vec_mask = sigma_toa_DX_mask * sigma_toa_mask;
	}
	if(incheck_scalar(15, param_og) == 1) ADREPORT(sigma_toa);

	//sigma_b0_ss, this is not extentable, so just declare it as a scalar
	PARAMETER(sigma_b0_ss);

	if(sigma_b0_ss < sigma_b0_ss_bound(0, 0) || sigma_b0_ss > sigma_b0_ss_bound(1, 0)) nll += Inf;


	if(incheck_scalar(16, param_og) == 1) ADREPORT(sigma_b0_ss);
  
	//settle with "D" as it will be used regardless type of model
	PARAMETER_VECTOR(D);

	for(int i = 0; i < (par_n_col(16, 0) + par_n_col(16, 1)); i++){
		if(D(i) < D_bound(0, i) || D(i) > D_bound(1, i)) nll += Inf;
	}

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
  
	//finally declare the latent variable "u"
	PARAMETER_VECTOR(u);
	
	//declare variables for all backtransformed parameters 
	Type g0_tem;
	Type sigma_tem;
	Type lambda0_tem;
	Type z_tem;
	Type shape_1_tem;
	Type shape_2_tem;
	Type shape_tem;
	Type scale_tem;
	Type b0_ss_tem;
	Type b1_ss_tem;
	Type b2_ss_tem;
	Type sigma_ss_tem;
	Type kappa_tem;
	Type alpha_tem;
	Type sigma_toa_tem;
	Type sigma_b0_ss_tem;
	Type D_tem;

	//mu = E(ss|x), since it is "session-mask-trap"
	//level data, use data_dist_theta's index
	vector<Type> mu(nrow_dx);
	mu.setZero();

///////////////////////////////////
    

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
			int index_data_full_D = lookup_data_full(is_animalID, s, 0, 1, 1, 0, 
				n_IDs_for_datafull, n_traps, 0);
			D_tem = D_vec_full[index_data_full_D] + D_vec_mask[index_data_mask];
			D_tem = area_unit * trans(D_tem, par_link(16));

			//if(D_tem < D_bound(0, 0) || D_tem > D_bound(1, 0)) nll += Inf;			

			p_dot(m - 1) = Type(1.0);

			for(int t = 1; t <= n_t; t++){
				//the arguments in the function of "int a", "vector<int> n_a/n_i_each_a",
				//are set to 0, and "int i" is set to 1 because currently we do not have
				//"ID-level" parameter extension
				int index_data_full = lookup_data_full(is_animalID, s, 0, 1, t,
					0, n_IDs_for_datafull, n_traps, 0);

				int index_data_dist_theta = lookup_data_dist_theta(s, t, m, n_traps, n_masks);
				
				if(detfn_index == 1){
					g0_tem = g0_vec_full(index_data_full) + g0_vec_mask(index_data_mask);
					g0_tem = trans(g0_tem, par_link(0));
					sigma_tem = sigma_vec_full(index_data_full) + sigma_vec_mask(index_data_mask);
					sigma_tem = trans(sigma_tem, par_link(1));

					//if(g0_tem < g0_bound(0, 0) || g0_tem > g0_bound(1, 0)) nll += Inf;
					//if(sigma_tem < sigma_bound(0, 0) || sigma_tem > sigma_bound(1, 0)) nll += Inf;

					detfn_param(0) = g0_tem;
					detfn_param(1) = sigma_tem;

				} else if(detfn_index == 2){
					lambda0_tem = lambda0_vec_full(index_data_full) + lambda0_vec_mask(index_data_mask);
					lambda0_tem = trans(lambda0_tem, par_link(2));
					sigma_tem = sigma_vec_full(index_data_full) + sigma_vec_mask(index_data_mask);
					sigma_tem = trans(sigma_tem, par_link(1));
					
					//if(lambda0_tem < lambda0_bound(0, 0) || lambda0_tem > lambda0_bound(1, 0)) nll += Inf;
					//if(sigma_tem < sigma_bound(0, 0) || sigma_tem > sigma_bound(1, 0)) nll += Inf;

					detfn_param(0) = lambda0_tem;
					detfn_param(1) = sigma_tem;

				} else if(detfn_index == 3){
					g0_tem = g0_vec_full(index_data_full) + g0_vec_mask(index_data_mask);
					g0_tem = trans(g0_tem, par_link(0));
					sigma_tem = sigma_vec_full(index_data_full) + sigma_vec_mask(index_data_mask);
					sigma_tem = trans(sigma_tem, par_link(1));
					z_tem = z_vec_full(index_data_full) + z_vec_mask(index_data_mask);
					z_tem = trans(z_tem, par_link(3));

					//if(g0_tem < g0_bound(0, 0) || g0_tem > g0_bound(1, 0)) nll += Inf;
					//if(sigma_tem < sigma_bound(0, 0) || sigma_tem > sigma_bound(1, 0)) nll += Inf;
					//if(z_tem < z_bound(0, 0) || z_tem > z_bound(1, 0)) nll += Inf;

					detfn_param(0) = g0_tem;
					detfn_param(1) = sigma_tem;	
					detfn_param(2) = z_tem;	

				} else if(detfn_index == 5){
					shape_1_tem = shape_1_vec_full(index_data_full) + shape_1_vec_mask(index_data_mask);
					shape_1_tem = trans(shape_1_tem, par_link(4));
					shape_2_tem = shape_2_vec_full(index_data_full) + shape_2_vec_mask(index_data_mask);
					shape_2_tem = trans(shape_2_tem, par_link(5));
					scale_tem = scale_vec_full(index_data_full) + scale_vec_mask(index_data_mask);
					scale_tem = trans(scale_tem, par_link(7));

					//if(shape_1_tem < shape_1_bound(0, 0) || shape_1_tem > shape_1_bound(1, 0)) nll += Inf;
					//if(shape_2_tem < shape_2_bound(0, 0) || shape_2_tem > shape_2_bound(1, 0)) nll += Inf;
					//if(scale_tem < scale_bound(0, 0) || scale_tem > scale_bound(1, 0)) nll += Inf;

					detfn_param(0) = shape_1_tem;
					detfn_param(1) = shape_2_tem;	
					detfn_param(2) = scale_tem;
					
				} else if(detfn_index == 4){
					shape_tem = shape_vec_full(index_data_full) + shape_vec_mask(index_data_mask);
					shape_tem = trans(shape_tem, par_link(6));
					scale_tem = scale_vec_full(index_data_full) + scale_vec_mask(index_data_mask);
					scale_tem = trans(scale_tem, par_link(7));

					//if(shape_tem < shape_bound(0, 0) || shape_tem > shape_bound(1, 0)) nll += Inf;
					//if(scale_tem < scale_bound(0, 0) || scale_tem > scale_bound(1, 0)) nll += Inf;

					detfn_param(0) = shape_tem;
					detfn_param(1) = scale_tem;	
				} else if(detfn_index == 6){
					b0_ss_tem = b0_ss_vec_full(index_data_full) + b0_ss_vec_mask(index_data_mask);
					b0_ss_tem = trans(b0_ss_tem, par_link(8));
					b1_ss_tem = b1_ss_vec_full(index_data_full) + b1_ss_vec_mask(index_data_mask);
					b1_ss_tem = trans(b1_ss_tem, par_link(9));

					//if(b0_ss_tem < b0_ss_bound(0, 0) || b0_ss_tem > b0_ss_bound(1, 0)) nll += Inf;
					//if(b1_ss_tem < b1_ss_bound(0, 0) || b1_ss_tem > b1_ss_bound(1, 0)) nll += Inf;

					detfn_param(0) = b0_ss_tem;
					detfn_param(1) = b1_ss_tem;

				} else if(detfn_index == 7){
					b0_ss_tem = b0_ss_vec_full(index_data_full) + b0_ss_vec_mask(index_data_mask);
					b0_ss_tem = trans(b0_ss_tem, par_link(8));
					b1_ss_tem = b1_ss_vec_full(index_data_full) + b1_ss_vec_mask(index_data_mask);
					b1_ss_tem = trans(b1_ss_tem, par_link(9));
					b2_ss_tem = b2_ss_vec_full(index_data_full) + b2_ss_vec_mask(index_data_mask);
					b2_ss_tem = trans(b2_ss_tem, par_link(10));

					//if(b0_ss_tem < b0_ss_bound(0, 0) || b0_ss_tem > b0_ss_bound(1, 0)) nll += Inf;
					//if(b1_ss_tem < b1_ss_bound(0, 0) || b1_ss_tem > b1_ss_bound(1, 0)) nll += Inf;
					//if(b2_ss_tem < b2_ss_bound(0, 0) || b2_ss_tem > b2_ss_bound(1, 0)) nll += Inf;

					detfn_param(0) = b0_ss_tem;
					detfn_param(1) = b1_ss_tem;
					detfn_param(2) = b2_ss_tem;
					detfn_param(3) = theta(index_data_dist_theta);

				} 
				//het is detfn_index == 8, not sure how to do it yet



				if(is_ss == 0){
					p_k(m - 1, t - 1) = (*detfn)(dx(index_data_dist_theta), detfn_param);
				} else if (is_ss_origin == 1){
					sigma_ss_tem = sigma_ss_vec_full(index_data_full) + sigma_ss_vec_mask(index_data_mask);
					sigma_ss_tem = trans(sigma_ss_tem, par_link(11));

					//if(sigma_ss_tem < sigma_ss_bound(0, 0) || sigma_ss_tem > sigma_ss_bound(1, 0)) nll += Inf;

					mu(index_data_dist_theta) = (*detfn)(dx(index_data_dist_theta), detfn_param);
					p_k(m - 1, t - 1) = 1 - pnorm((cutoff - mu(index_data_dist_theta)) / sigma_ss_tem);
				}
				std::cout <<  "dx for s=" << s << ", m="<< m << ", t="<< t <<": " << dx(index_data_dist_theta) << std::endl;
				std::cout <<  "p_k for s=" << s << ", m="<< m << ", t="<< t <<": " << p_k(m - 1, t - 1) << std::endl;
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
				int index_zi = lookup_n_detection(is_animalID, s, 0, i, 0, n_IDs, 0);
				Type Z_i = n_detection(index_zi);
				Type l_i = Type(0.0);
				for(int m = 1; m <= n_m; m++){
					int index_data_IDmask = lookup_data_IDmask(is_animalID, s, 0, i, m,
						0, n_IDs, 0, n_masks);
					int index_data_mask = lookup_data_mask(s, m, n_masks);
					//"D" is special, not related to trap nor ID, so we set id = 1 and t = 1
					//"sigma_toa" is not trap extendable nor ID either, so take this index as well
					int index_data_full_D = lookup_data_full(is_animalID, s, 0, 1, 1,
						0, n_IDs_for_datafull, n_traps, 0);
					Type D_tem = D_vec_full(index_data_full_D) + D_vec_mask(index_data_mask);
					D_tem = area_unit * trans(D_tem, par_link(16));
					//since we sum up original likelihood for each mask, so the inital value
					//should be 1 instead of 0
					Type fx = Type(1.0);
					Type fw = Type(1.0);
					
					
					
					//for fx=f(x|n;theta)
					//canceled out p_dot & lambda from original likelihood: 
					//fx = D_tem * p_dot(m - 1) / lambda_theta;
					fx = D_tem;
					
					//for fw=f(w_i|x,n;theta)
					for(int t = 1; t <= n_t; t++){
						int index_data_full = lookup_data_full(is_animalID, s, 0, i, t,
							0, n_IDs_for_datafull, n_traps, 0);
						
						if(is_ss == 0){
							fw *= pow(p_k(m - 1, t - 1), capt_bin(index_data_full)) * 
								pow((1 - p_k(m - 1, t - 1)), (1 - capt_bin(index_data_full)));
						} else if(is_ss_het == 0){
							//pow(p_k(), capt_bin()) term could be cancelled out with fy_ss
							fw *= pow((1 - p_k(m - 1, t - 1)), (1 - capt_bin(index_data_full)));
						}

					}
					//cancelled out p_dot from original likelihood: fw /= p_dot(m - 1);

					//the section below for 'fy' will be the same for all kinds of detfn excluding ss
					//make this part a function later
					//for fy=f(y_i|w_i,x,n;theta)
					Type fy = Type(1.0);
					
					Type fy_toa = Type(1.0);
					Type fy_bear = Type(1.0);
					Type fy_dist = Type(1.0);
					Type fy_ss = Type(1.0);

					//toa
					if(is_toa == 1){
						Type sigma_toa_tem = sigma_toa_vec_full(index_data_full_D) + 
							sigma_toa_vec_mask(index_data_mask);
						sigma_toa_tem = trans(sigma_toa_tem, par_link(14));

						//if(sigma_toa_tem < sigma_toa_bound(0, 0) || sigma_toa_tem > sigma_toa_bound(1, 0)) nll += Inf;

						Type toa_ssq_tem = toa_ssq(index_data_IDmask);
						fy_toa *= pow(sigma_toa_tem, (1 - Z_i)) * exp((-0.5) * toa_ssq_tem / pow(sigma_toa_tem, 2));
					}
					
					//bearing
					if(is_bearing == 1){
						for(int t = 1; t <= n_t; t++){
							int index_data_full = lookup_data_full(is_animalID, s, 0, i, t,
								0, n_IDs_for_datafull, n_traps, 0);
							int index_data_dist_theta = lookup_data_dist_theta(s, t, m, n_traps, n_masks);
							Type kappa_tem = kappa_vec_full(index_data_full) + kappa_vec_mask(index_data_mask);
							kappa_tem = trans(kappa_tem, par_link(12));

							//if(kappa_tem < kappa_bound(0, 0) || kappa_tem > kappa_bound(1, 0)) nll += Inf;

							fy_bear *= pow(exp(kappa_tem * cos(capt_bearing(index_data_full) - theta(index_data_dist_theta))) / 
								besselI(kappa_tem, Type(0)), capt_bin(index_data_full));
						}
					}

					
					//dist
					if(is_dist == 1){
						for(int t = 1; t <= n_t; t++){
							int index_data_full = lookup_data_full(is_animalID, s, 0, i, t,
																0, n_IDs_for_datafull, n_traps, 0);
							int index_data_dist_theta = lookup_data_dist_theta(s, t, m, n_traps, n_masks);
							Type alpha_tem = alpha_vec_full(index_data_full) + alpha_vec_mask(index_data_mask);
							alpha_tem = trans(alpha_tem, par_link(13));
							
							//if(alpha_tem < alpha_bound(0, 0) || alpha_tem > alpha_bound(1, 0)) nll += Inf;
							
							fy_dist *= pow(pow(pow(dx(index_data_dist_theta)/alpha_tem, alpha_tem) * exp(lgamma(alpha_tem)), (-1)) * pow(capt_dist(index_data_full), 
								(alpha_tem - 1)) * exp((-1) * alpha_tem * capt_dist(index_data_full) / dx(index_data_dist_theta)), capt_bin(index_data_full));
						}
					}

					//ss
					if(is_ss_origin == 1){
						for(int t = 1; t <= n_t; t++){
							int index_data_full = lookup_data_full(is_animalID, s, 0, i, t,
																0, n_IDs_for_datafull, n_traps, 0);
							int index_data_dist_theta = lookup_data_dist_theta(s, t, m, n_traps, n_masks);
							sigma_ss_tem = sigma_ss_vec_full(index_data_full) + sigma_ss_vec_mask(index_data_mask);
							sigma_ss_tem = trans(sigma_ss_tem, par_link(11));

							//if(sigma_ss_tem < sigma_ss_bound(0, 0) || sigma_ss_tem > sigma_ss_bound(1, 0)) nll += Inf;

							Type mu_tem = mu(index_data_dist_theta);
							
							fy_ss *= pow(dnorm(capt_ss(index_data_full), mu_tem, sigma_ss_tem), capt_bin(index_data_full));
						}
					}
					
					
					fy *= fy_toa * fy_bear * fy_dist * fy_ss;
					//end of the section for 'fy'

					//we sum up likelihood (not log-likelihood) of each mask 
					l_i += fw * fx * fy;
					//end for mask m
				}
			  
			  nll -= log(l_i);
			  //end for ID i
			}

			//end for if(n_i != 0)
		}
      //end for session s
    }

  
  return nll;
}

