/// @file distmodel.hpp
#ifndef distmodel_hpp
#define distmodel_hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type distmodel(objective_function<Type>* obj) {
    DATA_VECTOR(y); // data vector
    DATA_UPDATE(y); // data vector
    DATA_INTEGER(dclass);
    PARAMETER(mu);
    PARAMETER(sigma);
    PARAMETER(skew);
    PARAMETER(shape);
    PARAMETER(lambda);
    Type nll = 0.0;
    vector<Type> z = (y.array() - mu) * (Type(1.0)/sigma);
    REPORT(z);
    vector<Type>llh_vec = distfun::distlike(z, skew, shape, lambda, dclass);
    // scale to avoid problems in estimation
    llh_vec = llh_vec.array()/sigma;
    REPORT(llh_vec);
    vector<Type> tmp = llh_vec.array().log();
    REPORT(tmp);
    nll += Type(-1.0) * tmp.array().sum();
    return(nll);
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
