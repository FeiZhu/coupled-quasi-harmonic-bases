#include <cmath>
#include <ctime>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "nlopt.hpp"
#include "OptQNIPS.h"
#include "NonLinearEquation.h"
#include "CompoundConstraint.h"
#include "CoupledQuasiHarmonics.h"
using std::fstream;
using OPTPP::NLP;
using OPTPP::NLF0;
using OPTPP::NLF1;
using OPTPP::FDNLF1;
using OPTPP::OptQNIPS;
using OPTPP::NonLinearEquation;
using OPTPP::Constraint;
using OPTPP::CompoundConstraint;
using OPTPP::NLPFunction;
using OPTPP::NLPGradient;
using OPTPP::LineSearch;
using OPTPP::TrustRegion;
using OPTPP::TrustPDS;
using OPTPP::VanShanno;

//implemenation
int CoupledQuasiHarmonics::static_col=0;
int CoupledQuasiHarmonics::static_objective_type=0;
vector<double> CoupledQuasiHarmonics::static_x;

CoupledQuasiHarmonics::CoupledQuasiHarmonics(double **input_base_a, double **input_base_b, 
                                             double *input_eigenvalue_a, double *input_eigenvalue_b,
                                             int input_row_a, int input_row_b, int input_col,
                                             double *input_vertex_area_a, double *input_vertex_area_b, 
                                             double **input_F, double **input_G, int input_p)
    :base_a(input_base_a),base_b(input_base_b),eigenvalue_a(input_eigenvalue_a),eigenvalue_b(input_eigenvalue_b),
     row_a(input_row_a),row_b(input_row_b),col(input_col),vertex_area_a(input_vertex_area_a),vertex_area_b(input_vertex_area_b),
     F(input_F),G(input_G),p(input_p),temp_x(NULL),temp_f(0),temp_f_grad(NULL),F_transpose_base_a_A_minus_G_transpose_base_b_B(NULL),
     mu(0.132),f_eval_num(0),package_type(OPTPP),offterm_type(ORDERED),inner_product_type(AREA_WEIGHTED),objective_type(A_AND_B),
     scale_base(true),scale_corresponding_function(true),enable_output(true)
{
    static_col=col;
    static_objective_type=objective_type;
    //eigenfunctions have the same coupling strength by default
    coupling_scale.resize(col);
    for(int i=0;i<col;++i)
	coupling_scale[i]=1.0;
    //compute the mean value of base_a and base_b, for scaling
    vector<double> data;
    for(int i=0;i<row_a;++i)
	for(int j=0;j<col;++j)
	    data.push_back(fabs(base_a[i][j]));
    sort(data.begin(),data.end());
    base_a_mean=data[data.size()/2];
    data.clear();
    for(int i=0;i<row_b;++i)
	for(int j=0;j<col;++j)
	    data.push_back(fabs(base_b[i][j]));
    sort(data.begin(),data.end());
    base_b_mean=data[data.size()/2];
    //std::cout<<base_a_mean<<" "<<base_b_mean<<"\n";
}

CoupledQuasiHarmonics::~CoupledQuasiHarmonics()
{
    if(temp_x)
	delete[] temp_x;
    if(temp_f_grad)
	delete[] temp_f_grad;
    if(F_transpose_base_a_A_minus_G_transpose_base_b_B)
	delete[] F_transpose_base_a_A_minus_G_transpose_base_b_B;
}

int CoupledQuasiHarmonics::getCoupledBases(int new_col, double **new_base_a, double **new_base_b)
{
    //scale the bases before optimization
    if(scale_base)
    {
	double scale=base_a_mean/base_b_mean;
	for(int i=0;i<row_b;++i)
	{
	    for(int j=0;j<col;++j)
		base_b[i][j]*=scale;
	    vertex_area_b[i]/=scale*scale;
	}
    }
    //scale F and G with area of the region
    if(scale_corresponding_function)
    {
	for(int j=0;j<p;++j)
	{
	    double sum=0;
	    for(int i=0;i<row_a;++i)
		if(F[i][j]!=0)
		    sum+=vertex_area_a[i];
	    sum/=3.0;
	    for(int i=0;i<row_a;++i)
		F[i][j]/=sum;
	    sum=0;
	    for(int i=0;i<row_b;++i)
		if(G[i][j]!=0)
		    sum+=vertex_area_b[i];
	    sum/=3.0;
	    for(int i=0;i<row_b;++i)
		G[i][j]/=sum;
	}
    }
    //transformation matrix A & B are both of size (col x new_col)
    int num_unknown=2*new_col*col;
    int constraint_num=new_col*(new_col+1);
    if(objective_type==ONLY_A)
    {
	num_unknown=new_col*col;//only A
	constraint_num=new_col*(new_col+1)/2;
    }
    //expand A & B in col-order into x, A first, then B
    static_x.resize(num_unknown);
    getInitialGuess(static_x);//get initial guess
    temp_x=new double[num_unknown];
    F_transpose_base_a_A_minus_G_transpose_base_b_B=new double[p*new_col];//p x new_col size matrix, expanded in row order,used in objectiveGrad()
    if(package_type==NLOPT)
    {
	//nlopt::opt solver(nlopt::LN_COBYLA,num_unknown);//derivative-free method
	nlopt::opt solver(nlopt::AUGLAG,num_unknown);//augmented lagrangian
	solver.set_min_objective(nloptObjective,this);
	vector<double> tol(constraint_num,1.0e-6);
	solver.add_equality_mconstraint(nloptConstraint,this,tol);
	solver.set_xtol_rel(1e-6);
	solver.set_ftol_rel(1e-6);
	//solver.set_xtol_abs(1e-6);
	//solver.set_ftol_abs(1e-6);
	solver.set_maxeval(5000);
	nlopt::opt local_solver(nlopt::LD_MMA,num_unknown);
	local_solver.set_xtol_rel(1e-6);
	local_solver.set_ftol_rel(1e-6);
	//local_solver.set_xtol_abs(1e-6);
	//local_solver.set_ftol_abs(1e-6);
	local_solver.set_maxeval(5000);
	solver.set_local_optimizer(local_solver);
	double minf;
	nlopt::result result;
	result=solver.optimize(static_x,minf);
	if(result<0)
	    return -1;
	//work around for the NLopt bug
	if(result==nlopt::MAXEVAL_REACHED)
	{
	    for(int i=0;i<num_unknown;++i)
		static_x[i]=temp_x[i];
	    minf=temp_f;
	}
    }
    else if(package_type==OPTPP)
    {
	temp_f_grad=new double[num_unknown];
        //define equality constraint problem
        NLP *csp=new NLP(new NLF1(num_unknown,constraint_num,optppConstraint,initOptpp));
	ColumnVector c_rhs(constraint_num);
	for(int i=1;i<=c_rhs.size();++i)
	    c_rhs(i)=0;
	Constraint eqn=new NonLinearEquation(csp,c_rhs,constraint_num);
	CompoundConstraint *constraint=new CompoundConstraint(eqn); 
        //define the optimization problem
	NLF1 nlp(num_unknown,optppObjective,initOptpp,constraint,this);
        //define the solver algorithm
	OptQNIPS solver(&nlp);
	solver.setMeritFcn(VanShanno);
	solver.setSearchStrategy(LineSearch);
	solver.setMaxBacktrackIter(1.0e6);
	solver.setFcnTol(1e-8);
	solver.setConTol(1e-8);
	solver.setStepTol(1e-8);
	solver.setMaxIter(15000);
	solver.setMaxFeval(15000);
	solver.optimize();
	ColumnVector final_x=nlp.getXc();
	for(int i=0;i<num_unknown;++i)
	    static_x[i]=final_x(i+1);
	//solver.setOutputFile(std::cout);
	solver.printStatus("OPTPP: ");//The status is wirtten into .out file by opt++
	solver.cleanup();
    }
    else if(package_type==MATLAB)
    {
	//saveDataForMatlab("MatlabSolverInputs.m",static_x,row_a,row_b,col,eigenvalue_a,eigenvalue_b,mu,F,G,p,base_a,base_b,vertex_area_a,vertex_area_b);
	loadResultFromMatlab("result.txt",static_x);
    }
    else
    {
	std::cout<<"Error: No optimization package specified.\n";
	return -1;
    }
    //scale back F and G with area of the region
    if(scale_corresponding_function)
    {
	for(int j=0;j<p;++j)
	{
	    double sum=0;
	    for(int i=0;i<row_a;++i)
		if(F[i][j]!=0)
		    sum+=vertex_area_a[i];
	    sum/=3.0;
	    for(int i=0;i<row_a;++i)
		F[i][j]*=sum;
	    sum=0;
	    for(int i=0;i<row_b;++i)
		if(G[i][j]!=0)
		    sum+=vertex_area_b[i];
	    sum/=3.0;
	    for(int i=0;i<row_b;++i)
		G[i][j]*=sum;
	}
    }
    //scale back the bases after the optimization
    if(scale_base)
    {
	double scale=base_a_mean/base_b_mean;
	for(int i=0;i<row_b;++i)
	{
	    for(int j=0;j<col;++j)
		base_b[i][j]/=scale;
	    vertex_area_b[i]*=scale*scale;
	}
    }
    //compute the new bases
    for(int i=0;i<row_a;++i)
	for(int j=0;j<new_col;++j)
	{
	    new_base_a[i][j]=0.0;
	    int x_start_idx=j*col;
	    for(int k=0;k<col;++k)
		new_base_a[i][j]+=base_a[i][k]*static_x[k+x_start_idx];
	}
    for(int i=0;i<row_b;++i)
	for(int j=0;j<new_col;++j)
	{
	    if(objective_type==A_AND_B)
	    {
		new_base_b[i][j]=0.0;
		int x_start_idx=j*col+col*new_col;
		for(int k=0;k<col;++k)
		    new_base_b[i][j]+=base_b[i][k]*static_x[k+x_start_idx];
	    }
	    else
		new_base_b[i][j]=base_b[i][j];
	}
    if(enable_output)
    {
	std::cout<<"A: \n";
	for(int i=0;i<col;++i)
	{
	    for(int j=0;j<new_col;++j)
		std::cout<<static_x[j*col+i]<<" ";
	    std::cout<<"\n";
	}
	if(objective_type==A_AND_B)
	{
	    std::cout<<"B: \n";
	    for(int i=0;i<col;++i)
	    {
		for(int j=0;j<new_col;++j)
		    std::cout<<static_x[j*col+i+new_col*col]<<" ";
		std::cout<<"\n";
	    }
	}
    }

    return 1;
}

void CoupledQuasiHarmonics::getInitialGuess(vector<double> &x)
{
    int new_col=x.size()/(col*2);
    if(objective_type==ONLY_A)//optimize for only A
	new_col=x.size()/col;
    else//optimize for A and B
    {
	//set diagonal of A to 1, and off-diagonal to 0
	for(int i=0;i<col;++i)
	    for(int j=0;j<new_col;++j)
		x[j*col+i]=(i==j)?1:0;
    }
    int bias=col*new_col;//A_AND_B, set elements of B
    if(objective_type==ONLY_A)//ONLY_A, set elements of A
	bias=0; 
    for(int i=0;i<col;++i)
	for(int j=0;j<new_col;++j)
	{
	    if(i!=j)
		x[j*col+i+bias]=0;
	    else
	    {
		double minus_sqr=0,sum_sqr=0;
		for(int k=0;k<p;++k)
		{
		    double F_transpose_mult_base_a_k=0;
		    double G_transpose_mult_base_b_k=0;
		    for(int l=0;l<row_a;++l)
		    {
			double area_weight=(inner_product_type==AREA_WEIGHTED?vertex_area_a[l]:1.0);
			F_transpose_mult_base_a_k+=F[l][k]*base_a[l][i]*area_weight;
		    }
		    for(int l=0;l<row_b;++l)
		    {
			double area_weight=(inner_product_type==AREA_WEIGHTED?vertex_area_b[l]:1.0);
			G_transpose_mult_base_b_k+=G[l][k]*base_b[l][i]*area_weight;
		    }
		    minus_sqr+=(F_transpose_mult_base_a_k-G_transpose_mult_base_b_k)*(F_transpose_mult_base_a_k-G_transpose_mult_base_b_k);
		    sum_sqr+=(F_transpose_mult_base_a_k+G_transpose_mult_base_b_k)*(F_transpose_mult_base_a_k+G_transpose_mult_base_b_k);
		}
		x[j*col+i+bias]=(minus_sqr<=sum_sqr)?1:-1;
	    }
	}
}

void CoupledQuasiHarmonics::testGradients()
{
    //evaluate the gradient of objective
    int num_unknown=2*col*col;
    if(objective_type==ONLY_A)
	num_unknown=col*col;
    double *x=new double[num_unknown];
    double *f_grad=new double[num_unknown];
    F_transpose_base_a_A_minus_G_transpose_base_b_B=new double[p*col];//new_col=col
    srand((unsigned)time(0));
    int random_integer;
    int lowest=1, highest=10;
    int range=(highest-lowest)+1;
    //randomly generate initial x in range [0.1,1]
    for(int i=0;i<num_unknown;++i)
    {
	random_integer = lowest+rand()%range;
	x[i]=random_integer/10.0;
    }
    //analytical gradients
    objectiveGrad(num_unknown,x,f_grad);
    double *perturb=new double[num_unknown];//range: 1e-7~1e-6
    //randomly perturb the unknown and get numerical gradient using finite difference
    for(int i=0;i<num_unknown;++i)
    {
	random_integer = lowest+rand()%range;
	perturb[i]=random_integer/1.0e7;
	x[i]+=perturb[i];
    }
    double f1,f2;
    f2=objectiveVal(num_unknown,x);
    for(int i=0;i<num_unknown;++i)
        x[i]-=2.0*perturb[i];
    f1=objectiveVal(num_unknown,x);;
    double f_grad_times_dx=0;
    for(int i=0;i<num_unknown;++i)
        f_grad_times_dx+=f_grad[i]*2*perturb[i];
    std::cout<<"Objective, df, analytic: "<<setprecision(15)<<f_grad_times_dx<<", numerical: "<<setprecision(15)<<f2-f1;
    std::cout<<", rel_error="<<(f2-f1-f_grad_times_dx)/(fabs(f_grad_times_dx)>1e-20?fabs(f_grad_times_dx):1e-20)<<"\n";
    delete[] x;
    delete[] f_grad;
    delete[] perturb;
    delete[] F_transpose_base_a_A_minus_G_transpose_base_b_B;
    F_transpose_base_a_A_minus_G_transpose_base_b_B=NULL;
}

double CoupledQuasiHarmonics::objectiveVal(int n,const double *x)
{
    int new_col=n/(2*col);
    if(objective_type==ONLY_A)
	new_col=n/col;
    double *lambda_a=eigenvalue_a,*lambda_b=eigenvalue_b;
    double f=0;
    //off-diagonal terms
    for(int i=0;i<new_col;++i)
	for(int j=0;j<new_col;++j)
	{
	    double A_transopse_Lambda_A_ij=0,B_transpose_Lambda_B_ij=0;
	    for(int k=0;k<col;++k)
	    {
		A_transopse_Lambda_A_ij+=lambda_a[k]*x[k+i*col]*x[k+j*col];
		if(objective_type==A_AND_B)
		    B_transpose_Lambda_B_ij+=lambda_b[k]*x[k+i*col+col*new_col]*x[k+j*col+col*new_col];
	    }
	    if(offterm_type==NOT_ORDERED)//no order version
	    {
		if(i!=j)
		{
		    f+=A_transopse_Lambda_A_ij*A_transopse_Lambda_A_ij;
		    if(objective_type==A_AND_B)
			f+=B_transpose_Lambda_B_ij*B_transpose_Lambda_B_ij;
		}
	    }
	    else if(offterm_type==ORDERED)//ordered version
	    {
		if(i==j)
		{
		    A_transopse_Lambda_A_ij-=lambda_a[i];
		    if(objective_type==A_AND_B)
			B_transpose_Lambda_B_ij-=lambda_b[i];
		}
		f+=A_transopse_Lambda_A_ij*A_transopse_Lambda_A_ij;
		if(objective_type==A_AND_B)
		    f+=B_transpose_Lambda_B_ij*B_transpose_Lambda_B_ij;
	    }
	}
    double off_f=f;
    if(enable_output)
	std::cout<<"EvalNum: "<<++f_eval_num<<" Off-term: "<<setprecision(8)<<f<<" ";
    //coupling term
    for(int i=0;i<p;++i)
	for(int j=0;j<new_col;++j)
	{
	    double F_transpose_base_a_A_ij=0,G_transpose_base_b_B_ij=0;
	    for(int k=0;k<col;++k)
	    {
		for(int l=0;l<row_a;++l)
		{
		    double area_weight=(inner_product_type==AREA_WEIGHTED?vertex_area_a[l]:1.0);
		    F_transpose_base_a_A_ij+=F[l][i]*base_a[l][k]*x[j*col+k]*area_weight;
		}
		for(int l=0;l<row_b;++l)
		{
		    double area_weight=(inner_product_type==AREA_WEIGHTED?vertex_area_b[l]:1.0);
		    if(objective_type==A_AND_B)
			G_transpose_base_b_B_ij+=G[l][i]*base_b[l][k]*x[j*col+k+col*new_col]*area_weight;
		    else
			G_transpose_base_b_B_ij+=(j==k?G[l][i]*base_b[l][k]*area_weight:0);
		}
	    }
	    F_transpose_base_a_A_ij*=coupling_scale[j];
	    G_transpose_base_b_B_ij*=coupling_scale[j];
	    f+=mu*(F_transpose_base_a_A_ij-G_transpose_base_b_B_ij)*(F_transpose_base_a_A_ij-G_transpose_base_b_B_ij);
	}
    if(enable_output)
	std::cout<<"Couple-term: "<<setprecision(8)<<f-off_f<<" f: "<<setprecision(8)<<f<<"\n";
    return f;
}

void CoupledQuasiHarmonics::objectiveGrad(int n,const double *x, double *grad)
{
    int new_col=n/(2*col);
    if(objective_type==ONLY_A)
	new_col=n/col;
    double *lambda_a=eigenvalue_a,*lambda_b=eigenvalue_b;
    for(int i=0;i<p;++i)
	for(int j=0;j<new_col;++j)
	{
	    F_transpose_base_a_A_minus_G_transpose_base_b_B[i*new_col+j]=0;
	    for(int k=0;k<col;++k)
		for(int l=0;l<row_a;++l)
		{
		    double area_weight=(inner_product_type==AREA_WEIGHTED?vertex_area_a[l]:1.0);
		    F_transpose_base_a_A_minus_G_transpose_base_b_B[i*new_col+j]+=F[l][i]*base_a[l][k]*x[j*col+k]*area_weight;
		}
	    for(int k=0;k<col;++k)
		for(int l=0;l<row_b;++l)
		{
		    double area_weight=(inner_product_type==AREA_WEIGHTED?vertex_area_b[l]:1.0);
		    if(objective_type==A_AND_B)
			F_transpose_base_a_A_minus_G_transpose_base_b_B[i*new_col+j]-=G[l][i]*base_b[l][k]*x[j*col+k+col*new_col]*area_weight;
		    else
			F_transpose_base_a_A_minus_G_transpose_base_b_B[i*new_col+j]-=(j==k?G[l][i]*base_b[l][k]*area_weight:0);
		}
	    F_transpose_base_a_A_minus_G_transpose_base_b_B[i*new_col+j]*=coupling_scale[j];
	}
    for(int i=0;i<col;++i)
	for(int j=0;j<new_col;++j)
	{
	    //gradient with regard to A
	    grad[j*col+i]=0;
	    for(int k=0;k<col;++k)
		for(int l=0;l<new_col;++l)
		    grad[j*col+i]+=4*(lambda_a[i]*x[l*col+i]*x[l*col+k]*lambda_a[k]*x[j*col+k]);
	    if(offterm_type==NOT_ORDERED)//no order version
	    {	    
		for(int k=0;k<col;++k)
		    grad[j*col+i]-=4*x[j*col+k]*x[j*col+k]*lambda_a[k]*x[j*col+i]*lambda_a[i];
	    }
	    else if(offterm_type==ORDERED)//ordered version
	    {
		grad[j*col+i]-=4*lambda_a[i]*x[j*col+i]*lambda_a[j];
	    }
	    for(int k=0;k<p;++k)
		for(int l=0;l<row_a;++l)
		{
		    double area_weight=(inner_product_type==AREA_WEIGHTED?vertex_area_a[l]:1.0);
		    grad[j*col+i]+=2*mu*base_a[l][i]*F[l][k]*F_transpose_base_a_A_minus_G_transpose_base_b_B[k*new_col+j]*area_weight*coupling_scale[j];
		}
	    if(objective_type==A_AND_B)
	    {
		//gradient with regard to B
		int bias=col*new_col;
		grad[j*col+i+bias]=0;
		for(int k=0;k<col;++k)
		    for(int l=0;l<new_col;++l)
			grad[j*col+i+bias]+=4*(lambda_b[i]*x[l*col+i+bias]*x[l*col+k+bias]*lambda_b[k]*x[j*col+k+bias]);
		if(offterm_type==NOT_ORDERED)//no order version
		{
		    for(int k=0;k<col;++k)
			grad[j*col+i+bias]-=4*x[j*col+k+bias]*x[j*col+k+bias]*lambda_b[k]*x[j*col+i+bias]*lambda_b[i];	
		}
		else if(offterm_type==ORDERED)//ordered version
		{
		    grad[j*col+i+bias]-=4*lambda_b[i]*x[j*col+i+bias]*lambda_b[j];
		}
		for(int k=0;k<p;++k)
		    for(int l=0;l<row_b;++l)
		    {
			double area_weight=(inner_product_type==AREA_WEIGHTED?vertex_area_b[l]:1.0);
			grad[j*col+i+bias]-=2*mu*base_b[l][i]*G[l][k]*F_transpose_base_a_A_minus_G_transpose_base_b_B[k*new_col+j]*area_weight*coupling_scale[j];
		    }
	    }
	}
}

//constraint: A^T*A=I,B^T*B=I ==> A^T*A-I=0,B^T*B-I=0
//m: constraint num, n: variable num
void CoupledQuasiHarmonics::constraintVal(int m, int n, const double *x, double *c)
{
    ColumnVector x_vec(n),c_vec(m);
    for(int i=1;i<=n;++i)
	x_vec(i)=x[i-1];
    Matrix temp_grad(n,m);
    int temp_result;
    optppConstraint(NLPFunction,n,x_vec,c_vec,temp_grad,temp_result);
    for(int i=1;i<=m;++i)
	c[i-1]=c_vec(i);
}

//approximation of constraint gradient with central difference
void CoupledQuasiHarmonics::constraintGrad(int m, int n, const double *x, double *c_grad)
{
    ColumnVector x_vec(n),c_vec(m);
    for(int i=1;i<=n;++i)
	x_vec(i)=x[i-1];
    Matrix temp_grad(n,m);
    int temp_result;
    optppConstraint(NLPGradient,n,x_vec,c_vec,temp_grad,temp_result);
    for(int i=0;i<m;++i)
	for(int j=0;j<n;++j)
	    c_grad[i*n+j]=temp_grad(j+1,i+1);//column order
}

//save input for remote optimization using matlab
void CoupledQuasiHarmonics::saveDataForMatlab(const char *file_name, vector<double> &x0, int row_a, int row_b, int col,
                                              double *eigenvalue_a, double *eigenvalue_b, double mu, double **F, double **G,
                                              int p, double **base_a, double **base_b, double *vertex_area_a, double *vertex_area_b)
{
    fstream output_file(file_name,ios::out);
    if(!output_file)
    {
	std::cout<<"Error: Cannot open file "<<file_name<<"!\n";
	return;
    }
    output_file<<"% Input data for Coupled Quasi-harmonic Bases matlab implementation.\n";
    output_file<<"% Generated by exampleBasedDeformableSimulator.\n\n";
    output_file<<"function [Phi DX LambdaX Psi DY LambdaY F G mu Kprime K] = MatlabSolverInputs()\n";
    //Phi: base_a
    output_file<<"Phi=[ ";
    for(int i=0;i<row_a;++i)
	for(int j=0;j<col;++j)
	{
	    output_file<<base_a[i][j];
	    if(i!=row_a-1&&j==col-1)
		output_file<<"; ";
	    else
		output_file<<" ";
	}
    output_file<<"];\n";
    //DX: diagonal matrix with entries as vertex_area_a
    output_file<<"dx=[ ";
    for(int i=0;i<row_a;++i)
	output_file<<vertex_area_a[i]<<" ";
    output_file<<"];\n";
    output_file<<"DX=sparse(1:"<<row_a<<",1:"<<row_a<<",dx);\n";
    output_file<<"clear dx;\n";
    //LambdaX: eigenvalue_a
    output_file<<"LambdaX=[ ";
    for(int i=0;i<col;++i)
	output_file<<eigenvalue_a[i]<<" ";
    output_file<<"];\n";
    //Psi: base_b
    output_file<<"Psi=[ ";
    for(int i=0;i<row_b;++i)
	for(int j=0;j<col;++j)
	{
	    output_file<<base_b[i][j];
	    if(i!=row_b-1&&j==col-1)
		output_file<<"; ";
	    else
		output_file<<" ";
	}
    output_file<<"];\n";
    //DY: diagonal matrix with entries as vertex_area_b
    output_file<<"dy=[ ";
    for(int i=0;i<row_b;++i)
	output_file<<vertex_area_b[i]<<" ";
    output_file<<"];\n";
    output_file<<"DY=sparse(1:"<<row_b<<",1:"<<row_b<<",dy);\n";
    output_file<<"clear dy;\n";
    //LambdaY: eigenvalue_b
    output_file<<"LambdaY=[ ";
    for(int i=0;i<col;++i)
	output_file<<eigenvalue_b[i]<<" ";
    output_file<<"];\n";
    //F && G
    output_file<<"F = zeros( "<<row_a<<", "<<p<<");\n";
    for(int i=0;i<row_a;++i)
	for(int j=0;j<p;++j)
	    if(F[i][j]!=0)
		output_file<<"F("<<i+1<<","<<j+1<<")="<<F[i][j]<<";\n";
    output_file<<"G = zeros( "<<row_b<<", "<<p<<");\n";
    for(int i=0;i<row_b;++i)
	for(int j=0;j<p;++j)
	    if(G[i][j]!=0)
		output_file<<"G("<<i+1<<","<<j+1<<")="<<G[i][j]<<";\n";
    //mu
    output_file<<"mu = "<<mu<<";\n";
    //Kprime
    output_file<<"Kprime = "<<col<<";\n";
    //K
    output_file<<"K = "<<x0.size()/(2*col)<<";\n";
    output_file<<"end\n";
    output_file.close();
}

//load optimization result generated by matlab
void CoupledQuasiHarmonics::loadResultFromMatlab(const char *file_name, vector<double> &x)
{
    //result format: elements of A followed by elements of B
    int new_col=x.size()/(2*col);
    fstream input_file(file_name);
    if(!input_file)
    {
	std::cout<<"Error: Cannot open file "<<file_name<<"!\n";
	return;
    }
    for(int i=0;i<new_col;++i)
	for(int j=0;j<col;++j)
	    input_file>>x[j*new_col+i];
    for(int i=0;i<new_col;++i)
	for(int j=0;j<col;++j)
	    input_file>>x[j*new_col+i+col*new_col];
}

//****************************************************NLopt stuff**************************************************

double nloptObjective(unsigned int n, const double *x, double *grad, void *my_func_data)
{
    CoupledQuasiHarmonics *data=(CoupledQuasiHarmonics*)my_func_data;

    if(grad)
	data->objectiveGrad(n,x,grad);
    double f=data->objectiveVal(n,x);
    //workaround for NLopt bug: if the solver is terminated due to max evaluation times is reached, NLopt simply
    //returns initial guess. So we store the intemediate solution after every iteration
    memcpy(data->temp_x,x,sizeof(double)*n);
    data->temp_f=f;
    return f;
}

void nloptConstraint(unsigned int m, double *result, unsigned int n, const double *x, double *grad, void *my_func_data)
{
    CoupledQuasiHarmonics *data=(CoupledQuasiHarmonics*)my_func_data;
    if(grad)
	data->constraintGrad(m,n,x,grad);
    data->constraintVal(m,n,x,result);
}

//*****************************************************************************************************************

//****************************************************Opt++ stuff**************************************************
void initOptpp(int n, ColumnVector &x)
{
    //CoupledQuasiHarmonic::static_x stores the initial guess
    for(int i=0;i<n;++i)
	x(i+1)=CoupledQuasiHarmonics::static_x[i];
}

void optppObjective(int mode, int n, const ColumnVector &x, double &f_x, ColumnVector &grad_x, int &result, void *my_func_data)
{
    CoupledQuasiHarmonics *data=(CoupledQuasiHarmonics*)my_func_data;
    double *temp_x=data->temp_x;//temp_x is used as variable x here
    double *temp_grad=data->temp_f_grad;
    for(int i=0;i<n;++i)
	temp_x[i]=x(i+1);
    //function evaluation
    if(mode&NLPFunction)
    {
	f_x=data->objectiveVal(n,temp_x);
	result=NLPFunction;
    }
    //gradient evaluation
    if(mode&NLPGradient)
    {
	data->objectiveGrad(n,temp_x,temp_grad);
	for(int i=0;i<n;++i)
	    grad_x(i+1)=temp_grad[i];
	result=NLPGradient;
    }
}

//constraint, gradient approximated with central difference
//opt++ has bug in approximating gradient with finite difference, so I compute it myself
void optppConstraint(int mode, int n, const ColumnVector &x, ColumnVector &c_x, Matrix &c_grad_x, int &result)
{
    //constraint: A^T*A=I,B^T*B=I ==> A^T*A-I=0,B^T*B-I=0
    //Note: A^T*A and B^T*B are both symmetric matrices, inactive constraints need to be removed
    //      otherwise it leads to singular jacobian of optimization system
    int col=CoupledQuasiHarmonics::static_col;
    int new_col=n/(2*col);
    if(CoupledQuasiHarmonics::static_objective_type==CoupledQuasiHarmonics::ONLY_A)
	new_col=n/col;
    //function evaluation
    if(mode&NLPFunction)
    {
	int idx_1d=1;
	for(int i=1;i<=new_col;++i)
	    for(int j=1;j<=i;++j)
	    {
		c_x(idx_1d)=0;
		for(int k=1;k<=col;++k)
		    c_x(idx_1d)+=x((i-1)*col+k)*x((j-1)*col+k);
		if(i==j)
		    c_x(idx_1d)-=1.0;
		++idx_1d;
	    }
	if(CoupledQuasiHarmonics::static_objective_type==CoupledQuasiHarmonics::A_AND_B)
	{
	    int bias=new_col*col;
	    for(int i=1;i<=new_col;++i)
		for(int j=1;j<=i;++j)
		{
		    c_x(idx_1d)=0;
		    for(int k=1;k<=col;++k)
			c_x(idx_1d)+=x((i-1)*col+k+bias)*x((j-1)*col+k+bias);
		    if(i==j)
			c_x(idx_1d)-=1.0;
		    ++idx_1d;
		}
	}
	result=NLPFunction;
    }
    //gradient evaluation
    if(mode&NLPGradient)
    {
	ColumnVector temp_x=x;
	ColumnVector fplus(c_x.size()),fminus(c_x.size());
	Matrix c_grad_x_t(c_x.size(),n);
	Matrix temp_grad(n,c_x.size());
	int temp_result;
	double perturb=1.0e-8;
	for(int i=1;i<=n;++i)
	{
	    temp_x(i)+=perturb;
	    optppConstraint(NLPFunction,n,temp_x,fplus,temp_grad,temp_result);
	    temp_x(i)-=2*perturb;
	    optppConstraint(NLPFunction,n,temp_x,fminus,temp_grad,temp_result);
	    c_grad_x_t.Column(i)=(fplus-fminus)/(2*perturb);
	    temp_x(i)=x(i);
	}
	c_grad_x=c_grad_x_t.t();
	result=NLPGradient;
    }
}  
//*****************************************************************************************************************
