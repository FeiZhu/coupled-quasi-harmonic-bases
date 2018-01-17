//######################################################
// Implementation of EUROGRAPHICS 2013 paper:
// <Coupled Quasi-harmnoic Bases>
// @author: Fei Zhu, 01/14/2014
//######################################################

#ifndef COUPLED_QUASI_HARMONIC_BASES_H_
#define COUPLED_QUASI_HARMONIC_BASES_H_

#include <vector>
#include "NLF.h"
using std::vector;
using NEWMAT::ColumnVector;
using NEWMAT::Matrix;

//nlopt stuff
double nloptObjective(unsigned int n, const double *x, double *grad, void *my_func_data);
void nloptConstraint(unsigned int m, double *result, unsigned int n, const double *x, double *grad, void *my_func_data);
//opt++ stuff
void optppObjective(int mode, int n, const ColumnVector &x, double &f_x, ColumnVector &grad_x, int &result, void *my_func_data);
void optppConstraint(int mode, int n, const ColumnVector &x, ColumnVector &c_x, Matrix &c_grad_x, int &result);
void initOptpp(int n, ColumnVector &x);  

class CoupledQuasiHarmonics
{
public:
    //@base_a,base_b: pointer to the input eigenfunctions of shapes, size: num_vert X num_eigen, each column is an eigenfunction
    //@eigenvalue_a,eigenvalue_b: pointer to the input eigenvalues of shapes
    //@row_a,row_b: vertex numbers of the input shapes
    //@col: number of input eigenfunctions for each shape
    //@vertex_area_a,vertex_area_b: the vertex areas of the input shapes
    //@F,G: corresponding functions between shapes
    //@p: number of corresponding functions
    //@Note: all storage of pointers are allocated by caller
    CoupledQuasiHarmonics(double **input_base_a, double **input_base_b, double *input_eigenvalue_a, double *input_eigenvalue_b,int input_row_a, int input_row_b, int input_col,
                          double *input_vertex_area_a, double *input_vertex_area_b, double **input_F, double **input_G,int input_p);
    ~CoupledQuasiHarmonics();
    //@new_col: the number of coupled eigenfunctions for each shape, must be smaller than input 
    //@new_base_a,new_base_b: pointer to the new eigenfunctions, space allocated by the caller, each column is an eigenfunction
    //@return_value: if <0, failure
    int getCoupledBases(int new_col, double **new_base_a, double **new_base_b);
    void setMu(double new_mu){mu=new_mu;}
    double getMu(){return mu;}
    void setPackage(int package){package_type=package;}
    void setOfftermType(int offterm){offterm_type=(offterm>0?ORDERED:NOT_ORDERED);}
    void setInnerProductType(int inner_product){inner_product_type=(inner_product>0?AREA_WEIGHTED:STANDARD);}
    void setObjectiveType(int objective){objective_type=(objective>0?ONLY_A:A_AND_B);static_objective_type=objective_type;}//decide whether to optimize for A&B or only A
    void enableScaleBase(){scale_base=true;}//enable base_a && base_b to same scale
    void disableScaleBase(){scale_base=false;}
    void enableScaleCorrespondingFunction(){scale_corresponding_function=true;}//scale columns of corresponding function with area of that region
    void disableScaleCorrespondingFunction(){scale_corresponding_function=false;}
    void setCouplingScale(const vector<double> &input_coupling_scale){coupling_scale=input_coupling_scale;}//assume input_coupling_scale is of correct size (new_col)
    void enableOutput(){enable_output=true;}
    void disableOutput(){enable_output=false;}

    //***********************************************NLopt stuff********************************************************
    friend double nloptObjective(unsigned int n, const double *x, double *grad, void *my_func_data);
    friend void nloptConstraint(unsigned int m, double *result, unsigned int n, const double *x, double *grad, void *my_func_data);
    //******************************************************************************************************************

    //***********************************************Opt++ stuff********************************************************
    friend void optppObjective(int mode, int n, const ColumnVector &x, double &f_x, ColumnVector &grad_x, int &result, void *my_func_data);
    friend void optppConstraint(int mode, int n, const ColumnVector &x, ColumnVector &c_x, Matrix &c_grad_x, int &result);  
    friend void initOptpp(int n, ColumnVector &x);
    //******************************************************************************************************************

    void testGradients();//test if the derivation of gradients of objective function is correct
public:
    enum PACKAGES{
	NLOPT=0,
	OPTPP=1,
	MATLAB=2};
    enum OFFTERM{
	NOT_ORDERED=0,
	ORDERED=1};
    enum INNERPRODUCT{
	STANDARD=0,
	AREA_WEIGHTED=1};
    enum OBJECTIVE{
	A_AND_B=0,
	ONLY_A=1};
protected:
    void getInitialGuess(vector<double> &x);
    double objectiveVal(int n,const double *x);
    void objectiveGrad(int n,const double *x,double *grad);
    void constraintVal(int m,int n,const double *x,double *c);
    void constraintGrad(int m,int n,const double *x,double *c_grad);
protected://matlab stuff
    void saveDataForMatlab(const char *file_name,vector<double> &x0,int row_a,int row_b,int col,double *eigenvalue_a,double *eigenvalue_b,
                           double mu,double **F,double **G,int p,double **base_a,double **base_b,double *vertex_area_a,double *vertex_area_b);
    void loadResultFromMatlab(const char *file_name,vector<double> &x);
protected:
    double **base_a, **base_b, *eigenvalue_a, *eigenvalue_b,*vertex_area_a,*vertex_area_b,**F,**G;
    double base_a_mean,base_b_mean;
    int row_a,row_b,col,p;
    double mu;
    double *temp_x,temp_f,*temp_f_grad,*F_transpose_base_a_A_minus_G_transpose_base_b_B;//temp buffers
    int f_eval_num;
    int package_type;
    int offterm_type;
    int objective_type;
    bool scale_base;//scale base_a and base_b to the same scale
    bool scale_corresponding_function;// in case if binary functions representing corresponding regions are provided, 
                                      // F and G can be normalized by factor which approximately equals to the area of the region
    int inner_product_type;//decide whether inner product of F and Phi needs to be area weighted
    vector<double> coupling_scale;//set different coupling strength for different eigenfunctions, default all 1
    bool enable_output;
    //hack for opt++, since it has limitation passing custom data
    static int static_col;
    static int static_objective_type;
    static vector<double> static_x;
};

#endif
