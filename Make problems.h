	#include <iostream> 	
	#include <ilcplex/cplex.h>
///the following functions are used in Branch-and-Cut with postive rewards

	extern int make_master3
		(CPXCENVptr env,CPXCENVptr env_slp,CPXLPptr lp, CPXLPptr* slp, CPXLPptr* slp_v, double ***tran_prob, double **rhs_contin,double **rhs_stop, double **v,double **d, int *obj_ind);		
	
	extern int make_master2
		(CPXCENVptr env,CPXCENVptr env_slp,CPXLPptr lp, CPXLPptr* slp, CPXLPptr* slp_v, double ***tran_prob, double **rhs_contin,double **rhs_stop, double **v,double **d, int *obj_ind);		

	extern int make_slpvt3(CPXCENVptr env_slp, CPXLPptr slp_t, CPXLPptr slp_vt, double ***tran_prob, double **rhs_contin,double **rhs_stop, double *vt,double **d, int *obj_ind);
	extern int make_slpvt2(CPXCENVptr env_slp, CPXLPptr slp_t, CPXLPptr slp_vt, double ***tran_prob, double **rhs_contin,double **rhs_stop, double *vt,double **d, int *obj_ind);

///the following functions are used in CPLEX with postive rewards

	extern int parameters_calulator_2
		(CPXCENVptr env_slp, CPXLPptr* slp, CPXLPptr* slp_v, double ***tran_prob, double **rhs_contin,double **rhs_stop, double **v,double **d, int *obj_ind);

	extern int parameters_calulator_3
		(CPXCENVptr env_slp, CPXLPptr* slp, CPXLPptr* slp_v, double ***tran_prob, double **rhs_contin,double **rhs_stop, double **v,double **d, int *obj_ind);		
	
	extern int make_extensive_2(CPXCENVptr env_slp, CPXLPptr exten, double ***tran_prob, double **rhs_contin,double **rhs_stop, double **v,double **d, int *obj_ind);
	extern int make_extensive_3
		(CPXCENVptr env_slp, CPXLPptr exten, double ***tran_prob, double **rhs_contin,double **rhs_stop, double **v,double **d, int *obj_ind);


///the following functions are used in both Branch-and-cut and CPLEX with postive rewards

	extern int readdata3( double ***, double **,double **, int *, int *);

	extern int readdata2( double ***tran_prob, double **reward,double **stop_reward, int *obj_ind, int *instance_num);


///the following functions are used in Branch-and-Cut with negative rewards

	extern int make_master_negative(CPXCENVptr env,CPXCENVptr env_slp,CPXLPptr lp, CPXLPptr* slp, CPXLPptr* slp_v, double ***tran_prob, double **rhs_contin,double **rhs_stop, double **v,double **d, int *obj_ind);		

	extern int make_slpvt_negative(CPXCENVptr env_slp, CPXLPptr slp_t, CPXLPptr slp_vt, double ***tran_prob, double **rhs_contin,double **rhs_stop, double *vt,double **d, int *obj_ind);
		

///the following functions are used in CPLEX with negative rewards

	extern int parameters_calulator_negative
		(CPXCENVptr env_slp, CPXLPptr* slp, CPXLPptr* slp_v, double ***tran_prob, double **rhs_contin,double **rhs_stop, double **v,double **d, int *obj_ind);		

	
	extern int make_extensive_negative
		(CPXCENVptr env_slp, CPXLPptr exten, double ***tran_prob, double **rhs_contin,double **rhs_stop, double **v,double **d, int *obj_ind);


///the following function is used in both Branch-and-cut and CPLEX with negative rewards
	
	extern int readdata_negative( double ***tran_prob, double **reward,double **stop_reward, int *obj_ind, int *instance_num);