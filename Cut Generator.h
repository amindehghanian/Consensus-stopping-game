	#include <iostream> 	
	#include <ilcplex/cplex.h>

	#define BRANCH_AND_CUT 1
	#define CPLEX 2
	#define REWARD_POSITIVE 1
	#define REWARD_NEGATIVE -1

	const int num_patients=2;//number of patients	
	const int num_states=40; //number of states
	const int solution_method=BRANCH_AND_CUT; //Let solution_method be BRANCH_AND_CUT (or 1) for the Branch-and-Cut approach and CPLEX (or 2) for the CPLEX approach.
	const int reward_type=REWARD_POSITIVE;//Let reward_type be REWARD_POSITIVE (or 1) for positive (non-negative) rewards and REWARD_NEGATIVE (or -1) for negative rewards
	const int instance_beg=1; //This indicates the first instance number which should be solved be the solution method
	const int instance_end=30; //This indicates the last instance number which should be solved be the solution method
	const double time_limit=14400; //This indicates the time limiti for the selected solution method
	const int num_cores=1; //This indicates the number of threads is used to solve instances in parallel


	struct cutinfo {
	CPXENVptr env_slp;
	CPXLPptr slp_t;
	CPXLPptr* slp;
	CPXLPptr slp_vt;
	CPXLPptr* slp_v;	
	double **rhs_stop;
	double **rhs_contin;	
	double **v; //to store those big M's
	double *vt;
	double **v_upd;
	double *vt_upd;
	double **d; 	
	double *feas_h;			
	int obj_ind;
	int num_usercuts;
	bool *feas_h_status;
	int last_node;	
	};
	typedef struct cutinfo CUTINFO, *CUTINFOptr;

	/* The following structure will hold the information we need to 
	pass to the branch callback function */
	struct branchinfo {
	CPXENVptr env_slp;	//you need to be carfeul about this, since it may have been used at the same time in mycutcallback	
	CPXLPptr* slp;	
	double **rhs_stop;
	double *last_heur_solu;
	double *feas_h;	
	int obj_ind;
	bool *feas_h_status;	
	};
	typedef struct branchinfo BRANCHINFO, *BRANCHINFOptr;

	extern int CPXPUBLIC 
	mycutcallback (CPXCENVptr env, void *cbdata, int wherefrom,void *cbhandle, int *useraction_p);
	
	extern int CPXPUBLIC 
	usercutcallback (CPXCENVptr env, void *cbdata, int wherefrom,void *cbhandle, int *useraction_p);	

	extern int CPXPUBLIC 
	myincumbentcheck (CPXCENVptr env,void*cbdata,int wherefrom,void*cbhandle,double objval,double*x,int *isfeas_p,int *useraction_p);

	extern int CPXPUBLIC mybranchfunc(CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle, int type, int sos, //This is used with CPLEX 12.6
		int nodecnt, int bdcnt, const int *nodebeg, const int *indices, const char *lu, const double *bd, const double *nodeest, int *useraction_p);

	//extern int CPXPUBLIC mybranchfunc (CPXCENVptr env,void *cbdata,int wherefrom,void *cbhandle,int type,int sos, //This is used with CPLEX 12.4
	//int nodecnt,int bdcnt,const double *nodeest,const int *nodebeg,const int *indices,const char *lu,const int *bd,int *useraction_p);

	extern int CPXPUBLIC myheuristic (CPXCENVptr env,void *cbdata,int wherefrom, 
	void *cbhandle,double *objval_p,double *x,int *checkfeas_p,int *useraction_p);
	extern int v_update(CPXCENVptr env_slp, CPXLPptr* slp_v, double *lb, double *ub, double **v, double **rhs_stop, double **v_upd);

	extern int vt_update(CPXCENVptr env_slp, CPXLPptr slp_vt, double *lb, double *ub, double *vt,double **rhs_stop, double *vt_upd);
	

	extern int ub_update(CPXCENVptr env_slp, CPXLPptr* slp, double *lb, double *ub, double **v, double **rhs_stop, int *obj_ind);	

	extern int ubt_update(CPXCENVptr env_slp, CPXLPptr slp_t, double *lb, double *ub, double **v_upd, double *vt_upd, double **rhs_stop, int *obj_ind);

	extern bool feasibility(CPXCENVptr env_slp, CPXLPptr* slp, CPXLPptr* copy_slp, double *x, double **y,int **sol_state,double **rhs_stop);

	extern bool feasibility_strong(CPXCENVptr env_slp, CPXLPptr* slp, CPXLPptr* slp_v, double **y, double *lb, double *x,int **sol_state,double **rhs_stop);

	extern bool feasibility_check(CPXCENVptr env_slp, CPXLPptr* slp, CPXLPptr* copy_slp, double *x, double **solu,double **rhs_stop);	

	
	#define EPSILON 1E-4
	
	const double lambda=0.998829178;
	extern int counter; // We want "counter" to be a global variable
	extern int error_counter;
	
	

	extern int numz;
	extern int status;
	const int lb_local_search=0;
	const int ub_local_search=num_states;
		
		
	

