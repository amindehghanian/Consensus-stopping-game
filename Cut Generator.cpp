
#pragma region libraries 	
	#include <iostream> 	
	#include <ilcplex/cplex.h>	
	#include "Cut Generator.h"

	#include <ctype.h>
	#include <string.h>
	#include <math.h>
	#include <algorithm>   
	#include <stdio.h>// standard input output to use printer and other device 
	#include <stdlib.h>//C Standard General Utilities Library. This header defines several general purpose functions, including dynamic memory management, random number generation, communication with the environment, integer arthmetics, searching, sorting and converting.
	#include <iostream> // 
	#include <fstream> // to read or write in a file
	#include <iomanip>
	#include <ctime>
	#include <time.h>
	#include <sstream>
	#include <vector> 
	#include <algorithm>	
	using namespace std;	
#pragma endregion

	int CPXPUBLIC 
	mycutcallback (CPXCENVptr env, void *cbdata, int wherefrom,void *cbhandle, int *useraction_p){
	
#pragma region intializing data structures
	CUTINFOptr cutinfo = (CUTINFOptr) cbhandle;
	CPXCENVptr env_slp=cutinfo->env_slp;	
	CPXLPptr slp_t=cutinfo->slp_t;
	CPXLPptr* slp=cutinfo->slp;
	CPXLPptr slp_vt=cutinfo->slp_vt;
	CPXLPptr* slp_v=cutinfo->slp_v;
	double **v=cutinfo->v;
	double *vt=cutinfo->vt;
	double **rhs_contin=cutinfo->rhs_contin;
	double **rhs_stop=cutinfo->rhs_stop;
	double **v_upd= cutinfo->v_upd; // to store updated V's
	double *vt_upd= cutinfo->vt_upd; 
	double **d= cutinfo->d;	
	
	
	
	double *feas_h=cutinfo->feas_h;	
		
	int obj_ind=cutinfo->obj_ind;
	bool *feas_h_status=cutinfo->feas_h_status;		
	int *last_node=&cutinfo->last_node;
	int *num_usercuts=&cutinfo->num_usercuts;

	double   *x= new double [numz+num_patients];
	double   *lb= new double [numz];
	double   *ub= new double [numz+num_patients];
	int nz=0;// this variable counts the number of binary variables with the value 1
	double *cutval=new double[numz+num_patients];		
	double cutrhs;
	
	int *index= new int [numz+num_patients];
	for (int k=0;k<numz+num_patients;k++)
		index[k]=k;
	double **v_prim= new double *[num_patients];
	for (int k=0;k<num_patients;k++)
		v_prim[k]=new double [numz];

	*useraction_p = CPX_CALLBACK_DEFAULT; 	

	///I use the following pointers just in my first algorithm
	double **d_prim= new double *[num_patients];		
	double **y=new double *[ num_patients];
	double **pi=new double *[ num_patients];
	int **sol_state= new int *[num_patients];// to keep state of current solution  
	double *objval=new double [num_patients];
	for(int k=0; k<num_patients; k++){
	d_prim[k]=new double [numz];
	y[k]= new double [numz];
	pi[k]= new double [numz];
	sol_state[k]=new int [numz+4]; // before last one is to store the number of violated states with x(s)=1 while last one is for total number of violated states
	}
	for (int i=0;i<num_patients;i++)
	for (int j=0;j<numz+4;j++)
		sol_state[i][j]=0;

	for (int i=0;i <num_patients;i++)
	for (int j=0;j <numz;j++){
		v_prim[i][j]=v_upd[i][j]; 						
		d_prim[i][j]=d[i][j];
	}

	int beg=0;
	status = CPXgetcallbacknodex (env, cbdata, wherefrom, x,0, numz+num_patients-1); 
	if ( status ) {
		fprintf(stderr, "Failed to get node solution.\n");
		goto TERMINATE;
	}	
	
	int numiinf;
	status = CPXgetcallbacknodeinfo(env,cbdata,wherefrom,0,CPX_CALLBACK_INFO_NODE_SEQNUM,&numiinf);
	if ( status ) {
	fprintf(stderr, "Failed to get node solution.\n");
	goto TERMINATE;
	}	

	double gap;
	double best_int;
	double best_ub;
	status=CPXgetcallbackinfo(env, cbdata, wherefrom, CPX_CALLBACK_INFO_BEST_INTEGER, &best_int);
	status=CPXgetcallbackinfo(env, cbdata, wherefrom, CPX_CALLBACK_INFO_BEST_REMAINING, &best_ub);
	gap=(best_ub-best_int)/best_int;

	CPXLPptr copy_slp[num_patients];
	for (int i=0;i<num_patients;i++)
	{
		copy_slp[i] = CPXcloneprob (env_slp, slp[i], &status); // this creates a copy of slp problem 
		if ( status ) goto TERMINATE;
	}		
#pragma endregion
	//////////////////////////////////// checking integer feasibility of solution of the current node, and aggressively generating cuts		
	for(int i=0; i<numz; i++){				
		if (x[i] > (1.0 - EPSILON)){
			x[i] = 1.0;	
			nz++;
		}
		else
			x[i] = 0.0;			
	}



#pragma region Dynamically updating V 

	int node_level;	
	status = CPXgetcallbacknodeinfo(env,cbdata,wherefrom,0,CPX_CALLBACK_INFO_NODE_DEPTH,&node_level);
	status = CPXgetcallbacknodelb (env, cbdata, wherefrom,lb, 0, numz-1);
	status = CPXgetcallbacknodeub (env, cbdata, wherefrom,ub, 0, numz+num_patients-1);

	if( numiinf!=last_node[0]){		
		
		v_update(env_slp, slp_v, lb, ub, v, rhs_stop,v_upd);
		vt_update(env_slp, slp_vt, lb, ub, vt, rhs_stop,vt_upd);						
		ubt_update(env_slp,slp_t,lb,ub,v_upd,vt_upd,rhs_stop, &obj_ind);
		ub_update(env_slp,slp, lb, ub, v_upd, rhs_stop,&obj_ind);

		for(int i=0;i<num_patients;i++){				
				
		for(int j=0; j<numz+num_patients; j++)
				cutval[j]=0;		

		cutval[numz+i]=1;
		cutrhs=v_upd[i][obj_ind];				
		status = CPXcutcallbackaddlocal (env, cbdata, wherefrom,numz+num_patients, cutrhs, 'L',index, cutval);												
		if ( status ) goto TERMINATE;	
		num_usercuts[0]++;
		
		}

		for (int i=0;i <num_patients;i++)
		for (int j=0;j <numz;j++)
			v_prim[i][j]=v_upd[i][j];	

		for(int j=0; j<numz; j++)
				cutval[j]=0;	
		for(int j=0; j<num_patients; j++)
		cutval[numz+j]=1;
		cutrhs=vt_upd[obj_ind];						

		status = CPXcutcallbackaddlocal (env, cbdata, wherefrom,numz+num_patients, cutrhs, 'L',index, cutval);												
		if(status) goto TERMINATE;
		num_usercuts[0]++;
	}

#pragma endregion 


#pragma region	adding constraints to master problem in the root node
	if( numiinf!=last_node[0]&&numiinf==0){		
		double* init_solu= new double [numz];
		bool start_feas;
		for(int j=0;j<numz;j++){
			start_feas=false;
			for(int i=0;i<num_patients;i++)
				if(d[i][j]>rhs_stop[i][j]+EPSILON){
					start_feas=true;
					break;
				}

			if(start_feas==true)
				init_solu[j]=0;
			else
				init_solu[j]=1;			
		}

		CPXLPptr copy_slpd[num_patients];
		for (int i=0;i<num_patients;i++)
		{
			copy_slpd[i] = CPXcloneprob (env_slp, slp[i], &status); // this creates a copy of slp problem 
			if ( status ) goto TERMINATE;
		}
		
		feasibility(env_slp,slp,copy_slpd,init_solu,y,sol_state,rhs_stop);
		for(int i=0;i<num_patients;i++)
		if ( copy_slpd[i] != NULL ) {
			status = CPXfreeprob (env_slp, &copy_slpd[i]);
			if ( status ) {
				fprintf (stderr, "CPXfreeprob failed , error code %d.\n",
						status);
			}
		}
		
		int sol_state_sum=0;
		for(int i=0;i<num_patients;i++)
				sol_state_sum+=sol_state[i][numz+3];

		while(sol_state_sum>=1){		
			
			for(int j=0; j<numz; j++){	
				start_feas=true;
				for(int i=0;i<num_patients;i++)
				if(sol_state[i][j]!=1){
					start_feas=false;
					break;
				}

				if(start_feas==true) 
					init_solu[j]=1; 
				else
					init_solu[j]=0;
				}

			for (int i=0;i<num_patients;i++)
			{
				copy_slpd[i] = CPXcloneprob (env_slp, slp[i], &status); // this creates a copy of slp problem 
				if ( status ) goto TERMINATE;
			}
			feasibility(env_slp,slp,copy_slpd,init_solu,y,sol_state,rhs_stop);
			for(int i=0;i<num_patients;i++){		
				if(sol_state[i][numz+2]+sol_state[i][numz+3]>0){	
					for(int j=0; j<numz; j++){			
						if ((sol_state[i][j]==2) ||sol_state[i][j]==3)											
							cutval[j]=1;			
						else if(sol_state[i][j]==1) 			
							cutval[j]=sol_state[i][numz+2]+sol_state[i][numz+3];	 
						else
							cutval[j]=0;
					}	
					for (int j=numz;j<numz+num_patients;j++) 
						cutval[j]=0;

					cutrhs=(sol_state[i][numz+2]+sol_state[i][numz+3])*(sol_state[i][numz+1]);
		
					status = CPXcutcallbackaddlocal (env, cbdata, wherefrom,numz+num_patients, cutrhs, 'L',index, cutval);												
					}		
				}
					
			for(int i=0;i<num_patients;i++)
			if ( copy_slpd[i] != NULL ) {
				status = CPXfreeprob (env_slp, &copy_slpd[i]);
				if ( status ) {
					fprintf (stderr, "CPXfreeprob failed , error code %d.\n",
							status);
				}
			}

			sol_state_sum=0;
			for(int i=0;i<num_patients;i++)
				sol_state_sum+=sol_state[i][numz+3];
		
		}

		delete[] init_solu;
	}
#pragma endregion	



	///////////////////////////Add the my first algorithm cuts in the following/////////////////////////////		
	//double **v_x= new double *[num_patients];
	//for(int i=0;i<num_patients;i++)
	//	v_x[i]= new double[numz];
	//if(feasibility_strong(env_slp,slp,slp_v,v_x,lb,x,sol_state,rhs_stop)==false) goto ADD_FEASIBILITY_CUT;

	if(feasibility(env_slp, slp, copy_slp, x, y,sol_state,rhs_stop)==false) goto ADD_FEASIBILITY_CUT;	

# pragma region Optimality cut is added in the following		

	for(int k=0; k<num_patients; k++)
	{		
		status = CPXgetobjval (env_slp, copy_slp[k], &objval[k]);
		if ( status ) { fprintf (stderr, "Failed to obtain objval.\n");goto TERMINATE;}

		if(x[numz+k]>objval[k]+EPSILON){
				status = CPXgetpi (env_slp, copy_slp[k], pi[k], 0, numz-1); 
				if ( status ) {fprintf (stderr, "Failed to obtain dual variables.\n"); goto TERMINATE;}				
				cutrhs=0;
				for(int j=0; j<numz; j++)												   
				{
					if (x[j]==0)
					{	
						cutrhs +=pi[k][j] * rhs_contin[k][j];											
						if(rhs_stop[k][j]-d_prim[k][j]>0)
						cutval[j]=-pi[k][j] * (rhs_stop[k][j]-d_prim[k][j]);									
						else
						cutval[j]=-pi[k][j]*rhs_stop[k][j];																
					}					
					else if(x[j]==1 &&lb[j]<EPSILON)
					{	
						double infeasout;
						status = CPXgetrowinfeas (env_slp, slp[k], v_prim[k], &infeasout, j, j);
						infeasout= v_prim[k][j]-infeasout;
						cutrhs +=pi[k][j] * infeasout;
						cutval[j]=-pi[k][j] * (rhs_stop[k][j]-infeasout);							
					}
					else 
					{							
						cutval[j]=0;							
						cutrhs +=pi[k][j] *rhs_stop[k][j];
					}
				}
				for (int j=numz;j<numz+num_patients;j++) 
				cutval[j]=0;

				cutval[numz+k]=1;
							
				status = CPXcutcallbackaddlocal (env, cbdata, wherefrom, numz+num_patients, cutrhs, 'L',index, cutval);				
				if ( status )	goto TERMINATE;	
				num_usercuts[0]++;
				*useraction_p = CPX_CALLBACK_SET; 
				
			}
		}

#pragma endregion 

# pragma region adding the new point found to my heuristic
	double inc_obj=0;
	double feas_obj=0;
	for(int i=0;i<num_patients;i++){
		 inc_obj +=y[i][obj_ind];
		 feas_obj +=feas_h[numz+i];
	}
	if(inc_obj>feas_obj){
		for(int j=0;j<numz;j++)
			feas_h[j]=x[j];

		for(int i=0;i<num_patients;i++)
			feas_h[numz+i]=y[i][obj_ind];
		feas_h_status[0]=true;
	}

#pragma endregion 
	
	
 #pragma region 	stronger (newer) form of aggregated combinatorial optimality cut

	for(int j=0; j<numz; j++){
	if(sol_state[0][j]==0||sol_state[1][j]==0)
		cutval[j]=1;
	else
		cutval[j]=0;		
	}

	for (int k=numz;k<numz+num_patients;k++) 
		cutval[k]=0;
	cutrhs=1;
	
	status = CPXcutcallbackaddlocal (env, cbdata, wherefrom,numz+num_patients, cutrhs, 'G',index, cutval); // this cut should be local since it is redundant in the other nodes	  					
	if ( status ) 	goto TERMINATE;	

	num_usercuts[0]++;

#pragma endregion

ADD_FEASIBILITY_CUT: 


# pragma region In the following, extended feasibility cut is added, also it came out aggregated version is better than disaggregated one

	for(int i=0;i<num_patients;i++)
	if ( copy_slp[i] != NULL ) {
		status = CPXfreeprob (env_slp, &copy_slp[i]);
		if ( status ) fprintf (stderr, "CPXfreeprob failed , error code %d.\n",status);}
	
	for(int i=0;i<num_patients;i++){		
	if(sol_state[i][numz+2]+sol_state[i][numz+3]>0){		
		
		
		for(int j=0; j<numz; j++){			
			if ((sol_state[i][j]==2) ||sol_state[i][j]==3)											
				cutval[j]=1;			
			else if(sol_state[i][j]==1) 			
				cutval[j]=sol_state[i][numz+2]+sol_state[i][numz+3];	 
			else
				cutval[j]=0;
			}	

		for (int j=numz;j<numz+num_patients;j++) 
			cutval[j]=0;
		cutrhs=(sol_state[i][numz+2]+sol_state[i][numz+3])*(sol_state[i][numz+1]);		
		
		status = CPXcutcallbackaddlocal (env, cbdata, wherefrom,numz+num_patients, cutrhs, 'L',index, cutval);								
		if ( status ) goto TERMINATE;			
		num_usercuts[0]++;					

		}		
	}	
	
	//for(int i=0;i<num_patients;i++)
	//delete[] v_x[i];
	//delete[] v_x;

	*useraction_p = CPX_CALLBACK_SET; 
#pragma endregion

TERMINATE:
	
	counter++;	  

	if(last_node[0]!=numiinf)
	last_node[0]=numiinf;				
		
#pragma region deallocating dynamic memories
	  delete []objval;
	  delete []cutval;	  
	  delete[]index;
	  delete []x;
	  delete []lb;
	  delete []ub;	  	
	  for(int i=0;i<num_patients;i++){
		delete[]v_prim[i];
		delete []y[i];
		delete []pi[i];
		delete[] sol_state[i];
		delete[]d_prim[i];
	  }
	  delete[]v_prim;
	  delete[]d_prim;
	  delete []y;
	  delete []pi;
	  delete[] sol_state;
	  for(int i=0;i<num_patients;i++)
		if ( copy_slp[i] != NULL ) {
			status = CPXfreeprob (env_slp, &copy_slp[i]);
			if ( status ) {
				fprintf (stderr, "CPXfreeprob failed , error code %d.\n",status);
		}
	}
#pragma endregion

	return 0;
	}

	int CPXPUBLIC 
	usercutcallback (CPXCENVptr env, void *cbdata, int wherefrom,void *cbhandle, int *useraction_p){	
	
	CUTINFOptr cutinfo = (CUTINFOptr) cbhandle;
	CPXCENVptr env_slp=cutinfo->env_slp;		
	CPXLPptr slp_t=cutinfo->slp_t;
	CPXLPptr* slp=cutinfo->slp;		
	CPXLPptr slp_vt=cutinfo->slp_vt;
	CPXLPptr* slp_v=cutinfo->slp_v;				
	double *vt=cutinfo->vt;	
	double **v=cutinfo->v;	
	double **d=cutinfo->d;	
	double **rhs_stop=cutinfo->rhs_stop;
	double **v_upd= cutinfo->v_upd; 
	double *vt_upd= cutinfo->vt_upd; 
	int obj_ind=cutinfo->obj_ind;
	bool *feas_h_status=cutinfo->feas_h_status;		
	int *last_node=&cutinfo->last_node;
	int *num_usercuts=&cutinfo->num_usercuts;

	double   *x= new double [numz+num_patients];
	double   *lb= new double [numz];
	double   *ub= new double [numz];
	int nz=0;
	double *cutval=new double[numz+num_patients];		
	double cutrhs;
	
	int *index= new int [numz+num_patients];
	for (int k=0;k<numz+num_patients;k++)
		index[k]=k;


	*useraction_p = CPX_CALLBACK_DEFAULT; 	

	int numiinf;
	status = CPXgetcallbacknodeinfo(env,cbdata,wherefrom,0,CPX_CALLBACK_INFO_NODE_SEQNUM,&numiinf);
	if ( status ) goto TERMINATE;

	status = CPXgetcallbacknodex (env, cbdata, wherefrom, x,0, numz+num_patients-1); 
	if ( status ) {
		fprintf(stderr, "Failed to get node solution.\n");
		goto TERMINATE;
	}	
	
	double gap;
	double best_int;
	double best_ub;
	status=CPXgetcallbackinfo(env, cbdata, wherefrom, CPX_CALLBACK_INFO_BEST_INTEGER, &best_int);
	status=CPXgetcallbackinfo(env, cbdata, wherefrom, CPX_CALLBACK_INFO_BEST_REMAINING, &best_ub);
	gap=(best_ub-best_int)/best_int;
	//////////////////////////////////// checking integer feasibility of solution of the current node, and aggressively generating cuts	
	for(int i=0; i<numz; i++){				
		if (x[i] > (1.0 - EPSILON)){
			x[i] = 1.0;	
			nz++;
		}
		else
			x[i] = 0.0;			
	}
	
	///////////////////////////////Dynamically updating V's/////////////////////////////
	int node_level;	
	status = CPXgetcallbacknodeinfo(env,cbdata,wherefrom,0,CPX_CALLBACK_INFO_NODE_DEPTH,&node_level);
	status = CPXgetcallbacknodelb (env, cbdata, wherefrom,lb, 0, numz-1);
	status = CPXgetcallbacknodeub (env, cbdata, wherefrom,ub, 0, numz-1);

	if( numiinf!=last_node[0]){		
		
		v_update(env_slp, slp_v, lb, ub, v, rhs_stop,v_upd);
		vt_update(env_slp, slp_vt, lb, ub, vt, rhs_stop,vt_upd);						
		ubt_update(env_slp,slp_t,lb,ub,v_upd,vt_upd,rhs_stop,&obj_ind);
		ub_update(env_slp,slp, lb, ub, v_upd, rhs_stop,&obj_ind);

		for(int i=0;i<num_patients;i++){				
				
		for(int j=0; j<numz+num_patients; j++)
				cutval[j]=0;		

		cutval[numz+i]=1;
		cutrhs=v_upd[i][obj_ind];				
		status = CPXcutcallbackaddlocal (env, cbdata, wherefrom,numz+num_patients, cutrhs, 'L',index, cutval);												
		if ( status ) goto TERMINATE;	
		num_usercuts[0]++;
		}	

		for(int j=0; j<numz; j++)
				cutval[j]=0;	
		for(int j=0; j<num_patients; j++)
		cutval[numz+j]=1;
		cutrhs=vt_upd[obj_ind];						

		status = CPXcutcallbackaddlocal (env, cbdata, wherefrom,numz+num_patients, cutrhs, 'L',index, cutval);												
		if(status) goto TERMINATE;
		num_usercuts[0]++;

		*useraction_p = CPX_CALLBACK_SET; 
	}

TERMINATE:
	
	if(last_node[0]!=numiinf)
		last_node[0]=numiinf;				
		
	counter++;	
	  	 	  
	  delete []cutval;	  
	  delete[]index;
	  delete []x;
	  delete []lb;
	  delete []ub;	  	
	
	return 0;
	}

		
	int CPXPUBLIC 
	myincumbentcheck (CPXCENVptr env,void*cbdata,int wherefrom,void*cbhandle,double objval,double*x,int *isfeas_p,int *useraction_p)
	{

	CUTINFOptr cutinfo = (CUTINFOptr) cbhandle;
	CPXCENVptr env_slp=cutinfo->env_slp;	
	CPXLPptr* slp=cutinfo->slp;		
	double **v=cutinfo->v;	
	double **rhs_stop=cutinfo->rhs_stop;


	double **y=new double *[ num_patients];	
	double *objval_prime=new double [num_patients];
	for(int k=0; k<num_patients; k++){
	y[k]= new double [numz];	
	}

	CPXLPptr copy_slp[num_patients];
	for (int i=0;i<num_patients;i++)
	{
		copy_slp[i] = CPXcloneprob (env_slp, slp[i], &status); // this creates a copy of slp problem 
		if ( status ) {
		fprintf (stderr, "CPXcopylp failed to creates a copy of slp problem in node.\n");
		goto TERMINATE;
		}		
	}


	if(feasibility_check(env_slp,slp,copy_slp,x,y,rhs_stop)==false){			
		*isfeas_p=0;
		*useraction_p=CPX_CALLBACK_SET;
		goto TERMINATE;
						
	}

	for(int k=0; k<num_patients; k++)
	{	
			status = CPXgetobjval (env_slp, copy_slp[k], &objval_prime[k]);
			if ( status ) { 
			fprintf (stderr, "Failed to obtain objval_prime.\n");
			goto TERMINATE;
			}
			if(x[numz+k]>objval_prime[k]+EPSILON)
			{
				*isfeas_p=0;
				*useraction_p=CPX_CALLBACK_SET;
				goto TERMINATE;

			}	
	}	 

	*isfeas_p=1;
	*useraction_p=CPX_CALLBACK_SET;


TERMINATE:

	for(int i=0;i<num_patients;i++)
		if ( copy_slp[i] != NULL ) {
			status = CPXfreeprob (env_slp, &copy_slp[i]);
			if ( status ) {
				fprintf (stderr, "CPXfreeprob failed , error code %d.\n",
						status);
			}
		delete []y[i];
	}
	
	delete []y;	  
	delete []objval_prime;
	
	return 0;	
	}


	int CPXPUBLIC mybranchfunc(CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle, int type, int sos,
	int nodecnt, int bdcnt, const int *nodebeg, const int *indices, const char *lu, const double *bd, const double *nodeest, int *useraction_p) { //this branch function may be used with CPLEX 12.6
	BRANCHINFOptr branchinfo = (BRANCHINFOptr)cbhandle;
	CPXENVptr env_slp=branchinfo->env_slp;	
	CPXLPptr* slp=branchinfo->slp;	
	double **rhs_stop=branchinfo->rhs_stop;	
		
	int obj_ind=branchinfo->obj_ind;
	
	

	double   *lb= new double [numz];
	double   *ub= new double [numz];
	double   *x= new double [numz];
	double   **y= new double* [num_patients];
	for(int k=0; k<num_patients; k++)
	y[k]= new double [numz];
	bool   *fixed_var= new bool [numz];

	*useraction_p=CPX_CALLBACK_DEFAULT;

	int integral_feas;
	status = CPXgetcallbacknodeinfo(env,cbdata,wherefrom,0,CPX_CALLBACK_INFO_NODE_NIINF,&integral_feas);

	status = CPXgetcallbacknodelb (env, cbdata, wherefrom,lb, 0, numz-1);
	status = CPXgetcallbacknodeub (env, cbdata, wherefrom,ub, 0, numz-1);
	for (int i=0;i<numz;i++){
		if(ub[i]-lb[i] <=EPSILON){
			x[i]=ub[i];	
			fixed_var[i]=true;//fixed variabels
		}
		else{
			x[i]=0;
			fixed_var[i]=false;//non-fixed variabels
		}
	}



	for (int i=0;i<numz;i++){
		 if(x[i]<=EPSILON)
			x[i]=0;
		else if(x[i]>=1-EPSILON)
			x[i]=1;
		else{
			*useraction_p=CPX_CALLBACK_SET; //this prunes the current node if there is a non-integer fixed variable
			goto TERMINATE;
		}
	}
	
	CPXLPptr copy_slp[num_patients];
	for (int i=0;i<num_patients;i++){
		copy_slp[i] = CPXcloneprob (env_slp, slp[i], &status); 
		if ( status ) goto TERMINATE;		
	}

	if(feasibility_check(env_slp,slp,copy_slp,x,y,rhs_stop)==false){									
		*useraction_p=CPX_CALLBACK_SET;
		goto TERMINATE;
	}
		
	
	int nez=0;
	int nz_ind=0;

	bool f_fix;
	for (int j=0;j<numz;j++){
		f_fix=false;
		int i=0;
		while(i<num_patients && f_fix==false){
			if(y[i][j]>rhs_stop[i][j]+EPSILON) f_fix=true;
			i++;
		}
		if(fixed_var[j]==false&& f_fix==true) nez++;
	}

	int      seqnum1;
	double   objval;

	status = CPXgetcallbacknodeobjval (env, cbdata, wherefrom,&objval);
	if ( status ) goto TERMINATE;	

	if(nez==0){		
		goto LOCAL_SEARCH;
	}		
	char *varlu= new char [nez];
	double  *varbd= new double [nez];
	int  *index= new int [nez];

	
	for (int j=0;j<numz;j++){
		f_fix=false;
		int i=0;
		while(i<num_patients && f_fix==false){
			if(y[i][j]>rhs_stop[i][j]+EPSILON) f_fix=true;
			i++;
		}
		if(fixed_var[j]==false&& f_fix==true){
		varlu[nz_ind]='U';
		varbd[nz_ind]=0;
		index[nz_ind]=j;
		nz_ind++;
		}	
	}
	
	status = CPXbranchcallbackbranchbds(env, cbdata, wherefrom, nez, index, varlu, varbd, objval, NULL, &seqnum1);
	 if ( status )  goto TERMINATE;	 

	delete[]varlu;	
	delete[]varbd; 
	delete[]index;
	 *useraction_p=CPX_CALLBACK_SET;
	 goto TERMINATE;	 
	///////////////////////////////////adding dual reduction based on my optimality cut//////////////////
LOCAL_SEARCH:
	 int	bestj=-1;
	 int s1_obj=obj_ind/num_states;
	 int s2_obj=obj_ind-(s1_obj*num_states);


	 int *feas= new int[numz];
	 status = CPXgetcallbacknodeintfeas(env, cbdata, wherefrom,feas, 0, numz-1);
	 if ( status )  goto TERMINATE;	 
	 for(int j=0;j<numz;j++){
		 if(feas[j]==CPX_INTEGER_INFEASIBLE){
			int s1_j=j/num_states;
			int s2_j=j-(s1_j*num_states);
			if(min(s1_obj-s1_j,s2_obj-s2_j)>=lb_local_search&&max(s1_obj-s1_j,s2_obj-s2_j)<=ub_local_search){
				bestj=j;				
				break;
			}
		 }			
	 }

	 delete[] feas;
	// If there weren't any eligible variables, take default branch  
   if ( bestj < 0 ) {
      goto TERMINATE;
   }

   	 char *varlu_local= new char [1];
	 double  *varbd_local= new double [1];
	 int    seqnum2;
	 int    seqnum3;

	    // Up node  
	varlu_local[0] = 'L';
	varbd_local[0] =  1;
	status = CPXbranchcallbackbranchbds(env, cbdata, wherefrom, 1, &bestj, varlu_local, varbd_local, objval, NULL, &seqnum2);
	if ( status )  goto TERMINATE;
 
	/// Down node 
 
	varlu_local[0] = 'U';
	varbd_local[0] =0;
	status = CPXbranchcallbackbranchbds(env, cbdata, wherefrom, 1, &bestj, varlu_local, varbd_local, objval, NULL, &seqnum3);
	if ( status )  goto TERMINATE;
	
	delete[]varlu_local;
	delete[]varbd_local;
	*useraction_p = CPX_CALLBACK_SET;

TERMINATE:

	for(int i=0;i<num_patients;i++)
	if ( copy_slp[i] != NULL ) {
		status = CPXfreeprob (env_slp, &copy_slp[i]);
		if ( status ) {
			fprintf (stderr, "CPXfreeprob failed , error code %d.\n",
					status);
		}
	}
	
	for(int k=0; k<num_patients; k++)
	delete [] y[k];
	delete[]lb;
	delete[]ub;
	delete[]x;
	delete[]y;
	delete[]fixed_var;

	 return 0;
	 }


	
	 /*
	int CPXPUBLIC mybranchfunc (CPXCENVptr env,void *cbdata,int wherefrom,void *cbhandle,int type,int sos, //This is used with CPLEX 12.4
	int nodecnt,int bdcnt,const double *nodeest,const int *nodebeg,const int *indices,const char *lu,const int *bd,int *useraction_p){ 
	BRANCHINFOptr branchinfo = (BRANCHINFOptr)cbhandle;
	CPXENVptr env_slp=branchinfo->env_slp;	
	CPXLPptr* slp=branchinfo->slp;	
	double **rhs_stop=branchinfo->rhs_stop;	
		
	int obj_ind=branchinfo->obj_ind;
	
	

	double   *lb= new double [numz];
	double   *ub= new double [numz];
	double   *x= new double [numz];
	double   **y= new double* [num_patients];
	for(int k=0; k<num_patients; k++)
	y[k]= new double [numz];
	bool   *fixed_var= new bool [numz];

	*useraction_p=CPX_CALLBACK_DEFAULT;

	int integral_feas;
	status = CPXgetcallbacknodeinfo(env,cbdata,wherefrom,0,CPX_CALLBACK_INFO_NODE_NIINF,&integral_feas);

	status = CPXgetcallbacknodelb (env, cbdata, wherefrom,lb, 0, numz-1);
	status = CPXgetcallbacknodeub (env, cbdata, wherefrom,ub, 0, numz-1);
	for (int i=0;i<numz;i++){
		if(ub[i]-lb[i] <=EPSILON){
			x[i]=ub[i];	
			fixed_var[i]=true;//fixed variabels
		}
		else{
			x[i]=0;
			fixed_var[i]=false;//non-fixed variabels
		}
	}



	for (int i=0;i<numz;i++){
		 if(x[i]<=EPSILON)
			x[i]=0;
		else if(x[i]>=1-EPSILON)
			x[i]=1;
		else{
			*useraction_p=CPX_CALLBACK_SET; //this prunes the current node if there is a non-integer fixed variable
			goto TERMINATE;
		}
	}
	
	CPXLPptr copy_slp[num_patients];
	for (int i=0;i<num_patients;i++){
		copy_slp[i] = CPXcloneprob (env_slp, slp[i], &status); 
		if ( status ) goto TERMINATE;		
	}

	if(feasibility_check(env_slp,slp,copy_slp,x,y,rhs_stop)==false){									
		*useraction_p=CPX_CALLBACK_SET;
		goto TERMINATE;
	}
		
	
	int nez=0;
	int nz_ind=0;

	bool f_fix;
	for (int j=0;j<numz;j++){
		f_fix=false;
		int i=0;
		while(i<num_patients && f_fix==false){
			if(y[i][j]>rhs_stop[i][j]+EPSILON) f_fix=true;
			i++;
		}
		if(fixed_var[j]==false&& f_fix==true) nez++;
	}

	int      seqnum1;
	double   objval;

	status = CPXgetcallbacknodeobjval (env, cbdata, wherefrom,&objval);
	if ( status ) goto TERMINATE;	

	if(nez==0){		
		goto LOCAL_SEARCH;
	}		
	char *varlu= new char [nez];
	int  *varbd= new int [nez];
	int  *index= new int [nez];

	
	for (int j=0;j<numz;j++){
		f_fix=false;
		int i=0;
		while(i<num_patients && f_fix==false){
			if(y[i][j]>rhs_stop[i][j]+EPSILON) f_fix=true;
			i++;
		}
		if(fixed_var[j]==false&& f_fix==true){
		varlu[nz_ind]='U';
		varbd[nz_ind]=0;
		index[nz_ind]=j;
		nz_ind++;
		}	
	}
	
	status = CPXbranchcallbackbranchbds (env, cbdata, wherefrom,objval, nez, index, varlu, varbd,NULL, &seqnum1);
	if ( status )  goto TERMINATE;	 

	delete[]varlu;	
	delete[]varbd; 
	delete[]index;
	 *useraction_p=CPX_CALLBACK_SET;
	 goto TERMINATE;	 
	///////////////////////////////////adding dual reduction based on my optimality cut//////////////////
LOCAL_SEARCH:
	 int	bestj=-1;
	 int s1_obj=obj_ind/num_states;
	 int s2_obj=obj_ind-(s1_obj*num_states);


	 int *feas= new int[numz];
	 status = CPXgetcallbacknodeintfeas(env, cbdata, wherefrom,feas, 0, numz-1);
	 if ( status )  goto TERMINATE;	 
	 for(int j=0;j<numz;j++){
		 if(feas[j]==CPX_INTEGER_INFEASIBLE){
			int s1_j=j/num_states;
			int s2_j=j-(s1_j*num_states);
			if(min(s1_obj-s1_j,s2_obj-s2_j)>=lb_local_search&&max(s1_obj-s1_j,s2_obj-s2_j)<=ub_local_search){
				bestj=j;				
				break;
			}
		 }			
	 }

	 delete[] feas;
	// If there weren't any eligible variables, take default branch  
   if ( bestj < 0 ) {
      goto TERMINATE;
   }

   	 char *varlu_local= new char [1];
	 int  *varbd_local= new int [1];
	 int    seqnum2;
	 int    seqnum3;

	    // Up node  
	varlu_local[0] = 'L';
	varbd_local[0] =  1;
	status = CPXbranchcallbackbranchbds (env, cbdata, wherefrom,objval, 1, &bestj, varlu_local, varbd_local,NULL, &seqnum2);
	if ( status )  goto TERMINATE;
 
	/// Down node 
 
	varlu_local[0] = 'U';
	varbd_local[0] =0;
	status = CPXbranchcallbackbranchbds (env, cbdata, wherefrom,objval, 1, &bestj, varlu_local, varbd_local,NULL, &seqnum3);
	if ( status )  goto TERMINATE;
	
	delete[]varlu_local;
	delete[]varbd_local;
	*useraction_p = CPX_CALLBACK_SET;

TERMINATE:

	for(int i=0;i<num_patients;i++)
	if ( copy_slp[i] != NULL ) {
		status = CPXfreeprob (env_slp, &copy_slp[i]);
		if ( status ) {
			fprintf (stderr, "CPXfreeprob failed , error code %d.\n",
					status);
		}
	}
	
	for(int k=0; k<num_patients; k++)
	delete [] y[k];
	delete[]lb;
	delete[]ub;
	delete[]x;
	delete[]y;
	delete[]fixed_var;

	 return 0;
	 }
*/
	int CPXPUBLIC myheuristic (CPXCENVptr env,void *cbdata,int wherefrom, 
	void *cbhandle,double *objval_p,double *x,int *checkfeas_p,int *useraction_p){
	BRANCHINFOptr branchinfo = (BRANCHINFOptr) cbhandle;
	CPXENVptr env_slp=branchinfo->env_slp;	
	CPXLPptr* slp=branchinfo->slp;	
	double **rhs_stop=branchinfo->rhs_stop;
	double *last_heur_solu=branchinfo->last_heur_solu;
	double *feas_h=branchinfo->feas_h;
	
	bool *feas_h_status=branchinfo->feas_h_status;

	double **y=new double *[num_patients];	
	double *objval=new double [num_patients];
	for(int k=0; k<num_patients; k++)
	y[k]= new double [numz];

	CPXLPptr copy_slp[num_patients];
	for (int i=0;i<num_patients;i++)
	{
		copy_slp[i] = CPXcloneprob (env_slp, slp[i], &status); // this creates a copy of slp problem 
		if ( status ) {
		fprintf (stderr, "CPXcopylp failed to creates a copy of slp problem in node.\n");
		goto TERMINATE;
		}		
	}

	*useraction_p = CPX_CALLBACK_DEFAULT;

	int *feas= new int[numz];
	 status = CPXgetcallbacknodeintfeas(env, cbdata, wherefrom,feas, 0, numz-1);
	 if ( status )  goto TERMINATE;	 

	 for(int j=0;j<numz;j++){
		 if(feas[j]==CPX_INTEGER_INFEASIBLE){
			 x[j]=0.0;
		 }			
	 }
	
	 // to do not double check the same solution for next iteration
	 bool last_heuristic=true;
	 for(int j=0;j<numz;j++){
		 if(last_heur_solu[j]!=x[j])
			 last_heuristic=false;
	 }
	 if (last_heuristic==true) goto TERMINATE;

	for(int j=0;j<numz;j++)
	last_heur_solu[j]=x[j];
	
	if(feas_h_status[0]==true){
	for(int i=0;i<numz+num_patients;i++)
		x[i]=feas_h[i]; 
	}
	else 
		goto TERMINATE;

	feas_h_status[0]=false;
	

	double obj_modified=0;
	for(int j=0;j<num_patients;j++)
	obj_modified+=x[numz+j];
	*objval_p=obj_modified;

	*checkfeas_p=0;//Tells CPLEX does not need to check the solution for integer feasibility
	*useraction_p = CPX_CALLBACK_SET;
TERMINATE:
	
	for(int i=0;i<num_patients;i++)
	if ( copy_slp[i] != NULL ) {
		status = CPXfreeprob (env_slp, &copy_slp[i]);
		if ( status ) {
			fprintf (stderr, "CPXfreeprob failed , error code %d.\n",status);
		}
	}
	for(int k=0; k<num_patients; k++)
	delete[]y[k];

	delete[]y;
	delete[]objval;
	delete[]feas;

	return 0;
	}



	bool feasibility(CPXCENVptr env_slp, CPXLPptr* slp, CPXLPptr* copy_slp, double *x, double **y,int **sol_state,double **rhs_stop) {
	bool feasibility=false;
	int *cstat= new int [numz];			
	int *rstat= new int [numz];
	for(int j=0; j<numz; j++){
		cstat[j]=CPX_BASIC; rstat[j]=CPX_AT_LOWER;}
	for(int k=0; k<num_patients; k++)
	for(int j=0; j<numz+4; j++)
		sol_state[k][j]=0;

	for(int k=0; k<num_patients; k++){			
		for(int i=0; i<numz; i++)	
		{
				if (x[i]==1){						
							double value=rhs_stop[k][i];
							for(int j=0; j<numz; j++)
							status = CPXchgcoef (env_slp, copy_slp[k], i, j, 0);				
							status = CPXchgcoef (env_slp, copy_slp[k], i, i, 1);
							status = CPXchgrhs (env_slp, copy_slp[k], 1, &i, &value);		
				}
			}
			status = CPXcopybase (env_slp, copy_slp[k], cstat, rstat);
			status = CPXlpopt (env_slp, copy_slp[k]);
			if ( status ) goto TERMINATE;

			int lpstat = CPXgetstat (env_slp, copy_slp[k]);
			if (lpstat==CPX_STAT_OPTIMAL){																								

				status = CPXgetx (env_slp, copy_slp[k], y[k], 0, numz-1);
				if ( status ) goto TERMINATE;
				double infeasout;
				for(int j=0; j<numz; j++)												  
				{											
					if ( status ) goto TERMINATE;
					if (x[j]==0) 
					{																					
						if (y[k][j]>rhs_stop[k][j]){							
							sol_state[k][j]=2;
							sol_state[k][numz+2]++;							
						}
						else{
							sol_state[k][j]=0;							
							sol_state[k][numz]++;							
						}
							
					}
					else if(x[j]==1)
					{				
						status = CPXgetrowinfeas (env_slp, slp[k], y[k], &infeasout, j, j);
						if (y[k][j] - infeasout>rhs_stop[k][j]){
							sol_state[k][j]=3;													
							sol_state[k][numz+3]++;											
						}
						else{
							sol_state[k][j]=1;								
							sol_state[k][numz+1]++;							
						}
					}

				}	
			}
			else{
			cout<< "There should be some problem in the input data since the subproblems are not either feasible or finite\n";
			error_counter++;						
			goto TERMINATE;
			}
		}

	TERMINATE:
	delete[] cstat;
	delete[] rstat;

	int sum_sol_states=0;
	for(int i=0;i<num_patients;i++)
		sum_sol_states+=sol_state[i][numz+3];

	if(sum_sol_states==0)  feasibility=true;
		return feasibility;
	}

	
	bool feasibility_strong(CPXCENVptr env_slp, CPXLPptr* slp, CPXLPptr* slp_v, double **y, double *lb, double *x,int **sol_state,double **rhs_stop) {
		bool feasibility=false;
		int *cstat= new int [numz];			
		int *rstat= new int [2*numz];
		int *delstat=new int [2*numz];	
		vector <int> myrow;


		for(int k=0; k<num_patients; k++)
		for(int j=0; j<numz+4; j++)
			sol_state[k][j]=0;

		CPXLPptr copy_slp_v[num_patients];
		for (int i=0;i<num_patients;i++)
		{
			copy_slp_v[i] = CPXcloneprob (env_slp, slp_v[i], &status); 
			if ( status ) {
			fprintf (stderr, "CPXcopylp failed to creates a copy of slp problem in node.\n");
			goto TERMINATE;
			}		
		}
		for(int k=0; k<num_patients; k++){

			for(int j=0; j<numz; j++){
				cstat[j]=CPX_BASIC;
				rstat[j]=CPX_BASIC;
				rstat[numz+j]=CPX_BASIC;
				delstat[j]=0;delstat[numz+j]=0;
			}
			for(int j=0; j<numz; j++){
				if(x[j]==0){			
					delstat[numz+j]=1;
					rstat[j]=CPX_AT_LOWER;
					rstat[numz+j]=2;
				}		

				if(lb[j]==1){ 
					delstat[j]=1;
					rstat[numz+j]=CPX_AT_LOWER;
					rstat[j]=2;
				}
			}

			status = CPXdelsetrows (env_slp, copy_slp_v[k], delstat);
			for(int j=0;j<2*numz;j++){
			if(rstat[j]!=2)
			myrow.push_back(rstat[j]);}
			int *rstat_n= new int[myrow.size()];
			for(int j=0;j<myrow.size();j++)
				rstat_n[j]=myrow[j];

			status = CPXcopybase (env_slp, copy_slp_v[k], cstat, rstat);
			delete[] rstat_n;
			status = CPXlpopt (env_slp, copy_slp_v[k]); 
			if ( status ) goto TERMINATE;

				int lpstat = CPXgetstat (env_slp, copy_slp_v[k]);
				if (lpstat==CPX_STAT_OPTIMAL){																								

					status = CPXgetx (env_slp, copy_slp_v[k], y[k], 0, numz-1);
					if ( status ) goto TERMINATE;
					double infeasout;
					for(int j=0; j<numz; j++)												  
					{																	
						if (x[j]<1E-6) 
						{																					
							if (y[k][j]>rhs_stop[k][j]+EPSILON){							
								sol_state[k][j]=2;
								sol_state[k][numz+2]++;							
							}
							else{
								sol_state[k][j]=0;							
								sol_state[k][numz]++;							
							}
							
						}
						else if(x[j]>1-1E-6)
						{				
							status = CPXgetrowinfeas (env_slp, slp[k], y[k], &infeasout, j, j);
							if (y[k][j] - infeasout>rhs_stop[k][j]+EPSILON){
								sol_state[k][j]=3;													
								sol_state[k][numz+3]++;											
							}
							else{
								sol_state[k][j]=1;								
								sol_state[k][numz+1]++;							
							}
						}

					}	
				}
				else{
				cout<< "There should be some problem in the input data since the subproblems are not either feasible or finite\n";
				error_counter++;						
				goto TERMINATE;
				}
			}

	TERMINATE:

		for(int i=0;i<num_patients;i++)
		if ( copy_slp_v[i] != NULL ) {
			status = CPXfreeprob (env_slp, &copy_slp_v[i]);
			if ( status ) {
				fprintf (stderr, "CPXfreeprob failed , error code %d.\n",status);
			}			
		}	
		delete[] cstat;
		delete[] rstat;
		
		int sum_sol_states=0;
		for(int i=0;i<num_patients;i++)
			sum_sol_states+=sol_state[i][numz+3];

		if(sum_sol_states==0)  feasibility=true;
		return feasibility;
	}

	bool feasibility_check(CPXCENVptr env_slp,CPXLPptr* slp, CPXLPptr* copy_slp, double *x,double **solu,double **rhs_stop){
		bool feasibility;
		int *cstat= new int [numz];			
		int *rstat= new int [numz];
		for(int j=0; j<numz; j++){
		cstat[j]=CPX_BASIC; rstat[j]=CPX_AT_LOWER;}
		for(int k=0; k<num_patients; k++)
		for(int i=0; i<numz; i++)	
		{
				if (x[i]==1){						
							double value=rhs_stop[k][i];
							for(int j=0; j<numz; j++)
							status = CPXchgcoef (env_slp, copy_slp[k], i, j, 0);				
							status = CPXchgcoef (env_slp, copy_slp[k], i, i, 1);
							status = CPXchgrhs (env_slp, copy_slp[k], 1, &i, &value);		
				}
			}

		for(int k=0; k<num_patients; k++)
		{
			status = CPXcopybase (env_slp, copy_slp[k], cstat, rstat);
			status = CPXlpopt (env_slp, copy_slp[k]); 
			if ( status ) {
			fprintf (stderr, "Failed to optimize LP.\n");
			goto TERMINATE;
				}

			int lpstat = CPXgetstat (env_slp, copy_slp[k]);
			if (lpstat==CPX_STAT_OPTIMAL){																								

				status = CPXgetx (env_slp, copy_slp[k], solu[k], 0, numz-1);
				if ( status ) {
				fprintf (stderr, "Failed to get optimal solution of  the subproblem.\n");
				goto TERMINATE;
					}

				for(int j=0; j<numz; j++)												  
				{						
					if (x[j]==1) 
					{	
						double infeasout;
						status = CPXgetrowinfeas (env_slp, slp[k], solu[k], &infeasout, j, j);
						if ( status ) {
						fprintf (stderr, "Failed to get row infeasiblity.\n");
						goto TERMINATE;
							}							
						
						if (solu[k][j] - infeasout>rhs_stop[k][j]+EPSILON){
							feasibility=false;
							goto TERMINATE;
						}
					}

				}	
			}
			else{
			cout<< "There should be some problem in the input data since the subproblems are not either feasible or finite\n";
			error_counter++;			
			cout<<"lpstat is "<<lpstat<<endl;
			goto TERMINATE;
			}
		}
		feasibility=true;
TERMINATE:
		delete[] cstat;
		delete[] rstat;
	return feasibility;
	}


	int v_update(CPXCENVptr env_slp, CPXLPptr* slp_v, double *lb, double *ub, double **v, double **rhs_stop, double **v_upd){
	int *cstat= new int [numz];			
	int **rstat= new int *[num_patients];
	int **delstat= new int *[num_patients];
	vector <int> *myvector= new vector<int>[num_patients];	
	for (int i = 0; i < num_patients; i++){
		delstat[i]=new int [2*numz];	
		rstat[i]= new int [2*numz];		
	}
	for(int j=0; j<numz; j++)
		cstat[j]=CPX_BASIC;

	CPXLPptr copy_slp_v[num_patients];
	for (int i=0;i<num_patients;i++)
	{
		copy_slp_v[i] = CPXcloneprob (env_slp, slp_v[i], &status); 
		if ( status ) {
		fprintf (stderr, "CPXcopylp failed to creates a copy of slp problem in node.\n");
		goto TERMINATE;
		}		
	}

	for(int i=0;i<num_patients;i++)
	for(int j=0;j<2*numz;j++){
		rstat[i][j]=1;
		delstat[i][j]=0;
	}
	
	for(int k=0;k<num_patients;k++)
	for(int j=0;j<numz;j++){
		if(lb[j]>=1-EPSILON ||(v[k][j]==rhs_stop[k][j]&& ub[j]!=0)){			
				delstat[k][j]=1; rstat[k][j]=2; rstat[k][numz+j]=0;

		}		
		if(ub[j]<=EPSILON){	
				delstat[k][j+numz]=1; rstat[k][j]=0;rstat[k][numz+j]=2;
		}			
	}
		
	
	for (int k=0;k<num_patients;k++){
		status = CPXdelsetrows (env_slp, copy_slp_v[k], delstat[k]);
		if ( status ) 	goto TERMINATE;		
	}

	for(int i=0;i<num_patients;i++)
	for(int j=0;j<2*numz;j++){
		if(rstat[i][j]!=2)
		myvector[i].push_back(rstat[i][j]);}

	for (int k=0;k<num_patients;k++){
			int *rstat_n= new int[myvector[k].size()];
			for(int j=0;j<myvector[k].size();j++)
				rstat_n[j]=myvector[k][j];
			status = CPXcopybase (env_slp, copy_slp_v[k], cstat, rstat_n);			
			status = CPXlpopt (env_slp, copy_slp_v[k]);
			status = CPXgetx (env_slp, copy_slp_v[k], v_upd[k], 0, numz-1);
			delete[] rstat_n;
	}
	TERMINATE:

	for(int i=0;i<num_patients;i++)
	if ( copy_slp_v[i] != NULL ) {
		status = CPXfreeprob (env_slp, &copy_slp_v[i]);
		if ( status ) {
			fprintf (stderr, "CPXfreeprob failed , error code %d.\n",status);
		}
	}		
	for (int i = 0; i < num_patients; i++){
	delete[] delstat[i]; rstat[i];
	}
	delete[] delstat;
	delete[] cstat;
	delete[] rstat;
	delete[] myvector;
	return 0;
	}

int vt_update(CPXCENVptr env_slp, CPXLPptr slp_vt, double *lb, double *ub, double *vt,double **rhs_stop, double *vt_upd){
	int *cstat= new int [numz];			
	int *rstat= new int [2*numz];
	int *delstat= new int [2*numz];
	vector <int> myvector;

	for(int j=0; j<numz; j++)
		cstat[j]=CPX_BASIC;	

	CPXLPptr copy_slp_vt = CPXcloneprob (env_slp, slp_vt, &status); // this creates a copy of slp problem 
	if ( status ) goto TERMINATE;

	for(int j=0;j<2*numz;j++){
		rstat[j]=1;
		delstat[j]=0; 
	}
		
	double sum_vt;
	for(int j=0;j<numz;j++){
		sum_vt=0;
		for(int i=0;i<num_patients;i++)
			sum_vt+=rhs_stop[i][j];

		if(lb[j]>=1-EPSILON ||(vt[j]==sum_vt&& ub[j]!=0)){			
				delstat[j]=1; rstat[j]=2; rstat[numz+j]=0;

		}		
		else{			
				delstat[j]=0; 
		}

		if(ub[j]<=EPSILON){	
				delstat[j+numz]=1; rstat[j]=0;rstat[numz+j]=2;
		}			
		else			
			delstat[j+numz]=0;
		
	}
	
	status = CPXdelsetrows (env_slp, copy_slp_vt, delstat);
	if ( status ) 	goto TERMINATE;		
	
	for(int j=0;j<2*numz;j++){
		if(rstat[j]!=2)
		myvector.push_back(rstat[j]);}

	
	int *rstat_n= new int[myvector.size()];
	for(int j=0;j<myvector.size();j++)
		rstat_n[j]=myvector[j];
	status = CPXcopybase (env_slp, copy_slp_vt, cstat, rstat_n);			
	status = CPXlpopt (env_slp, copy_slp_vt);	
	status = CPXgetx (env_slp, copy_slp_vt, vt_upd, 0, numz-1);
	delete[] rstat_n;
	
	TERMINATE:

	if ( copy_slp_vt != NULL ) {
		status = CPXfreeprob (env_slp, &copy_slp_vt);
		if ( status ) {
			fprintf (stderr, "CPXfreeprob failed , error code %d.\n",status);
		}
	}			
	delete[] delstat;
	delete[] cstat;
	delete[] rstat;
	
	return 0;
	}

		
	int ub_update(CPXCENVptr env_slp, CPXLPptr* slp, double *lb, double *ub, double **v, double **rhs_stop,int *obj_ind){
	char sense;			
	int *indices= new int [numz];	
	char *lu= new char [numz];
	double *bd= new double [numz];	
	int *cstat= new int [numz];	
	
	for(int j=0; j<numz; j++){
		indices[j]=j;
		cstat[j]=CPX_BASIC;
	}
	CPXLPptr copy_slp[num_patients];
	for (int i=0;i<num_patients;i++){
		copy_slp[i] = CPXcloneprob (env_slp, slp[i], &status); 
		if ( status ) {fprintf (stderr, "CPXcopylp failed to creates a copy of slp problem in node.\n");
		goto TERMINATE;}		
	}

	sense='G';
	for(int i=0;i<num_patients;i++){
	for(int j=0; j<numz; j++){
		if(ub[j]==1)  status = CPXchgsense (env_slp, copy_slp[i], 1, &j, &sense);
		if(lb[j]==1) {bd[j]=rhs_stop[i][j]; lu[j]='B';}
		else {bd[j]=v[i][j]; lu[j]='U';}
		}
		 status = CPXchgbds (env_slp, copy_slp[i], numz, indices, lu, bd);
	}
	
	for (int i=0;i<num_patients;i++){						
			status = CPXcopybase (env_slp, copy_slp[i], cstat, NULL);			
			status = CPXlpopt (env_slp, copy_slp[i]);			
			status = CPXgetobjval (env_slp, copy_slp[i], &v[i][*obj_ind]);			
	}
	TERMINATE:

	for(int i=0;i<num_patients;i++)
	if ( copy_slp[i] != NULL ) {
		status = CPXfreeprob (env_slp, &copy_slp[i]);
		if ( status ) {fprintf (stderr, "CPXfreeprob failed , error code %d.\n",status);
		}
	}		

	delete[] indices;
	delete[] cstat;
	delete[] lu;
	delete[] bd;

	return 0;
	}
	
	int ubt_update(CPXCENVptr env_slp, CPXLPptr slp_t, double *lb, double *ub, double **v_upd, double *vt_upd, double **rhs_stop, int *obj_ind){
	char sense;			
	int *indices= new int [numz];	
	char *lu= new char [numz];
	double *bd= new double [numz];	
	int *cstat= new int [num_patients*numz];	
	
	for(int j=0; j<num_patients*numz; j++)		
		cstat[j]=CPX_BASIC;
	
	CPXLPptr copy_slp_t = CPXcloneprob (env_slp, slp_t, &status); 
	if ( status ) {fprintf (stderr, "CPXcopylp failed to creates a copy of slp problem in node.\n");goto TERMINATE;}		
	
	sense='G';
	for(int i=0;i<num_patients;i++){
	for(int j=0; j<numz; j++){
		int k;
		if(ub[j]==1){
			k=numz*i+j;
			status = CPXchgsense (env_slp, copy_slp_t, 1, &k, &sense);
		}
		if(lb[j]==1) {bd[j]=rhs_stop[i][j]; lu[j]='B';}
		else {bd[j]=v_upd[i][j]; lu[j]='U';}
		}

		for(int j=0; j<numz; j++)
		indices[j]=numz*i+j;
		 status = CPXchgbds (env_slp, copy_slp_t, numz, indices, lu, bd);
	}
	
	for(int j=0; j<numz; j++)
		indices[j]=numz*num_patients+j;
	status = CPXchgrhs (env_slp, copy_slp_t, numz, indices, vt_upd);
						
	status = CPXcopybase (env_slp, copy_slp_t, cstat, NULL);			
	status = CPXlpopt (env_slp, copy_slp_t);			
	status = CPXgetobjval (env_slp, copy_slp_t, &vt_upd[*obj_ind]);			
	


	TERMINATE:
	
	if ( copy_slp_t != NULL ) 
		status = CPXfreeprob (env_slp, &copy_slp_t);
	if ( status ) {fprintf (stderr, "CPXfreeprob failed , error code %d.\n",status);}
	

	delete[] indices;
	delete[] cstat;
	delete[] lu;
	delete[] bd;

	return 0;
	}
	
