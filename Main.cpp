/* This code is implementation of branch and cut algorithm for finding the best equilibrium of kidney exchange instances with three players. */ 	

	/* Bring in the declarations for the string and character functions, 
	malloc, and fabs */
	#include <ilcplex/cplex.h>
	#include <ctype.h>
	#include <string.h>
	#include <math.h>
	#include <stdio.h>// standard input output to use printer and other device 
	#include <stdlib.h>//C Standard General Utilities Library. This header defines several general purpose functions, including dynamic memory management, random number generation, communication with the environment, integer arthmetics, searching, sorting and converting.
	#include <iostream> // 
	#include <fstream> // to read or write in a file
	#include <iomanip>
	#include <ctime>
	#include <time.h>
	#include <sstream>
	#include <vector> 
	#include "Cut Generator.h"
	#include "Make problems.h"	
	#include <iomanip>
	
	using namespace std; 
			///////////////// it is better than we have the following variables ordered as patient_1 <=patient_2	
	int numz;	
	int counter=0;
	int error_counter=0;
	int status;
	
	void branch_and_cut();
	void cplex();
	void negative_reward_branch_and_cut();
	void negative_reward_cplex();

	int main(){
	if(solution_method==BRANCH_AND_CUT && reward_type==REWARD_POSITIVE) branch_and_cut();
	else if(solution_method==BRANCH_AND_CUT && reward_type==REWARD_NEGATIVE) negative_reward_branch_and_cut();
	else if(solution_method==CPLEX && reward_type==REWARD_POSITIVE) cplex();
	else if(solution_method==CPLEX && reward_type==REWARD_NEGATIVE) negative_reward_cplex();

	return 0;
}




	void branch_and_cut(){

#pragma omp parallel num_threads(num_cores)
	{
	int instance_num;
	#pragma omp for private(instance_num)
	for( instance_num=instance_beg;instance_num<=instance_end;instance_num++){

#pragma region initialization

	const double base=num_states;
	numz=pow (base, num_patients);	

	int obj_ind;// This specifies the initial state of the game for consideration of objective value			
	int num_usercuts=0;
	bool *feas_h_status= new bool [1];
	feas_h_status[0]=false;		
	int last_node=-1;	


	double diff;
	clock_t start, finish;
	start = clock();
	/////////////////////////start reading data///////////////////////////////////////
	double ***tran_prob= new double ** [num_patients];
	for (int k = 0; k < num_patients; k++)
	{
		tran_prob[k]= new double *[num_states];
	for (int i = 0; i < num_states; i++)
		tran_prob[k][i] = new double [num_states];
	}

		
		//taking reward matrix
	double **reward = new double * [num_patients];
	for (int i = 0; i < num_patients; i++)
	reward[i] = new double [num_states];

			//taking stopping reward matrix
	double **stop_reward = new double * [num_patients];
	for (int i = 0; i < num_patients; i++)
	stop_reward[i] = new double [num_states];

	double **rhs_contin = new double * [num_patients];
	double **rhs_stop = new double * [num_patients];
	double *vt= new double [numz];
	double **v= new double * [num_patients]; // to store V values
	double **d= new double * [num_patients]; // to store d values
	double **v_upd= new double * [num_patients]; // to store updated V's	
	double *vt_upd= new double [numz]; // to store updated V's		
	
	double *last_heur_solu= new double [numz]; // to store solution of last heuristic in order to do not double check the same solution for next iteration
	double *feas_h= new double [numz+num_patients]; // feasibility cut obtained in heuristic
			
	double *y=new double [numz+num_patients];	
	for (int i = 0; i < num_patients; i++){
		rhs_contin[i] = new double [numz];
		rhs_stop[i]= new double [numz];
		v[i]= new double [numz];
		d[i]= new double [numz];		
		v_upd[i]= new double [numz];
	}
	
				
	if( num_patients==2){
		readdata2( tran_prob, reward, stop_reward, &obj_ind, &instance_num );

		for(int i=0;i<num_patients;i++)
		for(int j1=0;j1<num_states;j1++)
		for(int j2=0;j2<num_states;j2++){
		if (j1!=num_states-1&& j2!=num_states-1)
			rhs_stop[i][j1*num_states+j2]=stop_reward[i][i*j2+(1-i)*j1];
		else
			rhs_stop[i][j1*num_states+j2]=0;
		rhs_contin[i][j1*num_states+j2]=reward[i][i*j2+(1-i)*j1];		
	}	

	}
	if( num_patients==3){
		readdata3( tran_prob, reward, stop_reward, &obj_ind, &instance_num );

		for(int j1=0;j1<num_states;j1++)
		for(int j2=0;j2<num_states;j2++)
		for(int j3=0;j3<num_states;j3++){
			if (j1!=num_states-1&& j2!=num_states-1&& j3!=num_states-1)
				rhs_stop[0][j1*num_states*num_states+j2*num_states+j3]=stop_reward[0][j1];
			else
				rhs_stop[0][j1*num_states*num_states+j2*num_states+j3]=0;
			rhs_contin[0][j1*num_states*num_states+j2*num_states+j3]=reward[0][j1];		
		}	

		for(int j1=0;j1<num_states;j1++)
		for(int j2=0;j2<num_states;j2++)
		for(int j3=0;j3<num_states;j3++){
			if (j1!=num_states-1&& j2!=num_states-1&& j3!=num_states-1)
				rhs_stop[1][j1*num_states*num_states+j2*num_states+j3]=stop_reward[1][j2];
			else
				rhs_stop[1][j1*num_states*num_states+j2*num_states+j3]=0;
			rhs_contin[1][j1*num_states*num_states+j2*num_states+j3]=reward[1][j2];		
		}	

		for(int j1=0;j1<num_states;j1++)
		for(int j2=0;j2<num_states;j2++)
		for(int j3=0;j3<num_states;j3++){
			if (j1!=num_states-1&& j2!=num_states-1&& j3!=num_states-1)
				rhs_stop[2][j1*num_states*num_states+j2*num_states+j3]=stop_reward[2][j3];
			else
				rhs_stop[2][j1*num_states*num_states+j2*num_states+j3]=0;
			rhs_contin[2][j1*num_states*num_states+j2*num_states+j3]=reward[2][j3];		
	}	

	}

	/////////////////////////finsih reading data///////////////////////////////////////



	CPXENVptr env;
	CPXENVptr env_slp;
	CPXLPptr lp;
	CPXLPptr slp[num_patients];
	CPXLPptr slp_vt;
	CPXLPptr slp_t;
	CPXLPptr slp_v[num_patients];// for calculation of V dynamically	
	
	env = CPXopenCPLEX (&status);
	if ( env == NULL ) {
	char  errmsg[1024];
	cout<<"Could not open CPLEX environment.\n";
	CPXgeterrorstring (env, status, errmsg);
	fprintf (stderr, "%s", errmsg);
	goto TERMINATE;
	}	
	status = CPXsetintparam (env, CPX_PARAM_THREADS, 1);
	//status = CPXsetintparam(env, CPX_PARAM_NODEFILEIND, 2);// CPLEX can use node files and write on hard disk
    //status = CPXsetdblparam(env, CPX_PARAM_WORKMEM, 2048); // the number of memory units available for CPLEX

	env_slp = CPXopenCPLEX (&status);
	if ( env_slp == NULL ) {
		char  errmsg[1024];
		cerr<<"Could not open CPLEX environment for subproblem.\n";
		CPXgeterrorstring (env_slp, status, errmsg);
		fprintf (stderr, "%s", errmsg);
		goto TERMINATE;
		}
	lp= CPXcreateprob (env, &status, "master problem");
	if ( lp == NULL ) {
	fprintf (stderr, "Failed to create master problem.\n");
	goto TERMINATE;
	}
	for(int i=0;i<num_patients;i++){
		slp[i]= CPXcreateprob (env_slp, &status, "subproblem1");
		if ( slp[i] == NULL ) {
		fprintf (stderr, "Failed to create subproblem.\n");
		goto TERMINATE;
		}
	}
	slp_vt= CPXcreateprob (env_slp, &status, "slp_vt");// we create a problem 
	if ( slp_vt == NULL ) {
	fprintf (stderr, "Failed to create subproblem.\n");
	goto TERMINATE;
	}
	for(int i=0;i<num_patients;i++){
		slp_v[i]= CPXcreateprob (env_slp, &status, "subproblem1");
		if ( slp_v[i] == NULL ) {
		fprintf (stderr, "Failed to create subproblem.\n");
		goto TERMINATE;
		}
	}
	
	slp_t= CPXcreateprob (env_slp, &status, "subproblem2");
	if ( slp_t == NULL ) {
	fprintf (stderr, "Failed to create master problem.\n");
	goto TERMINATE;
	}

	
	
	double mip_setup1; double mip_setup2; double diff_mip_setup;
	CPXgettime( env,  &mip_setup1);   	

	if(num_patients==2){ 
		make_master2(env, env_slp, lp,slp,slp_v,tran_prob, rhs_contin, rhs_stop,v,d, &obj_ind);
		make_slpvt2(env_slp,slp_t,slp_vt,tran_prob,rhs_contin,rhs_stop,vt,d, &obj_ind);
	}
	if(num_patients==3){
		make_master3(env, env_slp, lp,slp,slp_v,tran_prob, rhs_contin, rhs_stop,v,d, &obj_ind);
		make_slpvt3(env_slp,slp_t,slp_vt,tran_prob,rhs_contin,rhs_stop,vt,d, &obj_ind);
	}

	

	CPXgettime( env,  &mip_setup2);
	diff_mip_setup=mip_setup2-mip_setup1;
	
	
 
#pragma endregion
	
#pragma region Set parameters

	
	/// Assure linear mappings between the presolved and originalmodels 
	status = CPXsetintparam (env, CPX_PARAM_PRELINEAR, CPX_OFF);	
	

	status = CPXsetintparam (env, CPX_PARAM_SCRIND, CPX_ON);			

	status = CPXsetintparam (env, CPX_PARAM_MIPSEARCH, CPX_MIPSEARCH_TRADITIONAL);
	
	// Let MIP callbacks work on the original model 
	
	status = CPXsetintparam (env, CPX_PARAM_MIPCBREDLP, CPX_OFF);
	
	status = CPXsetintparam (env, CPX_PARAM_MIPDISPLAY, 2);	


	status = CPXsetintparam (env,  CPX_PARAM_SUBALG, CPX_ALG_PRIMAL);
		
	
	status = CPXsetintparam (env, CPX_PARAM_VARSEL , CPX_VARSEL_STRONG);
	
	status = CPXsetintparam (env, CPX_PARAM_NODESEL , CPX_NODESEL_BESTBOUND);

	CPXsetdblparam (env, CPX_PARAM_BTTOL , 0.0001);


	status = CPXsetintparam (env,  CPX_PARAM_HEURFREQ , -1);//Controls the frequency of applying heuristic for MIP though	

	status = CPXsetintparam (env, CPX_PARAM_REDUCE , 0 ); // This is ctitical. This way it can't use dual reduction.
	
	status = CPXsetintparam (env_slp,  CPX_PARAM_LPMETHOD  , CPX_ALG_DUAL ); 
#pragma endregion	

#pragma region	initializing heuristic
	
		for(int j=0;j<numz;j++)
		feas_h[j]=0;

		bool heur_bool=true;
		for(int i=0;i<num_patients;i++)
			if(rhs_stop[i][obj_ind] < d[i][obj_ind])heur_bool=false;

		if(heur_bool==true){
			feas_h[obj_ind]=1;

			for(int i=0;i<num_patients;i++)
				feas_h[numz+i]=rhs_stop[i][obj_ind];
		}
		else{
			for(int i=0;i<num_patients;i++)
				feas_h[numz+i]=d[i][obj_ind];			
		}
	
#pragma endregion

#pragma region	warm start
	if(num_patients==2){
		double *x= new double [2*(numz+num_patients)];
		double **solu= new double *[num_patients];
		int *indices=new int[2*(numz+num_patients)];
		CPXLPptr copy_slp[num_patients];
		for(int i=0;i<num_patients;i++){
			copy_slp[i] = CPXcloneprob (env_slp, slp[i], &status); 
			solu[i]= new double[numz];
		}
		int *beg=new int [2]; beg[0]=0;	beg[1]=numz+num_patients;
		int *effort_level= new int [2]; effort_level[0]=1;effort_level[1]=1;
		for(int i=0;i<num_states;i++)
		for(int j=0;j<num_states;j++){
			indices[i*num_states+j]=i*num_states+j;
			if(i!=num_states-1 && j!= num_states-1 && rhs_stop[0][i*num_states+j]==v[0][i*num_states+j]&&rhs_stop[1][i*num_states+j]==v[1][i*num_states+j])
				x[i*num_states+j]=1;						
			else
				x[i*num_states+j]=0;			
		}


		feasibility_check( env_slp, slp, copy_slp, x,solu,rhs_stop);
		for(int j=0;j<num_patients;j++){
			x[numz+j]=solu[j][obj_ind];
			indices[numz+j]=numz+j;
		}
		for(int j=0;j<numz+num_patients;j++){
			x[numz+num_patients+j]=feas_h[j];
			indices[numz+num_patients+j]=j;
		}
		status = CPXaddmipstarts (env, lp, 2, 2*(numz+num_patients), beg, indices,x, effort_level, NULL);
		if ( status ) goto TERMINATE;
		for(int i=0;i<num_patients;i++)
		if ( slp[i] != NULL ) {
			status = CPXfreeprob (env_slp, &copy_slp[i]);
			if ( status ) {
				fprintf (stderr, "CPXfreeprob failed , error code %d.\n",status);
			}
			delete[] solu[i];
		}
		delete[] beg;
		delete[] effort_level;
		delete[] x;
		delete[] indices;
		delete[] solu;
	}

	if(num_patients==3){

		double *x= new double [2*(numz+num_patients)];
		double **solu= new double *[num_patients];
		int *indices=new int[2*(numz+num_patients)];
		CPXLPptr copy_slp[num_patients];
		for(int i=0;i<num_patients;i++){
			copy_slp[i] = CPXcloneprob (env_slp, slp[i], &status); 
			solu[i]= new double[numz];
		}
		int *beg=new int [2]; beg[0]=0;	beg[1]=numz+num_patients;
		int *effort_level= new int [2]; effort_level[0]=1;effort_level[1]=1;
		for(int i1=0;i1<num_states;i1++)
		for(int i2=0;i2<num_states;i2++)
		for(int i3=0;i3<num_states;i3++){
			indices[i1*num_states*num_states+i2*num_states+i3]=i1*num_states*num_states+i2*num_states+i3;
			if(i1!=num_states-1 && i2!= num_states-1 &&i3!= num_states-1 && rhs_stop[0][i1*num_states*num_states+i2*num_states+i3]==v[0][i1*num_states*num_states+i2*num_states+i3]&&rhs_stop[1][i1*num_states*num_states+i2*num_states+i3]==v[1][i1*num_states*num_states+i2*num_states+i3]&&rhs_stop[2][i1*num_states*num_states+i2*num_states+i3]==v[2][i1*num_states*num_states+i2*num_states+i3])
				x[i1*num_states*num_states+i2*num_states+i3]=1;						
			else
				x[i1*num_states*num_states+i2*num_states+i3]=0;			
		}


		feasibility_check( env_slp, slp, copy_slp, x,solu,rhs_stop);
		for(int j=0;j<num_patients;j++){
			x[numz+j]=solu[j][obj_ind];
			indices[numz+j]=numz+j;
		}
		for(int j=0;j<numz+num_patients;j++){
			x[numz+num_patients+j]=feas_h[j];
			indices[numz+num_patients+j]=j;
		}
		status = CPXaddmipstarts (env, lp, 2, 2*(numz+num_patients), beg, indices,x, effort_level, NULL);
		if ( status ) goto TERMINATE;
		for(int i=0;i<num_patients;i++)
		if ( slp[i] != NULL ) {
			status = CPXfreeprob (env_slp, &copy_slp[i]);
			if ( status ) {
				fprintf (stderr, "CPXfreeprob failed , error code %d.\n",status);
			}
			delete[] solu[i];
		}
		delete[] beg;
		delete[] effort_level;
		delete[] x;
		delete[] indices;
		delete[] solu;
	}
#pragma endregion

	
#pragma region construct initialization
	CUTINFO cutinfo;

	cutinfo.env_slp=env_slp;			
	cutinfo.slp_t=slp_t;
	cutinfo.slp=slp;
	cutinfo.slp_vt=slp_vt;
	cutinfo.slp_v=slp_v;
	cutinfo.rhs_contin=rhs_contin;
	cutinfo.rhs_stop=rhs_stop;
	
	cutinfo.v=v;
	cutinfo.vt=vt;
	cutinfo.v_upd=v_upd;
	cutinfo.vt_upd=vt_upd;
	cutinfo.d=d;
	
	cutinfo.feas_h=feas_h;
			
	cutinfo.obj_ind=obj_ind;
	cutinfo.num_usercuts=num_usercuts;
	cutinfo.feas_h_status=feas_h_status;
	cutinfo.last_node=last_node;

	BRANCHINFO branchinfo;
	branchinfo.env_slp=env_slp;
	branchinfo.slp=slp;
	branchinfo.rhs_stop=rhs_stop;
	branchinfo.last_heur_solu=last_heur_solu;
	branchinfo.feas_h=feas_h;
			
	branchinfo.obj_ind=obj_ind;
	branchinfo.feas_h_status=feas_h_status;
#pragma endregion

#pragma region	Set up to use MIP callback 
		
	status = CPXsetbranchcallbackfunc (env, mybranchfunc, &branchinfo); 	
   if ( status )  goto TERMINATE;
      
	status = CPXsetusercutcallbackfunc(env, usercutcallback, &cutinfo);      
  if ( status )  goto TERMINATE;

   status = CPXsetlazyconstraintcallbackfunc (env, mycutcallback, &cutinfo); 
   if ( status )  goto TERMINATE;   

   status = CPXsetincumbentcallbackfunc (env, myincumbentcheck,&cutinfo);
   if ( status ) goto TERMINATE;
   
   status = CPXsetheuristiccallbackfunc (env, myheuristic, &branchinfo);
   if ( status )  goto TERMINATE;
#pragma endregion 

#pragma region solving MIP and getting the solution

    status = CPXsetdblparam  (env,  CPX_PARAM_TILIM , time_limit);// CPLEX time limit	
	//status = CPXsetintparam (env,  CPX_PARAM_NODELIM, 10);
   	double timestamp1; double timestamp2; double diffmip;
	CPXgettime( env,  &timestamp1);
   	status = CPXmipopt (env, lp); 
	if(status) goto TERMINATE;
	CPXgettime( env,  &timestamp2);
	diffmip=timestamp2-timestamp1;


   int solstat = CPXgetstat (env, lp);
   printf ("Solution status %d.\n", solstat);
    double objval=-1;
   status = CPXgetobjval (env, lp, &objval);
   if ( status ) {
      fprintf (stderr,"Failed to obtain objective value.\n");
      goto TERMINATE;
   }   

   status = CPXgetx (env, lp, y, 0, numz+num_patients-1);
   if ( status ) {
      fprintf (stderr, "Failed to obtain solution.\n");
      goto TERMINATE;
   }
   double best_ub=-1;
    status = CPXgetbestobjval (env, lp, &best_ub);
   // Write out the solution 

  

   printf ("Objective value %.10g\n", objval);
   cout<< "error_counter is "<< error_counter<<endl;
   cout<< "counter is "<< counter<<endl;
   
#pragma endregion


TERMINATE:
#pragma region saving the computational results in a file
	finish = clock();
    diff = (finish - start);
    printf("\n");
    printf("Running Time = %0.4f \n ",diff);

	double gap=-1;
	status = CPXgetmiprelgap (env, lp, &gap);   
	printf("optimality gap is %17.10g \n",gap);
	printf("best bound is %17.10g \n",best_ub);

	int nodecount= CPXgetnodecnt(  env,  lp );
	int cutcount;
	status= CPXgetnumcuts( env, lp, CPX_CUT_USER, &cutcount ); 
	printf("the number of nodes is %d \n", nodecount);
	printf("the number of user cuts is %d \n", cutinfo.num_usercuts);
	

	ofstream myfile;  
	myfile.open("summary.txt", fstream::app);
	myfile <<instance_num<<setw(12);
	myfile <<num_states<<setw(12);
	for(int j=numz;j<numz+num_patients;j++)
	myfile << y[j]<<setw(24);
	myfile << best_ub<<setw(24);
		
	myfile << cutinfo.num_usercuts<<setw(24);
	myfile << nodecount<<setw(24);
	myfile << objval<<setw(24);
	myfile << gap<<setw(24);
	myfile << diff_mip_setup<<setw(24);
	myfile <<  diffmip<<setw(24);
	myfile <<  "Branch-and-Cut"<<setw(24)<<endl;
	myfile.close();
#pragma endregion

#pragma region deallocate memory

	if ( lp != NULL ) {
		status = CPXfreeprob (env, &lp);
		if ( status ) {
			fprintf (stderr, "CPXfreeprob failed, error code1 %d.\n",status);
		}
	}

	
	 
	for(int i=0;i<num_patients;i++)
	if ( slp[i] != NULL ) {
		status = CPXfreeprob (env_slp, &slp[i]);
		if ( status ) {
			fprintf (stderr, "CPXfreeprob failed , error code %d.\n",status);
		}
	}
	
	if ( slp_t != NULL ) {
		status = CPXfreeprob (env_slp, &slp_t);
		if ( status ) {
			fprintf (stderr, "CPXfreeprob failed, error code2 %d.\n",status);
		}
	}

	if ( slp_vt != NULL ) {
		status = CPXfreeprob (env_slp, &slp_vt);
		if ( status ) {
			fprintf (stderr, "CPXfreeprob failed, error code2 %d.\n",status);
		}
	}
	for(int i=0;i<num_patients;i++)
	if ( slp_v[i] != NULL ) {
		status = CPXfreeprob (env_slp, &slp_v[i]);
		if ( status ) {
			fprintf (stderr, "CPXfreeprob failed , error code %d.\n",status);
		}
	}
			// Free the CPLEX environment, if necessary 
	if ( env_slp != NULL ) {
		status = CPXcloseCPLEX (&env_slp);

		if ( status ) {
			char errmsg[1024];
			fprintf (stderr, "Could not close CPLEX environment.\n");
			CPXgeterrorstring (env_slp, status, errmsg);
			fprintf (stderr, "%s", errmsg);			
		}
		}

	if ( env != NULL ) {
		status = CPXcloseCPLEX (&env);
		if ( status ) {
			char errmsg[1024];
			fprintf (stderr, "Could not close CPLEX environment.\n");
			CPXgeterrorstring (env, status, errmsg);
			fprintf (stderr, "%s", errmsg);
		}
	}

	for (int k = 0; k < num_patients; k++)
	{
	for (int i = 0; i < num_states; i++)
		delete []tran_prob[k][i];
	delete[] tran_prob[k];
	delete[] reward[k];
	delete[] stop_reward[k];	
	delete[] rhs_stop[k];	
	delete[] rhs_contin[k];
	delete[]v[k];
	delete[] d[k];
	delete[]v_upd[k];	
	}
	delete[] feas_h_status;
	delete[] tran_prob;
	delete[] reward;
	delete[] stop_reward;
	delete[] rhs_stop;	
	delete[] rhs_contin;
	delete[]vt;
	delete[]v;
	delete[] d;
	delete[]v_upd;
	delete[]vt_upd;
	
	delete[] last_heur_solu;
	delete[]feas_h;
			
	delete []y;
	
#pragma endregion 
	}
	}
}

	void cplex(){
		
	#pragma omp parallel num_threads(num_cores)
	{
		int instance_num;
		#pragma omp for private(instance_num)
		for( instance_num=instance_beg;instance_num<=instance_end;instance_num++){

	#pragma region initialization

		const double base=num_states;
		numz=pow (base, num_patients);	

		int obj_ind;// This specifies the initial state of the game for consideration of objective value			
	
		double diff;
		clock_t start, finish;
		start = clock();
		/////////////////////////start reading data///////////////////////////////////////
		double ***tran_prob= new double ** [num_patients];
		for (int k = 0; k < num_patients; k++)
		{
			tran_prob[k]= new double *[num_states];
		for (int i = 0; i < num_states; i++)
			tran_prob[k][i] = new double [num_states];
		}

		
			//taking reward matrix
		double **reward = new double * [num_patients];
		for (int i = 0; i < num_patients; i++)
		reward[i] = new double [num_states];

				//taking stopping reward matrix
		double **stop_reward = new double * [num_patients];
		for (int i = 0; i < num_patients; i++)
		stop_reward[i] = new double [num_states];

		double **rhs_contin = new double * [num_patients];
		double **rhs_stop = new double * [num_patients];
		double **v= new double * [num_patients]; // to store V values
		double **d= new double * [num_patients]; // to store d values
		
		double *y=new double [numz*(1+num_patients)];	
		for (int i = 0; i < num_patients; i++){
			rhs_contin[i] = new double [numz];
			rhs_stop[i]= new double [numz];
			v[i]= new double [numz];
			d[i]= new double [numz];		
	
		}
	
	
		if( num_patients==2){
			readdata2( tran_prob, reward, stop_reward, &obj_ind, &instance_num );

			for(int i=0;i<num_patients;i++)
			for(int j1=0;j1<num_states;j1++)
			for(int j2=0;j2<num_states;j2++){
			if (j1!=num_states-1&& j2!=num_states-1)
				rhs_stop[i][j1*num_states+j2]=stop_reward[i][i*j2+(1-i)*j1];
			else
				rhs_stop[i][j1*num_states+j2]=0;
			rhs_contin[i][j1*num_states+j2]=reward[i][i*j2+(1-i)*j1];		
		}	

		}
		if( num_patients==3){
			readdata3( tran_prob, reward, stop_reward, &obj_ind, &instance_num );

			for(int j1=0;j1<num_states;j1++)
			for(int j2=0;j2<num_states;j2++)
			for(int j3=0;j3<num_states;j3++){
				if (j1!=num_states-1&& j2!=num_states-1&& j3!=num_states-1)
					rhs_stop[0][j1*num_states*num_states+j2*num_states+j3]=stop_reward[0][j1];
				else
					rhs_stop[0][j1*num_states*num_states+j2*num_states+j3]=0;
				rhs_contin[0][j1*num_states*num_states+j2*num_states+j3]=reward[0][j1];		
			}	

			for(int j1=0;j1<num_states;j1++)
			for(int j2=0;j2<num_states;j2++)
			for(int j3=0;j3<num_states;j3++){
				if (j1!=num_states-1&& j2!=num_states-1&& j3!=num_states-1)
					rhs_stop[1][j1*num_states*num_states+j2*num_states+j3]=stop_reward[1][j2];
				else
					rhs_stop[1][j1*num_states*num_states+j2*num_states+j3]=0;
				rhs_contin[1][j1*num_states*num_states+j2*num_states+j3]=reward[1][j2];		
			}	

			for(int j1=0;j1<num_states;j1++)
			for(int j2=0;j2<num_states;j2++)
			for(int j3=0;j3<num_states;j3++){
				if (j1!=num_states-1&& j2!=num_states-1&& j3!=num_states-1)
					rhs_stop[2][j1*num_states*num_states+j2*num_states+j3]=stop_reward[2][j3];
				else
					rhs_stop[2][j1*num_states*num_states+j2*num_states+j3]=0;
				rhs_contin[2][j1*num_states*num_states+j2*num_states+j3]=reward[2][j3];		
		}	

		}		
	

		/////////////////////////finsih reading data///////////////////////////////////////

	
		CPXENVptr env_slp;	
		CPXLPptr slp[num_patients];		
		CPXLPptr slp_v[num_patients];
		CPXLPptr exten;


		env_slp = CPXopenCPLEX (&status);
		if ( env_slp == NULL ) {
			char  errmsg[1024];
			cerr<<"Could not open CPLEX environment for subproblem.\n";
			CPXgeterrorstring (env_slp, status, errmsg);
			fprintf (stderr, "%s", errmsg);
			goto TERMINATE;
			}

		for(int i=0;i<num_patients;i++){
			slp[i]= CPXcreateprob (env_slp, &status, "subproblem1");
			if ( slp[i] == NULL ) {
			fprintf (stderr, "Failed to create subproblem.\n");
			goto TERMINATE;
			}
		}
		for(int i=0;i<num_patients;i++){
			slp_v[i]= CPXcreateprob (env_slp, &status, "subproblem1");
			if ( slp_v[i] == NULL ) {
			fprintf (stderr, "Failed to create subproblem.\n");
			goto TERMINATE;
			}
		}
	
		exten= CPXcreateprob (env_slp, &status, "subproblem2");
		if ( exten == NULL ) {
		fprintf (stderr, "Failed to create master problem.\n");
		goto TERMINATE;
		}

		
		double mip_setup1; double mip_setup2; double diff_mip_setup;
		CPXgettime( env_slp,  &mip_setup1); 

		if(num_patients==2){
			parameters_calulator_2( env_slp,slp,slp_v,tran_prob, rhs_contin, rhs_stop,v,d, &obj_ind);	
			make_extensive_2(env_slp,exten,tran_prob,rhs_contin,rhs_stop,v,d, &obj_ind);
		}

		if(num_patients==3){
			parameters_calulator_3( env_slp, slp,slp_v,tran_prob, rhs_contin, rhs_stop,v,d, &obj_ind);	
			make_extensive_3(env_slp,exten,tran_prob,rhs_contin,rhs_stop,v,d, &obj_ind);
		}

		CPXgettime( env_slp,  &mip_setup2);
		diff_mip_setup=mip_setup2-mip_setup1;
	#pragma endregion


	#pragma region	warm start
		if(num_patients==2){
			double *x= new double [numz];
			double **solu= new double *[num_patients];
			double *w_x= new double [2*numz*(1+num_patients)];
			int *indices=new int[2*numz*(1+num_patients)];
			CPXLPptr copy_slp[num_patients];
			for(int i=0;i<num_patients;i++){
				copy_slp[i] = CPXcloneprob (env_slp, slp[i], &status); 
				solu[i]= new double[numz];
			}
			int *beg=new int [2]; beg[0]=0;	beg[1]=numz*(1+num_patients);
			int *effort_level= new int [2]; effort_level[0]=1;effort_level[1]=1;
			for(int i=0;i<num_states;i++)
			for(int j=0;j<num_states;j++){
				indices[i*num_states+j]=i*num_states+j;
				if(i!=num_states-1 && j!= num_states-1 && rhs_stop[0][i*num_states+j]==v[0][i*num_states+j]&&rhs_stop[1][i*num_states+j]==v[1][i*num_states+j])
					x[i*num_states+j]=1;						
				else
					x[i*num_states+j]=0;			
			}


			feasibility_check( env_slp, slp, copy_slp, x,solu,rhs_stop);
			for(int j=0;j<numz;j++)
				w_x[j]=x[j];
			for(int i=0;i<num_patients;i++)	
			for(int j=0;j<numz;j++)
				w_x[numz*(i+1)+j]=solu[i][j];		
	
			for(int j=0;j<numz*(1+num_patients);j++)
				indices[j]=j;

			for(int j=0;j<numz;j++)
				x[j]=0;
			x[obj_ind]=1;

			for(int i=0;i<num_patients;i++)
			if ( slp[i] != NULL ) {
				status = CPXfreeprob (env_slp, &copy_slp[i]);		
			}

			for(int i=0;i<num_patients;i++)
			copy_slp[i] = CPXcloneprob (env_slp, slp[i], &status); 

			feasibility_check( env_slp, slp, copy_slp, x,solu,rhs_stop);
			for(int j=0;j<numz;j++)
				w_x[numz*(1+num_patients)+j]=x[j];
			for(int i=0;i<num_patients;i++)	
			for(int j=0;j<numz;j++)
				w_x[numz*(1+num_patients)+numz*(i+1)+j]=solu[i][j];		

			for(int j=0;j<numz*(1+num_patients);j++)
			indices[numz*(1+num_patients)+j]=j;	
			status = CPXaddmipstarts (env_slp, exten, 2, 2*numz*(1+num_patients), beg, indices,w_x, effort_level, NULL);
			if ( status ) goto TERMINATE;
			for(int i=0;i<num_patients;i++)
			if ( slp[i] != NULL ) {
				status = CPXfreeprob (env_slp, &copy_slp[i]);
				if ( status ) {
					fprintf (stderr, "CPXfreeprob failed , error code %d.\n",status);
				}
				delete[] solu[i];
			}
			delete[] beg;
			delete[] effort_level; 
			delete[] x;
			delete[] w_x;
			delete[] indices;
			delete[] solu;
		}

		if(num_patients==3){

			double *x= new double [numz];
			double **solu= new double *[num_patients];
			double *w_x= new double [2*numz*(1+num_patients)];
			int *indices=new int[2*numz*(1+num_patients)];
			CPXLPptr copy_slp[num_patients];
			for(int i=0;i<num_patients;i++){
				copy_slp[i] = CPXcloneprob (env_slp, slp[i], &status); 
				solu[i]= new double[numz];
			}
			int *beg=new int [2]; beg[0]=0;	beg[1]=numz*(1+num_patients);
			int *effort_level= new int [2]; effort_level[0]=1;effort_level[1]=1;
			for(int i1=0;i1<num_states;i1++)
			for(int i2=0;i2<num_states;i2++)
			for(int i3=0;i3<num_states;i3++){
				indices[i1*num_states*num_states+i2*num_states+i3]=i1*num_states*num_states+i2*num_states+i3;
				if(i1!=num_states-1 && i2!= num_states-1 && i3!= num_states-1 && rhs_stop[0][i1*num_states*num_states+i2*num_states+i3]==v[0][i1*num_states*num_states+i2*num_states+i3]&&rhs_stop[1][i1*num_states*num_states+i2*num_states+i3]==v[1][i1*num_states*num_states+i2*num_states+i3]&&rhs_stop[2][i1*num_states*num_states+i2*num_states+i3]==v[2][i1*num_states*num_states+i2*num_states+i3])
					x[i1*num_states*num_states+i2*num_states+i3]=1;						
				else
					x[i1*num_states*num_states+i2*num_states+i3]=0;			
			}


			feasibility_check( env_slp, slp, copy_slp, x,solu,rhs_stop);
			for(int j=0;j<numz;j++)
				w_x[j]=x[j];
			for(int i=0;i<num_patients;i++)	
			for(int j=0;j<numz;j++)
				w_x[numz*(i+1)+j]=solu[i][j];		
	
			for(int j=0;j<numz*(1+num_patients);j++)
				indices[j]=j;

			for(int j=0;j<numz;j++)
				x[j]=0;
			x[obj_ind]=1;

			for(int i=0;i<num_patients;i++)
			if ( slp[i] != NULL ) {
				status = CPXfreeprob (env_slp, &copy_slp[i]);		
			}

			for(int i=0;i<num_patients;i++)
			copy_slp[i] = CPXcloneprob (env_slp, slp[i], &status); 

			feasibility_check( env_slp, slp, copy_slp, x,solu,rhs_stop);
			for(int j=0;j<numz;j++)
				w_x[numz*(1+num_patients)+j]=x[j];
			for(int i=0;i<num_patients;i++)	
			for(int j=0;j<numz;j++)
				w_x[numz*(1+num_patients)+numz*(i+1)+j]=solu[i][j];		

			for(int j=0;j<numz*(1+num_patients);j++)
			indices[numz*(1+num_patients)+j]=j;	
			status = CPXaddmipstarts (env_slp, exten, 2, 2*numz*(1+num_patients), beg, indices,w_x, effort_level, NULL);
			if ( status ) goto TERMINATE;
			for(int i=0;i<num_patients;i++)
			if ( slp[i] != NULL ) {
				status = CPXfreeprob (env_slp, &copy_slp[i]);
				if ( status ) {
					fprintf (stderr, "CPXfreeprob failed , error code %d.\n",status);
				}
				delete[] solu[i];
			}
			delete[] beg;
			delete[] effort_level;
			delete[] x;
			delete[] w_x;
			delete[] indices;
			delete[] solu;
		}
	#pragma endregion
	


	#pragma region solving MIP and getting the solution

		status = CPXsetintparam (env_slp, CPX_PARAM_SCRIND, CPX_ON);
		status = CPXsetintparam (env_slp, CPX_PARAM_THREADS, 1);	
		status = CPXsetdblparam  (env_slp,  CPX_PARAM_TILIM, time_limit);// CPLEX time limit	
		double timestamp1; double timestamp2; double diffmip;
		CPXgettime( env_slp,  &timestamp1);
   		status = CPXmipopt (env_slp, exten); // this is solver of our MILP problem 
		if(status) goto TERMINATE;
		CPXgettime( env_slp,  &timestamp2);
		diffmip=timestamp2-timestamp1;
   
   
   
	   double objval=-1;
	   status = CPXgetobjval (env_slp, exten, &objval);
	   if ( status ) {
		  fprintf (stderr,"Failed to obtain objective value.\n");
		  goto TERMINATE;
	   }   

	   status = CPXgetx (env_slp, exten, y, 0, numz*(1+num_patients)-1);
	   if ( status ) {
		  fprintf (stderr, "Failed to obtain solution.\n");
		  goto TERMINATE;
	   }
	   double best_ub=-1;
		status = CPXgetbestobjval (env_slp, exten, &best_ub);
	   // Write out the solution 


	   printf ("Objective value %.10g\n", objval);
   
	#pragma endregion

	TERMINATE:
	
	#pragma region saving the computational results in a file
		finish = clock();
		diff = (finish - start);
		printf("\n");
		printf("Running Time = %0.4f \n ",diff);

		double gap=-1;
		status = CPXgetmiprelgap (env_slp, exten, &gap);   
		printf("optimality gap is %17.10g \n",gap);
		printf("best bound is %17.10g \n",best_ub);

	
		int cutcount=0;
		for (int j=0;j<15;j++){
		int k;
		status= CPXgetnumcuts( env_slp, exten, j, &k ); 
		cutcount=cutcount+k;
		}
		int nodecount= CPXgetnodecnt(  env_slp,  exten );
		printf("the number of nodes is %d \n", nodecount);
		printf("the total number of  cuts is %d \n", cutcount);

		ofstream myfile;  

		myfile.open("summary.txt", fstream::app);
		myfile <<instance_num<<setw(12);
		myfile <<num_states<<setw(12);
	
		for( int j=0;j<num_patients;j++)
		myfile << y[numz*(j+1)+obj_ind]<<setw(24);
	
		myfile << best_ub<<setw(24); 
		
		myfile << cutcount<<setw(24);
		myfile << nodecount<<setw(24);	
		myfile << objval<<setw(24);
		myfile << gap<<setw(24);	
		myfile << diff_mip_setup<<setw(24);		
		myfile << diffmip<<setw(24);		
		myfile << "CPLEX"<<endl;
		myfile.close();
	#pragma endregion

		#pragma region deallocating memory


		if ( exten != NULL ) {
			status = CPXfreeprob (env_slp, &exten);
			if ( status ) {
				fprintf (stderr, "CPXfreeprob failed, error code2 %d.\n",status);
			}
		}

		
		for(int i=0;i<num_patients;i++)
		if ( slp[i] != NULL ) {
			status = CPXfreeprob (env_slp, &slp[i]);
			if ( status ) {
				fprintf (stderr, "CPXfreeprob failed , error code %d.\n",status);
			}
		}
	
		
		for(int i=0;i<num_patients;i++)
		if ( slp_v[i] != NULL ) {
			status = CPXfreeprob (env_slp, &slp_v[i]);
			if ( status ) {
				fprintf (stderr, "CPXfreeprob failed , error code %d.\n",status);
			}
		}
	
	
				// Free the CPLEX environment, if necessary 
		if ( env_slp != NULL ) {
			status = CPXcloseCPLEX (&env_slp);

			if ( status ) {
				char errmsg[1024];
				fprintf (stderr, "Could not close CPLEX environment.\n");
				CPXgeterrorstring (env_slp, status, errmsg);
				fprintf (stderr, "%s", errmsg);			
			}
			}

	

		for (int k = 0; k < num_patients; k++)
		{
			for (int i = 0; i < num_states; i++)
				delete []tran_prob[k][i];

			delete[] tran_prob[k];
			delete[] reward[k];
			delete[] stop_reward[k];	
			delete[] rhs_stop[k];	
			delete[] rhs_contin[k];
			delete[]v[k];
			delete[] d[k];

		}

		delete[] tran_prob;
		delete[] reward;
		delete[] stop_reward;
		delete[] rhs_stop;	
		delete[] rhs_contin;
		delete[]v;
		delete[] d;
		delete []y;
	
		#pragma endregion
		}
		}
	}

	void negative_reward_branch_and_cut(){
		#pragma omp parallel num_threads(num_cores)
		{
		int instance_num;
		#pragma omp for private(instance_num)
		for( instance_num=instance_beg;instance_num<=instance_end;instance_num++){
	#pragma region initialization
	
		const double base=num_states;
		numz=pow (base, num_patients);	

		int obj_ind;// This specifies the initial state of the game for consideration of objective value			
		int num_usercuts=0;
		bool *feas_h_status= new bool [1];
		feas_h_status[0]=false;		
		int last_node=-1;	

		double diff;
		clock_t start, finish;
		start = clock();
	#pragma region start reading data
		double ***tran_prob= new double ** [num_patients];
		for (int k = 0; k < num_patients; k++)
		{
			tran_prob[k]= new double *[num_states];
		for (int i = 0; i < num_states; i++)
			tran_prob[k][i] = new double [num_states];
		}

		
			//taking reward matrix
		double **reward = new double * [num_patients];
		for (int i = 0; i < num_patients; i++)
		reward[i] = new double [num_states];

				//taking stopping reward matrix
		double **stop_reward = new double * [num_patients];
		for (int i = 0; i < num_patients; i++)
		stop_reward[i] = new double [num_states];

		double **rhs_contin = new double * [num_patients];
		double **rhs_stop = new double * [num_patients];
		double *vt= new double [numz];
		double **v= new double * [num_patients]; // to store V values
		double **d= new double * [num_patients]; // to store d values
		double **v_upd= new double * [num_patients]; // to store updated V's	
		double *vt_upd= new double [numz]; // to store updated V's	
	
	
		double *last_heur_solu= new double [numz]; // to store solution of last heuristic in order to do not double check the same solution for next iteration
		double *feas_h= new double [numz+num_patients]; // feasibility cut obtained in heuristic
	
	
	
		double *y=new double [numz+num_patients];	
		for (int i = 0; i < num_patients; i++){
			rhs_contin[i] = new double [numz];
			rhs_stop[i]= new double [numz];
			v[i]= new double [numz];
			d[i]= new double [numz];		
			v_upd[i]= new double [numz];
		}
	
	
		readdata_negative( tran_prob, reward, stop_reward, &obj_ind, &instance_num );
	

		for(int i=0;i<num_patients;i++)
		for(int j1=0;j1<num_states;j1++)
		for(int j2=0;j2<num_states;j2++){
			if (j1!=num_states-1&& j2!=num_states-1)
				rhs_stop[i][j1*num_states+j2]=stop_reward[i][i*j2+(1-i)*j1];
			else
				rhs_stop[i][j1*num_states+j2]=0;
			rhs_contin[i][j1*num_states+j2]=reward[i][i*j2+(1-i)*j1];		
		}	

	

		#pragma endregion
	/////////////////////////finsih reading data///////////////////////////////////////


		CPXENVptr env;
		CPXENVptr env_slp;
		CPXLPptr lp;
		CPXLPptr slp[num_patients];
		CPXLPptr slp_vt;
		CPXLPptr slp_t;
		CPXLPptr slp_v[num_patients];// for calculation of V dynamically
	
		env = CPXopenCPLEX (&status);
		if ( env == NULL ) {
		char  errmsg[1024];
		cout<<"Could not open CPLEX environment.\n";
		CPXgeterrorstring (env, status, errmsg);
		fprintf (stderr, "%s", errmsg);
		goto TERMINATE;
		}	
		status = CPXsetintparam (env, CPX_PARAM_THREADS, 1);
		//status = CPXsetintparam(env, CPX_PARAM_NODEFILEIND, 2);// CPLEX can use node files and write on hard disk
		//status = CPXsetdblparam(env, CPX_PARAM_WORKMEM, 2048); // the number of memory units available for CPLEX

		env_slp = CPXopenCPLEX (&status);
		if ( env_slp == NULL ) {
			char  errmsg[1024];
			cerr<<"Could not open CPLEX environment for subproblem.\n";
			CPXgeterrorstring (env_slp, status, errmsg);
			fprintf (stderr, "%s", errmsg);
			goto TERMINATE;
			}
		lp= CPXcreateprob (env, &status, "master problem");
		if ( lp == NULL ) {
		fprintf (stderr, "Failed to create master problem.\n");
		goto TERMINATE;
		}
		for(int i=0;i<num_patients;i++){
			slp[i]= CPXcreateprob (env_slp, &status, "subproblem1");
			if ( slp[i] == NULL ) {
			fprintf (stderr, "Failed to create subproblem.\n");
			goto TERMINATE;
			}
		}
		slp_vt= CPXcreateprob (env_slp, &status, "slp_vt");// we create a problem 
		if ( slp_vt == NULL ) {
		fprintf (stderr, "Failed to create subproblem.\n");
		goto TERMINATE;
		}
		for(int i=0;i<num_patients;i++){
			slp_v[i]= CPXcreateprob (env_slp, &status, "subproblem1");
			if ( slp_v[i] == NULL ) {
			fprintf (stderr, "Failed to create subproblem.\n");
			goto TERMINATE;
			}
		}
	
		slp_t= CPXcreateprob (env_slp, &status, "subproblem2");
		if ( slp_t == NULL ) {
		fprintf (stderr, "Failed to create master problem.\n");
		goto TERMINATE;
		}

	
		double mip_setup1; double mip_setup2; double diff_mip_setup;
		CPXgettime( env,  &mip_setup1);   

		make_master_negative(env, env_slp, lp,slp,slp_v,tran_prob, rhs_contin, rhs_stop,v,d, &obj_ind);
	
		make_slpvt_negative(env_slp,slp_t,slp_vt,tran_prob,rhs_contin,rhs_stop,vt,d, &obj_ind);
	
		CPXgettime( env,  &mip_setup2);
		diff_mip_setup=mip_setup2-mip_setup1;

	#pragma endregion 
	
	#pragma region Set parameters

	
		/// Assure linear mappings between the presolved and originalmodels 
		status = CPXsetintparam (env, CPX_PARAM_PRELINEAR, CPX_OFF);	

		status = CPXsetintparam (env, CPX_PARAM_SCRIND, CPX_ON);		
	
		status = CPXsetintparam (env, CPX_PARAM_MIPSEARCH, CPX_MIPSEARCH_TRADITIONAL);
	
		// Let MIP callbacks work on the original model 
	
		status = CPXsetintparam (env, CPX_PARAM_MIPCBREDLP, CPX_OFF);
	
		status = CPXsetintparam (env, CPX_PARAM_MIPDISPLAY, 2);	


		status = CPXsetintparam (env,  CPX_PARAM_SUBALG, CPX_ALG_PRIMAL);
	
		status = CPXsetdblparam  (env, CPX_PARAM_CUTSFACTOR , 100000.0);
	
	
		status = CPXsetintparam (env, CPX_PARAM_VARSEL , CPX_VARSEL_STRONG);
	
		status = CPXsetintparam (env, CPX_PARAM_NODESEL , CPX_NODESEL_BESTBOUND);

		CPXsetdblparam (env, CPX_PARAM_BTTOL , 0.0001);

		status = CPXsetintparam (env,  CPX_PARAM_HEURFREQ , -1);

		status = CPXsetintparam (env, CPX_PARAM_REDUCE , 0 ); 
	
		status = CPXsetintparam (env_slp,  CPX_PARAM_LPMETHOD  , CPX_ALG_DUAL ); 
	#pragma endregion	

	#pragma region	initializing heuristic
		for(int j=0;j<numz;j++)
		feas_h[j]=0;

		if(rhs_stop[0][obj_ind] >= d[0][obj_ind] && rhs_stop[1][obj_ind] >= d[1][obj_ind]){
			feas_h[obj_ind]=1;
			feas_h[numz]=rhs_stop[0][obj_ind]; feas_h[numz+1]=rhs_stop[1][obj_ind];
		}
		else{
			feas_h[numz]=d[0][obj_ind]; feas_h[numz+1]=d[1][obj_ind];
		}
	#pragma endregion

	#pragma region	warm start
		double *x= new double [2*(numz+num_patients)];
		double **solu= new double *[num_patients];
		int *indices=new int[2*(numz+num_patients)];
		CPXLPptr copy_slp[num_patients];
		for(int i=0;i<num_patients;i++){
			copy_slp[i] = CPXcloneprob (env_slp, slp[i], &status); 
			solu[i]= new double[numz];
		}
		int *beg=new int [2]; beg[0]=0;	beg[1]=numz+num_patients;
		int *effort_level= new int [2]; effort_level[0]=1;effort_level[1]=1;
		for(int i=0;i<num_states;i++)
		for(int j=0;j<num_states;j++){
			indices[i*num_states+j]=i*num_states+j;
			if(i!=num_states-1 && j!= num_states-1 && rhs_stop[0][i*num_states+j]==v[0][i*num_states+j]&&rhs_stop[1][i*num_states+j]==v[1][i*num_states+j])
				x[i*num_states+j]=1;						
			else
				x[i*num_states+j]=0;			
		}


		feasibility_check( env_slp, slp, copy_slp, x,solu,rhs_stop);
		for(int j=0;j<num_patients;j++){
			x[numz+j]=solu[j][obj_ind];
			indices[numz+j]=numz+j;
		}
		for(int j=0;j<numz+num_patients;j++){
			x[numz+num_patients+j]=feas_h[j];
			indices[numz+num_patients+j]=j;
		}
		status = CPXaddmipstarts (env, lp, 2, 2*(numz+num_patients), beg, indices,x, effort_level, NULL);
		if ( status ) goto TERMINATE;
		for(int i=0;i<num_patients;i++)
		if ( slp[i] != NULL ) {
			status = CPXfreeprob (env_slp, &copy_slp[i]);
			if ( status ) {
				fprintf (stderr, "CPXfreeprob failed , error code %d.\n",status);
			}
			delete[] solu[i];
		}

		delete[] solu;
	#pragma endregion

	#pragma region construct initialization
		CUTINFO cutinfo;

		cutinfo.env_slp=env_slp;			
		cutinfo.slp_t=slp_t;
		cutinfo.slp=slp;
		cutinfo.slp_vt=slp_vt;
		cutinfo.slp_v=slp_v;
		cutinfo.rhs_contin=rhs_contin;
		cutinfo.rhs_stop=rhs_stop;
	
		cutinfo.v=v;
		cutinfo.vt=vt;
		cutinfo.v_upd=v_upd;
		cutinfo.vt_upd=vt_upd;
		cutinfo.d=d;
		cutinfo.feas_h=feas_h;	
		cutinfo.obj_ind=obj_ind;
		cutinfo.num_usercuts=num_usercuts;
		cutinfo.feas_h_status=feas_h_status;
		cutinfo.last_node=last_node;

		BRANCHINFO branchinfo;
		branchinfo.env_slp=env_slp;
		branchinfo.slp=slp;
		branchinfo.rhs_stop=rhs_stop;
		branchinfo.last_heur_solu=last_heur_solu;
		branchinfo.feas_h=feas_h;
		
		branchinfo.obj_ind=obj_ind;
		branchinfo.feas_h_status=feas_h_status;
	#pragma endregion

	#pragma region	Set up to use MIP callback 
		
		status = CPXsetbranchcallbackfunc (env, mybranchfunc, &branchinfo); 	
	   if ( status )  goto TERMINATE;
      
		status = CPXsetusercutcallbackfunc(env, usercutcallback, &cutinfo);      
	   if ( status )  goto TERMINATE;

	   status = CPXsetlazyconstraintcallbackfunc (env, mycutcallback, &cutinfo); 
	   if ( status )  goto TERMINATE;   

	   status = CPXsetincumbentcallbackfunc (env, myincumbentcheck,&cutinfo);
	   if ( status ) goto TERMINATE;
   
	   status = CPXsetheuristiccallbackfunc (env, myheuristic, &branchinfo);
	   if ( status )  goto TERMINATE;
	#pragma endregion   

	#pragma region solving MIP and getting the solution

		status = CPXsetdblparam  (env,  CPX_PARAM_TILIM , time_limit);// CPLEX time limit
		//status = CPXsetintparam (env,  CPX_PARAM_NODELIM, 10);
   		double timestamp1; double timestamp2; double diffmip;
		CPXgettime( env,  &timestamp1);
		status = CPXmipopt (env, lp); // this is solver of our MILP problem 
		if(status) goto TERMINATE;
		CPXgettime( env,  &timestamp2);
		diffmip=timestamp2-timestamp1;


	   int solstat = CPXgetstat (env, lp);
	   printf ("Solution status %d.\n", solstat);
		double objval=-1;
	   status = CPXgetobjval (env, lp, &objval);
	   if ( status ) {
		  fprintf (stderr,"Failed to obtain objective value.\n");
		  goto TERMINATE;
	   }   

	   status = CPXgetx (env, lp, y, 0, numz+num_patients-1);
	   if ( status ) {
		  fprintf (stderr, "Failed to obtain solution.\n");
		  goto TERMINATE;
	   }
	   double best_ub=-1;
		status = CPXgetbestobjval (env, lp, &best_ub);



	   printf ("Objective value %.10g\n", objval);
	   cout<< "error_counter is "<< error_counter<<endl;
	   cout<< "counter is "<< counter<<endl;
  
  
		double up_bo=0;
		double up_bo_cur;
		for(int i=0;i<num_patients;i++)
		for(int j1=0;j1<num_states;j1++)
		for(int j2=0;j2<num_states;j2++){
			if((j1!=num_states-1)&&(j2!=num_states)){
				up_bo_cur=rhs_stop[i][j1*num_states+j2]/d[i][j1*num_states+j2];
				if(up_bo_cur>up_bo)
					up_bo=up_bo_cur;
			}
		}
		up_bo=(d[0][obj_ind]+d[1][obj_ind])*up_bo;
		printf("upper bound is %17.10g \n",up_bo);
	#pragma endregion

	TERMINATE:
	#pragma region saving the computational results in a file
		finish = clock();
		diff = (finish - start);
		printf("\n");
		printf("Running Time = %0.4f \n ",diff);

		double gap=-1;
		status = CPXgetmiprelgap (env, lp, &gap);   
		printf("optimality gap is %17.10g \n",gap);
		printf("best bound is %17.10g \n",best_ub);

		int nodecount= CPXgetnodecnt(  env,  lp );
		int cutcount;
		status= CPXgetnumcuts( env, lp, CPX_CUT_USER, &cutcount ); 
		printf("the number of nodes is %d \n", nodecount);
		printf("the number of user cuts is %d \n", cutinfo.num_usercuts);
	

		ofstream myfile;  
		myfile.open("summary.txt", fstream::app);
		myfile <<instance_num<<setw(12);
		myfile <<num_states<<setw(12);
		for(int j=numz;j<numz+num_patients;j++)
		myfile << y[j]<<setw(24);
		myfile << best_ub<<setw(24);
	

		myfile << cutinfo.num_usercuts<<setw(24);
		myfile << nodecount<<setw(24);
		myfile << objval<<setw(24);
		myfile << gap<<setw(24);
		myfile << diff_mip_setup<<setw(24);
		myfile <<  diffmip<<setw(50);
		myfile <<  "Negative-reward-Branch-and-Cut"<<endl;
		myfile.close();
		#pragma endregion

		#pragma region deallocating memory
		if ( lp != NULL ) {
			status = CPXfreeprob (env, &lp);
			if ( status ) {
				fprintf (stderr, "CPXfreeprob failed, error code1 %d.\n",status);
			}
		}

	
	 
		for(int i=0;i<num_patients;i++)
		if ( slp[i] != NULL ) {
			status = CPXfreeprob (env_slp, &slp[i]);
			if ( status ) {
				fprintf (stderr, "CPXfreeprob failed , error code %d.\n",status);
			}
		}
	
		if ( slp_t != NULL ) {
			status = CPXfreeprob (env_slp, &slp_t);
			if ( status ) {
				fprintf (stderr, "CPXfreeprob failed, error code2 %d.\n",status);
			}
		}

		if ( slp_vt != NULL ) {
			status = CPXfreeprob (env_slp, &slp_vt);
			if ( status ) {
				fprintf (stderr, "CPXfreeprob failed, error code2 %d.\n",status);
			}
		}
		for(int i=0;i<num_patients;i++)
		if ( slp_v[i] != NULL ) {
			status = CPXfreeprob (env_slp, &slp_v[i]);
			if ( status ) {
				fprintf (stderr, "CPXfreeprob failed , error code %d.\n",status);
			}
		}
	
				// Free the CPLEX environment, if necessary 
		if ( env_slp != NULL ) {
			status = CPXcloseCPLEX (&env_slp);

			if ( status ) {
				char errmsg[1024];
				fprintf (stderr, "Could not close CPLEX environment.\n");
				CPXgeterrorstring (env_slp, status, errmsg);
				fprintf (stderr, "%s", errmsg);			
			}
			}

		if ( env != NULL ) {
			status = CPXcloseCPLEX (&env);
			if ( status ) {
				char errmsg[1024];
				fprintf (stderr, "Could not close CPLEX environment.\n");
				CPXgeterrorstring (env, status, errmsg);
				fprintf (stderr, "%s", errmsg);
			}
		}

		for (int k = 0; k < num_patients; k++)
		{
		for (int i = 0; i < num_states; i++)
			delete []tran_prob[k][i];

		delete[] tran_prob[k];
		delete[] reward[k];
		delete[] stop_reward[k];	
		delete[] rhs_stop[k];	
		delete[] rhs_contin[k];
		delete[]v[k];
		delete[] d[k];
		delete[]v_upd[k];	
		}
		delete[] feas_h_status;
		delete[] tran_prob;
		delete[] reward;
		delete[] stop_reward;
		delete[] rhs_stop;	
		delete[] rhs_contin;
		delete[]vt;
		delete[]v;
		delete[] d;
		delete[]v_upd;
		delete[]vt_upd;
	
		delete[] last_heur_solu;
		delete[]feas_h;	
		delete []y;
	
	
		delete[] beg;
		delete[] effort_level;
		delete[] x;
		delete[] indices;
		#pragma endregion
		}
		}
	}

	void negative_reward_cplex(){

	#pragma omp parallel num_threads(num_cores)
		{
		int instance_num;
		#pragma omp for private(instance_num)
		for( instance_num=instance_beg;instance_num<=instance_end;instance_num++){
	#pragma region initialization
		const double base=num_states;
		numz=pow (base, num_patients);	
			
		int obj_ind;// This specifies the initial state of the game for consideration of objective value			
	
		double diff;
		clock_t start, finish;
		start = clock();
	#pragma region start reading data
		double ***tran_prob= new double ** [num_patients];
		for (int k = 0; k < num_patients; k++)
		{
			tran_prob[k]= new double *[num_states];
		for (int i = 0; i < num_states; i++)
			tran_prob[k][i] = new double [num_states];
		}

		
			//taking reward matrix
		double **reward = new double * [num_patients];
		for (int i = 0; i < num_patients; i++)
		reward[i] = new double [num_states];

				//taking stopping reward matrix
		double **stop_reward = new double * [num_patients];
		for (int i = 0; i < num_patients; i++)
		stop_reward[i] = new double [num_states];

		double **rhs_contin = new double * [num_patients];
		double **rhs_stop = new double * [num_patients];
	
		double **v= new double * [num_patients]; // to store V values
		double **d= new double * [num_patients]; // to store d values
	
	
		double *y=new double [numz*(1+num_patients)];	
		for (int i = 0; i < num_patients; i++){
			rhs_contin[i] = new double [numz];
			rhs_stop[i]= new double [numz];
			v[i]= new double [numz];
			d[i]= new double [numz];	
		}
	
	
		readdata_negative( tran_prob, reward, stop_reward, &obj_ind, &instance_num );

		for(int i=0;i<num_patients;i++)
		for(int j1=0;j1<num_states;j1++)
		for(int j2=0;j2<num_states;j2++){
			if (j1!=num_states-1&& j2!=num_states-1)
				rhs_stop[i][j1*num_states+j2]=stop_reward[i][i*j2+(1-i)*j1];
			else
				rhs_stop[i][j1*num_states+j2]=0;
			rhs_contin[i][j1*num_states+j2]=reward[i][i*j2+(1-i)*j1];		
		}	

	

		#pragma endregion
	/////////////////////////finsih reading data///////////////////////////////////////


		CPXENVptr env_slp;
		CPXLPptr slp[num_patients];
		CPXLPptr slp_v[num_patients];
		CPXLPptr exten;

		env_slp = CPXopenCPLEX (&status);
		if ( env_slp == NULL ) {
			char  errmsg[1024];
			cerr<<"Could not open CPLEX environment for subproblem.\n";
			CPXgeterrorstring (env_slp, status, errmsg);
			fprintf (stderr, "%s", errmsg);
			goto TERMINATE;
			}


		for(int i=0;i<num_patients;i++){
			slp[i]= CPXcreateprob (env_slp, &status, "subproblem1");
			if ( slp[i] == NULL ) {
			fprintf (stderr, "Failed to create subproblem.\n");
			goto TERMINATE;
			}
		}
	
		for(int i=0;i<num_patients;i++){
			slp_v[i]= CPXcreateprob (env_slp, &status, "subproblem1");
			if ( slp_v[i] == NULL ) {
			fprintf (stderr, "Failed to create subproblem.\n");
			goto TERMINATE;
			}
		}
	
		exten= CPXcreateprob (env_slp, &status, "subproblem2");
		if ( exten == NULL ) {
		fprintf (stderr, "Failed to create master problem.\n");
		goto TERMINATE;
		}

	
	
		double mip_setup1; double mip_setup2; double diff_mip_setup;
		CPXgettime( env_slp,  &mip_setup1); 

		parameters_calulator_negative(env_slp,slp,slp_v,tran_prob, rhs_contin, rhs_stop,v,d, &obj_ind);	
		make_extensive_negative(env_slp,exten,tran_prob,rhs_contin,rhs_stop,v,d, &obj_ind);
	
		CPXgettime( env_slp,  &mip_setup2);
		diff_mip_setup=mip_setup2-mip_setup1;

	#pragma endregion 
	
	#pragma region	warm start
		double *x= new double [numz];
		double **solu= new double *[num_patients];
		double *w_x= new double [2*numz*(1+num_patients)];
		int *indices=new int[2*numz*(1+num_patients)];
		CPXLPptr copy_slp[num_patients];
		for(int i=0;i<num_patients;i++){
			copy_slp[i] = CPXcloneprob (env_slp, slp[i], &status); 
			solu[i]= new double[numz];
		}
		int *beg=new int [2]; beg[0]=0;	beg[1]=numz*(1+num_patients);
		int *effort_level= new int [2]; effort_level[0]=1;effort_level[1]=1;
		for(int i=0;i<num_states;i++)
		for(int j=0;j<num_states;j++){
			indices[i*num_states+j]=i*num_states+j;
			if(i!=num_states-1 && j!= num_states-1 && rhs_stop[0][i*num_states+j]==v[0][i*num_states+j]&&rhs_stop[1][i*num_states+j]==v[1][i*num_states+j])
				x[i*num_states+j]=1;						
			else
				x[i*num_states+j]=0;			
		}


		feasibility_check( env_slp, slp, copy_slp, x,solu,rhs_stop);
		for(int j=0;j<numz;j++)
			w_x[j]=x[j];
		for(int i=0;i<num_patients;i++)	
		for(int j=0;j<numz;j++)
			w_x[numz*(i+1)+j]=solu[i][j];		
	
		for(int j=0;j<numz*(1+num_patients);j++)
			indices[j]=j;

		for(int j=0;j<numz;j++)
			x[j]=0;
		x[obj_ind]=1;

		for(int i=0;i<num_patients;i++)
		if ( slp[i] != NULL ) {
			status = CPXfreeprob (env_slp, &copy_slp[i]);		
		}

		for(int i=0;i<num_patients;i++)
		copy_slp[i] = CPXcloneprob (env_slp, slp[i], &status); 

		feasibility_check( env_slp, slp, copy_slp, x,solu,rhs_stop);
		for(int j=0;j<numz;j++)
			w_x[numz*(1+num_patients)+j]=x[j];
		for(int i=0;i<num_patients;i++)	
		for(int j=0;j<numz;j++)
			w_x[numz*(1+num_patients)+numz*(i+1)+j]=solu[i][j];		

		for(int j=0;j<numz*(1+num_patients);j++)
		indices[numz*(1+num_patients)+j]=j;	
		status = CPXaddmipstarts (env_slp, exten, 2, 2*numz*(1+num_patients), beg, indices,w_x, effort_level, NULL);
		if ( status ) goto TERMINATE;
		for(int i=0;i<num_patients;i++)
		if ( slp[i] != NULL ) {
			status = CPXfreeprob (env_slp, &copy_slp[i]);
			if ( status ) {
				fprintf (stderr, "CPXfreeprob failed , error code %d.\n",status);
			}
			delete[] solu[i];
		}
		delete[] beg;
		delete[] effort_level;
		delete[] x;
		delete[] w_x;
		delete[] indices;
		delete[] solu;
	#pragma endregion

	#pragma region solving MIP and getting the solution

		status = CPXsetintparam (env_slp, CPX_PARAM_SCRIND, CPX_ON);
		status = CPXsetintparam (env_slp, CPX_PARAM_THREADS, 1);
		status = CPXsetdblparam  (env_slp,  CPX_PARAM_TILIM , time_limit);// CPLEX time limit
		double timestamp1; double timestamp2; double diffmip;
		CPXgettime( env_slp,  &timestamp1);
   		status = CPXmipopt (env_slp, exten); // this is solver of our MILP problem 
		if(status) goto TERMINATE;
		CPXgettime( env_slp,  &timestamp2);
		diffmip=timestamp2-timestamp1;
   
	   double objval=-1;
	   status = CPXgetobjval (env_slp, exten, &objval);
	   if ( status ) {
		  fprintf (stderr,"Failed to obtain objective value.\n");
		  goto TERMINATE;
	   }   

	   status = CPXgetx (env_slp, exten, y, 0, numz*(1+num_patients)-1);
	   if ( status ) {
		  fprintf (stderr, "Failed to obtain solution.\n");
		  goto TERMINATE;
	   }
	   double best_ub=-1;
		status = CPXgetbestobjval (env_slp, exten, &best_ub);


	   printf ("Objective value %.10g\n", objval);
	   cout<< "error_counter is "<< error_counter<<endl;
	   cout<< "counter is "<< counter<<endl;
	   if(v[0][obj_ind]==rhs_stop[0][obj_ind]&&v[1][obj_ind]==rhs_stop[1][obj_ind]) printf("easy instance");
	#pragma endregion

	TERMINATE:
		#pragma region saving the computational results in a file
		finish = clock();
		diff = (finish - start);
		printf("\n");
		printf("Running Time = %0.4f \n ",diff);

		double gap=-1;
		status = CPXgetmiprelgap (env_slp, exten, &gap);   
		printf("optimality gap is %17.10g \n",gap);
		printf("best bound is %17.10g \n",best_ub);

	
		int cutcount=0;
		for (int j=0;j<15;j++){
		int k;
		status= CPXgetnumcuts( env_slp, exten, j, &k ); 
		cutcount=cutcount+k;
		}
		int nodecount= CPXgetnodecnt(  env_slp,  exten );
		printf("the number of nodes is %d \n", nodecount);
		printf("the total number of  cuts is %d \n", cutcount);

		ofstream myfile;  
		myfile.open("summary.txt", fstream::app);
		myfile <<instance_num<<setw(12);
		myfile <<num_states<<setw(12);
	
		myfile << y[numz+obj_ind]<<setw(24);
		myfile << y[2*numz+obj_ind]<<setw(24);
		myfile << best_ub<<setw(24);
	

		myfile << cutcount<<setw(24);
		myfile << nodecount<<setw(24);	
		myfile << objval<<setw(24);
		myfile << gap<<setw(24);	
		myfile << diff_mip_setup<<setw(24);		
		myfile << diffmip<<setw(50);
		myfile << "Negative-reward-CPLEX"<<endl;
		myfile.close();
		#pragma endregion

		#pragma region deallocating memory
	

		if ( exten != NULL ) {
			status = CPXfreeprob (env_slp, &exten);
			if ( status ) {
				fprintf (stderr, "CPXfreeprob failed, error code2 %d.\n",status);
			}
		}


	 
		for(int i=0;i<num_patients;i++)
		if ( slp[i] != NULL ) {
			status = CPXfreeprob (env_slp, &slp[i]);
			if ( status ) {
				fprintf (stderr, "CPXfreeprob failed , error code %d.\n",status);
			}
		}
	
	
		for(int i=0;i<num_patients;i++)
		if ( slp_v[i] != NULL ) {
			status = CPXfreeprob (env_slp, &slp_v[i]);
			if ( status ) {
				fprintf (stderr, "CPXfreeprob failed , error code %d.\n",status);
			}
		}
	
				// Free the CPLEX environment, if necessary 
		if ( env_slp != NULL ) {
			status = CPXcloseCPLEX (&env_slp);

			if ( status ) {
				char errmsg[1024];
				fprintf (stderr, "Could not close CPLEX environment.\n");
				CPXgeterrorstring (env_slp, status, errmsg);
				fprintf (stderr, "%s", errmsg);			
			}
			}

	

		for (int k = 0; k < num_patients; k++)
		{
		for (int i = 0; i < num_states; i++)
			delete []tran_prob[k][i];

		delete[] tran_prob[k];
		delete[] reward[k];
		delete[] stop_reward[k];	
		delete[] rhs_stop[k];	
		delete[] rhs_contin[k];
		delete[]v[k];
		delete[] d[k];
		}
		delete[] tran_prob;
		delete[] reward;
		delete[] stop_reward;
		delete[] rhs_stop;	
		delete[] rhs_contin;
		delete[]v;
		delete[] d;
		delete []y;
		#pragma endregion
		}
		}

	}