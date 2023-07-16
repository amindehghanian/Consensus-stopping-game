#pragma region libraries
	#include <iostream> 	
	#include <ilcplex/cplex.h>	
	#include "Cut Generator.h"
	#include"Make problems.h"
	
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
	#include <algorithm>	
	using namespace std;	
#pragma endregion	



int make_master3
		(CPXCENVptr env,CPXCENVptr env_slp,CPXLPptr lp, CPXLPptr* slp, CPXLPptr* slp_v, double ***tran_prob, double **rhs_contin,double **rhs_stop, double **v,double **d, int *obj_ind)//making master problem	
	{		
#pragma region  start creating sub_problems for the  patients and calculating d_i(s)
	double *matval= new double [numz*numz];
	int *matind= new int [numz*numz];
	int *matbeg= new int [numz];
	int *matcnt= new int [numz];
	char *sense= new char [numz];
	double *obj= new double [numz];
	double *rhs= new double [numz];
	double *lb= new double [numz];
	double *ub= new double [numz];
	int *indices = new int[numz];
	

	for(int i1=0;i1<num_states;i1++)
	for(int i2=0;i2<num_states;i2++)
	for(int i3=0;i3<num_states;i3++)
	{
		for(int j1=0;j1<num_states;j1++)
		for(int j2=0;j2<num_states;j2++)
		for(int j3=0;j3<num_states;j3++)
		{
			matval[(i1*num_states*num_states+i2*num_states+i3)*numz+j1*num_states*num_states+j2*num_states+j3]=(-lambda)*tran_prob[0][j1][i1]*tran_prob[1][j2][i2]*tran_prob[2][j3][i3];
			matind[(i1*num_states*num_states+i2*num_states+i3)*numz+j1*num_states*num_states+j2*num_states+j3]=j1*num_states*num_states+j2*num_states+j3;

		}
		matval[(i1*num_states*num_states+i2*num_states+i3)*numz+i1*num_states*num_states+i2*num_states+i3]++;

		matbeg[i1*num_states*num_states+i2*num_states+i3]=numz*(i1*num_states*num_states+i2*num_states+i3);
		matcnt[i1*num_states*num_states+i2*num_states+i3]=numz;
		sense [i1*num_states*num_states+i2*num_states+i3]='E';
		obj[i1*num_states*num_states+i2*num_states+i3]=0;
		lb[i1*num_states*num_states+i2*num_states+i3]=0;
		ub[i1*num_states*num_states+i2*num_states+i3]=CPX_INFBOUND;
	}
	
	obj[*obj_ind]=1;
	

	for (int i=0;i<num_patients;i++){
		status = CPXcopylp (env_slp, slp[i], numz, numz, CPX_MAX, obj, rhs_contin[i],
						 sense, matbeg, matcnt, matind, matval, lb, ub, NULL);
		if ( status ) {
			fprintf (stderr, "CPXcopylp failed to populate subproblem 2.\n");
			goto TERMINATE;
		}
	
		status=CPXlpopt(env_slp,slp[i]);
		if ( status ) {
			fprintf (stderr, "CPXlpopt failed to optimizae subproblem 1.\n");
			goto TERMINATE;
		}

		status = CPXgetx (env_slp, slp[i], d[i], 0, numz-1);
		if ( status ) {
			fprintf (stderr, "CPXgetx failed to obtain d values of subproblem 1.\n");
			goto TERMINATE;
		}
	}

	
#pragma endregion	


#pragma region calculating v in the following	

	for(int k=0;k<num_patients;k++){

		slp_v[k] = CPXcloneprob (env_slp, slp[k], &status);
		if ( status ) {
		fprintf (stderr, "CPXcopylp failed to creates a copy of slp problem.\n");
		goto TERMINATE;
		}	
	
		for(int j=0;j<numz;j++)		
		{
			rhs[j]=rhs_stop[k][j];
			sense [j]='G';
			obj[j]=1;
			indices[j]=j;
			matbeg[j]=j;
		}

		for (int j=0;j<numz;j++)	
		if(d[0][j]>rhs_stop[0][j]+EPSILON || d[1][j]>rhs_stop[1][j]+EPSILON ||d[2][j]>rhs_stop[2][j]+EPSILON)
			rhs[j]=0;

		status = CPXchgobj (env_slp, slp_v[k], numz, indices, obj);
		if ( status ) {
		fprintf (stderr, "CPXcopylp failed to change obj fun.\n");
		goto TERMINATE;
		}

		status = CPXaddrows (env_slp, slp_v[k], 0, numz,numz, rhs,sense, matbeg, indices,
							obj, NULL, NULL);
		if ( status ) {
		fprintf (stderr, "CPXcopylp failed to add columns.\n");
		goto TERMINATE;
		}
	
		status = CPXchgsense (env_slp, slp_v[k], numz, indices, sense);
		CPXchgobjsen (env_slp, slp_v[k], CPX_MIN);

		status = CPXlpopt (env_slp, slp_v[k]);
		if ( status ) {
		fprintf (stderr, "CPXcopylp failed to optimize copy_slp_v.\n");
		goto TERMINATE;
		}	
 
		status = CPXgetx (env_slp, slp_v[k], v[k], 0, numz-1);
		if ( status ) {
		fprintf (stderr, "CPXcopylp failed to access solution values copy_slp_v.\n");
		goto TERMINATE;
		}

	}

#pragma endregion

#pragma region start creating master problems with column elimination for master problem with respect to d_i(S)

	char     *ctype = NULL;
	ctype= new char [numz];
	for(int i=0;i<numz;i++)
		ctype[i]='B';		
	for (int j1=0;j1<num_states;j1++)
	for (int j2=0;j2<num_states;j2++)
	for (int j3=0;j3<num_states;j3++)
	{
		if (j1!=num_states-1&& j2!=num_states-1 && j3!=num_states-1)
			ub[j1*num_states*num_states+j2*num_states+j3]=1;
		else	
			ub[j1*num_states*num_states+j2*num_states+j3]=0;
	}

		
	for (int i=0;i<numz;i++)	
	if(d[0][i]>rhs_stop[0][i]+EPSILON || d[1][i]>rhs_stop[1][i]+EPSILON ||d[2][i]>rhs_stop[2][i]+EPSILON)		
		ub[i]=0;
		
	status = CPXnewcols (env, lp, numz, NULL, NULL, ub, ctype, NULL);

	if ( status ) {
	fprintf (stderr, "Failed to create LP.\n");
	goto TERMINATE;
	}	
	double *z=new double [num_patients];
	double *upperbound=new double [num_patients];

	
	for (int i=0;i<num_patients;i++)
	{
		z[i]=1;
		upperbound[i]=v[i][*obj_ind];

	}



	status = CPXaddcols (env, lp, num_patients, 0, z, NULL, NULL, NULL, NULL, upperbound, NULL);
	if ( status ) {
	fprintf (stderr, "Failed to create LP.\n");
	goto TERMINATE;
	}	

	CPXchgobjsen (env, lp, CPX_MAX);

	delete[]z;
	delete[]ctype;
	delete[]upperbound;					
	//delete[] colname;	
	
#pragma endregion

TERMINATE:
#pragma region	deallocating dynamic memories
	delete [] matval;
	delete [] matind;
	delete [] matbeg;
	delete [] matcnt;
	delete [] sense;
	delete [] obj;
	delete [] rhs;
	delete [] lb;
	delete [] ub;
	delete []indices;
#pragma endregion				
		return 0;
	}
	


int make_master2
		(CPXCENVptr env,CPXCENVptr env_slp,CPXLPptr lp, CPXLPptr* slp, CPXLPptr* slp_v, double ***tran_prob, double **rhs_contin,double **rhs_stop, double **v,double **d, int *obj_ind)	//making master problem	
	{		
#pragma region  start creating sub_problems for the first patient and calculating d_1(s)
	double *matval= new double [numz*numz];
	int *matind= new int [numz*numz];
	int *matbeg= new int [numz];
	int *matcnt= new int [numz];
	char *sense= new char [numz];
	double *obj= new double [numz];
	double *rhs= new double [numz];
	double *lb= new double [numz];
	double *ub= new double [numz];
	int *indices = new int[numz];

	for(int i1=0;i1<num_states;i1++)
	for(int i2=0;i2<num_states;i2++)
	{
		for(int j1=0;j1<num_states;j1++)
		for(int j2=0;j2<num_states;j2++)
		{
			matval[(i1*num_states+i2)*numz+j1*num_states+j2]=(-lambda)*tran_prob[0][j1][i1]*tran_prob[1][j2][i2];
			matind[(i1*num_states+i2)*numz+j1*num_states+j2]=j1*num_states+j2;

		}
		matval[(i1*num_states+i2)*numz+i1*num_states+i2]++;

		matbeg[i1*num_states+i2]=numz*(i1*num_states+i2);
		matcnt[i1*num_states+i2]=numz;
		sense [i1*num_states+i2]='E';
		obj[i1*num_states+i2]=0;
		lb[i1*num_states+i2]=0;
		ub[i1*num_states+i2]=CPX_INFBOUND;
		indices[i1*num_states + i2] = i1*num_states + i2;

	}
	
	obj[*obj_ind]=1;

	status = CPXcopylp (env_slp, slp[0], numz, numz, CPX_MAX, obj, rhs_contin[0],
                     sense, matbeg, matcnt, matind, matval, lb,
                     ub, NULL);
	if ( status ) {
		fprintf (stderr, "CPXcopylp failed to populate subproblem 1.\n");
		goto TERMINATE;
	}
	status=CPXlpopt(env_slp,slp[0]);
	if ( status ) {
		fprintf (stderr, "CPXlpopt failed to optimizae subproblem 1.\n");
		goto TERMINATE;
	}

	status = CPXgetx (env_slp, slp[0], d[0], 0, numz-1);
	if ( status ) {
		fprintf (stderr, "CPXgetx failed to obtain d values of subproblem 1.\n");
		goto TERMINATE;
	}
#pragma endregion
	


#pragma region start creating sub_problems for the second patient and calculating d_2(s)
	status = CPXcopylp (env_slp, slp[1], numz, numz, CPX_MAX, obj, rhs_contin[1],
                     sense, matbeg, matcnt, matind, matval, lb,
                     ub, NULL);
	if ( status ) {
		fprintf (stderr, "CPXcopylp failed to populate subproblem 2.\n");
		goto TERMINATE;
	}
	
	status=CPXlpopt(env_slp,slp[1]);
	if ( status ) {
		fprintf (stderr, "CPXlpopt failed to optimizae subproblem 1.\n");
		goto TERMINATE;
	}

	status = CPXgetx (env_slp, slp[1], d[1], 0, numz-1);
	if ( status ) {
		fprintf (stderr, "CPXgetx failed to obtain d values of subproblem 1.\n");
		goto TERMINATE;
	}

#pragma endregion	


#pragma region calculating v in the following

	slp_v[0] = CPXcloneprob (env_slp, slp[0], &status);
	if ( status ) {
	fprintf (stderr, "CPXcopylp failed to creates a copy of slp problem.\n");
	goto TERMINATE;
	}
	

	for(int i1=0;i1<num_states;i1++)
	for(int i2=0;i2<num_states;i2++)
	{
		rhs[i1*num_states+i2]=rhs_stop[0][i1*num_states+i2];

		sense [i1*num_states+i2]='G';
		obj[i1*num_states+i2]=1;
		matbeg[i1*num_states+i2]=i1*num_states+i2;
	}

	for (int i=0;i<num_states;i++)
	for (int j=0;j<num_states;j++)
	{
		if (i==num_states-1|| j==num_states-1)
			rhs[i*num_states+j]=0;
	}
	for (int i=0;i<numz;i++)	
	if(d[0][i]>rhs_stop[0][i]+EPSILON || d[1][i]>rhs_stop[1][i]+EPSILON)
		rhs[i]=0; 

	status = CPXchgobj (env_slp, slp_v[0], numz, indices, obj);
	if ( status ) {
	fprintf (stderr, "CPXcopylp failed to change obj fun.\n");
	goto TERMINATE;
	}

	status = CPXaddrows (env_slp, slp_v[0], 0, numz,numz, rhs,sense, matbeg, indices,
                        obj, NULL, NULL);
	if ( status ) {
	fprintf (stderr, "CPXcopylp failed to add columns.\n");
	goto TERMINATE;
	}
	
	status = CPXchgsense (env_slp, slp_v[0], numz, indices, sense);
	CPXchgobjsen (env_slp, slp_v[0], CPX_MIN);

	status = CPXlpopt (env_slp, slp_v[0]);
	if ( status ) {
	fprintf (stderr, "CPXcopylp failed to optimize copy_slp_v.\n");
	goto TERMINATE;
	}	
 
	status = CPXgetx (env_slp, slp_v[0], v[0], 0, numz-1);
	if ( status ) {
	fprintf (stderr, "CPXcopylp failed to access solution values copy_slp_v.\n");
	goto TERMINATE;
	}

	//////////////calculatin of V2/////////
	slp_v[1] = CPXcloneprob (env_slp, slp_v[0], &status);
	if ( status ) {
	fprintf (stderr, "CPXcopylp failed to creates a copy of slp problem.\n");
	goto TERMINATE;
	}
	status = CPXchgrhs (env_slp, slp_v[1], numz, indices, rhs_contin[1]);
	if ( status ) {
	fprintf (stderr, "Failed to change righthand-side of the subproblem_V .\n");
	goto TERMINATE;
		}
	
	for (int j1=0; j1<num_states;j1++)
	for (int j2=0; j2<num_states;j2++)
	{
		rhs[j1*num_states+j2]=rhs_stop[1][j1*num_states+j2];
		indices[j1*num_states+j2]=numz+j1*num_states+j2;
	}
	
	for (int i=0;i<num_states;i++)
	for (int j=0;j<num_states;j++)
	{
		if (i==num_states-1|| j==num_states-1)
			rhs[i*num_states+j]=0;
	}

	for (int i=0;i<numz;i++)	
	if(d[0][i]>rhs_stop[0][i]+EPSILON || d[1][i]>rhs_stop[1][i]+EPSILON)
		rhs[i]=0; 

	status = CPXchgrhs (env_slp, slp_v[1], numz, indices, rhs);
	if ( status ) {
	fprintf (stderr, "Failed to change righthand-side of the subproblem_V .\n");
	goto TERMINATE;
		}


	status = CPXlpopt (env_slp, slp_v[1]);
	if ( status ) {
	fprintf (stderr, "CPXcopylp failed to optimize copy_slp_v.\n");
	goto TERMINATE;
	}

	status = CPXgetx (env_slp, slp_v[1], v[1], 0, numz-1);
	if ( status ) {
	fprintf (stderr, "CPXcopylp failed to access solution values copy_slp_v.\n");
	goto TERMINATE;
	}
	
#pragma endregion

#pragma region start creating the master problem

	char     *ctype = NULL;
	ctype= new char [numz];
	for(int i=0;i<numz;i++)
		ctype[i]='B';		
	for (int i=0;i<num_states;i++)
	for (int j=0;j<num_states;j++)
	{
		if (i!=num_states-1&& j!=num_states-1)
			ub[i*num_states+j]=1;
		else	
			ub[i*num_states+j]=0;
	}

		
	for (int i=0;i<numz;i++)	
	if(d[0][i]>rhs_stop[0][i]+EPSILON || d[1][i]>rhs_stop[1][i]+EPSILON)		
		ub[i]=0;
		

	status = CPXnewcols (env, lp, numz, NULL, NULL, ub, ctype, NULL);

	if ( status ) {
	fprintf (stderr, "Failed to create LP.\n");
	goto TERMINATE;
	}	
	double *z=new double [num_patients];
	double *upperbound=new double [num_patients];
	char *colname[num_patients];
	
	for (int i=0;i<num_patients;i++)
	{
		z[i]=1;
		upperbound[i]=v[i][*obj_ind];

	}

	colname[0]="z0";
	colname[1]="z1";

	status = CPXaddcols (env, lp, num_patients, 0, z, NULL, NULL, NULL, NULL, upperbound, colname);
	if ( status ) {
	fprintf (stderr, "Failed to create LP.\n");
	goto TERMINATE;
	}	

	CPXchgobjsen (env, lp, CPX_MAX);

	delete[]z;
	delete[]ctype;
	delete[]upperbound;					
	//delete[] colname;	
	
	
#pragma endregion


TERMINATE:

#pragma region	deallocating dynamic memories
	delete [] matval;
	delete [] matind;
	delete [] matbeg;
	delete [] matcnt;
	delete [] sense;
	delete [] obj;
	delete [] rhs;
	delete [] lb;
	delete [] ub;
	delete []indices;
#pragma endregion				
		return 0;
	}


	
int make_slpvt2
	(CPXCENVptr env_slp, CPXLPptr slp_t, CPXLPptr slp_vt, double ***tran_prob, double **rhs_contin,double **rhs_stop, double *vt,double **d, int *obj_ind){ //making slp_vt and slp_t
#pragma region	making slp_vt
	double *matval= new double [numz*numz];
	int *matind= new int [numz*numz];
	int *matbeg= new int [numz];
	int *matcnt= new int [numz];
	char *sense= new char [numz];
	double *obj= new double [numz];
	double *rhs= new double [numz];
	double *lb= new double [numz];
	double *ub= new double [numz];
	int * indices= new int[numz];

	for(int i1=0;i1<num_states;i1++)
	for(int i2=0;i2<num_states;i2++)
	{
		for(int j1=0;j1<num_states;j1++)
		for(int j2=0;j2<num_states;j2++)
		{
			matval[(i1*num_states+i2)*numz+j1*num_states+j2]=(-lambda)*tran_prob[0][j1][i1]*tran_prob[1][j2][i2];
			matind[(i1*num_states+i2)*numz+j1*num_states+j2]=j1*num_states+j2;

		}
		matval[(i1*num_states+i2)*numz+i1*num_states+i2]++;

		matbeg[i1*num_states+i2]=numz*(i1*num_states+i2);
		matcnt[i1*num_states+i2]=numz;
		rhs[i1*num_states+i2]=rhs_contin[0][i1*num_states+i2]+rhs_contin[1][i1*num_states+i2];
		sense [i1*num_states+i2]='G';
		obj[i1*num_states+i2]=1;
		lb[i1*num_states+i2]=0;
		ub[i1*num_states+i2]=CPX_INFBOUND;
	}
		

	status = CPXcopylp (env_slp, slp_vt, numz, numz, CPX_MIN, obj, rhs,sense, matbeg, matcnt, matind, matval, lb,ub, NULL);
	if ( status ) 	goto TERMINATE;
	
	
	for(int i1=0;i1<num_states;i1++)
	for(int i2=0;i2<num_states;i2++)
	{
		rhs[i1*num_states+i2]=rhs_stop[0][i1*num_states+i2]+rhs_stop[1][i1*num_states+i2];

		sense [i1*num_states+i2]='G';
		obj[i1*num_states+i2]=1;
		indices[i1*num_states+i2]=i1*num_states+i2;
		matbeg[i1*num_states+i2]=i1*num_states+i2;
	}

	for (int i=0;i<numz;i++)	
	if(d[0][i]>rhs_stop[0][i]+EPSILON || d[1][i]>rhs_stop[1][i]+EPSILON)
		rhs[i]=0; 
	status = CPXaddrows (env_slp, slp_vt, 0, numz,numz, rhs,sense, matbeg, indices,obj, NULL, NULL);
	status = CPXlpopt (env_slp, slp_vt);	
	status = CPXgetx (env_slp, slp_vt, vt, 0, numz-1);

#pragma endregion

	for(int i1=0;i1<num_states;i1++)
	for(int i2=0;i2<num_states;i2++)
	{		
		matbeg[i1*num_states+i2]=numz*(i1*num_states+i2);
		matcnt[i1*num_states+i2]=numz;
		rhs[i1*num_states+i2]=rhs_contin[0][i1*num_states+i2];
		sense [i1*num_states+i2]='E';
		obj[i1*num_states+i2]=0;
		lb[i1*num_states+i2]=0;
		ub[i1*num_states+i2]=CPX_INFBOUND;
	}
	obj[*obj_ind]=1;	

	status = CPXcopylp (env_slp, slp_t, numz, numz, CPX_MAX, obj, rhs,sense, matbeg, matcnt, matind, matval, lb,ub, NULL);
	if ( status ) 	goto TERMINATE;

	for(int i=0;i<numz;i++)	
		obj[i]=0;
	obj[*obj_ind]=1;
	status = CPXnewcols (env_slp, slp_t, numz, obj, NULL, NULL, NULL, NULL);
	if ( status ) goto TERMINATE;	

	for(int i1=0;i1<num_states;i1++)
	for(int i2=0;i2<num_states;i2++)
	{
		for(int j1=0;j1<num_states;j1++)
		for(int j2=0;j2<num_states;j2++)
		{
			matval[(i1*num_states+i2)*numz+j1*num_states+j2]=(-lambda)*tran_prob[0][i1][j1]*tran_prob[1][i2][j2];
			matind[(i1*num_states+i2)*numz+j1*num_states+j2]=numz+j1*num_states+j2;

		}
		matval[(i1*num_states+i2)*numz+i1*num_states+i2]++;

		matbeg[i1*num_states+i2]=numz*(i1*num_states+i2);
		sense [i1*num_states+i2]='E';				
		rhs [i1*num_states+i2]=rhs_contin[1][i1*num_states+i2];				
	}

	status = CPXaddrows (env_slp, slp_t, 0, numz,numz*numz, rhs,sense, matbeg, matind,matval, NULL, NULL);
	if ( status ) goto TERMINATE;

	double *val= new double [num_patients*numz];
	int *ind= new int [num_patients*numz];
	double cutrhs=100000;
	int beg=0;
	for(int j=0;j<num_patients*numz;j++)
		ind[j]=j;
	for(int j=0;j<numz;j++){
		for(int k=0;k<2*numz;k++)
			val[k]=0;
		val[j]=1, val[numz+j]=1;
		status=CPXaddrows(env_slp, slp_t, 0, 1, num_patients*numz, &cutrhs, "L", &beg, ind, val, NULL, NULL);
	}

	delete[] val;
	delete[] ind;
	

TERMINATE:
	delete [] matval;
	delete [] matind;
	delete [] matbeg;
	delete [] matcnt;
	delete [] sense;
	delete [] obj;
	delete [] rhs;
	delete [] lb;
	delete [] ub;
	delete []indices;

	return 0;
	}



	int make_slpvt3(CPXCENVptr env_slp, CPXLPptr slp_t, CPXLPptr slp_vt, double ***tran_prob, double **rhs_contin,double **rhs_stop, double *vt,double **d, int *obj_ind){
#pragma region	making slp_vt
	double *matval= new double [numz*numz];
	int *matind= new int [numz*numz];
	int *matbeg= new int [numz];
	int *matcnt= new int [numz];
	char *sense= new char [numz];
	double *obj= new double [numz];
	double *rhs= new double [numz];
	double *lb= new double [numz];
	double *ub= new double [numz];
	int *indices= new int[numz];

	for(int i1=0;i1<num_states;i1++)
	for(int i2=0;i2<num_states;i2++)
	for(int i3=0;i3<num_states;i3++)
	{
		for(int j1=0;j1<num_states;j1++)
		for(int j2=0;j2<num_states;j2++)
		for(int j3=0;j3<num_states;j3++)
		{
			matval[(i1*num_states*num_states+i2*num_states+i3)*numz+j1*num_states*num_states+j2*num_states+j3]=(-lambda)*tran_prob[0][j1][i1]*tran_prob[1][j2][i2]*tran_prob[2][j3][i3];
			matind[(i1*num_states*num_states+i2*num_states+i3)*numz+j1*num_states*num_states+j2*num_states+j3]=j1*num_states*num_states+j2*num_states+j3;

		}
		matval[(i1*num_states*num_states+i2*num_states+i3)*numz+i1*num_states*num_states+i2*num_states+i3]++;

		matbeg[i1*num_states*num_states+i2*num_states+i3]=numz*(i1*num_states*num_states+i2*num_states+i3);
		matcnt[i1*num_states*num_states+i2*num_states+i3]=numz;
		rhs[i1*num_states*num_states+i2*num_states+i3]=rhs_contin[0][i1*num_states*num_states+i2*num_states+i3]+rhs_contin[1][i1*num_states*num_states+i2*num_states+i3]+rhs_contin[2][i1*num_states*num_states+i2*num_states+i3];
		sense [i1*num_states*num_states+i2*num_states+i3]='G';
		obj[i1*num_states*num_states+i2*num_states+i3]=1;
		lb[i1*num_states*num_states+i2*num_states+i3]=0;
		ub[i1*num_states*num_states+i2*num_states+i3]=CPX_INFBOUND;
	}
		

	status = CPXcopylp (env_slp, slp_vt, numz, numz, CPX_MIN, obj, rhs,sense, matbeg, matcnt, matind, matval, lb,ub, NULL);
	if ( status ) 	goto TERMINATE;
	
	
	for(int i=0;i<numz;i++)	
	{
		rhs[i]=rhs_stop[0][i]+rhs_stop[1][i]+rhs_stop[2][i];
		sense [i]='G';
		obj[i]=1;
		indices[i]=i;
		matbeg[i]=i;
	}

	for (int i=0;i<numz;i++)	
	if(d[0][i]>rhs_stop[0][i]+EPSILON || d[1][i]>rhs_stop[1][i]+EPSILON ||d[2][i]>rhs_stop[2][i]+EPSILON)
		rhs[i]=0; 
	status = CPXaddrows (env_slp, slp_vt, 0, numz,numz, rhs,sense, matbeg, indices,obj, NULL, NULL);
	status = CPXlpopt (env_slp, slp_vt);	
	status = CPXgetx (env_slp, slp_vt, vt, 0, numz-1);
	 
#pragma endregion

	for(int i1=0;i1<num_states;i1++)
	for(int i2=0;i2<num_states;i2++)
	for(int i3=0;i3<num_states;i3++)
	{		
		matbeg[i1*num_states*num_states+i2*num_states+i3]=numz*(i1*num_states*num_states+i2*num_states+i3);
		matcnt[i1*num_states*num_states+i2*num_states+i3]=numz;
		rhs[i1*num_states*num_states+i2*num_states+i3]=rhs_contin[0][i1*num_states*num_states+i2*num_states+i3];
		sense [i1*num_states*num_states+i2*num_states+i3]='E';
		obj[i1*num_states*num_states+i2*num_states+i3]=0;
		lb[i1*num_states*num_states+i2*num_states+i3]=0;
		ub[i1*num_states*num_states+i2*num_states+i3]=CPX_INFBOUND;
	}
	obj[*obj_ind]=1;	

	status = CPXcopylp (env_slp, slp_t, numz, numz, CPX_MAX, obj, rhs,sense, matbeg, matcnt, matind, matval, lb,ub, NULL);
	if ( status ) 	goto TERMINATE;

	for(int k=1;k<num_patients;k++){

		for(int i=0;i<numz;i++)	
			obj[i]=0;
		obj[*obj_ind]=1;
		status = CPXnewcols (env_slp, slp_t, numz, obj, NULL, NULL, NULL, NULL);
		if ( status ) goto TERMINATE;	

		for(int i1=0;i1<num_states;i1++)
		for(int i2=0;i2<num_states;i2++)
		for(int i3=0;i3<num_states;i3++)
		{
			for(int j1=0;j1<num_states;j1++)
			for(int j2=0;j2<num_states;j2++)
			for(int j3=0;j3<num_states;j3++)
			{
				matval[(i1*num_states*num_states+i2*num_states+i3)*numz+j1*num_states*num_states+j2*num_states+j3]=(-lambda)*tran_prob[0][i1][j1]*tran_prob[1][i2][j2]*tran_prob[2][i3][j3];
				matind[(i1*num_states*num_states+i2*num_states+i3)*numz+j1*num_states*num_states+j2*num_states+j3]=numz+j1*num_states*num_states+j2*num_states+j3;

			}
			matval[(i1*num_states*num_states+i2*num_states+i3)*numz+i1*num_states*num_states+i2*num_states+i3]++;

			matbeg[i1*num_states*num_states+i2*num_states+i3]=numz*(i1*num_states*num_states+i2*num_states+i3);
			sense [i1*num_states*num_states+i2*num_states+i3]='E';				
			rhs [i1*num_states*num_states+i2*num_states+i3]=rhs_contin[k][i1*num_states*num_states+i2*num_states+i3];				
		}

		status = CPXaddrows (env_slp, slp_t, 0, numz,numz*numz, rhs,sense, matbeg, matind,matval, NULL, NULL);
		if ( status ) goto TERMINATE;

	}


	double *val= new double [num_patients*numz];
	int *ind= new int [num_patients*numz];
	double cutrhs=100000;
	int beg=0;
	for(int j=0;j<num_patients*numz;j++)
		ind[j]=j;
	for(int j=0;j<numz;j++){
		for(int k=0;k<num_patients*numz;k++)
			val[k]=0;
		val[j]=1, val[numz+j]=1; val[2*numz+j]=1;
		status=CPXaddrows(env_slp, slp_t, 0, 1, num_patients*numz, &cutrhs, "L", &beg, ind, val, NULL, NULL);
	}

	delete[] val;
	delete[] ind;

TERMINATE:
	delete [] matval;
	delete [] matind;
	delete [] matbeg;
	delete [] matcnt;
	delete [] sense;
	delete [] obj;
	delete [] rhs;
	delete [] lb;
	delete [] ub;
	delete []indices;

	return 0;
	}



	int parameters_calulator_3
		(CPXCENVptr env_slp, CPXLPptr* slp, CPXLPptr* slp_v, double ***tran_prob, double **rhs_contin,double **rhs_stop, double **v,double **d, int *obj_ind)		//calculating parameters of 3 player instances	
	{		
#pragma region  start creating sub_problems for the  patients and calculating d_i(s)
	double *matval= new double [numz*numz];
	int *matind= new int [numz*numz];
	int *matbeg= new int [numz];
	int *matcnt= new int [numz];
	char *sense= new char [numz];
	double *obj= new double [numz];
	double *rhs= new double [numz];
	double *lb= new double [numz];
	double *ub= new double [numz];
	int *indices = new int[numz];

	for(int i1=0;i1<num_states;i1++)
	for(int i2=0;i2<num_states;i2++)
	for(int i3=0;i3<num_states;i3++)
	{
		for(int j1=0;j1<num_states;j1++)
		for(int j2=0;j2<num_states;j2++)
		for(int j3=0;j3<num_states;j3++)
		{
			matval[(i1*num_states*num_states+i2*num_states+i3)*numz+j1*num_states*num_states+j2*num_states+j3]=(-lambda)*tran_prob[0][j1][i1]*tran_prob[1][j2][i2]*tran_prob[2][j3][i3];
			matind[(i1*num_states*num_states+i2*num_states+i3)*numz+j1*num_states*num_states+j2*num_states+j3]=j1*num_states*num_states+j2*num_states+j3;

		}
		matval[(i1*num_states*num_states+i2*num_states+i3)*numz+i1*num_states*num_states+i2*num_states+i3]++;

		matbeg[i1*num_states*num_states+i2*num_states+i3]=numz*(i1*num_states*num_states+i2*num_states+i3);
		matcnt[i1*num_states*num_states+i2*num_states+i3]=numz;
		sense [i1*num_states*num_states+i2*num_states+i3]='E';
		obj[i1*num_states*num_states+i2*num_states+i3]=0;
		lb[i1*num_states*num_states+i2*num_states+i3]=0;
		ub[i1*num_states*num_states+i2*num_states+i3]=CPX_INFBOUND;
	}
	
	obj[*obj_ind]=1;
	

	for (int i=0;i<num_patients;i++){
		status = CPXcopylp (env_slp, slp[i], numz, numz, CPX_MAX, obj, rhs_contin[i],
						 sense, matbeg, matcnt, matind, matval, lb, ub, NULL);
		if ( status ) {
			fprintf (stderr, "CPXcopylp failed to populate subproblem 2.\n");
			goto TERMINATE;
		}
	
		status=CPXlpopt(env_slp,slp[i]);
		if ( status ) {
			fprintf (stderr, "CPXlpopt failed to optimizae subproblem 1.\n");
			goto TERMINATE;
		}

		status = CPXgetx (env_slp, slp[i], d[i], 0, numz-1);
		if ( status ) {
			fprintf (stderr, "CPXgetx failed to obtain d values of subproblem 1.\n");
			goto TERMINATE;
		}
	}

	
#pragma endregion	


#pragma region calculating v in the following

	for(int k=0;k<num_patients;k++){

		slp_v[k] = CPXcloneprob (env_slp, slp[k], &status);
		if ( status ) {
		fprintf (stderr, "CPXcopylp failed to creates a copy of slp problem.\n");
		goto TERMINATE;
		}	
	
		for(int j=0;j<numz;j++)		
		{
			rhs[j]=rhs_stop[k][j];
			sense [j]='G';
			obj[j]=1;
			indices[j]=j;
			matbeg[j]=j;
		}
		
		status = CPXchgobj (env_slp, slp_v[k], numz, indices, obj);
		if ( status ) {
		fprintf (stderr, "CPXcopylp failed to change obj fun.\n");
		goto TERMINATE;
		}

		status = CPXaddrows (env_slp, slp_v[k], 0, numz,numz, rhs,sense, matbeg, indices,
							obj, NULL, NULL);
		if ( status ) {
		fprintf (stderr, "CPXcopylp failed to add columns.\n");
		goto TERMINATE;
		}
	
		status = CPXchgsense (env_slp, slp_v[k], numz, indices, sense);
		CPXchgobjsen (env_slp, slp_v[k], CPX_MIN);

		status = CPXlpopt (env_slp, slp_v[k]);
		if ( status ) {
		fprintf (stderr, "CPXcopylp failed to optimize copy_slp_v.\n");
		goto TERMINATE;
		}	
 
		status = CPXgetx (env_slp, slp_v[k], v[k], 0, numz-1);
		if ( status ) {
		fprintf (stderr, "CPXcopylp failed to access solution values copy_slp_v.\n");
		goto TERMINATE;
		}

	}

#pragma endregion
	
TERMINATE:
#pragma region	deallocating dynamic memories
	delete [] matval;
	delete [] matind;
	delete [] matbeg;
	delete [] matcnt;
	delete [] sense;
	delete [] obj;
	delete [] rhs;
	delete [] lb;
	delete [] ub;
	delete []indices;
#pragma endregion				
		return 0;
	}


	
	int parameters_calulator_2
		(CPXCENVptr env_slp, CPXLPptr* slp, CPXLPptr* slp_v, double ***tran_prob, double **rhs_contin,double **rhs_stop, double **v,double **d, int *obj_ind)		//calculating parameters of 2 player instances	
	{		
#pragma region  start creating sub_problems for the first patient and calculating d_1(s)
	double *matval= new double [numz*numz];
	int *matind= new int [numz*numz];
	int *matbeg= new int [numz];
	int *matcnt= new int [numz];
	char *sense= new char [numz];
	double *obj= new double [numz];
	double *rhs= new double [numz];
	double *lb= new double [numz];
	double *ub= new double [numz];
	int *indices = new int[numz];

	for(int i1=0;i1<num_states;i1++)
	for(int i2=0;i2<num_states;i2++)
	{
		for(int j1=0;j1<num_states;j1++)
		for(int j2=0;j2<num_states;j2++)
		{
			matval[(i1*num_states+i2)*numz+j1*num_states+j2]=(-lambda)*tran_prob[0][j1][i1]*tran_prob[1][j2][i2];
			matind[(i1*num_states+i2)*numz+j1*num_states+j2]=j1*num_states+j2;

		}
		matval[(i1*num_states+i2)*numz+i1*num_states+i2]++;

		matbeg[i1*num_states+i2]=numz*(i1*num_states+i2);
		matcnt[i1*num_states+i2]=numz;
		sense [i1*num_states+i2]='E';
		obj[i1*num_states+i2]=0;
		lb[i1*num_states+i2]=0;
		ub[i1*num_states+i2]=CPX_INFBOUND;
	}
	
	obj[*obj_ind]=1;

	status = CPXcopylp (env_slp, slp[0], numz, numz, CPX_MAX, obj, rhs_contin[0],
                     sense, matbeg, matcnt, matind, matval, lb,
                     ub, NULL);
	if ( status ) {
		fprintf (stderr, "CPXcopylp failed to populate subproblem 1.\n");
		goto TERMINATE;
	}
	status=CPXlpopt(env_slp,slp[0]);
	if ( status ) {
		fprintf (stderr, "CPXlpopt failed to optimizae subproblem 1.\n");
		goto TERMINATE;
	}

	status = CPXgetx (env_slp, slp[0], d[0], 0, numz-1);
	if ( status ) {
		fprintf (stderr, "CPXgetx failed to obtain d values of subproblem 1.\n");
		goto TERMINATE;
	}
#pragma endregion
	


#pragma region start creating sub_problems for the second patient and calculating d_2(s)
	status = CPXcopylp (env_slp, slp[1], numz, numz, CPX_MAX, obj, rhs_contin[1],
                     sense, matbeg, matcnt, matind, matval, lb,
                     ub, NULL);
	if ( status ) {
		fprintf (stderr, "CPXcopylp failed to populate subproblem 2.\n");
		goto TERMINATE;
	}
	
	status=CPXlpopt(env_slp,slp[1]);
	if ( status ) {
		fprintf (stderr, "CPXlpopt failed to optimizae subproblem 1.\n");
		goto TERMINATE;
	}

	status = CPXgetx (env_slp, slp[1], d[1], 0, numz-1);
	if ( status ) {
		fprintf (stderr, "CPXgetx failed to obtain d values of subproblem 1.\n");
		goto TERMINATE;
	}

#pragma endregion	


#pragma region calculating v in the following

	slp_v[0] = CPXcloneprob (env_slp, slp[0], &status);
	if ( status ) {
	fprintf (stderr, "CPXcopylp failed to creates a copy of slp problem.\n");
	goto TERMINATE;
	}
	


	for(int i1=0;i1<num_states;i1++)
	for(int i2=0;i2<num_states;i2++)
	{
		rhs[i1*num_states+i2]=rhs_stop[0][i1*num_states+i2];

		sense [i1*num_states+i2]='G';
		obj[i1*num_states+i2]=1;
		indices[i1*num_states+i2]=i1*num_states+i2;
		matbeg[i1*num_states+i2]=i1*num_states+i2;
	}


	status = CPXchgobj (env_slp, slp_v[0], numz, indices, obj);
	if ( status ) {
	fprintf (stderr, "CPXcopylp failed to change obj fun.\n");
	goto TERMINATE;
	}

	status = CPXaddrows (env_slp, slp_v[0], 0, numz,numz, rhs,sense, matbeg, indices,
                        obj, NULL, NULL);
	if ( status ) {
	fprintf (stderr, "CPXcopylp failed to add columns.\n");
	goto TERMINATE;
	}
	
	status = CPXchgsense (env_slp, slp_v[0], numz, indices, sense);
	CPXchgobjsen (env_slp, slp_v[0], CPX_MIN);

	status = CPXlpopt (env_slp, slp_v[0]);
	if ( status ) {
	fprintf (stderr, "CPXcopylp failed to optimize copy_slp_v.\n");
	goto TERMINATE;
	}	
 
	status = CPXgetx (env_slp, slp_v[0], v[0], 0, numz-1);
	if ( status ) {
	fprintf (stderr, "CPXcopylp failed to access solution values copy_slp_v.\n");
	goto TERMINATE;
	}

	//////////////calculatin of V2/////////
	slp_v[1] = CPXcloneprob (env_slp, slp_v[0], &status);
	if ( status ) {
	fprintf (stderr, "CPXcopylp failed to creates a copy of slp problem.\n");
	goto TERMINATE;
	}
	status = CPXchgrhs (env_slp, slp_v[1], numz, indices, rhs_contin[1]);
	if ( status ) {
	fprintf (stderr, "Failed to change righthand-side of the subproblem_V .\n");
	goto TERMINATE;
		}

	for (int j1=0; j1<num_states;j1++)
	for (int j2=0; j2<num_states;j2++)
	{
		rhs[j1*num_states+j2]=rhs_stop[1][j1*num_states+j2];
		indices[j1*num_states+j2]=numz+j1*num_states+j2;
	}
	
	status = CPXchgrhs (env_slp, slp_v[1], numz, indices, rhs);
	if ( status ) {
	fprintf (stderr, "Failed to change righthand-side of the subproblem_V .\n");
	goto TERMINATE;
		}


	status = CPXlpopt (env_slp, slp_v[1]);
	if ( status ) {
	fprintf (stderr, "CPXcopylp failed to optimize copy_slp_v.\n");
	goto TERMINATE;
	}

	status = CPXgetx (env_slp, slp_v[1], v[1], 0, numz-1);
	if ( status ) {
	fprintf (stderr, "CPXcopylp failed to access solution values copy_slp_v.\n");
	goto TERMINATE;
	}
	
#pragma endregion

	
TERMINATE:
#pragma region	deallocating dynamic memories
	delete [] matval;
	delete [] matind;
	delete [] matbeg;
	delete [] matcnt;
	delete [] sense;
	delete [] obj;
	delete [] rhs;
	delete [] lb;
	delete [] ub;
	delete []indices;
#pragma endregion				
		return 0;
	}


	int make_extensive_2(CPXCENVptr env_slp, CPXLPptr exten, double ***tran_prob, double **rhs_contin,double **rhs_stop, double **v,double **d, int *obj_ind){//constructing the extensive form problem for two player instances
	double *matval= new double [numz*numz];
	int *matind= new int [numz*numz];
	int *matbeg= new int [numz];
	int *matcnt= new int [numz];
	char *sense= new char [numz];
	double *obj= new double [numz];
	double *rhs= new double [numz];
	double *lb= new double [numz];
	double *ub= new double [numz];
	int * indices= new int[numz];

	CPXchgobjsen (env_slp, exten, CPX_MAX);	
	char     *ctype = NULL;
	ctype= new char [numz];
	for(int i=0;i<numz;i++)
		ctype[i]='B';		
	for (int i=0;i<num_states;i++)
	for (int j=0;j<num_states;j++)
	{
		if (i!=num_states-1&& j!=num_states-1)
			ub[i*num_states+j]=1;
		else	
			ub[i*num_states+j]=0;
	}

		
	status = CPXnewcols (env_slp, exten, numz, NULL, NULL, ub, ctype, NULL);
	if ( status ) goto TERMINATE;	

	/////adding the constraint set (1) for the first patient
	for(int k=0;k<num_patients;k++){
	for(int i=0;i<numz;i++)	
		obj[i]=0;
	obj[*obj_ind]=1;
	status = CPXnewcols (env_slp, exten, numz, obj, NULL, NULL, NULL, NULL);
	if ( status ) goto TERMINATE;	
	for(int i1=0;i1<num_states;i1++)
	for(int i2=0;i2<num_states;i2++)
	{
		for(int j1=0;j1<num_states;j1++)
		for(int j2=0;j2<num_states;j2++)
		{
			matval[(i1*num_states+i2)*numz+j1*num_states+j2]=(-lambda)*tran_prob[0][i1][j1]*tran_prob[1][i2][j2];
			matind[(i1*num_states+i2)*numz+j1*num_states+j2]=(k+1)*numz+j1*num_states+j2;

		}
		matval[(i1*num_states+i2)*numz+i1*num_states+i2]++;

		matbeg[i1*num_states+i2]=numz*(i1*num_states+i2);
		sense [i1*num_states+i2]='G';				
		rhs [i1*num_states+i2]=rhs_contin[k][i1*num_states+i2];				
	}

	status = CPXaddrows (env_slp, exten, 0, numz,numz*numz, rhs,sense, matbeg, matind,matval, NULL, NULL);
	if ( status ) goto TERMINATE;

	///adding the constraint set (3) for the first patient
	for(int i=0;i<numz;i++)	
	{
		rhs[i]=d[k][i];
		sense [i]='G';
		obj[i]=1;
		indices[i]=(k+1)*numz+i;
		matbeg[i]=i;
	}

	status = CPXaddrows (env_slp, exten, 0, numz,numz, rhs,sense, matbeg, indices,obj, NULL, NULL);

	
	///////adding the constraint set (2) for the first patient
	for(int i=0;i<numz;i++)	
	{
		for(int j=0;j<numz;j++)				
			matind[i*numz+j]=(k+1)*numz+j;

		matbeg[i]=numz*i;
		sense [i]='L';
		indices[i]=i;
		rhs[i]=rhs_contin[k][i];
	}			
			
	status = CPXaddrows (env_slp, exten, 0, numz,numz*numz, rhs,sense, matbeg, matind,matval, NULL, NULL);
	if ( status ) goto TERMINATE;

			   //adding the constraint set (4) for the first patient
	for(int i=0;i<numz;i++)	{
		matbeg[i]=i;
		obj[i]=1;		
		sense [i]='L';		
		indices[i]=(k+1)*numz+i;
	}
	
	status = CPXaddrows (env_slp, exten, 0, numz,numz, v[k],sense, matbeg, indices, obj, NULL, NULL);
	if ( status ) goto TERMINATE;

	}

		 
	for(int k=0;k<num_patients;k++){
		///////constraint set (1)

		///////constraint set (3)
		for (int i=0;i<numz;i++)				
			rhs[i]=-rhs_stop[k][i]+d[k][i];
		
		for(int i=0;i<numz;i++)
		status = CPXchgcoef (env_slp, exten, (4*k+1)*numz+i,i, rhs[i]);
				
		 ///////constraint set (2)
		for (int i=0;i<numz;i++)		
			rhs[i]=-rhs_stop[k][i]+d[k][i];
		
		for(int i=0;i<numz;i++)
		status = CPXchgcoef (env_slp, exten, (4*k+2)*numz+i, i, rhs[i]);

				///////constraint set (4)
		for (int i=0;i<numz;i++)				
			rhs[i]=v[k][i]-rhs_stop[k][i];
		
		for(int i=0;i<numz;i++)
		status = CPXchgcoef (env_slp, exten, (4*k+3)*numz+i, i, rhs[i]);
	}
	
	
	
	TERMINATE:
	delete [] matval;
	delete [] matind;
	delete [] matbeg;
	delete [] matcnt;
	delete [] sense;
	delete [] obj;
	delete [] rhs;
	delete [] lb;
	delete [] ub;
	delete []indices;

	return 0;
	}


	int make_extensive_3(CPXCENVptr env_slp, CPXLPptr exten, double ***tran_prob, double **rhs_contin,double **rhs_stop, double **v,double **d, int *obj_ind){ //constructing the extensive form problem for three player instances
	double *matval= new double [numz*numz];
	int *matind= new int [numz*numz];
	int *matbeg= new int [numz];
	int *matcnt= new int [numz];
	char *sense= new char [numz];
	double *obj= new double [numz];
	double *rhs= new double [numz];
	double *lb= new double [numz];
	double *ub= new double [numz];
	int * indices= new int[numz];

	CPXchgobjsen (env_slp, exten, CPX_MAX);	
	char     *ctype = NULL;
	ctype= new char [numz];
	for(int i=0;i<numz;i++)
		ctype[i]='B';		
	for (int i1=0;i1<num_states;i1++)
	for (int i2=0;i2<num_states;i2++)
	for (int i3=0;i3<num_states;i3++)
	{
		if (i1!=num_states-1&&	i2!=num_states-1&&i3!=num_states-1)
			ub[i1*num_states*num_states+i2*num_states+i3]=1;
		else	
			ub[i1*num_states*num_states+i2*num_states+i3]=0;
	}

		
	status = CPXnewcols (env_slp, exten, numz, NULL, NULL, ub, ctype, NULL);
	if ( status ) goto TERMINATE;	

	/////adding the constraint set (1) for the first patient
	for(int k=0;k<num_patients;k++){
	for(int i=0;i<numz;i++)	
		obj[i]=0;
	obj[*obj_ind]=1;
	status = CPXnewcols (env_slp, exten, numz, obj, NULL, NULL, NULL, NULL);
	if ( status ) goto TERMINATE;	
	for(int i1=0;i1<num_states;i1++)
	for(int i2=0;i2<num_states;i2++)
	for(int i3=0;i3<num_states;i3++)
	{
		for(int j1=0;j1<num_states;j1++)
		for(int j2=0;j2<num_states;j2++)
		for(int j3=0;j3<num_states;j3++)
		{
			matval[(i1*num_states*num_states+i2*num_states+i3)*numz+j1*num_states*num_states+j2*num_states+j3]=(-lambda)*tran_prob[0][i1][j1]*tran_prob[1][i2][j2]*tran_prob[2][i3][j3];
			matind[(i1*num_states*num_states+i2*num_states+i3)*numz+j1*num_states*num_states+j2*num_states+j3]=(k+1)*numz+j1*num_states*num_states+j2*num_states+j3;

		}
		matval[(i1*num_states*num_states+i2*num_states+i3)*numz+i1*num_states*num_states+i2*num_states+i3]++;

		matbeg[i1*num_states*num_states+i2*num_states+i3]=numz*(i1*num_states*num_states+i2*num_states+i3);
		sense [i1*num_states*num_states+i2*num_states+i3]='G';				
		rhs [i1*num_states*num_states+i2*num_states+i3]=rhs_contin[k][i1*num_states*num_states+i2*num_states+i3];				
	}

	status = CPXaddrows (env_slp, exten, 0, numz,numz*numz, rhs,sense, matbeg, matind,matval, NULL, NULL);
	if ( status ) goto TERMINATE;
	
	///adding the constraint set (3) for the first patient
	for(int i=0;i<numz;i++)	
	{
		rhs[i]=d[k][i];
		sense [i]='G';
		obj[i]=1;
		indices[i]=(k+1)*numz+i;
		matbeg[i]=i;
	}

	status = CPXaddrows (env_slp, exten, 0, numz,numz, rhs,sense, matbeg, indices,obj, NULL, NULL);

	
	///////adding the constraint set (2) for the first patient
	for(int i=0;i<numz;i++)	
	{
		for(int j=0;j<numz;j++)				
			matind[i*numz+j]=(k+1)*numz+j;

		matbeg[i]=numz*i;
		sense [i]='L';
		indices[i]=i;
		rhs[i]=rhs_contin[k][i];
	}			
			
	status = CPXaddrows (env_slp, exten, 0, numz,numz*numz, rhs,sense, matbeg, matind,matval, NULL, NULL);
	if ( status ) goto TERMINATE;

			   //adding the constraint set (4) for the first patient
	for(int i=0;i<numz;i++)	{
		matbeg[i]=i;
		obj[i]=1;		
		sense [i]='L';		
		indices[i]=(k+1)*numz+i;
	}
	
	status = CPXaddrows (env_slp, exten, 0, numz,numz, v[k],sense, matbeg, indices, obj, NULL, NULL);
	if ( status ) goto TERMINATE;
	
	}
 
	for(int k=0;k<num_patients;k++){
		///////constraint set (1)

		///////constraint set (3)
		for (int i=0;i<numz;i++)				
			rhs[i]=-rhs_stop[k][i]+d[k][i];
		
		for(int i=0;i<numz;i++)
		status = CPXchgcoef (env_slp, exten, (4*k+1)*numz+i,i, rhs[i]);
				
		 ///////constraint set (2)
		for (int i=0;i<numz;i++)		
			rhs[i]=-rhs_stop[k][i]+d[k][i];
		
		for(int i=0;i<numz;i++)
		status = CPXchgcoef (env_slp, exten, (4*k+2)*numz+i, i, rhs[i]);

				///////constraint set (4)
		for (int i=0;i<numz;i++)				
			rhs[i]=v[k][i]-rhs_stop[k][i];
		
		for(int i=0;i<numz;i++)
		status = CPXchgcoef (env_slp, exten, (4*k+3)*numz+i, i, rhs[i]);
	}
	
	//status=CPXwriteprob(env_slp,exten,"exten.lp",NULL);
	
	TERMINATE:
	delete [] matval;
	delete [] matind;
	delete [] matbeg;
	delete [] matcnt;
	delete [] sense;
	delete [] obj;
	delete [] rhs;
	delete [] lb;
	delete [] ub;
	delete []indices;

	return 0;
	}


int readdata3( double ***tran_prob, double **reward,double **stop_reward, int *obj_ind, int *instance_num) //read data for three player instances
	{
		int patient_1; // Determines which lines of  the text file for ptobability matrix and immediate rewards should be read for patient 1
		int patient_2; // Determines which lines of  the text file for ptobability matrix and immediate rewards should be read for patient 2
		int patient_3; // Determines which lines of  the text file for ptobability matrix and immediate rewards should be read for patient 2
		//reading the initial state
		int row_counter=1;
		stringstream ss_nstates;
		ss_nstates << num_states;
		string str1 = ss_nstates.str();		
		string str2="GFRs-"+str1+".txt";

		ifstream ifile(str2);
		string line;
		while ( getline(ifile, line)){
			if(row_counter==*instance_num){
				vector < int > data;
				int value;
				istringstream iss(line);
				while (iss >> value)
				{
					data.push_back(value);
				}				
					*obj_ind=data[0]*num_states*num_states+data[1]*num_states+data[2];
			}
			row_counter ++;
		}
		ifile.close();

		// reading the types of patients and stopping rewards
		row_counter=1;
		vector < vector < double > > info;
		str2="Inputs-"+str1+".txt";

		ifstream sfile(str2);
		while ( getline(sfile, line))
		{			
			if(row_counter==*instance_num){
				vector < double > data;
				double value;
				istringstream iss(line);
				while (iss >> value)
				{
					data.push_back(value);
				}
				info.push_back(data);
			}
			
			row_counter += 1;
		}

		patient_1=2*info[0][0]+info[0][1];
		patient_2=2*info[0][4]+info[0][5];
		patient_3=2*info[0][8]+info[0][9];

		for(int i=0;i<num_states-2;i++)
			stop_reward[0][i]=info[0][2];
		stop_reward[0][num_states-2]=info[0][3];
		stop_reward[0][num_states-1]=0;

		for(int i=0;i<num_states-2;i++)
			stop_reward[1][i]=info[0][6];
		stop_reward[1][num_states-2]=info[0][7];
		stop_reward[1][num_states-1]=0;

		for(int i=0;i<num_states-2;i++)
			stop_reward[2][i]=info[0][10];
		stop_reward[2][num_states-2]=info[0][11];
		stop_reward[2][num_states-1]=0;
		sfile.close();

		row_counter=0;		
		info.clear();
		str2="Rewards-"+str1+".txt";

		ifstream file(str2);		
		while ( getline(file, line) )
		{
			if(row_counter==patient_1){
				vector < double > data;
				double value;
				istringstream iss(line);
				while (iss >> value)
				{
					data.push_back(value);
				}
					info.push_back(data);
			}
			if(row_counter==patient_2){
				vector < double > data;
				double value;
				istringstream iss(line);
				while (iss >> value)
				{
					data.push_back(value);
				}
					info.push_back(data);
			}
			if(row_counter==patient_3){
				vector < double > data;
				double value;
				istringstream iss(line);
				while (iss >> value)
				{
					data.push_back(value);
				}
					info.push_back(data);
			}
			row_counter += 1;
		}
		for ( vector < vector < double > > :: size_type i = 0, size = info.size(); i < size; ++i)
		{
			for ( vector < double > :: size_type j = 0, length = info[i].size(); j < length; ++j)
			{	
				reward[i][j]=info[i][j];					
			}			
		}
		file.close();

		row_counter=0;
		info.clear();
		str2="Probabilities-"+str1+".txt";

		ifstream pfile(str2);
		while ( getline(pfile, line))
		{		
				vector < double > data;
				double value;
				istringstream iss(line);
				while (iss >> value)
				{
					data.push_back(value);
				}
				info.push_back(data);		
			row_counter += 1;
		}
	
		
		for(int j=0;j<num_states;j++)
		{	
			for(int k=0;k<num_states;k++)
			{
				tran_prob[0][j][k]=info[(patient_1)* num_states +j][k];
				tran_prob[1][j][k]=info[(patient_2)* num_states +j][k];					
				tran_prob[2][j][k]=info[(patient_3)* num_states +j][k];	
			}
		}

		pfile.close();
		info.clear();
		
		return 0;
	}
	

int readdata2( double ***tran_prob, double **reward,double **stop_reward, int *obj_ind, int *instance_num)// read data for two player instances
	{

		int patient_1;
		int patient_2;
		//reading the initial state
		int row_counter=1;		
		stringstream ss_nstates;
		ss_nstates << num_states;
		string str1 = ss_nstates.str();		
		string str2="GFRs-"+str1+".txt";
				
		ifstream ifile(str2);
		string line;
		while ( getline(ifile, line)){
			if(row_counter==*instance_num){
				vector < int > data;
				int value;
				istringstream iss(line);
				while (iss >> value)
				{
					data.push_back(value);
				}				
					*obj_ind=data[0]*num_states+data[1];
			}
			row_counter ++;
		}
		ifile.close();
		// reading the types of patients and stopping rewards
		row_counter=1;
		vector < vector < double > > info;
		str2="Inputs-"+str1+".txt";
		ifstream sfile(str2);
		while ( getline(sfile, line))
		{			
			if(row_counter==*instance_num){
				vector < double > data;
				double value;
				istringstream iss(line);
				while (iss >> value)
				{
					data.push_back(value);
				}
				info.push_back(data);
			}
			
			row_counter += 1;
		}

		patient_1=2*info[0][0]+info[0][1];
		patient_2=2*info[0][4]+info[0][5];

		for(int i=0;i<num_states-2;i++)
			stop_reward[0][i]=info[0][2];
		stop_reward[0][num_states-2]=info[0][3];
		stop_reward[0][num_states-1]=0;

		for(int i=0;i<num_states-2;i++)
			stop_reward[1][i]=info[0][6];
		stop_reward[1][num_states-2]=info[0][7];
		stop_reward[1][num_states-1]=0;
		sfile.close();

		row_counter=0;		
		info.clear();
		str2="Rewards-"+str1+".txt";

		ifstream file(str2);		
		while ( getline(file, line) )
		{
			if(row_counter==patient_1){
				vector < double > data;
				double value;
				istringstream iss(line);
				while (iss >> value)
				{
					data.push_back(value);
				}
					info.push_back(data);
			}
			if(row_counter==patient_2){
				vector < double > data;
				double value;
				istringstream iss(line);
				while (iss >> value)
				{
					data.push_back(value);
				}
					info.push_back(data);
			}
			row_counter += 1;
		}
		for ( vector < vector < double > > :: size_type i = 0, size = info.size(); i < size; ++i)
		{
			for ( vector < double > :: size_type j = 0, length = info[i].size(); j < length; ++j)
			{	
				reward[i][j]=info[i][j];					
			}			
		}
		file.close();

		row_counter=0;
		info.clear();
		str2="Probabilities-"+str1+".txt";

		ifstream pfile(str2);
		while ( getline(pfile, line))
		{		
				vector < double > data;
				double value;
				istringstream iss(line);
				while (iss >> value)
				{
					data.push_back(value);
				}
				info.push_back(data);		
			row_counter += 1;
		}
	
		
		for(int j=0;j<num_states;j++)
		{	
			for(int k=0;k<num_states;k++)
			{
				tran_prob[0][j][k]=info[(patient_1)* num_states +j][k];
				tran_prob[1][j][k]=info[(patient_2)* num_states +j][k];					
			}
		}

		pfile.close();
		info.clear();
		
		return 0;
	}




int make_master_negative(CPXCENVptr env,CPXCENVptr env_slp,CPXLPptr lp, CPXLPptr* slp, CPXLPptr* slp_v, double ***tran_prob, double **rhs_contin,double **rhs_stop, double **v,double **d, int *obj_ind)		//construct the master problem for negative reward instances
	{		
#pragma region  start creating sub_problems for the first patient and calculating d_1(s)
	double *matval= new double [numz*numz];
	int *matind= new int [numz*numz];
	int *matbeg= new int [numz];
	int *matcnt= new int [numz];
	char *sense= new char [numz];
	double *obj= new double [numz];
	double *rhs= new double [numz];
	double *lb= new double [numz];
	double *ub= new double [numz];
	int *indices = new int[numz];

	for(int i1=0;i1<num_states;i1++)
	for(int i2=0;i2<num_states;i2++)
	{
		for(int j1=0;j1<num_states;j1++)
		for(int j2=0;j2<num_states;j2++)
		{
			matval[(i1*num_states+i2)*numz+j1*num_states+j2]=(-lambda)*tran_prob[0][j1][i1]*tran_prob[1][j2][i2];
			matind[(i1*num_states+i2)*numz+j1*num_states+j2]=j1*num_states+j2;

		}
		matval[(i1*num_states+i2)*numz+i1*num_states+i2]++;

		matbeg[i1*num_states+i2]=numz*(i1*num_states+i2);
		matcnt[i1*num_states+i2]=numz;
		sense [i1*num_states+i2]='E';
		obj[i1*num_states+i2]=0;
		lb[i1*num_states+i2]=-CPX_INFBOUND;;
		ub[i1*num_states+i2]=CPX_INFBOUND;
	}
	
	obj[*obj_ind]=1;

	status = CPXcopylp (env_slp, slp[0], numz, numz, CPX_MAX, obj, rhs_contin[0],
                     sense, matbeg, matcnt, matind, matval, lb,
                     ub, NULL);
	if ( status ) {
		fprintf (stderr, "CPXcopylp failed to populate subproblem 1.\n");
		goto TERMINATE;
	}
	status=CPXlpopt(env_slp,slp[0]);
	if ( status ) {
		fprintf (stderr, "CPXlpopt failed to optimizae subproblem 1.\n");
		goto TERMINATE;
	}

	status = CPXgetx (env_slp, slp[0], d[0], 0, numz-1);
	if ( status ) {
		fprintf (stderr, "CPXgetx failed to obtain d values of subproblem 1.\n");
		goto TERMINATE;
	}
#pragma endregion
	


#pragma region start creating sub_problems for the second patient and calculating d_2(s)
	status = CPXcopylp (env_slp, slp[1], numz, numz, CPX_MAX, obj, rhs_contin[1],
                     sense, matbeg, matcnt, matind, matval, lb,
                     ub, NULL);
	if ( status ) {
		fprintf (stderr, "CPXcopylp failed to populate subproblem 2.\n");
		goto TERMINATE;
	}
	
	status=CPXlpopt(env_slp,slp[1]);
	if ( status ) {
		fprintf (stderr, "CPXlpopt failed to optimizae subproblem 1.\n");
		goto TERMINATE;
	}

	status = CPXgetx (env_slp, slp[1], d[1], 0, numz-1);
	if ( status ) {
		fprintf (stderr, "CPXgetx failed to obtain d values of subproblem 1.\n");
		goto TERMINATE;
	}

#pragma endregion	


#pragma region calculating v in the following

	slp_v[0] = CPXcloneprob (env_slp, slp[0], &status);
	if ( status ) {
	fprintf (stderr, "CPXcopylp failed to creates a copy of slp problem.\n");
	goto TERMINATE;
	}
	


	for(int i1=0;i1<num_states;i1++)
	for(int i2=0;i2<num_states;i2++)
	{
		rhs[i1*num_states+i2]=rhs_stop[0][i1*num_states+i2];

		sense [i1*num_states+i2]='G';
		obj[i1*num_states+i2]=1;
		indices[i1*num_states+i2]=i1*num_states+i2;
		matbeg[i1*num_states+i2]=i1*num_states+i2;
	}

	for (int i=0;i<num_states;i++)
	for (int j=0;j<num_states;j++)
	{
		if (i==num_states-1|| j==num_states-1)
			rhs[i*num_states+j]=-CPX_INFBOUND;
	}
	for (int i=0;i<numz;i++)	
	if(d[0][i]>rhs_stop[0][i]+EPSILON || d[1][i]>rhs_stop[1][i]+EPSILON)
		rhs[i]=-CPX_INFBOUND; //this implements the column elimination idea for calulating V

	status = CPXchgobj (env_slp, slp_v[0], numz, indices, obj);
	if ( status ) {
	fprintf (stderr, "CPXcopylp failed to change obj fun.\n");
	goto TERMINATE;
	}

	status = CPXaddrows (env_slp, slp_v[0], 0, numz,numz, rhs,sense, matbeg, indices,
                        obj, NULL, NULL);
	if ( status ) {
	fprintf (stderr, "CPXcopylp failed to add columns.\n");
	goto TERMINATE;
	}
	
	status = CPXchgsense (env_slp, slp_v[0], numz, indices, sense);
	CPXchgobjsen (env_slp, slp_v[0], CPX_MIN);

	status = CPXlpopt (env_slp, slp_v[0]);
	if ( status ) {
	fprintf (stderr, "CPXcopylp failed to optimize copy_slp_v.\n");
	goto TERMINATE;
	}	
 
	status = CPXgetx (env_slp, slp_v[0], v[0], 0, numz-1);
	if ( status ) {
	fprintf (stderr, "CPXcopylp failed to access solution values copy_slp_v.\n");
	goto TERMINATE;
	}

	//////////////calculatin of V2/////////
	slp_v[1] = CPXcloneprob (env_slp, slp_v[0], &status);
	if ( status ) {
	fprintf (stderr, "CPXcopylp failed to creates a copy of slp problem.\n");
	goto TERMINATE;
	}
	status = CPXchgrhs (env_slp, slp_v[1], numz, indices, rhs_contin[1]);
	if ( status ) {
	fprintf (stderr, "Failed to change righthand-side of the subproblem_V .\n");
	goto TERMINATE;
		}

	for (int j1=0; j1<num_states;j1++)
	for (int j2=0; j2<num_states;j2++)
	{
		rhs[j1*num_states+j2]=rhs_stop[1][j1*num_states+j2];
		indices[j1*num_states+j2]=numz+j1*num_states+j2;
	}
	
	for (int i=0;i<num_states;i++)
	for (int j=0;j<num_states;j++)
	{
		if (i==num_states-1|| j==num_states-1)
			rhs[i*num_states+j]=-CPX_INFBOUND;
	}

	for (int i=0;i<numz;i++)	
	if(d[0][i]>rhs_stop[0][i]+EPSILON || d[1][i]>rhs_stop[1][i]+EPSILON)
		rhs[i]=-CPX_INFBOUND; //this implements the column elimination idea for calulating V

	status = CPXchgrhs (env_slp, slp_v[1], numz, indices, rhs);
	if ( status ) {
	fprintf (stderr, "Failed to change righthand-side of the subproblem_V .\n");
	goto TERMINATE;
		}


	status = CPXlpopt (env_slp, slp_v[1]);
	if ( status ) {
	fprintf (stderr, "CPXcopylp failed to optimize copy_slp_v.\n");
	goto TERMINATE;
	}

	status = CPXgetx (env_slp, slp_v[1], v[1], 0, numz-1);
	if ( status ) {
	fprintf (stderr, "CPXcopylp failed to access solution values copy_slp_v.\n");
	goto TERMINATE;
	}
	
#pragma endregion

#pragma region start creating master problems 

	char     *ctype = NULL;
	ctype= new char [numz];
	for(int i=0;i<numz;i++)
		ctype[i]='B';		
	for (int i=0;i<num_states;i++)
	for (int j=0;j<num_states;j++)
	{
		if (i!=num_states-1&& j!=num_states-1)
			ub[i*num_states+j]=1;
		else	
			ub[i*num_states+j]=0;
	}

		
	for (int i=0;i<numz;i++)	
	if(d[0][i]>rhs_stop[0][i]+EPSILON || d[1][i]>rhs_stop[1][i]+EPSILON)		
		ub[i]=0;
		

	status = CPXnewcols (env, lp, numz, NULL, NULL, ub, ctype, NULL);

	if ( status ) {
	fprintf (stderr, "Failed to create LP.\n");
	goto TERMINATE;
	}	
	double *z=new double [num_patients];
	double *upperbound=new double [num_patients];
	double *lowerbound=new double [num_patients];
	
	
	for (int i=0;i<num_patients;i++)
	{
		z[i]=1;
		upperbound[i]=v[i][*obj_ind];
		lowerbound[i]=-CPX_INFBOUND;
	}

	
	status = CPXaddcols (env, lp, num_patients, 0, z, NULL, NULL, NULL, lowerbound, upperbound, NULL);
	if ( status ) {
	fprintf (stderr, "Failed to create LP.\n");
	goto TERMINATE;
	}	

	CPXchgobjsen (env, lp, CPX_MAX);

	delete[]z;
	delete[]ctype;
	delete[]upperbound;					
	delete[]lowerbound;					
	
	
#pragma endregion
	
TERMINATE:

#pragma region	deallocating dynamic memories
	delete [] matval;
	delete [] matind;
	delete [] matbeg;
	delete [] matcnt;
	delete [] sense;
	delete [] obj;
	delete [] rhs;
	delete [] lb;
	delete [] ub;
	delete []indices;
#pragma endregion				
		return 0;
	}



int make_slpvt_negative(CPXCENVptr env_slp, CPXLPptr slp_t, CPXLPptr slp_vt, double ***tran_prob, double **rhs_contin,double **rhs_stop, double *vt,double **d, int *obj_ind){ //construct the slp_vt and slp_t for negative reward instances
#pragma region	making slp_vt
	double *matval= new double [numz*numz];
	int *matind= new int [numz*numz];
	int *matbeg= new int [numz];
	int *matcnt= new int [numz];
	char *sense= new char [numz];
	double *obj= new double [numz];
	double *rhs= new double [numz];
	double *lb= new double [numz];
	double *ub= new double [numz];
	int * indices= new int[numz];

	for(int i1=0;i1<num_states;i1++)
	for(int i2=0;i2<num_states;i2++)
	{
		for(int j1=0;j1<num_states;j1++)
		for(int j2=0;j2<num_states;j2++)
		{
			matval[(i1*num_states+i2)*numz+j1*num_states+j2]=(-lambda)*tran_prob[0][j1][i1]*tran_prob[1][j2][i2];
			matind[(i1*num_states+i2)*numz+j1*num_states+j2]=j1*num_states+j2;

		}
		matval[(i1*num_states+i2)*numz+i1*num_states+i2]++;

		matbeg[i1*num_states+i2]=numz*(i1*num_states+i2);
		matcnt[i1*num_states+i2]=numz;
		rhs[i1*num_states+i2]=rhs_contin[0][i1*num_states+i2]+rhs_contin[1][i1*num_states+i2];
		sense [i1*num_states+i2]='G';
		obj[i1*num_states+i2]=1;
		lb[i1*num_states+i2]=-CPX_INFBOUND;
		ub[i1*num_states+i2]=CPX_INFBOUND;
	}
		

	status = CPXcopylp (env_slp, slp_vt, numz, numz, CPX_MIN, obj, rhs,sense, matbeg, matcnt, matind, matval, lb,ub, NULL);
	if ( status ) 	goto TERMINATE;
	
	
	for(int i1=0;i1<num_states;i1++)
	for(int i2=0;i2<num_states;i2++)
	{
		rhs[i1*num_states+i2]=rhs_stop[0][i1*num_states+i2]+rhs_stop[1][i1*num_states+i2];

		sense [i1*num_states+i2]='G';
		obj[i1*num_states+i2]=1;
		indices[i1*num_states+i2]=i1*num_states+i2;
		matbeg[i1*num_states+i2]=i1*num_states+i2;
	}

	for (int i=0;i<numz;i++)	
	if(d[0][i]>rhs_stop[0][i]+EPSILON || d[1][i]>rhs_stop[1][i]+EPSILON)
		rhs[i]=-CPX_INFBOUND; 
	status = CPXaddrows (env_slp, slp_vt, 0, numz,numz, rhs,sense, matbeg, indices,obj, NULL, NULL);
	status = CPXlpopt (env_slp, slp_vt);	
	status = CPXgetx (env_slp, slp_vt, vt, 0, numz-1);

#pragma endregion

	for(int i1=0;i1<num_states;i1++)
	for(int i2=0;i2<num_states;i2++)
	{		
		matbeg[i1*num_states+i2]=numz*(i1*num_states+i2);
		matcnt[i1*num_states+i2]=numz;
		rhs[i1*num_states+i2]=rhs_contin[0][i1*num_states+i2];
		sense [i1*num_states+i2]='E';
		obj[i1*num_states+i2]=0;
		lb[i1*num_states+i2]=-CPX_INFBOUND;
		ub[i1*num_states+i2]=CPX_INFBOUND;
	}
	obj[*obj_ind]=1;	

	status = CPXcopylp (env_slp, slp_t, numz, numz, CPX_MAX, obj, rhs,sense, matbeg, matcnt, matind, matval, lb,ub, NULL);
	if ( status ) 	goto TERMINATE;

	for(int i=0;i<numz;i++)	
		obj[i]=0;
	obj[*obj_ind]=1;
	status = CPXnewcols (env_slp, slp_t, numz, obj, NULL, NULL, NULL, NULL);
	if ( status ) goto TERMINATE;	

	for(int i1=0;i1<num_states;i1++)
	for(int i2=0;i2<num_states;i2++)
	{
		for(int j1=0;j1<num_states;j1++)
		for(int j2=0;j2<num_states;j2++)
		{
			matval[(i1*num_states+i2)*numz+j1*num_states+j2]=(-lambda)*tran_prob[0][i1][j1]*tran_prob[1][i2][j2];
			matind[(i1*num_states+i2)*numz+j1*num_states+j2]=numz+j1*num_states+j2;

		}
		matval[(i1*num_states+i2)*numz+i1*num_states+i2]++;

		matbeg[i1*num_states+i2]=numz*(i1*num_states+i2);
		sense [i1*num_states+i2]='E';				
		rhs [i1*num_states+i2]=rhs_contin[1][i1*num_states+i2];				
	}

	status = CPXaddrows (env_slp, slp_t, 0, numz,numz*numz, rhs,sense, matbeg, matind,matval, NULL, NULL);
	if ( status ) goto TERMINATE;

	double *val= new double [num_patients*numz];
	int *ind= new int [num_patients*numz];
	double cutrhs=100000;
	int beg=0;
	for(int j=0;j<num_patients*numz;j++)
		ind[j]=j;
	for(int j=0;j<numz;j++){
		for(int k=0;k<2*numz;k++)
			val[k]=0;
		val[j]=1, val[numz+j]=1;
		status=CPXaddrows(env_slp, slp_t, 0, 1, num_patients*numz, &cutrhs, "L", &beg, ind, val, NULL, NULL);
	}

	delete[] val;
	delete[] ind;
	

TERMINATE:
	delete [] matval;
	delete [] matind;
	delete [] matbeg;
	delete [] matcnt;
	delete [] sense;
	delete [] obj;
	delete [] rhs;
	delete [] lb;
	delete [] ub;
	delete []indices;

	return 0;
	}



int readdata_negative( double ***tran_prob, double **reward,double **stop_reward, int *obj_ind, int *instance_num) //read the data for negative reward instances
	{

		int patient_1;
		int patient_2;
		//reading the initial state
		int row_counter=1;		
		stringstream ss_nstates;
		ss_nstates << num_states;
		string str1 = ss_nstates.str();		
		string str2="GFRs-Negative-"+str1+".txt";
				
		ifstream ifile(str2);		
		string line;
		while ( getline(ifile, line)){
			if(row_counter==*instance_num){
				vector < int > data;
				int value;
				istringstream iss(line);
				while (iss >> value)
				{
					data.push_back(value);
				}				
					*obj_ind=data[0]*num_states+data[1];
			}
			row_counter ++;
		}

		// reading the types of patients and stopping rewards
		row_counter=1;
		vector < vector < double > > info;
		str2="Inputs-Negative-"+str1+".txt";
		ifstream sfile(str2);		
		while ( getline(sfile, line))
		{			
			if(row_counter==*instance_num){
				vector < double > data;
				double value;
				istringstream iss(line);
				while (iss >> value)
				{
					data.push_back(value);
				}
				info.push_back(data);
			}
			
			row_counter += 1;
		}

		patient_1=2*info[0][0]+info[0][1];
		patient_2=2*info[0][num_states+1]+info[0][num_states+2];

		for(int i=0;i<num_states-1;i++)
			stop_reward[0][i]=info[0][2+i];

		stop_reward[0][num_states-1]=0;

		for(int i=0;i<num_states-1;i++)
			stop_reward[1][i]=info[0][num_states+3+i];		
		stop_reward[1][num_states-1]=0;
		sfile.close();

		row_counter=0;		
		info.clear();
		str2="Rewards-Negative-"+str1+".txt";
		ifstream file(str2);
		while ( getline(file, line) )
		{
			if(row_counter==patient_1){
				vector < double > data;
				double value;
				istringstream iss(line);
				while (iss >> value)
				{
					data.push_back(value);
				}
					info.push_back(data);
			}
			if(row_counter==patient_2){
				vector < double > data;
				double value;
				istringstream iss(line);
				while (iss >> value)
				{
					data.push_back(value);
				}
					info.push_back(data);
			}
			row_counter += 1;
		}
		for ( vector < vector < double > > :: size_type i = 0, size = info.size(); i < size; ++i)
		{
			for ( vector < double > :: size_type j = 0, length = info[i].size(); j < length; ++j)
			{	
				reward[i][j]=info[i][j];					
			}			
		}
		file.close();

		row_counter=0;
		info.clear();
		str2="Probabilities-Negative-"+str1+".txt";

		ifstream pfile(str2);
		while ( getline(pfile, line))
		{		
				vector < double > data;
				double value;
				istringstream iss(line);
				while (iss >> value)
				{
					data.push_back(value);
				}
				info.push_back(data);		
			row_counter += 1;
		}
	
		
		for(int j=0;j<num_states;j++)
		{	
			for(int k=0;k<num_states;k++)
			{
				tran_prob[0][j][k]=info[(patient_1)* num_states +j][k];
				tran_prob[1][j][k]=info[(patient_2)* num_states +j][k];					
			}
		}

		pfile.close();
		info.clear();
		
		return 0;
	}	

int parameters_calulator_negative(CPXCENVptr env_slp, CPXLPptr* slp, CPXLPptr* slp_v, double ***tran_prob, double **rhs_contin,double **rhs_stop, double **v,double **d, int *obj_ind)	//calculating parameters of extensive form problem for negative-reward instances	
	{		
#pragma region  start creating sub_problems for the first patient and calculating d_1(s)
	double *matval= new double [numz*numz];
	int *matind= new int [numz*numz];
	int *matbeg= new int [numz];
	int *matcnt= new int [numz];
	char *sense= new char [numz];
	double *obj= new double [numz];
	double *rhs= new double [numz];
	double *lb= new double [numz];
	double *ub= new double [numz];
	int *indices = new int[numz];

	for(int i1=0;i1<num_states;i1++)
	for(int i2=0;i2<num_states;i2++)
	{
		for(int j1=0;j1<num_states;j1++)
		for(int j2=0;j2<num_states;j2++)
		{
			matval[(i1*num_states+i2)*numz+j1*num_states+j2]=(-lambda)*tran_prob[0][j1][i1]*tran_prob[1][j2][i2];
			matind[(i1*num_states+i2)*numz+j1*num_states+j2]=j1*num_states+j2;

		}
		matval[(i1*num_states+i2)*numz+i1*num_states+i2]++;

		matbeg[i1*num_states+i2]=numz*(i1*num_states+i2);
		matcnt[i1*num_states+i2]=numz;
		sense [i1*num_states+i2]='E';
		obj[i1*num_states+i2]=0;
		lb[i1*num_states+i2]=-CPX_INFBOUND;
		ub[i1*num_states+i2]=CPX_INFBOUND;
	}
	
	obj[*obj_ind]=1;

	status = CPXcopylp (env_slp, slp[0], numz, numz, CPX_MAX, obj, rhs_contin[0],
                     sense, matbeg, matcnt, matind, matval, lb,
                     ub, NULL);
	if ( status ) {
		fprintf (stderr, "CPXcopylp failed to populate subproblem 1.\n");
		goto TERMINATE;
	}
	status=CPXlpopt(env_slp,slp[0]);
	if ( status ) {
		fprintf (stderr, "CPXlpopt failed to optimizae subproblem 1.\n");
		goto TERMINATE;
	}

	status = CPXgetx (env_slp, slp[0], d[0], 0, numz-1);
	if ( status ) {
		fprintf (stderr, "CPXgetx failed to obtain d values of subproblem 1.\n");
		goto TERMINATE;
	}
#pragma endregion
	


#pragma region start creating sub_problems for the second patient and calculating d_2(s)
	status = CPXcopylp (env_slp, slp[1], numz, numz, CPX_MAX, obj, rhs_contin[1],
                     sense, matbeg, matcnt, matind, matval, lb,
                     ub, NULL);
	if ( status ) {
		fprintf (stderr, "CPXcopylp failed to populate subproblem 2.\n");
		goto TERMINATE;
	}
	
	status=CPXlpopt(env_slp,slp[1]);
	if ( status ) {
		fprintf (stderr, "CPXlpopt failed to optimizae subproblem 1.\n");
		goto TERMINATE;
	}

	status = CPXgetx (env_slp, slp[1], d[1], 0, numz-1);
	if ( status ) {
		fprintf (stderr, "CPXgetx failed to obtain d values of subproblem 1.\n");
		goto TERMINATE;
	}

#pragma endregion	


#pragma region calculating v in the following

	slp_v[0] = CPXcloneprob (env_slp, slp[0], &status);
	if ( status ) {
	fprintf (stderr, "CPXcopylp failed to creates a copy of slp problem.\n");
	goto TERMINATE;
	}
	


	for(int i1=0;i1<num_states;i1++)
	for(int i2=0;i2<num_states;i2++)
	{
		rhs[i1*num_states+i2]=rhs_stop[0][i1*num_states+i2];

		sense [i1*num_states+i2]='G';
		obj[i1*num_states+i2]=1;
		indices[i1*num_states+i2]=i1*num_states+i2;
		matbeg[i1*num_states+i2]=i1*num_states+i2;
	}

	for (int i=0;i<num_states;i++)
	for (int j=0;j<num_states;j++)
	{
		if (i==num_states-1|| j==num_states-1)
			rhs[i*num_states+j]=-CPX_INFBOUND;
	}
	
	status = CPXchgobj (env_slp, slp_v[0], numz, indices, obj);
	if ( status ) {
	fprintf (stderr, "CPXcopylp failed to change obj fun.\n");
	goto TERMINATE;
	}

	status = CPXaddrows (env_slp, slp_v[0], 0, numz,numz, rhs,sense, matbeg, indices,
                        obj, NULL, NULL);
	if ( status ) {
	fprintf (stderr, "CPXcopylp failed to add columns.\n");
	goto TERMINATE;
	}
	
	status = CPXchgsense (env_slp, slp_v[0], numz, indices, sense);
	CPXchgobjsen (env_slp, slp_v[0], CPX_MIN);

	status = CPXlpopt (env_slp, slp_v[0]);
	if ( status ) {
	fprintf (stderr, "CPXcopylp failed to optimize copy_slp_v.\n");
	goto TERMINATE;
	}	
 
	status = CPXgetx (env_slp, slp_v[0], v[0], 0, numz-1);
	if ( status ) {
	fprintf (stderr, "CPXcopylp failed to access solution values copy_slp_v.\n");
	goto TERMINATE;
	}

	//////////////calculatin of V2/////////
	slp_v[1] = CPXcloneprob (env_slp, slp_v[0], &status);
	if ( status ) {
	fprintf (stderr, "CPXcopylp failed to creates a copy of slp problem.\n");
	goto TERMINATE;
	}
	status = CPXchgrhs (env_slp, slp_v[1], numz, indices, rhs_contin[1]);
	if ( status ) {
	fprintf (stderr, "Failed to change righthand-side of the subproblem_V .\n");
	goto TERMINATE;
		}

	for (int j1=0; j1<num_states;j1++)
	for (int j2=0; j2<num_states;j2++)
	{
		rhs[j1*num_states+j2]=rhs_stop[1][j1*num_states+j2];
		indices[j1*num_states+j2]=numz+j1*num_states+j2;
	}

	for (int i=0;i<num_states;i++)
	for (int j=0;j<num_states;j++)
	{
		if (i==num_states-1|| j==num_states-1)
			rhs[i*num_states+j]=-CPX_INFBOUND;
	}
	
	status = CPXchgrhs (env_slp, slp_v[1], numz, indices, rhs);
	if ( status ) {
	fprintf (stderr, "Failed to change righthand-side of the subproblem_V .\n");
	goto TERMINATE;
		}


	status = CPXlpopt (env_slp, slp_v[1]);
	if ( status ) {
	fprintf (stderr, "CPXcopylp failed to optimize copy_slp_v.\n");
	goto TERMINATE;
	}

	status = CPXgetx (env_slp, slp_v[1], v[1], 0, numz-1);
	if ( status ) {
	fprintf (stderr, "CPXcopylp failed to access solution values copy_slp_v.\n");
	goto TERMINATE;
	}
	
#pragma endregion


TERMINATE:
#pragma region	deallocating dynamic memories
	delete [] matval;
	delete [] matind;
	delete [] matbeg;
	delete [] matcnt;
	delete [] sense;
	delete [] obj;
	delete [] rhs;
	delete [] lb;
	delete [] ub;
	delete []indices;
#pragma endregion				
		return 0;
	}


int make_extensive_negative(CPXCENVptr env_slp, CPXLPptr exten, double ***tran_prob, double **rhs_contin,double **rhs_stop, double **v,double **d, int *obj_ind){////constructing the extensive form problem for negative reward instances
	double *matval= new double [numz*numz];
	int *matind= new int [numz*numz];
	int *matbeg= new int [numz];
	int *matcnt= new int [numz];
	char *sense= new char [numz];
	double *obj= new double [numz];
	double *rhs= new double [numz];
	double *lb= new double [numz];
	double *ub= new double [numz];
	int * indices= new int[numz];

	CPXchgobjsen (env_slp, exten, CPX_MAX);	
	char     *ctype = NULL;
	ctype= new char [numz];
	for(int i=0;i<numz;i++)
		ctype[i]='B';		
	for (int i=0;i<num_states;i++)
	for (int j=0;j<num_states;j++)
	{
		if (i!=num_states-1&& j!=num_states-1)
			ub[i*num_states+j]=1;
		else	
			ub[i*num_states+j]=0;
	}

		
	
		
	status = CPXnewcols (env_slp, exten, numz, NULL, NULL, ub, ctype, NULL);
	if ( status ) goto TERMINATE;	

	/////adding the constraint set (1) for the first patient
	for(int k=0;k<num_patients;k++){
	for(int i=0;i<numz;i++)	{
		obj[i]=0;
		lb[i]=-CPX_INFBOUND;
	}
	obj[*obj_ind]=1;
	status = CPXnewcols (env_slp, exten, numz, obj, lb, NULL, NULL, NULL);
	if ( status ) goto TERMINATE;	
	for(int i1=0;i1<num_states;i1++)
	for(int i2=0;i2<num_states;i2++)
	{
		for(int j1=0;j1<num_states;j1++)
		for(int j2=0;j2<num_states;j2++)
		{
			matval[(i1*num_states+i2)*numz+j1*num_states+j2]=(-lambda)*tran_prob[0][i1][j1]*tran_prob[1][i2][j2];
			matind[(i1*num_states+i2)*numz+j1*num_states+j2]=(k+1)*numz+j1*num_states+j2;

		}
		matval[(i1*num_states+i2)*numz+i1*num_states+i2]++;

		matbeg[i1*num_states+i2]=numz*(i1*num_states+i2);
		sense [i1*num_states+i2]='G';				
		rhs [i1*num_states+i2]=rhs_contin[k][i1*num_states+i2];				
	}

	status = CPXaddrows (env_slp, exten, 0, numz,numz*numz, rhs,sense, matbeg, matind,matval, NULL, NULL);
	if ( status ) goto TERMINATE;

	///adding the constraint set (3) for the first patient
	for(int i=0;i<numz;i++)	
	{
		rhs[i]=d[k][i];
		sense [i]='G';
		obj[i]=1;
		indices[i]=(k+1)*numz+i;
		matbeg[i]=i;
	}

	status = CPXaddrows (env_slp, exten, 0, numz,numz, rhs,sense, matbeg, indices,obj, NULL, NULL);

	
	///////adding the constraint set (2) for the first patient
	for(int i=0;i<numz;i++)	
	{
		for(int j=0;j<numz;j++)				
			matind[i*numz+j]=(k+1)*numz+j;

		matbeg[i]=numz*i;
		sense [i]='L';
		indices[i]=i;
		rhs[i]=rhs_contin[k][i];
	}			
			
	status = CPXaddrows (env_slp, exten, 0, numz,numz*numz, rhs,sense, matbeg, matind,matval, NULL, NULL);
	if ( status ) goto TERMINATE;

			   //adding the constraint set (4) for the first patient
	for(int i=0;i<numz;i++)	{
		matbeg[i]=i;
		obj[i]=1;		
		sense [i]='L';		
		indices[i]=(k+1)*numz+i;
	}
	
	status = CPXaddrows (env_slp, exten, 0, numz,numz, v[k],sense, matbeg, indices, obj, NULL, NULL);
	if ( status ) goto TERMINATE;

	}

		 
	for(int k=0;k<num_patients;k++){
		///////constraint set (1)

		///////constraint set (3)
		for (int i=0;i<numz;i++)				
			rhs[i]=-rhs_stop[k][i]+d[k][i];
		
		for(int i=0;i<numz;i++)
		status = CPXchgcoef (env_slp, exten, (4*k+1)*numz+i,i, rhs[i]);
				
		 ///////constraint set (2)
		for (int i=0;i<numz;i++)		
			rhs[i]=-rhs_stop[k][i]+d[k][i];
		
		for(int i=0;i<numz;i++)
		status = CPXchgcoef (env_slp, exten, (4*k+2)*numz+i, i, rhs[i]);

				///////constraint set (4)
		for (int i=0;i<numz;i++)				
			rhs[i]=v[k][i]-rhs_stop[k][i];
		
		for(int i=0;i<numz;i++)
		status = CPXchgcoef (env_slp, exten, (4*k+3)*numz+i, i, rhs[i]);
	}
	
	
	
	TERMINATE:
	delete [] matval;
	delete [] matind;
	delete [] matbeg;
	delete [] matcnt;
	delete [] sense;
	delete [] obj;
	delete [] rhs;
	delete [] lb;
	delete [] ub;
	delete []indices;

	return 0;
	}

