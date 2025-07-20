/*********************************************
 * OPL 20.1.0.0 Model
 * Author: Godswill Ezeorah
 * Creation Date: 28 Jun 2021 at 12:07:53
 *********************************************/
 int n=4;
int mac=2;
int M=1000;
range job=1..n;
range stage=1..mac;
range stage2=2..mac;
int p[job][stage]=[[2,3],[4,5],[6,2],[5,3]];
int d[job]=[5,10,12,10];
dvar boolean x[job][job];
dvar float+ c[job][stage];
dvar float+ T[job][stage];
dvar float+ TT;
minimize  TT;
subject to {
 c1: forall(j in job) c[j][1]>=p[j][1];
 c2: forall(j in job, i2 in stage2,i in stage) c[j][i2]>=c[j][i]+p[j][i];
 c3: forall(i in stage, j in job, k in job )
 	if (j<k)
 	c[j][i]>=c[k][i]+p[j][i]-M*(1-x[j][k]);
 c4: forall(i in stage,j in job, k in job)
 	if (j<k)
 	c[k][i]>=c[j][i]+p[k][i]-M*(x[j][k]);
 c5: forall(j in job, i in stage) T[j][i]>=c[j][i]-d[j];
 c6: TT==sum(j in job, i in stage) T[j][i];
}