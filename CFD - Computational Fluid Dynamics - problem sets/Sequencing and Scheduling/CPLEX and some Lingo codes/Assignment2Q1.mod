/*********************************************
 * OPL 20.1.0.0 Model
 * Author: ezeor
 * Creation Date: 1 Jul 2021 at 14:48:31
 *********************************************/

/*We exclude worker 5 from this optimization since we know that packaging(28sec) and
  checkout(10sec) must be done last, which both sum up to 38sec, which is not optimal
  so we consider only the first 4 workers*/
int workers=4;
int tasks=6;
int min2hr=60*60;
range n=1..tasks;
range m=1..workers;

int p[n]=[9,9,21,6,6,22];
int add[m]=[18,0,0,28]; /*compulsory added task that must be performed in worker 1 and 4*/

dvar boolean x[n][m];
dvar int+ cycle_time;
dvar int+ Wcyc_time[n][m];
minimize cycle_time;

subject to {
  c1: forall(i in n) sum(j in m)x[i][j]==1;
  c2: forall(i in n,j in m) Wcyc_time[i][j]==x[i][j]*p[i];
  c3: forall(j in m) cycle_time>=sum(i in n)Wcyc_time[i][j]+add[j];
}

float cap=round(min2hr/cycle_time); /*Since each Burritos must be a complete product*/
execute
{
  writeln("Number of Workers = ",workers+1);
  writeln("Minimum Cycle time = ",cycle_time);
  write("Maximum Capacity = ", cap, " Burritos/hr");  
}