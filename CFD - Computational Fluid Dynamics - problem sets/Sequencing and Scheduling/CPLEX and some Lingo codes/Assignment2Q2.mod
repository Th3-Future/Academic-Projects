/*********************************************
 * OPL 20.1.0.0 Model
 * Author: ezeor
 * Creation Date: 1 Jul 2021 at 16:10:28
 *********************************************/

/*For this optimization we reduced the number of workers to 3, since the last 
  (i.e., 4th worker) will have a cycle time of 38sec for packaging(28sec) 
  and checkout(10sec)*/
int workers=3;
int tasks=6;
int min2hr=60*60;
int worker4=38;

range n=1..tasks;
range m=1..workers;
range cnt=1..2;
int p[n]=[9,9,21,6,6,22];
int add[m]=[18,0,0]; /*compulsory added task that must be performed in worker 1*/

dvar boolean x[n][m];
dvar int+ cycle_time;
dvar int+ Wcyc_time[n][m];
minimize cycle_time;

subject to 
{
  c1: forall(i in n) sum(j in m)x[i][j]==1;
  c2: forall(i in n,j in m) Wcyc_time[i][j]==x[i][j]*p[i];
  c3: forall(j in m) cycle_time>=sum(i in n)Wcyc_time[i][j]+add[j];
}

int cyc_time[cnt]=[cycle_time,worker4];
/*Below line, compares the optimized cycle time for the 1st three stations and station 4*/
int actual_cyc_time=max(i in cnt)cyc_time[i]; 
float cap=round(min2hr/actual_cyc_time); /*Since each Burritos must be a complete product*/
execute
{
  writeln("Number of Workers = ",workers+1);
  writeln("Minimum Cycle time = ",actual_cyc_time);
  write("Maximum Capacity = ", cap, " Burritos/hr");  
}