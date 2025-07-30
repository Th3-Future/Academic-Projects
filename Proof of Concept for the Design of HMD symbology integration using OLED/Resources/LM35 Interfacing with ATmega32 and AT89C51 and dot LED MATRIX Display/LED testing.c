#include<reg51.h>
sbit l1=P2^0;
sbit l2=P2^1;
sbit l3=P2^2;
sbit s1=P3^0;
sbit s2=P3^1;
sbit s3=P3^2;

void main()
{
while(1)
{
	P2=0x00;
if(s1==0)
	{l1=1;
		while(s1==0);
	}
if(s2==0)
	{l2=1;
		while(s2==0);	
	}
if(s3==0)
	{l3=1;
		while(s3==0);
	}	
	l1=l2=l3=0;
}
}