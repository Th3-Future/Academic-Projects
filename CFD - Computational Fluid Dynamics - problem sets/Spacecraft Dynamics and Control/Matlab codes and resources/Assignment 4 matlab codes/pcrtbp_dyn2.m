function dX=pcrtbp_dyn2(~,X)

%mass ratio
mu=0.01213;

%extract states
x=X(1,1);
y=X(2,1);
dx=X(3,1);
dy=X(4,1);

%distances to primaries
rs1=sqrt((x+mu)^2+y^2);
rs2=sqrt((x-1+mu)^2+y^2);

%dynamics

dX=-[dx;
    dy;
    2*dy+x-(1-mu)*(x+mu)/(rs1^3)-mu*(x+mu-1)/(rs2^3);
    -2*dx+y-(1-mu)*y/rs1^3-mu*y/rs2^3];

end
