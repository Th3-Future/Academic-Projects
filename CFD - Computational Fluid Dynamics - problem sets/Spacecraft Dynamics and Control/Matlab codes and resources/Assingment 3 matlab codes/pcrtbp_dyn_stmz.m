function dX=pcrtbp_dyn_stmz(~,X)

%mass ratio
mu=0.01213;

%extract states
x=X(1,1);
y=X(2,1);
dx=X(3,1);
dy=X(4,1);
z=0;
%extract columns of state transition matrix
X1=X(5:6,1);
X2=X(7:8,1);


%distances to primaries
rs1=sqrt((x+mu)^2+y^2+z^2);
rs2=sqrt((x-1+mu)^2+y^2+z^2);

%position from primaries
rs1v=[x+mu;y;z];
rs2v=[x+mu-1;y;z];


%linearization of z-equation about trajectory


F=[0 1;-(1-mu)/rs1^3-mu/rs2^3 0];
%dynamics

dX=[dx;
    dy;
    2*dy+x-(1-mu)*(x+mu)/(rs1^3)-mu*(x+mu-1)/(rs2^3);
    -2*dx+y-(1-mu)*y/rs1^3-mu*y/rs2^3;
    F*X1;
    F*X2];

end
