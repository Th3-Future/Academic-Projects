function dX=pcrtbp_dyn_stm(~,X)

%mass ratio
mu=0.01213;

%extract states
x=X(1,1);
y=X(2,1);
dx=X(3,1);
dy=X(4,1);

%extract columns of state transition matrix
X1=X(5:8,1);
X2=X(9:12,1);
X3=X(13:16,1);
X4=X(17:20,1);

%distances to primaries
rs1=sqrt((x+mu)^2+y^2);
rs2=sqrt((x-1+mu)^2+y^2);

%position from primaries
rs1v=[x+mu;y];
rs2v=[x+mu-1;y];


%linearization about trajectory
omx=[0 -1;1 0];

F=[zeros(2,2) eye(2);
    -omx*omx-(1-mu)*(eye(2)-3*(rs1v*rs1v.')/rs1^2)/rs1^3-mu*(eye(2)-3*(rs2v*rs2v.')/rs2^2)/rs2^3 -2*omx];

%dynamics

dX=[dx;
    dy;
    2*dy+x-(1-mu)*(x+mu)/(rs1^3)-mu*(x+mu-1)/(rs2^3);
    -2*dx+y-(1-mu)*y/rs1^3-mu*y/rs2^3;
    F*X1;
    F*X2;
    F*X3;
    F*X4];

end
