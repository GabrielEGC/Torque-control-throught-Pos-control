syms I1 I2 m1 m2 l1 l2 lc2 lc1 g
syms q1 q2 dq1 dq2
syms kp kd t
s1=sin(q1);s2=sin(q2);s12=sin(q1+q2);c2=cos(q2);
q=[q1;q2];dq=[dq1;dq2];

M=[I1+I2+m2*l1^2+2*m2*l1*lc2*c2, I2+m2*l1*lc2*c2;I2+m2*l1*lc2*c2, I2];
C=[-2*m2*l1*lc2*s2*dq2, -m2*l1*lc2*s2*dq2; m2*l1*lc2*s2*dq2, 0];
Tg=[-m1*g*lc1*s1-m2*g*(l1*s1+lc2*s12);-m2*g*lc2*s12];
phiqdq=C*dq+Tg;

xw=l1*sin(q1)+l2*sin(q1+q2);
dxwdq=transpose(jacobian(xw,q));
Hxw=jacobian(dxwdq,q);

yw=-l1*cos(q1)-l2*cos(q1+q2);
dywdq=transpose(jacobian(yw,q));
dyw=transpose(dywdq)*dq;
Hyw=jacobian(dywdq,q);

ywd=-0.6+0.1*sin(2*t);
dywd=0.2*cos(2*t);
ddywd=-0.4*sin(2*t);
fxd=5; %(?)

A=transpose(dxwdq)*(M^-1)*dxwdq;
A2=transpose(dxwdq)*(M^-1);
varphi=transpose(dq)*Hxw*dq-transpose(dxwdq)*(M^-1)*phiqdq;

psi1=transpose(dq)*Hyw*dq+transpose(dywdq)*(M^-1)*(-dxwdq*(A^-1)*varphi-phiqdq);
psi2=transpose(dywdq)*(M^-1)*(eye(2)-dxwdq*(A^-1)*A2);



rho1=[psi1;-(A^-1)*varphi];
rho2=[psi2;-(A^-1)*A2];
v=[ddywd-kd*(dyw-dywd)-kp*(yw-ywd);fxd];

u=(rho2^-1)*(v-rho1);%[0;0];%
%u=min(max(u,-39),40);
ddq=(M^-1)*(u-phiqdq+dxwdq*fxd);
dxdt=[dq;ddq];