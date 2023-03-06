function [dP,S_normal,Ts_mag]=FSP(tau,aphi,friction,Pp,Sv,strike_,dip_)

gamma=( sqrt(friction^2+1) +friction)^2;
tau=tau/norm(tau); 
e=eig(tau);
phi=(e(2)-e(3))/(e(1)-e(3));

    if aphi>1 && aphi<=2 
        S3= (Sv/phi - Pp + gamma*Pp) / (1/phi -1 + gamma);
        S1=(Sv-S3)/phi + S3;
        S2=Sv;
    end

    if aphi<=1  %tau(3,3)<=tau(2,2) && tau(3,3)<=tau(1,1)
        S1=Sv;
        S3=Pp+(S1-Pp)/gamma;
        S2=(1-phi)*S3+phi*S1; 
    end

    if aphi>2 %tau(3,3)>=tau(2,2) && tau(3,3)>=tau(1,1)
        S3=Sv;
        S1=Pp+gamma*(S3-Pp);
        S2=(1-phi)*S3+phi*S1; 
    end

    D=S1-S3;diag=(S1+S2+S3)/3;
    [v,d]=eig(tau);
    dd=d(3,3)-d(1,1);
    t=tau*-D/dd; 
    t(1,1)=t(1,1)+diag;
    t(2,2)=t(2,2)+diag; 
    t(3,3)=t(3,3)+diag;
    
    for i_mechanism=1:length(strike_)
        strike=strike_(i_mechanism);
        dip=dip_(i_mechanism);
        n1 = -sin(dip*pi/180).*sin(strike*pi/180);
        n2 =  sin(dip*pi/180).*cos(strike*pi/180);
        n3 = -cos(dip*pi/180);
        N=[n1;n2;n3;];
        T=-t*N;  
        S_normal(i_mechanism,1) = dot(T,N);B = cross(T,N); 
        Ts = cross(N,B); 
        Ts_mag(i_mechanism,1) = sqrt(Ts(1)^2 + Ts(2)^2 + Ts(3)^2);
    end
    
   
    S_n=S_normal;
    S_normal=S_normal+Pp; %%%% effective normal if at hydrostatic
    critP=-(S_n-Ts_mag/friction);%S_n*friction;
    dP=Pp-critP;
    
    
    
    
%     mdp=min(dP);
%     dP_all(:,i)=dP-min(dP);
%     ShearProp=100*Ts_mag./Shear_crit;
%     ShearProp_all(:,i)=ShearProp;
%     
%         LabettecritP=-(LabetteS_normal-LabetteTs_mag/friction);%S_n*friction;
%     critP=-(LabetteS_normal-LabetteTs_mag/friction);%S_n*friction;
%     LabettedP(i)=Pp-critP-mdp;