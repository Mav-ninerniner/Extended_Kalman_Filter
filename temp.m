Noise=[1 1];
UL_ob=4.1;
UL_pre=4.2;
P_old=[1 1];
SOC_pre=1;
Up_pre=1.2;
Q = Noise(1);   
R = Noise(2); 
n=1;  
alp=0.04;
beta=2;                                                       
lamda=1.5;
%% ---------------Weight Prediction---------------
 Wm=[lamda/(n+lamda),0.5/(n+lamda)+zeros(1,2*n)]  ;                            
 Wc=Wm;
 Wc(1)=Wc(1)+(1-alp^2+beta);

%% --------------Initialization--------------     
X_pre = [SOC_pre-0.009; Up_pre];
Xsigma=X_pre;
pk=sqrt((n+lamda)*P_old);

%% ---------------Calculate Sigma Points------------
sigma1 = Xsigma+pk;
sigma2 = Xsigma-pk;
sigma=[Xsigma sigma1 sigma2];

%% -------------State Prediction----------------
sxk=0;     
for ks=1:2*n+1
    sxk=Wm(ks)*sigma(ks)+sxk;      
end

spk=0;    
for kp=1:2*n+1
    spk=Wc(kp)*(sigma(kp)-sxk)*(sigma(kp)-sxk)'+spk;  
end
spk=spk+Q;

%% ------Measurement Update--------
syk=0;     
for ky=1:2*n+1
    syk=syk+Wm(ky)*UL_pre;     
end

pyy=0;     
for kpy=1:2*n+1
    pyy=Wc(kpy)*(UL_pre-syk)*(UL_pre-syk)'+pyy;    
end
pyy=pyy+R;

%% --------------Cross-correlation between sigma points --------------
pxy=0;     
for kxy=1:2*n+1
    pxy=Wc(kxy)*(sigma(kxy)-sxk)*(UL_pre-syk)'+pxy;  
end

%% ------------Kalman gain-----------
kgs=pxy/pyy; 

%% ------------State and error covariance matrix update------
X_upd=X_pre+kgs*(UL_ob-syk);    
P_update=spk-kgs*pyy*kgs';          

%% --------------Output--------------
SOC_upd = X_upd;
P_upd = P_update;
