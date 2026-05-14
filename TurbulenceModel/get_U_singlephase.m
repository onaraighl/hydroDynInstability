function [zG,uG,Re_star,ReG,tauG,kG]=get_U_singlephase(ReG)

nn=2;

kappa=0.4;
S_guess=ReG/sqrt(2)
S=fzero(@ff_int,S_guess)

Re_star=S;
% logzG=-12:0.01:log(1);
% zG=exp(logzG);

zG=0:0.001:1;

for i=1:length(zG)
    zz=zG(i);
    if(zz==0)
        uG(i)=0;
    else
        uG(i)=quad(@(x) vel_G(x),0,zG(i));
    end
    duG(i)=vel_G(zG(i));
end

Re_av=ReG*(sum(uG)/length(uG))

Avd=exp(6.3)*(S^(-2));
psi1=1-exp(-zG.^nn/Avd);
psi2=1-exp(-(1-zG).^nn/Avd);
gg0=zG.*(1-zG);
rr0=1/abs(1-(ReG*ReG/(S*S)));
gg1=(zG.*zG.*zG+(rr0^(5/2))*((1-zG).^3))./(rr0*rr0*(1-zG).*(1-zG)+rr0*zG.*(1-zG)+zG.*zG);
gg=gg0.*gg1;

ct=0.55;
kG=(1/(ct^2))*psi1.*psi2*((S/ReG))^2;
tauG=kappa*(S/ReG)*gg.*psi1.*psi2.*duG;

% logzG=-12:0.01:log(1);
% zG=exp(logzG);
% 
% for i=1:length(zG)
%     zz=zG(i);
%     if(zz==0)
%         uG(i)=0;
%     else
%         uG(i)=quad(@(x) vel_G(x),0,zG(i));
%     end
% end
% 
% dv=1/S;
% zwall=zG/dv;
% uwall=uG*(ReG/S);
% ulin=zwall;
% 
% for i=1:length(zwall)
%     if(uwall(i)>=0.95*ulin(i))
%         xx=zwall(i);
%     end
% end
% 
% xx



%**************************************************************************

function y=ff_int(S)
    yG1=quad(@(x) vel_G1(x,S),0,1);
    y=yG1;
end

%**************************************************************************

function y=vel_G1(x,S)
    
    Avd=exp(6.3)*(S^(-2));
    
    rr=1/abs(1-(ReG*ReG/(S*S)));
    
    gg0=x.*(1-x);
    gg1=(x.*x.*x+(rr^(5/2))*((1-x).^3))./(rr*rr*(1-x).*(1-x)+rr*x.*(1-x)+x.*x);
    gg=gg0.*gg1;
    
    psi1=1-exp(-x.^nn/Avd);
    psi2=1-exp(-(1-x).^nn/Avd);
    y=(S*S/ReG)*((1-(ReG^2/S^2)*x)./(1+kappa*(S/sqrt(rr))*gg.*psi1.*psi2));
    
end

%**************************************************************************

function y=vel_G(x)
    
    Avd=exp(6.3)*(S^(-2));
    
    rr=1/abs(1-(ReG*ReG/(S*S)));
    
    gg0=x.*(1-x);
    gg1=(x.*x.*x+(rr^(5/2))*((1-x).^3))./(rr*rr*(1-x).*(1-x)+rr*x.*(1-x)+x.*x);
    gg=gg0.*gg1;
    
    psi1=1-exp(-x.^nn/Avd);
    psi2=1-exp(-(1-x).^nn/Avd);
    y=(S*S/ReG)*(((1-(ReG^2/S^2)*x))./(1+kappa*(S/sqrt(rr))*gg.*psi1.*psi2));
end

%**************************************************************************

end