%% load various fit types

%fit of normalized dot product of two state vectors of reaction A->B with
%A(0)=0 and B(0)=0 and time offset tOffset before reaction starts with time between
%state vectors dt

%alpha dot product of A species vector
%beta dot product of A species and B species vectors
%gamma dot product of B species vector
%tOffset time until reaction starts
%tau reaction rate constant
%dt time between the two state vectors

tempString=['beta.*exp(1).^((t+dt>tOffset).*tau.^(-1).*(dt+t+(-1).*tOffset)).*(alpha.* '...
'exp(1).^(2.*(t>tOffset).*tau.^(-1).*((-1).*t+tOffset))+(-2).*beta.*exp(1) '...
'.^((t>tOffset).*tau.^(-1).*((-1).*t+tOffset)).*((-1)+exp(1).^((t>tOffset).* '...
'tau.^(-1).*((-1).*t+tOffset)))+((-1)+exp(1).^((t>tOffset).*tau.^(-1).*(( '...
'-1).*t+tOffset))).^2.*gamma).^(-1/2).*(alpha+((-1)+exp(1).^( '...
'(t+dt>tOffset).*tau.^(-1).*(dt+t+(-1).*tOffset))).*(2.*beta+((-1)+exp(1) '...
'.^((t+dt>tOffset).*tau.^(-1).*(dt+t+(-1).*tOffset))).*gamma)).^(-1).*(exp( '...
'1).^((-2).*(t+dt>tOffset).*tau.^(-1).*(dt+t+(-1).*tOffset)).*(alpha+((-1)+ '...
'exp(1).^((t+dt>tOffset).*tau.^(-1).*(dt+t+(-1).*tOffset))).*(2.*beta+((-1) '...
'+exp(1).^((t+dt>tOffset).*tau.^(-1).*(dt+t+(-1).*tOffset))).*gamma))).^( '...
'1/2)+alpha.*exp(1).^((t+dt>tOffset).*tau.^(-1).*(dt+t+(-1).*tOffset)+ '...
'(t>tOffset).*tau.^(-1).*((-1).*t+tOffset)).*(alpha.*exp(1).^(2.*(t>tOffset).* '...
'tau.^(-1).*((-1).*t+tOffset))+(-2).*beta.*exp(1).^((t>tOffset).*tau.^(-1) '...
'.*((-1).*t+tOffset)).*((-1)+exp(1).^((t>tOffset).*tau.^(-1).*((-1).*t+ '...
'tOffset)))+((-1)+exp(1).^((t>tOffset).*tau.^(-1).*((-1).*t+tOffset))) '...
'.^2.*gamma).^(-1/2).*(alpha+((-1)+exp(1).^((t+dt>tOffset).*tau.^(-1).*(dt+ '...
't+(-1).*tOffset))).*(2.*beta+((-1)+exp(1).^((t+dt>tOffset).*tau.^(-1).*( '...
'dt+t+(-1).*tOffset))).*gamma)).^(-1).*(exp(1).^((-2).*(t+dt>tOffset).* '...
'tau.^(-1).*(dt+t+(-1).*tOffset)).*(alpha+((-1)+exp(1).^((t+dt>tOffset).* '...
'tau.^(-1).*(dt+t+(-1).*tOffset))).*(2.*beta+((-1)+exp(1).^( '...
'(t+dt>tOffset).*tau.^(-1).*(dt+t+(-1).*tOffset))).*gamma))).^(1/2)+(-2).* '...
'beta.*exp(1).^((t+dt>tOffset).*tau.^(-1).*(dt+t+(-1).*tOffset)+(t>tOffset).* '...
'tau.^(-1).*((-1).*t+tOffset)).*(alpha.*exp(1).^(2.*(t>tOffset).*tau.^(-1) '...
'.*((-1).*t+tOffset))+(-2).*beta.*exp(1).^((t>tOffset).*tau.^(-1).*((-1).* '...
't+tOffset)).*((-1)+exp(1).^((t>tOffset).*tau.^(-1).*((-1).*t+tOffset)))+( '...
'(-1)+exp(1).^((t>tOffset).*tau.^(-1).*((-1).*t+tOffset))).^2.*gamma).^( '...
'-1/2).*(alpha+((-1)+exp(1).^((t+dt>tOffset).*tau.^(-1).*(dt+t+(-1).* '...
'tOffset))).*(2.*beta+((-1)+exp(1).^((t+dt>tOffset).*tau.^(-1).*(dt+t+(-1) '...
'.*tOffset))).*gamma)).^(-1).*(exp(1).^((-2).*(t+dt>tOffset).*tau.^(-1).*( '...
'dt+t+(-1).*tOffset)).*(alpha+((-1)+exp(1).^((t+dt>tOffset).*tau.^(-1).*( '...
'dt+t+(-1).*tOffset))).*(2.*beta+((-1)+exp(1).^((t+dt>tOffset).*tau.^(-1).* '...
'(dt+t+(-1).*tOffset))).*gamma))).^(1/2)+beta.*exp(1).^(2.*(t+dt>tOffset).* '...
'tau.^(-1).*(dt+t+(-1).*tOffset)+(t>tOffset).*tau.^(-1).*((-1).*t+tOffset) '...
').*(alpha.*exp(1).^(2.*(t>tOffset).*tau.^(-1).*((-1).*t+tOffset))+(-2).* '...
'beta.*exp(1).^((t>tOffset).*tau.^(-1).*((-1).*t+tOffset)).*((-1)+exp(1) '...
'.^((t>tOffset).*tau.^(-1).*((-1).*t+tOffset)))+((-1)+exp(1).^((t>tOffset).* '...
'tau.^(-1).*((-1).*t+tOffset))).^2.*gamma).^(-1/2).*(alpha+((-1)+ '...
'exp(1).^((t+dt>tOffset).*tau.^(-1).*(dt+t+(-1).*tOffset))).*(2.*beta+((-1) '...
'+exp(1).^((t+dt>tOffset).*tau.^(-1).*(dt+t+(-1).*tOffset))).*gamma)).^(-1) '...
'.*(exp(1).^((-2).*(t+dt>tOffset).*tau.^(-1).*(dt+t+(-1).*tOffset)).*( '...
'alpha+((-1)+exp(1).^((t+dt>tOffset).*tau.^(-1).*(dt+t+(-1).*tOffset))).*( '...
'2.*beta+((-1)+exp(1).^((t+dt>tOffset).*tau.^(-1).*(dt+t+(-1).*tOffset))).* '...
'gamma))).^(1/2)+(-1).*exp(1).^((t+dt>tOffset).*tau.^(-1).*(dt+t+(-1).* '...
'tOffset)).*gamma.*(alpha.*exp(1).^(2.*(t>tOffset).*tau.^(-1).*((-1).*t+ '...
'tOffset))+(-2).*beta.*exp(1).^((t>tOffset).*tau.^(-1).*((-1).*t+tOffset)) '...
'.*((-1)+exp(1).^((t>tOffset).*tau.^(-1).*((-1).*t+tOffset)))+((-1)+exp(1) '...
'.^((t>tOffset).*tau.^(-1).*((-1).*t+tOffset))).^2.*gamma).^(-1/2).*( '...
'alpha+((-1)+exp(1).^((t+dt>tOffset).*tau.^(-1).*(dt+t+(-1).*tOffset))).*( '...
'2.*beta+((-1)+exp(1).^((t+dt>tOffset).*tau.^(-1).*(dt+t+(-1).*tOffset))).* '...
'gamma)).^(-1).*(exp(1).^((-2).*(t+dt>tOffset).*tau.^(-1).*(dt+t+(-1).* '...
'tOffset)).*(alpha+((-1)+exp(1).^((t+dt>tOffset).*tau.^(-1).*(dt+t+(-1).* '...
'tOffset))).*(2.*beta+((-1)+exp(1).^((t+dt>tOffset).*tau.^(-1).*(dt+t+(-1) '...
'.*tOffset))).*gamma))).^(1/2)+exp(1).^(2.*(t+dt>tOffset).*tau.^(-1).*(dt+ '...
't+(-1).*tOffset)).*gamma.*(alpha.*exp(1).^(2.*(t>tOffset).*tau.^(-1).*(( '...
'-1).*t+tOffset))+(-2).*beta.*exp(1).^((t>tOffset).*tau.^(-1).*((-1).*t+ '...
'tOffset)).*((-1)+exp(1).^((t>tOffset).*tau.^(-1).*((-1).*t+tOffset)))+(( '...
'-1)+exp(1).^((t>tOffset).*tau.^(-1).*((-1).*t+tOffset))).^2.*gamma).^( '...
'-1/2).*(alpha+((-1)+exp(1).^((t+dt>tOffset).*tau.^(-1).*(dt+t+(-1).* '...
'tOffset))).*(2.*beta+((-1)+exp(1).^((t+dt>tOffset).*tau.^(-1).*(dt+t+(-1) '...
'.*tOffset))).*gamma)).^(-1).*(exp(1).^((-2).*(t+dt>tOffset).*tau.^(-1).*( '...
'dt+t+(-1).*tOffset)).*(alpha+((-1)+exp(1).^((t+dt>tOffset).*tau.^(-1).*( '...
'dt+t+(-1).*tOffset))).*(2.*beta+((-1)+exp(1).^((t+dt>tOffset).*tau.^(-1).* '...
'(dt+t+(-1).*tOffset))).*gamma))).^(1/2)+exp(1).^((t+dt>tOffset).*tau.^(-1) '...
'.*(dt+t+(-1).*tOffset)+(t>tOffset).*tau.^(-1).*((-1).*t+tOffset)).* '...
'gamma.*(alpha.*exp(1).^(2.*(t>tOffset).*tau.^(-1).*((-1).*t+tOffset))+( '...
'-2).*beta.*exp(1).^((t>tOffset).*tau.^(-1).*((-1).*t+tOffset)).*((-1)+ '...
'exp(1).^((t>tOffset).*tau.^(-1).*((-1).*t+tOffset)))+((-1)+exp(1).^( '...
'(t>tOffset).*tau.^(-1).*((-1).*t+tOffset))).^2.*gamma).^(-1/2).*(alpha+(( '...
'-1)+exp(1).^((t+dt>tOffset).*tau.^(-1).*(dt+t+(-1).*tOffset))).*(2.*beta+( '...
'(-1)+exp(1).^((t+dt>tOffset).*tau.^(-1).*(dt+t+(-1).*tOffset))).*gamma)) '...
'.^(-1).*(exp(1).^((-2).*(t+dt>tOffset).*tau.^(-1).*(dt+t+(-1).*tOffset)).* '...
'(alpha+((-1)+exp(1).^((t+dt>tOffset).*tau.^(-1).*(dt+t+(-1).*tOffset))).*( '...
'2.*beta+((-1)+exp(1).^((t+dt>tOffset).*tau.^(-1).*(dt+t+(-1).*tOffset))).* '...
'gamma))).^(1/2)+(-1).*exp(1).^(2.*(t+dt>tOffset).*tau.^(-1).*(dt+t+(-1).* '...
'tOffset)+(t>tOffset).*tau.^(-1).*((-1).*t+tOffset)).*gamma.*(alpha.*exp( '...
'1).^(2.*(t>tOffset).*tau.^(-1).*((-1).*t+tOffset))+(-2).*beta.*exp(1).^( '...
'(t>tOffset).*tau.^(-1).*((-1).*t+tOffset)).*((-1)+exp(1).^((t>tOffset).*tau.^( '...
'-1).*((-1).*t+tOffset)))+((-1)+exp(1).^((t>tOffset).*tau.^(-1).*((-1).*t+ '...
'tOffset))).^2.*gamma).^(-1/2).*(alpha+((-1)+exp(1).^((t+dt>tOffset).* '...
'tau.^(-1).*(dt+t+(-1).*tOffset))).*(2.*beta+((-1)+exp(1).^( '...
'(t+dt>tOffset).*tau.^(-1).*(dt+t+(-1).*tOffset))).*gamma)).^(-1).*(exp(1) '...
'.^((-2).*(t+dt>tOffset).*tau.^(-1).*(dt+t+(-1).*tOffset)).*(alpha+((-1)+ '...
'exp(1).^((t+dt>tOffset).*tau.^(-1).*(dt+t+(-1).*tOffset))).*(2.*beta+((-1) '...
'+exp(1).^((t+dt>tOffset).*tau.^(-1).*(dt+t+(-1).*tOffset))).*gamma))).^( '...
'1/2)'];

normDotProductFit = fittype( tempString, 'independent', 't', 'dependent', 'y','problem','dt' );

%fit of mean variable x of reaction A->B with A(0)=1 and B(0)=0 and time offset tOffset before reaction starts
%x value of A is xa, x value of B is xb

%tOffset time until reaction starts
%tau reaction rate constant
%xa x value species A
%xb x value species B

meanAtoBFit = fittype( '(exp((tOffset-t).*(t>tOffset)./tau).*(xa-xb)+xb)', 'independent', 't', 'dependent', 'y' );

%fit of gaussian of width sigma convolved with step function from stepStart
%to stepEnd with height stepHeight

%sigma width of gaussian
%stepEnd end of step function
%stepHeight height of step function
%stepStart start of step function

gaussConvolveStepFit = fittype( '(1/2)*stepHeight*(-erf((stepStart-x)/(sqrt(2)*sigma))+erf((stepEnd-x)/(sqrt(2)*sigma)))', 'independent', 'x', 'dependent', 'y' );


%fit of intensity standard deviation of reaction A->B with A(0)=a0 and B(0)=0 and time offset tOffset before reaction starts
%two species vector consisting of two step functions of width na and nb

%a0 starting concentration A
%alpha y variance of background intensity
%width of A
%width of B
%total width of vector
%tOffset time until reaction starts
%tau reaction rate constant

yStandardDeviationFit = fittype( 'sqrt((-(na+nb-ntot)/(ntot.^2)+2*(a0-1)/ntot+((na-a0*ntot*exp((tOffset-x)*(x>tOffset)/tau)).^2)/(na*(ntot.^2))+((nb+a0*ntot*(exp((tOffset-x)*(x>tOffset)/tau)-1)).^2)/(nb*(ntot.^2))+((a0-1).^2)*(1+alpha*(1/(ntot-nb-na)-1)))/(ntot-1))', 'independent', 'x', 'dependent', 'y' );



