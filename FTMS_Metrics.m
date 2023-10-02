function varargout = FTMS_Metrics(filename,varargin)
%% Description 

% This code takes a formula list and calculates various number- and abundance-weighed 
% metrics. They are exported in an .xlsx file. The code can be integrated
% in other codes.


% Example: 
% FTMS_Metrics('Sample 1_Final.xlsx')

%% Copyright and License Notice: 

% Copyright © 2022 Old Dominion University Research Foundation, Norfolk VA, USA
% All rights reserved.

% This file is part of the Toolbox for Environmental Research (TEnvR). Please cite the toolbox as follows: 
% Goranov, A. I., Sleighter, R. L., Yordanov, D. A., and Hatcher, P. G. (2023): 
% TEnvR: MATLAB-Based Toolbox for Environmental Research, Analytical Methods, doi: 10.1039/d3ay00750b.

% TEnvR is free software for non-commercial use: you can redistribute it and/or modify 
% %it under the terms of the GNU General Public License as published by the Free Software Foundation, 
% either version 3 of the License, or (at your option) any later version. 
% Users wishing to use TEnvR commercially must obtain a commercial license from Old Dominion University Research Foundation. 

% TEnvR is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
% of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

% You should have received a copy of the GNU General Public License along with TEnvR
% (located in the Supplementary files directory of the toolbox). If not, see <https://www.gnu.org/licenses/>.

% Please contact the corresponding authors Drs. Aleksandar Goranov (aleksandar.i.goranov@gmail.com) 
% and Patrick Hatcher (phatcher@odu.edu) with any questions or concerns.

%% Code
format=FTMS_ConfigurationToolbox;

if nargin==1 %FTMS_Metrics(filename)
    copyfile(filename, [filename(1:end-5) '_Metrics.xlsx'],'f')
    Data=xlsread(filename);
    Int=Data(:,format.Column_Magnitude);
    RelInt=Int./sum(Int);
    Data=[Data,RelInt];
elseif nargin==2 %FTMS_Metrics(filename,sheet)
    copyfile(filename, [filename(1:end-5) '_Metrics.xlsx'],'f')
    sheet=varargin{1};
    Data=xlsread(filename,sheet);
    Int=Data(:,format.Column_Magnitude);
    RelInt=Int./sum(Int);
    Data=[Data,RelInt];
elseif nargin==4 %FTMS_Metrics(filename,sheet,loc1,loc2)
    copyfile(filename, [filename(1:end-5) '_Metrics.xlsx'],'f')
    sheet=varargin{1};
    location1=varargin{2};
    location2=varargin{3};
    Data=xlsread(filename,sheet);
    Int=Data(:,format.Column_Magnitude);
    RelInt=Int./sum(Int);
    Data=[Data,RelInt];

% For codes with relative magnitude already calculated
elseif nargin==5 % FTMS_Metrics(filename,sheet,loc1,loc2,true) 
    sheet=varargin{1};
    location1=varargin{2};
    location2=varargin{3};
    Data=xlsread(filename,sheet);
    RelInt=Data(:,format.Column_RelMagnitude);
    Data=[Data,RelInt];
else
    error('Error! You have entered too many variables for this function!')
end

C=Data(:,format.Column_C);
H=Data(:,format.Column_H);
O=Data(:,format.Column_O);
N=Data(:,format.Column_N);
S=Data(:,format.Column_S);
P=Data(:,format.Column_P);
Cl=Data(:,format.Column_E);
Mass=Data(:,format.Column_ExactMass);
OC=Data(:,format.Column_OC);
HC=Data(:,format.Column_HC);
DBE=Data(:,format.Column_DBE);      % Already corrected for presence of E
DBEC=Data(:,format.Column_DBEC);    % Already corrected for presence of E
AImod=Data(:,format.Column_AImod);  % Already corrected for presence of E

% Finding condensed molecules (AImod >=0.67, C >=15)
BlackCarbon=Data(AImod >= 0.67 & C >= 15,format.Column_mz:format.Column_RelMagnitude);
BlackCarbon_DBE=BlackCarbon(:,format.Column_DBE);
BlackCarbon_RelInt=BlackCarbon(:,format.Column_RelMagnitude);

% Calculating parameters
NC=N./C;    
ClC=Cl./C;
HN=H./N;
ON=O./N;
HCl=H./Cl;
OCl=O./Cl;
NCl=N./Cl;
HS=H./S;
HP=H./P;
OS=O./S;
OP=O./P;
NS=N./S;
PS=P./S;
DBEH=DBE./H;
DBEO=DBE./O;
DBEminusO=DBE-O;
Ring=(BlackCarbon_DBE-3.1374)/2.436;

if format.Heteroelement_Halogen
    NOSC=4-(((4*C)+H-(3*N)-(2*O)-(2*S)+(5*P)-Cl)./C);
elseif format.Heteroelement_N15
	NOSC=4-(((4*C)+H-(3*N)-(3*E)-(2*O)-(2*S)+(5*P))./C);
else
	NOSC=4-(((4*C)+H-(3*N)-(2*O)-(2*S)+(5*P))./C);
end

% Intensity-weighing all parameters
NCw=NC.*RelInt;               
ClCw=ClC.*RelInt;
HNw=HN.*RelInt;
ONw=ON.*RelInt;
HClw=HCl.*RelInt;
OClw=OCl.*RelInt;
NClw=NCl.*RelInt;
HSw=HS.*RelInt;
HPw=HP.*RelInt;
OSw=OS.*RelInt;
OPw=OP.*RelInt;
NSw=NS.*RelInt;
PSw=PS.*RelInt;
DBEHw=DBEH.*RelInt;
DBEOw=DBEO.*RelInt;
DBEminusOw=DBEminusO.*RelInt;
Ringw=Ring.*BlackCarbon_RelInt;
NOSCw=NOSC.*RelInt;
Cw=C.*RelInt;
Hw=H.*RelInt;
Ow=O.*RelInt;
Nw=N.*RelInt;
Clw=Cl.*RelInt;
Sw=S.*RelInt;
Pw=P.*RelInt;
Massw=Mass.*RelInt;
OCw=OC.*RelInt;
HCw=HC.*RelInt;
DBEw=DBE.*RelInt;
DBECw=DBEC.*RelInt;
AImodw=AImod.*RelInt;

% Averages of parameters
C_AVG=mean(nonzeros(C(isfinite(C))));
H_AVG=mean(nonzeros(H(isfinite(H))));
O_AVG=mean(nonzeros(O(isfinite(O))));
N_AVG=mean(nonzeros(N(isfinite(N))));
Cl_AVG=mean(nonzeros(Cl(isfinite(Cl))));
S_AVG=mean(nonzeros(S(isfinite(S))));
P_AVG=mean(nonzeros(P(isfinite(P))));
OC_AVG=mean(nonzeros(OC(isfinite(OC))));
HC_AVG=mean(nonzeros(HC(isfinite(HC))));
NC_AVG=mean(nonzeros(NC(isfinite(NC))));
ClC_AVG=mean(nonzeros(ClC(isfinite(ClC))));
HN_AVG=mean(nonzeros(HN(isfinite(HN))));
ON_AVG=mean(nonzeros(ON(isfinite(ON))));
HCl_AVG=mean(nonzeros(HCl(isfinite(HCl))));
OCl_AVG=mean(nonzeros(OCl(isfinite(OCl))));
NCl_AVG=mean(nonzeros(NCl(isfinite(NCl))));
HS_AVG=mean(nonzeros(HS(isfinite(HS))));
HP_AVG=mean(nonzeros(HP(isfinite(HP))));
OS_AVG=mean(nonzeros(OS(isfinite(OS))));
OP_AVG=mean(nonzeros(OP(isfinite(OP))));
NS_AVG=mean(nonzeros(NS(isfinite(NS))));
PS_AVG=mean(nonzeros(PS(isfinite(PS))));

Mass_AVG=mean(nonzeros(Mass(isfinite(Mass))));
DBE_AVG=mean(nonzeros(DBE(isfinite(DBE))));
DBEC_AVG=mean(nonzeros(DBEC(isfinite(DBEC))));
DBEH_AVG=mean(nonzeros(DBEH(isfinite(DBEH))));
DBEO_AVG=mean(nonzeros(DBEO(isfinite(DBEO))));
DBEminusO_AVG=mean(nonzeros(DBEminusO(isfinite(DBEminusO))));
AImod_AVG=mean(nonzeros(AImod(isfinite(AImod))));
Ring_AVG=mean(nonzeros(Ring(isfinite(Ring))));
NOSC_AVG=mean(nonzeros(NOSC(isfinite(NOSC))));

% Standard Deviations
C_STD=std(nonzeros(C(isfinite(C))));
H_STD=std(nonzeros(H(isfinite(H))));
O_STD=std(nonzeros(O(isfinite(O))));
N_STD=std(nonzeros(N(isfinite(N))));
Cl_STD=std(nonzeros(Cl(isfinite(Cl))));
S_STD=std(nonzeros(S(isfinite(S))));
P_STD=std(nonzeros(P(isfinite(P))));
OC_STD=std(nonzeros(OC(isfinite(OC))));
HC_STD=std(nonzeros(HC(isfinite(HC))));
NC_STD=std(nonzeros(NC(isfinite(NC))));
ClC_STD=std(nonzeros(ClC(isfinite(ClC))));
HN_STD=std(nonzeros(HN(isfinite(HN))));
ON_STD=std(nonzeros(ON(isfinite(ON))));
HCl_STD=std(nonzeros(HCl(isfinite(HCl))));
OCl_STD=std(nonzeros(OCl(isfinite(OCl))));
NCl_STD=std(nonzeros(NCl(isfinite(NCl))));
HS_STD=std(nonzeros(HS(isfinite(HS))));
HP_STD=std(nonzeros(HP(isfinite(HP))));
OS_STD=std(nonzeros(OS(isfinite(OS))));
OP_STD=std(nonzeros(OP(isfinite(OP))));
NS_STD=std(nonzeros(NS(isfinite(NS))));
PS_STD=std(nonzeros(PS(isfinite(PS))));

Mass_STD=std(nonzeros(Mass(isfinite(Mass))));
DBE_STD=std(nonzeros(DBE(isfinite(DBE))));
DBEC_STD=std(nonzeros(DBEC(isfinite(DBEC))));
DBEH_STD=std(nonzeros(DBEH(isfinite(DBEH))));
DBEO_STD=std(nonzeros(DBEO(isfinite(DBEO))));
DBEminusO_STD=std(nonzeros(DBEminusO(isfinite(DBEminusO))));
AImod_STD=std(nonzeros(AImod(isfinite(AImod))));
Ring_STD=std(nonzeros(Ring(isfinite(Ring))));
NOSC_STD=std(nonzeros(NOSC(isfinite(NOSC))));

% Averages of intensity-weighed parameters
Cw_AVG=sum(Cw(isfinite(Cw)));
Hw_AVG=sum(Hw(isfinite(Hw)));
Ow_AVG=sum(Ow(isfinite(Ow)));
Nw_AVG=sum(Nw(isfinite(Nw)));
Clw_AVG=sum(Clw(isfinite(Clw)));
Sw_AVG=mean(Sw(isfinite(Sw)));
Pw_AVG=mean(Pw(isfinite(Pw)));
OCw_AVG=sum(OCw(isfinite(OCw)));
HCw_AVG=sum(HCw(isfinite(HCw)));
NCw_AVG=sum(NCw(isfinite(NCw)));
ClCw_AVG=sum(ClCw(isfinite(ClCw)));
HNw_AVG=sum(HNw(isfinite(HNw)));
ONw_AVG=sum(ONw(isfinite(ONw)));
HClw_AVG=sum(HClw(isfinite(HClw)));
OClw_AVG=sum(OClw(isfinite(OClw)));
NClw_AVG=sum(NClw(isfinite(NClw)));
HSw_AVG=sum(HSw(isfinite(HSw)));
HPw_AVG=sum(HPw(isfinite(HPw)));
OSw_AVG=sum(OSw(isfinite(OSw)));
OPw_AVG=sum(OPw(isfinite(OPw)));
NSw_AVG=sum(NSw(isfinite(NSw)));
PSw_AVG=sum(PSw(isfinite(PSw)));

Massw_AVG=sum(Massw(isfinite(Massw)));
DBEw_AVG=sum(DBEw(isfinite(DBEw)));
DBECw_AVG=sum(DBECw(isfinite(DBECw)));
DBEHw_AVG=sum(DBEHw(isfinite(DBEHw)));
DBEOw_AVG=sum(DBEOw(isfinite(DBEOw)));
DBEminusOw_AVG=sum(DBEminusOw(isfinite(DBEminusOw)));
AImodw_AVG=sum(AImodw(isfinite(AImodw)));
Ringw_AVG=sum(Ringw(isfinite(Ringw)));
NOSCw_AVG=sum(NOSCw(isfinite(NOSCw)));

titleAVGs={'C' 'H' 'O' 'N' 'S' 'P' 'Cl' 'O/C' 'H/C' 'N/C' 'Cl/C' 'H/N' 'O/N' 'H/S' 'H/P' 'O/S' 'O/P' 'N/S' 'P/S' 'H/E' 'O/E' 'N/E' 'ExactMass' 'DBE' 'DBE/C' 'DBE/H' 'DBE/O' 'DBE-O' 'AImod' 'Ring' 'NOSC'};
AVGs=[C_AVG,H_AVG,O_AVG,N_AVG,S_AVG,P_AVG,Cl_AVG,OC_AVG,HC_AVG,NC_AVG,ClC_AVG,HN_AVG,ON_AVG,HS_AVG,HP_AVG,OS_AVG,OP_AVG,NS_AVG,PS_AVG,HCl_AVG,OCl_AVG,NCl_AVG,Mass_AVG,DBE_AVG,DBEC_AVG,DBEH_AVG,DBEO_AVG,DBEminusO_AVG,AImod_AVG,Ring_AVG,NOSC_AVG];
STDs=[C_STD,H_STD,O_STD,N_STD,S_STD,P_STD,Cl_STD,OC_STD,HC_STD,NC_STD,ClC_STD,HN_STD,ON_STD,HS_STD,HP_STD,OS_STD,OP_STD,NS_STD,PS_STD,HCl_STD,OCl_STD,NCl_STD,Mass_STD,DBE_STD,DBEC_STD,DBEH_STD,DBEO_STD,DBEminusO_STD,AImod_STD,Ring_STD,NOSC_STD];

titleAVGsw={'Cw' 'Hw' 'Ow' 'Nw' 'Sw' 'P' 'Clw' 'O/Cw' 'H/Cw' 'N/Cw' 'Cl/Cw' 'H/Nw' 'O/Nw' 'H/Sw' 'H/Pw' 'O/Sw' 'O/Pw' 'N/Sw' 'P/Sw' 'H/Ew' 'O/Ew' 'N/Ew' 'ExactMassw' 'DBEw' 'DBE/Cw' 'DBE/Hw' 'DBE/Ow' 'DBE-Ow' 'AImodw' 'Ringw' 'NOSCw'};
AVGsw=[Cw_AVG,Hw_AVG,Ow_AVG,Nw_AVG,Sw_AVG,Pw_AVG,Clw_AVG,OCw_AVG,HCw_AVG,NCw_AVG,ClCw_AVG,HNw_AVG,ONw_AVG,HSw_AVG,HPw_AVG,OSw_AVG,OPw_AVG,NSw_AVG,PSw_AVG,HClw_AVG,OClw_AVG,NClw_AVG,Massw_AVG,DBEw_AVG,DBECw_AVG,DBEHw_AVG,DBEOw_AVG,DBEminusOw_AVG,AImodw_AVG,Ringw_AVG,NOSCw_AVG];

output1=[titleAVGs;num2cell(AVGs);num2cell(STDs)]; output1=[{' ';'AVG (num)';'STDEV'},output1];
output2=[titleAVGsw;num2cell(AVGsw)]; output2=[{' ';'AVG (magn)'},output2];
warning('off','MATLAB:xlswrite:AddSheet')

if nargout == 0
    if nargin==1
        xlswrite([filename(1:end-5) '_Metrics.xlsx'],output1,1,'U7'); 
        xlswrite([filename(1:end-5) '_Metrics.xlsx'],output2,1,'U11');
    elseif nargin==2
        xlswrite([filename(1:end-5) '_Metrics.xlsx'],output1,sheet,'U7'); 
        xlswrite([filename(1:end-5) '_Metrics.xlsx'],output2,sheet,'U11');
    elseif nargin==4
        xlswrite([filename(1:end-5) '_Metrics.xlsx'],output1,sheet,location1);
        xlswrite([filename(1:end-5) '_Metrics.xlsx'],output2,sheet,location2);
    else
        xlswrite(filename,output1,sheet,location1);
        xlswrite(filename,output2,sheet,location2);
    end
elseif nargout == 3
    varargout{1}=AVGs;
    varargout{2}=STDs;
    varargout{3}=AVGsw;
end

if length(varargin) ~=4 
    disp(['Finished calculating metrics on ' char(filename)])
end

end

% Common Errors:

% Error: No magnitude-weighted average values coming out, only number-weighed
% Caused by: You have an empty line separating your data into two, thus the 
% empty line is assigned as NaN, thus the error in weighting. 