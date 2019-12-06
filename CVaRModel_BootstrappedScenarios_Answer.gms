$TITLE Value at Risk and Conditional Value at Risk models
* VaR_CVaR.gms: Value at Risk and Conditional Value at Risk models.
$eolcom //
option optcr=0, reslim=120;

SET
         scen
         Asset
;
ALIAS(Asset,i);
ALIAS(scen, s);

PARAMETER
         RetScen(i, s)
;

$GDXIN Scenarios01
$LOAD Asset, scen, RetScen
$GDXIN


SCALARS
        Budget        'Nominal investment budget'
        alpha         'Confidence level'
        Lambda        'Risk aversion parameter'
        CVaRLim       'CVaR limit'
        ExpRetLim     'Expected Return Limit'
;

Budget = 100.0;
alpha  = 0.95;
CVaRLim = Budget*10;
ExpRetLim = -100;

PARAMETERS
        pr(s)       'Scenario probability'
        P(i,s)      'Final values'
        EP(i)       'Expected final values'
;

pr(s) = 1.0 / CARD(s);


P(i,s) = 1 + RetScen( i, s );


EP(i) = SUM(s, pr(s) * P(i,s));


POSITIVE VARIABLES
         x(i)            'Holdings of assets in monetary units (not proportions)'
         VaRDev(s)       'Measures of the deviation from the VaR'
;

VARIABLES
         losses(s)       'The scenario loss function'
         VaR             'The alpha Value-at-Risk'
         CVaR            'The alpha Conditional Value-at-Risk'
         ExpectedReturn  'Expected return of the portfolio'
         obj             'objective function value'
;



EQUATIONS
         BudgetCon       'Equation defining the budget constraint'
         ReturnCon       'Equation defining the portfolio expected return'
         LossDefCon(s)   'Equation defining the losses'
         VaRDevCon(s)    'Equation defining the VaRDev variable'
         CVaRDefCon      'Equation defining the CVaR'
         ObjectivFunc    'lambda formulation of the MeanCVaR model'
         CVaRLimCon      'Constraint limiting the CVaR'
         ReturnLimCon    'Constraint on a minimum expected return'
;



*--Objective------

*--s.t.-----------
BudgetCon ..             sum(i, x(i)) =E= Budget;

ReturnCon ..             ExpectedReturn =E= sum(i, EP(i)*x(i));

LossDefCon(s) ..         Losses(s) =E= -1*sum(i, P(i, s)*x(i) );

VaRDevCon(s) ..          VaRDev(s) =G= Losses(s) - VaR;

CVaRDefCon ..            CVaR =E= VaR + (sum(s, pr(s)*VarDev(s) ) )/(1 - alpha);

ObjectivFunc ..          Obj =E= (1-lambda)*ExpectedReturn - lambda*CVaR;

CVaRLimCon ..            CVaR =L= CVaRLim;

ReturnLimCon ..          ExpectedReturn =G= ExpRetLim;



*--Models-----------



//Let's build an equal weight portfolio first,
//by fixing the X values to be eqaully weighted:
X.fx(i) = Budget/card(i);
display X.l;


MODEL CVaRModel 'The Conditional Value-at-Risk Model' /BudgetCon, ReturnCon, LossDefCon, VaRDevCon,CVaRDefCon, ObjectivFunc, CVaRLimCon, ReturnLimCon/;



*------------CVaR----------------------
//We need both terms in the objective function to be under control in order for the model to calculate mean return and CVaR values correctly:
lambda = 0.999;
SOLVE CVaRModel Maximizing OBJ Using LP;


DISPLAY X.l, ExpectedReturn.l, VaR.L, CVaR.L;


Parameters
    MuTarget,
    CVaRTarget,
    BestCase,
    WorstCase
;

//Now we use the expected return and the CVaR of the equal weight strategy as target benchmarks in the CVaR model
MuTarget = ExpectedReturn.l;
CVaRTarget = CVaR.L;
display MuTarget, CVaRTarget;

parameter PortValue(s);
PortValue(s) = sum(i, P(i, s)*x.l(i) );

parameter SummaryReport(*,*);
SummaryReport(s,'PortValue') = PortValue(s);

BestCase = Smax(s, PortValue(s));
WorstCase = Smin(s, PortValue(s));
display BestCase, WorstCase;



*EXECUTE_UNLOAD 'SummaryReport.gdx', SummaryReport;
*EXECUTE 'gdxxrw.exe SummaryReport.gdx O=PortValue.xls par=SummaryReport rng=sheet1!a1' ;

//The next two lines are used to free the X variable again
X.lo(i) = 0;
X.up(i) = Budget*10;


*------------------Eff Front with equidistant CVaR--------------------*
set pp /pp0*pp10/;

PARAMETERS
         RunningCVaR(pp)          'Optimal level of portfolio CVaR'
         RunningReturn(pp)        'Portfolio return'
         RunningAllocation(pp,i)  'Optimal asset allocation'
         MaxCVar
         MinCVar
;

*first we find the biggest possible average, hence maximum CVaR:
lambda = 0.001;
SOLVE CVaRModel Maximizing OBJ Using LP;
MaxCVar = CVaR.l;
RunningCVaR(pp)$(ord(pp)=1) = CVaR.l;
RunningReturn(pp)$(ord(pp)=1)  = ExpectedReturn.l;
RunningAllocation(pp,i)$(ord(pp)=1)     = x.l(i);
display ExpectedReturn.l;
*Then we find the lowest possible variance:
lambda = 0.999;
SOLVE CVaRModel Maximizing OBJ Using LP;
MinCVar = CVaR.l;
RunningCVaR(pp)$(ord(pp)=card(pp)) = CVaR.l;
RunningReturn(pp)$(ord(pp)=card(pp))  = ExpectedReturn.l;
RunningAllocation(pp,i)$(ord(pp)=card(pp))     = x.l(i);
display ExpectedReturn.l;
display MaxCVar, MinCVar;
$exit
*Then we find the equidistant variances in between
lambda = 0;
CVarLim = MaxCVar;
loop(pp$(ord(pp)>1 and ord(pp)<card(pp)),
CVarLim = CVarLim - (MaxCVar-MinCVar)/(card(pp)-1);
SOLVE CVaRModel Maximizing OBJ Using LP;
*display VarLim, PortVariance.l;
RunningCVaR(pp) = CVaR.l;
RunningReturn(pp)  = ExpectedReturn.l;
RunningAllocation(pp,i)     = x.l(i);
);

display RunningCVaR, RunningReturn, RunningAllocation;


parameter SummaryReport2(*,*);
* Store results by rows
SummaryReport2(i,pp) = RunningAllocation(pp,i);
SummaryReport2('CVaR',pp) = RunningCVaR(pp);
SummaryReport2('Return',pp) = RunningReturn(pp);


DISPLAY SummaryReport2;
* Write SummaryReport into an Excel file

EXECUTE_UNLOAD 'Summary2.gdx', SummaryReport2;
EXECUTE 'gdxxrw.exe Summary2.gdx O=MeanCVaRFrontier.xls par=SummaryReport2 rng=sheet1!a1' ;


//Let's minimize CVaR with the target return from the equal weight portfolio
ExpRetLim = Mutarget;

Lambda = 1;
CVaRLim = CVaRLim+100;
SOLVE CVaRModel Maximizing OBJ Using LP;
display X.l, ExpectedReturn.l, CVaR.l;


//Let's maximize expected return with the target CVaR from the equal weight portfolio
CVaRLim = CVaRtarget;


Lambda = 0;
ExpRetLim = -100;
SOLVE CVaRModel Maximizing OBJ Using LP;
display X.l, ExpectedReturn.l, CVaR.l;


















