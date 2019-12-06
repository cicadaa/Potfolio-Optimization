$TITLE Value at Risk and Conditional Value at Risk models
* VaR_CVaR.gms: Value at Risk and Conditional Value at Risk models.
$eolcom //
option optcr=0, reslim=120;

SET
         scen
         Asset
         w
         PWD
         PD
         Date
         p
;
ALIAS(Asset,i);
ALIAS(scen, s);
ALIAS(Date, t);


PARAMETER
         RetScen(i, s)
         AllRetScen(p,i,s)
        
         
;

$GDXIN RScen
$LOAD Date,PWD,PD,p,w, scen, Asset, AllRetScen
$GDXIN


*$exit
SCALARS
        Budget        'Nominal investment budget'
        alpha         'Confidence level'
        Lambda        'Risk aversion parameter'
        CVaRLim       'CVaR limit'
        ExpRetLim     'Expected Return Limit'
;

Budget = 100000.0;
alpha  = 0.95;
CVaRLim = Budget*10;
ExpRetLim = -100000;

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
         x_0(i)          'Holdings for the flat cost regime'
         x_1(i)          'Holdings for the linear cost regime'
         buy(i)          'Amount to buy'
         sell(i)         'Amount to be sold';
         
BINARY VARIABLE
    Y(i)          'Indicator variable for assets included in the portfolio';

PARAMETER
    xlow(i)    'lower bound for active variables' ;


SCALARS
  FlatCost 'Normalized fixed transaction cost' / 20 /
  PropCost 'Normalized proportional transaction cost' / 0.001 /;
  

VARIABLES
         losses(s)       'The scenario loss function'
         VaR             'The alpha Value-at-Risk'
         CVaR            'The alpha Conditional Value-at-Risk'
         ExpectedReturn  'Expected return of the portfolio'
         obj             'objective function value'
         PortReturn
;


xlow(i) = 0.0;
x.UP(i) = 100000.0;//incase short sell
x_0.UP(i) = 20000;

EQUATIONS
         BudgetCon       'Equation defining the budget constraint'
         ReturnCon       'Equation defining the portfolio expected return'
         LossDefCon(s)   'Equation defining the losses'
         VaRDevCon(s)    'Equation defining the VaRDev variable'
         CVaRDefCon      'Equation defining the CVaR'
         ObjectivFunc    'lambda formulation of the MeanCVaR model'
         CVaRLimCon      'Constraint limiting the CVaR'
         ReturnLimCon    'Constraint on a minimum expected return'
         
         HoldingCon(i)        'Constraint defining the holdings'
         ReturnDefWithCost    'Equation defining the portfolio return with cost'
         FlatCostBounds(i)    'Upper bounds for flat transaction fee'
         LinCostBounds(i)     'Upper bonds for linear transaction fee'
         
         BalanceCon(i)    'Same sell and buy amount';
         



BudgetCon ..             sum(i, x(i)) =E= Budget;

ReturnCon ..             ExpectedReturn =E= sum(i, EP(i)*x(i));

LossDefCon(s) ..         Losses(s) =E= -1*sum(i, P(i, s)*x(i) );

VaRDevCon(s) ..          VaRDev(s) =G= Losses(s) - VaR;

CVaRDefCon ..            CVaR =E= VaR + (sum(s, pr(s)*VarDev(s) ) )/(1 - alpha);

ReturnDefWithCost..      PortReturn =e= SUM(i, ( EP(i)*x_0(i) - FlatCost*Y(i) ) ) +
                          SUM(i, (EP(i) - PropCost)*x_1(i));

ObjectivFunc ..          Obj =E= (1-lambda)*PortReturn - lambda*CVaR;

CVaRLimCon ..            CVaR =L= CVaRLim;

ReturnLimCon ..          ExpectedReturn =G= ExpRetLim;


HoldingCon(i)..           x(i) =e= x_0(i) + x_1(i); 

FlatCostBounds(i)..       x_0(i) =l= x_0.UP(i) * Y(i);

LinCostBounds(i)..        x_1(i) =l= Y(i);


BalanceCon(i)..          x(i) - buy(i) + sell(i) =e=  InitHold(i);




*--Models-----------



//Let's build an equal weight portfolio first,
//by fixing the X values to be eqaully weighted:
X.fx(i) = Budget/card(i);

MODEL CVaRModel 'The Conditional Value-at-Risk Model' /BudgetCon, ReturnCon, LossDefCon, VaRDevCon,CVaRDefCon, ObjectivFunc, CVaRLimCon, ReturnLimCon
 HoldingCon,ReturnDefWithCost,FlatCostBounds,LinCostBounds/;



*------------CVaR----------------------
//We need both terms in the objective function to be under control in order for the model to calculate mean return and CVaR values correctly:
lambda = 0.999;
SOLVE CVaRModel Maximizing OBJ Using MIP;


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

//The next two lines are used to free the X variable again
X.lo(i) = 0;
X.up(i) = Budget*10;


//Let's maximize expected return with the target CVaR from the equal weight portfolio
CVaRLim = CVaRtarget;
Lambda = 0;
ExpRetLim = -100000;
SOLVE CVaRModel Maximizing OBJ Using MIP; 
display X.l, ExpectedReturn.l, CVaR.l,x_0.l;
$exit

*init holding=================================================================================
Parameter
        InitHold(i) 'initial hold of assets';



MODEL RevisionCVaRModel 'The Conditional Value-at-Risk Model' /BalanceCon, ReturnCon, LossDefCon, VaRDevCon,CVaRDefCon, ObjectivFunc, CVaRLimCon, ReturnLimCon
 HoldingCon,ReturnDefWithCost,FlatCostBounds,LinCostBounds/;


loop(p$ord(p)>1,
    RetScen( i, s ) = AllRetScen(p,i,s);
    initHold(i) = X.l * AllRetScen(p-1，i，s);
    CVaRLim = SUM(i,InitHold(i))*10;
*max-----------------------
    Lambda = 1;
    CVaRLim = CVaRLim+100000;
    SOLVE RevisionCVaRModel Maximizing OBJ Using MIP;
    x(p)ExpectedReturn.l, CVaR.l;
*min-----------------------
    
    
);
















