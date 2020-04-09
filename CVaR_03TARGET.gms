$TITLE Value at Risk and Conditional Value at Risk models
* VaR_CVaR.gms: Value at Risk and Conditional Value at Risk models.
$eolcom //
option optcr=0, reslim=120;

SET
         scen
         Asset
         p
         w
         t
         PWD(p,w,t)
         PD(p,t)
;
ALIAS(Asset,i);
ALIAS(scen, s);


PARAMETER
         RetScen(i, s)
         AllRetScen(p,i,s)
         IndexReturns(i,t)        
;


$GDXIN RScen
$LOAD p,w,t, scen,PWD, Asset, AllRetScen, IndexReturns
$GDXIN
*display IndexReturns;
*$exit
scalar counter;
counter=336;
loop(p,
    loop(w,
        loop(t$(ORD(t)=counter),
                PD(p,t)=yes;
               
        );
        counter=counter+1;
    );

);


RetScen(i, s) = SUM(p$(ord(p) eq 1),AllRetScen(p,i,s));

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
        FP(i,s)      'Final values'
        EP(i)       'Expected final values'
        
        InitHold(i)
;

pr(s) = 1.0 / CARD(s);


FP(i,s) = 1 + RetScen( i, s );


EP(i) = SUM(s, pr(s) * FP(i,s));


POSITIVE VARIABLES
         x(i)            'Holdings of assets in monetary units (not proportions)'
         VaRDev(s)       'Measures of the deviation from the VaR'
         x_0(i)          'Holdings for the flat cost regime'
         x_1(i)          'Holdings for the linear cost regime'
         ;

BINARY VARIABLE
    Y(i)          'Indicator variable for assets included in the portfolio';

PARAMETER
    xlow(i)    'lower bound for active variables' ;
    

// In case short sales are allowed these bounds must be set properly.
xlow(i) = 0.0;
*x.UP(i) = 100000.0;;

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
         
        ReturnCon02
        LossDefCon02(s) 
         
   

;



*--Objective------

BudgetCon ..             sum(i, x(i)) =E= Budget;

ReturnCon ..             ExpectedReturn =E= sum(i, EP(i)*x(i));
ReturnCon02 ..           ExpectedReturn =E= sum(i, EP(i)*InitHold(i));

LossDefCon(s) ..         Losses(s) =E= -1*sum(i, FP(i, s)*x(i) );
LossDefCon02(s) ..       Losses(s) =E= -1*sum(i, FP(i, s)*InitHold(i) );


VaRDevCon(s) ..          VaRDev(s) =G= Losses(s) - VaR;

CVaRDefCon ..            CVaR =E= VaR + (sum(s, pr(s)*VarDev(s) ) )/(1 - alpha);

ObjectivFunc ..          Obj =E= (1-lambda)*ExpectedReturn - lambda*CVaR;

CVaRLimCon ..            CVaR =L= CVaRLim;

ReturnLimCon ..          ExpectedReturn =G= ExpRetLim;




*--Models-----------



//Let's build an equal weight portfolio first,
//by fixing the X values to be eqaully weighted:
X.fx(i)= Budget/card(i);
display X.l;


MODEL CVaRModel 'The Conditional Value-at-Risk Model' /BudgetCon, ReturnCon, LossDefCon, VaRDevCon,CVaRDefCon, ObjectivFunc, CVaRLimCon, ReturnLimCon/;

MODEL CVaRModel02 'The Conditional Value-at-Risk Model' / ReturnCon, LossDefCon, VaRDevCon,CVaRDefCon, ObjectivFunc, CVaRLimCon, ReturnLimCon/;


*------------CVaR----------------------
//We need both terms in the objective function to be under control in order for the model to calculate mean return and CVaR values correctly:
lambda = 0.5;
SOLVE CVaRModel Maximizing OBJ Using MIP;


Parameters
        MuTarget(p),
        CVaRTarget(p),
        BestCase(p),
        WorstCase(p),
        AverageCase(p),
        PortValue(s),
        Revision(p,i),
        Tot(p)
        Totx(p) 
;
    PortValue(s) = sum(i, FP(i, s)*x.l(i) );
    BestCase(p)$(ord(p) eq 1) = Smax(s, PortValue(s));
    WorstCase(p)$(ord(p) eq 1) = Smin(s, PortValue(s));
//Now we use the expected return and the CVaR of the equal weight strategy as target benchmarks in the CVaR model
AverageCase(p)$(ord(p) eq 1) = ExpectedReturn.l;
MuTarget(p)$(ord(p) eq 1) = ExpectedReturn.l/Budget;
CVaRTarget(p)$(ord(p) eq 1) = CVaR.L/Budget;




*--s.t.---------------------------------------------------


loop(p$(ord(p)>1 ),

    pr(s) = 1.0 / CARD(s);

*   Retscen(i,s) = AllRetScen(p,i,s)$(ord(p) 1 );
    InitHold(i) = x.l(i)*PROD(t$PD(p,t),(1+IndexReturns(i,t)));
    X.fx(i)=  InitHold(i);
    FP(i,s) = 1 + RetScen( i, s );
    EP(i) = SUM(s, pr(s) * FP(i,s));


    lambda = 0.5;
    SOLVE CVaRModel02 Maximizing obj Using LP;

    PortValue(s) = sum(i, FP(i, s)*x.l(i));
    BestCase(p) = Smax(s, PortValue(s));
    WorstCase(p) = Smin(s, PortValue(s));
    
    AverageCase(p) = ExpectedReturn.l;
    MuTarget(p) = ExpectedReturn.l/SUM(i,InitHold(i));
    CVaRTarget(p) = CVaR.L/SUM(i,InitHold(i));
 

    );

EXECUTE_UNLOAD 'Targets_Fixed.gdx',p,MuTarget,CVaRTarget,BestCase,AverageCase,WorstCase;
$exit 














