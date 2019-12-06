$TITLE V1/N strategies

$eolcom //
option optcr=0, reslim=120;

SET
         scen
         Asset
         period;         
ALIAS(Asset,i);
ALIAS(scen, s);
ALIAS(period, t);




PARAMETER
         RetScen(t, i, s)
         RetScenrol1(i, s)
         target(t);

$GDXIN RollingScenarios01
$LOAD period, scen, Asset, RetScen
$GDXIN

SCALARS
        Budget        'Nominal investment budget'
        alpha         'Confidence level'
        Lambda        'Risk aversion parameter'
        CVaRLim       'CVaR limit'
        ExpRetLim     'Expected Return Limit';

Budget = 100000.0;
alpha  = 0.95;
CVaRLim = Budget*10;
ExpRetLim = -100;

PARAMETERS
        pr(s)       'Scenario probability'
        P(i,s)      'Final values'
        EP(i)       'Expected final values';
        
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
    BudgetCon ..             sum(i, x(i)) =E= Budget;

    ReturnCon ..             ExpectedReturn =E= sum(i, EP(i)*x(i));

    LossDefCon(s) ..         Losses(s) =E= -1*sum(i, P(i, s)*x(i) );

    VaRDevCon(s) ..          VaRDev(s) =G= Losses(s) - VaR;

    CVaRDefCon ..            CVaR =E= VaR + (sum(s, pr(s)*VarDev(s) ) )/(1 - alpha);

    ObjectivFunc ..          Obj =E= (1-lambda)*ExpectedReturn - lambda*CVaR;

    CVaRLimCon ..            CVaR =L= CVaRLim;

    ReturnLimCon ..          ExpectedReturn =G= ExpRetLim;
 MODEL CVaRModel 'The Conditional Value-at-Risk Model' /BudgetCon, ReturnCon, LossDefCon,
 VaRDevCon,CVaRDefCon, ObjectivFunc, CVaRLimCon, ReturnLimCon/;
    
loop(t,
    RetScenrol1(i, s) = RetScen(t, i, s);
    pr(s) = 1.0 / CARD(s);
    P(i,s) = 1 + RetScenrol1(i, s) ;
    EP(i) = SUM(s, pr(s) * P(i,s));

    X.fx(i) = Budget/card(i);

   
    lambda = 0.999;

    SOLVE CVaRModel Maximizing OBJ Using LP;
    target(t) = CVaR.l;
    );
  


EXECUTE_UNLOAD 'Targets.gdx', target;
*EXECUTE 'gdxxrw.exe Targets.gdx O=MeanCVaRFrontier.xls par=SummaryReport2 rng=sheet1!a1' ;
