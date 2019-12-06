$TITLE Value at Risk and Conditional Value at Risk models
* VaR_CVaR.gms: Value at Risk and Conditional Value at Risk models.
$eolcom //
*option optcr=0, reslim=120;

SET
         scen
         Asset
         p;
         
ALIAS(Asset,i);
ALIAS(scen, s);

PARAMETER
         RetScen(i, s)
         AllRetScen(p,i,s)
         MuTarget(p)
         CVaRTarget(p)
         ERMAX(p)
         CVaRMAX(p);


$GDXIN RScen
$LOAD p, scen, Asset, AllRetScen
$GDXIN

$GDXIN Targets
$LOAD MuTarget,CVaRTarget
$GDXIN

SCALARS
              
        alpha         'Confidence level'
        Lambda        'Risk aversion parameter'
        CVaRLim       'CVaR limit'
        ExpRetLim     'Expected Return Limit'
        Budget ;


alpha  = 0.95;

Budget = 100000;

CVaRLim = Budget*10;

ExpRetLim = -100000;

RetScen(i, s) = SUM(p$(ord(p) eq 1),AllRetScen(p,i,s));


PARAMETERS
        pr(s)       'Scenario probability'
        FP(i,s)      'Final values'
        EP(i)       'Expected final values'
;

pr(s) = 1.0 / CARD(s);

FP(i,s) = 1 + RetScen( i, s );

EP(i) = SUM(s, pr(s) * FP(i,s));


POSITIVE VARIABLES
         
         x(i)            'Holdings of assets in monetary units (not proportions)'
         VaRDev(s)       'Measures of the deviation from the VaR'
         
         x_0(i)          'Holdings for the flat cost regime'
         x_1(i)          'Holdings for the linear cost regime'
         
         buy(i)
         sell(i)
         InitHold(i);

BINARY VARIABLE
    Y(i)          'Indicator variable for assets included in the portfolio'

    Yb(i) Indicator variable for assets to be purchased
    Ys(i) Indicator variable for assets to be sold;



SCALARS
  FlatCost 'Normalized fixed transaction cost' / 20 /
  PropCost 'Normalized proportional transaction cost' / 0.001 /;
  

VARIABLES
         losses(s)       'The scenario loss function'
         VaR             'The alpha Value-at-Risk'
         CVaR            'The alpha Conditional Value-at-Risk'
         ExpectedReturn  'Expected return of the portfolio'
         obj             'objective function value'
         PortReturn;

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
         
         BalanceCon(i)
         BinBuyLimits(i)
         BinSellLimits(i)
         InCon
         
;



*--Objective------

BudgetCon ..             sum(i, x(i)) =E= Budget;
InCon ..                 sum(i, x(i)) =E= SUM(i,InitHold.l(i));



ReturnCon ..             ExpectedReturn =E= sum(i, EP(i)*x(i));

LossDefCon(s) ..         Losses(s) =E= -1*sum(i, FP(i, s)*x(i) );

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


BalanceCon(i)..           x(i) - buy(i) + sell(i) =e=  x.l(i) ;

BinBuyLimits(i)..         buy(i)  =l=  Yb(i);

BinSellLimits(i)..        sell(i)  =l=  Ys(i);




MODEL CVaRModel 'The Conditional Value-at-Risk Model' /BudgetCon, ReturnCon, LossDefCon, VaRDevCon,CVaRDefCon, ObjectivFunc, CVaRLimCon, ReturnLimCon
 HoldingCon,ReturnDefWithCost,FlatCostBounds,LinCostBounds/;

//Let's maximize expected return with the target CVaR from the equal weight portfolio
Lambda = 0;
CVaRLim = SUM(p$(ord(p) eq 1),CVaRtarget(p));
ExpRetLim = -SUM(p$(ord(p) eq 1),MuTarget(p))
SOLVE CVaRModel Maximizing OBJ Using RMIP;
display X.l;
 
ERMAX(p)$(ord(p) eq 1) = ExpectedReturn.l;
CVaRMAX(p)$(ord(p) eq 1) = CVaR.l;


Parameters

    PortValue(s)
    Revision(p,i)
     
;


*--s.t.---------------------------------------------------


      

MODEL RevisionCVaRModel 'The Conditional Value-at-Risk Model' /BalanceCon,InCon,ReturnCon,LossDefCon,VaRDevCon,CVaRDefCon,
ObjectivFunc,CVaRLimCon,ReturnLimCon,HoldingCon,ReturnDefWithCost,FlatCostBounds,LinCostBounds /;
//Retscen(i,s) = SUM(p$(ord(p) eq 2),AllRetScen(p,i,s) );
//InitHold.l(i) = x.l(i)*SUM(s, pr(s) * (SUM(p$(ord(p) eq 2),AllRetScen(p,i,s) ) +1));

//ExpRetLim = -SUM(i,InitHold.l(i));
//CVaRLim = SUM(p$(ord(p) eq 2), CVaRtarget(p));

   

//Lambda = 0;
//SOLVE RevisionCVaRModel Maximizing OBJ Using RMIP;
//DISPLAY sell.l,buy.l;


loop(p$(ord(p)>1 and ord(p)<3),   
    Retscen(i,s) = AllRetScen(p,i,s);
    InitHold.l(i) = x.l(i)*SUM(s, pr(s) * (AllRetScen(p-1,i,s)+1));
*    Budget = SUM(i,InitHold.l(i));

    ExpRetLim = -SUM(i,InitHold.l(i));
    CVaRLim = CVaRtarget(p); 

    Lambda = 0.001;
    SOLVE RevisionCVaRModel Maximizing OBJ Using RMIP;
    ERMAX(p) = ExpectedReturn.l;
    CVaRMAX(p) = CVaR.l;
    DISPLAY sell.l,buy.l,x.l;
    );
    


//EXECUTE_UNLOAD 'CVaR_MAX03.gdx', ERMAX,CVaRMAX;














