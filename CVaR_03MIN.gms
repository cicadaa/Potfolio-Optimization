$TITLE Value at Risk and Conditional Value at Risk models
* VaR_CVaR.gms: Value at Risk and Conditional Value at Risk models.
$eolcom //
*option optcr=0, reslim=120;

SET
         scen
         Asset
         Date
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




*display IndexReturns;

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
PARAMETER
         RetScen(i, s)
         AllRetScen(p,i,s)
         MuTarget(p)
         CVaRTarget(p)
         ERMIN(p)
         CVaRMIN(p)
         AverageCase(p);


$GDXIN Targets_new
$LOAD MuTarget,CVaRTarget,AverageCase
$GDXIN


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

SCALARS
              
        alpha         'Confidence level'
        Lambda        'Risk aversion parameter'
        CVaRLim       'CVaR limit'
        ExpRetLim     'Expected Return Limit'
        Budget ;


alpha  = 0.95;

Budget = 100000;

CVaRLim = Budget*10;

ExpRetLim = 100000;

RetScen(i, s) = SUM(p$(ord(p) eq 1),AllRetScen(p,i,s));


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
         
         buy(i)
         sell(i)
         ;

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
*         InCon
         CVaRLimCon01
         HoldingCon01(i)
         SellBuy
;



*--Objective------

BudgetCon ..             sum(i, x(i)) =E= Budget;
*InCon ..                sum(i, x(i)) =E= SUM(i,InitHold(i));
ReturnCon ..             ExpectedReturn =E= sum(i, EP(i)*x(i));

LossDefCon(s) ..         Losses(s) =E= -1*sum(i, FP(i, s)*x(i) );

VaRDevCon(s) ..          VaRDev(s) =G= Losses(s) - VaR;

CVaRDefCon ..            CVaR =E= VaR + (sum(s, pr(s)*VarDev(s) ) )/(1 - alpha);

ReturnDefWithCost..      PortReturn =e= SUM(i, ( EP(i)*x_0(i) - FlatCost*Y(i) ) ) +
                          SUM(i, (EP(i) - PropCost)*x_1(i));

ObjectivFunc ..          Obj =E= (1-lambda)*PortReturn - lambda*CVaR;

CVaRLimCon ..             CVaR =L= CVaRLim;

ReturnLimCon ..           ExpectedReturn =G= ExpRetLim;


HoldingCon(i)..           x(i) =e= x_0(i) + x_1(i);
*HoldingCon01(i)..         x(i) =e= x_0(i) + x_1(i); 

FlatCostBounds(i)..       x_0(i) =l= x_0.UP(i) * Y(i);

LinCostBounds(i)..        x_1(i) =l= Y(i);


BalanceCon(i)..           x(i) - buy(i) + sell(i) =e=  InitHold(i) ;

SellBuy..                 SUM(i,buy(i)) =E=  SUM(i,sell(i));

BinBuyLimits(i)..         buy(i)  =l=  Yb(i);

BinSellLimits(i)..        sell(i)  =l=  Ys(i);




MODEL CVaRModel 'The Conditional Value-at-Risk Model' /BudgetCon,ReturnCon, LossDefCon, VaRDevCon,CVaRDefCon, ObjectivFunc, ReturnLimCon,HoldingCon,ReturnDefWithCost,FlatCostBounds,LinCostBounds/;

//Let's maximize expected return with the target CVaR from the equal weight portfolio

Lambda = 0.99;

*ExpRetLim = SUM(p$(ord(p) eq 1),AverageCase(p));
ExpRetLim = SUM(p$(ord(p) eq 1),MuTarget(p))*Budget;

*Budget;
*SOLVE CVaRModel Minimizing CVaR Using MIP;
SOLVE CVaRModel Maximizing OBJ Using MIP;
ERMIN(p)$(ord(p) eq 1) = ExpectedReturn.l;
CVaRMIN(p)$(ord(p) eq 1) = CVaR.l;
DISPLAY x.l,ExpRetLim;
$exit
Parameters

    PortValue(s)
    Revision(p,i)
    Tot(p) 
;

*--s.t.---------------------------------------------------

MODEL RevisionCVaRModel 'The Conditional Value-at-Risk Model' /BalanceCon,SellBuy, ReturnCon,LossDefCon,VaRDevCon,CVaRDefCon,
ObjectivFunc,ReturnLimCon,HoldingCon,ReturnDefWithCost,FlatCostBounds,LinCostBounds /;
//Retscen(i,s) = SUM(p$(ord(p) eq 2),AllRetScen(p,i,s) );
//InitHold.l(i) = x.l(i)*SUM(s, pr(s) * (SUM(p$(ord(p) eq 2),AllRetScen(p,i,s) ) +1));

//ExpRetLim = -SUM(p$(ord(p) eq 2), MuTarget(p));
//CVaRLim = SUM(p$(ord(p) eq 2), CVaRtarget(p));

   

//Lambda = 0;
//SOLVE RevisionCVaRModel Maximizing OBJ Using MIP;



loop(p$(ord(p)>1),
*and ord(p)<20),
    pr(s) = 1.0 / CARD(s);
    Retscen(i,s) = AllRetScen(p,i,s);
    InitHold(i) = x.l(i)*PROD(t$PD(p,t),(1+IndexReturns(i,t)));
         
    Tot(p) = SUM(i,InitHold(i))/SUM(i,x.l(i)); 

    ExpRetLim = MuTarget(p)*SUM(i,InitHold(i));
    FP(i,s) = 1 + RetScen( i, s );
    EP(i) = SUM(s, pr(s) * FP(i,s));

  
    Lambda = 0.99;
*    SOLVE RevisionCVaRModel Minimizing CVaR Using MIP;
    SOLVE RevisionCVaRModel Maximizing OBJ Using MIP;
    ERMIN(p) = ExpectedReturn.l;
    CVaRMIN(p) = CVaR.l;
    DISPLAY sell.l,buy.l,x.l,fp;
    );
DISPLAY Tot;
SCALAR
muti;
muti = PROD(p$(ord(p)>1 and ord(p)<186),Tot(p));
DISPLAY muti; 
EXECUTE_UNLOAD 'CVaR_MIN_d.gdx', ERMIN,CVaRMIN;














