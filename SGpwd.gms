$TITLE sccenario generation
* VaR_CVaR.gms: Value at Risk and Conditional Value at Risk models.
$eolcom //
option optcr=0, reslim=120;


SET Asset Underlying asset;

SET Date Time periods;

ALIAS(Asset,i);
ALIAS(Date,t);
Parameter
    IndexReturns(i,t) 'Index return';

$GDXIN Data
$LOAD Asset, Date, IndexReturns
$GDXIN

display i,t, IndexReturns;

    
set
       
         scen /s1*s250/
         w    /w1*w4/ 
         p /p1*p186/
         PWD(p,w,t)
         PD(p,t)
;


// outcomment next line, if you want to reseed the random number generator at each new run of the file
*execseed=gmillisec(jnow)

SCALAR
        randnum
        start
        end;

Parameter RetScenWeeks(i,scen,w);
Parameter AllRetScen(p, i, scen);
Parameter RollingScen(t);
scalar counter;
counter=336;

loop(p,
    loop(w,
        loop(t$(ORD(t)=counter),
                PD(p,t)=yes;
                PWD(p,w,t)=YES;
        );
        counter=counter+1;
    );

);
display PWD,PD;
*$exit


loop(p,
    start = 1+4*(ord(p)-1);
    end = 336+4*(ord(p)-1);
    loop(scen,
        loop(w,

                     randnum = uniformint(start, end);
                     RetScenWeeks(i, scen,  w) = SUM(t$(ord(t)=randnum), IndexReturns(i, t)) ;
            );
    );
    AllRetScen(p, i,scen) = prod(w, (1 + RetScenWeeks(i, scen,w ))) - 1;
);


EXECUTE_UNLOAD 'RScen.gdx',p,w,t, scen, Date,PWD, Asset, AllRetScen, IndexReturns;
