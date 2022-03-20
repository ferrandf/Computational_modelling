function Fn=numericalNormalFluxPureConvection(uL,uR,velocityDotNormal)

if velocityDotNormal>0
    Fn=velocityDotNormal*uL;
else
    Fn=velocityDotNormal*uR;
end

