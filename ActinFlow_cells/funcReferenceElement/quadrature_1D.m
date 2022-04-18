function [z, w] = quadrature_1D(ngaus)
%
% [z, w] = quadrature_1D(ngaus)

if ngaus == 1
    z = 0;
    w = 2;
elseif ngaus == 2
    pos1 = 1/sqrt(3);
    z = [-pos1; pos1];
    w = [ 1 1 ];
elseif ngaus == 3
    pos1 = sqrt(3/5);
    z = [-pos1; 0; pos1];
    w = [ 5/9   8/9   5/9 ];
elseif ngaus == 4
    pos1 = sqrt(525+70*sqrt(30))/35;
    pos2 = sqrt(525-70*sqrt(30))/35;
    z = [-pos1; -pos2; pos2; pos1];
    w1 = sqrt(30)*(3*sqrt(30)-5)/180;
    w2 = sqrt(30)*(3*sqrt(30)+5)/180;
    w = [w1   w2   w2   w1];
elseif ngaus == 5
    r70 = sqrt(70);
    pos1 = sqrt(245+14*r70)/21;
    pos2 = sqrt(245-14*r70)/21;
    z = [-pos1; - pos2; 0; pos2; pos1];
    w1 = (7+5*r70)*3*r70/(100*(35+2*r70));
    w2 = -(-7+5*r70)*3*r70/(100*(-35+2*r70));
    w0 = 128/225;
    w = [w1,w2,w0,w2,w1];
elseif ngaus == 6
    z = [0.23861918608319690863050172168066;     0.66120938646626451366139959501991;
         0.93246951420315202781230155449406;    -0.23861918608319690863050172168066;
        -0.66120938646626451366139959501991;    -0.93246951420315202781230155449406];
    w = [0.46791393457269104738987034398801      0.36076157304813860756983351383812, ...
         0.17132449237917034504029614217260      0.46791393457269104738987034398891, ...
         0.36076157304813860756983351383816      0.17132449237917034504029614217271 ];
elseif ngaus == 7
    z = [0.
         0.40584515137739716690660641207692;     0.74153118559939443986386477328078
         0.94910791234275852452618968404784;    -0.40584515137739716690660641207692
        -0.74153118559939443986386477328078;    -0.94910791234275852452618968404784 ];
    w = [0.4179591836734693877551020408166, ...
         0.38183005050511894495036977548841      0.27970539148927666790146777142377, ...
         0.12948496616886969327061143267904      0.38183005050511894495036977548964, ...
         0.27970539148927666790146777142336      0.12948496616886969327061143267912 ];
elseif ngaus == 8
    z = [0.18343464249564980493947614236027;     0.52553240991632898581773904918921
         0.79666647741362673959155393647586;     0.96028985649753623168356086856950
        -0.18343464249564980493947614236027;    -0.52553240991632898581773904918921
        -0.79666647741362673959155393647586;    -0.96028985649753623168356086856950];
    w = [0.3626837833783619829651504492780       0.31370664587788728733796220198797, ...
         0.22238103445337447054435599442573      0.10122853629037625915253135431028, ...
         0.3626837833783619829651504492834       0.31370664587788728733796220198807, ...
         0.22238103445337447054435599442632      0.10122853629037625915253135431015];
elseif ngaus == 9
    z = [0.
         .32425342340380892903853801464336;     .61337143270059039730870203934149
         .83603110732663579429942978806972;     .96816023950762608983557620290365
        -.32425342340380892903853801464336;    -.61337143270059039730870203934149
        -.83603110732663579429942978806972;    -.96816023950762608983557620290365];
    w = [.3302393550012597631645250692903;
         .3123470770400028400686304065887;      .2606106964029354623187428694188
         .18064816069485740405847203124168;    0.81274388361574411971892158110806e-1
        .3123470770400028400686304065836;       .2606106964029354623187428694150
        .18064816069485740405847203124263;     0.81274388361574411971892158110938e-1 ]';
elseif ngaus == 10
    z = [.14887433898163121088482600112972;      .43339539412924719079926594316579
         .67940956829902440623432736511487;      .86506336668898451073209668842349
         .97390652851717172007796401208445;     -.14887433898163121088482600112972
        -.43339539412924719079926594316579;     -.67940956829902440623432736511487
        -.86506336668898451073209668842349;     -.97390652851717172007796401208445];
    w = [.2955242247147528701738929946601;       .26926671930999635509122692156867
         .21908636251598204399553493422796;      .14945134915058059314577633966048
        0.66671344308688137593568809894898e-1;   .2955242247147528701738929946484
        .26926671930999635509122692157323;       .21908636251598204399553493422877
        .14945134915058059314577633965578;      0.66671344308688137593568809893298e-1]';
else
    error('unavailable quadrature')
end


