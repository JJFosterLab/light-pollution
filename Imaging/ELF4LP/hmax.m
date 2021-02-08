function hm = hmax(angles, rmsid, epsilon, xstep, chains, symAx)
%find the half max of a rotational image difference function
%% set up variables
    crv = fit(angles',rmsid,'smoothingspline');
    if ~exist('rmsid','var')
        warning('Please specify a response variable to fit')
    end
    if ~exist('epsilon','var')
        epsilon = 10^-5; %10e-5 takes ≈1/50 seconds
    end
    if ~exist('xstep','var')
        xstep = min(diff(angles))*10^-1; %smallest x step to take
    end
    if ~exist('chains','var')
        chains = 20; %test multiple start points
    end
    if ~exist('symAx','var')
        symAx = 0; %look either side of 0
    end
%% set up search region
%     x1 = -min(diff(angles));%    -epsilon*10;
%     x2 = min(diff(angles));%   epsilon*10;
    %Search from _maximum_ difference not miniumum
    left_x = angles(angles<symAx); right_x = angles(angles>symAx);
    left_rmsid = rmsid(angles<symAx); right_rmsid = rmsid(angles>symAx);
    left_max_x = left_x(left_rmsid == max(left_rmsid));
    right_max_x = right_x(right_rmsid == max(right_rmsid));
    x1 = left_max_x+xstep;%    -just right of negative maxima
    x2 = right_max_x-xstep;%    -just left of positive maxima
    %starting points for chains randomly distributed 
%     x1start = x1 + rand(1,chains)*abs(x1-symAx);
%     x2start = x2 - rand(1,chains)*abs(x2-symAx);
%     x1start = x1 + rand(1,4)*abs(x1-symAx)/chains;
%     x2start = x2 + rand(1,4)*abs(x2-symAx)/chains;
%% set up function to minimise
%     fx = 10^16;%starting value should be big
    inputs = struct();
    inputs.curve = crv;
    inputs.rms_imdiff = rmsid;
    f_of_x = @(xx, inputs) abs(inputs.curve(xx)./max(inputs.rms_imdiff) - 0.5);
    %use a built in optimiser, Newton-Rhaphson is not precise enough
    fo = optimset('tolX',xstep,'tolFun',epsilon, 'FunValCheck','on');
    xx1 = nan(chains,1); out1 = nan(chains,1); exit1 = zeros(chains,1);
    xx2 = nan(chains,1); out2 = nan(chains,1); exit2 = zeros(chains,1);
%% search for local minima in multiple chains
found_flag1 = 0;found_flag2 = 0;
thresh = 0.001;
while (~found_flag1 || ~found_flag2) %can this  be faster?
    if  ~found_flag1 
        x1start = -180*rand(1,chains);
        for (i = 1:chains)
    %         [xx1(i),out1(i)] = fminunc(@(xx) f_of_x(xx, inputs), x1start(i), fo);%left of 0°
    %         [xx2(i),out2(i)] = fminunc(@(xx) f_of_x(xx, inputs), x2start(i), fo);%right of 0°
            [xx1(i),out1(i),exit1(i)] = fminbnd(@(xx) f_of_x(xx, inputs),...
                x1start(i), symAx, fo);%left of 0°
                if (min(out1) < thresh)%are any of the answers close?
                    found_flag1 = 1;
                end
        end
    end
    if  ~found_flag2 
        x2start = 180*rand(1,chains);
        for (i = 1:chains)
            [xx2(i),out2(i),exit2(i)] = fminbnd(@(xx) f_of_x(xx, inputs),...
                symAx, x2start(i), fo);%right of 0°
                if (min(out2) < thresh)%are any of the answers close?
                    found_flag2 = 1;
                end
        end
    end        
end
%     [xx1,out1,exit1]
%     [xx2,out2,exit2]
%% choose the best values from chain output
%     xx1=xx1(exit1);out1=out1(exit1);%this wasn't helpful
%     xx2=xx2(exit2);out2=out2(exit2);%this wasn't helpful
    xx1=xx1(out1<thresh);out1=out1(out1<thresh);
    xx2=xx2(out2<thresh);out2=out2(out2<thresh);
    x1 = xx1(out1 == min(out1));
    x2 = xx2(out2 == min(out2));
    if(length(x1) > 1)
        warning('Multiple values for left half-max found, taking median')
        x1 = median(x1);
    end
    if(length(x2) > 1)
        warning('Multiple values for right half-max found, taking median')
        x2 = median(x2);
    end
%     if (x1 < x2)
%         warning('Halfmax order reversed')
%         xtmp = x1;
%         x1 = x2;
%         x2 = xtmp;
%     end
    %even better, use GlobalSearch 
    %if you can afford the global optimation toolbox
%                 gs = GlobalSearch('FunctionTolerance');
%                 problem1 = createOptimProblem('fminunc','x1',x1,...
%                     'objective',@(x1) f_of_x(x1, inputs),...
%                     'lb',[-180],'ub',[0]);
%                 x1 = run(gs,problem1);
%     x1 = fminsearch(@(x1) f_of_x(x1, inputs), x1, fo);%left of 0°
%     x2 = fminsearch(@(x2) f_of_x(x2, inputs), x2, fo);%right of 0°
    %defunct Newton-Raphson version
%     while fx > epsilon
%     fx =  f_of_x(x2);
%     fx_der = (f_of_x(x2+epsilon) -  f_of_x(x2) )/epsilon;
%         if(fx > epsilon)
%             x2 = x2-fx/fx_der;
%         end
%     end
%     fx = 10^16;
%     while fx > epsilon
%         fx =  f_of_x(x1);
%         fx_der = (f_of_x(x1+epsilon) -  f_of_x(x1) )/epsilon;
%         if(fx > epsilon)
%             x1 = x1-fx/fx_der;
%         end
%     end
    hm = [x1,x2];
end