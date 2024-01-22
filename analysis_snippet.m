%% analysis snippet

load('Data/ET-pilot-full-fixed-20231102T143522-0')
% make sure psignifit installed and on path

%% reconfigure response data
% 90 is forward, 270 is backward
pcorrect = NaN(length(pa.speed),length(pa.direction));
data = NaN(length(pa.speed), 3, length(pa.direction));
for speed = 1:length(pa.speed)
    idx_speed = find(pa.fullFactorial(:,3) == pa.speed(speed));
    for direction = 1:length(pa.direction)
        idx_dir = find(pa.fullFactorial(:,4) == pa.direction(direction));
        idx = intersect(idx_speed, idx_dir);
        pcorrect(speed,direction) = sum(eq(pa.LR(idx), pa.LRresponse(idx)))/pa.nRepeats;
        data(speed, :, direction) = [pa.speed(speed), sum(eq(pa.LR(idx), pa.LRresponse(idx))), pa.nRepeats];
    end
end

%% fit to psychometric function
% options for psignifit
options             = struct;   % initialize as an empty struct
options.sigmoidName = 'weibull';   
options.expType     = '2AFC';   % choose 2-AFC as the paradigm of the experiment
                                % this sets the guessing rate to .5 and
                                % fits the rest of the parameters
options.fixedPars = NaN(5,1);                                
options.fixedPars(5) = 0;       % fix eta (dispersion) to zero

% reshape matrix to add distance to constraint/surround variable in first
% column
C=permute(data,[1 3 2]);
C = reshape(C,[],size(data,2),1);

% calculate distance to constraint and distance to the surround
% last argument determines the length of the constraint, based on range of
% depth (m)
[dconst, dsurr] = DistanceToConstraint(ds, pa, 0.05);

% add distance to the constraint (dconst) or distance to surround (dsurr)
% to the first column of data matrix to indicate "stimulus intensity"
a = dconst;
C(:,1) = a(:);

% plot with psychometric curve fit
result = psignifit(C,options);
figure, plotPsych(result, options);
thresh = exp(result.Fit(1));

