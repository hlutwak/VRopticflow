%% Analysis code
% get behavior files, conccatenate same trials, run psignifit

% add psignifit toolbox
% addpath('/Users/hopelutwak/Documents/MATLAB/psignifit')
addpath(genpath('/Users/hopelutwak/Documents/MATLAB/'))

addpath(genpath('/Applications/Psychtoolbox'))
% addpath('/Users/hopelutwak/Documents/GitHub/VRopticflow/Analysis')
addpath('/Users/hlutwak/Documents/GitHub/VRopticflow/Analysis')

% addpath(genpath('C:\Users\hlutw\Documents\MATLAB\psignifit-master'))
% assign data folder
% dataFolder = '/Users/hopelutwak/Documents/GitHub/VRopticflow/Data';
dataFolder = '/Users/hlutwak/Documents/GitHub/VRopticflow/Data';

figFolder = '/Users/hlutwak/Documents/GitHub/VRopticflow/Figures';

% analysisFolder = '/Users/hopelutwak/Documents/GitHub/VRopticflow/Analysis';
analysisFolder = '/Users/hlutwak/Documents/GitHub/VRopticflow/Analysis';

% names of files
S = dir(fullfile(dataFolder,'*.mat'));

%
% which subjects data to analyze
subjects = ["MG"]; %,"DL", "PL","MG", "SM", "IK", "JO", "KZ"]; %"MP","DL", "PL","MG", "SM", "IK", "JO", "KZ","IG"
data_set = 1; % 1 = 1.5 deg tolerance
% all: "PL", "MP", "SM", "JL", "IK", "JO", "KZ", "IG"
% all with good eyetracking trials: subjects = ["MP","DL","PL", "MG", "SM", "IK", "JO", "KZ","IG"];

stims =  ["full-1", "full-2"]; %, "monocular-2"]; %add "copy" to have pa.good_trials, and/or dconst and dsurround based on vertical eye movements
% % ["full-1", "full-2"];["dots-1", "dots-2"] ["monocular-1", "monocular-2"]
ideal_eye = 1; % use measurements of data_const and data_surr based on ideal eye movements, otherwise use eyetracking vertical movements
% depth_range = .05; % additive
depth_est = 0;
depth_range = 1.05; % multiplicative
data_const = [];
data_surr= [];

% loop over all subjects

for s  = 1:length(subjects)

    data = [];
    %load appropriate files
    count = 0;
    for f = 1:length(S)
        subj = contains(S(f).name,subjects(s));
        sti = contains(S(f).name,stims) && contains(S(f).name, 'eyetracking');
        %add "eyetrakcing" to have pa.good_trials, and/or dconst and dsurround based on vertical eye movements

        if subj && sti
            load(fullfile(dataFolder,S(f).name));
            display(fullfile(dataFolder,S(f).name));

            if ideal_eye
                n_conditions = length(pa.speed)*length(pa.direction);
                conditions = fullfact([numel(pa.speed), numel(pa.direction)]);
                data_session = [];
                for cond = 1:n_conditions
                    idx_speed = find(pa.fullFactorial(:,3) == pa.speed(conditions(cond,1)));
                    idx_direction = find(pa.fullFactorial(:,4) == pa.direction(conditions(cond,2)));
                    idx = intersect(idx_speed, idx_direction);
                    if isfield(pa, 'goodTrials')
                        idx = intersect(idx, pa.goodTrials{data_set});
                    end
                    data_session(cond,:) = [nan(1) nan(1) pa.speed(conditions(cond,1)), rad2deg(pa.direction(conditions(cond,2))), sum(eq(pa.LR(idx), pa.LRresponse(idx))), length(idx)];

                end
%                 pa.speed = [0.3];
%                 pa.direction = deg2rad(90);
                [dconst, dsurr] = DistanceToConstraint(ds, pa, depth_range, depth_est);
                data_session(:,1) = dconst(:);
                data_session(:,2) = dsurr(:);
                data = [data; data_session];
                data_const = [data(:,1) data(:,end-1:end)]; % to surround data_const = [data(:,1) data(:,end-1:end)];
                data_surr = [data(:,2) data(:,end-1:end)];
                count = count+1;


            else
                data_const = [data_const; pa.data_const];
                data_surr = [data_surr; pa.data_surr];
%                 options.poolxTol = 0.005;
%                 options.nblocks = 24;

            end
        end

    end
    disp([pa.subjectName, ' goodtrials = ', num2str(sum(data_const(:,end))/(length(data_const)*20)*100), '%']) %pa.nRepeats=20 except MP in one environment


    options             = struct;   % initialize as an empty struct
    options.sigmoidName = 'weibull'; %'weibull';
    options.expType     = '2AFC';   % choose 2-AFC as the paradigm of the experiment
    % this sets the guessing rate to .5 and
    % fits the rest of the parameters
    options.fixedPars = NaN(5,1);
    if ideal_eye
%         4 speeds, 6 directions
        options.dataColor = [255,153,255; 255,102,255; 255,51,255; 204,0,204;
            255,153,153; 255,102,102; 255,51,51; 204,0,0;
            255,204,153; 255,178,102; 255, 153, 51; 204,102,0;
            204,255,153; 178,255,102; 153,255,51; 102,204,0;
            153,255,255; 102,255,255; 51,255,255; 0,204,204;
            153,153,255; 102,102,255; 51,51,255; 0,0,204]/255;
% options.poolxTol = 0.0025;
    else
        options.poolxTol = 0.005;

    end

    %     options.fixedPars(3) = .01; % fix lapse rate at 1%
%         options.fixedPars(5) = 0;       % fix eta (dispersion) to zero

    result_const = psignifit(data_const,options);

    figure, plotPsych(result_const, options);
    title([subjects(s),' distance to constraint, depth range = ', num2str(depth_range), 'depth estimate = ', num2str(depth_est)])
    
    if isfield(pa, 'goodTrials')
        figname = [subjects(s)+'_const_'+stims(1)+'_multiplicative'+'.eps'];
    else
        figname = [subjects(s)+'_const_'+stims(1)+'.eps'];
    end
    saveas(gcf, fullfile(figFolder, figname), 'epsc')

%     options.fixedPars(3) = result_const.Fit(3); % fix lapse rate at calculated for constraint
    result_surr = psignifit(data_surr,options);


    figure, plotPsych(result_surr, options);
    title([subjects(s)+ ' distance to surround'])
    if isfield(pa, 'good_trials')
        figname = [subjects(s)+'_surr_'+stims(1)+'_multiplicative'+'.eps'];
    else
        figname = [subjects(s)+'_surr_'+stims(1)+'.eps'];
    end
    saveas(gcf, fullfile(figFolder, figname), 'epsc')


    display([subjects(s) + ' const dev = '+   num2str(result_const.deviance)])
    display([subjects(s) + ' surr dev = ' + num2str(result_surr.deviance)])

end

%% iterate over different values of distance to const
% additive distance
% distances  = logspace(-2, 1.5, 20);
% distances  = logspace(-1.5, -.5, 10);

% multiplicative
depth_ranges = linspace(1.0001, 1.3, 10); % 0.5% to 100%
% percentages = round((distances-1)*100,2);

% depth est
depth_estimates = linspace(-.4,.4, 11); %+.35m to hit the ground plane


% window of distances cube could physically be distances  = linspace(.05, .6, 10);
% constraint_length_opt distances = linspace(0.025, 0.15, 10)
dev = zeros(length(depth_estimates),length(depth_ranges));

tic
for depth = 1:length(depth_estimates)
    for range = 1:length(depth_ranges) 
        [dconst, dsurr] = DistanceToConstraint(ds, pa, depth_ranges(range), depth_estimates(depth));
        data_const(:,1) = repmat(dconst(:), count, 1);
        % run psignifit
        result = psignifit(data_const,options);
    %     figure, plotPsych(result, options);
    %     title(num2str(distances(d)))
    %     set(gca, 'FontSize', 16)
    %     thresh = exp(result.Fit(1));
        dev(depth,range) = result.deviance;
    %     if mod(d,5) == 0
            display([depth_estimates(depth),depth_ranges(range)])
    %     end
    end
end
toc

% imagesc(dev)

% [val,idx] = min(dev);
% display(['dev = ' num2str(val),' distance = ' num2str(depth_ranges(idx))])

xx = meshgrid(depth_estimates, depth_ranges)';
yy = meshgrid(depth_ranges, depth_estimates);
zz = dev;
figure, surf(xx,yy,zz)

[min_val,idx]=min(dev(:))
[row,col]=ind2sub(size(dev),idx);
depth_ranges(col)
depth_estimates(row)

title([subjects stims(1) 'opt depth range = ' num2str(depth_ranges(col)) 'opt depth estimate = ', num2str(depth_estimates(row)), 'deviance = ' num2str(min_val)])

% % 
% figure, semilogx(distances,dev,'o-','LineWidth', 5)
% set(gca, 'FontSize', 16)

% figure, hold on, plot(depth_ranges,dev,'o-','LineWidth', 5)
% set(gca, 'FontSize', 16)
% ticks = 1:5:length(depth_ranges);

% xticks(distances(ticks))
% pticks = split(num2str(percentages));
% xticklabels(pticks(ticks))
% 
% xticks([1.1, 1.2, 1.3])
% xticklabels({"1.1", "1.2", "1.3"})
% 
% hold on, plot(distances, dev, 'o-','LineWidth', 5)

% if idx == length(distances)
%     distances_ext = linspace(distances(end)+mean(diff(distances)), (distances(end))+10*mean(diff(distances)), 10);
%     dev_ext = zeros(1,length(distances_ext));
%     
%     for d = 1:length(distances_ext)
%         [dconst, dsurr] = DistanceToConstraint(ds, pa, distances_ext(d));
%         data_const(:,1) = repmat(dconst(:), count, 1);
%         
%         % run psignifit
%         result = psignifit(data_const,options);
%         figure, plotPsych(result, options);
%         %     thresh = exp(result.Fit(1));
%         dev_ext(d) = result.deviance;
%     end
%     
%     dev = [dev dev_ext];
%     distances = [distances distances_ext];
% elseif idx == 1
%     distances_prev = linspace(max(0,distances(1)-10*mean(diff(distances))), distances(1)-mean(diff(distances)), 10);
%     dev_prev = zeros(1,length(distances_prev));
%     
%     for d = 1:length(distances_prev)
%         [dconst, dsurr] = DistanceToConstraint(ds, pa, distances_prev(d));
%         data_const(:,1) = repmat(dconst(:), count, 1);
%         
%         % run psignifit
%         result = psignifit(data_const,options);
%         figure, plotPsych(result, options);
%         %     thresh = exp(result.Fit(1));
%         dev_prev(d) = result.deviance;
%     end
%     
%     dev = [dev_prev dev];
%     distances = [distances_prev distances];
% end

% subj distances deviances
% save dev and distances
% if contains(stims, "monocular")
%     save( fullfile(subjects(s), distances.mat) "distances"


%% deviance differences
x = categorical({'constraint', 'surround'});
conds = categorical({'full', 'dots', 'monocular'});
conds = reordercats(conds, {'full', 'dots', 'monocular'});
subjects = ["MP","DL","PL", "MG", "SM", "JL", "IK", "JO", "KZ","IG"];

const.full = [16, 37, 52, 24, 29, 33, 38, 61, 21, 40];
surr.full = [312,137,184,195,219,225, 112, 230,143,139];

const.dots = [71	93	29	83	62	91	85	26	120	38];
surr.dots = [233	269	106	144	134	215	203	301	264	110];

const.mono = [41	29	44	18	42	29	40	47	39	24];
surr.mono = [203	208	253	150	250	225	153	221	217	157];

% with eyetracking (no JL, IG)
subjects = ["MP","DL","PL", "MG", "SM", "IK", "JO", "KZ"]; %"IG"

const.full = [15.6351, 5.8561, 13.427, 20.5445, 12.5146, 17.7019, 4.0937, 5.3757, 5.861];
surr.full = [68.7766, 42.1598, 36.401, 28.4526, 58.2991, 42.5696, 45.2817, 22.1162, 23.7202];

const.dots = [55.4234	70.218	9.7616	52.2366	24.4333		53.953	26.908	22.4821 3.8244];
surr.dots = [125.7537	122.5027	37.2422	38.6214	39.2498		73.7502	93.9493	100.7684 16.9668];

const.mono = [4.7275	12.8345	5.1145	7.4692	4.8169		11.1978	10.7408	8.3045	10.7063];
surr.mono = [52.7765	65.9968	30.8971	60.341	35.1009		53.887	34.6784	67.3172	16.7235];

% redoing eyetracking to update constraint calculation (no added rotation
% before to constraint)
subjects = ["MP","DL","PL", "MG", "SM", "IK", "JO", "KZ"]; %"IG"

const.full = [22.3472	1.4375	53.8067	6.0718	3.9189		10.908	5.3659	0.84972	];
surr.full = [209.1479	83.329	71.928	73.1586	142.0699		26.9868	114.8144	77.3615	];

const.dots = [47.8804	38.1428	0.37073	30.3843	42.5853		12.5412	58.3317	94.6148	];
surr.dots = [97.8805	215.0529	42.0954	56.2199	73.4183		98.8093	153.8852	188.0491];

const.mono = [31.4526	12.8345	9.2975	25.7822	2.1776		0.69794	4.1503	17.516	];
surr.mono = [131.869	64.6475	142.8303	58.7146	127.65		16.9672	113.381	104.641	];


% with just good trials, no eyetracking, seeded cubes no JL or IG
subjects = ["MP","DL","PL", "MG", "SM", "IK", "JO", "KZ"];

const.full= [15.6395	28.7171	35.2079	24.8349	29.7911		37.9016	44.5579	19.4831];
surr.full = [235.6173	137.8694	103.4511	133.819	153.7909		50.0077	159.5812	72.6018];

const.dots = [63.897	86.9814	53.2302	67.3524	60.8678		69.4851	26.3628	71.8504];
surr.dots = [202.006	206.9385	207.1306	93.8245	216.0456		128.5636	237.5121	148.3069];

const.mono = [27.1673	29.6859	46.3152	24.0491	40.472		32.7438	41.0679	29.8143];
surr.mono = [120.5473	138.4028	231.3806	79.1458	224.4793		52.4259	151.5077	150.3862];


% multiplicative range, sedded cubes no JL or IG

const.full = [15.9151	30.2447	35.3319	24.5125	30.1324	38.3559	46.3583	18.7115];
surr.full = [235.6138	137.748	103.4843	133.8528	153.9385	50.0077	159.6656	72.7322];
const.dots = [67.0991	96.9531	57.3236	69.7838	66.7175	72.1655	28.2833	76.8247];
surr.dots = [223.7979	234.2082	216.3863	124.0553	229.7325	143.1377	249.8042	164.549];
const.mono = [27.6668	29.6249	49.3243	24.8407	44.1944	32.2815	44.8319	29.4757];
surr.mono = [120.5473	138.4324	231.4612	79.1445	224.579	52.4371	151.5863	150.4013];

% with eyetracking 0,75 thresh (no JL, IG)
subjects = ["MP","DL","PL", "MG", "SM", "IK", "JO", "KZ"]; %"IG"

const.full = [14.5646	33.0393	12.8474	31.9263	28.7374		4.0394	2.9203	22.5231];
surr.full = [132.9029	40.0209	22.2642	26.2903	40.658		19.3174	41.2065	26.0243];

const.dots = [34.5412	0.83839	0.53432	34.5412	20.8139		12.067	63.3388	2.3312];
surr.dots = [68.5071	61.2749	9.9836	68.5071	23.7643		39.7869	100.545	5.3749];

const.mono = [48.492	6.7729	67.751	1.2891	2.1174		1.3429	1.6376	17.9469];
surr.mono = [76.9227	42.8086	39.8211	8.0531	83.7478		12.7231	7.3378	35.1828];


% % const vs surround
% fig = figure();
% for d = 1:length(const)
%     hold on, plot(x, [const(d), surr(d)],'.-','MarkerSize',20,'LineWidth', 2)
% end
% title('Monocular')
% set(gca, 'FontSize', 16)
% 
% ylim([0. 350])

c = lines(length(subjects));

% const for each stimulus
x = categorical(["full" "monocular" "dots"]);
x = reordercats(x,string(x));

fig = figure();
for range = 1:length(const.full)
%     hold on, plot(conds, [const.full(d)', const.dots(d)', const.mono(d)'], '.-','MarkerSize',20,'LineWidth', 2)
%     hold on, plot([conds(1) conds(3)], [const.full(d)', const.mono(d)'], '.-','MarkerSize',20,'LineWidth', 2) 
    hold on, plot(x, [surr.full(range)', surr.mono(range)', surr.dots(range)'], '.-', 'color', c(range,:),'MarkerSize',20,'LineWidth', 2) 
end
set(gca, 'FontSize', 16)
hold on, plot(x, mean([surr.full', surr.mono', surr.dots']), '.-', 'color', [.25 .25 .25],'MarkerSize',30,'LineWidth', 5)
legend([subjects, "mean"])

%   figname = "deviance_comparison";
%  saveas(gcf, fullfile(figFolder, figname), 'epsc')

fig = figure();
for range = 1:length(const.full)
%     hold on, plot(conds, [const.full(d)', const.dots(d)', const.mono(d)'], '.-','MarkerSize',20,'LineWidth', 2)
%     hold on, plot([conds(1) conds(3)], [const.full(d)', const.mono(d)'], '.-','MarkerSize',20,'LineWidth', 2) 
    hold on, plot(x, [const.full(range)', const.mono(range)', const.dots(range)'], '.-', 'color', c(range,:),'MarkerSize',20,'LineWidth', 2) 

end
set(gca, 'FontSize', 16)
hold on, plot(x, mean([const.full', const.mono', const.dots']), '.-', 'color', [.25 .25 .25],'MarkerSize',30,'LineWidth', 5)
legend([subjects, "mean"])
ylim([0 100])
figname = "deviance_combined";
% saveas(gcf, fullfile(figFolder, figname), 'epsc')

% scatter const vs surround
figure
hold on, scatter(surr.full, const.full, 50, 'filled')
hold on, scatter(surr.dots, const.dots, 50, 'filled')
hold on, scatter(surr.mono, const.mono, 50, 'filled')
axis equal
xlim([0 350])
ylim([0 150])
xlabel('surround')

hold on
plot([min([xlim ylim]) max([xlim ylim])], [min([xlim ylim]) max([xlim ylim])], '--k')


%% plot optimal depth ranges
% const for each stimulus
x = categorical(["full" "monocular" "dots"]);
x = reordercats(x,string(x));

subjects = ["MP","DL","PL", "MG", "SM","IK", "JO", "KZ"];
c = lines(length(subjects));

%additive
full = [0.05	0.1274	0.05	0.05	0.127		0.19	0.015];
dots = [0.3	    0.2976	0.29	NaN	    0.45		0.12	0.45];
monocular = [0.29	0.2976	0.12	0.19	0.127		0.12	0.05];

% full = (3*[1.0311	1.0828	1.0311	1.0001	1.0725	1.0518	1.0932	1.0104]-3)*100;
% dots = (3*[1.2069	1.2276	1.2586	1.1035	1.1966	1.1759	1.2173	1.1966]-3)*100;
% monocular = (3*[1.2483	1.0828	1.0725	1.0621	1.0828	1.0001	1.0828	1.0311]-3)*100;

% 
full = [1.0311	1.0828	1.0311	1.0001	1.0725	1.0518	1.0932	1.0104];
dots = [1.2069	1.2276	1.2586	1.1035	1.1966	1.1759	1.2173	1.1966];
monocular = [1.2483	1.0828	1.0725	1.0621	1.0828	1.0001	1.0828	1.0311];

% depth est
full = [0	0	0.04	0	0.04	0.08	-0.04	0.04];
dots = [0.16	0.16	-0.04	0.2	0	0.16	0.04	0.16];
monocular = [0.04	0.04	-0.04	0.04	-0.04	0.12	-0.04	0.08];

% depth est and range
full = [0	-0.08	0	0	0.16	0.08	0	0;
1.0001	1.1334	1.0334	1.0001	1.1667	1.0667	1.1001	1.0001];

dots = [0.24	0.24	0	0.4	0.08	0.32	0.16	0.32;
1.1667	1.1667	1.2334	1.3	1.2667	1.2334	1.2	1.2];

monocular = [-0.08	0.08	-0.16	0	-0.08	0.08	-0.08	0.08;
1.2	1.1001	1.1334	1.0667	1.1001	1.0001	1.1001	1.0667];




fig = figure();
for range = 1:length(full)
%     hold on, plot(conds, [const.full(d)', const.dots(d)', const.mono(d)'], '.-','MarkerSize',20,'LineWidth', 2)
%     hold on, plot([conds(1) conds(3)], [const.full(d)', const.mono(d)'], '.-','MarkerSize',20,'LineWidth', 2) 
    hold on, plot(x, [full(1,range) monocular(1,range) dots(1,range)], '.-', 'color', c(range,:),'MarkerSize',20,'LineWidth', 2) 
end
set(gca, 'FontSize', 16)
hold on, plot(x, nanmean([full(1,:)' monocular(1,:)' dots(1,:)']), '.-', 'color', [.25 .25 .25],'MarkerSize',30,'LineWidth', 5)
legend([subjects, "mean"])

y = nanmean([full(1,:)' monocular(1,:)' dots(1,:)']);
err = std([full(1,:)' monocular(1,:)' dots(1,:)'])/sqrt(length(full));
errorbar(x,nanmean([full(1,:)' monocular(1,:)' dots(1,:)']), err, 'LineWidth', 2)


fig = figure();
for range = 1:length(full)
%     hold on, plot(conds, [const.full(d)', const.dots(d)', const.mono(d)'], '.-','MarkerSize',20,'LineWidth', 2)
%     hold on, plot([conds(1) conds(3)], [const.full(d)', const.mono(d)'], '.-','MarkerSize',20,'LineWidth', 2) 
    hold on, plot(x, [full(2,range) monocular(2,range) dots(2,range)], '.-', 'color', c(range,:),'MarkerSize',20,'LineWidth', 2) 
end
set(gca, 'FontSize', 16)
hold on, plot(x, nanmean([full(2,:)' monocular(2,:)' dots(2,:)']), '.-', 'color', [.25 .25 .25],'MarkerSize',30,'LineWidth', 5)
legend([subjects, "mean"])

y = nanmean([full(2,:)' monocular(2,:)' dots(2,:)']);
err = std([full(2,:)' monocular(2,:)' dots(2,:)'])/sqrt(length(full));
errorbar(x,nanmean([full(2,:)' monocular(2,:)' dots(2,:)']), err, 'LineWidth', 2)




figname = "optimal_depthrange";

%errorbar plot
bias = [full(1,:)' monocular(1,:)' dots(1,:)'];
variance = [full(2,:)' monocular(2,:)' dots(2,:)'];
est  = pa.objectdist+bias;
ypos = est.*variance-est;
yneg = est./variance-est;

figure, errorbar(x, est(6,:), yneg(6,:), ypos(6,:))

meanest = mean(est);
meanvar = mean(variance);
ymeanpos = meanest.*meanvar-pa.objectdist;
ymeanneg = meanest./meanvar-pa.objectdist;
figure
hold on, errorbar(x, meanest, ymeanneg, ymeanpos)

%%
%load constraint_legnth_opt.mat
% replace correct array in mat file

save (fullfile(analysisFolder, "constraint_length_opt.mat"), "constraint_length_opt")

%% full vs dots vs monocular
figure, hold on, plot(depth_ranges, dev, 'LineWidth', 2)

figure
for s = 1:3
    hold on, plot(constraint_length_opt(8).distances{s}, constraint_length_opt(8).dev{s}, 'LineWidth', 2)
    
end
legend('full', 'dots', 'monocular')
set(gca, 'FontSize', 16)

%% eye movements improve resposne?


for s  = 1:length(subjects)

    data = [];
    %load appropriate files
    count = 0;
    for f = 1:length(S)
        subj = contains(S(f).name,subjects(s));
        sti = contains(S(f).name,stims) && contains(S(f).name, 'eyetracking');
        %add "eyetrakcing" to have pa.good_trials, and/or dconst and dsurround based on vertical eye movements

        if subj && sti
            load(fullfile(dataFolder,S(f).name));
            display(fullfile(dataFolder,S(f).name));

            if ideal_eye
                n_conditions = length(pa.speed)*length(pa.direction);
                conditions = fullfact([numel(pa.speed), numel(pa.direction)]);
                data_session = [];
                for cond = 1:n_conditions
                    idx_speed = find(pa.fullFactorial(:,3) == pa.speed(conditions(cond,1)));
                    idx_direction = find(pa.fullFactorial(:,4) == pa.direction(conditions(cond,2)));
                    idx = intersect(idx_speed, idx_direction);
                    if isfield(pa, 'good_trials')
                        idx = intersect(idx, pa.goodTrials{data_set});
                    end
                    data_session(cond,:) = [nan(1) nan(1) pa.speed(conditions(cond,1)), rad2deg(pa.direction(conditions(cond,2))), sum(eq(pa.LR(idx), pa.LRresponse(idx))), length(idx)];

                end
                [dconst, dsurr] = DistanceToConstraint(ds, pa, depth_range);
                data_session(:,1) = dconst(:);
                data_session(:,2) = dsurr(:);
                data = [data; data_session];
                data_const = [data(:,1) data(:,end-1:end)]; % to surround data_const = [data(:,1) data(:,end-1:end)];
                data_surr = [data(:,2) data(:,end-1:end)];
                count = count+1;


            else
                data_const = [data_const; pa.data_const];
                data_surr = [data_surr; pa.data_surr];
%                 options.poolxTol = 0.005;
%                 options.nblocks = 24;

            end
        end

    end
end
