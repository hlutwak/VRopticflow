%% Analysis code
% get behavior files, conccatenate same trials, run psignifit

% add psignifit toolbox
addpath('/Users/hopelutwak/Documents/MATLAB/psignifit')
addpath(genpath('/Applications/Psychtoolbox'))
addpath('/Users/hopelutwak/Documents/GitHub/VRopticflow/Analysis')
% addpath(genpath('C:\Users\hlutw\Documents\MATLAB\psignifit-master'))
% assign data folder
dataFolder = '/Users/hopelutwak/Documents/GitHub/VRopticflow/Data';
% dataFolder = 'cC:\Users\hlutw\Documents\GitHub\VRopticflow\Data';
figFolder = '/Users/hopelutwak/Documents/GitHub/VRopticflow/Figures';
analysisFolder = '/Users/hopelutwak/Documents/GitHub/VRopticflow/Analysis';
% names of files
S = dir(fullfile(dataFolder,'*.mat'));

% which subjects data to analyze
subjects = ["DL"]; %,"DL", "PL","MG", "SM", "IK", "JO", "KZ","IG"

% all: "PL", "MP", "SM", "JL", "IK", "JO", "KZ", "IG"
% all with good eyetracking trials: subjects = ["MP","DL","PL", "MG", "SM", "IK", "JO", "KZ","IG"];

stims = ["full-1", "full-2"]; %add "copy" to have pa.good_trials, and/or dconst and dsurround based on vertical eye movements
% % ["full-1", "full-2"];["dots-1", "dots-2"] ["monocular-1", "monocular-2"]
ideal_eye = 1; % use measurements of data_const and data_surr based on ideal eye movements, otherwise use eyetracking vertical movements
depth_range = .05; % additive
depth_range = 1.01; % multiplicative
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
                    if isfield(pa, 'good_trials')
                        idx = intersect(idx, pa.good_trials);
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
    disp([pa.subjectName, ' goodtrials = ', num2str(sum(data_const(:,end))/(length(data_const)*20)*100), '%']) %pa.nRepeats=20 except MP in one environment


    options             = struct;   % initialize as an empty struct
    options.sigmoidName = 'weibull';
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
    title([subjects(s),' distance to constraint, depth range = ', num2str(depth_range)])
    
    if isfield(pa, 'good_trials')
        figname = [subjects(s)+'_const_'+stims(1)+'_goodtrials'+'.eps'];
    else
        figname = [subjects(s)+'_const_'+stims(1)+'.eps'];
    end
%     saveas(gcf, fullfile(figFolder, figname), 'epsc')

    options.fixedPars(3) = result_const.Fit(3); % fix lapse rate at calculated for constraint
    result_surr = psignifit(data_surr,options);


    figure, plotPsych(result_surr, options);
    title([subjects(s)+ ' distance to surround'])
    if isfield(pa, 'good_trials')
        figname = [subjects(s)+'_surr_'+stims(1)+'_goodtrials'+'.eps'];
    else
        figname = [subjects(s)+'_surr_'+stims(1)+'.eps'];
    end
%     saveas(gcf, fullfile(figFolder, figname), 'epsc')


    display([subjects(s) + ' const dev = '+   num2str(result_const.deviance)])
    display([subjects(s) + ' surr dev = ' + num2str(result_surr.deviance)])

end

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

% with eyetracking (no JL)
subjects = ["MP","DL","PL", "MG", "SM", "IK", "JO", "KZ","IG"];

const.full = [15.6351, 5.8561, 13.427, 20.5445, 12.5146, 17.7019, 4.0937, 5.3757, 5.861];
surr.full = [68.7766, 42.1598, 36.401, 28.4526, 58.2991, 42.5696, 45.2817, 22.1162, 23.7202];

const.dots = [55.4234	70.218	9.7616	52.2366	24.4333		53.953	26.908	22.4821 3.8244];
surr.dots = [125.7537	122.5027	37.2422	38.6214	39.2498		73.7502	93.9493	100.7684 16.9668];

const.mono = [4.7275	12.8345	5.1145	7.4692	4.8169		11.1978	10.7408	8.3045	10.7063];
surr.mono = [52.7765	65.9968	30.8971	60.341	35.1009		53.887	34.6784	67.3172	16.7235];

% with just good trials, seeded cubes no JL or IG
subjects = ["MP","DL","PL", "MG", "SM", "IK", "JO", "KZ"];

const.full= [15.6395	28.7171	35.2079	24.8349	29.7911		37.9016	44.5579	19.4831];
surr.full = [235.6173	137.8694	103.4511	133.819	153.7909		50.0077	159.5812	72.6018];

const.dots = [63.897	86.9814	53.2302	67.3524	60.8678		69.4851	26.3628	71.8504];
surr.dots = [202.006	206.9385	207.1306	93.8245	216.0456		128.5636	237.5121	148.3069];

const.mono = [27.1673	29.6859	46.3152	24.0491	40.472		32.7438	41.0679	29.8143];
surr.mono = [120.5473	138.4028	231.3806	79.1458	224.4793		52.4259	151.5077	150.3862];

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
for d = 1:length(const.full)
%     hold on, plot(conds, [const.full(d)', const.dots(d)', const.mono(d)'], '.-','MarkerSize',20,'LineWidth', 2)
%     hold on, plot([conds(1) conds(3)], [const.full(d)', const.mono(d)'], '.-','MarkerSize',20,'LineWidth', 2) 
    hold on, plot(x, [surr.full(d)', surr.mono(d)', surr.dots(d)'], '.-', 'color', c(d,:),'MarkerSize',20,'LineWidth', 2) 
end
set(gca, 'FontSize', 16)
hold on, plot(x, mean([surr.full', surr.mono', surr.dots']), '.-', 'color', [.25 .25 .25],'MarkerSize',30,'LineWidth', 5)
legend([subjects, "mean"])

%   figname = "deviance_comparison";
%  saveas(gcf, fullfile(figFolder, figname), 'epsc')

fig = figure();
for d = 1:length(const.full)
%     hold on, plot(conds, [const.full(d)', const.dots(d)', const.mono(d)'], '.-','MarkerSize',20,'LineWidth', 2)
%     hold on, plot([conds(1) conds(3)], [const.full(d)', const.mono(d)'], '.-','MarkerSize',20,'LineWidth', 2) 
    hold on, plot(x, [const.full(d)', const.mono(d)', const.dots(d)'], '.-', 'color', c(d,:),'MarkerSize',20,'LineWidth', 2) 

end
set(gca, 'FontSize', 16)
hold on, plot(x, mean([const.full', const.mono', const.dots']), '.-', 'color', [.25 .25 .25],'MarkerSize',30,'LineWidth', 5)
legend([subjects, "mean"])
ylim([0 100])
figname = "deviance_combined";
% saveas(gcf, fullfile(figFolder, figname), 'epsc')

% scatter const vs surround
figure, scatter(const.dots, surr.dots, 50, 'filled')
axis equal
xlim([0 100])
ylim([0 350])

hold on
plot([min([xlim ylim]) max([xlim ylim])], [min([xlim ylim]) max([xlim ylim])], '--k')


%% plot optimal depth ranges
% const for each stimulus
x = categorical(["full" "monocular" "dots"]);
x = reordercats(x,string(x));

subjects = ["MP","DL","PL", "MG", "SM", "JO", "KZ"];
c = lines(length(subjects));


full = [0.05	0.1274	0.05	0.05	0.127		0.19	0.015];
dots = [0.3	    0.2976	0.29	NaN	    0.45		0.12	0.45];
monocular = [0.29	0.2976	0.12	0.19	0.127		0.12	0.05];

fig = figure();
for d = 1:length(full)
%     hold on, plot(conds, [const.full(d)', const.dots(d)', const.mono(d)'], '.-','MarkerSize',20,'LineWidth', 2)
%     hold on, plot([conds(1) conds(3)], [const.full(d)', const.mono(d)'], '.-','MarkerSize',20,'LineWidth', 2) 
    hold on, plot(x, [full(d) monocular(d) dots(d)], '.-', 'color', c(d,:),'MarkerSize',20,'LineWidth', 2) 
end
set(gca, 'FontSize', 16)
hold on, plot(x, nanmean([full' monocular' dots']), '.-', 'color', [.25 .25 .25],'MarkerSize',30,'LineWidth', 5)
legend([subjects, "mean"])


%% iterate over different values of distance to const
% additive distance
% distances  = logspace(-2, 1.5, 20);
% distances  = logspace(-1.5, -.5, 10);

% multiplicative
distances = linspace(1.0001, 1.3, 30); % 0.5% to 100%
percentages = round((distances-1)*100,2);


% window of distances cube could physically be distances  = linspace(.05, .6, 10);
% constraint_length_opt distances = linspace(0.025, 0.15, 10)
dev = zeros(1,length(distances));

tic
for d = 1 :length(distances) 
    [dconst, dsurr] = DistanceToConstraint(ds, pa, distances(d));
    data_const(:,1) = repmat(dconst(:), count, 1);
    % run psignifit
    result = psignifit(data_const,options);
    figure, plotPsych(result, options);
    title(num2str(distances(d)))
    set(gca, 'FontSize', 16)
%     thresh = exp(result.Fit(1));
    dev(d) = result.deviance;
%     if mod(d,5) == 0
        display(d)
%     end
end

[val,idx] = min(dev);
toc

% 
figure, semilogx(distances,dev,'o-','LineWidth', 5)
set(gca, 'FontSize', 16)

figure, plot(distances,dev,'o-','LineWidth', 5)
set(gca, 'FontSize', 16)
ticks = 1:5:length(distances);
xticks(distances(ticks))
pticks = split(num2str(percentages));
xticklabels(pticks(ticks))


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
%%
%load constraint_legnth_opt.mat
% replace correct array in mat file

save (fullfile(analysisFolder, "constraint_length_opt.mat"), "constraint_length_opt")

%% full vs dots vs monocular
figure, hold on, plot(distances, dev, 'LineWidth', 2)

figure
for s = 1:3
    hold on, plot(constraint_length_opt(8).distances{s}, constraint_length_opt(8).dev{s}, 'LineWidth', 2)
    
end
legend('full', 'dots', 'monocular')
set(gca, 'FontSize', 16)
