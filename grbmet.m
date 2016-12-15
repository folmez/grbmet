function [escape_times, pos_vec, time_vec] = grbmet(varargin)
% [escape_times, pos_vec, time_vec] = GRBMET(varargin) simulates regulated
% (or reflected) Brownian Motion and generates escape times at a given
% point.
%
% The parameter M in this code represents the location of a boundary 
% control. Such a Brownian Motion is called "regulated" or "reflected".
% When M=inf, the model is the standard Brownian Motion. When M=finite,
% Browninan motion will be absorbed at -M when it attempts to go below it.
%
% e.g.
%   grbmet('M', 5,   'nr_bouts', 1000);
%   grbmet('M', inf, 'nr_bouts', 1000);   
%
% If you find our method useful, please cite the following paper: 
% [1] RBM
%
%   1. If you want to change the escape distance, try:
%       grbmet('a', 3);
%
%   2. If you want to record position and time vectors of each Brownian
%   Motion escape bout, try:
%       grbmet('record_position_and_time', 1);
%
%   3. If you want to plot trajectories of each Brownian Motion bout, try:
%       grbmet('plot_trajectories', 1);
%
% Version 1.0   (2016 December)
% Version 0.0	(2016 July)
% Copyright (C) Fatih Olmez (folmez@gmail.com)
% grbmet comes with ABSOLUTELY NO WARRANTY
%

% Input parameters
M = 10;
a = 1;
nr_bouts = 1e3;
display_results = 1;
record_position_and_time = 0;
plot_trajectories = 0;
dt = 0.0001;

i=1;
while i<=length(varargin),
    switch varargin{i},
        case 'M',                   M = varargin{i+1};
        case 'a',                   a = varargin{i+1};
        case 'nr_bouts',           	nr_bouts = varargin{i+1};
        case 'display_results',     display_results = varargin{i+1};
        case 'record_position_and_time'
            record_position_and_time = varargin{i+1};
        case 'plot_trajectories'
            plot_trajectories = varargin{i+1};
            record_position_and_time = 1;
        otherwise,
            display(varargin{i});
            error('Unexpected inputs!!!');
    end
    i = i+2;
end

% Model
sqrt_dt = sqrt(dt); % simulation time step
escape_times = zeros(nr_bouts,1);
initial_pt = 0;
nr_jumps_guess = 1e6;
if display_results
    fprintf('Sim# \tEscape time\tStart\tEnd\t#jumps\tSimTime\n');
end
pos_vec{nr_bouts} = [];
time_vec{nr_bouts} = [];
tSIM = tic;
for j = 1:nr_bouts
    % Initialize
    t = zeros(nr_jumps_guess, 1);   % time vector
    X = zeros(nr_jumps_guess, 1);   % position vector
    total_time = 0;            
    X(1) = initial_pt;      
    
    % Run Brownian Motion
    current_stop = initial_pt;
    ind_ct = 2;
    while true
        % Jump to next point
        step = sqrt_dt*randn;
        next_stop = current_stop + step;
        
        % Control at -M
        if next_stop < -M, next_stop = -M; end
        
        % Update
        current_stop = next_stop; total_time = total_time + dt;
        
        % Record state and time if desired
        if record_position_and_time
            if ind_ct>length(X)
                X(end+1:end+nr_jumps_guess) = zeros(nr_jumps_guess,1 );
                t(end+1:end+nr_jumps_guess) = zeros(nr_jumps_guess,1 );
            end
            X(ind_ct) = current_stop;
            t(ind_ct) = t(ind_ct-1) + dt;
        end
        
        % Jump count
        ind_ct = ind_ct + 1;
        
        % Complete the process if the escape goal is achieved
        if current_stop>a, last_pt = current_stop; break; end        
    end
    
    % Print bout summary
    if display_results
        fprintf('%i\t%6.2f\t\t%1.2f\t%1.2f\t%1.1e\t%3.2fm\n', ...
            j, total_time, initial_pt, last_pt, ind_ct-1, toc(tSIM)/60);
    end
    
    % Remove unused entries from the position and time vectors
    t(ind_ct:end) = []; X(ind_ct:end) = [];
        
    % Plot trajectories
    if record_position_and_time && plot_trajectories
        figure, plot(t, X); hold on;
        plot([min(t) max(t)], [-M -M], 'k');
        plot([min(t) max(t)], [a a], 'g');
        ylim([-M-0.5 a+0.5]);
        title('Press a key to close this figure and continue ...');
        legend('Brownian motion', 'Reflecting boundary', 'Target', ...
            'Location', 'SouthWest');
        pause; close;
    end
    
    % Record position, time vectors and the bout duration        
    pos_vec{j} = X; time_vec{j} = t; escape_times(j) = total_time;
end

% Plot approximate PDF of bout durations
nr_bins = min(round(nr_bouts*0.1), 5e1);
etbins = logspace(log10(min(escape_times)), ...
    log10(max(escape_times)+1e-10), nr_bins);
[etn, etx] = hist(escape_times, etbins);
bin_lengths = [etx(2)-etx(1), etx(3:end)-etx(1:end-2), ...
    etx(end)-etx(end-1)]*0.5;
etn = (etn./bin_lengths)/nr_bouts;
nis = find(etn);
figure, loglog(etx(nis), etn(nis), 'b.' , 'MarkerSize', 10);

% Calculate theoretical PDF: The formula used below was derived from the
% CDF formula from the following paper: Dybiec, Bartlomiej, Ewa
% Gudowska-Nowak, and Peter Hanggi. "Levy-Brownian motion on finite
% intervals: Mean first passage time analysis." Physical Review E 73.4
% (2006): 046104.
etPDFx = logspace(log10(min(escape_times)), log10(max(escape_times)), 100);
etPDFn = zeros(100, 1);
if isinf(M)
    etPDFn = a/sqrt(2*pi)*etPDFx.^(-1.5).*exp((-a^2)./(2*etPDFx));
else
    x0 = a;
    L = M+a;
    sigma = 1/sqrt(2);
    jj = 2*(0:1e4)+1;
    for k=1:length(etPDFn)
        etPDFn(k) = (4/pi) * sum( (1./jj) .* sin(jj*pi/2) .* ...
            cos(jj*pi/2*(L-x0)/L) .* ...
            (exp(-etPDFx(k)*(jj*pi/2*sigma/L).^2)) .* ...
            ((jj*pi/2*sigma/L).^2) );
    end
end

% Plot theoretical PDF
hold on; plot(etPDFx, etPDFn, 'k');
axis tight;
legend('Escape time approx, PDF', 'Escape time theoretical PDF');

end