function von_mises_slider_demo
    % ==============
    % 1. initialize parameters
    % ==============
    % range of theta
    theta = linspace(-180, 180, 361);

    % initialize parameters
    A_init    = 1;    % amplitude
    k_init    = 1;    % tuning width
    phi_init  = 0;    % theta pref

    % calculate initial curve
    y_init = von_mises_variant(theta, A_init, k_init, phi_init);

    % ==============
    % 2. plotting
    % ==============
    fig = figure('Name','Von Mises Variant Slider Demo',...
                 'Position',[400, 300, 700, 450]);

    % plot initial curve
    hPlot = plot(theta, y_init, 'LineWidth', 2);
    hold on; grid on;
    xlim([-180, 180]);
    ylim([0, 10]);

    % show current value of parameters
    titleStr = sprintf('A = %.2f, k = %.2f, \\phi = %.2f', A_init, k_init, phi_init);
    title(titleStr, 'FontSize', 12);

    % ==========================
    % 3. Adding slide bar (uicontrol)
    % ==========================
    % 3.1 Amplitude slide
    uicontrol('Style','text', ...
              'Position',[70 60 100 20], ...
              'String','Amplitude', ...
              'BackgroundColor',[0.8 0.8 0.8]);
    sliderA = uicontrol('Style','slider', ...
                        'Min', 0, 'Max', 10, ...
                        'Value', A_init, ...
                        'Position',[70 40 100 20], ...
                        'Callback',@(src,~)updatePlot());

    % 3.2 tuning width slide
    uicontrol('Style','text', ...
              'Position',[250 60 100 20], ...
              'String','tuning width', ...
              'BackgroundColor',[0.8 0.8 0.8]);
    sliderK = uicontrol('Style','slider', ...
                        'Min', 0, 'Max', 10, ...
                        'Value', k_init, ...
                        'Position',[250 40 100 20], ...
                        'Callback',@(src,~)updatePlot());

    % 3.3 preferred angle slide
    uicontrol('Style','text', ...
              'Position',[430 60 100 20], ...
              'String','preferred angle', ...
              'BackgroundColor',[0.8 0.8 0.8]);
    sliderPhi = uicontrol('Style','slider', ...
                          'Min', -180, 'Max', 180, ...
                          'Value', phi_init, ...
                          'Position',[430 40 100 20], ...
                          'Callback',@(src,~)updatePlot());

    % ================
    % 4. call back function
    % ================
    function updatePlot()
        % get current parameters
        A_val   = sliderA.Value;
        k_val   = sliderK.Value;
        phi_val = sliderPhi.Value;

        % calculate new distribution
        y_new = von_mises_variant(theta, A_val, k_val, phi_val);

        % update curve
        set(hPlot, 'YData', y_new);

        % refresh title
        titleStr = sprintf('A = %.2f, k = %.2f, \\phi = %.2f', A_val, k_val, phi_val);
        title(titleStr, 'FontSize', 12);

        drawnow;
    end

end

% ========================================
%  Function: von_mises_variant
%  M(θ) = A * exp{k [cos(2(θ - φ)) - 1]}
% ========================================
function y = von_mises_variant(theta, A, k, phi)
    y = A * exp( k .* ( cosd( 2*(theta - phi) ) - 1 ) );
end
