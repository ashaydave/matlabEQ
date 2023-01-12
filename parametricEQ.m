classdef parametricEQ < audioPlugin
    properties
        HighShelf_freq = 10000;
        HighShelf_gain = 0;
        HighMid_freq = 5000;
        HighMid_gain = 0;
        LowMid_freq = 500;
        LowMid_gain = 0;
        LowShelf_freq = 100;
        LowShelf_gain = 0;

        HPF_freq = 50;

        fs = 44100;
        fn = 22050;

        BYPASS = 'off';
       
    end

    %% Plugin UI
    properties (Constant)
        PluginInterface = audioPluginInterface(...
            audioPluginParameter('HighShelf_freq',...
            'DisplayName','High Shelf Freq',...
            'Label', 'Hz',...
            'Mapping',{'log', 2500, 20000}), ...
            audioPluginParameter('HighShelf_gain',...
            'DisplayName','High Shelf Gain',...
            'Label', 'dB',...
            'Mapping',{'lin', -12, 12}), ...
            audioPluginParameter('HighMid_freq',...
            'DisplayName','High Mid Freq',...
            'Label', 'Hz',...
            'Mapping',{'log', 800, 12500}),...
            audioPluginParameter('HighMid_gain',...
            'DisplayName','High Mid Gain',...
            'Label', 'dB',...
            'Mapping',{'lin', -12, 12}),...
            audioPluginParameter('LowMid_freq',...
            'DisplayName','Low Mid Freq',...
            'Label', 'Hz',...
            'Mapping',{'log', 255, 500}), ...
             audioPluginParameter('LowMid_gain',...
            'DisplayName','Low Mid Gain',...
            'Label', 'dB',...
            'Mapping',{'lin', -12, 12}), ...
             audioPluginParameter('LowShelf_freq',...
            'DisplayName','Low Shelf Freq',...
            'Label', 'Hz',...
            'Mapping',{'log', 75, 225}), ...
             audioPluginParameter('LowShelf_gain',...
            'DisplayName','Low Shelf Gain',...
            'Label', 'dB',...
            'Mapping',{'lin', -12, 12}), ...
             audioPluginParameter('HPF_freq',...
            'DisplayName','DC Block Freq',...
            'Label', 'Hz',...
            'Mapping',{'log', 30, 400}), ...
             audioPluginParameter('BYPASS',...
            'DisplayName','Bypass',...
            'Label', 'Hz',...
            'Mapping',{'enum', 'off', 'on'}));
    end
    %% Bi Quads data structs
    properties (Access = private)
        filter_HS = struct('w', [0 0; 0 0], 'a0', 1, 'a1', 0, 'a2', 0, 'b0', 1, 'b1', 0, 'b2', 0); %canonical form which uses the fewest number of delays and two zeroes for left and right channel
        filter_HMF = struct ('w', [0 0; 0 0], 'a0', 1, 'a1', 0, 'a2', 0, 'b0', 1, 'b1', 0, 'b2', 0);
        filter_LMF = struct ('w', [0 0; 0 0], 'a0', 1, 'a1', 0, 'a2', 0, 'b0', 1, 'b1', 0, 'b2', 0);
        filter_LS = struct ('w', [0 0; 0 0], 'a0', 1, 'a1', 0, 'a2', 0, 'b0', 1, 'b1', 0, 'b2', 0);
        filter_HPF = struct ('w', [0 0; 0 0], 'a0', 1, 'a1', 0, 'a2', 0, 'b0', 1, 'b1', 0, 'b2', 0);

    end
    %% functions

    methods
        function out = process(plugin,in) %flow process
            out = zeros(size(in)); %out vector but filled with zeros

            for ch = 1:min(size(in)) %1 dimension of input will be (1) two because stereo and 1 dimension of input will be (2) 512 for eg, buffer size/cache size of the daw
                x = in(:,ch);

                [y1, plugin.filter_HS.w(:,ch)] = processBiquad(x, plugin.filter_HS, ch);
                [y2, plugin.filter_HMF.w(:,ch)] = processBiquad(y1, plugin.filter_HMF, ch);
                [y3, plugin.filter_LMF.w(:,ch)] = processBiquad(y2, plugin.filter_LMF, ch);
                [y4, plugin.filter_LS.w(:,ch)] = processBiquad(y3, plugin.filter_LS, ch);
                [y5, plugin.filter_HPF.w(:,ch)] = processBiquad(y4, plugin.filter_HPF, ch);

                if strcmp(plugin.BYPASS, 'on') %if Bypass is on, then output the input as it is (no filtering)
                    out(:, ch) = x;
                else
                    out(:,ch) = y5;
                 
                end
            end
        end
%% Reset the plugin
        function reset(plugin)
            plugin.fs = getSampleRate(plugin); %check if sample rate has changed
            plugin.fn = plugin.fs/2; %set the Nyquist frequency
           
            %set filter states to 0
            plugin.filter_HS.w = [0 0; 0 0];
            plugin.filter_HMF.w = [0 0; 0 0];
            plugin.filter_LMF.w = [0 0; 0 0];
            plugin.filter_LS.w = [0 0; 0 0];
            plugin.filter_HPF.w = [0 0; 0 0];
        end
       
 %% High Shelf

        function set.HighShelf_freq(plugin, val) 
            plugin.HighShelf_freq = val; %update value of high shelf freq when changed
            update_HighShelf(plugin);
        end
        
        function set.HighShelf_gain(plugin, val)
            plugin.HighShelf_gain = val;
            update_HighShelf(plugin);
        end
        
        function update_HighShelf(plugin) %Using BiQuad equations (Robert Bristow Johnson)
            Q=0.5;
            f0=plugin.HighShelf_freq;
            gain = plugin.HighShelf_gain;
            w0=2*pi*f0/plugin.fs;
            alpha=sin(w0)/(2*Q);
            A=sqrt(db2mag(gain));
            
            plugin.filter_HS.a0 =    A*( (A+1) + (A-1)*cos(w0) + 2*sqrt(A)*alpha );
            plugin.filter_HS.a1 = -2*A*( (A-1) + (A+1)*cos(w0)                   );
            plugin.filter_HS.a2 =    A*( (A+1) + (A-1)*cos(w0) - 2*sqrt(A)*alpha );
            plugin.filter_HS.b0 =        (A+1) - (A-1)*cos(w0) + 2*sqrt(A)*alpha;
            plugin.filter_HS.b1 =    2*( (A-1) - (A+1)*cos(w0)                   );
            plugin.filter_HS.b2 =        (A+1) - (A-1)*cos(w0) - 2*sqrt(A)*alpha;
            
        end

%% High Mid
        function set.HighMid_freq(plugin, val)
            plugin.HighMid_freq = val;
            update_HighMid(plugin);
        end
        
        function set.HighMid_gain(plugin, val)
            plugin.HighMid_gain = val;
            update_HighMid(plugin);
        end
        
        function update_HighMid(plugin)
            Q=0.5;
            f0=plugin.HighMid_freq;
            gain = plugin.HighMid_gain;
            w0=2*pi*f0/plugin.fs;
            alpha=sin(w0)/(2*Q);
            A=sqrt(db2mag(gain));
            
            plugin.filter_HMF.a0 =   1 + alpha*A;
            plugin.filter_HMF.a1 =  -2*cos(w0);
            plugin.filter_HMF.a2 =   1 - alpha*A;
            plugin.filter_HMF.b0 =   1 + alpha/A;
            plugin.filter_HMF.b1 =  -2*cos(w0);
            plugin.filter_HMF.b2 =   1 - alpha/A;
        end

%% Low Mid
        function set.LowMid_freq(plugin, val)
            plugin.LowMid_freq = val;
            update_LowMid(plugin);
        end
        
        function set.LowMid_gain(plugin, val)
            plugin.LowMid_gain = val;
            update_LowMid(plugin);
        end
        
        function update_LowMid(plugin)
            Q=0.5;
            f0=plugin.LowMid_freq;
            gain = plugin.LowMid_gain;
            w0=2*pi*f0/plugin.fs;
            alpha=sin(w0)/(2*Q);
            A=sqrt(db2mag(gain));
            
            plugin.filter_LMF.a0 =   1 + alpha*A;
            plugin.filter_LMF.a1 =  -2*cos(w0);
            plugin.filter_LMF.a2 =   1 - alpha*A;
            plugin.filter_LMF.b0 =   1 + alpha/A;
            plugin.filter_LMF.b1 =  -2*cos(w0);
            plugin.filter_LMF.b2 =   1 - alpha/A;
        end
%% Low Shelf

        function set.LowShelf_freq(plugin, val)
            plugin.LowShelf_freq = val;
            update_LowShelf(plugin);
        end
        
        function set.LowShelf_gain(plugin, val)
            plugin.LowShelf_gain = val;
            update_LowShelf(plugin);
        end
        
        function update_LowShelf(plugin)
            Q=0.5;
            f0=plugin.LowShelf_freq;
            gain = plugin.LowShelf_gain;
            w0=2*pi*f0/plugin.fs;
            alpha=sin(w0)/(2*Q);
            A=sqrt(db2mag(gain));
            
            plugin.filter_LS.a0 =   A * (  (A+1)  -  (A-1)*cos(w0) + 2*sqrt(A)*alpha);
            plugin.filter_LS.a1 = 2*A * (  (A-1)  -  (A+1)*cos(w0)                  );
            plugin.filter_LS.a2 =   A * (  (A+1)  -  (A-1)*cos(w0) - 2*sqrt(A)*alpha);
            plugin.filter_LS.b0 =       (  (A+1)  +  (A-1)*cos(w0) + 2*sqrt(A)*alpha);
            plugin.filter_LS.b1 =  -2 * (  (A-1)  +  (A+1)*cos(w0)                  );
            plugin.filter_LS.b2 =       (  (A+1)  +  (A-1)*cos(w0) - 2*sqrt(A)*alpha);
        end
%% HPF

         function set.HPF_freq(plugin, val)
            plugin.HPF_freq = val;
            update_HPF(plugin);
            
        end


        function update_HPF(plugin)
            f0=plugin.HPF_freq;
            Q = 0.5;
            w0=2*pi*f0/plugin.fs;
            alpha=sin(w0)/(2*Q);

            
            plugin.filter_HPF.a0 =  (1 + cos(w0))/2;
            plugin.filter_HPF.a1 = -(1 + cos(w0));
            plugin.filter_HPF.a2 =  (1 + cos(w0))/2;
            plugin.filter_HPF.b0 =   1 + alpha;
            plugin.filter_HPF.b1 =  -2*cos(w0);
            plugin.filter_HPF.b2 =   1 - alpha;
        end
        
        function set.BYPASS(plugin, val)
            plugin.BYPASS = val;
        end
    end
end

