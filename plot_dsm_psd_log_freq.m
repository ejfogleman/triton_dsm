function plot_dsm_psd_log_freq( datvar, delta, Fs, d, pltres )
% plot_dsm_psd_log_freq.m
%	Plot psd in dB; scale based on spectrum; log x axis
%       Provides 2x scaling to give V^2/Hz for + freq.
%
%	Usage:
%	plot_dsm_psd_log_freq( datvar, delta, Fs, d, pltres )
%
%	Where:
%	datvar(:,1) = normalized frequency vector (cycles/sample)
%	datvar(:,2) = psd (codes^2/(cycles/sample))
%	datvar(:,3) = res psd (codes^2/(cycles/sample))
% 	delta = lsb size (V)
% 	Fs = sample rate (Hz)
%       d = dynamic range of vertical axis
% 	pltres:  0=>plot both, 1=>plot sig only, 2=>plot res onlya

datlen = length( datvar(:,1) );

if ( Fs < 10e6 )
	fa = datvar( 1:datlen/2, 1 )*Fs/1000;  % scale in kHz; truncate to [0,0.5] 
else
	fa =  datvar( 1:datlen/2, 1 )*Fs/1e6;	% scale in MHz
end

dBSxx = 10*log10( 2*(delta^2/Fs)*datvar( 1:datlen/2, 2 ) );  % this scaling gives V^2/Hz
dBSrr = 10*log10( 2*(delta^2/Fs)*datvar( 1:datlen/2, 3 ) );  % this scaling gives V^2/Hz
h = axes;

if ( pltres == 2 )
	% res only
	semilogx( fa, dBSrr, 'b-' );
elseif ( pltres == 1 )
	% sig only
	semilogx( fa, dBSxx, 'b-' );
else
	% sig and res
	semilogx( fa, dBSrr, 'r-', fa, dBSxx, 'b-' );
end	

set( h, 'FontSize', 15 );

% find maximum over the band of interest for scaling y-axis
dBSxx_max = max( dBSxx );
dBSxx_max = 10*( ceil( dBSxx_max/10 ) );

if ( Fs < 10e6 )
	axis( [ 0.5e-5*(Fs/1000) 0.5*(Fs/1000) -d+dBSxx_max dBSxx_max ] );
xlabel('Frequency (kHz)');
else
	axis( [ 0.5e-5*(Fs/1e6) 0.5*(Fs/1e6) -d+dBSxx_max dBSxx_max ] );
	xlabel('Frequency (MHz)');
end
	
grid; 

ylabel('PSD (dBV/Hz)');


