function plot_dsm_psd_linear_freq( datvar, delta, Fs, d, OSR, pltres )
% plot_dsm_psd_linear_freq.m
%	Plot psd in dB; scale based on spectrum; linear x axis
%       Provides 2x scaling to give V^2/Hz for + freq
%
%	Usage:
%	plot_dsm_psd_linear_freq( datvar, delta, Fs, d, OSR, pltres )
%
%	Where:
% 	datvar(:,1) = normalized frequency vector (cycles/sample)
%	datvar(:,2) = psd (codes^2/(cycles/sample))
%	datvar(:,3) = res psd (codes^2/(cycles/sample))
% 	delta = lsb size (V)
% 	Fs = sample rate (Hz)
%       OSR = oversampling ratio 
%       d = dynamic range of vertical axis
%	pltres:  0=>plot both, 1=>plot sig only, 2=>plot res only

if ( Fs < 10e6 )
	fa = datvar(:,1) * Fs/1000;  % make x axis scale kHz
else
	fa = datvar(:,1) * Fs/1e6;	% make x axis scale MHz 
end

dBSxx = 10*log10( 2*(delta^2/Fs)*datvar(:,2) );  % this scaling gives V^2/Hz
dBSrr = 10*log10( 2*(delta^2/Fs)*datvar(:,3) );  % this scaling gives V^2/Hz
clf;
h = axes;

if ( pltres == 2 )
	% res only
	plot( fa, dBSrr, 'b-' );
elseif ( pltres == 1 )
	% sig only
	plot( fa, dBSxx, 'b-' );
else
	% sig and res
	plot( fa, dBSrr, 'r-', fa, dBSxx, 'b-' );
end	

set( h, 'FontSize', 15 );

% find maximum over the band of interest for scaling y-axis
dBSxx_max = max( dBSxx( 1 : ( length( dBSxx ) / 2 / OSR ) ) ); 
dBSxx_max = 10*( ceil( dBSxx_max/10 ) );

% xtick_vector = [ 0 0.1 0.2 0.3 0.4 0.5 ] / OSR;

if ( Fs < 10e6 )
	axis( [ 0 0.5*(Fs/1000)/OSR -d+dBSxx_max dBSxx_max ] );
	% set( gca, 'xtick', xtick_vector ); 
	xlabel('Frequency (kHz)');
else
	axis( [ 0 0.5*(Fs)/1e6/OSR -d+dBSxx_max dBSxx_max ] );
	xlabel('Frequency (MHz)');
end
grid; 

ylabel('PSD (dBV/Hz)');


