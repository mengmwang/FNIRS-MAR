function h=mar_algo()
    h.diff = @cal_sigdiff;
    h.rb = @reduce_basis;
    h.estimate = @estimate_iter;
    h.ls = @ls_estimate;
    h.overlapadd = @overlapadd;
	h.artefact = @add_artefact;
end

function yt = cal_sigdiff(y,dis)
% function used to calculate signal difference
% Input: dis: displacement
y_1 = y((dis+1):end,:);
y_0 = y(1:(end-dis),:);
yt = y_1 - y_0;
end

function [basis_r,basis_rt] = reduce_basis(y, coeff)
% function used to calculate reduced basis functions (Br and Br_tilde)
N = length(y);
dct_mat = dctmtx(N);
basis_mat = dct_mat';

basis_r = zeros(length(basis_mat),coeff);
basis_rt = zeros(N-1,coeff);

theta = basis_mat'*y;

[~,ind] = sort(abs(theta));
index = flipud(ind);
princ = index(1:coeff); 

for c = 1:1:coeff
    basis_r(:,c) = basis_mat(:,princ(c));   
end
for k = 1:1:N-1
    basis_rt(k,:) = basis_r(k+1,:) - basis_r(k,:);
end
end

function param_hat = estimate_iter(H,y,param_0,alpha)
% function implementation by NS
% H: regressor matrix
% y: observation vector
% param_0: initial estimation pinv(H'*H)*H'*y
% alpha: tuning parameter (scalar between 0 and 1). alpha=0.1

param_hat = param_0;
for i = 1:100
    w = (y - H*param_hat);
%     figure; plot(w)
    w = exp(-alpha*w.^2);
    w = w/sum(w);
    w = diag(w);
    param_hat = pinv(H'*w*H)*H'*w*y;
end
end

function lse = ls_estimate(X,y)
% least square estimator of X*b = y.
% uses pseudo-inverse of X.
%lse = X\y;
lse = pinv(X'*X)*X'*y;
end

function [x,zo]=overlapadd(f,win,inc)
%OVERLAPADD join overlapping frames together X=(F,WIN,INC)
%
% Usage for frequency-domain processing:
%       S=...;                              % input signal
%       OV=2;                               % overlap factor of 2 (4 is also often used)
%       INC=20;                             % set frame increment in samples
%       NW=INC*OV;                          % DFT window length
%       W=sqrt(hamming(NW,'periodic'));     % omit sqrt if OV=4
%       W=W/sqrt(sum(W(1:INC:NW).^2));      % normalize window
%       F=rfft(enframe(S,W,INC),NW,2);      % do STFT: one row per time frame, +ve frequencies only
%       ... process frames ...
%       X=overlapadd(irfft(F,NW,2),W,INC);  % reconstitute the time waveform (omit "X=" to plot waveform)
%
% Inputs:  F(NR,NW) contains the frames to be added together, one
%                   frame per row.
%          WIN(NW)  contains a window function to multiply each frame.
%                   WIN may be omitted to use a default rectangular window
%                   If processing the input in chunks, WIN should be replaced by
%                   ZI on the second and subsequent calls where ZI is the saved
%                   output state from the previous call.
%          INC      gives the time increment (in samples) between
%                   succesive frames [default = NW].
%
% Outputs: X(N,1) is the output signal. The number of output samples is N=NW+(NR-1)*INC.
%          ZO     Contains the saved state to allow a long signal
%                 to be processed in chunks. In this case X will contain only N=NR*INC
%                 output samples.
%
%	   Copyright (C) Mike Brookes 2009
%      Version: $Id: overlapadd.m 2470 2012-11-02 15:27:24Z dmb $
%
%   VOICEBOX is a MATLAB toolbox for speech processing.
%   Home page: http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You can obtain a copy of the GNU General Public License from
%   http://www.gnu.org/copyleft/gpl.html or by writing to
%   Free Software Foundation, Inc.,675 Mass Ave, Cambridge, MA 02139, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nr,nf]=size(f);            % number of frames and frame length
if nargin<2
    win=nf;                 % default increment
end
if isstruct(win)
    w=win.w;
    if ~numel(w) && length(w)~=nf
        error('window length does not match frames size');
    end
    inc=win.inc;
    xx=win.xx;
else
    if nargin<3
        inc=nf;
    end
    if numel(win)==1 && win==fix(win) && nargin<3       % win has been omitted
        inc=win;
        w=[];
    else
        w=win(:).';
        if length(w)~=nf
            error('window length does not match frames size');
        end
        if all(w==1)
            w=[];
        end
    end
    xx=[];      % partial output from previous call is null
end
nb=ceil(nf/inc);        % number of overlap buffers
no=nf+(nr-1)*inc;       % buffer length
z=zeros(no,nb);                      % space for overlapped output speech
if numel(w)
    z(repmat(1:nf,nr,1)+repmat((0:nr-1)'*inc+rem((0:nr-1)',nb)*no,1,nf))=f.*repmat(w,nr,1);
else
    z(repmat(1:nf,nr,1)+repmat((0:nr-1)'*inc+rem((0:nr-1)',nb)*no,1,nf))=f;
end
x=sum(z,2);
if ~isempty(xx)
    x(1:length(xx))=x(1:length(xx))+xx;     % add on leftovers from previous call
end
if nargout>1            % check if we want to preserve the state
    mo=inc*nr;          % completed output samples
    if no<mo
        x(mo,1)=0;
        zo.xx=[];
    else
        zo.xx=x(mo+1:end);
        zo.w=w;
        zo.inc=inc;
        x=x(1:mo);
    end
elseif ~nargout
    if isempty(xx)
        k1=nf-inc;  % dubious samples at start
    else
        k1=0;
    end
    k2=nf-inc;      % dubious samples at end
    plot(1+(0:nr-1)*inc,x(1+(0:nr-1)*inc),'>r',nf+(0:nr-1)*inc,x(nf+(0:nr-1)*inc),'<r', ...
        1:k1+1,x(1:k1+1),':b',k1+1:no-k2,x(k1+1:end-k2),'-b',no-k2:no,x(no-k2:no),':b');
    xlabel('Sample Number');
    title(sprintf('%d frames of %d samples with %.0f%% overlap = %d samples',nr,nf,100*(1-inc/nf),no));
end
end

% tools add simulated artefacts
function y = add_artefact(x,n1,n2,s1,s2)
% generate added motion artefacts
% x: data
% n1/n2: no of type1/type2 artefacts
% s1/s2: scaling factors of type1/type2 artefacts
% y: signal with added motion artefacts

[total_t, noc] = size(x);
ma_t1 = zeros(total_t, noc);
ma_type_1 = zeros(total_t, noc);
ma_type_2 = zeros(total_t, noc);

t_1 = randi([1 total_t], 1, n1); % the starting point of motion artifacts
t_2 = randi([1 total_t], 1, n2);

r_1 = ones(1, n1);
r_2 = ones(1, n2);

for m = 1:n1
    rd_1(1,m) = randi([1 total_t-t_1(1,m)]);
    ma_t1(t_1(1,m):t_1(1,m)+rd_1(1,m),1) = r_1(1, m);
    ma_type_1 = ma_type_1 + ma_t1;
    ma_t1 = zeros(total_t, noc);
end

for n = 1:n2
    rd_2(1,n) = randi([1 total_t-t_2(1,n)]);
    ma_type_2(t_2(1,n),1) = r_2(1, n);
    ma_type_2(t_2(1,n)+rd_2(1,n),1) = r_2(1, n);
end

% change the amplitude
m1 = ma_type_1 * s1; 
m2 = ma_type_2 * s2;

y=x+m1+m2;
end