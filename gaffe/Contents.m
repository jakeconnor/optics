% The gaffe toolbox.
% Version 0.0.13 (20100527) 27-May-2010
%
% Author: Edward Grace <ej.grace@imperial.ac.uk>
%
% The GAFFE toolbox is a collection of routines for the Generalised
% Adaptive Fast Fourier Evolution (GAFFE) of N-dimensional, potentially 
% nonlinear partial differential equations.
%
% Given a generalised, possibly nonlinear, N dimensional wave
% equation that can be written in terms of linear and nonlinear
% operators with respect to a single evolutionary parameter this
% toolkit allows one to evolve this system adapting step-size and mesh-size 
% in a fully adaptive manner.
%
% Solver related functions
%  evolve        - Evolve a given field according to the operators.
%  evolvemeshset - Set mesh resizing convention for evolve.
%  fftspace      - Make coordinates for FFT ordered signals.
%
% Demos
%  gaffe_demo_linear   - A simple example of linear 1D gaussian dispersion.
%  gaffe_demo_soliton  - Follow a 1D simulation of a second-order soliton.
%  gaffe_demo_steepen  - Simulate nonlinear self-steepening in 1D.
%  gaffe_demo_gaussian - A 2D example of nonlinear self-focusing.
%  
% Operators
%  DefaultDiffraction    - Canonical diffraction operator in 2D.
%  DefaultDispersion     - Canonical dispersion (1D diffraction).
%  DefaultIdentity       - No effect, useful for disabling effects.
%  DefaultKerr           - Canonical (Kerr) self phase modulation.
%  DefaultSelfSteepening - Canonical 1D self-steepening operator.
%  DefaultDiffractionFibich - Diffraction operator for case study.  
%
% Solver support functions
%  structmerge           - Merge structures.
%
% Mathematical support functions
%  fftnpad               - Padding of N dimensional FFT ordered data.
%  fibonacci             - Generalised Fibonacci function.
%  ifibonacci            - Inverse Fibonacci function.
%  isfibonacci           - Determine if a number is a Fibonacci number.
%  nearest_humble        - Determine nearest humble number.
%  l2norm                - L2 norm of an N dimensional field.
%  width                 - Return width of a field along given dimension.
% 
% Gaussian beam functions
%  rtm_distance          - Generate a ray transfer matrix for distance.
%  rtm_lens              - Generate a ray transfer matrix for a lens.
%  wr2iq                 - Form an inverse Q parameter from W & R.
%  iq2iq                 - Transform an inverse Q parameter with a RTM.
%  iq2w0                 - Determine the waist of a gaussian.
%  iq2wr                 - Given an inverse Q return the spot size and R.
%  iq2z0                 - Return the distance to the waist for a given iq.
%  zr                    - Return the Rayleigh range for beam waist.
%  gaussian              - Evaluate a Gaussian given an inverse Q.
%
%
%gaffe - noun; 1. A clumsy social mistake.
%              2. A blunder.
%              3. Gaston Lagaffe - http://www.gastonlagaffe.com/


% $Author: graceej $ $Date: 2010/05/23 13:22:51 $
% $Revision: 1.9 $
% Copyright (c) 2010, Edward J. Grace
% All rights reserved.
 
% Redistribution and use in source and binary forms, with or 
% without modification, are permitted provided that the following 
% conditions are met:
 
%     * Redistributions of source code must retain the above 
%       copyright notice, this list of conditions and the following 
%       disclaimer.
%     * Redistributions in binary form must reproduce the above 
%       copyright notice, this list of conditions and the following 
%       disclaimer in the documentation and/or other materials 
%       provided with the distribution.
%     * Neither the name of the Imperial College London nor the 
%       names of its contributors may be used to endorse or 
%       promote products derived  this software without specific 
%       prior written permission.
 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND 
% CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, 
% INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS 
% BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED 
% TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
% DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON 
% ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR 
% TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF 
% THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF 
% SUCH DAMAGE.

