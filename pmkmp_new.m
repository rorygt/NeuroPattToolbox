function map=pmkmp_new(scheme, n)
% PMKMP_NEW Makes perceptually balanced colormaps with rainbow-like colors
%
% ARGUMENTS:
%        scheme -- Colormap scheme. Can be linear schemes 'isol' or 
%                  'cubicl', or circular schemes 'ostwald', 'ostwald_o', or
%                  'edge'. Default 'ostwald_o'.
%        n      -- Number of elements in colormap, cannot be greater than
%                  256. Default is the number of colormap elements in the
%                  active figure.
%
% OUTPUT: 
%        map    -- Nx3 colormap defined by RGB triplets.
%
% USAGE:
%{
    schemes = {'isol', 'cubicl', 'ostwald', 'ostwald_o', 'edge'}
    linearData = peaks(100);
    circularData = mod(repmat(linspace(0,10,50), 50, 1), 2*pi);
    data = {linearData, circularData};
    for idata = 1:length(data)
        for ischeme = 1:length(schemes)
            figure(idata)
            h = subplot(2,3,ischeme);
            imagesc(data{idata})
            colormap(h, pmkmp_new(schemes{ischeme}))
            colorbar
            title(schemes{ischeme})
        end
    end
%}
%
% This function primarily implements the PMKMP function written by Matteo
% Niccoli on the MATLAB File Exchange, but has been modified by Paul
% Martin to add a couple of new colormap schemes: 
%    'ostwald'   Colormap based on Ostwald hue trianges, suitable for use
%                with circular data
%    'ostwald_o' Variant of ostwald map to give less distinct but more
%                visually appealing colours (blue-pink rather than
%                purple-green)
%
%
% _________________________________________
% ORIGINAL DOCUMENTATION FOR PMKMP FUNCTION
% _________________________________________
% JUSTIFICATION: rainbow, or spectrum color schemes are considered a poor
% choice for scientific data display by many in the scientific community
% (see for example reference 1 and 2) in that they introduce artifacts 
% that mislead the viewer. "The rainbow color map appears as if it�s separated
% into bands of almost constant hue, with sharp transitions between hues. 
% Viewers perceive these sharp transitions as sharp transitions in the data,
% even when this is not the casein how regularly spaced (interval) data are
% displayed (quoted from reference 1). This submission is intended to share
% the results of my work to create more perceptually balanced, 
% rainbow-like color maps. Please see output arguments section for descriptions.
%
%
%   arguments: (input)
%   scheme - can be one of the following strings:
%     'IsoL'       Lab-based isoluminant rainbow with constant luminance L*=60
%                  For interval data displayed with external lighting
%
%     'CubicL'	   Lab-based rainbow scheme with cubic-law luminance (default)
%                  For interval data displayed without external lighting
%
%     'Edge'       Diverging Black-blue-cyan-white-yellow-red-black scheme
%                  For ratio data (ordered, constant scale, natural zero)  
%
%   n - scalar specifying number of points in the colorbar. Maximum n=256
%      If n is not specified, the size of the colormap is determined by the
%      current figure. If no figure exists, MATLAB creates one.
%
%
%   arguments: (output)
%   map - colormap of the chosen scheme
%   - IsoL is based on work in paper 2 in the reference section.
%     In both this paper and in several others this is indicated as the
%     best for displaying interval data with external lighting.
%     This is so as to allow the lighting to provide the shading to
%     highlight the details of interest. If lighting is combined with a
%     colormap that has its own luminance function associated - even as
%     simple as a linear increase this will confuse the viewer. The only 
%     difference from the paper is that I changed the value of constant 
%     luminance to L*=60 to make it brighter that the authors' example.
%
%   - CubicL too is based on some of the ideas in paper 2 in the 
%      reference section but rather than using a linearly increasing
%      L* function such as the one used by those authors, I am
%      using a compressive or cubic law function for the increase in 
%      L*.  L* ranges between 31 and 90 in the violet to yellowish 
%      portion of the colormap, then decreases to about 80 to get 
%      to the red (please refer to figure L_a_b_PlotsCubicL.png).
%      The choice to start at 31 was a matter of taste. 
%      I like having violet instead of black at the cold end of the
%      colormap. The latter choice was so as to have red and not
%      white at the warm end  of the colorbar, which is also a 
%      matter of personal taste. As a result,  there is an inversion in 
%      the L* trend, but I believe because it is a smooth one that
%      this is an acceptable compromise and the resulting
%      colormap is much of an improvement over the standard 
%      rainbow or spectrum schemes, which  typically have at least 3 sharp 
%      L* inversions. Please run CompareLabPlotsUsingColorspace.m or see
%      figures: L_plot_for_CubicL_colormap.png, L_plot_for_jet_colormap.png,
%      and L_plot_for_spectrum_colormap.png for a demonstration
%
%   - Edge is based on the Goethe Edge Colors described in the book in 
%     reference 3. In practice the colormap resembles a cold color map attached
%     to a warm color map. But the science behind it is rigorous and the
%     experimental work is based on is very intriguing to me: an alternative
%     to the Newtonian spectrum. This is not perceptually balanced in a
%     strict sense but because it does not have green it is perceptually
%     improved in a harmonious sense (refer to paper reference 10 for a review
%     of the concept of harmony in color visualization).
%
%
%   Example: 128-color rainbow with cubic-law luminance (default)
%     %  load mandrill;
%     %  imagesc(X);
%     %  colormap(pmkmp(128));
%     %  colorbar;
%   See files examples.m, examples1.m, and example2.m for more examples 
%   See files MakeLabPlotUsingColorspace.m and CompareLabPlotsUsingColorspace.m 
%   for some demonstrations
%
%
%   See also: JET, HSV, GRAY, HOT, COOL, BONE, COPPER, PINK, FLAG, PRISM,
%             COLORMAP, RGBPLOT
% 
%
%   Other submissions of interest
% 
%     Haxby color map
%     www.mathworks.com/matlabcentral/fileexchange/25690-haxby-color-map
% 
%     Colormap and colorbar utilities
%     www.mathworks.com/matlabcentral/fileexchange/24371-colormap-and-color
%     bar-utilities-sep-2009
% 
%     Lutbar
%     www.mathworks.com/matlabcentral/fileexchange/9137-lutbar-a-pedestrian-colormap-toolbarcontextmenu-creator
% 
%     usercolormap
%     www.mathworks.com/matlabcentral/fileexchange/7144-usercolormap
% 
%     freezeColors
%     www.mathworks.com/matlabcentral/fileexchange/7943
%
%
%     Bipolar Colormap
%     www.mathworks.com/matlabcentral/fileexchange/26026
%
%     colorGray
%     www.mathworks.com/matlabcentral/fileexchange/12804-colorgray
%
%     mrgb2gray
%     www.mathworks.com/matlabcentral/fileexchange/5855-mrgb2gray
%
%     CMRmap
%     www.mathworks.com/matlabcentral/fileexchange/2662-cmrmap-m
%
%     real2rgb & colormaps
%     www.mathworks.com/matlabcentral/fileexchange/23342-real2rgb-colormaps
%
%   Acknowledgements
% 
%     For input to do this researchI was inspired by: 
%     ColorSpiral - http://bsp.pdx.edu/Software/ColorSpiral.m
%     Despite an erroneous assumption about conversion/equivalence to 
%     grayscale (which CMRmap achieves correctly) the main idea is ingenious
%     and the code is well written. It also got me interested in perceptual
%     colormaps. See reference 5 for paper
%     
%     For function architecture and code syntax I was inspired by:
%     Light Bartlein Color Maps 
%     www.mathworks.com/matlabcentral/fileexchange/17555
%     (and comments posted therein)
% 
%     For idea on world topgraphy in examples.m I was inspired by:
%     Cold color map
%     www.mathworks.cn/matlabcentral/fileexchange/23865-cold-colormap
%
%     To generate the spectrum in examples1.m I used:
%     Spectral and XYZ Color Functions
%     www.mathworks.com/matlabcentral/fileexchange/7021-spectral-and-xyz-color-functions
%     
%     For Lab=>RGB conversions I used:
%     Colorspace transforamtions
%     www.mathworks.com/matlabcentral/fileexchange/28790-colorspace-transformations
%
%
%     For the figures in example 2 I used:
%     Shaded pseudo color
%     http://www.mathworks.cn/matlabcentral/fileexchange/14157-shaded-pseudo-color
%
%     For plots in CompareLabPlotsUsingColorspace.m I used:
%     cline
%     http://www.mathworks.cn/matlabcentral/fileexchange/14677-cline
%
%     For some ideas in general on working in Lab space:
%     Color scale
%     www.mathworks.com/matlabcentral/fileexchange/11037
%     http://blogs.mathworks.com/steve/2006/05/09/a-lab-based-uniform-color-scale/
%
%     A great way to learn more about improved colormaps and making colormaps:
%     MakeColorMap
%     www.mathworks.com/matlabcentral/fileexchange/17552
%     blogs.mathworks.com/videos/2007/11/15/practical-example-algorithm-development-for-making-colormaps/
%
%
%  References
% 
%     1)  Borland, D. and Taylor, R. M. II (2007) - Rainbow Color Map (Still) 
%         Considered Harmful
%         IEEE Computer Graphics and Applications, Volume 27, Issue 2
%         Pdf paper included in submission
% 
%     2)  Kindlmann, G. Reinhard, E. and Creem, S. Face-based Luminance Matching
%         for Perceptual Colormap Generation
%         IEEE - Proceedings of the conference on Visualization '02
%         www.cs.utah.edu/~gk/papers/vis02/FaceLumin.pdf
% 
%     3)  Koenderink, J. J. (2010) - Color for the Sciences
%         MIT press, Cambridge, Massachusset
% 
%     4)  Light, A. and Bartlein, P.J. (2004) - The end of the rainbow? 
%         Color schemes for improved data graphics.
%         EOS Transactions of the American Geophysical Union 85 (40)
%         Reprint of Article with Comments and Reply
%         http://geography.uoregon.edu/datagraphics/EOS/Light-and-Bartlein.pdf
% 
%     5)  McNames, J. (2006) An effective color scale for simultaneous color
%         and gray-scale publications
%         IEEE Signal Processing Magazine, Volume 23, Issue1
%         http://bsp.pdx.edu/Publications/2006/SPM_McNames.pdf
%
%     6)  Rheingans, P.L. (2000), Task-based Color Scale Design
%         28th AIPR Workshop: 3D Visualization for Data Exploration and Decision Making
%         www.cs.umbc.edu/~rheingan/pubs/scales.pdf.gz
% 
%     7)  Rogowitz, B.E. and  Kalvin, A.D. (2001) - The "Which Blair project":
%         a quick visual method for evaluating perceptual color maps. 
%         IEEE - Proceedings of the conference on Visualization �01
%         www.research.ibm.com/visualanalysis/papers/WhichBlair-Viz01Rogowitz_Kalvin._final.pdf
% 
%     8)  Rogowitz, B.E. and  Kalvin, A.D. - Why Should Engineers and Scientists
%         Be Worried About Color?
%         www.research.ibm.com/people/l/lloydt/color/color.HTM
% 
%     9)  Rogowitz, B.E. and  Kalvin, A.D. - How NOT to Lie with Visualization
%         www.research.ibm.com/dx/proceedings/pravda/truevis.htm
%
%     10) Wang, L. and Mueller,K (2008) - Harmonic Colormaps for Volume Visualization
%         IEEE/ EG Symposium on Volume and Point-Based Graphics
%         http://www.cs.sunysb.edu/~mueller/papers/vg08_final.pdf
%
%     11) Wyszecki, G. and Stiles W. S. (2000) - Color Science: Concepts and 
%         Methods, Quantitative Data and Formulae, 2nd Edition, John Wiley and Sons
% 
%
%  Author: Matteo Niccoli
%  e-mail address: matteo@matteoniccoli.com
%  Release: 1.02
%  Release date: October 13 2010


% error checking, defaults, valid schemes
%narginchk(0,2)
%nargoutchk(0,1)


if nargin<1
    scheme = 'ostwald_o';
end
if nargin<2
    n = size(get(gcf,'colormap'),1);
end
if n>256
    error('Maximum number of 256 points for colormap exceeded');
end

switch lower(scheme)
    case 'isol'
        baseMap = IsoL;
    case 'cubicl'
        baseMap = CubicL;
    case 'edge'
        baseMap = Edge;
    case 'ostwald'
        baseMap = Ostwald;
    case 'ostwald_o'
        baseMap = Ostwald_o;
    otherwise
        error(['Invalid scheme ' scheme])
end
idx1 = linspace(1,n,size(baseMap,1));
idx2 = 1:1:n;
warning off
map = interp1(idx1,baseMap,idx2,'cubic');
warning on

% Clean up dirty zeros
map(map < 0) = 0;
 end

function baseMap = Edge
baseMap =    [0 0 0;
              0 0 1;
              0 1 1;
              1 1 1;
              1 1 0;
              1 0 0
              0 0 0];
end

function baseMap = Ostwald

t  = 1 - logspace(-1,0,10)';
t2 = [t; flipud(t(2:end))];
t3 = [t2, circshift(t2,5),circshift(t2,-5)];
baseMap =    [t3; t3(1,:)];
end

function baseMap = Ostwald_o

t  = (.2:.1:1)';
t2 = [t; flipud(t(2:end))];
t3 = [t2, circshift(t2,10),circshift(t2,-10)];
baseMap =    [t3; t3(1,:)];
end

function baseMap = IsoL
baseMap =   [0.9102    0.2236    0.8997
             0.4027    0.3711    1.0000
             0.0422    0.5904    0.5899
             0.0386    0.6206    0.0201
             0.5441    0.5428    0.0110
             1.0000    0.2288    0.1631];
end
 
function baseMap = CubicL
 baseMap =  [0.4706         0    0.5216;
             0.5137    0.0527    0.7096;
             0.4942    0.2507    0.8781;
             0.4296    0.3858    0.9922;
             0.3691    0.5172    0.9495;
             0.2963    0.6191    0.8515;
             0.2199    0.7134    0.7225;
             0.2643    0.7836    0.5756;
             0.3094    0.8388    0.4248;
             0.3623    0.8917    0.2858;
             0.5200    0.9210    0.3137;
             0.6800    0.9255    0.3386;
             0.8000    0.9255    0.3529;
             0.8706    0.8549    0.3608;
             0.9514    0.7466    0.3686;
             0.9765    0.5887    0.3569];
end
