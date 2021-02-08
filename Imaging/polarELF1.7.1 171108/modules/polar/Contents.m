% Contents of folder "POLAR"
% Contains functions related to calculating polarisataion light fields. See polar.m for the main calculation routine and examples of plotting functions.
%
% Dependencies: No toolbox functions required.
%
% Files:
%   polar                           - Example script
%   polar_main1_HDRscenes           - Step 1, unwarps and calibrates images
%   polar_main2_filter              - Step 2, filters images to 2/4/8/16 degrees
%   polar_main3_stokes              - Step 3, calculatea Stokes' vectors and Int/AoP/DoLP matrices
%   polar_main4_stitch              - Step 4, stitches together different compass directions into a single image
%   polar_aop3d                     - Calculates 3d Polarisation vectors
%   polar_gridsphere                - Produces a grid that regularly covers a sphere; based on FileExchange GridSphere
%   im_stokes                       - Calculates Stokes vectors for a single image stack
%
% Suibfolder plot:
%   polar_plot_int                  - Plots intensity image.
%   polar_plot_dolp                 - Plots degree of linear polarisation as gray scale.
%   polar_plot_aop                  - Plots angle of polarisation as colour.
%   polar_plot_int3d                - Plots intensity image on a sphere.
%   polar_plot_dolp3d               - Plots degree of linear polarisation as colour on a sphere.
%   polar_plot_aop3d                - Plots angle of polarisation as vectors on a sphere.
%   polar_plot_sphere               - Plots a monochromatic sphere.
