function main(varargin)
%main - The interface of PlatEMO.
%
%   main() displays the GUI of PlatEMO.
%
%   main('-Name',Value,'-Name',Value,...) runs one algorithm on a problem
%   with the specified parameter setting.
%
% All the acceptable properties:
%   '-N'            <positive integer>  population size
%   '-M'            <positive integer>  number of objectives
%   '-D'            <positive integer>  number of variables
%	'-algorithm'    <function handle>   algorithm function
%	'-problem'      <function handle>   problem function
%	'-evaluation'   <positive integer>  maximum number of evaluations
%   '-run'          <positive integer>  run number
%   '-save'         <integer>           number of saved populations
%   '-outputFcn'	<function handle>   function invoked after each generation
%
%   Example:
%       main()
%
%   displays the GUI of PlatEMO.
%
%       main('-algorithm',@ARMOEA,'-problem',@DTLZ2,'-N',200,'-M',10)
%
%   runs AR-MOEA on 10-objective DTLZ2 with a population size of 200.
%
%       main('-algorithm',{@KnEA,0.4},'-problem',{@WFG4,6})
%
%   runs KnEA on WFG4, and sets the parameters in KnEA and WFG4.
%
%       for i = 1 : 10
%           main('-algorithm',@RVEA,'-problem',@LSMOP1,'-run',i,'-save',5)
%       end
%
%   runs RVEA on LSMOP1 for 10 times, and each time saves 5 populations to
%   a file in PlatEMO/Data/RVEA.

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

   cd(fileparts(mfilename('fullpath')));
   addpath(genpath(cd));
   clear;
   close all;
   warning('off')
   
   %Global = GLOBAL('-algorithm',@EAGO,'-problem',@LSMOP9,'-N',276,'-M',10,'-D',500);


   
   for i = 1:7
       problem = sprintf("%s%d",'@DTLZ',i);


       
       Global = GLOBAL('-algorithm',@EAGO,'-problem',eval(problem),'-N',276,'-M',10,'-D',1000);
       

       Global.Start();
       close all
   end
	%Global.Start();


