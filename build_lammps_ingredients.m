% Build LAMMPS inputs from ONE previously saved simulation. No growth is run.
% Edit inputFolder, then run this script. A dumps subfolder is also accepted.
addpath(fullfile(fileparts(mfilename('fullpath')),'SRC'));

inputFolder = fullfile(fileparts(mfilename('fullpath')),'runs_pkg_cleaned', ...
    'single_xA_0p0500_N_6714_20261026_122029_480'); % '' opens a folder picker

builderOptions = struct();
builderOptions.outputDir = ''; % '' = put the three LAMMPS files in the simulation folder
builderOptions.sourceFile = ''; % '' = auto; or select a specific saved MAT file
builderOptions.allowLatestSnapshot = false; % true permits an unfinished run's latest network dump
builderOptions.exportConfig = struct();
builderOptions.exportConfig.lammpsKLD = 0.414;
% Optional: lammpsNRho, lammpsRhoRange, lammpsZThickness, lammpsPeriodicPolicy.
% b and L come from the saved simulation. Loops + connected + dangling Kuhn
% segments count toward the radius; only distinct-node bonds go to LAMMPS.
% Existing LAMMPS files are preserved. Set outputDir to a fresh folder to rebuild.

if isempty(inputFolder)
    inputFolder = uigetdir(pwd,'Select one saved simulation folder');
    if isequal(inputFolder,0)
        disp('LAMMPS ingredients build cancelled.');
        return;
    end
end

lammpsBuild = buildLammpsIngredients(inputFolder,builderOptions);
disp(lammpsBuild);
