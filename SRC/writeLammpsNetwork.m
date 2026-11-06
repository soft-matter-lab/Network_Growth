function info = writeLammpsNetwork(nodes, bonds, activeEnds, cfg)
% Export a relaxed state in the NetworkGen bpm/poly + local/density formats.
% Writes network.data, bond.table, network.localdensity.table in cfg.outputDir.
% Direct calls always export; saveLammps controls only the automatic run hook.
% Kuhn count defaults to completed chains + loops + dangling chains.
% No free reservoir particles or explicit dangling tips are exported.
cfg = mergeConfig(defaultNetworkAConfig(),cfg);
if isempty(cfg.outputDir)
    error('NetworkGrowth:LammpsOutputDir', 'Set outputDir for LAMMPS export.');
end
validateattributes(nodes,{'numeric'},{'real','finite','2d','nonempty','ncols',3});
validateattributes(bonds,{'numeric'},{'real','finite','2d','ncols',5});
validateattributes(activeEnds,{'numeric'},{'real','finite','2d','ncols',4});
validateattributes(cfg.L,{'numeric'},{'real','finite','scalar','positive'});
validateattributes(cfg.b,{'numeric'},{'real','finite','scalar','positive'});
validateattributes(cfg.lammpsKLD,{'numeric'},{'real','finite','scalar','nonnegative'});
validateattributes(cfg.lammpsNRho,{'numeric'},{'real','finite','scalar','integer','>=',2});
validateattributes(cfg.lammpsRhoRange,{'numeric'},{'real','finite','vector','numel',2,'nonnegative'});
validateattributes(cfg.lammpsZThickness,{'numeric'},{'real','finite','scalar','positive'});
validateattributes(cfg.lammpsIncludeDanglingKuhn,{'logical','numeric'},{'scalar','binary'});
policy = validatestring(cfg.lammpsPeriodicPolicy,{'warn','error'});
if cfg.lammpsRhoRange(2) <= cfg.lammpsRhoRange(1)
    error('NetworkGrowth:LammpsDensityRange', 'lammpsRhoRange must increase.');
end
nAtoms = size(nodes,1); % actual formed crosslinkers, never the free A reservoir
if ~isempty(bonds)
    validateattributes(bonds(:,1:2),{'numeric'},{'integer','>=',1,'<=',nAtoms});
    validateattributes(bonds(:,3),{'numeric'},{'integer','positive'});
    validateattributes(bonds(:,4:5),{'numeric'},{'integer'});
end
if ~isempty(activeEnds)
    validateattributes(activeEnds(:,3),{'numeric'},{'integer','>=',1,'<=',nAtoms});
    validateattributes(activeEnds(:,4),{'numeric'},{'integer','nonnegative'});
end

loopRows = find(bonds(:,1) == bonds(:,2));
exportRows = find(bonds(:,1) ~= bonds(:,2));
exportBonds = bonds(exportRows,:);
completedKuhn = sum(bonds(:,3)); % includes self-loops before export filtering
danglingKuhn = sum(activeEnds(:,4)); % bare arms contribute zero; free B is absent
totalKuhn = completedKuhn + double(cfg.lammpsIncludeDanglingKuhn)*danglingKuhn;
areaPerCrosslinker = totalKuhn/nAtoms * pi*cfg.b^2/4;
rRepulsive = sqrt(areaPerCrosslinker/pi);
info = struct('status','pending','dataFile','','bondTableFile','', ...
    'localDensityFile','','atomCount',nAtoms,'bondCount',numel(exportRows), ...
    'exportedBondRows',exportRows,'omittedSelfLoopRows',loopRows, ...
    'completedKuhn',completedKuhn,'danglingKuhn',danglingKuhn, ...
    'includeDanglingKuhn',logical(cfg.lammpsIncludeDanglingKuhn), ...
    'totalKuhn',totalKuhn,'kuhnPerCrosslinker',totalKuhn/nAtoms, ...
    'areaPerCrosslinker',areaPerCrosslinker,'rRepulsive',rRepulsive, ...
    'RLower',0.8*rRepulsive,'RUpper',4*rRepulsive,'pairRepulsionCutoff',2*rRepulsive, ...
    'rho0',0.8*4^2,'kLD',cfg.lammpsKLD,'nRho',cfg.lammpsNRho, ...
    'rhoRange',cfg.lammpsRhoRange,'nonMinimumImageBondRows',zeros(0,1));

% bpm/poly rejects an empty bond table. There is no transferable network in
% this case; report the skipped export without failing a valid jammed run.
if isempty(exportRows)
    info.status = 'skipped_no_interjunction_bonds';
    warning('NetworkGrowth:LammpsEmptyNetwork', ...
        'LAMMPS export skipped: no bonds between distinct junctions in %s.',cfg.outputDir);
    return;
end

% BPM initializes bonds using the closest periodic image. A data file cannot
% encode the independent edge image offsets of a non-minimum-image strand.
xy = mod(nodes(:,1:2),cfg.L);
p = exportBonds(:,1); q = exportBonds(:,2);
originalVector = nodes(q,1:2) + exportBonds(:,4:5)*cfg.L - nodes(p,1:2);
minimumVector = xy(q,:) - xy(p,:);
minimumVector = minimumVector - cfg.L*round(minimumVector/cfg.L);
mismatch = any(abs(originalVector-minimumVector) > 1e-9*cfg.L,2);
info.nonMinimumImageBondRows = exportRows(mismatch);
if any(mismatch)
    message = sprintf(['%d strands use non-minimum periodic images. LAMMPS bpm/poly ' ...
        'will use their minimum-image vectors; see nonMinimumImageBondRows in results.lammps.'],sum(mismatch));
    if strcmp(policy,'error')
        error('NetworkGrowth:LammpsPeriodicImages','%s',message);
    end
    warning('NetworkGrowth:LammpsPeriodicImages','%s',message);
end
if any(vecnorm(minimumVector,2,2) < 1e-10)
    error('NetworkGrowth:LammpsZeroLengthBond', ...
        'Distinct bonded junctions coincide periodically; bpm/poly requires nonzero reference lengths.');
end
if ~isempty(loopRows)
    warning('NetworkGrowth:LammpsSelfLoops', ...
        ['Omitting %d self-loops from LAMMPS bonds (atom-to-itself bonds are invalid). ' ...
         'Their Kuhn segments remain in the density-radius total.'],numel(loopRows));
end

names = {'network.data','bond.table','network.localdensity.table'};
files = cellfun(@(name)fullfile(cfg.outputDir,name),names,'UniformOutput',false);
for k = 1:numel(files)
    if isfile(files{k})
        error('NetworkGrowth:ExistingLammpsFiles','Refusing to overwrite %s.',files{k});
    end
end
if ~isfolder(cfg.outputDir), mkdir(cfg.outputDir); end
% Stage all three before making any visible as completed output.
partial = cellfun(@(file)[file '.part'],files,'UniformOutput',false);
cleanup = onCleanup(@()removePartialFiles(partial));
writeData(partial{1},xy,exportBonds,cfg,info);
writeBonds(partial{2},exportBonds,cfg.b);
writeDensity(partial{3},cfg,info);
for k = 1:numel(files)
    [ok,msg] = movefile(partial{k},files{k});
    if ~ok, error('NetworkGrowth:LammpsWrite','%s',msg); end
end
info.dataFile = files{1};
info.bondTableFile = files{2};
info.localDensityFile = files{3};
info.status = 'written';
if any(mismatch), info.status = 'written_with_image_changes'; end
fprintf('Saved LAMMPS files in %s\n',cfg.outputDir);
fprintf('LAMMPS: N_Kuhn = %g, formed crosslinkers = %d, r_repulsive = %.10g, rc = %.10g\n', ...
    info.totalKuhn,nAtoms,info.rRepulsive,info.pairRepulsionCutoff);
end

function writeData(file,xy,bonds,cfg,info)
fid = openOutput(file);
cleanup = onCleanup(@()fclose(fid));
fprintf(fid,'Network Growth saved state; use dimension 2 and periodic x/y\n\n');
fprintf(fid,'%d atoms\n%d bonds\n1 atom types\n1 bond types\n\n',size(xy,1),size(bonds,1));
fprintf(fid,'0 %.17g xlo xhi\n0 %.17g ylo yhi\n',cfg.L,cfg.L);
fprintf(fid,'%.17g %.17g zlo zhi\n\n',-cfg.lammpsZThickness/2,cfg.lammpsZThickness/2);
fprintf(fid,'# r_repulsive %.17g; pair rc %.17g; N_Kuhn %.17g\n', ...
    info.rRepulsive,info.pairRepulsionCutoff,info.totalKuhn);
fprintf(fid,'# omitted_self_loops %d; non_minimum_image_bonds %d\n\n', ...
    numel(info.omittedSelfLoopRows),numel(info.nonMinimumImageBondRows));
fprintf(fid,'Atoms # bpm/sphere\n\n');
% NetworkGen convention: molecule=1, type=1, diameter=1, density=1, z=0.
fprintf(fid,'%d 1 1 1 1 %.17g %.17g 0\n',[(1:size(xy,1))',xy]');
fprintf(fid,'\nBonds\n\n');
fprintf(fid,'%d 1 %d %d\n',[(1:size(bonds,1))',bonds(:,1:2)]');
checkWrite(fid,file);
end

function writeBonds(file,bonds,b)
fid = openOutput(file);
cleanup = onCleanup(@()fclose(fid));
fprintf(fid,'# Chain stats: bondID atom1 atom2 N_Kuhn b (bpm/poly)\n\nKEY\nN %d\n\n',size(bonds,1));
fprintf(fid,'%d %d %d %d %.17g\n',[(1:size(bonds,1))',bonds(:,1:3),repmat(b,size(bonds,1),1)]');
checkWrite(fid,file);
end

function writeDensity(file,cfg,info)
fid = openOutput(file);
cleanup = onCleanup(@()fclose(fid));
fprintf(fid,'# Network Growth local/density: F(rho) = kLD*(rho-rho0)^2\n');
fprintf(fid,'# r_repulsive %.17g; kLD %.17g; rho0 %.17g\n',info.rRepulsive,info.kLD,info.rho0);
fprintf(fid,'1 %d\n\n',cfg.lammpsNRho);
fprintf(fid,'%.17g %.17g\n1\n1\n',info.RLower,info.RUpper);
lo = cfg.lammpsRhoRange(1); hi = cfg.lammpsRhoRange(2);
fprintf(fid,'%.17g %.17g %.17g\n',lo,hi,(hi-lo)/(cfg.lammpsNRho-1));
% LAMMPS reads exactly N_rho values and computes drho=(max-min)/(N_rho-1).
% The reference generator used N_rho+1 samples, inconsistent with its header.
rho = linspace(lo,hi,cfg.lammpsNRho)';
fprintf(fid,'%.17g\n',cfg.lammpsKLD*(rho-info.rho0).^2);
fprintf(fid,'\n');
checkWrite(fid,file);
end

function fid = openOutput(file)
[fid,msg] = fopen(file,'w');
if fid < 0, error('NetworkGrowth:LammpsWrite','Cannot open %s: %s',file,msg); end
end

function checkWrite(fid,file)
[msg,err] = ferror(fid);
if err ~= 0, error('NetworkGrowth:LammpsWrite','Cannot write %s: %s',file,msg); end
end

function removePartialFiles(files)
for k = 1:numel(files)
    if isfile(files{k}), delete(files{k}); end
end
end
