function [nodes, rmsSpeed, it] = equilibrateNetwork(nodes, bonds, activeEnds, cfg, rEq, viz)
% Mechanical relaxation of the current connectivity.
% Supports two schemes selected by cfg.eqIntegrator:
%   - "damped2nd" (default): m*x'' + gamma*x' = F
%   - "overdamped":          x' = mobility*F

if nargin < 6, viz = []; end

N = size(nodes,1);
it = 0;
rmsSpeed = 0;
if N < 1 || isempty(bonds)
    return;
end

integrator = lower(strtrim(cfg.eqIntegrator));
if isempty(integrator)
    integrator = 'damped2nd';
end

vState = zeros(N,2); % second-order state (ignored in overdamped mode)

for it = 1:cfg.maxEq
    F = networkForces(nodes, bonds, cfg.b, cfg.kT, cfg.kRep, rEq, cfg.L, ...
        cfg.repelBondedPairs, cfg.useRepulsionBinning, cfg.repulsionFullPairThreshold);

    switch integrator
        case 'overdamped'
            v = cfg.mobility * F;
            speed = vecnorm(v, 2, 2);
            rmsSpeed = sqrt(mean(speed.^2));
            if rmsSpeed < cfg.tolV
                break;
            end
            vmax = max(speed);
            dt = min(cfg.dtEq, cfg.maxMove / max(vmax, 1e-14));
            nodes(:,1:2) = nodes(:,1:2) + v * dt;

        otherwise
            % Damped second-order step:
            %   v_{n+1} = v_n + dt*(F/m - gamma*v_n)
            %   x_{n+1} = x_n + dt*v_{n+1}
            massEq = max(cfg.massEq, 1e-12);
            gammaEq = max(cfg.gammaEq, 0);
            accelFree = F / massEq;

            vTrial = dampedVelocityUpdate(vState, accelFree, cfg.dtEq, gammaEq);
            vmax = max(vecnorm(vTrial, 2, 2));
            dt = min(cfg.dtEq, cfg.maxMove / max(vmax, 1e-14));
            vState = dampedVelocityUpdate(vState, accelFree, dt, gammaEq);

            % If inertia starts to oppose the current descent direction, zero it.
            if sum(vState(:) .* F(:)) < 0
                vState(:) = 0;
            end

            nodes(:,1:2) = nodes(:,1:2) + dt * vState;

            speed = vecnorm(vState, 2, 2);
            rmsSpeed = sqrt(mean(speed.^2));
            if rmsSpeed < cfg.tolV
                break;
            end
    end

    if cfg.showLivePlot && ~isempty(viz) && mod(it, cfg.eqPlotEvery) == 0
        drawMotion(nodes, bonds, cfg.L, activeEnds, viz, cfg.fA, cfg.b);
        drawnow;
    end
end

% Final speed metric at returned coordinates.
F = networkForces(nodes, bonds, cfg.b, cfg.kT, cfg.kRep, rEq, cfg.L, ...
    cfg.repelBondedPairs, cfg.useRepulsionBinning, cfg.repulsionFullPairThreshold);
if strcmp(integrator, 'overdamped')
    v = cfg.mobility * F;
else
    v = vState;
end
rmsSpeed = sqrt(mean(vecnorm(v, 2, 2).^2));
end

function vNext = dampedVelocityUpdate(vNow, accelFree, dt, gammaEq)
if gammaEq > 0
    damp = exp(-gammaEq * dt);
    vNext = damp .* vNow + (1 - damp) / gammaEq .* accelFree;
else
    vNext = vNow + dt * accelFree;
end
end
