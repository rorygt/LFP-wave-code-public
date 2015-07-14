function [phi, v0, v_direction] = orderParameter(vx, vy)
% Function to calculate order parameter, from equation 1 in Collective
% Motion by Vicsek & Zafeiris.
% VX and VY must be 10x10xtime vectors giving the x and y components of the
% vector field.

n = numel(vx(:,:,1));

% Calculate average velocity magnitude
v0 = squeeze(mean(nanmean(sqrt(vx.^2 + vy.^2))));

% Calculate average normalized velocity magnitude
vx_sum = squeeze(sum(nansum(vx)));
vy_sum = squeeze(sum(nansum(vy)));
phi = 1./(n*v0) .* sqrt(vx_sum.^2 + vy_sum.^2);

% Calculate average velocity direction
if nargout > 2
    v_direction = atan2(vy_sum, vx_sum);
end

end