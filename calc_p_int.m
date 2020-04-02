
function [p_int]= calc_p_int(Q,W)

% Calculatig p_con: the number of rating intersections within observed data range
% as a percentage of the total number of rating intersections

% Experimenting the possibility of using estimated p_con as a priori of
% AMHG strength.

% J. Wang June 2, 2015





Q_dimension = size(Q);
n_sections = Q_dimension(2); %the number of cross sections along this river segment
fake_b_loga = zeros(n_sections, 2); %initiate a matrix to store rating coefficients
  
minW=10e6;
maxW=0;
% Loop through each cross section:
%figure;
for i_section = 1:n_sections
    Q_this_section = Q{i_section};
    W_this_section = W{i_section};
    
    if min(log10(W_this_section)) < minW
    minW=min(log10(W_this_section));
    end
     if max(log10(W_this_section)) > maxW
    maxW=max(log10(W_this_section));
    end

    
    %re order W_this_section
    [W_this_section, sort_ind] = sort(W_this_section);
    Q_this_section = Q_this_section(sort_ind);
    
    %Re-scale log(Q) to be [0, 10]
    fQ_min = 0; 
    fQ_max = 10;
    log_W_this_section = log10(W_this_section);
    fW_min = mean(log_W_this_section(1:2));
    fW_max = mean(log_W_this_section((length(log_W_this_section)-1):length(log_W_this_section))); 
     fQ = (log_W_this_section-fW_min).*(fQ_max-fQ_min)./(fW_max-fW_min); 


    %Rescaled pseudo AHG rating curves:
    X = [ones(size(fQ)) fQ];
    fake_coef = polyfit(fQ,log_W_this_section,1);
    fake_b_loga(i_section,:)= fake_coef; % store [b, log10(a)] for this cross-sectional rating



end


fake_b = fake_b_loga(:,1); 
fake_loga = fake_b_loga(:,2); 


%Rescaled pseudo AHG rating curves
n_line = size(fake_b_loga);
fake_XY_intersections = [];
for i_line = 1:n_line(1)
    for j_line = 1:n_line(1)
        if j_line>i_line
            p1 = [fake_b(i_line), fake_loga(i_line)]; p2 = [fake_b(j_line), fake_loga(j_line)];
            %calculate this intersection
            x_intersect = fzero(@(x) polyval(p1-p2,x),3);
            y_intersect = polyval(p1,x_intersect);
            fake_XY_intersections=[fake_XY_intersections, [x_intersect;y_intersect]];
            %hold on; plot(x_intersect, y_intersect, 'rO', 'MarkerSize',2)
        end
    end
end
fake_logQ_intersected = fake_XY_intersections(1,:); fake_logW_intersected = fake_XY_intersections(2,:);

%
%check the number of intersections within observation ranges:
[fake_within_range_indices]= ...
    find( (fake_logQ_intersected >= fQ_min & fake_logQ_intersected <= fQ_max));% & (fake_logW_intersected > minW & fake_logW_intersected < maxW ));
p_int = length(fake_within_range_indices)/length(fake_logQ_intersected);
