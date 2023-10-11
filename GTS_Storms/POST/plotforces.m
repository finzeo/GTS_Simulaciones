t=csvread("user_data_forces.csv"); plot(t(:,3)); ylim([0 1]); grid; set(gca, 'ytick', 0:0.05:1);
