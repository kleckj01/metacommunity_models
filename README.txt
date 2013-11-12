Metacommunity_models

The function metacommunity_input specifies parameter values and saves them in a file called input.mat

The function metacommunity_simulation loads parameters from the input.mat file and runs the simulations

The function plot_dynamics plots the temporal dynamics in individual patches



The following commands run everything:

metacommunity_input("in_test.mat")
load("in_test.mat")
metacommunity_simulation("in_test.mat","out_test.mat")
load("out_test.mat")
plot_dynamics(t,x)

