#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <cmath>
#include <matplotlibcpp.h>
namespace plt = matplotlibcpp;
using namespace std;
// Constants
const double R = 8.314; // Gas constant in J/(mol K)
const double learning_rate = 1e-6; // Learning rate for optimization
// Function to calculate work of sintering Q
vector<double> work_of_sintering(const vector<double>& t, double T, double Q) {
vector<double> integrand(t.size());
vector<double> result(t.size(), 0.0);
for (size_t i = 0; i < t.size(); ++i) {
integrand[i] = exp(-Q / (R * T)) / T;
}
// Numerical integration using cumulative trapezoidal rule
for (size_t i = 1; i < t.size(); ++i) {
result[i] = result[i - 1] + 0.5 * (integrand[i] + integrand[i - 1]) * (t[i] - t[i - 1]);
}
return result;
}
// Sigmoid function to model the densification curve
double sigmoid(double Theta, double a, double b, double rho_i, double rho_f) {
// Add a small value to avoid log(0) issues
double Theta_safe = Theta == 0 ? 1e-8 : Theta;
return rho_f + (rho_i - rho_f) / (1 + exp((log(Theta_safe) - a) / b));
}
// Mean squared error loss function
double mse_loss(const vector<double>& predicted, const vector<double>& actual) {
double error = 0.0;
for (size_t i = 0; i < predicted.size(); ++i) {
error += pow(predicted[i] - actual[i], 2);
}
return error / predicted.size();
}
// Gradient descent to optimize parameters Q, a, and b
void gradient_descent(const vector<double>& t, const vector<double>& T, const
vector<double>& rho_exp, vector<double>& params, int epochs) {
double Q = params[0], a = params[1], b = params[2];
for (int epoch = 0; epoch < epochs; ++epoch) {
vector<double> Theta_fitted = work_of_sintering(t, T[0], Q); // Using first temperature as
example
vector<double> rho_fitted(t.size());
for (size_t i = 0; i < t.size(); ++i) {
rho_fitted[i] = sigmoid(Theta_fitted[i], a, b, 0.5, 0.9);
}
// Calculate gradients for Q, a, b using finite differences
double dQ = 0.0, da = 0.0, db = 0.0;
for (size_t i = 0; i < t.size(); ++i) {
double error = rho_fitted[i] - rho_exp[i];
dQ += error * (-Theta_fitted[i] * exp(-Q / (R * T[0])) / (R * T[0] * T[0]));
da += error * (rho_fitted[i] - 0.9) * log(Theta_fitted[i]) / (b * pow(1 +
exp((log(Theta_fitted[i]) - a) / b), 2));
db += error * (rho_fitted[i] - 0.9) * (log(Theta_fitted[i]) - a) / (b * b * pow(1 +
exp((log(Theta_fitted[i]) - a) / b), 2));
}
// Update parameters
Q -= learning_rate * dQ;
a -= learning_rate * da;
b -= learning_rate * db;
// Print loss for each epoch (optional for tracking)
if (epoch % 100 == 0) {
double loss = mse_loss(rho_fitted, rho_exp);
cout << "Epoch " << epoch << ", Loss: " << loss << endl;
}
}
params[0] = Q;
params[1] = a;
params[2] = b;
}
int main() {
// Experimental data (user input or loaded from files)
vector<double> t_data = {0, 100, 200, 300, 400, 500}; // Time in seconds
vector<double> T_data = {1200, 1250, 1300, 1350, 1400, 1450}; // Temperature in Kelvin
vector<double> rho_data = {0.5, 0.6, 0.7, 0.8, 0.85, 0.91}; // Relative Density
// Initial guesses for parameters [Q, a, b]
vector<double> p0 = {40000, 0.1, 0.1};
// Perform gradient descent to fit the curve and get optimal parameters
int epochs = 10000; // Number of optimization iterations
gradient_descent(t_data, T_data, rho_data, p0, epochs);
// Calculate fitted relative density using the optimized parameters
vector<double> Theta_fitted = work_of_sintering(t_data, T_data[0], p0[0]); // Using first
temperature as an example
vector<double> rho_fitted(t_data.size());
for (size_t i = 0; i < Theta_fitted.size(); ++i) {
rho_fitted[i] = sigmoid(Theta_fitted[i], p0[1], p0[2], 0.5, 0.9);
}
// Plot the experimental vs fitted data
plt::plot(t_data, rho_data, "o", {{"label", "Experimental Data"}});
plt::plot(t_data, rho_fitted, "-o", {{"label", "Fitted MSC"}, {"color", "blue"}});
plt::xlabel("Time (seconds)");
plt::ylabel("Relative Density");
plt::legend();
plt::title("Master Sintering Curve (MSC)");
plt::grid(true);
plt::show();
// For debugging: Plot work of sintering (Theta) over time
plt::plot(t_data, Theta_fitted, "o-", {{"label", "Work of Sintering (Theta)"}, {"color", "green"}});
plt::xlabel("Time (seconds)");
plt::ylabel("Work of Sintering (Theta)");
plt::legend();
plt::title("Work of Sintering vs Time");
plt::grid(true);
plt::show();
// Output the optimized parameters
cout << "Optimized Activation Energy (Q): " << p0[0] / 1000 << " kJ/mol" << endl;
cout << "Optimized a: " << p0[1] << endl;
cout << "Optimized b: " << p0[2] << endl;
return 0;
}