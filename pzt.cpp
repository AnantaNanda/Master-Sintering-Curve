#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "matplotlibcpp.h" // Ensure this file is in your project directory
namespace plt = matplotlibcpp;
using namespace std;
// Function to load data from a text file
void loadData(const string &filename, vector<double> &t, vector<double> &v, vector<double> &i) {
ifstream infile(filename);
double time, voltage, current;
while (infile >> time >> voltage >> current) {
t.push_back(time);
v.push_back(voltage);
i.push_back(current);
}
infile.close();
}
int main() {
vector<double> t, v, i;
// Load data from the PZT.txt file
loadData("PZT.txt", t, v, i);
// Power calculation
vector<double> power(t.size());
for (size_t j = 0; j < t.size(); ++j) {
power[j] = v[j] * i[j]; // Power in watts
}
// Material and geometric properties
double density = 8000; // Density of sample in kg/m^3
double dia = 6e-3; // Diameter of sample in meters
double h = 5.5e-3; // Height of sample in meters
double area = 3.14 * pow(dia / 2.0, 2); // Cross-sectional area of sample
// Current density, electric field, and conductivity
vector<double> cd(t.size()), E(t.size()), sigma(t.size()), ln_sigma(t.size());
for (size_t j = 0; j < t.size(); ++j) {
cd[j] = i[j] / area; // Current density
E[j] = v[j] / h; // Electric field
sigma[j] = cd[j] / E[j]; // Conductivity
ln_sigma[j] = log(sigma[j]); // Log of conductivity
}
// Volume and power per unit mass
double vs = 3.14 * pow(dia / 2.0, 2) * h; // Volume of the sample
vector<double> power_m(t.size());
for (size_t j = 0; j < t.size(); ++j) {
power_m[j] = power[j] / (vs * density); // Power per unit mass
}
// Inconel properties and Joule heating calculations
double R_inconel = 1.31e-6; // Resistivity of inconel in ohm-m
double density_inconel = 8113; // Density of inconel in kg/m^3
double h_inconel = 3.9e-3; // Height of inconel (bottom)
double dia_inconel = 6.05e-3; // Diameter of inconel (bottom)
vector<double> joule_inconel(t.size()), joule_inconel_m(t.size());
for (size_t j = 0; j < t.size(); ++j) {
joule_inconel[j] = i[j] * i[j] * R_inconel * h_inconel / (3.14 * pow(dia_inconel / 2.0, 2));
joule_inconel_m[j] = joule_inconel[j] / (h_inconel * 3.14 * pow(dia_inconel / 2.0, 2) * density_inconel);
}
// Plotting the data using matplotlibcpp
plt::figure_size(1200, 780); // Set figure size
// Plot power per unit mass for the sample (3YSZ)
plt::named_plot("Power per unit mass (3YSZ)", t, power_m);
// Plot joule heating per unit mass in inconel
plt::named_plot("Joule heating in Inconel", t, joule_inconel_m);
// Add labels and title
plt::xlabel("Time (seconds)");
plt::ylabel("Power per unit mass (W/kg)");
plt::title("Power per Unit Mass vs Time");
plt::legend(); // Add legend
plt::grid(true); // Add grid to the plot
// Save the plot as a high-resolution image
plt::save("FS_P_vs_time.png");
// Show the plot
plt::show();
cout << "Plotting complete." << endl;
return 0;
}