#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>

using namespace std;

const double epsilon = 1.0; //ok
const double delta = 0.1; //ok
const int nx = 150, ny = 100; //ok
const double xmax = delta * nx, ymax = delta * ny; //ok
const double V1 = 10.0, V2 = 0.0; //ok
const double TOL = 1e-8; //ok

// sprawdzone, dobrze
double rho1(double x, double y) {
    double sigma_x = 0.1 * xmax;
    double sigma_y = 0.1 * ymax;
    return exp(-pow((x - 0.35 * xmax) / sigma_x, 2) - pow((y - 0.5 * ymax) / sigma_y, 2));
}

// sprawdzone, dobrze
double rho2(double x, double y) {
    double sigma_x = 0.1 * xmax;
    double sigma_y = 0.1 * ymax;
    return -exp(-pow((x - 0.65 * xmax) / sigma_x, 2) - pow((y - 0.5 * ymax) / sigma_y, 2));
}

double warunek_stopu(const vector<vector<double>> &V, const vector<vector<double>> &rho) {
    double S = 0.0;
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            S += pow(delta, 2) * ((0.5 * pow((V[i+1][j] - V[i][j]) / delta, 2)) + (0.5 * pow((V[i][j+1] - V[i][j]) / delta, 2)) - rho[i][j] * V[i][j]);
        }
    }
    return S;
}

void relaksacja_globalna(const vector<vector<double>> &rho, double omega, vector<vector<double>> &V, vector<double> &S_values) {
    vector<vector<double>> V_new(nx + 1, vector<double>(ny + 1, 0.0));

    for (int i = 0; i <= nx; i++) {
        V[i][0] = V1;
        V[i][ny] = V2;
        V_new[i][0] = V1;
        V_new[i][ny] = V2;
    }

    double S_prev = 1e10;
    int count = 0;
    while (true) {
        
        // Punkty wewnetrzne
        for (int i = 1; i < nx; ++i) {
            for (int j = 1; j < ny; ++j) {
                V_new[i][j] =  1.0 / 4.0 * (V[i + 1][j] + V[i - 1][j] + V[i][j + 1] + V[i][j - 1] +
                                             delta * delta / epsilon * rho[i][j]);
            }
        }

        // Neumann
        for (int j = 1; j < ny; ++j) {
            V_new[0][j] = V_new[1][j];
            V_new[nx][j] = V_new[nx - 1][j];
        }

        // Mieszanie rozwiazan
        // Potencjalnie zle zakresy petli
        for (int i = 0; i < nx + 1; ++i) {
            for (int j = 0; j < ny + 1; ++j) {
                V[i][j] = (1 - omega) * V[i][j] + omega * V_new[i][j];
            }
        }


        double S = warunek_stopu(V_new, rho);
        S_values.push_back(S);

        // sprawdzanie warunku stopu
        if (fabs(S - S_prev) / fabs(S_prev) < TOL) {
            cout << "S = " << S << ", S_prev = " << S_prev << endl;
            cout << "Warunek stopu spelniony" << endl;
            cout << "Count = " << count << endl;
            break;
        }
        S_prev = S;
        // V = V_new;
        count++;
    }
}

void relaksacja_lokalna(const vector<vector<double>> &rho, double omega, vector<vector<double>> &V, vector<double> &S_values) {

    for (int i = 0; i <= nx; i++) {
        V[i][0] = V1;
        V[i][ny] = V2;
    }

    double S_prev = 1e10;
    int count = 0;
    while (true) {
        // Wewnetrzne punkty
        for (int i = 1; i < nx; ++i) {
            for (int j = 1; j < ny; ++j) {
                V[i][j] = (1 - omega) * V[i][j] +
                          omega / 4.0 * (V[i + 1][j] + V[i - 1][j] + V[i][j + 1] + V[i][j - 1] +
                                         delta * delta / epsilon * rho[i][j]);
            }
        }

        // Neumann
        for (int j = 1; j < ny; ++j) {
            V[0][j] = V[1][j];
            V[nx][j] = V[nx - 1][j];
        }


        double S = warunek_stopu(V, rho);
        S_values.push_back(S);

        // Sprawdzanie warunku stopu
        if (abs(S - S_prev) / S_prev < TOL) {
            cout << "S = " << S << ", S_prev = " << S_prev << endl;
            cout << "Warunek stopu spelniony" << endl;
            cout << "Count = " << count << endl;
            break;
        }
        S_prev = S;
        count++;
    }
}

int main() {
    vector<vector<double>> rho(nx + 1, vector<double>(ny + 1, 0.0));
    for (int i = 0; i <= nx; ++i) {
        for (int j = 0; j <= ny; ++j) {
            double x = i * delta;
            double y = j * delta;
            rho[i][j] = rho1(x, y) + rho2(x, y);
        }
    }

    vector<vector<double>> V_global_06(nx + 1, vector<double>(ny + 1, 0.0));
    vector<vector<double>> V_global_10(nx + 1, vector<double>(ny + 1, 0.0));
    vector<vector<double>> V_local_10(nx + 1, vector<double>(ny + 1, 0.0));
    vector<vector<double>> V_local_14(nx + 1, vector<double>(ny + 1, 0.0));
    vector<vector<double>> V_local_18(nx + 1, vector<double>(ny + 1, 0.0));
    vector<vector<double>> V_local_19(nx + 1, vector<double>(ny + 1, 0.0));
    vector<double> S_global_06, S_global_10, S_local_10, S_local_14, S_local_18, S_local_19;

    relaksacja_globalna(rho, 0.6, V_global_06, S_global_06);
    relaksacja_globalna(rho, 1.0, V_global_10, S_global_10);

    relaksacja_lokalna(rho, 1.0, V_local_10, S_local_10);
    relaksacja_lokalna(rho, 1.4, V_local_14, S_local_14);
    relaksacja_lokalna(rho, 1.8, V_local_18, S_local_18);
    relaksacja_lokalna(rho, 1.9, V_local_19, S_local_19);

    
    ofstream potential_file_06("potential_V_06.csv");
    for (int i = 0; i <= nx; ++i) {
        for (int j = 0; j <= ny; ++j) {
            potential_file_06 << V_global_06[i][j];
            if (j < ny) potential_file_06 << ",";
        }
        potential_file_06 << "\n";
    }
    potential_file_06.close();

    ofstream potential_file_10("potential_V_10.csv");
    for (int i = 0; i <= nx; ++i) {
        for (int j = 0; j <= ny; ++j) {
            potential_file_10 << V_global_10[i][j];
            if (j < ny) potential_file_10 << ",";
        }
        potential_file_10 << "\n";
    }
    potential_file_10.close();

    ofstream theta_file_06("theta_06.csv");
    for (int i = 1; i < nx; ++i) {
        for (int j = 1; j < ny; ++j) {
            double gradient_V = (V_global_06[i + 1][j] - 2 * V_global_06[i][j] + V_global_06[i - 1][j]) / pow(delta, 2) +
                                (V_global_06[i][j + 1] - 2 * V_global_06[i][j] + V_global_06[i][j - 1]) / pow(delta, 2);
            double theta = gradient_V + rho[i][j] / epsilon;
            theta_file_06 << theta;
            if (j < ny - 1) theta_file_06 << ",";
        }
        theta_file_06 << "\n";
    }
    theta_file_06.close();

    ofstream theta_file_10("theta_10.csv");
    for (int i = 1; i < nx; ++i) {
        for (int j = 1; j < ny; ++j) {
            double gradient_V = (V_global_10[i + 1][j] - 2 * V_global_10[i][j] + V_global_10[i - 1][j]) / pow(delta, 2) +
                                (V_global_10[i][j + 1] - 2 * V_global_10[i][j] + V_global_10[i][j - 1]) / pow(delta, 2);

            double theta = gradient_V + rho[i][j] / epsilon;
            theta_file_10 << theta;
            if (j < ny - 1) theta_file_10 << ",";
        }
        theta_file_10 << "\n";
    }
    theta_file_10.close();

    ofstream outfile_global_06("global_relaxation_06.csv");
    outfile_global_06 << setprecision(10);
    outfile_global_06 << "Iteration,S_global_06\n";
    for (size_t i = 0; i < S_global_06.size(); ++i) {
        outfile_global_06 << i << ",";
        outfile_global_06 << (i < S_global_06.size() ? S_global_06[i] : 0) << "\n";
    }
    outfile_global_06.close();

    ofstream outfile_global_10("global_relaxation_10.csv");
    outfile_global_10 << setprecision(10);
    outfile_global_10 << "Iteration,S_global_10\n";
    for (size_t i = 0; i < S_global_10.size(); ++i) {
        outfile_global_10 << i << ",";
        outfile_global_10 << (i < S_global_10.size() ? S_global_10[i] : 0) << "\n";
    }
    outfile_global_10.close();

    ofstream outfile_local_10("local_relaxation_10.csv");
    outfile_local_10 << setprecision(10);
    outfile_local_10 << "Iteration,S_local_10\n";
    for (size_t i = 0; i < S_local_10.size(); ++i) {
        outfile_local_10 << i << ",";
        outfile_local_10 << (i < S_local_10.size() ? S_local_10[i] : 0) << "\n";
    }
    outfile_local_10.close();

    ofstream outfile_local_14("local_relaxation_14.csv");
    outfile_local_14 << setprecision(10);
    outfile_local_14 << "Iteration,S_local_14\n";
    for (size_t i = 0; i < S_local_14.size(); ++i) {
        outfile_local_14 << i << ",";
        outfile_local_14 << (i < S_local_14.size() ? S_local_14[i] : 0) << "\n";
    }
    outfile_local_14.close();

    ofstream outfile_local_18("local_relaxation_18.csv");
    outfile_local_18 << setprecision(10);
    outfile_local_18 << "Iteration,S_local_18\n";
    for (size_t i = 0; i < S_local_18.size(); ++i) {
        outfile_local_18 << i << ",";
        outfile_local_18 << (i < S_local_18.size() ? S_local_18[i] : 0) << "\n";
    }
    outfile_local_18.close();

    ofstream outfile_local_19("local_relaxation_19.csv");
    outfile_local_19 << setprecision(10);
    outfile_local_19 << "Iteration,S_local_19\n";
    for (size_t i = 0; i < S_local_19.size(); ++i) {
        outfile_local_19 << i << ",";
        outfile_local_19 << (i < S_local_19.size() ? S_local_19[i] : 0) << "\n";
    }
    outfile_local_19.close();

    return 0;
}
