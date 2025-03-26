#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <string>

// Constants
const double DELTA = 0.01;
const double RHO = 1.0;
const double MU = 1.0;
const int NX = 200, NY = 90;
const int J1 = 55, I1 = 50;
const int IT_MAX = 20000;

double q_wy(double Qwe) {
    return Qwe * (std::pow(NY * DELTA, 3) - std::pow(J1 * DELTA, 3) - 3 * J1 * DELTA * std::pow(NY * DELTA, 2) + 3 * NY * DELTA * std::pow(J1 * DELTA, 2)) / std::pow(NY * DELTA, 3);
}


void inicjalizacja(std::vector<std::vector<double>> &vec) {
    for (int i = 0; i <= NX; i++) {
        for (int j = 0; j <= NY; j++) {
            vec[i][j] = 0.0;
        }
    }
}


// WB dla psi
void wb_psi(std::vector<std::vector<double>> &psi, double Qwe, double Qwy) {
    double y;

    // Brzeg A
    for (int j = J1; j <= NY; j++) {
        y = j * DELTA;
        psi[0][j] = (Qwe / (2.0 * MU)) * (std::pow(y, 3) / 3.0 - std::pow(y, 2) / 2.0 * (J1 + NY) * DELTA + y * J1 * NY * DELTA * DELTA);
    }

    // Brzeg C
    for (int j = 0; j <= NY; j++) {
        y = j * DELTA;
        psi[NX][j] = (Qwy / (2 * MU)) * (std::pow(y, 3) / 3.0 - std::pow(y, 2) / 2.0 * NY * DELTA) +
                     (Qwe * std::pow(J1 * DELTA, 2) * (-J1 * DELTA + 3 * NY * DELTA)) / (12.0 * MU);
    }

    // Brzeg B
    for (int i = 1; i < NX; i++) {
        psi[i][NY] = psi[0][NY];
    }

    // Brzeg D
    for (int i = I1; i < NX; i++) {
        psi[i][0] = psi[0][J1];
    }

    // Brzeg E
    for (int j = 1; j <= J1; j++) {
        psi[I1][j] = psi[0][J1];
    }

    // Brzeg F
    for (int i = 1; i <= I1; i++) {
        psi[i][J1] = psi[0][J1];
    }
}

// WB dla zeta
void wb_zeta(std::vector<std::vector<double>> &zeta, std::vector<std::vector<double>> &psi, double Qwe, double Qwy) {
    double DELTA_2_INV = 2.0 / (DELTA * DELTA);
    double y;

    // Brzeg A
    for (int j = J1; j <= NY; j++) {
        y = j * DELTA;
        zeta[0][j] = (Qwe / (2.0 * MU)) * (2.0 * y - J1 * DELTA - NY * DELTA);
    }

    // Brzeg C
    for (int j = 0; j <= NY; j++) {
        y = j * DELTA;
        zeta[NX][j] = (Qwy / (2.0 * MU)) * (2.0 * y - NY * DELTA);
    }

    // Brzeg B
    for (int i = 1; i < NX; i++) {
        zeta[i][NY] = DELTA_2_INV * (psi[i][NY - 1] - psi[i][NY]);
    }

    // Brzeg D
    for (int i = I1 + 1; i < NX; i++) {
        zeta[i][0] = DELTA_2_INV * (psi[i][1] - psi[i][0]);
    }

    // Brzeg E
    for (int j = 1; j < J1; j++) {
        zeta[I1][j] = DELTA_2_INV * (psi[I1 + 1][j] - psi[I1][j]);
    }

    // Brzeg F
    for (int i = 1; i <= I1; i++) {
        zeta[i][J1] = DELTA_2_INV * (psi[i][J1 + 1] - psi[i][J1]);
    }

    // Wierzchołek E/F
    zeta[I1][J1] = 0.5 * (zeta[I1 - 1][J1] + zeta[I1][J1 - 1]);
}

double compute_gamma(const std::vector<std::vector<double>>& psi, const std::vector<std::vector<double>>& zeta) {
    double result = 0.0;
    for(int i = 1; i < NX; ++i) 
    {
        result += psi[i + 1][J1 + 2] + psi[i - 1][J1 + 2] + psi[i][J1 + 3] + psi[i][J1 + 1] - 4.0 * psi[i][J1 + 2] - DELTA * DELTA * zeta[i][J1 + 2];
    }
    return result;
}

// Algorytm relaksacji
void algorytm_relaksacji(std::vector<std::vector<double>> &psi, std::vector<std::vector<double>> &zeta, double Qwe, double Qwy) {
    wb_psi(psi, Qwe, Qwy);
    double omega;

    for (int it = 1; it <= IT_MAX; it++) {
        omega = (it >= 2000) ? 1.0 : 0.0;

        for(int i = 1; i < NX; ++i)
        {
            for(int j = 1; j < NY; ++j) 
            {
                if (!(i == 0 || i == NX || j == 0 || j == NY || (i == I1 && j <= J1) || (j == J1 && i <= I1)))
                {
                    psi[i][j] = 0.25 * (psi[i + 1][j] + psi[i - 1][j] + psi[i][j + 1] + psi[i][j - 1] - DELTA * DELTA * zeta[i][j]);
                    zeta[i][j] = 0.25 * (zeta[i + 1][j] + zeta[i - 1][j] + zeta[i][j + 1] + zeta[i][j - 1]) - omega * (RHO / (16.0 * MU)) * ((psi[i][j + 1] - psi[i][j - 1]) * (zeta[i + 1][j] - zeta[i - 1][j]) - (psi[i + 1][j] - psi[i - 1][j]) * (zeta[i][j + 1] - zeta[i][j - 1]));
                }
            }
        }
        wb_zeta(zeta, psi, Qwe, Qwy);

        // Check error
        double gamma = compute_gamma(psi, zeta);

        std::cout << "Iteration " << it << ", Gamma = " << gamma << std::endl;
    }
}


void save_to_file(std::vector<std::vector<double>> &psi, std::vector<std::vector<double>> &zeta, std::string title) {
    std::ofstream file(title);
    for (int i = 1; i < NX; i++) {
        for (int j = 1; j < NY; j++) {
            file << DELTA * i << " " << DELTA * j << " " << psi[i][j] << " " << zeta[i][j] << "\n";
        }
    }
    file.close();
}

void predkosc(std::vector<std::vector<double>> &u, std::vector<std::vector<double>> &v, std::vector<std::vector<double>> &psi) {
    for (int i = 1; i < NX; i++) {
        for (int j = 1; j < NY; j++) {
            // u(x, y) = dpsi / dy
            u[i][j] = (psi[i][j + 1] - psi[i][j - 1]) / (2.0 * DELTA);
            // v(x, y) = -dpsi / dy
            v[i][j] = -(psi[i + 1][j] - psi[i - 1][j]) / (2.0 * DELTA);
        }
    }
}

int main() {
    std::vector<double> Q_values = {-1000.0, -4000.0, 4000.0};
    std::vector<std::vector<double>> psi(NX + 1, std::vector<double>(NY + 1, 0.0));
    std::vector<std::vector<double>> zeta(NX + 1, std::vector<double>(NY + 1, 0.0));
    std::vector<std::vector<double>> u(NX + 1, std::vector<double>(NY + 1, 0.0));
    std::vector<std::vector<double>> v(NX + 1, std::vector<double>(NY + 1, 0.0));

    for (double Qwe : Q_values) {
        inicjalizacja(psi);
        inicjalizacja(zeta);
        inicjalizacja(u);
        inicjalizacja(v);
        double Qwy = q_wy(Qwe);
        algorytm_relaksacji(psi, zeta, Qwe, Qwy);
        save_to_file(psi, zeta, "psi_zeta_Q" + std::to_string(static_cast<int>(Qwe)) + ".txt");
        predkosc(u, v, psi);
        save_to_file(u, v, "u_v_Q" + std::to_string(static_cast<int>(Qwe)) + ".txt");
    }

    return 0;
}
