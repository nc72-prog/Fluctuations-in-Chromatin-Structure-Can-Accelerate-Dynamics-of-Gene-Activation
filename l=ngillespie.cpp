
#include <iostream>
#include <random>
#include <cmath>
#include <fstream>
#include <vector>

int main() {

    int L    = 1000;            // DNA length
    int m    = L / 2;           // target position
    int n    = 10;              // number of sites that open at once
    double u     = 1e3;         // sliding rate
    double koff  = 10;          // TF dissociation rate
    double kon   = 1e2;         // TF binding base rate
    int num_trials = 50;       // trials per (alpha, beta) grid point

    int nAlpha = 15;
    int nBeta  = 15;
    double alphaMin = 1e-2, alphaMax = 1e3;
    double betaMin  = 1e-2, betaMax  = 1e3;

    std::vector<double> alphas(nAlpha), betas(nBeta);

    for (int i = 0; i < nAlpha; i++) {
        double logVal = log10(alphaMin)
            + i * (log10(alphaMax) - log10(alphaMin)) / (nAlpha - 1);
        alphas[i] = pow(10.0, logVal);
    }

    for (int i = 0; i < nBeta; i++) {
        double logVal = log10(betaMin)
            + i * (log10(betaMax) - log10(betaMin)) / (nBeta - 1);
        betas[i] = pow(10.0, logVal);
    }

    std::ofstream dataFile("gillespie_heatmap.dat");
    dataFile << "# beta\talpha\tT_sim\n";

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    for (double beta : betas) {
        for (double alpha : alphas) {

            double T_total = 0.0;

            for (int trial = 0; trial < num_trials; trial++) {

                bool dna_open = false;
                int open_start = 0, open_end = 0;
                bool tf_bound = false;
                int position = 0;
                double t = 0.0;

                while (true) {

                    double r_total = 0.0;
                    double alpha_open = alpha / n;
                    double kon_site = kon / n;

                    if (!dna_open) {
                        r_total = alpha_open;
                    } else if (!tf_bound) {
                        r_total = n * kon_site + beta;
                    } else {
                        r_total = koff + beta;
                        if (position > open_start) r_total += u;
                        if (position < open_end)   r_total += u;
                    }

                    double r1;
                    do {
                        r1 = dis(gen);
                    } while (r1 == 0.0);

                    double r2 = dis(gen);
                    t += -log(r1) / r_total;

                    if (!dna_open) {

                        int max_start = L - n + 1;
                        open_start = (int)(dis(gen) * max_start) + 1;
                        open_end = open_start + n - 1;

                        dna_open = true;
                        tf_bound = false;

                    } else if (!tf_bound) {

                        if (r2 < (n * kon_site) / r_total) {
                            tf_bound = true;
                            position = open_start + (int)(dis(gen) * n);

                            if (position == m) break;

                        } else {
                            dna_open = false;
                        }

                    } else {

                        double r_left  = (position > open_start) ? u : 0.0;
                        double r_right = (position < open_end) ? u : 0.0;

                        double cum = r_left / r_total;

                        if (r2 < cum) {

                            position--;
                            if (position == m) break;

                        } else {

                            cum += r_right / r_total;

                            if (r2 < cum) {

                                position++;
                                if (position == m) break;

                            } else {

                                cum += koff / r_total;

                                if (r2 < cum) {
                                    tf_bound = false;
                                } else {
                                    dna_open = false;
                                    tf_bound = false;
                                    position = 0;
                                }
                            }
                        }
                    }
                }

                T_total += t;
            }

            double T_sim = T_total / num_trials;
            dataFile << beta << "\t" << alpha << "\t" << T_sim << "\n";
        }

        std::cout << "beta=" << beta << " done\n";
    }

    dataFile.close();
    std::cout << "Done. Data written to gillespie_heatmap.dat\n";

    return 0;
}
