#include <iostream>
#include <fstream>
#include <cmath>
#include <array>
#include <sys/stat.h>  
#include <sys/types.h> 
#include <direct.h>   
#include <sstream>  // For std::ostringstream
#include <iomanip>  // For setw, setfill
#include <string>   // For string operations
using namespace std;

const int nx = 64;        // Grid size
const int ny = 64;
const int npop = 9;           // D2Q9 directions


const int nsteps = 20;
const int noutput = 1;

// LBM velocities and weights
const int cx[npop] = { 0, 1, 0, -1, 0, 1, -1, -1, 1 };
const int cy[npop] = { 0, 0, 1, 0, -1, 1, 1, -1, -1 };
const double weights[npop] = { 4.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0 };

// Simulation parameters 
double tau_water = 1.0;    // Relaxation time (water)
double tau_air = 1.0;      // Relaxation time (air)
double rho_water0 = 1.0;    // Water reference density
double rho_air0 = 0.1;  // Air reference density
int radius = 10;
double G = -0.08;
const double cs2 = 1.0 / 3.0; // Speed of sound squared


//External force
double g = 0.01; 


// Field arrays
double rho_water[nx * ny];
double rho_air[nx * ny];

double u1_water[nx * ny];
double u2_water[nx * ny];

double u1_air[nx * ny];
double u2_air[nx * ny];

double f_water_mem[nx * ny * npop];
double f_air_mem[nx * ny * npop];

double feq_water[nx * ny * npop];
double feq_air[nx * ny * npop];

//Mass conservation
double initial_total_mass_water = 0.0;
double initial_total_mass_air = 0.0;

void saveToBinary(
    const std::string& sim_dir,
    const std::string& id,
    const std::string& var_name,
    int step,
    double* data,
    size_t size
) {
    std::ostringstream filename;
    filename << sim_dir << "/" << id << "_" << var_name
             << std::setw(6) << std::setfill('0') << step << ".bin";

    std::ofstream file(filename.str(), std::ios::binary);
    if (!file) {
        std::cerr << "Erro ao abrir " << filename.str() << " para escrita." << std::endl;
        return;
    }

    file.write(reinterpret_cast<const char*>(data), size * sizeof(double));
    file.close();
}

void generateSimulationInfoFile(
    const std::string& SIM_DIR,
    const std::string& SIM_ID,
    const std::string& VELOCITY_SET,
    const int NX,
    const int NY,
    const int NSTEPS,
    const int MACRO_SAVE,
    const double TAU,
    const double MLUPS
) {
    std::string INFO_FILE = SIM_DIR + "/" + SIM_ID + "_info.txt";
    try {
        std::ofstream file(INFO_FILE);

        if (!file.is_open()) {
            std::cerr << "Erro ao abrir o arquivo: " << INFO_FILE << std::endl;
            return;
        }

        file << "---------------------------- SIMULATION INFORMATION ----------------------------\n"
             << "                           Simulation ID: " << SIM_ID << '\n'
             << "                           Velocity set: " << VELOCITY_SET << '\n'
             << "                           Precision: double\n"
             << "                           NX: " << NX << '\n'
             << "                           NY: " << NY << '\n'
             << "                           Tau: " << TAU << '\n'
             << "                           Save steps: " << MACRO_SAVE << '\n'
             << "                           Nsteps: " << NSTEPS << '\n'
             << "                           MLUPS: " << MLUPS << '\n'
             << "--------------------------------------------------------------------------------\n";

        file.close();
        std::cout << "Arquivo de informações gerado em: " << INFO_FILE << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Erro ao gerar arquivo de informações: " << e.what() << std::endl;
    }
}

//#############################################################


int main(int argc, char** argv) 
{
    if (argc == 3) {
        tau_water = atof(argv[1]);
        radius = atoi(argv[2]);
    }

    // Pointers for distributions
    double* f_water = f_water_mem;          // Water distribution
    double* f_air = f_air_mem;  // Air distribution
    double* f_water_tmp = new double[nx * ny * npop] ();    // Temporary for water
    double* f_air_tmp = new double[nx * ny * npop] ();// Temporary for air

    // Initialization of fields - Bubble setup
    for (int idx = 0; idx < nx * ny; ++idx)
    {
        int ix = idx % nx;
        int iy = idx / nx;
        double dx = ix - nx / 2.0;
        double dy = iy - ny / 2.0;
        double r = sqrt(dx * dx + dy * dy);

        // Initialize densities
        if (r <= radius) 
        {
            // Inside bubble: air dominant
            rho_water[idx] = 0.00001 * rho_water0;
            rho_air[idx] = rho_air0;
        }
        else 
        {
            // Outside bubble: water dominant
            rho_water[idx] = rho_water0;
            rho_air[idx] = 0.00001 * rho_air0;
        }

        u1_water[idx] = 0.0;
        u2_water[idx] = 0.0;
        u1_air[idx] = 0.0;
        u2_air[idx] = 0.0;

        // Initialize distributions to equilibrium
        for (int k = 0; k < npop; ++k) 
        {
            f_water[9 * idx + k] = weights[k] * rho_water[idx];
            f_air[9 * idx + k] = weights[k] * rho_air[idx];
        }
    }


    //Initial mass
    for (int idx = 0; idx < nx * ny; ++idx)
    {
        initial_total_mass_water += rho_water[idx];
        initial_total_mass_air += rho_air[idx];
    }



    std::time_t start = std::time(nullptr);
    // LBM LOOP
    for (int t = 0; t < nsteps; t++)
    {

        //Compute macroscopic densities
        for (int idx = 0; idx < nx * ny; idx++)
        {
            rho_water[idx] = 0.0;
            rho_air[idx] = 0.0;
            for (int k = 0; k < npop; k++)
            {
                rho_water[idx] += f_water[9 * idx + k];
                rho_air[idx] += f_air[9 * idx + k];

                // Clamp densities to avoid numerical instability
                //rho_water[idx] = fmax(rho_water[idx], 1e-10);
                //rho_air[idx] = fmax(rho_air[idx], 1e-10);
            }
        }

        //Mass checking 
        double current_total_mass_water = 0.0; //Reser masses 
        double current_total_mass_air = 0.0;


        for (int idx = 0; idx < nx * ny; idx++) 
        {
            current_total_mass_water += rho_water[idx];
            current_total_mass_air += rho_air[idx];
        }

        double mass_error_water = fabs(current_total_mass_water - initial_total_mass_water) / initial_total_mass_water;
        double mass_error_air = fabs(current_total_mass_air - initial_total_mass_air) / initial_total_mass_air;

        if (t % 1 == 0)
        {
            cout << "Step " << t
                << ": Water mass error = " << mass_error_water * 100 << "%, "
                << "Air mass error = " << mass_error_air * 100 << "%" << endl;
        }




        //Collision and streaming 
        for (int iY = 0; iY < ny; iY++)
        {
            for (int iX = 0; iX < nx; iX++)
            {
                int idx = iY * nx + iX; //Current cell index

                //SHAN-CHEN MODEL
                double Fx_water = 0.0, Fy_water = 0.0;
                double Fx_air = 0.0, Fy_air = 0.0;

                double fx_water = 0.0, fy_water = 0.0;
                double fx_air = 0.0, fy_air = 0.0;

                for (int k = 0; k < npop; k++)
                {
                    int iX2 = (iX + cx[k] + nx) % nx; // Neighbor x-position
                    int iY2 = (iY + cy[k] + ny) % ny; // Neighbor y-position
                    int idx2 = iY2 * nx + iX2; // Neighbor cell index

                    
                    //FIRST WAY
                    // water force from air density
                    fx_water += weights[k] * cx[k] * (1 - std::exp(-rho_air[idx2]));
                    fy_water += weights[k] * cy[k] * (1 - std::exp(-rho_air[idx2]));

                    // air force from water density
                    fx_air += weights[k] * cx[k] * (1 - std::exp(-rho_water[idx2] ));
                    fy_air += weights[k] * cy[k] * (1 - std::exp(-rho_water[idx2] ));


                    //SECOND WAY 
                    //Fx_water += weights[k] * cx[k] * rho_air[idx2];
                    //Fy_water += weights[k] * cy[k] * rho_air[idx2];

                    // Air force from water density
                    //Fx_air += weights[k] * cx[k] * rho_water[idx2];
                    //Fy_air += weights[k] * cy[k] * rho_water[idx2];



                }

                //FIRST WAY
                Fx_water = -G * (1 - std::exp(-rho_water[idx])) * fx_water;
                Fy_water = -G * (1 - std::exp(-rho_water[idx] )) * fy_water;
                Fx_air = -G * (1 - std::exp(-rho_air[idx] )) * fx_air;
                Fy_air = -G * (1 - std::exp(-rho_air[idx] )) * fy_air;


                //SECOND WAY
                //Fx_water = -G * rho_water[idx] * Fx_water;
                //Fy_water = -G * rho_water[idx] * Fy_water;
                //Fx_air = -G * rho_air[idx] * Fx_air;
                //Fy_air = -G * rho_air[idx] * Fy_air;


                //ADD external force 

                //VELOCITIES CALCULATIONS

                //Water velocity 

           
                u1_water[idx] = (f_water[9 * idx + 1] - f_water[9 * idx + 3] + f_water[9 * idx + 5]
                    - f_water[9 * idx + 6] - f_water[9 * idx + 7] + f_water[9 * idx + 8]) / rho_water[idx];
                u2_water[idx] = (f_water[9 * idx + 2] - f_water[9 * idx + 4] + f_water[9 * idx + 5]
                    + f_water[9 * idx + 6] - f_water[9 * idx + 7] - f_water[9 * idx + 8]) / rho_water[idx];

                // Add force contribution (half-way bounce-back) - Guo's forcing correction
                u1_water[idx] += 0.5 * (Fx_water ) / rho_water[idx];
                u2_water[idx] += 0.5 * (Fy_water ) / rho_water[idx];
      

                // Air velocity
             
                u1_air[idx] = (f_air[9 * idx + 1] - f_air[9 * idx + 3] + f_air[9 * idx + 5] - f_air[9 * idx + 6] - f_air[9 * idx + 7] + f_air[9 * idx + 8]) / rho_air[idx];
                u2_air[idx] = (f_air[9 * idx + 2] - f_air[9 * idx + 4] + f_air[9 * idx + 5] + f_air[9 * idx + 6] - f_air[9 * idx + 7] - f_air[9 * idx + 8]) / rho_air[idx];


                // Add force contribution
                u1_air[idx] += 0.5 * (Fx_air ) / rho_air[idx];
                u2_air[idx] += 0.5 * (Fy_air ) / rho_air[idx];
         



                //GUO FORCING SCHEME
                double u1w = u1_water[idx];
                double u2w = u2_water[idx];
                double u1a = u1_air[idx];
                double u2a = u2_air[idx];

                double fguo_water[npop];
                double fguo_air[npop];

                for (int k = 0; k < npop; ++k)
                {
                    fguo_water[k] = weights[k] * (1.0 - 0.5 / tau_water) * (
                        (3.0 * (cx[k] - u1w) + 9.0 * cx[k] * (cx[k] * u1w + cy[k] * u2w)) * Fx_water
                        + (3.0 * (cy[k] - u2w) + 9.0 * cy[k] * (cx[k] * u1w + cy[k] * u2w)) * Fy_water
                        );

                    fguo_air[k] = weights[k] * (1.0 - 0.5 / tau_air) * (
                        (3.0 * (cx[k] - u1a) + 9.0 * cx[k] * (cx[k] * u1a + cy[k] * u2a)) * Fx_air
                        + (3.0 * (cy[k] - u2a) + 9.0 * cy[k] * (cx[k] * u1a + cy[k] * u2a)) * Fy_air
                        );

                }

                //FORCE CHECKING
                //double sumF_water = 0.0;
                //double sumF_air = 0.0;

                //for (int k = 0; k < npop; ++k) 
                //{
                //    sumF_water += fguo_water[k];
                //    sumF_air += fguo_air[k];
                //}

                //std::cout << "Check cell " << idx
                //    << ": Sum fguo_water = " << sumF_water
                //    << ", Sum fguo_air = " << sumF_air << std::endl;




                //PURE BGK - "JUST FOR CHECKING"
                //double rhow = rho_water[idx];
                //double usq_water = u1w * u1w + u2w * u2w;

                //feq_water[0] = weights[0] * rhow * (1.0 - 1.5 * usq_water);
                //feq_water[1] = weights[1] * rhow * (1.0 + 3 * u1w + 4.5 * u1w * u1w - 1.5 * usq_water);
                //feq_water[2] = weights[2] * rhow * (1.0 + 3 * u2w + 4.5 * u2w * u2w - 1.5 * usq_water);
                //feq_water[3] = weights[3] * rhow * (1.0 - 3 * u1w + 4.5 * u1w * u1w - 1.5 * usq_water);
                //feq_water[4] = weights[4] * rhow * (1.0 - 3 * u2w + 4.5 * u2w * u2w - 1.5 * usq_water);
                //feq_water[5] = weights[5] * rhow * (1.0 + 3 * (u1w + u2w) + 4.5 * (u1w + u2w) * (u1w + u2w) - 1.5 * usq_water);
                //feq_water[6] = weights[6] * rhow * (1.0 + 3 * (-u1w + u2w) + 4.5 * (-u1w + u2w) * (-u1w + u2w) - 1.5 * usq_water);
                //feq_water[7] = weights[7] * rhow * (1.0 + 3 * (-u1w - u2w) + 4.5 * (u1w + u2w) * (u1w + u2w) - 1.5 * usq_water);
                //feq_water[8] = weights[8] * rhow * (1.0 + 3 * (u1w - u2w) + 4.5 * (u1w - u2w) * (u1w - u2w) - 1.5 * usq_water);

                //double rhoa = rho_air[idx];
                //double usq_air = u1a * u1a + u2a * u2a;

                //feq_air[0] = weights[0] * rhoa * (1.0 - 1.5 * usq_air);
                //feq_air[1] = weights[1] * rhoa * (1.0 + 3 * u1a + 4.5 * u1a * u1a - 1.5 * usq_air);
                //feq_air[2] = weights[2] * rhoa * (1.0 + 3 * u2a + 4.5 * u2a * u2a - 1.5 * usq_air);
                //feq_air[3] = weights[3] * rhoa * (1.0 - 3 * u1a + 4.5 * u1a * u1a - 1.5 * usq_air);
                //feq_air[4] = weights[4] * rhoa * (1.0 - 3 * u2a + 4.5 * u2a * u2a - 1.5 * usq_air);
                //feq_air[5] = weights[5] * rhoa * (1.0 + 3 * (u1a + u2a) + 4.5 * (u1a + u2a) * (u1a + u2a) - 1.5 * usq_air);
                //feq_air[6] = weights[6] * rhoa * (1.0 + 3 * (-u1a + u2a) + 4.5 * (-u1a + u2a) * (-u1a + u2a) - 1.5 * usq_air);
                //feq_air[7] = weights[7] * rhoa * (1.0 + 3 * (-u1a - u2a) + 4.5 * (u1a + u2a) * (u1a + u2a) - 1.5 * usq_air);
                //feq_air[8] = weights[8] * rhoa * (1.0 + 3 * (u1a - u2a) + 4.5 * (u1a - u2a) * (u1a - u2a) - 1.5 * usq_air);

                //for (int k = 0; k < npop; ++k)
                //{
                //    f_water[9 * idx + k] = f_water[9 * idx + k]
                //        - (f_water[9 * idx + k] - feq_water[k]) / tau_water
                //        + fguo_water[k];

                //    f_air[9 * idx + k] = f_air[9 * idx + k]
                //        - (f_air[9 * idx + k] - feq_air[k]) / tau_air
                //        + fguo_air[k];


                //}

                ////STREAMING 

                //for (int k = 0; k < npop; ++k)
                //{
                //    int iX2 = (iX + cx[k] + nx) % nx;
                //    int iY2 = (iY + cy[k] + ny) % ny;
                //    int idx2 = iY2 * nx + iX2;

                //    f_water_tmp[9 * idx2 + k] = f_water[9 * idx + k];
                //    f_air_tmp[9 * idx2 + k] = f_air[9 * idx + k];
                //}



                //REGULARIZED BGK COLLISION --------- WATER

                double rhow = rho_water[idx];
                double usq_water = u1w * u1w + u2w * u2w;

                feq_water[0] = weights[0] * rhow * (1.0 - 1.5 * usq_water);
                feq_water[1] = weights[1] * rhow * (1.0 + 3 * u1w + 4.5 * u1w * u1w - 1.5 * usq_water);
                feq_water[2] = weights[2] * rhow * (1.0 + 3 * u2w + 4.5 * u2w * u2w - 1.5 * usq_water);
                feq_water[3] = weights[3] * rhow * (1.0 - 3 * u1w + 4.5 * u1w * u1w - 1.5 * usq_water);
                feq_water[4] = weights[4] * rhow * (1.0 - 3 * u2w + 4.5 * u2w * u2w - 1.5 * usq_water);
                feq_water[5] = weights[5] * rhow * (1.0 + 3 * (u1w + u2w) + 4.5 * (u1w + u2w) * (u1w + u2w) - 1.5 * usq_water);
                feq_water[6] = weights[6] * rhow * (1.0 + 3 * (-u1w + u2w) + 4.5 * (-u1w + u2w) * (-u1w + u2w) - 1.5 * usq_water);
                feq_water[7] = weights[7] * rhow * (1.0 + 3 * (-u1w - u2w) + 4.5 * (u1w + u2w) * (u1w + u2w) - 1.5 * usq_water);
                feq_water[8] = weights[8] * rhow * (1.0 + 3 * (u1w - u2w) + 4.5 * (u1w - u2w) * (u1w - u2w) - 1.5 * usq_water);

                // Compute full non-equilibrium stress tensor
                double Pi_xx_w = 0.0, Pi_xy_w = 0.0, Pi_yy_w = 0.0;
                for (int k = 0; k < npop; ++k)
                {
                    double f_neq_w = f_water[9 * idx + k] - feq_water[k];
                    Pi_xx_w += cx[k] * cx[k] * f_neq_w;   // c_x * c_x
                    Pi_xy_w += cx[k] * cy[k] * f_neq_w;   // c_x * c_y
                    Pi_yy_w += cy[k] * cy[k] * f_neq_w;   // c_y * c_y
                }

                // Remove trace to make it traceless
                const double trace_w = (Pi_xx_w + Pi_yy_w) ;
                Pi_xx_w -= trace_w * 0.5 ;
                Pi_yy_w -= trace_w * 0.5 ;

                // Reconstruct regularized non-equilibrium
                double f_reg_water[npop];
                for (int k = 0; k < npop; ++k)
                {
                    double Hab_Pi_ab_w = (cx[k] * cx[k] - cs2) * Pi_xx_w
                        + 2.0 * cx[k] * cy[k] * Pi_xy_w
                        + (cy[k] * cy[k] - cs2) * Pi_yy_w;

                    f_reg_water[k] = feq_water[k]
                        + (1.0 - 1.0 / tau_water)
                        * weights[k] * Hab_Pi_ab_w / (2.0 * cs2 * cs2);
                    f_reg_water[k] += fguo_water[k];

                }



                //REGULARIZED BGK COLLISION --------- AIR

                double rhoa = rho_air[idx];
                double usq_air = u1a * u1a + u2a * u2a;

                // Compute equilibrium
                feq_air[0] = weights[0] * rhoa * (1.0 - 1.5 * usq_air);
                feq_air[1] = weights[1] * rhoa * (1.0 + 3 * u1a + 4.5 * u1a * u1a - 1.5 * usq_air);
                feq_air[2] = weights[2] * rhoa * (1.0 + 3 * u2a + 4.5 * u2a * u2a - 1.5 * usq_air);
                feq_air[3] = weights[3] * rhoa * (1.0 - 3 * u1a + 4.5 * u1a * u1a - 1.5 * usq_air);
                feq_air[4] = weights[4] * rhoa * (1.0 - 3 * u2a + 4.5 * u2a * u2a - 1.5 * usq_air);
                feq_air[5] = weights[5] * rhoa * (1.0 + 3 * (u1a + u2a) + 4.5 * (u1a + u2a) * (u1a + u2a) - 1.5 * usq_air);
                feq_air[6] = weights[6] * rhoa * (1.0 + 3 * (-u1a + u2a) + 4.5 * (-u1a + u2a) * (-u1a + u2a) - 1.5 * usq_air);
                feq_air[7] = weights[7] * rhoa * (1.0 + 3 * (-u1a - u2a) + 4.5 * (u1a + u2a) * (u1a + u2a) - 1.5 * usq_air);
                feq_air[8] = weights[8] * rhoa * (1.0 + 3 * (u1a - u2a) + 4.5 * (u1a - u2a) * (u1a - u2a) - 1.5 * usq_air);

                // Compute full non-equilibrium stress tensor
                double Pi_xx_a = 0.0, Pi_xy_a = 0.0, Pi_yy_a = 0.0;
                for (int k = 0; k < npop; ++k)
                {
                    double f_neq_a = f_air[9 * idx + k] - feq_air[k];
                    Pi_xx_a += cx[k] * cx[k] * f_neq_a;   // c_x * c_x
                    Pi_xy_a += cx[k] * cy[k] * f_neq_a;   // c_x * c_y
                    Pi_yy_a += cy[k] * cy[k] * f_neq_a;   // c_y * c_y
                }

                // Remove trace to make it traceless
                const double trace_a = (Pi_xx_a + Pi_yy_a);
                Pi_xx_a -= trace_a * 0.5;
                Pi_yy_a -= trace_a * 0.5;

                // Reconstruct regularized non-equilibrium
                double f_reg_air[npop];
                for (int k = 0; k < npop; ++k)
                {
                    double Hab_Pi_ab_a = (cx[k] * cx[k] - cs2) * Pi_xx_a
                        + 2.0 * cx[k] * cy[k] * Pi_xy_a
                        + (cy[k] * cy[k] - cs2) * Pi_yy_a;

                    f_reg_air[k] = feq_air[k]
                        + (1.0 - 1.0 / tau_air)
                        * weights[k] * Hab_Pi_ab_a / (2.0 * cs2 * cs2);
                    f_reg_air[k] += fguo_air[k];

                }


                //STREAMING 

                for (int k = 0; k < npop; ++k)
                {
                    int iX2 = (iX + cx[k] + nx) % nx;
                    int iY2 = (iY + cy[k] + ny) % ny;
                    int idx2 = iY2 * nx + iX2;

                    // Water streaming
                    f_water_tmp[9 * idx2 + k] = f_reg_water[k];

                    // Air streaming
                    f_air_tmp[9 * idx2 + k] = f_reg_air[k];
                }


            }
        }

        // Swap pointers
        std::swap(f_water, f_water_tmp);
        std::swap(f_air, f_air_tmp);


        if (t % noutput == 0) {
            // diretorio de execução (000 como um id)
            std::string SIM_DIR = "000";
            std::string ID = "000";

            // criar pasta de simulação se não existir
            _mkdir(SIM_DIR.c_str());

            saveToBinary(SIM_DIR, ID, "rho_water", t, rho_water, nx * ny);
            saveToBinary(SIM_DIR, ID, "rho_air",   t, rho_air,   nx * ny);
            saveToBinary(SIM_DIR, ID, "u1_water",  t, u1_water,  nx * ny);
            saveToBinary(SIM_DIR, ID, "u1_air",    t, u1_air,    nx * ny);

            std::cout << "Salvo passo " << t << " em " << SIM_DIR << std::endl;
        }

    }

    std::time_t finish = std::time(nullptr);
    std::cout << "Simulation finished in " << (finish - start)
        << " seconds" << std::endl;

    
generateSimulationInfoFile(
        "000/",        
        "000",        
        "D2Q9",     
        nx, ny,
        nsteps,
        noutput,
        tau_water,
        000.0 // placeholder para MLUPS
    );

    // Cleanup
    delete[] f_water_tmp;
    delete[] f_air_tmp;
    return 0;

}