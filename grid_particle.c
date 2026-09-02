#include "stdio.h"
#include "stdlib.h"
#include <math.h>
#include "time.h" // For random number seed

typedef struct {
    double x;           // x-coordinate
    double y;           // y-coordinate
    double v[2];        // Velocity [vx, vy]
    double r;           // Radius of the particle
    double theta;       // Angle of rotation
    double omega;       // Angular velocity
    double m;           // Mass of the particle
    double mom_I;       // Moment of inertia
} particle;

// Function to generate a random double in a given range
double random_double(double min, double max) {
    return min + (rand() / (double)RAND_MAX) * (max - min);
}

// Initialize grid-based particle placement
void initial_condition_r(int npoints_max, double radius, int *npoints, particle p[])
{
    double grid_multiplier = 10.0;
    double grid_min_x = -0.2 * grid_multiplier;  // Left bound of the grid
    double grid_max_x = 0.2 * grid_multiplier;   // Right bound of the grid
    double grid_min_y = 0.0 * grid_multiplier;   // Bottom bound of the grid
    double grid_max_y = 0.4 * grid_multiplier;
    double grid_size;    // Size of each grid cell
    int grid_cols,grid_rows;
    int k = 0;
    
    if (4*radius >= (grid_max_x - grid_min_x)) {
        printf("radius too large for x extend of box\n");
        exit(1);
    }
    if (4*radius >= (grid_max_y - grid_min_y)) {
        printf("radius too large for y extend of box\n");
        exit(1);
    }

    
    grid_size = 4.0 * radius;
    grid_cols = (int)( (grid_max_x - grid_min_x)/grid_size );
    grid_rows = (int)( (grid_max_y - grid_min_y)/grid_size );
    printf("cells % d x %d\n",grid_cols, grid_rows);
    
    for (int row = 0; row < grid_rows && k < npoints_max; row++) {
        for (int col = 0; col < grid_cols && k < npoints_max; col++) {
            double cell_min_x = grid_min_x + col * grid_size;
            double cell_max_x = cell_min_x + grid_size;
            double cell_min_y = grid_min_y + row * grid_size;
            double cell_max_y = cell_min_y + grid_size;

            // Randomly position the particle within the cell
            p[k].x = random_double(cell_min_x + radius, cell_max_x - radius);
            p[k].y = random_double(cell_min_y + radius, cell_max_y - radius);
            
            // Centered particles
            //p[k].x = cell_min_x + grid_size/2;
            //p[k].y = cell_min_y + grid_size/2;
            
            // Assign random velocities in range [-10, 10]
            p[k].v[0] = random_double(0.0, 0.0);
            p[k].v[1] = random_double(0.0, 0.0);

            // Set particle properties
            p[k].r = radius;  // Radius
            p[k].m = 0.2;   // Mass
            p[k].theta = random_double(0, 2 * M_PI);  // Random orientation
            p[k].omega = random_double(-1.0, 1.0);    // Random angular velocity
            p[k].mom_I = (1.0 / 2.0) * p[k].m * p[k].r * p[k].r; // Moment of inertia

            k++;
        }
    }
    *npoints = k;
}








void H_interaction(int npoints, const particle Y[], particle h[])
{
    // Normal Particle Interactions
    int i,j;
    double distance, delta, normal[2], F_n[2];
    double k_n = 100000;
    double c_n = 5000.0;
    
    for (i=0; i<npoints; i++) {
        // get particle i Y[i]
        for (j=0; j<npoints; j++) {
            
            if (i == j) {
                continue;
            }
            distance = sqrt((Y[j].x - Y[i].x)*(Y[j].x - Y[i].x) + (Y[j].y - Y[i].y)*(Y[j].y - Y[i].y));
            delta = (Y[i].r + Y[j].r) - distance;
            normal[0] = (Y[j].x - Y[i].x) /distance;
            normal[1] = (Y[j].y - Y[i].y) /distance;
            
            // chcek for overlap between i and j
            if (delta > 0) {
                //printf("in contact i %d | j \%d\n",i,j);
                //F_n[0] = 0;// fix me later
                //F_n[1] = k_n * overlap - c_n * Y[k].v[1];;

                // compute normal force from j on i
                double relative[2] = {Y[j].v[0] - Y[i].v[0] , Y[j].v[1] - Y[i].v[1]  };
                
                relative[0] += -1 * Y[j].omega * (Y[j].y - Y[i].y)/2.0;
                relative[1] +=      Y[j].omega * (Y[j].x - Y[i].x)/2.0;
                
                relative[0] +=        Y[i].omega * (Y[i].y - Y[j].y)/2.0;
                relative[1] += -1  *  Y[i].omega * (Y[i].x - Y[j].x)/2.0;
                
                double v_dot_n  = relative[0] * normal[0] + relative[1] * normal[1];
                double v_n[2];
                v_n[0] = v_dot_n * normal[0];
                v_n[1] = v_dot_n * normal[1];
               
                
                F_n[0] = -1 * (k_n * delta * normal[0] - c_n * v_n[0]);
                F_n[1] = -1 * (k_n * delta * normal[1] - c_n * v_n[1]);
                //printf("delta %+1.10e\n",delta);
                //printf("F_n[0] %+1.10e\n",F_n[0]);
                //printf("F_n[1] %+1.10e\n",F_n[1]);
                //printf("normal[0] %+1.10e\n",normal[0]);
                //printf("normal[1] %+1.10e\n",normal[1]);
                //printf("---\n");
                // add force to i from j
                h[i].v[0] += F_n[0] / Y[i].m;
                h[i].v[1] += F_n[1] / Y[i].m;

                // add force to j from i
                //xh[j].v[0] -= F_n[0] / Y[j].m;
                //h[j].v[1] -= F_n[1] / Y[j].m;

            }
        }
    }
}

// Contact with walls
void H_contact(int npoints, const particle Y[], double time, particle h[])
{
    double grid_multiplier = 10.0;
    int k;
    //Ceiling Boundary
    double y_max = 0.4 * grid_multiplier; //Ceiling Boundary
    // Wall Boundaries
    double x_min = -0.2 * grid_multiplier; // Left boundary
    double x_max = .2 * grid_multiplier;  // Right boundary
    

    for (k=0; k < npoints; k++){
        h[k].v[0] = 0.0;
        h[k].v[1] = 0.0;
        h[k].omega = 0.0;
    }
    
    // gravity and wall forces only
    for (k=0; k < npoints; k++){
        double g = 9.8; // Gravity
        double m = Y[k].m; // Mass of ball
        double mom_I =  Y[k].mom_I;
        double dt = .00010;
        const double *v = Y[k].v;
        //double r = 0.01;
        double k_n  = 50000.0; // Spring constant
        double k_t = 5;
        double c_n = 10000.0; // damping cpnstant
        double c_t = 50000.0;
        double F_ext[] = { 0.0, 0.0 };
        double torque_ext = 0;
        double F_n[2];
        double F_t[2];
        double r[2]= {0.0, -Y[k].r};
        double overlap = (Y[k].y - Y[k].r);
        double u_s = 0.8; //Friction Threshold
        double V_n[2];
        double V_t[2];
        double N_wall[] = {0.0, -1.0};
        double disp_t[2];
        
        
        //Compute Vectors
        
        V_n[0] = (v[0]*N_wall[0] + v[1]* N_wall[1]) * N_wall[0];
        V_n[1] = (v[0]*N_wall[0] + v[1]* N_wall[1]) * N_wall[1];
        
        V_t[0] = (v[0] - V_n[0]);
        V_t[1] = (v[1] - V_n[1]);
        
        disp_t[0] = V_t[0] * dt;
        disp_t[1] = V_t[1] * dt;
        
        //Computer Displacements
        
        //Floor Boundary
        // Check if the radius overlaps with ground
        F_ext[1] = -Y[k].m * g;
        if (overlap <= 0){
            
            overlap = overlap * -1;
            
            // Force Normal
            F_n[0] = 0;// fix me later
            F_n[1] = k_n * overlap - c_n * Y[k].v[1];;
            double Fn_magnitude = sqrt(F_n[0] * F_n[0] + F_n[1] * F_n[1]);
            
            
            // Force Tangential
            F_t[0] = -1 * k_t * disp_t[0] - c_t * Y[k].v[0];
            F_t[1] = 0;// fix me later
            double Ft_magnitude = sqrt(F_t[0] * F_t[0] + F_t[1] * F_t[1]);
            
            if (Ft_magnitude > Fn_magnitude * u_s){
                F_t[0] = F_t[0] * u_s * Fn_magnitude/Ft_magnitude;
                F_t[1] = F_t[1] * u_s * Fn_magnitude/Ft_magnitude;
            }
            // add forces into F_ext
            F_ext[0] += (F_n[0] + F_t[0]); // Force X
            F_ext[1] += (F_n[1] + F_t[1]); // Force Y
            
            // r cross ft
            double r_cross = (r[0] * F_t[1]) - (r[1] * F_t[0]);
            //r_cross = r[1] * F_t[0];

            // sum torques into torque_ext
            torque_ext += r_cross;
            
        }
        //Ceiling Boundar
        // Check for ceiling collision
        if (Y[k].y + Y[k].r >= y_max) {
            double overlap = (Y[k].y + Y[k].r) - y_max;
            overlap = overlap;

            // Force Normal
            F_n[0] = 0;
            F_n[1] = -k_n * overlap - c_n * Y[k].v[1];
            double Fn_magnitude = fabs(F_n[1]);
           
            F_t[0] = -k_t * disp_t[0] - c_t * Y[k].v[0];
            F_t[1] = 0;
            double Ft_magnitude = fabs(F_t[0]);

            if (Ft_magnitude > Fn_magnitude * u_s) {
                F_t[0] = F_t[0] * u_s * Fn_magnitude / Ft_magnitude;
            }

            // Add forces into F_ext
            F_ext[0] += F_t[0];
            F_ext[1] += F_n[1];

            // Calculate torque (r cross F_t)
            double r_cross = (r[0] * F_t[1]) - (r[1] * F_t[0]);
            // Sum torque into torque_ext
            torque_ext += r_cross;
        }

        

        // Check for left boundary collision
        if (Y[k].x - Y[k].r <= x_min) {
            double overlap = x_min - (Y[k].x - Y[k].r); // Overlap with the left wall
            overlap = overlap;

            // Force Normal
            F_n[0] = k_n * overlap - c_n * Y[k].v[0];
            F_n[1] = 0;
            double Fn_magnitude = fabs(F_n[0]);

            // Force Tangential
            F_t[0] = 0;
            F_t[1] = -1 * k_t * disp_t[1] - c_t * Y[k].v[1];
            double Ft_magnitude = fabs(F_t[1]);

            // Friction threshold
            if (Ft_magnitude > Fn_magnitude * u_s) {
                F_t[1] = F_t[1] * u_s * Fn_magnitude / Ft_magnitude;
            }

            // Add forces into F_ext
            F_ext[0] += F_n[0]; // Normal force in x
            F_ext[1] += F_t[1]; // Tangential force in y
        }

        // Right Wall
        if (Y[k].x + Y[k].r >= x_max) {
            double overlap = (Y[k].x + Y[k].r) - x_max;
            overlap = overlap;

            // Force Normal
            F_n[0] = -k_n * overlap - c_n * Y[k].v[0];
            F_n[1] = 0;
            double Fn_magnitude = F_n[0];

            // Force Tangential
            F_t[0] = 0;
            F_t[1] = -1 * k_t * disp_t[1] - c_t * Y[k].v[1];
            double Ft_magnitude = F_t[1];

            // Friction threshold
            if (Ft_magnitude > Fn_magnitude * u_s) {
                F_t[1] = F_t[1] * u_s * Fn_magnitude / Ft_magnitude;
            }

            // Add forces into F_ext
            F_ext[0] += F_n[0]; // Normal force in x
            F_ext[1] += F_t[1]; // Tangential force in y
        }

        
        
        
        //Partivle update
        h[k].v[0]  += F_ext[0] / Y[k].m ; // Velocity in the x diretion
        h[k].v[1]  += F_ext[1] / Y[k].m; // Adds acceleratoin?
        h[k].omega += torque_ext/ Y[k].mom_I;
        
        //Position update
        h[k].x = Y[k].v[0];
        h[k].y = Y[k].v[1];
        h[k].theta = Y[k].omega;
        
    }
    
    // particle-particle forces
    // h[i] += force_from_particle_j
    H_interaction(npoints, Y, h);
    
}





int main(int nargs, char *args[])
{
    srand(time(NULL));
    double dt = 0.00001;
    int count = 0;
    double time = 0;
    FILE *fp;
    int steps = 0;
    

    double time_max = 1;
    // RHS = rate of change
    particle *ynew, *yold;
    particle *rhs;
    int k, npoints_max = 100, npoints = 1000;
    
    ynew = (particle*)malloc(sizeof(particle)*npoints);
    yold = (particle*)malloc(sizeof(particle)*npoints);
    rhs = (particle*)malloc(sizeof(particle)*npoints);

   // initial_condition(npoints, yold);
    printf("init np_max %d np %d\n", npoints_max, npoints);
    initial_condition_r(npoints_max, 0.0125 * 10.0, &npoints, yold);
    printf("finalized init np_max %d np %d\n", npoints_max, npoints);
    //exit(1);

    fp = fopen("multi_trajectory.txt", "w"); // "w" => write, "r" => read
    
    fprintf(fp, "# time radius x y vx vy theta omega\n");

    //Time Integrstor
    while (time <= time_max){

        time += dt;
        
        H_contact(npoints, yold, time, rhs);
        
        // Position changes
        // Logic => change in position / time * time =  change in position
        for (k=0; k<npoints; k++) {
            ynew[k].x = yold[k].x + dt * rhs[k].x;
            ynew[k].y = yold[k].y + dt * rhs[k].y;
            
            ynew[k].v[0] = yold[k].v[0] + dt * rhs[k].v[0];
            ynew[k].v[1] = yold[k].v[1] + dt * rhs[k].v[1];
            
            ynew[k].theta = yold[k].theta + dt * rhs[k].theta;
            ynew[k].omega = yold[k].omega + dt * rhs[k].omega;
            
            ynew[k].m = yold[k].m;
            ynew[k].mom_I = yold[k].mom_I;
            
            ynew[k].r = yold[k].r;
        }

        if (steps % 1000 == 0) {
            for (k=0; k<npoints; k++) {
                fprintf(fp, "%lf %lf %lf %lf %lf %lf %lf %lf\n", time, ynew[k].r, ynew[k].x, ynew[k].y, ynew[k].v[0], ynew[k].v[1], ynew[k].theta, ynew[k].omega);
            }
        }
        
        for (k=0; k<npoints; k++) {
            yold[k] = ynew[k];
        }
        
        steps++;
    }

    fclose(fp);
    return 0;
    
}
